/*!
 * \file crack_generator.cpp
 * Author: Stephan Grein
 * Created on: 09/18/2018
 */

#include "crack_generator.h"
#include "lib_grid/lib_grid.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"
#include "lib_grid/algorithms/remeshing/delaunay_triangulation.h"
#include "lib_grid/refinement/regular_refinement.h"
#include <cmath>

#define UG_ENABLE_WARNINGS

namespace ug {
	namespace crack_generator {
		////////////////////////////////////////////////////////////////////////////////
		/// BuildCompleteCrack
		////////////////////////////////////////////////////////////////////////////////
		void BuildCompleteCrack
		(
			number crackInnerLength=0.2,
			number innerThickness=0.1,
			number crackOuterLength=2.0,
			number angle = 10
		)
		{
			/// grid management
			Grid g;
			SubsetHandler sh(g);
			Selector sel(g);
			sh.set_default_subset_index(0);

			AInt aInt;
			g.attach_to_vertices(aPosition);
			g.attach_to_vertices(aInt);

			Grid::VertexAttachmentAccessor<AInt> aaIntVertex(g, aInt);
			Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

			/// create crack tip and base and connect with edges
			Vertex* crackTipVtx = *g.create<RegularVertex>();
			aaPos[crackTipVtx] = vector3(0, 0, 0); /// crack tip

			vector3 crackBaseTop, crackBaseBottom;
			crackBaseTop.x() = -innerThickness * cos(deg_to_rad(angle));
			crackBaseTop.y() = -innerThickness * sin(deg_to_rad(angle));
			crackBaseTop.z() = 0;

			crackBaseBottom.x() = -innerThickness * cos(-deg_to_rad(angle));
			crackBaseBottom.y() = -innerThickness * sin(-deg_to_rad(angle));
			crackBaseBottom.z() = 0;

			Vertex* crackBaseTopVtx = *g.create<RegularVertex>();
			Vertex* crackBaseBottomVtx = *g.create<RegularVertex>();
			aaPos[crackBaseTopVtx] = crackBaseTop;
			aaPos[crackBaseBottomVtx] = crackBaseBottom;

			*g.create<RegularEdge>(EdgeDescriptor(crackTipVtx, crackBaseTopVtx));
			*g.create<RegularEdge>(EdgeDescriptor(crackTipVtx, crackBaseBottomVtx));
			sh.set_default_subset_index(1);

			crackBaseTop.x() = -crackInnerLength * cos(deg_to_rad(angle));
			crackBaseTop.y() = -crackInnerLength * sin(deg_to_rad(angle));
			crackBaseTop.z() = 0;
			crackBaseBottom.x() = -crackInnerLength * cos(-deg_to_rad(angle));
			crackBaseBottom.y() = -crackInnerLength * sin(-deg_to_rad(angle));
			crackBaseBottom.z() = 0;
			Vertex* v3 = *g.create<RegularVertex>();
			Vertex* v4 = *g.create<RegularVertex>();
			aaPos[v3] = crackBaseTop;
			aaPos[v4] = crackBaseBottom;
			*g.create<RegularEdge>(EdgeDescriptor(crackBaseTopVtx, v3));
			*g.create<RegularEdge>(EdgeDescriptor(crackBaseBottomVtx, v4));

			sh.set_default_subset_index(2);
			crackBaseTop.x() = -crackOuterLength * cos(deg_to_rad(angle));
			crackBaseTop.y() = -crackOuterLength * sin(deg_to_rad(angle));
			crackBaseTop.z() = 0;
			crackBaseBottom.x() = -crackOuterLength * cos(-deg_to_rad(angle));
			crackBaseBottom.y() = -crackOuterLength * sin(-deg_to_rad(angle));
			crackBaseBottom.z() = 0;

			Vertex* v5 = *g.create<RegularVertex>();
			Vertex* v6 = *g.create<RegularVertex>();
			aaPos[v5] = crackBaseTop;
			aaPos[v6] = crackBaseBottom;
			*g.create<RegularEdge>(EdgeDescriptor(v3, v5));
			*g.create<RegularEdge>(EdgeDescriptor(v4, v6));

			AssignSubsetColors(sh);
			SaveGridToFile(g, sh, "crack_generator_step_1.ugx");
			sh.set_default_subset_index(3);

			/// outermost square
			number outerDistance = VecDistance(crackBaseTop, crackBaseBottom);
			vector3 centerOuter = crackBaseTop;
			centerOuter.y() = centerOuter.y() + outerDistance / 2;
			number squareOuterDiameter = VecDistance(centerOuter, aaPos[crackTipVtx]);
			number outerHeight = squareOuterDiameter;
			outerHeight = 2*squareOuterDiameter - outerDistance;

			vector3 topLeft, bottomLeft, topRight, bottomRight;
			topLeft = crackBaseTop;
			topLeft.y() = topLeft.y() - outerHeight/2;

			bottomLeft = crackBaseBottom;
			bottomLeft.y() = bottomLeft.y() + outerHeight/2;
			Vertex* bottomLeftVtx = *g.create<RegularVertex>();
			Vertex* bottomRightVtx = *g.create<RegularVertex>();
			aaPos[bottomLeftVtx] = topLeft;
			aaPos[bottomRightVtx] = bottomLeft;
			*g.create<RegularEdge>(EdgeDescriptor(v5, bottomLeftVtx));
			*g.create<RegularEdge>(EdgeDescriptor(v6, bottomRightVtx));

			topRight = topLeft;
			topRight.x() = topRight.x() + 2.0*squareOuterDiameter;
			bottomRight = bottomLeft;
			bottomRight.x() = bottomRight.x() + 2.0*squareOuterDiameter;

			sh.set_default_subset_index(7);
			Vertex* topRightVtx = *g.create<RegularVertex>();
			Vertex* v10 = *g.create<RegularVertex>();
			aaPos[topRightVtx] = topRight;
			aaPos[v10] = bottomRight;

			sh.set_default_subset_index(5);
			*g.create<RegularEdge>(EdgeDescriptor(bottomLeftVtx, topRightVtx));
			sh.set_default_subset_index(6);
			*g.create<RegularEdge>(EdgeDescriptor(bottomRightVtx, v10));
			sh.set_default_subset_index(7);
			*g.create<RegularEdge>(EdgeDescriptor(topRightVtx, v10));

			AssignSubsetColors(sh);
			sh.set_default_subset_index(0);
			SaveGridToFile(g, sh, "crack_generator_step_2.ugx");

			/// innermost square
			number innerDistance = VecDistance(aaPos[crackBaseTopVtx], aaPos[crackBaseBottomVtx]);

			vector3 centerInner = aaPos[crackBaseTopVtx];
			centerInner.y() = centerInner.y() + innerDistance / 2;
			number squareInnerDiameter = VecDistance(centerInner, aaPos[crackTipVtx]);
			number innerHeight = squareInnerDiameter;
			innerHeight = 2*squareInnerDiameter - innerDistance;
			UG_COND_THROW(std::signbit(innerHeight), "inner height cannot be negative");

			topLeft = aaPos[crackBaseTopVtx];
			topLeft.y() = topLeft.y() - innerHeight/2;
			bottomLeft = aaPos[crackBaseBottomVtx];
			bottomLeft.y() = bottomLeft.y() + innerHeight/2;
			Vertex* v11 = *g.create<RegularVertex>();
			Vertex* v12 = *g.create<RegularVertex>();
			aaPos[v11] = topLeft;
			aaPos[v12] = bottomLeft;
			*g.create<RegularEdge>(EdgeDescriptor(crackBaseTopVtx, v11));
			*g.create<RegularEdge>(EdgeDescriptor(crackBaseBottomVtx, v12));

			topRight = topLeft;
			topRight.x() = topRight.x() + 2*squareInnerDiameter;
			bottomRight = bottomLeft;
			bottomRight.x() = bottomRight.x() + 2*squareInnerDiameter;
			Vertex* v13 = *g.create<RegularVertex>();
			Vertex* v14 = *g.create<RegularVertex>();
			aaPos[v13] = topRight;
			aaPos[v14] = bottomRight;
			*g.create<RegularEdge>(EdgeDescriptor(v11, v13));
			*g.create<RegularEdge>(EdgeDescriptor(v12, v14));
			*g.create<RegularEdge>(EdgeDescriptor(v13, v14));

			AssignSubsetColors(sh);
			SaveGridToFile(g, sh, "crack_generator_step_3.ugx");

			/// middle square (Refine this square)
			sh.set_default_subset_index(1);
			number middleDistance = VecDistance(aaPos[v3], aaPos[v4]);

			vector3 centerMiddle = aaPos[v3];
			centerMiddle.y() = centerMiddle.y() + middleDistance / 2;
			number squareMiddleDiameter = VecDistance(centerMiddle, aaPos[crackTipVtx]);
			number middleHeight = squareMiddleDiameter;
			middleHeight = 2*squareMiddleDiameter - middleDistance;
			UG_COND_THROW(std::signbit(middleHeight), "middle height cannot be negative");

			topLeft = aaPos[v3];
			topLeft.y() = topLeft.y() - middleHeight/2;
			bottomLeft = aaPos[v4];
			bottomLeft.y() = bottomLeft.y() + middleHeight/2;
			Vertex* v15 = *g.create<RegularVertex>();
			Vertex* v16 = *g.create<RegularVertex>();
			aaPos[v15] = topLeft;
			aaPos[v16] = bottomLeft;
			*g.create<RegularEdge>(EdgeDescriptor(v3, v15));
			*g.create<RegularEdge>(EdgeDescriptor(v4, v16));

			topRight = topLeft;
			topRight.x() = topRight.x() + 2.0*squareMiddleDiameter;
			bottomRight = bottomLeft;
			bottomRight.x() = bottomRight.x() + 2.0*squareMiddleDiameter;
			Vertex* v17 = *g.create<RegularVertex>();
			Vertex* v18 = *g.create<RegularVertex>();
			aaPos[v17] = topRight;
			aaPos[v18] = bottomRight;
			*g.create<RegularEdge>(EdgeDescriptor(v15, v17));
			*g.create<RegularEdge>(EdgeDescriptor(v16, v18));
			*g.create<RegularEdge>(EdgeDescriptor(v17, v18));

			/// Rename subsets
			EraseEmptySubsets(sh);
			sh.subset_info(0).name = "Inner square";
			sh.subset_info(1).name = "Middle square";
			sh.subset_info(2).name = "Outer square";
			sh.subset_info(3).name = "Left boundary";
			sh.subset_info(4).name = "Back boundary";
			sh.subset_info(5).name = "Front boundary";
			sh.subset_info(6).name = "Right boundary";
			AssignSubsetColors(sh);
			SaveGridToFile(g, sh, "crack_generator_step_4.ugx");

			/// Triangulate bottom surface
			for (int i = 0; i < sh.num_subsets(); i++) {
				SelectSubsetElements<Edge>(sel, sh, i, true);
			}
			TriangleFill_SweepLine(g, sel.edges_begin(), sel.edges_end(), aPosition, aInt, &sh, 7);
			SelectSubsetElements<Face>(sel, sh, 7, true);
			QualityGridGeneration(g, sel.faces_begin(), sel.faces_end(), aaPos, 10);
			AssignSubsetColors(sh);
			SaveGridToFile(g, sh, "crack_generator_step_5.ugx");

			/// Extrude towards top
			vector3 normal = ug::vector3(0, 0, 2*squareOuterDiameter);
			std::vector<Edge*> edges;
			for (int i = 0; i < sh.num_subsets(); i++) {
				SelectSubsetElements<Edge>(sel, sh, i, true);
			}
			edges.assign(sel.edges_begin(), sel.edges_end());
			Extrude(g, NULL, &edges, NULL, normal, aaPos, EO_CREATE_FACES, NULL);
			sh.subset_info(7).name = "Bottom boundary";

			/// Triangulate top surface
			TriangleFill_SweepLine(g, edges.begin(), edges.end(), aPosition, aInt, &sh, 8);
			QualityGridGeneration(g, sh.begin<Face>(8), sh.end<Face>(8), aaPos, 30.0);
			sh.subset_info(8).name = "Top boundary";
			AssignSubsetColors(sh);
			SaveGridToFile(g, sh, "crack_generator_step_6.ugx");

			/// Tetrahedralize whole grid
			Tetrahedralize(g, 5, false, false, aPosition, 1);
			AssignSubsetColors(sh);
			SaveGridToFile(g, sh, "crack_generator_step_7.ugx");
		}

	////////////////////////////////////////////////////////////////////////////////
	/// CREATE_RECT
	////////////////////////////////////////////////////////////////////////////////
	void create_rect
	(
		vector3 bottomLeft,
		vector3 bottomRight,
		vector3 leftMDLayer,
		vector3 rightMDLayer,
		vector3 topLeft,
		vector3 topRight,
		Grid& g,
		SubsetHandler& sh,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		AInt& aInt,
		number h_r_0,
		Selector& sel,
		number depth,
		size_t si_offset,
		std::vector<Vertex*>& verts
	)
	{
		std::stringstream step;

		sh.set_default_subset_index(si_offset);
		Vertex* bottomLeftVertex = *g.create<RegularVertex>();
		aaPos[bottomLeftVertex] = bottomLeft;
		Vertex* bottomRightVertex = *g.create<RegularVertex>();
		aaPos[bottomRightVertex] = bottomRight;
		Edge* e1 = *g.create<RegularEdge>(EdgeDescriptor(bottomLeftVertex, bottomRightVertex));
		AssignSubsetColors(sh);
		step << "crack_generator_simple_step_" << si_offset+1 << ".ugx";
		SaveGridToFile(g, sh, step.str().c_str());
		step.str(""); step.clear();

		Vertex* leftMDLayerVertex = *g.create<RegularVertex>();
		aaPos[leftMDLayerVertex] = leftMDLayer;
		Vertex* rightMDLayerVertex = *g.create<RegularVertex>();
		aaPos[rightMDLayerVertex] = rightMDLayer;
		Edge* e2 = *g.create<RegularEdge>(EdgeDescriptor(bottomLeftVertex, leftMDLayerVertex));
		Edge* e3 = *g.create<RegularEdge>(EdgeDescriptor(bottomRightVertex, rightMDLayerVertex));
		AssignSubsetColors(sh);
		step << "crack_generator_simple_step_" << si_offset+2 << ".ugx";
		SaveGridToFile(g, sh, step.str().c_str());
		step.str(""); step.clear();

		sh.set_default_subset_index(1+si_offset);
		Vertex* topLeftVertex = *g.create<RegularVertex>();
		aaPos[topLeftVertex] = topLeft;
		Vertex* topRightVertex = *g.create<RegularVertex>();
		aaPos[topRightVertex] = topRight;
		Edge* e4 = *g.create<RegularEdge>(EdgeDescriptor(leftMDLayerVertex, topLeftVertex));
		AssignSubsetColors(sh);
		step << "crack_generator_simple_step_" << si_offset+3 << ".ugx";
		SaveGridToFile(g, sh, step.str().c_str());
		step.str(""); step.clear();

		Edge* e5 = *g.create<RegularEdge>(EdgeDescriptor(rightMDLayerVertex, topRightVertex));
		Edge* e6 = *g.create<RegularEdge>(EdgeDescriptor(topLeftVertex, topRightVertex));
		number leftMDLayertoRightMDLayer = VecDistance(leftMDLayer, rightMDLayer);
		number pos = 0;
		std::vector<Vertex*> vertices;
		vertices.push_back(topLeftVertex);
		/// REFINE HORIZONTAL
		while (pos < leftMDLayertoRightMDLayer) {
			Vertex* v3 = *g.create<RegularVertex>();
			vector3 temp3 = topLeft;
			vector3 dir;
			VecSubtract(dir, rightMDLayer, leftMDLayer);
			VecNormalize(dir, dir);
			VecScaleAdd(temp3, 1, temp3, pos, dir);
			aaPos[v3] = temp3;
			vertices.push_back(v3);
			pos+=h_r_0;
		}
		vertices.push_back(topRightVertex);
		for (size_t i = 0; i < vertices.size()-1; i++) {
			*g.create<RegularEdge>(EdgeDescriptor(vertices[i], vertices[i+1]));
		}
		vertices.clear();

		/// REFINE VERTICAL
		number leftMDLayerToTopLeft = VecDistance(leftMDLayer, topLeft);
		pos = 0;
		sel.clear();
		std::vector<Vertex*> vertices2;
		std::vector<Vertex*> vertices3;
		vertices.push_back(leftMDLayerVertex);
		vertices2.push_back(rightMDLayerVertex);
		while (pos < leftMDLayerToTopLeft) {
			pos+=h_r_0;
			Vertex* v1 = *g.create<RegularVertex>();
			Vertex* v2 = *g.create<RegularVertex>();
			vector3 temp1 = leftMDLayer;
			vector3 temp2 = rightMDLayer;
			vector3 dir;
			VecSubtract(dir, topLeft, leftMDLayer);
			VecNormalize(dir, dir);
			VecScaleAdd(temp1, 1, temp1, pos, dir);
			VecScaleAdd(temp2, 1, temp2, pos, dir);
			aaPos[v1] = temp1;
			aaPos[v2] = temp2;
			number pos2=0;
			vertices3.push_back(v1);
			while (pos2 < leftMDLayertoRightMDLayer) {
				vector3 temp3;
				vector3 dir;
				VecSubtract(dir, rightMDLayer, leftMDLayer);
				VecNormalize(dir, dir);
				Vertex* v3 = *g.create<RegularVertex>();
				VecScaleAdd(temp3, 1, temp1, pos2, dir);
				aaPos[v3] = temp3;
				vertices3.push_back(v3);
				pos2+=h_r_0;
			}
			vertices3.push_back(v2);

			for (size_t i = 0; i < vertices3.size()-1; i++) {
				*g.create<RegularEdge>(EdgeDescriptor(vertices3[i], vertices3[i+1]));
			}

			vertices3.clear();
			vertices.push_back(v1);
			vertices2.push_back(v2);
		}
		vertices.push_back(topLeftVertex);
		vertices2.push_back(topRightVertex);


		for (size_t i = 0; i < vertices.size()-1; i++) {
			*g.create<RegularEdge>(EdgeDescriptor(vertices[i], vertices[i+1]));
			*g.create<RegularEdge>(EdgeDescriptor(vertices2[i], vertices2[i+1]));
		}
		vertices.clear();
		vertices2.clear();

		sh.set_default_subset_index(si_offset);
		Edge* e7 = *g.create<RegularEdge>(EdgeDescriptor(leftMDLayerVertex, rightMDLayerVertex));

		pos = 0;
		vertices.push_back(leftMDLayerVertex);
		vertices2.push_back(bottomLeftVertex);
		while (pos < leftMDLayertoRightMDLayer) {
			Vertex* v1 = *g.create<RegularVertex>();
			Vertex* v2 = *g.create<RegularVertex>();
			vector3 temp1 = leftMDLayer;
			vector3 temp2 = bottomLeft;
			vector3 dir;
			VecSubtract(dir, rightMDLayer, leftMDLayer);
			VecNormalize(dir, dir);
			VecScaleAdd(temp1, 1, temp1, pos, dir);
			VecScaleAdd(temp2, 1, temp2, pos, dir);
			aaPos[v1] = temp1;
			aaPos[v2] = temp2;
			vertices.push_back(v1);
			vertices2.push_back(v2);
			pos+=h_r_0;
		}

		/// REFINE BD DOMAIN VERTICAL AND HORIZONTAL
		pos = 0;
		number distance = VecDistance(leftMDLayer, bottomLeft);
		std::vector<Vertex*> vertices4;
		std::vector<Vertex*> vertices5;
		std::vector<Vertex*> vertices6;
		vertices4.push_back(bottomLeftVertex);
		vertices5.push_back(bottomRightVertex);
		while (pos < distance) {
			Vertex* v3 = *g.create<RegularVertex>();
			Vertex* v4 = *g.create<RegularVertex>();
			vector3 temp3 = bottomLeft;
			vector3 temp4 = bottomRight;
			vector3 dir;
			VecSubtract(dir, leftMDLayer, bottomLeft);
			VecNormalize(dir, dir);
			VecScaleAdd(temp3, 1, temp3, pos, dir);
			VecScaleAdd(temp4, 1, temp4, pos, dir);
			aaPos[v3] = temp3;
			aaPos[v4] = temp4;
			vertices4.push_back(v3);
			vertices5.push_back(v4);
			pos+=h_r_0;
			number dist = VecDistance(leftMDLayer, rightMDLayer);
			number pos2 = 0;
			vertices6.push_back(v3);
			while (pos2 < dist) {
				Vertex* v5 = *g.create<RegularVertex>();
				vector3 temp5 = temp3;
				vector3 dir;
				VecSubtract(dir, rightMDLayer, leftMDLayer);
				VecNormalize(dir, dir);
				VecScaleAdd(temp5, 1, temp5, pos2, dir);
				aaPos[v5] = temp5;
				pos2+=h_r_0;
				vertices6.push_back(v5);
			}
			vertices6.push_back(v4);

			for (size_t i = 0; i < vertices6.size()-1; i++) {
				*g.create<RegularEdge>(EdgeDescriptor(vertices6[i], vertices6[i+1]));
			}
			vertices6.clear();
		}
		vertices4.push_back(leftMDLayerVertex);
		vertices5.push_back(rightMDLayerVertex);

		for (size_t i = 0; i < vertices4.size()-1; i++) {
			*g.create<RegularEdge>(EdgeDescriptor(vertices4[i], vertices4[i+1]));
			*g.create<RegularEdge>(EdgeDescriptor(vertices5[i], vertices5[i+1]));
		}

		vertices.push_back(rightMDLayerVertex);
		vertices2.push_back(bottomRightVertex);

		for (size_t i = 0; i < vertices.size()-1; i++) {
			*g.create<RegularEdge>(EdgeDescriptor(vertices[i], vertices[i+1]));
		}
		vertices.clear();
		vertices2.clear();
		sel.clear();
		sh.set_default_subset_index(si_offset);

		/// erase the old coarse edges
		g.erase(e1); g.erase(e2); g.erase(e3);
		g.erase(e4); g.erase(e5); g.erase(e6);
		g.erase(e7);

		AssignSubsetColors(sh);
		step << "crack_generator_simple_step_" << si_offset+4 << ".ugx";
		SaveGridToFile(g, sh, step.str().c_str());
		step.str(""); step.clear();

		verts.push_back(bottomLeftVertex);
		verts.push_back(bottomRightVertex);
	}

	////////////////////////////////////////////////////////////////////////////////
	/// BuildSimpleCrack
	////////////////////////////////////////////////////////////////////////////////
	void BuildSimpleCrack
	(
		number height,
		number width,
		number depth,
		number thickness,
		number spacing,
		number r_0,
		number h
	)
	{
		/// check user input
		UG_COND_THROW(thickness >= height || thickness >= width || thickness >= depth,
				"Thickness of bridging domain layers can't be larger then height of whole geometry.");
		if (fabs(fmod(width, h*r_0)) > SMALL || fabs(fmod(height, h*r_0)) > SMALL || fabs(fmod(depth, h*r_0)) > SMALL) {
			UG_LOGN("Width, height, or depth not evenly divisable by h*r_0, expect on borders non uniform spacing.")
		}

		/// grid management
		Grid g;
	    SubsetHandler sh(g);
	    sh.set_default_subset_index(0);
	    Selector sel(g);
	    AInt aInt;
	    g.attach_to_vertices(aPosition);
	    g.attach_to_vertices(aInt);
	    Grid::VertexAttachmentAccessor<AInt> aaIntVertex(g, aInt);
	    Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);

	    //// first (upper) rectangle
   		vector3 bottomLeft = ug::vector3(0, 0, 0);
		vector3 bottomRight = ug::vector3(width, 0, 0);
		vector3 leftMDLayer = ug::vector3(0, thickness, 0);
		vector3 rightMDLayer = ug::vector3(width, thickness, 0);
		vector3 topLeft = ug::vector3(0, height, 0);
		vector3 topRight = ug::vector3(width, height, 0);
		std::vector<std::pair<vector3, ug::vector3> > boxes;
		boxes.push_back(std::make_pair(leftMDLayer, topRight));
		boxes.push_back(std::make_pair(bottomLeft, rightMDLayer));
		size_t si_offset = 0;
		std::vector<Vertex*> verts;
		create_rect(bottomLeft, bottomRight, leftMDLayer, rightMDLayer, topLeft, topRight, g, sh, aaPos, aInt, h*r_0, sel, depth, si_offset, verts);
		si_offset = 3;

		/// Second (lower) rectangle
   		bottomLeft = vector3(0, -spacing, 0);
		bottomRight = vector3(width, -spacing, 0);
		leftMDLayer = vector3(0, -spacing-thickness, 0);
		rightMDLayer = vector3(width, -spacing-thickness, 0);
		topLeft = vector3(0, -spacing-height, 0);
		topRight = vector3(width, -spacing-height, 0);
		boxes.push_back(std::make_pair(topLeft, rightMDLayer));
		boxes.push_back(std::make_pair(leftMDLayer, bottomRight));
		create_rect(bottomLeft, bottomRight, leftMDLayer, rightMDLayer, topLeft, topRight, g, sh, aaPos, aInt, h*r_0, sel, depth, si_offset, verts);

		/// Connect the lower and upper rectangle
		sh.set_default_subset_index(2*si_offset);
		Edge* e1 = *g.create<RegularEdge>(EdgeDescriptor(verts[0], verts[2]));
		Edge* e2 = *g.create<RegularEdge>(EdgeDescriptor(verts[1], verts[3]));
		boxes.push_back(std::make_pair(vector3(0, -spacing, 0), ug::vector3(width, 0, 0)));

		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_8.ugx");
	    RemoveDoubles<3>(g, g.begin<Vertex>(), g.end<Vertex>(), aaPos, 0.0001);

		/// Triangulate bottom
	    /// TODO: Triangulate bottom manually by hand to achieve optimal uniform grid orientation
		for (int i = 0; i < sh.num_subsets(); i++) {
			SelectSubsetElements<Edge>(sel, sh, i, true);
		}
		TriangleFill_SweepLine(g, sel.edges_begin(), sel.edges_end(), aPosition, aInt, &sh, sh.num_subsets());
		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_9.ugx");

		/// Reassign the elements in the layers to subsets -> start beyond the current subsets, thus this ordering is the same as below
		size_t siFaces = sh.num_subsets()-1;
		for (size_t i = 0; i < boxes.size(); i++) {
			sel.clear();
			SelectSubsetElements<Face>(sel, sh, siFaces, true);
			Selector::traits<Face>::iterator fit = sel.faces_begin();
			Selector::traits<Face>::iterator fit_end = sel.faces_end();
			vector3 min, max;
			min = boxes[i].first;
			max = boxes[i].second;
			size_t si = sh.num_subsets();
			Selector sel2(g);
			for (; fit != fit_end; ++fit) {
				if(BoxBoundProbe(CalculateCenter(*fit, aaPos), min, max)) {
					sel2.select(*fit);
				}
			}
			CloseSelection(sel2);
			AssignSelectionToSubset(sel2, sh, si);
			sel2.clear();
		}
		EraseEmptySubsets(sh);
		AssignSubsetColors(sh);
		sel.clear();
		UG_LOGN("Triangulate bottom surface...")
		/// Retriangulate all
		for (int i = 0; i < sh.num_subsets(); i++) {
			SelectSubsetElements<Face>(sel, sh, i, true);
		}
		QualityGridGeneration(g, sel.faces_begin(), sel.faces_end(), aaPos, 30);
		SaveGridToFile(g, sh, "crack_generator_simple_step_10.ugx");

		/// Extrude all in steps to ensure uniformity we use h*r_0
		vector3 normal = ug::vector3(0, 0, depth);
		std::vector<Edge*> edges;
		for (int i = 0; i < sh.num_subsets(); i++) {
			SelectSubsetElements<Edge>(sel, sh, i, true);
		}
		edges.assign(sel.edges_begin(), sel.edges_end());
		VecScale(normal, normal, 0.5/h*r_0);
		number totalLength = normal.z();
		UG_LOGN("Extruding...")
		while (totalLength < depth) {
			Extrude(g, NULL, &edges, NULL, normal, aaPos, EO_CREATE_FACES, NULL);
			totalLength += normal.z();
		}
		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_11.ugx");

		/// Triangulate top
		UG_LOGN("Triangulate top surface...")
		TriangleFill_SweepLine(g, edges.begin(), edges.end(), aPosition, aInt, &sh, sh.num_subsets());
		AssignSubsetColors(sh);
		EraseEmptySubsets(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_12.ugx");

		siFaces = sh.num_subsets()-1;
		/// Reassign the elements in the layers to subsets (uses ordering from above)
		for (size_t i = 0; i < boxes.size(); i++) {
			sel.clear();
			SelectSubsetElements<Face>(sel, sh, siFaces, true);
			Selector::traits<Face>::iterator fit = sel.faces_begin();
			Selector::traits<Face>::iterator fit_end = sel.faces_end();
			vector3 min, max;
			min = boxes[i].first;
			min.z() = depth;
			max = boxes[i].second;
			max.z() = depth;
			Selector sel2(g);
			for (; fit != fit_end; ++fit) {
				if(BoxBoundProbe(CalculateCenter(*fit, aaPos), min, max)) {
					sel2.select(*fit);
				}
			}
			CloseSelection(sel2);
			AssignSelectionToSubset(sel2, sh, i);
			sel2.clear();
		}
		EraseEmptySubsets(sh);
		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_13.ugx");

		UG_LOGN("Tetrahedralize...")
		/// Tetrahedralize whole grid (TODO: If we don't pre-refine Tetgen breaks
		/// down and disrespects the boundaries somehow)
		Tetrahedralize(g, 5, false, true, aPosition, 1);

		EraseEmptySubsets(sh);
		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_14.ugx");

		/// Reassign the elements in the layers to subsets (uses ordering from above)
		for (size_t i = 0; i < boxes.size(); i++) {
			sel.clear();
			SelectSubsetElements<Volume>(sel, sh, siFaces, true);
			Selector::traits<Volume>::iterator fit = sel.volumes_begin();
			Selector::traits<Volume>::iterator fit_end = sel.volumes_end();
			vector3 min, max;
			min = boxes[i].first;
			min.z() = -depth;
			max = boxes[i].second;
			max.z() = depth;
			Selector sel2(g);
			for (; fit != fit_end; ++fit) {
				if(BoxBoundProbe(CalculateCenter(*fit, aaPos), min, max)) {
					sel2.select(*fit);
				}
			}
			CloseSelection(sel2);
			AssignSelectionToSubset(sel2, sh, i);
			sel2.clear();
		}

		/// Set subset names
		sh.subset_info(0).name = "FE1";
		sh.subset_info(1).name = "BD1";
		sh.subset_info(2).name = "FE2";
		sh.subset_info(3).name = "BD2";
		sh.subset_info(4).name = "MD";

		/// Reassign top boundary
		sel.clear();
		SelectSubsetElements<Face>(sel, sh, 0, true);
		Selector::traits<Face>::iterator fit = sel.faces_begin();
		Selector::traits<Face>::iterator fit_end = sel.faces_end();
		vector3 min, max;
		min = vector3(0, height, 0);
		max = vector3(width, height, depth);
		Selector sel2(g);
		for (; fit != fit_end; ++fit) {
			if(BoxBoundProbe(CalculateCenter(*fit, aaPos), min, max)) {
				sel2.select(*fit);
			}
		}
		CloseSelection(sel2);
		AssignSelectionToSubset(sel2, sh, sh.num_subsets());
		sel2.clear();

		/// Reassign bottom boundary
		sel.clear();
		SelectSubsetElements<Face>(sel, sh, 2, true);
		fit = sel.faces_begin();
		fit_end = sel.faces_end();
		min = vector3(0, -height-spacing, 0);
		max = vector3(width, -height-spacing, depth);
		for (; fit != fit_end; ++fit) {
			if(BoxBoundProbe(CalculateCenter(*fit, aaPos), min, max)) {
				sel2.select(*fit);
			}
		}

		CloseSelection(sel2);
		AssignSelectionToSubset(sel2, sh, sh.num_subsets());
		sel2.clear();

		EraseEmptySubsets(sh);
		sh.subset_info(5).name = "Top";
		sh.subset_info(6).name = "Bottom";
		SaveGridToFile(g, sh, "crack_generator_simple_step_15.ugx");
		sel.clear();

		/// Save final grid after optimization
		UG_LOGN("Writing final grid...")
		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_final.ugx");

		UG_COND_THROW(sh.num_subsets() != 7, "Number of subsets not seven (7). "
				"Something must have gone wrong. Use final grid with care!")
		}
	}
}
