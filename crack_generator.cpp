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
		/// BUILDCOMPLETECRACK
		////////////////////////////////////////////////////////////////////////////////
		/// TODO: Pre-refine the inner MD square
		/// TODO: Improve triangulation and tetrahedralization
		void BuildCompleteCrack
		(
			number crackInnerLength=0.2,
			number innerThickness=0.1,
			number crackOuterLength=2.0,
			number angle = 10
		) {
		Grid g;
	    SubsetHandler sh(g);
	    sh.set_default_subset_index(0);
	    g.attach_to_vertices(aPosition);
	    AInt aInt;
	    g.attach_to_vertices(aInt);
	    Grid::VertexAttachmentAccessor<AInt> aaIntVertex(g, aInt);

	    Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
	    Selector sel(g);

		Vertex* startVertex = *g.create<RegularVertex>();
		aaPos[startVertex] = ug::vector3(0, 0, 0); /// crack tip

		ug::vector3 endPoint1, endPoint2;
		endPoint1.x() = -innerThickness * cos(deg_to_rad(angle));
		endPoint1.y() = -innerThickness * sin(deg_to_rad(angle));
		endPoint1.z() = 0;

		endPoint2.x() = -innerThickness * cos(-deg_to_rad(angle));
		endPoint2.y() = -innerThickness * sin(-deg_to_rad(angle));
		endPoint2.z() = 0;

		Vertex* v1 = *g.create<RegularVertex>();
		Vertex* v2 = *g.create<RegularVertex>();
		aaPos[v1] = endPoint1;
		aaPos[v2] = endPoint2;

		Edge* e1 = *g.create<RegularEdge>(EdgeDescriptor(startVertex, v1));
		Edge* e2 = *g.create<RegularEdge>(EdgeDescriptor(startVertex, v2));

	    sh.set_default_subset_index(1);

		endPoint1.x() = -crackInnerLength * cos(deg_to_rad(angle));
		endPoint1.y() = -crackInnerLength * sin(deg_to_rad(angle));
		endPoint1.z() = 0;

		endPoint2.x() = -crackInnerLength * cos(-deg_to_rad(angle));
		endPoint2.y() = -crackInnerLength * sin(-deg_to_rad(angle));
		endPoint2.z() = 0;

		Vertex* v3 = *g.create<RegularVertex>();
		Vertex* v4 = *g.create<RegularVertex>();
		aaPos[v3] = endPoint1;
		aaPos[v4] = endPoint2;

		Edge* e3 = *g.create<RegularEdge>(EdgeDescriptor(v1, v3));
		Edge* e4 = *g.create<RegularEdge>(EdgeDescriptor(v2, v4));

	    sh.set_default_subset_index(2);

		endPoint1.x() = -crackOuterLength * cos(deg_to_rad(angle));
		endPoint1.y() = -crackOuterLength * sin(deg_to_rad(angle));
		endPoint1.z() = 0;

		endPoint2.x() = -crackOuterLength * cos(-deg_to_rad(angle));
		endPoint2.y() = -crackOuterLength * sin(-deg_to_rad(angle));
		endPoint2.z() = 0;

		Vertex* v5 = *g.create<RegularVertex>();
		Vertex* v6 = *g.create<RegularVertex>();
		aaPos[v5] = endPoint1;
		aaPos[v6] = endPoint2;

		Edge* e5 = *g.create<RegularEdge>(EdgeDescriptor(v3, v5));
		Edge* e6 = *g.create<RegularEdge>(EdgeDescriptor(v4, v6));

		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_step_1.ugx");
	    sh.set_default_subset_index(3);

	    /// outermost square
		number outerDistance = VecDistance(endPoint1, endPoint2);

		ug::vector3 centerOuter = endPoint1;
		centerOuter.y() = centerOuter.y() + outerDistance / 2;
		number squareOuterDiameter = VecDistance(centerOuter, aaPos[startVertex]);
		number outerHeight = squareOuterDiameter;
		outerHeight = 2*squareOuterDiameter - outerDistance;

		ug::vector3 topLeft, bottomLeft, topRight, bottomRight;
		topLeft = endPoint1;
		topLeft.y() = topLeft.y() - outerHeight/2;

		bottomLeft = endPoint2;
		bottomLeft.y() = bottomLeft.y() + outerHeight/2;

		Vertex* v7 = *g.create<RegularVertex>();
		Vertex* v8 = *g.create<RegularVertex>();
		aaPos[v7] = topLeft;
		aaPos[v8] = bottomLeft;

		Edge* e7 = *g.create<RegularEdge>(EdgeDescriptor(v5, v7));
		Edge* e8 = *g.create<RegularEdge>(EdgeDescriptor(v6, v8));

		topRight = topLeft;
		topRight.x() = topRight.x() + 2.0*squareOuterDiameter;

		bottomRight = bottomLeft;
		bottomRight.x() = bottomRight.x() + 2.0*squareOuterDiameter;

	    sh.set_default_subset_index(7);
		Vertex* v9 = *g.create<RegularVertex>();
		Vertex* v10 = *g.create<RegularVertex>();
		aaPos[v9] = topRight;
		aaPos[v10] = bottomRight;

	    sh.set_default_subset_index(5);
		Edge* e9 = *g.create<RegularEdge>(EdgeDescriptor(v7, v9));
	    sh.set_default_subset_index(6);
		Edge* e10 = *g.create<RegularEdge>(EdgeDescriptor(v8, v10));
	    sh.set_default_subset_index(7);
		Edge* e11 = *g.create<RegularEdge>(EdgeDescriptor(v9, v10));

		AssignSubsetColors(sh);
	    sh.set_default_subset_index(0);
		SaveGridToFile(g, sh, "crack_generator_step_2.ugx");

		/// innermost square
		number innerDistance = VecDistance(aaPos[v1], aaPos[v2]);

		ug::vector3 centerInner = aaPos[v1];
		centerInner.y() = centerInner.y() + innerDistance / 2;
		number squareInnerDiameter = VecDistance(centerInner, aaPos[startVertex]);
		number innerHeight = squareInnerDiameter;
		innerHeight = 2*squareInnerDiameter - innerDistance;

		UG_COND_THROW(std::signbit(innerHeight), "inner height cannot be negative");

		topLeft = aaPos[v1];
		topLeft.y() = topLeft.y() - innerHeight/2;
		bottomLeft = aaPos[v2];
		bottomLeft.y() = bottomLeft.y() + innerHeight/2;

		Vertex* v11 = *g.create<RegularVertex>();
		Vertex* v12 = *g.create<RegularVertex>();
		aaPos[v11] = topLeft;
		aaPos[v12] = bottomLeft;

		Edge* e12 = *g.create<RegularEdge>(EdgeDescriptor(v1, v11));
		Edge* e13 = *g.create<RegularEdge>(EdgeDescriptor(v2, v12));

		topRight = topLeft;
		topRight.x() = topRight.x() + 2*squareInnerDiameter;
		bottomRight = bottomLeft;
		bottomRight.x() = bottomRight.x() + 2*squareInnerDiameter;

		Vertex* v13 = *g.create<RegularVertex>();
		Vertex* v14 = *g.create<RegularVertex>();
		aaPos[v13] = topRight;
		aaPos[v14] = bottomRight;

		Edge* e14 = *g.create<RegularEdge>(EdgeDescriptor(v11, v13));
		Edge* e15 = *g.create<RegularEdge>(EdgeDescriptor(v12, v14));
		Edge* e16 = *g.create<RegularEdge>(EdgeDescriptor(v13, v14));

		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_step_3.ugx");

		/// middle square
	    sh.set_default_subset_index(1);
		number middleDistance = VecDistance(aaPos[v3], aaPos[v4]);

		ug::vector3 centerMiddle = aaPos[v3];
		centerMiddle.y() = centerMiddle.y() + middleDistance / 2;
		number squareMiddleDiameter = VecDistance(centerMiddle, aaPos[startVertex]);
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

		Edge* e17 = *g.create<RegularEdge>(EdgeDescriptor(v3, v15));
		Edge* e18 = *g.create<RegularEdge>(EdgeDescriptor(v4, v16));

		topRight = topLeft;
		topRight.x() = topRight.x() + 2.0*squareMiddleDiameter;
		bottomRight = bottomLeft;
		bottomRight.x() = bottomRight.x() + 2.0*squareMiddleDiameter;

		Vertex* v17 = *g.create<RegularVertex>();
		Vertex* v18 = *g.create<RegularVertex>();
		aaPos[v17] = topRight;
		aaPos[v18] = bottomRight;

		Edge* e19 = *g.create<RegularEdge>(EdgeDescriptor(v15, v17));
		Edge* e20 = *g.create<RegularEdge>(EdgeDescriptor(v16, v18));
		Edge* e21 = *g.create<RegularEdge>(EdgeDescriptor(v17, v18));

		EraseEmptySubsets(sh); /// necessary since previously assigned outerSquare is empty because this compromises the boundary

		/// Rename subsets
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
		for (size_t i = 0; i < sh.num_subsets(); i++) {
			SelectSubsetElements<ug::Edge>(sel, sh, i, true);
		}
		TriangleFill_SweepLine(g, sel.edges_begin(), sel.edges_end(), aPosition, aInt, &sh, 7);
		SelectSubsetElements<ug::Face>(sel, sh, 7, true);
		QualityGridGeneration(g, sel.faces_begin(), sel.faces_end(), aaPos, 10);
		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_step_5.ugx");

		/// Extrude towards top
		ug::vector3 normal = ug::vector3(0, 0, 2*squareOuterDiameter);
		std::vector<Edge*> edges;
		for (size_t i = 0; i < sh.num_subsets(); i++) {
				SelectSubsetElements<ug::Edge>(sel, sh, i, true);
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
		ug::vector3 bottomLeft,
		ug::vector3 bottomRight,
		ug::vector3 leftMDLayer,
		ug::vector3 rightMDLayer,
		ug::vector3 topLeft,
		ug::vector3 topRight,
		Grid& g,
		SubsetHandler& sh,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		AInt& aInt,
		size_t refinements,
		Selector& sel,
		number depth,
		size_t si_offset,
		std::vector<ug::Vertex*>& verts
	) {
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

		sh.set_default_subset_index(si_offset);
		Edge* e7 = *g.create<RegularEdge>(EdgeDescriptor(leftMDLayerVertex, rightMDLayerVertex));

		/// This is pre-refinement of the Bridging Domains...
		/// TODO: refine as long as lattice constant not reached! -> lattice constant given and length of subset => #refs = length of subset / lattice constant = number of refinements
		/// TODO Then go from bottom left corner up to top left corner -> introduce vertices on both sides in regular spacing subset_length / #refs -> create edge -> refine edge until lattice constant met.
		sel.clear();
		sel.select(e1);
		sel.select(e7);
		for (size_t i = 0; i < refinements; i++) {
			Refine(g, sel);
		}
		sel.clear();
		sh.set_default_subset_index(si_offset);

		AssignSubsetColors(sh);
		step << "crack_generator_simple_step_" << si_offset+4 << ".ugx";
		SaveGridToFile(g, sh, step.str().c_str());
		step.str(""); step.clear();

		verts.push_back(bottomLeftVertex);
		verts.push_back(bottomRightVertex);
		}

		////////////////////////////////////////////////////////////////////////////////
		/// BUILDSIMPLECRACK
		////////////////////////////////////////////////////////////////////////////////
		void BuildSimpleCrack
		(
			number height,
			number width,
			number depth,
			number thickness,
			number spacing,
			size_t preRefinements,
			size_t postRefinements
		) {
		/// TODO: refine not for number of refinemnts, but input lattice constant r_0, and refine until met
		///       --> refine all subsets to have uniform mesh:
		/// Algorithm:
		/// 1. Create line from bottomLeft to bottomRight
		/// 2. Create line from bottomLeft to leftMDLayer and leftMDLayer to topLeft
		/// 3. Create line from bottomRight to rightMDLayer and rightMDLayer to topRight
		/// 4. Create line from rightMDLayer to leftMDLayer and topRight to topLeft
		UG_COND_THROW(thickness==height, "Thickness can't be the same as height.");

		Grid g;
	    SubsetHandler sh(g);
	    sh.set_default_subset_index(0);
	    g.attach_to_vertices(aPosition);
	    AInt aInt;
	    g.attach_to_vertices(aInt);
	    Grid::VertexAttachmentAccessor<AInt> aaIntVertex(g, aInt);

	    Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
	    Selector sel(g);

	    //// First rectangle
   		ug::vector3 bottomLeft = ug::vector3(0, 0, 0);
		ug::vector3 bottomRight = ug::vector3(width, 0, 0);
		ug::vector3 leftMDLayer = ug::vector3(0, thickness, 0);
		ug::vector3 rightMDLayer = ug::vector3(width, thickness, 0);
		ug::vector3 topLeft = ug::vector3(0, height, 0);
		ug::vector3 topRight = ug::vector3(width, height, 0);
		std::vector<std::pair<ug::vector3, ug::vector3> > boxes;
		boxes.push_back(std::make_pair(leftMDLayer, topRight));
		boxes.push_back(std::make_pair(bottomLeft, rightMDLayer));
		size_t si_offset = 0;
		std::vector<ug::Vertex*> verts;
		create_rect(bottomLeft, bottomRight, leftMDLayer, rightMDLayer, topLeft, topRight, g, sh, aaPos, aInt, preRefinements, sel, depth, si_offset, verts);
		si_offset = 3;

		/// Second rectangle
   		bottomLeft = ug::vector3(0, -spacing, 0);
		bottomRight = ug::vector3(width, -spacing, 0);
		leftMDLayer = ug::vector3(0, -spacing-thickness, 0);
		rightMDLayer = ug::vector3(width, -spacing-thickness, 0);
		topLeft = ug::vector3(0, -spacing-height, 0);
		topRight = ug::vector3(width, -spacing-height, 0);
		boxes.push_back(std::make_pair(topLeft, rightMDLayer));
		boxes.push_back(std::make_pair(leftMDLayer, bottomRight));
		create_rect(bottomLeft, bottomRight, leftMDLayer, rightMDLayer, topLeft, topRight, g, sh, aaPos, aInt, preRefinements, sel, depth, si_offset, verts);

		/// Connect the two rectangles
		sh.set_default_subset_index(2*si_offset);
		ug::Edge* e1 = *g.create<RegularEdge>(EdgeDescriptor(verts[0], verts[2]));
		ug::Edge* e2 = *g.create<RegularEdge>(EdgeDescriptor(verts[1], verts[3]));
		boxes.push_back(std::make_pair(ug::vector3(0, -spacing, 0), ug::vector3(width, 0, 0)));

		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_9.ugx");

		/// Triangulate bottom
		for (size_t i = 0; i < sh.num_subsets(); i++) {
			SelectSubsetElements<ug::Edge>(sel, sh, i, true);
		}
		TriangleFill_SweepLine(g, sel.edges_begin(), sel.edges_end(), aPosition, aInt, &sh, sh.num_subsets());

		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_10.ugx");

		/// Reassign the elements in the layers to subsets -> start beyond the current subsets, thus this ordering is the same as below
		size_t siFaces = sh.num_subsets()-1;
		for (size_t i = 0; i < boxes.size(); i++) {
			sel.clear();
			SelectSubsetElements<Face>(sel, sh, siFaces, true);
			Selector::traits<Face>::iterator fit = sel.faces_begin();
			Selector::traits<Face>::iterator fit_end = sel.faces_end();
			ug::vector3 min, max;
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
		/// Retriangulate all
		for (size_t i = 0; i < sh.num_subsets(); i++) {
			SelectSubsetElements<ug::Face>(sel, sh, i, true);
		}
		QualityGridGeneration(g, sel.faces_begin(), sel.faces_end(), aaPos, 30);
		SaveGridToFile(g, sh, "crack_generator_simple_step_11.ugx");

		/// Extrude all in steps to ensure uniformity
		ug::vector3 normal = ug::vector3(0, 0, depth);
		std::vector<Edge*> edges;
		for (size_t i = 0; i < sh.num_subsets(); i++) {
			SelectSubsetElements<ug::Edge>(sel, sh, i, true);
		}
		edges.assign(sel.edges_begin(), sel.edges_end());
		VecScale(normal, normal, 1.0/((number)(2*(preRefinements+1))));
		number totalLength = normal.z();
		while (totalLength < depth) {
			std::cout << "Extrude! length: " << totalLength << std::endl;
			Extrude(g, NULL, &edges, NULL, normal, aaPos, EO_CREATE_FACES, NULL);
			totalLength += normal.z();
		}
		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_12.ugx");

		/// Triangulate top
		TriangleFill_SweepLine(g, edges.begin(), edges.end(), aPosition, aInt, &sh, sh.num_subsets());
		AssignSubsetColors(sh);
		EraseEmptySubsets(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_13.ugx");

		siFaces = sh.num_subsets()-1;
		/// Reassign the elements in the layers to subsets (uses ordering from above)
		for (size_t i = 0; i < boxes.size(); i++) {
			sel.clear();
			SelectSubsetElements<Face>(sel, sh, siFaces, true);
			Selector::traits<Face>::iterator fit = sel.faces_begin();
			Selector::traits<Face>::iterator fit_end = sel.faces_end();
			ug::vector3 min, max;
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
		SaveGridToFile(g, sh, "crack_generator_simple_step_14.ugx");

		/// Tetrahedralize whole grid (Note: If we don't pre-refine Tetgen brakes down and disrespects the boundaries somehow)
		Tetrahedralize(g, 5, false, true, aPosition, 1);

		EraseEmptySubsets(sh);
		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_15.ugx");

		/// Reassign the elements in the layers to subsets (uses ordering from above)
		for (size_t i = 0; i < boxes.size(); i++) {
			sel.clear();
			SelectSubsetElements<Volume>(sel, sh, siFaces, true);
			Selector::traits<Volume>::iterator fit = sel.volumes_begin();
			Selector::traits<Volume>::iterator fit_end = sel.volumes_end();
			ug::vector3 min, max;
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
		ug::vector3 min, max;
		min = ug::vector3(0, height, 0);
		max = ug::vector3(width, height, depth);
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
		min = ug::vector3(0, -height-spacing, 0);
		max = ug::vector3(width, -height-spacing, depth);
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
		SaveGridToFile(g, sh, "crack_generator_simple_step_16.ugx");

		/// Refine BD subsets (This is post-refinement) TODO probably not necessary!
		std::vector<std::string> bdDomains;
		bdDomains.push_back("BD1");
		bdDomains.push_back("BD2");
		sel.clear();
		for (std::vector<std::string>::const_iterator it = bdDomains.begin(); it != bdDomains.end(); ++it) {
			SelectSubsetElements<ug::Face>(sel, sh, sh.get_subset_index(it->c_str()), true);
			SelectSubsetElements<ug::Volume>(sel, sh, sh.get_subset_index(it->c_str()), true);
			SelectSubsetElements<ug::Edge>(sel, sh, sh.get_subset_index(it->c_str()), true);
			SelectSubsetElements<ug::Vertex>(sel, sh, sh.get_subset_index(it->c_str()), true);
			CloseSelection(sel);
			for (size_t i = 0; i < postRefinements; i++) {
				Refine(g, sel);
			}
			sel.clear();
		}
		sel.clear();

		/// Save final grid after optimization
		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_simple_step_final.ugx");

		UG_COND_THROW(sh.num_subsets() != 7, "Number of subsets not seven (7). "
				"Something must have gone wrong. Use final grid with care!")
		}
	}
}
