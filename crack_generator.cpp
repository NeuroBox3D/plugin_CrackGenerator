/*!
 * \file crack_generator.cpp
 *
 */

#include "crack_generator.h"
#include "lib_grid/lib_grid.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"
#include <cmath>

#define UG_ENABLE_WARNINGS

namespace ug {
	/// Note/TODO: We could pre-refine the inner squares!
	/// TODO: Improve triangulation and tetrahedralization!
	void BuildCrack
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
		aaPos[startVertex] = ug::vector3(0, 0, 0);

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
	    sh.set_default_subset_index(0);

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

		Vertex* v9 = *g.create<RegularVertex>();
		Vertex* v10 = *g.create<RegularVertex>();
		aaPos[v9] = topRight;
		aaPos[v10] = bottomRight;

		Edge* e9 = *g.create<RegularEdge>(EdgeDescriptor(v7, v9));
		Edge* e10 = *g.create<RegularEdge>(EdgeDescriptor(v8, v10));
		Edge* e11 = *g.create<RegularEdge>(EdgeDescriptor(v9, v10));

		AssignSubsetColors(sh);
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

		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_step_4.ugx");

		/// Triangulate bottom surface
		SelectSubsetElements<ug::Edge>(sel, sh, 0, true);
		SelectSubsetElements<ug::Edge>(sel, sh, 1, true);
		SelectSubsetElements<ug::Edge>(sel, sh, 2, true);
		TriangleFill_SweepLine(g, sel.edges_begin(), sel.edges_end(), aPosition, aInt, &sh, 3);
		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_step_5.ugx");

		/// Extrude towards top
		ug::vector3 normal = ug::vector3(0, 0, 2*squareOuterDiameter);
		std::vector<Edge*> edges;
		SelectSubsetElements<ug::Edge>(sel, sh, 0, true);
		SelectSubsetElements<ug::Edge>(sel, sh, 1, true);
		SelectSubsetElements<ug::Edge>(sel, sh, 2, true);
		SelectSubsetElements<ug::Edge>(sel, sh, 3, true);
		edges.assign(sel.edges_begin(), sel.edges_end());
		Extrude(g, NULL, &edges, NULL, normal, aaPos, EO_CREATE_FACES, NULL);

		/// Triangulate top surface
		TriangleFill_SweepLine(g, edges.begin(), edges.end(), aPosition, aInt, &sh, 4);

		AssignSubsetColors(sh);
		sh.subset_info(0).name = "Inner square";
		sh.subset_info(1).name = "Middle square";
		sh.subset_info(2).name = "Outer square";
		sh.subset_info(3).name = "Bottom surface";
		sh.subset_info(4).name = "Top surface";
		SaveGridToFile(g, sh, "crack_generator_step_6.ugx");

		/// Tetrahedralize whole grid
		Tetrahedralize(g, 5, false, false, aPosition, 1);
		sh.subset_info(5).name = "Volumes";
		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_step_7.ugx");
  }
}
