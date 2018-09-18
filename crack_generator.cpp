/*
 * dendrite_generator.cpp
 *
 */

#include "crack_generator.h"
#include "lib_grid/lib_grid.h"
#include "lib_grid/algorithms/remove_duplicates_util.h"
#include "bridge/domain_bridges/selection_bridge.h"
#include <algorithm>
#include "../ProMesh/mesh.h"
#include "../ProMesh/tools/grid_generation_tools.h"
#include "../ProMesh/tools/remeshing_tools.h"
#include "../ProMesh/tools/selection_tools.h"
#include "../ProMesh/tools/new_tools.h"
#include "../ProMesh/tools/subset_tools.h"
#include "../ProMesh/tools/refinement_tools.h"
#include "../ProMesh/tools/coordinate_transform_tools.h"
#include "../ProMesh/tools/topology_tools.h"

#define UG_ENABLE_WARNINGS

using namespace ug::promesh;

namespace ug {
	void BuildCrack() {
		Grid g;
	    SubsetHandler sh(g);
	    sh.set_default_subset_index(0);
	    g.attach_to_vertices(aPosition);

	    Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
	    Selector sel(g);

		number crackInnerLength = 0.2;
		number crackOuterLength = 2.0;
		number innerThickness = 0.1;
		number angle = 10;
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
		number outerHeight = 5 - outerDistance;
		std::cout << "outerHeight: " << outerHeight << std::endl;

		ug::vector3 topLeft, bottomLeft, topRight, bottomRight;
		topLeft = endPoint1;
		std::cout << "topLeft: " << topLeft << std::endl;
		topLeft.y() = topLeft.y() - outerHeight / 2;
		std::cout << "topLeft: " << topLeft << std::endl;

		bottomLeft = endPoint2;
		bottomLeft.y() = bottomLeft.y() + outerHeight / 2;

		Vertex* v7 = *g.create<RegularVertex>();
		Vertex* v8 = *g.create<RegularVertex>();
		aaPos[v7] = topLeft;
		aaPos[v8] = bottomLeft;

		Edge* e7 = *g.create<RegularEdge>(EdgeDescriptor(v5, v7));
		Edge* e8 = *g.create<RegularEdge>(EdgeDescriptor(v6, v8));

		topRight = topLeft;
		topRight.x() = topRight.x() + outerHeight+outerDistance;

		bottomRight = bottomLeft;
		bottomRight.x() = bottomRight.x() + outerHeight+outerDistance;

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
		number innerHeight = 0.25 - innerDistance;
		UG_COND_THROW(signbit(innerHeight), "inner height cannot be negative");

		topLeft = aaPos[v1];
		topLeft.y() = topLeft.y() - innerHeight / 2;
		bottomLeft = aaPos[v2];
		bottomLeft.y() = bottomLeft.y() + innerHeight / 2;

		Vertex* v11 = *g.create<RegularVertex>();
		Vertex* v12 = *g.create<RegularVertex>();
		aaPos[v11] = topLeft;
		aaPos[v12] = bottomLeft;

		Edge* e12 = *g.create<RegularEdge>(EdgeDescriptor(v1, v11));
		Edge* e13 = *g.create<RegularEdge>(EdgeDescriptor(v2, v12));

		topRight = topLeft;
		topRight.x() = topRight.x() + innerHeight+innerDistance;
		bottomRight = bottomLeft;
		bottomRight.x() = bottomRight.x() + innerHeight+innerDistance;

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
		number middleHeight = 0.50 - middleDistance;
		UG_COND_THROW(signbit(middleHeight), "middle height cannot be negative");

		topLeft = aaPos[v3];
		topLeft.y() = topLeft.y() - middleHeight / 2;
		bottomLeft = aaPos[v4];
		bottomLeft.y() = bottomLeft.y() + middleHeight / 2;

		Vertex* v15 = *g.create<RegularVertex>();
		Vertex* v16 = *g.create<RegularVertex>();
		aaPos[v15] = topLeft;
		aaPos[v16] = bottomLeft;

		Edge* e17 = *g.create<RegularEdge>(EdgeDescriptor(v3, v15));
		Edge* e18 = *g.create<RegularEdge>(EdgeDescriptor(v4, v16));

		topRight = topLeft;
		topRight.x() = topRight.x() + middleHeight+middleDistance;
		bottomRight = bottomLeft;
		bottomRight.x() = bottomRight.x() + middleHeight+middleDistance;

		Vertex* v17 = *g.create<RegularVertex>();
		Vertex* v18 = *g.create<RegularVertex>();
		aaPos[v17] = topRight;
		aaPos[v18] = bottomRight;

		Edge* e19 = *g.create<RegularEdge>(EdgeDescriptor(v15, v17));
		Edge* e20 = *g.create<RegularEdge>(EdgeDescriptor(v16, v18));
		Edge* e21 = *g.create<RegularEdge>(EdgeDescriptor(v17, v18));

		AssignSubsetColors(sh);
		SaveGridToFile(g, sh, "crack_generator_step_4.ugx");
  }
}
