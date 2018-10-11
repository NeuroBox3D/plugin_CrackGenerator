/*!
 * \file crack_generator.h
 * Author: Stephan Grein
 */
#ifndef UG__PLUGINS__CRACK_GENERATOR__CRACK_GENERATOR_H
#define  UG__PLUGINS__CRACK_GENERATOR__CRACK_GENERATOR_H

#include <common/types.h>

namespace ug {
	namespace crack_generator {
		/*!
		 * \brief builds a complex crack geometry
		 * \param[in] crackInnerLength
		 * \param[in] innerThickness
		 * \param[in] crackOuterLength
		 * \param[in] angle
		 */
		void BuildCompleteCrack
		(
			number crackInnerLength,
			number innerThickness,
			number crackOuterLength,
			number angle
		);

		/*!
		 * \brief builds a simple crack geometry
		 * \param[in] height of cuboid
		 * \param[in] width of cuboid
		 * \param[in] depth of cuboid
		 * \param[in] thickness of one briding domain
		 * \param[in] spacing the size of the MD domain
		 * \param[in] h finess of grid
		 * \param[in] r_0 lattice constant
		 */
		void BuildSimpleCrack
		(
			number height,
			number width,
			number depth,
			number thickness,
			number spacing,
			number h,
			number r_0
		);
	}
}

#endif /// UG__PLUGINS__CRACK_GENERATOR__CRACK_GENERATOR_H
