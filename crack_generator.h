/*!
 * \file crack_generator.h
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
		 * \param[in] height
		 * \param[in] width
		 * \param[in] depth
		 * \param[in] thickness
		 * \param[in] spacing
		 * \param[in] preRefinement
		 * \param[in] postRefinement
		 */
		void BuildSimpleCrack
		(
			number height,
			number width,
			number depth,
			number thickness,
			number spacing,
			size_t preRefinement,
			size_t postRefinement
		);
	}
}

#endif /// UG__PLUGINS__CRACK_GENERATOR__CRACK_GENERATOR_H
