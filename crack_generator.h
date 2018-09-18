/*!
 * \file crack_generator.h
 */
#ifndef UG__PLUGINS__CRACK_GENERATOR__CRACK_GENERATOR_H
#define  UG__PLUGINS__CRACK_GENERATOR__CRACK_GENERATOR_H

#include <common/types.h>

namespace ug {
	/*!
	 * \brief builds a crack geometry
	 * \param[in] crackInnerLength
	 * \param[in] innerThickness
	 * \param[in] crackOuterLength
	 * \param[in] angle
	 */
	void BuildCrack(number crackInnerLength, number innerThickness, number crackOuterLength, number angle);
}

#endif /// UG__PLUGINS__CRACK_GENERATOR__CRACK_GENERATOR_H
