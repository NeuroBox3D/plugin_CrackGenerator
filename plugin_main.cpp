/*!
 * \file plugin_main.cpp
 * Registry for crack generator functions.
 */
#include "registry/registry.h"
#include "common/ug_config.h"
#include "common/error.h"
#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"
#include "crack_generator.h"
#include <string>

using namespace std;

extern "C" UG_API void
InitUGPlugin_CrackGenerator(ug::bridge::Registry* reg, string parentGroup)
{
  string grp(parentGroup); grp.append("CrackGenerator/");
  reg->add_function("BuildCompleteCrack", &ug::BuildCompleteCrack, "", "crackInnerLength#innerThickness#crackOuterLength#angle (degree)", grp);
  reg->add_function("BuildSimpleCrack", &ug::BuildSimpleCrack, "", "height#width#depth#thickness#refinement#spacing", grp);
}
