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
  reg->add_function("BuildCrack", &ug::BuildCrack, "", "crackInnerLength#innerThickness#crackOuterLength#angle (degree)", grp);
}
