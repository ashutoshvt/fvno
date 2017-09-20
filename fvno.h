#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/wavefunction.h"
namespace psi { namespace fvno{
void gs_mp2_density_vv(SharedWavefunction ref_wfn, std::shared_ptr<PSIO> psio);
}}
