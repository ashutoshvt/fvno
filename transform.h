#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/wavefunction.h"
namespace psi { namespace fvno{
void transform_to_mo(SharedWavefunction ref_wfn, std::shared_ptr<PSIO> psio);
}}
