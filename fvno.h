#ifndef FVNO_H
#define FVNO_H

#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/liboptions/liboptions.h"

namespace psi { namespace fvno{
class FVNO{
public:
    FVNO(SharedWavefunction ref_wfn, std::shared_ptr<PSIO> psio, Options& options);
    ~FVNO();
    int occ_;
    int vir_;
    int nmo_;
    int nso_;
    int nirrep_;
    int frz_vno_;
    SharedWavefunction ref_wfn_;
    std::shared_ptr<PSIO> psio_;
    SharedVector epsilon_;
    SharedMatrix C_;
    SharedMatrix F_so_;
    SharedMatrix F_mo_;
    SharedMatrix gs_density_;
    void gs_mp2_density_vv();
    void truncate_VNOs();
    void semicanonicalize_VNOs();
};
}}
#endif

