#ifndef FVNO_H
#define FVNO_H

#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libtrans/integraltransform.h"

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
    int on_cutoff_;
    SharedWavefunction ref_wfn_;
    std::shared_ptr<PSIO> psio_;
    IntegralTransform * ints_;
    SharedVector epsilon_;
    SharedVector ONs_;
    SharedMatrix C_;
    SharedMatrix F_so_;
    SharedMatrix F_mo_;
    SharedMatrix VNOs_;
    SharedMatrix gs_density_;
    void transform_mo_mp2();
    void gs_mp2_density_vv();
    void truncate_VNOs();
    void semicanonicalize_VNOs();

};
}}
#endif

