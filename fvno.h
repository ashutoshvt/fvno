#ifndef FVNO_H
#define FVNO_H

#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
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
    int my_frozen_uocc_;
    double on_cutoff_;
    bool print_;
    bool vno_basis_;
    bool semicanonical_basis_;
    bool canonical_basis_;
    bool prop_correction_;
    int pert_density_order_;
    Options options_;
    SharedWavefunction ref_wfn_;
    std::shared_ptr<PSIO> psio_;
    IntegralTransform * ints_;
    std::shared_ptr<MintsHelper> mints_;
    std::vector<SharedMatrix> dipole_;
    std::vector<SharedMatrix> angmom_;
    std::vector<std::string> cart_;
    SharedVector epsilon_;
    SharedMatrix C_;
    SharedMatrix F_so_;
    SharedMatrix F_mo_;
    SharedMatrix gs_density_;
    SharedMatrix pert_density_singles_;
    SharedMatrix pert_density_doubles_;
    SharedMatrix combined_density_;
    SharedMatrix VNOs_;
    SharedVector ONs_;
    SharedMatrix SCan_VMO_;
    void transform_mo_mp2();
    void preppert(Options& options);
    //void sort_pert(std::string pert, std::string cart, SharedMatrix pert_mat);
    void gs_mp2_density_vv();
    void pert_density_vv(Options & options);
    void final_density(Options& options);
    void truncate_VNOs(Options& options);
    void final_basis();

};
}}
#endif

