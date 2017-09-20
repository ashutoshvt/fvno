#include "fvno.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/psifiles.h"
#include "psi4/libpsi4util/PsiOutStream.h" 
#define ID(x) ints.DPD_ID(x)

namespace psi { namespace fvno{

    FVNO::FVNO(SharedWavefunction ref_wfn, std::shared_ptr<PSIO> psio, Options& options){

        psio_ = psio;
        ref_wfn_ = ref_wfn;
        nirrep_  = ref_wfn_->nirrep();
        nso_ = ref_wfn_->nso();
        nmo_ = ref_wfn_->nmo();
        occ_ = ref_wfn_->doccpi()[0]; // hard-coded for symmetry C1
        vir_ = nmo_-occ_;
        frz_vno_ = options.get_int("FRZ_VNO");
        outfile->Printf("\n\t Number of frozen VNOs: %d\n",frz_vno_);

        epsilon_ = SharedVector(new Vector(" orbital energies",nmo_));
        epsilon_->copy(ref_wfn_->epsilon_a()->clone());
        C_ = SharedMatrix(new Matrix("Canonical MO Coefficients", nso_, nmo_));
        C_->copy(ref_wfn_->Ca()->clone());
        F_so_ = SharedMatrix(new Matrix("Fock Matrix SO basis", nso_, nso_));
        F_so_->copy(ref_wfn_->Fa()->clone());
        F_mo_ = SharedMatrix(new Matrix("Fock Matrix Canonical MO basis", nmo_, nmo_));
        F_mo_->transform(F_so_, C_);
        gs_density_ = SharedMatrix(new Matrix("ground state mp2 density MO basis (vir-vir)", vir_, vir_));
    }

    FVNO::~FVNO()
    {}

    void FVNO::gs_mp2_density_vv(){

        dpdbuf4 D,D1;
        dpdfile2 Density;

        std::vector<std::shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::occ);
        spaces.push_back(MOSpace::vir);
        IntegralTransform ints(ref_wfn_, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly);
        ints.set_keep_iwl_so_ints(true);
        ints.transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
        /* Use the IntegralTransform object's DPD instance, for convenience.
           This is also the reason that I needed to combine the transformation
           code with the construction of density code inside a single function 
           as when ints object gets out of scope, the dpd object associated with ints 
           goes out of scope as well.
        */
        dpd_set_default(ints.get_dpd_id());

        dpdbuf4 K;
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "D <ij|ab>");
        global_dpd_->buf4_close(&K);

        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
        global_dpd_->buf4_copy(&D, PSIF_LIBTRANS_DPD , "D' <ij|ab>");
        global_dpd_->buf4_close(&D);

        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "D' <ij|ab>");
        global_dpd_->buf4_mat_irrep_init(&D, 0); /* here and below: hard-coded for symmetry C1 (irrep 0) */
        global_dpd_->buf4_mat_irrep_rd(&D, 0);
            for(int i=0,ij=0;i<occ_;i++)
            for(int j=0;j<occ_;j++,ij++){
            for(int a=0,ab=0;a<vir_;a++)
            for(int b=0;b<vir_;b++,ab++)
           D.matrix[0][ij][ab] /= (epsilon_->get(i) + epsilon_->get(j) - F_mo_->get(a+occ_,a+occ_) - F_mo_->get(b+occ_,b+occ_));
         }
        global_dpd_->buf4_mat_irrep_wrt(&D, 0);
        global_dpd_->buf4_mat_irrep_close(&D, 0);
        global_dpd_->buf4_close(&D);


        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
        global_dpd_->buf4_scmcopy(&D, PSIF_LIBTRANS_DPD, "D 2<ij|ab> - <ij|ba>", 2);
        global_dpd_->buf4_sort_axpy(&D, PSIF_LIBTRANS_DPD, pqsr, 0, 5, "D 2<ij|ab> - <ij|ba>", -1);
        global_dpd_->buf4_close(&D);
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
        global_dpd_->buf4_mat_irrep_init(&D, 0);
        global_dpd_->buf4_mat_irrep_rd(&D, 0);
           for(int i=0,ij=0;i<occ_;i++)
           for(int j=0;j<occ_;j++,ij++){
           for(int a=0,ab=0;a<vir_;a++)
           for(int b=0;b<vir_;b++,ab++)
              D.matrix[0][ij][ab] /= (epsilon_->get(i) + epsilon_->get(j) - F_mo_->get(a+occ_,a+occ_) - F_mo_->get(b+occ_,b+occ_));
          }
        global_dpd_->buf4_mat_irrep_wrt(&D, 0);
        global_dpd_->buf4_mat_irrep_close(&D, 0);
        global_dpd_->buf4_close(&D);

        global_dpd_->file2_init(&Density,PSIF_LIBTRANS_DPD,0,1,1,"d(a,b)");
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "D' <ij|ab>");
        global_dpd_->buf4_init(&D1, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
        global_dpd_->contract442(&D1,&D,&Density,2,2,2.0,0);
        global_dpd_->buf4_close(&D);
        global_dpd_->buf4_close(&D1);
        global_dpd_->file2_close(&Density);

        global_dpd_->file2_init(&Density,PSIF_LIBTRANS_DPD,0,1,1,"d(a,b)");
        global_dpd_->file2_mat_init(&Density);
        global_dpd_->file2_mat_rd(&Density);
         for(int a=0;a<vir_;a++)
            for(int b=0;b<vir_;b++)
              gs_density_->set(a,b,Density.matrix[0][a][b]);
        global_dpd_->file2_mat_close(&Density);
        global_dpd_->file2_close(&Density);

    }

    void FVNO::truncate_VNOs(){
    
    /* CVMO : Canonical virtual MOs
       VNO: Virtual Natural orbitals
       ONs: Occupation numbers 
    */

    SharedMatrix eigenvectors(new Matrix("Eigen-vectors (CVMO->VNO transformation)", vir_, vir_));
    SharedVector eigenvalues(new Vector("Eigen-values (ONs)", vir_));
    SharedVector zero_vec(new Vector("zero vector", vir_));
    gs_density_->diagonalize(eigenvectors, eigenvalues, descending);
    eigenvectors->print();
    eigenvalues->print();
    for(int a=0;a<frz_vno_;a++)
      eigenvectors->set_column(0, vir_-a-1,zero_vec);
    eigenvectors->print();   

    }

    void FVNO::semicanonicalize_VNOs(){
        int a;
    }

}}
