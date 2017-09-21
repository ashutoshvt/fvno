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
#define ID(x) ints_->DPD_ID(x)

namespace psi { namespace fvno{

    FVNO::FVNO(SharedWavefunction ref_wfn, std::shared_ptr<PSIO> psio, Options& options){

        psio_ = psio;
        ref_wfn_ = ref_wfn;
        options_ = options;
        nirrep_  = ref_wfn_->nirrep();
        nso_ = ref_wfn_->nso();
        nmo_ = ref_wfn_->nmo();
        occ_ = ref_wfn_->doccpi()[0]; // hard-coded for symmetry C1
        vir_ = nmo_-occ_;
        frz_vno_ = options.get_int("FRZ_VNO");
        epsilon_ = SharedVector(new Vector(" orbital energies",nmo_));
        epsilon_->copy(ref_wfn_->epsilon_a()->clone());
        C_ = SharedMatrix(new Matrix("Canonical MO Coefficients", nso_, nmo_));
        C_->copy(ref_wfn_->Ca()->clone());
        F_so_ = SharedMatrix(new Matrix("Fock Matrix SO basis", nso_, nso_));
        F_so_->copy(ref_wfn_->Fa()->clone());
        F_mo_ = SharedMatrix(new Matrix("Fock Matrix Canonical MO basis", nmo_, nmo_));
        F_mo_->transform(F_so_, C_);
        gs_density_ = SharedMatrix(new Matrix("ground state mp2 density MO basis (vir-vir)", vir_, vir_));
        VNOs_ = SharedMatrix(new Matrix("Eigen-vectors (CVMO->VNO transformation)", vir_, vir_));
        ONs_ = SharedVector(new Vector("Eigen-values (Occupation numbers)", vir_));
        SCan_VMO_ = SharedMatrix(new Matrix(" Semicanonical basis (CVMO->VNO->SCMO transformation)", nmo_, vir_));
    }

    FVNO::~FVNO()
    {}

    void FVNO::transform_mo_mp2(){

        std::vector<std::shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::occ);
        spaces.push_back(MOSpace::vir);
        ints_ = new IntegralTransform(ref_wfn_, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly);
        ints_->set_keep_iwl_so_ints(true);
        ints_->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);

        /* Use the IntegralTransform object's DPD instance, for convenience.
           This is also the reason that I thought that I needed to combine the transformation
           code with the construction of density code inside a single function 
           as when ints object gets out of scope, the dpd object associated with ints 
           goes out of scope as well. Well, I can always make ints a member of the class
           to avoid this problem.
        */
        dpd_set_default(ints_->get_dpd_id());

        dpdbuf4 K, Dijab;
        dpdfile2 Dia;
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "D <ij|ab>");
        global_dpd_->buf4_close(&K);

        global_dpd_->buf4_init(&Dijab, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "Dijab");
        global_dpd_->buf4_mat_irrep_init(&Dijab, 0); /* here and below: hard-coded for symmetry C1 (irrep 0) */
        for(int i=0,ij=0;i<occ_;i++)
          for(int j=0;j<occ_;j++,ij++){
            for(int a=0,ab=0;a<vir_;a++)
              for(int b=0;b<vir_;b++,ab++)
                Dijab.matrix[0][ij][ab] = 1.0/(epsilon_->get(i) + epsilon_->get(j) - F_mo_->get(a+occ_,a+occ_) - F_mo_->get(b+occ_,b+occ_));
         }
        global_dpd_->buf4_mat_irrep_wrt(&Dijab, 0);
        global_dpd_->buf4_mat_irrep_close(&Dijab, 0);
        global_dpd_->buf4_close(&Dijab);


        global_dpd_->file2_init(&Dia, PSIF_LIBTRANS_DPD, 0, 0, 1, "Dia");
        global_dpd_->file2_mat_init(&Dia);
        global_dpd_->file2_mat_rd(&Dia);
         for(int i=0; i<occ_; i++)
           for(int a=0; a<vir_; a++)
             Dia.matrix[0][i][a] = 1.0/(epsilon_->get(i) - F_mo_->get(a+occ_,a+occ_));
        global_dpd_->file2_mat_close(&Dia);
        global_dpd_->file2_close(&Dia);
        

        }

        
    void FVNO::gs_mp2_density_vv(){

        dpdbuf4 D, D1, Dia, Dijab;
        dpdfile2 Density;
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
        global_dpd_->buf4_copy(&D, PSIF_LIBTRANS_DPD , "tIjAb(MP2)");
        global_dpd_->buf4_close(&D);

        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "tIjAb(MP2)");
        global_dpd_->buf4_init(&Dijab, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "Dijab");
        global_dpd_->buf4_dirprd(&Dijab, &D);
        global_dpd_->buf4_close(&D);


        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
        global_dpd_->buf4_scmcopy(&D, PSIF_LIBTRANS_DPD, "2 tIjAb - tIjBa (MP2)", 2);
        global_dpd_->buf4_sort_axpy(&D, PSIF_LIBTRANS_DPD, pqsr, 0, 5, "2 tIjAb - tIjBa (MP2)", -1);
        global_dpd_->buf4_close(&D);

        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa (MP2)");
        global_dpd_->buf4_dirprd(&Dijab, &D);
        global_dpd_->buf4_close(&D);
        global_dpd_->buf4_close(&Dijab);

        global_dpd_->file2_init(&Density, PSIF_LIBTRANS_DPD, 0, 1, 1, "d(a,b)");
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "tIjAb(MP2)");
        global_dpd_->buf4_init(&D1, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa (MP2)");
        global_dpd_->contract442(&D1, &D, &Density, 2, 2, 1.0, 0); // don't know why it was 2.0 instead of 1.0
        global_dpd_->buf4_close(&D);
        global_dpd_->buf4_close(&D1);
        global_dpd_->file2_close(&Density);

    }



    void FVNO::preppert(){

       /* Right now hard-coded just for dipole perturbations */

       cart_.push_back("X"); cart_.push_back("Y"); cart_.push_back("Z");
       mints_ = std::shared_ptr<MintsHelper>(new MintsHelper(ref_wfn_->basisset(), options_, 0));
       dipole_ = mints_->ao_dipole();
       for(int p=0; p<3; p++){
         dipole_[p]->transform(C_);
         sort_pert("mu", cart_[p], dipole_[p]); // hard-coded for dipole perturbations
        }   
    }

    void FVNO::sort_pert(std::string pert, std::string cart, SharedMatrix pert_mat){
   
       dpdfile2 f, Dia; 
   
       std::string lbl = pert + cart + "_IJ";
       global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 0, lbl.c_str());
       global_dpd_->file2_mat_init(&f);
       for(int i=0; i<occ_; i++)
         for(int j=0; j<occ_; j++) 
           f.matrix[0][i][j] = pert_mat->get(i,j);
       global_dpd_->file2_mat_wrt(&f);
       global_dpd_->file2_mat_close(&f);
       global_dpd_->file2_close(&f);

       lbl = pert + cart + "_AB";
       global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 1, 1, lbl.c_str());
       global_dpd_->file2_mat_init(&f);
       for(int a=0; a<vir_; a++) 
         for(int b=0; b<vir_; b++) 
           f.matrix[0][a][b] = pert_mat->get(a+occ_,b+occ_); 
       global_dpd_->file2_mat_wrt(&f);
       global_dpd_->file2_mat_close(&f);
       global_dpd_->file2_close(&f);

       lbl = pert + cart + "_IA";
       global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl.c_str());
       global_dpd_->file2_mat_init(&f);
       for(int i=0; i<occ_; i++)
         for(int a=0; a<vir_; a++) 
           f.matrix[0][i][a] = pert_mat->get(i, a+occ_);
       global_dpd_->file2_mat_wrt(&f);
       global_dpd_->file2_mat_close(&f);
       global_dpd_->file2_close(&f);
    
       /* need only mu_ia/D_ia actually */
       global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl.c_str());
       global_dpd_->file2_init(&Dia, PSIF_LIBTRANS_DPD, 0, 0, 1, "Dia");
       global_dpd_->file2_dirprd(&Dia, &f);
       global_dpd_->file2_close(&f);
       global_dpd_->file2_close(&Dia);
    }


    void FVNO::pert_density_vv(){

       dpdfile2 f, fc, D;
       std::string lbl; 

       /* lets try to finish it tonight 
          singles contribution to pert_density 
       */  
 
       global_dpd_->file2_init(&D, PSIF_LIBTRANS_DPD, 0, 1, 1, "d(a,b)");
       for(int p=0; p<3; p++){
         lbl = "mu" + cart_[p] + "_IA"; // hard-coded for "mu"
         global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl.c_str());
         global_dpd_->file2_copy(&f, PSIF_LIBTRANS_DPD, "tmp");
         global_dpd_->file2_init(&fc, PSIF_LIBTRANS_DPD, 0, 0, 1, "tmp");
         global_dpd_->contract222(&f, &fc, &D, 1, 1, 1, 1);
         global_dpd_->file2_close(&f);
         global_dpd_->file2_close(&fc);
       }
       global_dpd_->file2_close(&D);
       psio_->close(PSIF_LIBTRANS_DPD, 1);

    }

    void FVNO::truncate_VNOs(){
    
        /* [CVMO : Canonical virtual MOs] [VNO: Virtual Natural orbitals]  [ONs: Occupation numbers] 
           I would like to gather the density matrix from the dpd structute here inside this function 
           insread of gs_mp2_density_vv()
        */

        bool loaded = gs_density_->load(psio_, PSIF_LIBTRANS_DPD, "d(a,b)", vir_);

        SharedVector zero_vec(new Vector("zero vector", vir_));
        gs_density_->diagonalize(VNOs_, ONs_, descending);
        VNOs_->print();
        ONs_->print();
        
        if (frz_vno_){
          for(int a=0;a<frz_vno_;a++)
            VNOs_->set_column(0, vir_-a-1,zero_vec);
        }
        else if (on_cutoff_){
          int count = 0;
          for(int a=0; a<vir_; a++){
            if (fabs(ONs_->get(a)) > on_cutoff_)
                count++;   
            else break;
        } 
          frz_vno_ = count;
          outfile->Printf("\n\t Cutoff for ONs: %20.10lf\n",on_cutoff_);

          for(int a=0; a<frz_vno_; a++)
            VNOs_->set_column(0, vir_-a-1,zero_vec);
        }

          outfile->Printf("\n\t Number of frozen VNOs: %d\n",frz_vno_);
          VNOs_->print();   

    }

    void FVNO::semicanonicalize_VNOs(){
    
        SharedMatrix F_v(new Matrix("vir-vir block of fock", vir_, vir_));            
        SharedMatrix C_v(new Matrix("vir columns of MO C ", nmo_, vir_));            
        SharedVector SCanVOrbEnergies(new Vector("SemiCanonical Virtual Orbital energies", vir_));            
        SharedMatrix VNOs_VSCan(new Matrix("VNOs to VSCan", vir_, vir_));            
        SharedMatrix tmp_m(new Matrix("temporary matrix", nmo_, vir_));            

        for(int a=0; a<vir_; a++)  
          for(int b=0; b<vir_; b++)  
            F_v->set(a, b, F_mo_->get(a+occ_, b+occ_));

        for(int p=0; p<nmo_; p++)  
          for(int a=0; a<vir_; a++)  
            C_v->set(p, a, C_->get(p, a+occ_));  

        F_v->transform(VNOs_);
        F_v->diagonalize(VNOs_VSCan, SCanVOrbEnergies, descending); 
        tmp_m->gemm(0, 0, 1, C_v, VNOs_, 0); 
        SCan_VMO_->gemm(0, 0, 1, tmp_m, VNOs_VSCan, 0); 
      
        for (int p=0; p<nmo_;  p++)
          for (int a=0; a<vir_; a++)
             C_->set(p, a+occ_, SCan_VMO_->get(p,a));
        C_->print();
    }

}}
