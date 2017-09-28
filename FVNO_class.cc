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

        outfile->Printf("\nConstructor starting:\n");
        psio_ = psio;
        ref_wfn_ = ref_wfn;
        options_ = options;
        nirrep_  = ref_wfn_->nirrep();
        nso_ = ref_wfn_->nso();
        nmo_ = ref_wfn_->nmo();
        occ_ = ref_wfn_->doccpi()[0]; // hard-coded for symmetry C1
        vir_ = nmo_-occ_;
        print_ = false;
        frz_vno_ = options.get_int("FRZ_VNO");
        my_frozen_uocc_ = options.get_int("MY_FROZEN_UOCC");
        on_cutoff_ = options.get_double("ON_CUTOFF");
        vno_basis_ = options.get_bool("VNO_BASIS"); 
        //outfile->Printf("\nVNO_BASIS: %d\n", vno_basis_);
        semicanonical_basis_ = options.get_bool("SEMICANONICAL_BASIS"); 
        canonical_basis_ = options.get_bool("CANONICAL_BASIS"); 
        prop_correction_ = options.get_bool("PROP_CORRECTION"); 
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
        SCan_VMO_ = SharedMatrix(new Matrix(" Semicanonical basis (CVMO->VNO->SCVMO transformation)", nmo_, vir_));
        outfile->Printf("\nConstructor over:\n");
    }

    FVNO::~FVNO()
    {   
    }

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
         for(int i=0; i<occ_; i++)
           for(int a=0; a<vir_; a++)
             Dia.matrix[0][i][a] = 1.0/(epsilon_->get(i) - F_mo_->get(a+occ_,a+occ_));
        global_dpd_->file2_mat_wrt(&Dia);
        global_dpd_->file2_mat_close(&Dia);
        global_dpd_->file2_close(&Dia);
        

        }

        
    void FVNO::gs_mp2_density_vv(){

        dpdbuf4 D, D1, Dia, Dijab;
        dpdfile2 Density;
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
        global_dpd_->buf4_init(&Dijab, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "Dijab");
        global_dpd_->buf4_dirprd(&Dijab, &D);
        global_dpd_->buf4_copy(&D, PSIF_LIBTRANS_DPD , "tIjAb (MP2)");
        global_dpd_->buf4_close(&Dijab);
        global_dpd_->buf4_close(&D);

        if (print_){
        SharedMatrix T2(new Matrix("T2 in matrix form (vir-vir)", vir_, vir_));
        T2->zero();
        SharedVector T2_vec(new Vector("T2 in vector form (vir)", vir_));
        T2_vec->zero();

        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "tIjAb (MP2)");
        global_dpd_->buf4_mat_irrep_init(&D, 0);
        global_dpd_->buf4_mat_irrep_rd(&D, 0);
        for(int a=0,ab=0; a<vir_; a++)
          for(int b =0; b<vir_; b++,ab++)
            for(int i=0,ij=0; i<occ_;i++ )
              for(int j=0; j<occ_;j++,ij++)
                T2->add(a,b, fabs(D.matrix[0][ij][ab])) ;
        global_dpd_->buf4_mat_irrep_close(&D,0);
        global_dpd_->buf4_close(&D);

        FILE * ptr = fopen("T2.dat","w");
        for(int a=0; a<vir_; a++){
          for(int b=0; b<vir_; b++){
            fprintf(ptr,"%20.15lf", T2->get(a,b));
            }
            fprintf(ptr,"\n");
        }

        for(int a=0; a<vir_; a++)
          for(int b =0; b<vir_; b++)
            T2_vec->add(a, fabs(T2->get(a,b))) ;
          T2_vec->print();

        }  

        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "tIjAb (MP2)");
        global_dpd_->buf4_scmcopy(&D, PSIF_LIBTRANS_DPD, "2 tIjAb - tIjBa (MP2)", 2);
        global_dpd_->buf4_sort_axpy(&D, PSIF_LIBTRANS_DPD, pqsr, 0, 5, "2 tIjAb - tIjBa (MP2)", -1);
        global_dpd_->buf4_close(&D);

        global_dpd_->file2_init(&Density, PSIF_LIBTRANS_DPD, 0, 1, 1, "d(a,b)");
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "tIjAb (MP2)");
        global_dpd_->buf4_init(&D1, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa (MP2)");
        global_dpd_->contract442(&D1, &D, &Density, 2, 2, 2.0, 0); // don't know why it was 2.0 instead of 1.0, lets play with this
        global_dpd_->buf4_close(&D);
        global_dpd_->buf4_close(&D1);
        global_dpd_->file2_close(&Density);

        global_dpd_->file2_init(&Density, PSIF_LIBTRANS_DPD, 0, 1, 1,"d(a,b)");
        global_dpd_->file2_mat_init(&Density);
        global_dpd_->file2_mat_rd(&Density);
        for(int a=0; a<vir_; a++)
          for(int b=0; b<vir_; b++)
            gs_density_->set(a,b,Density.matrix[0][a][b]);
        global_dpd_->file2_mat_close(&Density);
        global_dpd_->file2_close(&Density);

       combined_density_ = SharedMatrix(new Matrix("combined density MO basis (vir-vir)", vir_, vir_));
       combined_density_->copy(gs_density_->clone());
        //gs_density_->print();

    }


    /* Eventually, prepert should take a perturbation as an argument.
       Right now its hard-coded just for dipole perturbations */


    void FVNO::preppert(Options& options){  

       cart_.push_back("X"); cart_.push_back("Y"); cart_.push_back("Z");
       std::string lbl, lbl1 ;
       dpdfile2 Dia, f, f1;

       if (options.get_bool("CONSTRUCT_X1")){
       mints_ = std::shared_ptr<MintsHelper>(new MintsHelper(ref_wfn_->basisset(), options, 0));
       dipole_ = mints_->ao_dipole();
       for(int p=0; p<3; p++){
         dipole_[p]->transform(C_);
         /* 
         Sorting the perturbation matrices into occupied and
         virtual blocks in dpdfile2 structures.
         */
         lbl = "Mu" + cart_[p] + "_IJ";
         global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 0, lbl.c_str());
         global_dpd_->file2_mat_init(&f);
         for(int i=0; i<occ_; i++)
           for(int j=0; j<occ_; j++)
             f.matrix[0][i][j] = dipole_[p]->get(i,j);
         global_dpd_->file2_mat_wrt(&f);
         global_dpd_->file2_mat_close(&f);
         global_dpd_->file2_close(&f);

         
         /*global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 0, lbl.c_str());
         global_dpd_->file2_mat_init(&f);
         global_dpd_->file2_mat_rd(&f);
         global_dpd_->file2_mat_print(&f, "outfile");
         global_dpd_->file2_close(&f);*/

         //outfile->Printf("\nMU_X over:\n");

         lbl = "Mu" + cart_[p] + "_AB";
         global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 1, 1, lbl.c_str());
         global_dpd_->file2_mat_init(&f);
         for(int a=0; a<vir_; a++)
           for(int b=0; b<vir_; b++)
             f.matrix[0][a][b] = dipole_[p]->get(a+occ_,b+occ_);
         global_dpd_->file2_mat_wrt(&f);
         global_dpd_->file2_mat_close(&f);
         global_dpd_->file2_close(&f);

         lbl = "Mu" + cart_[p] + "_IA";
         global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl.c_str());
         global_dpd_->file2_mat_init(&f);
         for(int i=0; i<occ_; i++)
           for(int a=0; a<vir_; a++)
             f.matrix[0][i][a] = dipole_[p]->get(i, a+occ_);
         global_dpd_->file2_mat_wrt(&f);
         global_dpd_->file2_mat_close(&f);
         global_dpd_->file2_close(&f);
   
         lbl = "Mu" + cart_[p] + "_IA";
         global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl.c_str());
         lbl = "X_Mu_" + cart_[p] + "_IA";
         global_dpd_->file2_copy(&f, PSIF_LIBTRANS_DPD, lbl.c_str());
         global_dpd_->file2_close(&f);
          

         /* need only mu_ia/D_ia actually */
         global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl.c_str());
         global_dpd_->file2_init(&Dia, PSIF_LIBTRANS_DPD, 0, 0, 1, "Dia");
         global_dpd_->file2_dirprd(&Dia, &f);
         global_dpd_->file2_close(&f);
         global_dpd_->file2_close(&Dia);

         if (print_){
         global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl.c_str());
         global_dpd_->file2_mat_init(&f);
         global_dpd_->file2_mat_rd(&f);

         FILE * ptr = fopen(lbl.c_str(),"w"); 
         for(int i=0; i<occ_; i++){
           for(int a=0; a<vir_; a++){
             fprintf(ptr,"%20.15lf", f.matrix[0][i][a]);
             //fprintf(ptr,"%20.15lf", dipole_[p]->get(i,a+occ_));
            }
           fprintf(ptr,"\n");
         }
         //global_dpd_->file2_mat_print(&f, "outfile");
         global_dpd_->file2_close(&f);
            }
        }
      }

       if (options.get_bool("READ_X1")){
         outfile->Printf("\nInside read_X1\n");
         dpdfile2 f; 
         std::string lbl, lbl1;
         for(int p=0; p<3; p++){
           lbl = "X_Mu_" + cart_[p] + "_IA"; // expects the fole to read from with this name
           global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl.c_str());
           global_dpd_->file2_mat_init(&f);
           FILE *ptr;
           double value;
           ptr = fopen(lbl.c_str(), "r");
           for(int i=0; i<occ_; i++)
             for(int a=0; a<vir_; a++){
               fscanf(ptr, "%lf", &f.matrix[0][i][a]);
               }
             fclose(ptr);
           global_dpd_->file2_mat_wrt(&f);
           global_dpd_->file2_mat_close(&f);
           global_dpd_->file2_close(&f);        
           }
         outfile->Printf("\nREAD_X1 OVER\n");
          }

       if(prop_correction_) {
         double four_X1_Mu1 = 0;
         for(int p=0; p<3; p++){
           lbl = "Mu" + cart_[p] + "_IA";
           lbl1 = "X_Mu_" + cart_[p] + "_IA";
           global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl.c_str());
           global_dpd_->file2_init(&f1, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl1.c_str());
           four_X1_Mu1 += 4.0 * global_dpd_->file2_dot(&f, &f1);
           global_dpd_->file2_close(&f);
           global_dpd_->file2_close(&f1);
           }
         outfile->Printf("\n 4X1Mu1 -> %20.15lf\n", four_X1_Mu1);
          }
        }

         /*if(prop_correction_) {

           double value = 0;
           for(int p=0; p<3; p++){

             // *** Mu * X1 *** /

             lbl = "Mu" + cart_[p] + "_IA";
             lbl1 = "X_Mu_" + cart_[p] + "_IA";
             global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl.c_str());
             global_dpd_->file2_init(&f1, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl1.c_str());
             value += 4.0 * global_dpd_->file2_dot(&f, &f1);
             global_dpd_->file2_close(&f);
             global_dpd_->file2_close(&f1);


             / *** L2 * MuBAR * X1 + L2 * MuBAR * X2 ***/

            /*double value1 = 0;

            dpdfile2 xc, mu, X1;
            global_dpd_->file2_init(&xc, PSIF_LIBTRANS_DPD, 0, 0, 0, "XC_IJ");
            global_dpd_->file2_init(&mu, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl);
            global_dpd_->file2_init(&X1, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl1);
            global_dpd_->contract222(&X1, &mu, &xc, 0, 0, 2.0, 0.0);
            global_dpd_->file2_close(&X1);
            global_dpd_->file2_close(&mu);
            global_dpd_->file2_close(&xc);

            global_dpd_->file2_init(&lt, PSIF_LIBTRANS_DPD, 0, 0, 0, "Lt_IJ");
            global_dpd_->file2_init(&xc, PSIF_LIBTRANS_DPD, 0, 0, 0, "XC_IJ");
            value1 += global_dpd_->file2_dot(&lt, &xc);
            global_dpd_->file2_close(&xc);
            global_dpd_->file2_close(&lt);

            global_dpd_->buf4_init(&z2, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
            global_dpd_->buf4_scm(&z2, 0);
            global_dpd_->buf4_close(&z2);

            sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
            global_dpd_->file2_init(&X1, PSIF_LIBTRANS_DPD, irrep_x, 0, 1, lbl);

            global_dpd_->buf4_init(&z2, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");

            global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
            sprintf(lbl, "%sBAR_MbIj", pert_c);
            global_dpd_->buf4_init(&mu2, PSIF_CC_LR, irrep_c, 10, 0, 10, 0, 0, lbl);
            global_dpd_->contract244(&X1, &mu2, &Z, 0, 0, 1, 1, 0);
            global_dpd_->buf4_close(&mu2);
            global_dpd_->buf4_axpy(&Z, &z2, -1);
            global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
            global_dpd_->buf4_close(&Z);
            global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
            global_dpd_->buf4_axpy(&Z, &z2, -1);
            global_dpd_->buf4_close(&Z);

            global_dpd_->file2_close(&X1);

            sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega);
            global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

            global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");

            sprintf(lbl, "%sBAR_AE", pert_c);
            global_dpd_->file2_init(&mu1, PSIF_LIBTRANS_DPD, irrep_c, 1, 1, lbl);
            global_dpd_->contract424(&X2, &mu1, &Z, 3, 1, 0, 1, 0);
            global_dpd_->file2_close(&mu1);
            global_dpd_->buf4_axpy(&Z, &z2, 1);
            global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
            global_dpd_->buf4_close(&Z);
            global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
            global_dpd_->buf4_axpy(&Z, &z2, 1);
            global_dpd_->buf4_close(&Z);

            global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");

            sprintf(lbl, "%sBAR_MI", pert_c);
            global_dpd_->file2_init(&mu1, PSIF_LIBTRANS_DPD, irrep_c, 0, 0, lbl);
            global_dpd_->contract244(&mu1, &X2, &Z, 0, 0, 0, 1, 0);
            global_dpd_->file2_close(&mu1);
            global_dpd_->buf4_axpy(&Z, &z2, -1);
            global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
            global_dpd_->buf4_close(&Z);
            global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
            global_dpd_->buf4_axpy(&Z, &z2, -1);
            global_dpd_->buf4_close(&Z);

            global_dpd_->buf4_close(&X2);

            global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
            polar += global_dpd_->buf4_dot(&l2, &z2);
            global_dpd_->buf4_close(&l2);

            global_dpd_->buf4_close(&z2);

           }
         outfile->Printf("\n 4X1Mu1 -> %20.15lf\n", value);
          }  */

    void FVNO::pert_density_vv(Options & options){

       dpdfile2 f, fc, D;
       dpdbuf4 fbar, fbar1, t2, Dijab;
       std::string lbl; 
       std::string lbl1; 

       /* 
          Singles contribution to pert_density 
       */  

       pert_density_singles_ = SharedMatrix(new Matrix("singles perturbed density MO basis (vir-vir)", vir_, vir_));
       pert_density_doubles_ = SharedMatrix(new Matrix("doubles perturbed density MO basis (vir-vir)", vir_, vir_));
 
       outfile->Printf("\nConstructing singles contribution to pert_density\n");
       for(int p=0; p<3; p++){
         lbl = "X_Mu_" + cart_[p] + "_IA"; // hard-coded for "mu"
         lbl1 = "tmp_mu" + cart_[p] + "_IA"; 
         global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl.c_str());
         global_dpd_->file2_copy(&f, PSIF_LIBTRANS_DPD, lbl1.c_str());
         global_dpd_->file2_close(&f);
         }
        

       global_dpd_->file2_init(&D, PSIF_LIBTRANS_DPD, 0, 1, 1, "d_pert_singles(a,b)");
       for(int p=0; p<3; p++){
         lbl = "X_Mu_" + cart_[p] + "_IA"; // hard-coded for "mu"
         lbl1 = "tmp_mu" + cart_[p] + "_IA"; 
         global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl.c_str());
         global_dpd_->file2_init(&fc, PSIF_LIBTRANS_DPD, 0, 0, 1, lbl1.c_str());
         global_dpd_->contract222(&f, &fc, &D, 1, 1, 1, 1);
         global_dpd_->file2_close(&f);
         global_dpd_->file2_close(&fc);
       }
       global_dpd_->file2_close(&D);

        
        global_dpd_->file2_init(&D, PSIF_LIBTRANS_DPD, 0, 1, 1,"d_pert_singles(a,b)");
        global_dpd_->file2_mat_init(&D);
        global_dpd_->file2_mat_rd(&D);
        for(int a=0; a<vir_; a++)
          for(int b=0; b<vir_; b++)
            pert_density_singles_->set(a,b,D.matrix[0][a][b]);
        global_dpd_->file2_mat_close(&D);
        global_dpd_->file2_close(&D);


        //pert_density_singles_->print();


        /* Doubles contribution to pert_density */
         

       for(int p=0; p<3; p++){
         lbl = "Mu" + cart_[p] + "_IjAb";
         global_dpd_->buf4_init(&fbar, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, lbl.c_str());

         lbl = "Mu" + cart_[p] + "_AB"; 
         global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 1, 1, lbl.c_str());
         global_dpd_->buf4_init(&t2, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "tIjAb (MP2)");
         global_dpd_->contract424(&t2, &f, &fbar, 3, 1, 0, 1, 0);
         global_dpd_->contract244(&f, &t2, &fbar, 1, 2, 1, 1, 1);
         global_dpd_->buf4_close(&t2);
         global_dpd_->file2_close(&f);

         lbl = "Mu" + cart_[p] + "_IJ"; 
         global_dpd_->file2_init(&f, PSIF_LIBTRANS_DPD, 0, 0, 0, lbl.c_str());
         global_dpd_->buf4_init(&t2, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "tIjAb (MP2)");
         global_dpd_->contract424(&t2, &f, &fbar, 1, 0, 1, -1, 1);
         global_dpd_->contract244(&f, &t2, &fbar, 0, 0, 0, -1, 1);
         global_dpd_->buf4_close(&t2);
         global_dpd_->file2_close(&f);

         global_dpd_->buf4_init(&Dijab, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "Dijab");
         global_dpd_->buf4_dirprd(&Dijab, &fbar);
         global_dpd_->buf4_close(&Dijab);

         lbl = "Mu" + cart_[p] + "_XIjAb (MP2)";
         global_dpd_->buf4_copy(&fbar, PSIF_LIBTRANS_DPD , lbl.c_str());
         global_dpd_->buf4_close(&fbar);


         global_dpd_->buf4_init(&fbar, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, lbl.c_str());
         lbl = "2 Mu" + cart_[p] + "_XIjAb (MP2)" + " - Mu" + cart_[p] + "_XIjBa (MP2)";
         global_dpd_->buf4_scmcopy(&fbar, PSIF_LIBTRANS_DPD, lbl.c_str(), 2);
         global_dpd_->buf4_sort_axpy(&fbar, PSIF_LIBTRANS_DPD, pqsr, 0, 5, lbl.c_str(), -1);
         global_dpd_->buf4_close(&fbar);

         pert_density_order_ = options.get_bool("PERT_DENSITY_ORDER") ;

         if (pert_density_order_ == 1){

         global_dpd_->buf4_init(&t2, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa (MP2)");
         global_dpd_->file2_init(&D, PSIF_LIBTRANS_DPD, 0, 1, 1, "d_pert_doubles(a,b)");
         lbl = "Mu" + cart_[p] + "_XIjAb (MP2)";
         global_dpd_->buf4_init(&fbar, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, lbl.c_str());
         global_dpd_->contract442(&t2, &fbar, &D, 2, 2, 2.0, 1);
         global_dpd_->buf4_close(&fbar);
         global_dpd_->buf4_close(&t2);
         global_dpd_->file2_close(&D);

         }

         if (pert_density_order_ == 2){

         global_dpd_->file2_init(&D, PSIF_LIBTRANS_DPD, 0, 1, 1, "d_pert_doubles(a,b)");
         lbl = "Mu" + cart_[p] + "_XIjAb (MP2)";
         lbl1 = "2 Mu" + cart_[p] + "_XIjAb (MP2)" + " - Mu" + cart_[p] + "_XIjBa (MP2)";
         global_dpd_->buf4_init(&fbar, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, lbl.c_str());
         global_dpd_->buf4_init(&fbar1, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, lbl1.c_str());
         global_dpd_->contract442(&fbar1, &fbar, &D, 2, 2, 2.0, 1); 
         global_dpd_->buf4_close(&fbar);
         global_dpd_->buf4_close(&fbar1);
         global_dpd_->file2_close(&D);

         }

         }

         global_dpd_->file2_init(&D, PSIF_LIBTRANS_DPD, 0, 1, 1,"d_pert_doubles(a,b)");
         global_dpd_->file2_mat_init(&D);
         global_dpd_->file2_mat_rd(&D);
         for(int a=0; a<vir_; a++)
           for(int b=0; b<vir_; b++)
             pert_density_doubles_->set(a,b,D.matrix[0][a][b]);
         global_dpd_->file2_mat_close(&D);
         global_dpd_->file2_close(&D);

        if (print_){

        SharedMatrix X2(new Matrix("X2 in matrix form (vir-vir)", vir_, vir_));

       for(int p=0; p<3; p++){
        X2->zero();
        lbl = "Mu" + cart_[p] + "_XIjAb (MP2)";
        global_dpd_->buf4_init(&t2, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, lbl.c_str());
        global_dpd_->buf4_mat_irrep_init(&t2, 0);
        global_dpd_->buf4_mat_irrep_rd(&t2, 0);
        for(int a=0,ab=0; a<vir_; a++)
          for(int b =0; b<vir_; b++,ab++)
            for(int i=0,ij=0; i<occ_;i++ )
              for(int j=0; j<occ_;j++,ij++)
                    X2->add(a,b, fabs(t2.matrix[0][ij][ab])) ;
        global_dpd_->buf4_mat_irrep_close(&t2,0);
        global_dpd_->buf4_close(&t2);

        FILE * ptr = fopen(lbl.c_str(),"w");
        for(int a=0; a<vir_; a++){
          for(int b=0; b<vir_; b++){
            fprintf(ptr,"%20.15lf", X2->get(a,b));
            }
            fprintf(ptr,"\n");
            }
        }
      }
    }

    void FVNO::final_density(Options& options){

        /* Combine the ground state density with singles and doubles perturbed densities */

       bool only_gs_density = options.get_bool("ONLY_GS_DENSITY") ;
       bool combined_pert_density = options.get_bool("COMBINED_PERT_DENSITY") ;
       bool only_pert_density = options.get_bool("ONLY_PERT_DENSITY") ;
       bool pert_singles = options.get_bool("PERT_DENSITY_SINGLES") ;
       bool pert_doubles = options.get_bool("PERT_DENSITY_DOUBLES") ;

       if (only_gs_density) {
          combined_density_ = SharedMatrix(new Matrix("combined density MO basis (vir-vir)", vir_, vir_));
          combined_density_->copy(gs_density_->clone());
        }
       else if (combined_pert_density) {
          combined_density_ = SharedMatrix(new Matrix("combined density MO basis (vir-vir)", vir_, vir_));
          combined_density_->copy(gs_density_->clone());
          if (pert_singles) {combined_density_->add(pert_density_singles_); outfile->Printf("\n PERT_SINGLES\n");}
          if (pert_doubles) {combined_density_->add(pert_density_doubles_); outfile->Printf("\n PERT_DOUBLES\n");}
         combined_density_->scale(0.5);
       } 

       else if (only_pert_density) {
          combined_density_ = SharedMatrix(new Matrix("combined density MO basis (vir-vir)", vir_, vir_));
          combined_density_->zero();
          double value = 0;
          if (pert_singles) {
            /* Just experimenting //
            pert_density_singles_->zero();
            for(int a=0; a<vir_; a++)
              for(int b=0; b<vir_; b++)
                for(int i=0; i<occ_; i++){
                  value = dipole_[2]->get(i,a+occ_) * dipole_[2]->get(i,b+occ_); // hard-coded for muZ
                  value /= ((epsilon_->get(i) - F_mo_->get(a+occ_,a+occ_)) * (epsilon_->get(i) - F_mo_->get(b+occ_,b+occ_)));
                  combined_density_->add(a, b, value);
                }
                  //combined_density_->print();*/
                combined_density_->add(pert_density_singles_); outfile->Printf("\n PERT_SINGLES\n");
            }
          if (pert_doubles) {combined_density_->add(pert_density_doubles_); outfile->Printf("\n PERT_DOUBLES\n");}
          }
        }

    void FVNO::truncate_VNOs(){
    
        /* [CVMO : Canonical virtual MOs] [VNO: Virtual Natural orbitals]  [ONs: Occupation numbers] 
           I would like to gather the density matrix from the dpd structute here inside this function 
           insread of gs_mp2_density_vv()
        */

        //bool loaded = gs_density_->load(psio_, PSIF_LIBTRANS_DPD, "d(a,b)", vir_); 

        SharedVector zero_vec(new Vector("zero vector", vir_));
        combined_density_->diagonalize(VNOs_, ONs_, descending);
        VNOs_->print();
        ONs_->print();
        zero_vec->zero();;
        
        /* sort the VNOs based on their absolute values of eigenvalues */

         if (pert_density_order_ == 1){

         SharedMatrix No_temp(new Matrix("temporary matrix (full NO) ",vir_,vir_));
         SharedVector indx(new Vector("index for occupation numbers",vir_));


         for(int a=0; a<vir_; a++)
           ONs_->set(a,fabs(ONs_->get(a)));


        double X;
        int IX;

        for(int a=0; a<vir_;a++)
          indx->set(a,a);

        for(int a=0; a<vir_; a++)
          { for(int b=0; b<vir_-a-1; b++)
            { if (ONs_->get(b) < ONs_->get(b+1))
                { X = ONs_->get(b+1);
                  ONs_->set(b+1, ONs_->get(b));
                  ONs_->set(b,X) ;
                  IX = indx->get(b+1);
                  indx->set(b+1,indx->get(b));
                  indx->set(b,IX) ;
              }
            }
          }

        indx->print();

        for (int a=0; a<vir_; a++)
           for (int b=0; b<vir_; b++)
               No_temp->set(a,b,VNOs_->get(a,indx->get(b)));

        for (int a=0; a<vir_; a++)
           for (int b=0; b<vir_; b++)
             VNOs_->set(a,b,No_temp->get(a,b));


          ONs_->print();   
          VNOs_->print();   

        }

        /* sorting procedure ends*/

        if (frz_vno_){
          for(int a=0; a<frz_vno_; a++)
            VNOs_->set_column(0, vir_-a-1, zero_vec);
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

    void FVNO::final_basis(){

        SharedMatrix C_v(new Matrix("vir columns of MO C ", nmo_, vir_));            
        SharedMatrix tmp_m(new Matrix("temporary matrix", nmo_, vir_));            

        if (vno_basis_){
           for(int p=0; p<nmo_; p++)  
             for(int a=0; a<vir_; a++)  
               C_v->set(p, a, C_->get(p, a+occ_));  

           tmp_m->gemm(0, 0, 1, C_v, VNOs_, 0);   

        for (int p=0; p<nmo_;  p++)
          for (int a=0; a<vir_; a++)
             C_->set(p, a+occ_, tmp_m->get(p,a));
        
        SharedMatrix tmp1_m(new Matrix("temporary matrix", vir_, vir_));            
        SharedMatrix tmp2_m(new Matrix("temporary matrix", vir_, vir_));            
        tmp2_m->copy(VNOs_->clone());
        tmp1_m->gemm(1, 0, 1, VNOs_, tmp2_m, 0);   
        //tmp1_m->print();

        }

        else if (semicanonical_basis_) {
        SharedMatrix F_v(new Matrix("vir-vir block of fock", vir_, vir_));            
        SharedVector SCanVOrbEnergies(new Vector("SemiCanonical Virtual Orbital energies", vir_));            
        SharedMatrix VNOs_VSCan(new Matrix("VNOs to VSCan", vir_, vir_));            

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
        }

        else if (canonical_basis_) {
        
           SharedVector zero_vec(new Vector("zero vector", nmo_));
           zero_vec->zero();

             for(int a=0; a<my_frozen_uocc_; a++)
               C_->set_column(0, vir_-a-1, zero_vec);
        }

        //C_->print();
    }

}}
