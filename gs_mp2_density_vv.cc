#include "psi4/libdpd/dpd.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/wavefunction.h"
#include "fvno.h"

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

namespace psi { namespace fvno{

SharedMatrix gs_mp2_density_vv(SharedWavefunction ref_wfn, std::shared_ptr<PSIO> psio)
{
     int nirrep  = ref_wfn->nirrep();
     int nmo= ref_wfn->nmo();
     int occ = ref_wfn->doccpi()[0]; // hard-coded for symmetry C1
     int vir = nmo-occ;
     dpdbuf4 D,D1;
     dpdfile2 Density;


     SharedVector epsilon (new Vector(" orbital energies",nmo));
     SharedMatrix C(new Matrix("Canonical MO Coefficients", nmo, nmo));
     SharedMatrix F(new Matrix("Fock Matrix Canonical MO basis", nmo, nmo));
     SharedMatrix gs_density(new Matrix("ground state mp2 density MO basis (vir-vir)", vir, vir));
     epsilon->copy(ref_wfn->epsilon_a()->clone());
     F->copy(ref_wfn->Fa()->clone());
     C->copy(ref_wfn->Ca()->clone());
     F->transform(C);

     std::vector<std::shared_ptr<MOSpace> > spaces;
     spaces.push_back(MOSpace::occ);
     spaces.push_back(MOSpace::vir);
     IntegralTransform ints(ref_wfn, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly);
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
     psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
     global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "D <ij|ab>");
     global_dpd_->buf4_close(&K);

     global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
     global_dpd_->buf4_copy(&D, PSIF_LIBTRANS_DPD , "D' <ij|ab>");
     global_dpd_->buf4_close(&D);

     global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, 0, 5, 0, 5, 0, "D' <ij|ab>");
     global_dpd_->buf4_mat_irrep_init(&D, 0);
     global_dpd_->buf4_mat_irrep_rd(&D, 0);
         for(int i=0,ij=0;i<occ;i++)
         for(int j=0;j<occ;j++,ij++){
         for(int a=0,ab=0;a<vir;a++)
         for(int b=0;b<vir;b++,ab++)
        D.matrix[0][ij][ab] /= (epsilon->get(i) + epsilon->get(j) - F->get(a+occ,a+occ) - F->get(b+occ,b+occ));
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
        for(int i=0,ij=0;i<occ;i++)
        for(int j=0;j<occ;j++,ij++){
        for(int a=0,ab=0;a<vir;a++)
        for(int b=0;b<vir;b++,ab++)
           D.matrix[0][ij][ab] /= (epsilon->get(i) + epsilon->get(j) - F->get(a+occ,a+occ) - F->get(b+occ,b+occ));
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
      for(int a=0;a<vir;a++)
         for(int b=0;b<vir;b++)
           gs_density->set(a,b,Density.matrix[0][a][b]);
     global_dpd_->file2_mat_close(&Density);
     global_dpd_->file2_close(&Density);

    return gs_density;

 }

}} // end namespace 
