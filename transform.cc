#include "psi4/libdpd/dpd.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/wavefunction.h"
#include"transform.h"
// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

namespace psi { namespace fvno{

void transform_to_mo(SharedWavefunction ref_wfn, std::shared_ptr<PSIO> psio)
{
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
    IntegralTransform ints(ref_wfn, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly);
    ints.set_keep_iwl_so_ints(true);
    ints.transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());

    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    psio->open(PSIF_CC_DINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_sort(&K, PSIF_CC_DINTS, prqs, ID("[O,O]"), ID("[V,V]"), "D <ij|ab>");
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, 1);
    psio->close(PSIF_CC_DINTS, 1);

 }

}} // end namespace 
