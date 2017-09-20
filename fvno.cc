/*
 * @BEGIN LICENSE
 *
 * fvno by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "fvno.h"

namespace psi{ namespace fvno{


extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "FVNO" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /* The number of virtual NOs to be removed */
        options.add_int("FRZ_VNO", 0);
    }

    return true;
}


extern "C"
SharedWavefunction fvno(SharedWavefunction ref_wfn, Options& options)
{
    /*
     * This plugin shows a simple way of obtaining MO basis integrals, directly from a DPD buffer.  It is also
     * possible to generate integrals with labels (IWL) formatted files, but that's not shown here.
     */
    int print = options.get_int("PRINT");
    int frz_vno = options.get_int("FRZ_VNO");

    // Grab the global (default) PSIO object, for file I/O
    std::shared_ptr<PSIO> psio(_default_psio_lib_);

    // Have the reference (SCF) wavefunction, ref_wfn
    if(!ref_wfn) throw PSIEXCEPTION("SCF has not been run yet!");

     std::shared_ptr<FVNO> FVNO_obj(new FVNO(ref_wfn, psio, options)); 

     /* Transform the ao integrals to the MO basis <ij|ab> type integrals 
       and construct the virtual-virtual block of ground state MP2 
       second order reduced density matrix. In spin orbitals,
       D(a,b) = 0.5 * \sum_ijc t^{ac}_{ij} * t^{bc}_{ij}
       In spin adpated form: 
       D(a,b) = \sum_ijc (2.0 * t^{ac}_{ij} - t^{ca}{ij}) * t^{bc}_{ij}
    */

     FVNO_obj->gs_mp2_density_vv();

    /* Truncation of virtual NOs based on a given cutoff
       or number of frozen virtual NOs specified in the input
    */

     FVNO_obj->truncate_VNOs();

    /* Semi-canonicalize the VNO basis to speed up the convergence
       of CC amplitude equations.
    */

     FVNO_obj->semicanonicalize_VNOs();


/* 
    Other steps to follow below:
    1. Create a pertubed density based on guesses of X1s and X2s
    2. Write down the expression for the perturbed demnsity as well
    3. Don't forget to analyze the sparsity pattern of singles and
       doubles amplitudes in ccenergy, cclambda and ccresponse modules 
       in all these kinds of basis.
    4. add an option of localizing the occupied orbitals
    5. add an option for the canonical scheme where high energy 
       virtuals should be removed.
    6. add an option for the finite difference scheme so that the field
       can be added here after the HF stage. 
    7. maybe have different functions for the below. Include only the logic here
       in the main file.
    A. FVNO
    B. localization of occupied orbitals
    C. FCMO
    D. FVNO ++
    E. FIN_DIFF
*/
    psio->close(PSIF_LIBTRANS_DPD, 1);
    return ref_wfn;
}

}} // End Namespaces
