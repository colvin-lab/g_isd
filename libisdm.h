/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

/* Modified 2012-09-26
 * 
 * Joshua L. Phillips - jphillips@lanl.gov
 * T-6/CNLS - Los Alamos National Laboratory
 *
 * This code was originally taken from the file src/tools/gmx_gyrate.c
 * from GROMACS release 4.5.5 and modified for performing shape
 * analysis as described in:
 *
 * Dima, R. I., & Thirumalai, D. (2004). Asymmetry in the shapes of
 * folded and denatured states of proteins. The Journal of Physical
 * Chemistry B, 108(21), 6564-6570. doi:10.1021/jp037128y   
 */ 


#include <math.h>
#include <string.h>

#include <gromacs/statutil.h>
#include <gromacs/sysstuff.h>
#include <gromacs/typedefs.h>
#include <gromacs/smalloc.h>
#include <gromacs/macros.h>
#include <gromacs/vec.h>
#include <gromacs/pbc.h>
#include <gromacs/copyrite.h>
#include <gromacs/futil.h>
#include <gromacs/statutil.h>
#include <gromacs/index.h>
#include <gromacs/mshift.h>
#include <gromacs/xvgr.h>
#include <gromacs/rmpbc.h>
#include <gromacs/txtdump.h>
#include <gromacs/tpxio.h>
#include <gromacs/gstat.h>
#include <gromacs/gmx_ana.h>



extern real calc_drms(int iatoms, rvec frame[], rvec rframe[], real drms[]);

extern real calc_rg(int iatoms, rvec frame[]);

extern real calc_ang(int iatoms, rvec frame[], rvec rframe[], real ang[]);

extern real calc_dih(int iatoms, rvec frame[], rvec rframe[], real dihs[]);

extern real calc_angdih(int iatoms, rvec frame[], rvec rframe[], 
                        real angdih[]);

extern real calc_phipsi(int iatoms, rvec frame[], rvec rframe[], 
                        real phipsi[]);

extern real calc_msd(int iatoms, rvec x[], rvec xref[], real msd[]);

extern real calc_scaled_msd(int iatoms, rvec x[], rvec xref[], real msd[]);

extern real calc_mirror_msd(int iatoms, rvec frame[], rvec rframe[], 
                            real mirror_msd[]);

extern real calc_corr_atoms(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_corr_angs(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_mammoth(int iatoms, rvec frame[], rvec rframe[], int *rnum);

extern real calc_esa(int iatoms, rvec frame[], rvec rframe[]);

extern real call_ISDM(int iatoms, rvec cframe[], rvec rframe[], real diff[], 
                        const char *ISDM);