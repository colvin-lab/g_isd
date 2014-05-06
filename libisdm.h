/* 
 * 
 * Tim Connolly - tconnolly@ucmerced.edu
 * Copyright (c) 2014, Regents of the University of California
 * Released under BSD 2-Clause License (see "LICENSE" file)
 * 
 * This library of functions implements various comparison ISDMs
 * intended to quanitfy the difference between two input structures
 * from a trajectory. The gromacs libraries are required and
 * some functions are performed in libmammothmod.c and libesa.c.
 * 
 * See help output for citation information.
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



extern real calc_drms(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_sdrms(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_rg(int iatoms, rvec frame[]);

extern real calc_ang(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_dih(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_angdih(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_ang2(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_dih2(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_angdih2(int iatoms, rvec frame[], rvec rframe[], double ang_multi);

extern real calc_phipsi(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_msd(int iatoms, rvec x[], rvec xref[]);

extern real calc_srms(int iatoms, rvec x[], rvec xref[]);

extern real calc_mmsd(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_corr_atoms(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_corr_angs(int iatoms, rvec frame[], rvec rframe[]);

extern real calc_mammoth(int iatoms, rvec frame[], rvec rframe[], int *rnum);

extern real calc_esa(int iatoms, rvec frame[], rvec rframe[]);

extern real call_ISDM(int iatoms, rvec cframe[], rvec rframe[], const char *ISDM);