/* 
 * 
 * Tim Connolly - tconnolly@ucmerced.edu
 * Copyright (c) 2014, Regents of the University of California
 * Released under BSD 2-Clause License (see "LICENSE" file)
 * 
 * This code was modified from the file src/tools/gmx_gyrate.c
 */


#include <math.h>
#include <string.h>

#include "libisdm.h"
#include "eigensolver.h"

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



void mat_mult_mat(real* mat1, real* mat2, int m, int n, int o, real* out, gmx_bool bMP)
{
  /* Assume a is an array of doubles m by n in dimensions.
   * Assume b is an array of doubles n by o in dimensions.
   * Out should point to enough memory to store m by o doubles.
   */
  int i;
  #pragma omp parallel for schedule(dynamic) if (bMP)
  for (i = 0; i < m; i++) {
    int j, k;
    for (k = 0; k < o; k++) {
      out[(i * m) + k] = 0;
      for (j = 0; j < n; j++) {
        out[(i * m) + k] += mat1[(i * m) + j] * mat2[(j * n) + k];
      }
    }
  }
}



void scl_mult_mat(real scl, real* mat, int m, int n, real* out)
{
  /* Multiply all elements of matrix mat by scalar scl.
   * 
   * This should still work even if mat and out point to the same thing.
   */
  int i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      out[(i * m) + j] = scl * mat[(i * m) + j];
    }
  }
}



void mat_transpose(real* mat, int m,int n, real* out)
{
  /* Transpose m by n array of reals a into n by m array of reals out.
   */
  int i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      out[(j * n) + i] = mat[(i * m) + j];
    }
  }
}



void calc_EISD(real** MDS, int nframes, int d, real** EISD)
{
  /* Calculate the approximate ISD from the euclidean dimensionally reduced 
   * coordinates in MDS. Uses the number of dimensions specified by d. The 
   * number of structures is specified by nframes, and EISD should be a 
   * matrix of nframes x nframes.
   */
  int  i, j, k;
  real kEISD;
  for (i = 0; i < nframes; i++)
  {
    for (j = 0; j < nframes; j++)
    {
      // Same structure.
      if (i == j)
      {
        EISD[i][j] = 0.0;
        continue;
      }
      
      // Different structures.
      EISD[i][j] = 0.0;
      for (k = 0; k < d; k++)
      {
        kEISD        = MDS[i][k] - MDS[j][k];
        EISD[i][j]  += kEISD * kEISD;
      }
      EISD[i][j] = sqrt(EISD[i][j]);
    }
  }
}



real calc_rcc(real** ISD, real** EISD, int nframes)
{
  int i, j;
  int N = nframes * (nframes - 1) / 2;
  double sCov, sISD, sEISD, sISD2, sEISD2, mISD, mEISD, vISD, vEISD;
  sCov = 0.0; sISD = 0.0; sEISD = 0.0; sISD2 = 0.0; sEISD2 = 0.0;
  
  // Variances, means, and the means of squares.
  for (i = 0; i < (nframes - 1); i++)
  {
    for (j = (i + 1); j < nframes; j++)
    {
      sISD   += ISD[i][j];
      sEISD  += EISD[i][j];
      sISD2  += ISD[i][j]  * ISD[i][j];
      sEISD2 += EISD[i][j] * EISD[i][j];
    }
  }
  mISD  = sISD  / N;
  mEISD = sEISD / N;
  vISD  = (sISD2  / N) - (mISD  * mISD);
  vEISD = (sEISD2 / N) - (mEISD * mEISD);
  
  // Covariance.
  for (i = 0; i < (nframes - 1); i++)
  {
    for (j = (i + 1); j < nframes; j++)
    {
      sCov += (ISD[i][j] - mISD) * (EISD[i][j] - mEISD);
    }
  }
  
  // Correlation coefficient, R.
  return (sCov / N) / (sqrt(vISD) * sqrt(vEISD));
}



int gmx_isdcmds(int argc,char *argv[])
{
  const char *desc[] = {
    "[TT]g_isdcmds[tt] implements classical multi-dimensional scaling ",
    "by first calculating the matrix of inter-structure distances (ISD). ",
    "The default ISDM if one is not, chosen by the user is RMSD. Only ",
    "one ISDM can be chosen at a time. The -xpm option is required. An ",
    "upper threshold for the ISD can be specified with the setmax option."
  };
  
  
  
  static gmx_bool bANG=FALSE, bDIH=FALSE, bANGDIH=FALSE, bDRMS=FALSE;
  static gmx_bool bPHIPSI=FALSE, bSRMS=FALSE, bPCOR=FALSE, bMAMMOTH=FALSE;
  static gmx_bool bACOR=FALSE, bESA=FALSE, bRMSD=FALSE, bMIR=FALSE;
  static gmx_bool bRG=FALSE, bSRG=FALSE, bE2E=FALSE, bSE2E=FALSE;
  static gmx_bool bANG2=FALSE, bDIH2=FALSE, bANGDIH2=FALSE, bANGDIH2G=FALSE;
  static gmx_bool bRROT=FALSE, bSDRMS=FALSE, bPHIPSI2=FALSE;
  static int nt          = -1;
  static real setmax     = -1.0;
  static real rcutoff    =  1.1;
  static real noisefloor =  0.0;
  static gmx_bool bNoise = FALSE;
  static gmx_bool bMP    = FALSE;
  t_pargs pa[] = {
    { "-ang", FALSE, etBOOL, {&bANG},
    "ISDM: Mean cosine of difference of backbone angles for each "
    "set of three atoms. Assumes only CA atoms." },
    { "-dih", FALSE, etBOOL, {&bDIH},
    "ISDM: Mean cosine of difference of backbone dihedrals for "
    "each set of four atoms. Assumes only CA atoms." },
    { "-angdih", FALSE, etBOOL, {&bANGDIH},
    "ISDM: Geometric mean of ang and dih ISDMs." },
    { "-ang2", FALSE, etBOOL, {&bANG2},
    "ISDM: Attempts to euclideanize -ang." },
    { "-dih2", FALSE, etBOOL, {&bDIH2},
    "ISDM: Attempts to euclideanize -dih." },
    { "-angdih2", FALSE, etBOOL, {&bANGDIH2},
    "ISDM: Attempts to euclideanize -angdih." },
    { "-angdih2g", FALSE, etBOOL, {&bANGDIH2G},
    "ISDM: Attempts to euclideanize -angdih. Geometric mean." },
    { "-phipsi", FALSE, etBOOL, {&bPHIPSI},
    "ISDM: Mean cosine of difference of phi and psi angles. "
    "Assumes only backbone atoms." },
    { "-phipsi2", FALSE, etBOOL, {&bPHIPSI2},
    "ISDM: Attempts to euclideanize -phipsi." },
    { "-drms", FALSE, etBOOL, {&bDRMS},
    "ISDM: Mean difference of the paired distances matrix for all "
    "atoms. Distance RMS(D)." },
    { "-sdrms", FALSE, etBOOL, {&bSDRMS},
    "ISDM: Mean difference of the paired distances matrix for all "
    "atoms scaled by 2 * geometric mean of Rg. Scaled distance "
    "RMS(D)." },
    { "-rg", FALSE, etBOOL, {&bRG},
    "ISDM: Calculates difference in Rg. Only compares size. " },
    { "-srg", FALSE, etBOOL, {&bSRG},
    "ISDM: Calculates difference in Rg scaled by mean Rg. " },
    { "-e2e", FALSE, etBOOL, {&bE2E},
    "ISDM: Calculates difference in end-to-end distance. " },
    { "-se2e", FALSE, etBOOL, {&bSE2E},
    "ISDM: Calculates difference in end-to-end distance scaled "
    "by (2 * Rg). " },
    { "-mir", FALSE, etBOOL, {&bMIR},
    "ISDM: RMSD with the mirror of the reference structure. " },
    { "-rrot", FALSE, etBOOL, {&bRROT},
    "ISDM: RMSD with random rotation of reference structure. " },
    { "-srms", FALSE, etBOOL, {&bSRMS},
    "ISDM: Scaled RMSD. RMSD between the structure and reference "
    "divided by the RMSD between the structure and mirror of the "
    "reference created by multiplying the coordinates by the "
    "negative identity matrix." },
    { "-rmsd", FALSE, etBOOL, {&bRMSD},
    "ISDM: Standard RMSD." },
    { "-pcor", FALSE, etBOOL, {&bPCOR},
    "ISDM: Position correlation. Correlation coefficient of the "
    "positions is computed after alignment. Only positive "
    "correlation is considered. Negative correlations are set to "
    "zero." },
    { "-acor", FALSE, etBOOL, {&bACOR},
    "ISDM: Angle correlation. Correlation coefficient of the "
    "backbone angles (see ang ISDM) is computed. "
    "Only positive correlation is considered. Negative correlations "
    "are set to zero." },
    { "-mammoth", FALSE, etBOOL, {&bMAMMOTH},
    "ISDM: MAMMOTH (MAtching Molecular Models Obtained from "
    "Theory). Compares segments of residues chosen by sequence "
    "alignment. Attempts to focus on correct secondary structure "
    "moreso than tertiary structure. Source code modified for "
    "compatibility. For this ISDM, please cite: \n\n"
    "Ortiz, AR, Strauss, CE, Olmea, O (2002). MAMMOTH "
    "(Matching molecular models obtained from theory): An automated "
    "method for model comparison. Protein Sci. 11 (11), 2606–2621.\n"},
    { "-esa", FALSE, etBOOL, {&bESA},
    "ISDM: Elastic shape analysis. Based on image analysis. "
    "Warps structure onto the reference structure. Original source "
    "code ported from Matlab to C. For this ISDM, please cite: \n\n"
    "Liu W, Srivastava A, Zhang J (2011) A Mathematical Framework "
    "for Protein Structure Comparison. PLoS Comput Biol 7(2): "
    "e1001075.\n\nAssume only CA atoms." },
    { "-mp", FALSE, etBOOL, {&bMP},
      "Use OpenMP commands for parallel processing. "},
    { "-nt", FALSE, etINT, {&nt},
      "Limit the maximum number of threads for parallel processing. "},
    { "-noise", FALSE, etBOOL, {&bNoise},
    "If this flag is set, additional information is sent to "
    "stdout. The tool calculates the number of positive eigenvalues "
    "and the number of positive eigenvalues that can be accounted "
    "by two sources of noise. (1) Algorithmic noise based on the "
    "negative eigenvalues, (2) thermal noise based on the expected "
    "variation of folded proteins, and (3) the combined noise. "
    "An estimate of thermal noise can be set manually with the "
    "option -noisefloor." },
    { "-setmax", FALSE, etREAL, {&setmax},
    "Set maximum value to threshold the xpm file. Must be greater "
    "than the average inter-structure distance." },
    { "-rcutoff", FALSE, etREAL, {&rcutoff},
    "Set cutoff value for the correlation coefficient. Only applies "
    "if the -rcc output is set. The correlation coefficient (R) "
    "will be calculated for each dimensional until rcutoff is "
    "reached. The value should be between 0 and 1." },
    { "-noisefloor", FALSE, etREAL, {&noisefloor},
    "Only applies if the -noise option is set. Manually sets the "
    "the estimate of thermal noise used by the dimensionality "
    "estimator." },
  };
  
  
  
  FILE       *out;
  t_trxstatus *status;
  t_topology top;
  int        ePBC;
  rvec       *x, **frames;
  real       *nweights, *iweights, abscoor, maxcoor;
  real       *diff, **ISDmat, *P2, *J, *P2J, *B, *BT, *E, *V, *MDSa;
  real       **Va, **MDS, **EISD, *EISDm, Rcc, sumne, cumpe;
  double     *avgdiff, *maxdiff, avgISD, maxISD;
  matrix     box;
  real       t, xpm_max, pi = 3.14159265358979;
  int        *maxframe, *rnum, maxcoori;
  int        i, j, k, m, n, p, np, d, iatoms, natoms, nframes, nframes2;
  int        percent_calcs, finished_calcs, noptions;
  gmx_bool   bDFLT, bFit, bISD, bMDS, bEig, bVec, bRcc, bMRg, bDRg, bPy, bM;
  char       buf[256];
  char       *ISDM, *grpname, title[256], title2[256], *rname;
  atom_id    *index;
  output_env_t oenv;
  gmx_rmpbc_t  gpbc=NULL;
  const char *leg[]  = { "D" }; 
  #define NLEG asize(leg) 
  t_filenm fnm[] = {
    { efTRX, "-f",   NULL,       ffREAD }, 
    { efTPS, NULL,   NULL,       ffREAD },
    { efNDX, NULL,   NULL,       ffOPTRD },
    { efXVG, "-eig", "eigvals",  ffOPTWR },
    { efXVG, "-rcc", "corrcoef", ffOPTWR },
    { efXVG, "-mrg", "mrgcorr",  ffOPTWR },
    { efXVG, "-drg", "drgcorr",  ffOPTWR },
    { efDAT, "-vec", "eigvecs",  ffOPTWR },
    { efDAT, "-isd", "isdcsv",   ffOPTWR },
    { efDAT, "-mds", "mdscsv",   ffOPTWR },
    { efDAT, "-py",  "mayapy",   ffOPTWR },
    { efDAT, "-m",   "disp6D",   ffOPTWR },
  }; 
  #define NFILE asize(fnm)
  int npargs;
  
  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  
  // Lots of black magic with this one. The oenv is used by many things.
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
                    NFILE,fnm,npargs,pa,asize(desc),desc,0,NULL,&oenv);
  
  // If there are no options at command line, do default behavior.
  bDFLT = !(bANG || bDIH || bANGDIH || bPHIPSI || bDRMS || bSRMS || bRMSD || 
  bPCOR || bACOR || bMAMMOTH || bESA || bRG || bSRG || bE2E || 
  bSE2E || bMIR || bRROT || bSDRMS || bANG2 || bDIH2 || 
  bANGDIH2 || bPHIPSI2 || bANGDIH2G);
  
  bFit  =  (bDFLT || bRMSD || bMIR || bSRMS || bPCOR);
  
#ifdef _OPENMP
  if (nt > 0)
  {
    omp_set_num_threads(nt);
  }
#endif
  
  /* Reads the tpr file. Outputs a ton of info.
   * 
   * I think this is the line that forces you to have a -s at prompt.
   */
  read_tps_conf(ftp2fn(efTPS, NFILE, fnm), title, &top, &ePBC, &x, NULL, box, TRUE);
  
  // Asks you to choose a selection of atoms at prompt.
  get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &iatoms, &index, &grpname);
  
  // For error checking.
  noptions = 0;
  // Check which ISDM will be used. Default is RMSD.
  if (bDFLT || bRMSD)
  {
    fprintf(stderr,"\nUsing RMSD as ISDM.\n");
    ISDM = "RMSD";
    noptions++;
  }
  
  if (bANG)
  {
    fprintf(stderr,"\nUsing backbone angles as ISDM.\n");
    ISDM = "ANG";
    noptions++;
  }
  
  if (bDIH)
  {
    fprintf(stderr,"\nUsing backbone dihedrals as ISDM.\n");
    ISDM = "DIH";
    noptions++;
  }
  
  if (bANG2)
  {
    fprintf(stderr,"\nUsing backbone angles as ISDM.\n");
    ISDM = "ANG2";
    noptions++;
  }
  
  if (bDIH2)
  {
    fprintf(stderr,"\nUsing backbone dihedrals as ISDM.\n");
    ISDM = "DIH2";
    noptions++;
  }
  
  if (bANGDIH)
  {
    fprintf(stderr,"\nUsing geometric mean of angles and dihedrals as ISDM.\n");
    ISDM = "ANGDIH";
    noptions++;
  }
  
  if (bANGDIH2)
  {
    fprintf(stderr,"\nUsing geometric mean of angles and dihedrals as ISDM.\n");
    ISDM = "ANGDIH2";
    noptions++;
  }
  
  if (bANGDIH2G)
  {
    fprintf(stderr,"\nUsing geometric mean of angles and dihedrals as ISDM.\n");
    ISDM = "ANGDIH2G";
    noptions++;
  }
  
  if (bPHIPSI)
  {
    fprintf(stderr,"\nUsing phi and psi angles as ISDM.\n");
    ISDM = "PHIPSI";
    noptions++;
  }
  
  if (bPHIPSI2)
  {
    fprintf(stderr,"\nUsing phi and psi angles as ISDM.\n");
    ISDM = "PHIPSI2";
    noptions++;
  }
  
  if (bDRMS)
  {
    fprintf(stderr,"\nUsing distance RMS as ISDM.\n");
    ISDM = "DRMS";
    noptions++;
  }
  
  if (bSDRMS)
  {
    fprintf(stderr,"\nUsing scaled distance RMS as ISDM.\n");
    ISDM = "SDRMS";
    noptions++;
  }
  
  if (bRG)
  {
    fprintf(stderr,"\nUsing Rg difference as ISDM.\n");
    ISDM = "RG";
    noptions++;
  }
  
  if (bSRG)
  {
    fprintf(stderr,"\nUsing scaled Rg difference as ISDM.\n");
    ISDM = "SRG";
    noptions++;
  }
  
  if (bE2E)
  {
    fprintf(stderr,"\nUsing end-to-end distance as ISDM.\n");
    ISDM = "E2E";
    noptions++;
  }
  
  if (bSE2E)
  {
    fprintf(stderr,"\nUsing scaled end-to-end distance as ISDM.\n");
    ISDM = "SE2E";
    noptions++;
  }
  
  if (bMIR)
  {
    fprintf(stderr,"\nUsing mirrored RMSD as ISDM.\n");
    ISDM = "MIR";
    noptions++;
  }
  
  if (bSRMS)
  {
    fprintf(stderr,"\nUsing scaled RMSD as ISDM.\n");
    ISDM = "SRMS";
    noptions++;
  }
  
  if (bPCOR)
  {
    fprintf(stderr,"\nUsing position correlation as ISDM.\n");
    ISDM = "PCOR";
    noptions++;
  }
  
  if (bACOR)
  {
    fprintf(stderr,"\nUsing backbone angle correlation as ISDM.\n");
    ISDM = "ACOR";
    noptions++;
  }
  
  if (bRROT)
  {
    fprintf(stderr,"\nUsing RMSD with random rotation as ISDM.\n");
    noptions++;
    
    // Additional stuff for option.
    srand(time(NULL));
  }
  
  if (bMAMMOTH)
  {
    fprintf(stderr,"\nUsing MAMMOTH comparison as ISDM.\n");
    noptions++;
    
    // Additional stuff for option.
    snew(rnum,iatoms);
    //printf(stderr,"\nOutput sequence (tool).\n\n");
    for (i = 0; i < iatoms; i++)
    {
      rname = *(top.atoms.resinfo[top.atoms.atom[index[i]].resind].name);
      // Convert to integers.
      if (!(strcmp(rname, "ALA")))
      {
        rnum[i] = 0;
      }
      else if (!(strcmp(rname, "CYS")))
      {
        rnum[i] = 1;
      }
      else if (!(strcmp(rname, "ASP")))
      {
        rnum[i] = 2;
      }
      else if (!(strcmp(rname, "GLU")))
      {
        rnum[i] = 3;
      }
      else if (!(strcmp(rname, "PHE")))
      {
        rnum[i] = 4;
      }
      else if (!(strcmp(rname, "GLY")))
      {
        rnum[i] = 5;
      }
      else if (!(strcmp(rname, "HIS")) || !(strcmp(rname, "HID")) || 
        !(strcmp(rname, "HIE")) || !(strcmp(rname, "HIP")) || 
        !(strcmp(rname, "HSD")) || !(strcmp(rname, "HSE")) || 
        !(strcmp(rname, "HSP")))
      {
        rnum[i] = 6;
      }
      else if (!(strcmp(rname, "ILE")))
      {
        rnum[i] = 7;
      }
      else if (!(strcmp(rname, "LYS")))
      {
        rnum[i] = 8;
      }
      else if (!(strcmp(rname, "LEU")))
      {
        rnum[i] = 9;
      }
      else if (!(strcmp(rname, "MET")))
      {
        rnum[i] = 10;
      }
      else if (!(strcmp(rname, "ASN")))
      {
        rnum[i] = 11;
      }
      else if (!(strcmp(rname, "PRO")))
      {
        rnum[i] = 12;
      }
      else if (!(strcmp(rname, "GLN")))
      {
        rnum[i] = 13;
      }
      else if (!(strcmp(rname, "ARG")))
      {
        rnum[i] = 14;
      }
      else if (!(strcmp(rname, "SER")))
      {
        rnum[i] = 15;
      }
      else if (!(strcmp(rname, "THR")))
      {
        rnum[i] = 16;
      }
      else if (!(strcmp(rname, "VAL")))
      {
        rnum[i] = 17;
      }
      else if (!(strcmp(rname, "TRP")))
      {
        rnum[i] = 18;
      }
      else if (!(strcmp(rname, "TYR")))
      {
        rnum[i] = 19;
      }
      else
      {
        rnum[i] = 20;
      }
    }
  }
  
  if (bESA)
  {
    fprintf(stderr,"\nUsing ESA comparison as ISDM.\n"
    "For this ISDM, please cite: \n\n"
    "Liu W, Srivastava A, Zhang J (2011) A Mathematical Framework "
    "for Protein Structure Comparison. PLoS Comput Biol 7(2): "
    "e1001075.\n" );
    noptions++;
  }
  
  // Throw an error if multiple -ISDM options were given by the user.
  if (noptions > 1)
  {
    gmx_fatal(FARGS,"\nThis tool only supports using one optional ISDM at a time.\n");
  }
  
  // Check for error on -setmax before doing the calculations.
  if (setmax != -1.0)
  {
    if (setmax <= 0.0)
    {
      gmx_fatal(FARGS,"\nThe argument for -setmax must be greater than 0.\n");
    }
  }
  
  
  /* Opens trj. Reads first frame. Returns status. Allocates mem for x.
   * 
   * Not sure which argument determines which atoms to pull info for.
   */
  printf("\nCounting the number of frames.\n");
  natoms=read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
  
  // Now that we have iatoms, allocate memory for other arrays.
  snew(nweights, natoms);
  snew(iweights, iatoms);
  snew(diff, iatoms);
  
  // Initialize nweights to zeros.
  for (i=0; i < natoms; i++)
  {
    nweights[i] = 0;
  }
  // Makes an array of weights. Necessary for reset_x.
  for (i=0; i < iatoms; i++)
  {
    // Give a value for the weights.
    nweights[(int)index[i]] = 1;
    iweights[i] = 1;
    // While we're at it, initialize diff to zeros.
    diff[i] = 0;
  }
  
  // Output which files?
  bEig = opt2bSet("-eig", NFILE, fnm);
  bRcc = opt2bSet("-rcc", NFILE, fnm);
  bMRg = opt2bSet("-mrg", NFILE, fnm);
  bDRg = opt2bSet("-drg", NFILE, fnm);
  bVec = opt2bSet("-vec", NFILE, fnm);
  bISD = opt2bSet("-isd", NFILE, fnm);
  bMDS = opt2bSet("-mds", NFILE, fnm);
  bPy  = opt2bSet("-py",  NFILE, fnm);
  bM   = opt2bSet("-m",   NFILE, fnm);
  
  nframes = 0;
  do
  {
    /* This loop doesn't do anything.
     * 
     * It's just the most reliable way to find the number of frames.
     */
    nframes++;
  } while(read_next_x(oenv, status, &t, natoms, x, box));
  // Close the trajectory.
  close_trj(status);
  // Throw an error if there aren't enough frames.
  if (nframes < 2)
  {
    gmx_fatal(FARGS, "\nThe trajectory must have at least 2 frames.\n");
  }
  
  
  
  // Create an array to hold all frames.
  snew(frames,  nframes);
  // Create arrays based on nframes.
  snew(maxdiff, nframes);
  snew(avgdiff, nframes);
  snew(ISDmat,  nframes);
  for (i = 0; i < nframes; i++)
  {
    avgdiff[i] = 0;
    snew(ISDmat[i], nframes);
  }
  nframes2 = nframes * nframes;
  
  /* Opens trj. Reads first frame. Returns status. Allocates mem for x.
   * 
   * Not sure which argument determines which atoms to pull info for.
   */
  printf("\nStoring trajectory to memory.\n");
  natoms=read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
  
  // Initialize index to keep track of current frame.
  i = 0;
  // This is for removing periodic boundary conditions.
  gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms,box);
  do
  {
    // Set aside new memory to store this frame.
    snew(frames[i], iatoms);
    // Removes periodic boundary conditions from x.
    gmx_rmpbc(gpbc, natoms, box, x);
    // Centers x. The NULL arguments are necessary to fit based on subset.
    reset_x(natoms, NULL, natoms, NULL, x, nweights);
    // Saves the current frame into frames.
    for (j=0; j<iatoms; j++)
    {
      copy_rvec(x[(int)index[j]], frames[i][j]);
    }
    // Increment frame index.
    i++;
  } while(read_next_x(oenv, status, &t, natoms, x, box));
  // Close the trajectory.
  close_trj(status);
  // Closes the thing that removes periodic boundary conditions.
  gmx_rmpbc_done(gpbc);
  
  // Initialize to 0.
  maxISD = 0.0;
  avgISD = 0.0;
  
  
  
  /* Main calculation loop.
   */
  printf("\nCalculating inter-structure distances. \n");
  
  /* Originally this was designed to only loop through each pair of i and j
   * one time to save half of the calculations. Eventually it became 
   * impractical to make sure that each ISDM was symmetrical, so now the
   * algorithm takes the performance hit in favor of accuracy and simplicity.
   */
  percent_calcs  = 1;
  finished_calcs = 0;
  
  // Loop through reference frames.
  #pragma omp parallel for schedule(dynamic) if (bMP)
  for (i = 0; i < nframes; i++)
  {
    // Some memory required by each thread.
    real   ISD;
    double dISD;
    matrix rrot, rrotx, rroty, rrotz;
    rvec *iframe, *jframe, *cframe, *rframe, rrot_xyz, xold;
    if (bRROT)
    {
      snew(iframe,iatoms);
      // Use up the first few random numbers that usually aren't random.
      rrot_xyz[0] = (real)rand();
      rrot_xyz[1] = (real)rand();
      rrot_xyz[2] = (real)rand();
    }
    if (bFit)
    {
      snew(jframe,iatoms);
    }
    
    // Loop through fitting frames.
    for (j = 0; j < nframes; j++)
    {
      /* In this section, we'll put calls to all of the ISDMs.
       * 
       * Each should have its own if statement, so it is only executed
       * if that option is specified at the command line.
       * 
       * This function doesn't use the output stored in diff.
       */
      
      // Skip for i == j (comparing structure with self).
      if (i == j)
      {
        ISDmat[i][j] = 0;
        continue;
      }
      
      // Copy ith frame.
      if (bRROT)
      {
        // Make a copy of the ith frame.
        copy_rvecn(frames[i], iframe, 0, iatoms);
        rframe = iframe;
      }
      else
      {
        rframe = frames[i];
      }
      
      // Fit the jth frame.
      if (bFit)
      {
        // Need to make a copy of the fit frame or bad stuff will happen.
        copy_rvecn(frames[j], jframe, 0, iatoms);
        // Aligns jframe to current reference frame.
        do_fit(iatoms, iweights, frames[i], jframe);
        cframe = jframe;
      }
      else
      {
        cframe = frames[j];
      }
      
      // Calls most ISDM options.
      if (bDFLT || bRMSD || bSRMS || bRG || bSRG || bE2E || bSE2E || 
        bMIR || bANG || bDIH || bANGDIH || bPHIPSI || bDRMS || 
        bSDRMS || bPCOR || bACOR || bANG2 || bDIH2 || bANGDIH2 || 
        bPHIPSI2 || bANGDIH2G)
      {
        ISD = call_ISDM(iatoms, cframe, rframe, ISDM);
      }
      
      // RMSD with random rotation. User gives -rrot option.
      if (bRROT)
      {
        // Solve for three random numbers.
        for (k = 0; k < 3; k++)
        {
          rrot_xyz[k] = 2.0 * pi * ((real)rand() / RAND_MAX) - pi;
        }
        // Create x, y, z rotation matrices and multiply.
        clear_mat(rrotx);
        clear_mat(rroty);
        clear_mat(rrotz);
        /*      Rx = rrotx[rows][cols]
         * 
         *      |   1.0   |   0.0   |   0.0   |
         * Rx = |   0.0   |  cos(x) |  sin(x) |
         *      |   0.0   | -sin(x) |  cos(x) |
         */
        rrotx[0][0] = 1.0;
        rrotx[1][1] = cos(rrot_xyz[0]);
        rrotx[2][2] = rrotx[1][1];
        rrotx[1][2] = sin(rrot_xyz[0]);
        rrotx[2][1] = -1.0 * rrotx[1][2];
        /*      Ry = rroty[rows][cols]
         * 
         *      |  cos(x) |   0.0   | -sin(x) |
         * Ry = |   0.0   |   1.0   |   0.0   |
         *      |  sin(x) |   0.0   |  cos(x) |
         */
        rroty[1][1] = 1.0;
        rroty[0][0] = cos(rrot_xyz[1]);
        rroty[2][2] = rroty[0][0];
        rroty[2][0] = sin(rrot_xyz[1]);
        rroty[0][2] = -1.0 * rroty[2][0];
        /*      Rz = rrotz[rows][cols]
         * 
         *      |  cos(x) |  sin(x) |   0.0   |
         * Rz = | -sin(x) |  cos(x) |   0.0   |
         *      |   0.0   |   0.0   |   1.0   |
         */
        rrotz[2][2] = 1.0;
        rrotz[0][0] = cos(rrot_xyz[2]);
        rrotz[1][1] = rrotz[0][0];
        rrotz[0][1] = sin(rrot_xyz[2]);
        rrotz[1][0] = -1.0 * rrotz[0][1];
        // Multiply rotation matrices.
        mmul(rrotx, rroty, rrot);
        copy_mat(rrot, rrotx);
        mmul(rrotx, rrotz, rrot);
        // Apply random rotation.
        for (k = 0; k < iatoms; k++)
        {
          for (m = 0; m < 3; m++)
          {
            xold[m] = rframe[k][m];
          }
          for (m = 0; m < 3; m++)
          {
            rframe[k][m] = 0;
            for (n = 0; n < 3; n++)
            {
              rframe[k][m] += rrot[m][n] * xold[n];
            }
          }
        }
        // Calculate RMSD after rotation.
        ISD = sqrt(calc_msd(iatoms, cframe, rframe));
      }
      
      // MAMMOTH. User gives -mammoth option.
      if (bMAMMOTH)
      {
        // Calculate MAMMOTH comparison.
        ISD = calc_mammoth(iatoms, cframe, rframe, rnum);
      }
      
      // ESA.
      if (bESA)
      {
        // Calculate ESA comparison.
        ISD = calc_esa(iatoms, cframe, rframe);
      }
      
      // Use doubles instead of reals for the summations.
      dISD = (double)ISD;
      // Add difference to the difference matrix.
      ISDmat[i][j] = ISD;
      // Update the max and avg difference for scaling.
      if (dISD > maxdiff[i])
      {
        maxdiff[i] = dISD;
      }
      avgdiff[i] += dISD;
      
      // Debugging.
      //printf("On the %i th loop. \n",j);
    }
    
    // Average difference for each frame.
    avgdiff[i] /= (nframes - 1);
    
    // Update progress output. OpenMP critical section.
    #pragma omp critical
    {
      finished_calcs += nframes;
      while ((double)(finished_calcs) / nframes2 >= (double)percent_calcs / 100)
      {
        fprintf(stderr, "Approximately %i percent complete. \r", percent_calcs);
        percent_calcs++;
      }
    } // End of OpenMP critical section.
    
    // Free memory used in parallel section.
    if (bRROT)
    {
      sfree(iframe);
    }
    if (bFit)
    {
      sfree(jframe);
    }
  } // End of OpenMP parallel for loop.
  fprintf(stderr, "\n\n\n");
  
  // Find the final average of differences.
  for (i = 0; i < nframes; i++)
  {
    avgISD += avgdiff[i];
  }
  avgISD /= nframes;
  
  
  if (bISD)
  {
    // Opens the output file.
    out = opt2FILE("-isd", NFILE, fnm, "w");
    
    // Write output.
    for (i = 0; i < nframes; i++)
    {
      fprintf(out, "%12.8f", ISDmat[i][0]);
      for (j = 1; j < nframes; j++)
      {
        fprintf(out, ",%12.8f", ISDmat[i][j]);
      }
      fprintf(out, "\n");
    }
    
    // Close the output file.
    ffclose(out);
  }
  
  
  /* Implements Torgerson's classical multi-dimensional scaling (CMDS) 
   * algorithm.
   * 
   * 1) Convert the ISD matrix to the squared proximities matrix (P2).
   * 
   * 2) Perform double centering on P2.
   *        [B = (-1/2) * J * P2 * J, where J = I - (UnitMatrix / nframes)]
   * 
   * 3) Solve for the eigenvalues and eigenvectors.
   * 
   * 4) Keep only the dimensions corresponding to positive eigenvalues. 
   *    The rest are imaginary dimensions.
   *        [The requirement to keep the dimension here is that both the 
   *         eigenvalue and the root of the eigenvalue must be greater than 
   *         zero to rule out dimensions below the precision limit.]
   * 
   * 5) Sign convention. This may not be necessary.
   */
  fprintf(stderr, "Performing MDS.\n");
  // Allocate memory.
  snew(J,   nframes * nframes);
  snew(P2,  nframes * nframes);
  snew(P2J, nframes * nframes);
  snew(V,   nframes * nframes);
  snew(E,   nframes);
  snew(Va,  nframes);
  for (i = 0; i < nframes; i++)
  {
    Va[i] = &V[i * nframes];
  }
  
  // Step 1.
  fprintf(stderr, "MDS step 1 of 5. \r");
  for (i = 0; i < nframes; i++)
  {
    for (j = 0; j < nframes; j++)
    {
      P2[(i * nframes) + j] = ISDmat[i][j] * ISDmat[i][j];
    }
  }
  
  // Step 2.
  fprintf(stderr, "MDS step 2 of 5. \r");
  // Constructs J.
  for (i = 0; i < nframes; i++)
  {
    for (j = 0; j < nframes; j++)
    {
      if (i == j)
      {
        J[(i * nframes) + j] = 1.0 - (1.0 / nframes);
      }
      else
      {
        J[(i * nframes) + j] = -1.0 / nframes;
      }
    }
  }
  // Solve for B.
  mat_mult_mat(P2,   J, nframes, nframes, nframes, P2J, bMP);
  B = P2; // Finished with the memory in P2. Reuse it to store B.
  scl_mult_mat(-0.5, J, nframes, nframes, J);
  mat_mult_mat(J,  P2J, nframes, nframes, nframes, B, bMP);
  
  // Step 3.
  fprintf(stderr, "MDS step 3 of 5. \r");
  // Fix assymetry in B caused by precision limits.
  BT = J; // Finished with the memory in J. Reuse it to store BT.
  mat_transpose(B, nframes, nframes, BT);
  for (i = 0; i < nframes; i++)
  {
    for (j = 0; j < nframes; j++)
    {
      B[(i * nframes) + j] = (B[(i * nframes) + j] + 
      BT[(i * nframes) + j]) / 2.0;
    }
  }
  // Call the eigensolver which uses a lapack backend.
  // E and V are sorted ascending by eigensolver.
  eigensolver(B, nframes, 0, nframes, E, V);
  
  // Step 4.
  fprintf(stderr, "MDS step 4 of 5. \r");
  // Find eigenvalues > 0.0.
  for (i = 0; i < nframes; i++)
  {
    if (E[i] > 0.0)
    {
      if (sqrt(E[i]) > 0.0)
      {
        p  = i;
        np = nframes - p;
        break;
      }
    }
    if (i == (nframes - 1))
    {
      gmx_fatal(FARGS,"\nThere are zero positive eigenvalues.\n");
    }
  }
  // Save coordinates in reduced dimensions.
  snew(MDSa,  nframes * np);
  snew(MDS,   nframes);
  for (i = 0; i < nframes; i++)
  {
    MDS[i] = &MDSa[np * i];
  }
  for (i = 0; i < nframes; i++)
  {
    for (j = 0; j < np; j++)
    {
      MDS[i][j] = sqrt(E[nframes - j - 1]) * Va[nframes - j - 1][i];
    }
  }
  
  // Step 5.
  fprintf(stderr, "MDS step 5 of 5. \r");
  for (j = 0; j < np; j++)
  {
    maxcoor = -1.0;
    for (i = 0; i < nframes; i++)
    {
      abscoor = abs(MDS[i][j]);
      if (abscoor > maxcoor)
      {
        maxcoor  = abscoor;
        maxcoori = i;
      }
    }
    
    if (MDS[maxcoori][j] < 0.0)
    {
      for (i = 0; i < nframes; i++)
      {
        MDS[i][j] *= -1.0;
      }
    }
  }
  
  // Output the eigenvectors.
  if (bVec)
  {
    // Opens the output file.
    out = opt2FILE("-vec", NFILE, fnm, "w");
    
    // Write output.
    for (j = 0; j < nframes; j++)
    {
      fprintf(out, "%15.6e", Va[0][j]);
      for (i = 1; i < nframes; i++)
      {
        fprintf(out, ",%15.6e", Va[i][j]);
      }
      fprintf(out, "\n");
    }
    
    // Close the output file.
    ffclose(out);
  }
  
  // Release memory.
  sfree(J);
  sfree(P2);
  sfree(P2J);
  sfree(V);
  sfree(Va);
  fprintf(stderr, "\nClassical MDS Complete. \n\n");
  
  
  // Output the eigenvalues.
  if (bEig)
  {
    // Opens the output file.
    out = xvgropen(opt2fn("-eig", NFILE, fnm), 
                   "MDS Eigenvalues", 
                   "Dimension", 
                   "Eigenvalue", 
                   oenv);
    
    // Write output.
    for (i = 1; i <= nframes; i++)
    {
      // Print in reversed order.
      fprintf(out, "%-6i %15.8f \n", i, E[nframes - i]);
    }
    
    // Close the output file.
    ffclose(out);
  }
  
  // Estimate the number of dimensions explained by noise.
  if (bNoise)
  {
    printf("\n\n");
    printf("Positive eigenvalues correspond to real dimensions. ");
    printf("Negative eigenvalues correspond to imaginary dimensions.\n\n");
    
    // Sum the positive and negative eigenvalues.
    cumpe = 0.0;
    for (i = p; i < nframes; i++)
    {
      cumpe += E[i];
    }
    sumne = 0.0;
    for (i = 0; i < p; i++)
    {
      sumne += E[i];
    }
    printf("Sum of positive eigenvalues: %12.6f \n", cumpe);
    printf("Sum of negative eigenvalues: %12.6f \n", sumne);
    
    
    // Output the explained noise.
    printf("%-6i eigenvalues are positive.\n", np);
    printf("%-6i eigenvalues are zero or negative.\n", nframes - np);
    cumpe = 0.0;
    for (i = p; i < nframes; i++)
    {
      cumpe += E[i];
      if (cumpe > abs(sumne))
      {
        break;
      }
    }
    printf("%-6i positive eigenvalues can be explained by negative "
    "eigenvalues.\n", i - p);
    cumpe = 0.0;
    for (i = p; i < nframes; i++)
    {
      cumpe += E[i];
      if (cumpe > noisefloor)
      {
        break;
      }
    }
    printf("%-6i positive eigenvalues can be explained by estimated "
    "thermal noise.\n", i - p);
    cumpe = 0.0;
    for (i = p; i < nframes; i++)
    {
      cumpe += E[i];
      if (cumpe > (abs(sumne) + noisefloor))
      {
        break;
      }
    }
    printf("%-6i positive eigenvalues can be explained by estimated "
    "thermal noise and algorithmic noise combined.\n\n", i - p);
  }
  
  // Output dimensionally reduced coordinates.
  if (bMDS)
  {
    // Opens the output file.
    out = opt2FILE("-mds", NFILE, fnm, "w");
    
    // Write output.
    for (i = 0; i < nframes; i++)
    {
      fprintf(out, "%12.8f", MDS[i][0]);
      for (j = 1; j < np; j++)
      {
        fprintf(out, ",%12.8f", MDS[i][j]);
      }
      fprintf(out, "\n");
    }
    
    // Close the output file.
    ffclose(out);
  }
  
  // Allocates memory to store the approximated ISD.
  if (bRcc || bMRg || bDRg || bPy || bM)
  {
    snew(EISD,  nframes);
    snew(EISDm, nframes * nframes);
    for (i = 0; i < nframes; i++)
    {
      EISD[i] = &EISDm[nframes * i];
    }
  }
  
  // Reduced dimensional visualization.
  if (bPy)
  {
    // Calculate accuracy of the displayed results.
    calc_EISD(MDS, nframes, 6, EISD);
    Rcc = calc_rcc(ISDmat, EISD, nframes);
    fprintf(stdout, "The accuracy for 6D MDS is R = %8.4f.\n\n", Rcc);
    
    // Opens the output file.
    out = opt2FILE("-py", NFILE, fnm, "w");
    
    // Python script header (py).
    fprintf(out, "# Plots MDS output in 6 dimensions:\n");
    fprintf(out, "# x, y, z, r, g, b\n\n");
    
    // Import modules (py).
    fprintf(out, "from mayavi import mlab\n");
    fprintf(out, "import numpy as np\n\n");
    
    // Save data to numpy array (py).
    fprintf(out, "# Save data to numpy array.\n");
    fprintf(out, "MDS = np.array([[%8.4f", MDS[0][0]);
    for (j = 1; j < 6; j++)
    {
      fprintf(out, ",%8.4f", MDS[0][j]);
    }
    fprintf(out, "]");
    for (i = 1; i < nframes; i++)
    {
      fprintf(out, ",\n                [%8.4f", MDS[i][0]);
      for (j = 1; j < 6; j++)
      {
        fprintf(out, ",%8.4f", MDS[i][j]);
      }
      fprintf(out, "]");
    }
    fprintf(out, "])\n\n");
    
    // Calculate box center and range, center at zero (py).
    fprintf(out, "# Calculate box center and range.\n");
    fprintf(out, "bctr = np.mean(MDS, 0)\n");
    fprintf(out, "MDS  = np.subtract(MDS, bctr)\n");
    fprintf(out, "bmin = np.min(MDS) #- Rbead\n");
    fprintf(out, "bmax = np.max(MDS) #+ Rbead\n\n");
    
    // Split MDS by dimensions. Recenter and rescale rgb dimensions (py).
    fprintf(out, "# Split MDS by dimensions. Recenter to 0.5.\n");
    fprintf(out, "xyz, rgb = np.hsplit(MDS, 2)\n");
    fprintf(out, "color_sf = 0.8 / (bmax - bmin)\n");
    fprintf(out, "rgb = np.add(np.multiply(rgb, color_sf), 0.5)\n");
    fprintf(out, "s = np.array([0.01])\n");
    fprintf(out, "s = s[0]\n\n");
    
    /*
     *        // Display first coordinate and set up figure (py).
     *        fprintf(out, "# Display first coordinate and set up figure.\n");
     *        fprintf(out, "x = xyz[0, 0]\n");
     *        fprintf(out, "y = xyz[0, 1]\n");
     *        fprintf(out, "z = xyz[0, 2]\n");
     *        fprintf(out, "r = rgb[0, 0]\n");
     *        fprintf(out, "g = rgb[0, 1]\n");
     *        fprintf(out, "b = rgb[0, 2]\n");
     *        fprintf(out, "if r > 1.0:\n    r = 1.0\n");
     *        fprintf(out, "if r < 0.0:\n    r = 0.0\n");
     *        fprintf(out, "if g > 1.0:\n    g = 1.0\n");
     *        fprintf(out, "if g < 0.0:\n    g = 0.0\n");
     *        fprintf(out, "if b > 1.0:\n    b = 1.0\n");
     *        fprintf(out, "if b < 0.0:\n    b = 0.0\n");
     *        fprintf(out, "mlab.points3d(x, y, z, color=(r, g, b), ");
     *        fprintf(out, "extent=[bmin, bmax, bmin, bmax, bmin, bmax])\n\n");
     */
    
    // Display coordinates (py).
    fprintf(out, "# Display coordinates.\n");
    fprintf(out, "for i in range(0, %i):\n", nframes);
    fprintf(out, "    x = xyz[i, 0]\n");
    fprintf(out, "    y = xyz[i, 1]\n");
    fprintf(out, "    z = xyz[i, 2]\n");
    fprintf(out, "    r = rgb[i, 0]\n");
    fprintf(out, "    g = rgb[i, 1]\n");
    fprintf(out, "    b = rgb[i, 2]\n");
    fprintf(out, "    if r > 1.0:\n        r = 1.0\n");
    fprintf(out, "    if r < 0.0:\n        r = 0.0\n");
    fprintf(out, "    if g > 1.0:\n        g = 1.0\n");
    fprintf(out, "    if g < 0.0:\n        g = 0.0\n");
    fprintf(out, "    if b > 1.0:\n        b = 1.0\n");
    fprintf(out, "    if b < 0.0:\n        b = 0.0\n");
    fprintf(out, "    mlab.points3d(x, y, z, s, color=(r, g, b), scale_factor=1)\n\n");
    
    // Close the output file.
    ffclose(out);
  }
  
  // Reduced dimensional visualization.
  if (bM)
  {
    // Calculate accuracy of the displayed results.
    calc_EISD(MDS, nframes, 6, EISD);
    Rcc = calc_rcc(ISDmat, EISD, nframes);
    fprintf(stdout, "The accuracy for 6D MDS is R = %8.4f.\n\n", Rcc);
    
    // Opens the output file.
    out = opt2FILE("-m", NFILE, fnm, "w");
    
    // Octave function and comments.
    fprintf(out, "function disp6D(varargin)\n");
    fprintf(out, "%% function disp6D(varargin)\n");
    fprintf(out, "%%\n");
    fprintf(out, "%% 'PauseTime': Pause between frames (numeric).\n");
    fprintf(out, "%% 'TimeStep' : Time per frame (numeric).\n");
    fprintf(out, "%% 'Radius'   : Sphere size (numeric).\n");
    fprintf(out, "%% 'Res'      : Sphere resolution (numeric).\n");
    fprintf(out, "%% 'Title'    : Figure title (char).\n");
    fprintf(out, "%% 'GIFName'  : Create animated GIF (char).\n");
    fprintf(out, "%% 'GIFStep'  : Frames per image (numeric).\n");
    fprintf(out, "%% 'ShowLine' : Connect spheres (logical).\n");
    fprintf(out, "%%\n%% Defaults   : \n");
    fprintf(out, "%% No pause, 1.0 ps time step, radius auto, sphere \n");
    fprintf(out, "%% resolution 6, no title, no GIF, 1.0 frame GIF step, \n");
    fprintf(out, "%% no line.\n");
    fprintf(out, "%%\n%% Plots MDS output in 6 dimensions:\n");
    fprintf(out, "%% x, y, z, r, g, b\n\n");
    
    // Defaults.
    fprintf(out, "%% Set defaults.\n");
    fprintf(out, "defPauseTime = -1.0;\n");
    fprintf(out, "defTimeStep  = -1.0;\n");
    fprintf(out, "defRadius    = -1.0;\n");
    fprintf(out, "defRes       = -1.0;\n");
    fprintf(out, "defTitle     = '';\n");
    fprintf(out, "defGIFName   = '';\n");
    fprintf(out, "defGIFStep   = -1.0;\n");
    fprintf(out, "defShowLine  = false;\n\n");
    fprintf(out, "%% Initialize parser.\n");
    fprintf(out, "p = inputParser;\n");
    fprintf(out, "addOptional(p, 'PauseTime', defPauseTime, @isnumeric);\n");
    fprintf(out, "addOptional(p, 'TimeStep',  defTimeStep,  @isnumeric);\n");
    fprintf(out, "addOptional(p, 'Radius',    defRadius,    @isnumeric);\n");
    fprintf(out, "addOptional(p, 'Res',       defRes,       @isnumeric);\n");
    fprintf(out, "addOptional(p, 'Title',     defTitle,     @ischar);\n");
    fprintf(out, "addOptional(p, 'GIFName',   defGIFName,   @ischar);\n");
    fprintf(out, "addOptional(p, 'GIFStep',   defGIFStep,   @isnumeric);\n");
    fprintf(out, "addOptional(p, 'ShowLine',  defShowLine);\n\n");
    fprintf(out, "parse(p, varargin{:});\n\n");
    fprintf(out, "fnm = p.Results.GIFName;");
    
    
    // Save data to matrix.
    fprintf(out, "%% Save data to matrix called MDS.\n");
    fprintf(out, "MDS = [%8.4f", MDS[0][0]);
    for (j = 1; j < 6; j++)
    {
      fprintf(out, ", %8.4f", MDS[0][j]);
    }
    for (i = 1; i < nframes; i++)
    {
      fprintf(out, ";\n       %8.4f", MDS[i][0]);
      for (j = 1; j < 6; j++)
      {
        fprintf(out, ", %8.4f", MDS[i][j]);
      }
    }
    fprintf(out, "];\n\n");
    
    // Accuracy of MDS.
    fprintf(out, "%% Print correlation coefficient of MDS and ISD.\n");
    fprintf(out, "fprintf('The accuracy of MDS is: %%8.4f', %8.4f)\n", Rcc);
    
    // Set bead radius. Calculate box center and range.
    fprintf(out, "%% Calculate plot limits.\n");
    fprintf(out, "bsize = max(max(MDS)) - min(min(MDS));\n");
    fprintf(out, "if (p.Results.Radius < 0.0)\n  R = 0.01 * bsize;\n");
    fprintf(out, "else\n  R = p.Results.Radius;\nend\n");
    fprintf(out, "bctr  = mean(MDS);\n");
    fprintf(out, "bmin  = min(min(MDS)) - R;\n");
    fprintf(out, "bmax  = max(max(MDS)) + R;\n\n");
    
    // Split MDS by dimensions. Recenter and rescale rgb dimensions.
    fprintf(out, "%% Split MDS by dimensions. Recenter to 0.5.\n");
    fprintf(out, "x = MDS(:,1);\n");
    fprintf(out, "y = MDS(:,2);\n");
    fprintf(out, "z = MDS(:,3);\n");
    fprintf(out, "r = MDS(:,4);\n");
    fprintf(out, "g = MDS(:,5);\n");
    fprintf(out, "b = MDS(:,6);\n\n");
    fprintf(out, "%% Rescale colors.\n");
    fprintf(out, "color_sf = 0.8 / bsize;\n");
    fprintf(out, "r = (r - bctr(4)) * color_sf + 0.5;\n");
    fprintf(out, "g = (g - bctr(5)) * color_sf + 0.5;\n");
    fprintf(out, "b = (b - bctr(6)) * color_sf + 0.5;\n\n");
    
    // Setup main figure.
    fprintf(out, "%% Setup figure.\n");
    fprintf(out, "n = %6i;\n", nframes);
    fprintf(out, "if (p.Results.Res < 0.0)\n  Res = 6;\nelse\n  ");
    fprintf(out, "Res = p.Results.Res;\nend\n");
    fprintf(out, "[Sx, Sy, Sz] = sphere(Res);\n");
    fprintf(out, "Sx = R * Sx; Sy = R * Sy; Sz = R * Sz;\n");
    fprintf(out, "figure;\n");
    fprintf(out, "axis([bmin, bmax, bmin, bmax, bmin, bmax]);\n\n");
    fprintf(out, "%% The following line can be commented in Matlab.\n");
    fprintf(out, "axis('equal');\n");
    fprintf(out, "%% Uncomment the following line in Matlab.\n");
    fprintf(out, "%%axis('vis3d');\n\n");
    
    // Display coordinates.
    fprintf(out, "%% Display 6D coordinates.\n");
    fprintf(out, "hold on;\n");
    fprintf(out, "for i = 1:n\n  ");
    fprintf(out, "c = [r(i), g(i), b(i)];\n  ");
    fprintf(out, "h = surf(Sx + x(i), Sy + y(i), Sz + z(i));\n  ");
    fprintf(out, "set(h, 'FaceColor', c, 'EdgeColor', 'none');\n\n  ");
    
    // Figure and file options.
    fprintf(out, "if (p.Results.PauseTime > 0.0)\n    ");
    fprintf(out, "pause(p.Results.PauseTime);\n  end\n  ");
    fprintf(out, "if (p.Results.TimeStep > 0.0)\n    ");
    fprintf(out, "iTime  = num2str(i * p.Results.TimeStep);\n    ");
    fprintf(out, "iTime  = strcat(iTime,' ns');\n  ");
    fprintf(out, "else\n    ");
    fprintf(out, "iTime  = '';\n  end\n  ");
    fprintf(out, "iTitle = strcat(p.Results.Title,' ',iTime);\n  ");
    fprintf(out, "title(iTitle);\n\n  ");
    fprintf(out, "if (~strcmp(fnm,''))\n    ");
    fprintf(out, "if (i == 1)\n      ");
    fprintf(out, "f  = getframe(gcf);\n      ");
    fprintf(out, "im = frame2im(f);\n      ");
    fprintf(out, "[imind,cm] = rgb2ind(im,256);\n      ");
    fprintf(out, "imwrite(imind,cm,fnm,'gif','Loopcount',inf);\n      ");
    fprintf(out, "n  = 0;\n      continue\n    end\n    ");
    fprintf(out, "n  = n + 1;\n    ");
    fprintf(out, "if (n >= p.Results.GIFStep)\n      ");
    fprintf(out, "n  = 0;\n      ");
    fprintf(out, "f  = getframe(gcf);\n      ");
    fprintf(out, "im = frame2im(f);\n      ");
    fprintf(out, "[imind,cm] = rgb2ind(im,256);\n      ");
    fprintf(out, "imwrite(imind,cm,fnm,'gif','WriteMode','append');\n    ");
    fprintf(out, "end\n  end\n");
    fprintf(out, "end\n\n");
    
    // Line option.
    fprintf(out, "%% Optionally draw a line to show the time component.\n");
    fprintf(out, "if (p.Results.ShowLine)\n    c = [0.75, 0.75, 0.75];\n    ");
    fprintf(out, "plot3(x, y, z, 'LineWidth', 1, 'Color', c);\n");
    fprintf(out, "end\nhold off;\n\n");
    
    // Close .m script.
    fprintf(out, "end");
    // Close the output file.
    ffclose(out);
  }
  
  // Tests the accuracy of the dimensionally reduced coordinates.
  if (bRcc)
  {
    // Error checking for rcutoff.
    if (rcutoff < 0.0)
    {
      gmx_fatal(FARGS,"\nThe argument for -rcutoff must be greater "
      "than or equal to 0.\n");
    }
    fprintf(stderr, "\nCalculating accuracy of dimensionality reduction."
    "\n");
    
    // Opens the output file.
    out = xvgropen(opt2fn("-rcc", NFILE, fnm), 
                   "Accuracy of Dimensionality Reduction", 
                   "Dimension", 
                   "Correlation Coefficient, R", 
                   oenv);
    
    for (d = 1; d <= np; d++)
    {
      // Calculate correlation coefficient.
      calc_EISD(MDS, nframes, d, EISD);
      Rcc = calc_rcc(ISDmat, EISD, nframes);
      
      // Write to file.
      fprintf(out, "%-6i %12.8f \n", i, Rcc);
      
      if (Rcc > rcutoff)
      {
        break;
      }
    }
    
    // Close file.
    ffclose(out);
    
    printf("\nThe rcutoff is: %12.8f \n", rcutoff);
    printf("The final correlation coefficient is: %12.8f \n", Rcc);
    printf("The estimated dimensionality is: %-6i \n", d);
  }
  
  // Tests correlation between ISD and Rg.
  if (bMRg)
  {
    fprintf(stderr, "\nCalculating correlation of ISD with Rg.\n");
    
    // Opens the output file.
    out = xvgropen(opt2fn("-mrg", NFILE, fnm), 
                   "Correlation of Rg with ISD", 
                   "Radius of Gyration, Rg (nm)", 
                   "ISD", 
                   oenv);
    
    /* Main calculation loop.
     */
    printf("\nCalculating Rg matrix. \n");
    
    percent_calcs = 1;
    
    // Loop through reference frames.
    for (i = 0; i < nframes; i++)
    {
      // Loop through fitting frames.
      for (j = 0; j < nframes; j++)
      {
        // Skip for i == j (comparing structure with self).
        if (i == j)
        {
          EISD[i][j] = 0.0;
          continue;
        }
        
        EISD[i][j] = call_ISDM(iatoms, frames[j], frames[i], "MRG");
        
        // Write to file.
        fprintf(out, "%12.8f %12.8f \n", EISD[i][j], ISDmat[i][j]);
      }
      
      // Update progress output.
      while ((double)(i+1)/nframes >= (double)percent_calcs/100)
      {
        fprintf(stderr, "Approximately %i percent complete. \r", 
                percent_calcs);
        percent_calcs++;
      }
    }
    
    Rcc = calc_rcc(ISDmat, EISD, nframes);
    printf("The Rg vs ISD correlation coefficient is: %12.8f \n", Rcc);
    // Close file.
    ffclose(out);
  }
  
  // Tests correlation between ISD and Rg difference.
  if (bDRg)
  {
    fprintf(stderr, "\nCalculating correlation of ISD with Rg difference."
    "\n");
    
    // Opens the output file.
    out = xvgropen(opt2fn("-drg", NFILE, fnm), 
                   "Correlation of Rg difference with ISD", 
                   "Radius of Gyration Difference (nm)", 
                   "ISD", 
                   oenv);
    
    /* Main calculation loop.
     */
    printf("\nCalculating Rg difference matrix. \n");
    
    percent_calcs = 1;
    
    // Loop through reference frames.
    for (i = 0; i < nframes; i++)
    {
      // Loop through fitting frames.
      for (j = 0; j < nframes; j++)
      {
        // Skip for i == j (comparing structure with self).
        if (i == j)
        {
          EISD[i][j] = 0;
          continue;
        }
        
        EISD[i][j] = call_ISDM(iatoms, frames[j], frames[i], "RG");
        
        // Write to file.
        fprintf(out, "%12.8f %12.8f \n", EISD[i][j], ISDmat[i][j]);
      }
      
      // Update progress output.
      while ((double)(i+1)/nframes >= (double)percent_calcs/100)
      {
        fprintf(stderr, "Approximately %i percent complete. \r", 
                percent_calcs);
        fflush(stderr);
        percent_calcs++;
      }
    }
    
    Rcc = calc_rcc(ISDmat, EISD, nframes);
    printf("The Rg difference vs ISD correlation coefficient is: "
    "%12.8f \n", Rcc);
    // Close file.
    ffclose(out);
  }
  
  
  // Closing.
  thanx(stderr);
  return 0;
}
