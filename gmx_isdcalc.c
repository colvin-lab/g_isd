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



void write_mat_levels(FILE *out,int n_x, int n_y,int *nlevels,real lo,real hi,t_rgb rlo,t_rgb rhi);
void write_mat_levels3(FILE *out,int n_x,int n_y,int *nlevels,real lo,real mid,real hi,t_rgb rlo,t_rgb rmid,t_rgb rhi);
void write_mat_pixels(FILE *out,int n_x,int n_y,real **matrix,real lo,real hi,int nlevels);
void write_mat_pixels3(FILE *out,int n_x,int n_y,real **matrix,real lo,real mid,real hi,int nlevels);



int gmx_isdcalc(int argc,char *argv[])
{
    const char *desc[] = {
        "[TT]g_isdcalc[tt] implements a list of measures designed to ",
        "compare two frames from a trajectory and quantify their difference. ",
        "This tool analyzes the input by comparing each pair of frames in ",
        "the trajectory and recording some statistics. The quantified ",
        "difference between two structures is referred to as the inter-",
        "structure distance (ISD), and the methods used to measure the ",
        "ISD are referred to as ISDMs. ",
        "If no output options are chosen, the mean overall ISD, the overall ",
        "maximum ISD, and the maximally different pair of frames are output ",
        "to stdout. ",
        "The mean and maximum interstructure distance (ISD) for each frame ",
        "can be output with the -mean and -max options. Variance ",
        "calculations are optional, and the -var option outputs variance ",
        "of ISD for each frame to a file and the overall variance to stdout. ",
        "Optionally, a subset of frames can be used as reference structures ",
        "(still compared to all other frames) with the bf and ef ",
        "options. ",
        "The default ISDM if one is not chosen by the user is RMSD, and only ",
        "one ISDM can be chosen at a time. ",
        "For targeted and folding MD simulations, the trajectory can be ",
        "compared against a reference structure (taken from the topology ",
        "file) with the -ref output option. If -ref is the only output ",
        "chosen, the tool will skip the calculations for mean ISD and should ",
        "complete in less time. "
    };
    
    static gmx_bool bANG=FALSE, bDIH=FALSE, bANGDIH=FALSE, bDRMS=FALSE;
    static gmx_bool bPHIPSI=FALSE, bSRMS=FALSE, bPCOR=FALSE, bMAMMOTH=FALSE;
    static gmx_bool bACOR=FALSE, bESA=FALSE, bRMSD=FALSE, bMIR=FALSE;
    static gmx_bool bRG=FALSE, bSRG=FALSE, bE2E=FALSE, bANGDIH2G=FALSE;
    static gmx_bool bANG2=FALSE, bDIH2=FALSE, bANGDIH2=FALSE, bSE2E=FALSE;
    static gmx_bool bRROT=FALSE, bSDRMS=FALSE, bPHIPSI2=FALSE, bGMRG=FALSE;
    static gmx_bool bMRMS=FALSE;
    static int      user_bf = -1, user_ef = -1, user_td = -1, user_tdsens = -1;
    static real     setmax = -1.0;
    t_pargs pa[] = {
        { "-ang", FALSE, etBOOL, {&bANG},
            "ISDM: Mean cosine of difference of backbone angles for each "
            "set of three atoms. Assumes only CA atoms." },
        { "-dih", FALSE, etBOOL, {&bDIH},
            "ISDM: Mean cosine of difference of backbone dihedrals for "
            "each set of four atoms. Assumes only CA atoms." },
        { "-angdih", FALSE, etBOOL, {&bANGDIH},
            "ISDM: Geometric mean of ang and dih measures." },
        { "-ang2", FALSE, etBOOL, {&bANG2},
            "ISDM: Attempts to euclideanize -ang." },
        { "-dih2", FALSE, etBOOL, {&bDIH2},
            "ISDM: Attempts to euclideanize -dih." },
        { "-angdih2", FALSE, etBOOL, {&bANGDIH2},
            "ISDM: Attempts to euclideanize -angdih." },
        { "-angdih2g", FALSE, etBOOL, {&bANGDIH2G},
            "ISDM: Geometric mean based on -angdih2." },
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
        { "-gmrg", FALSE, etBOOL, {&bGMRG},
            "ISDM: Calculates geometric mean of Rg. " },
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
        { "-mrms", FALSE, etBOOL, {&bMRMS},
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
            "backbone angles (see ang measure) is computed. "
            "Only positive correlation is considered. Negative correlations "
            "are set to zero." },
//         { "-mammoth", FALSE, etBOOL, {&bMAMMOTH},
//             "ISDM: MAMMOTH (MAtching Molecular Models Obtained from "
//             "Theory). Compares segments of residues chosen by sequence "
//             "alignment. Attempts to focus on correct secondary structure "
//             "moreso than tertiary structure. Source code modified for "
//             "compatibility. For this measure, please cite: \n\n"
//             "Ortiz, AR, Strauss, CE, Olmea, O (2002). MAMMOTH "
//             "(Matching molecular models obtained from theory): An automated "
//             "method for model comparison. Protein Sci. 11 (11), 2606–2621.\n"},
        { "-esa", FALSE, etBOOL, {&bESA},
            "ISDM: Elastic shape analysis. Based on image analysis. "
            "Warps structure onto the reference structure. Original source "
            "code ported from Matlab to C. Assumes only CA atoms. "
            "For this measure, please cite: \n\n"
            "Liu W, Srivastava A, Zhang J (2011) A Mathematical Framework "
            "for Protein Structure Comparison. PLoS Comput Biol 7(2): "
            "e1001075.\n\n" },
        { "-bf", FALSE, etINT, {&user_bf},
            "Compare range of frames from bf to ef to all other " 
            "frames. The bf and ef options are applied after the b, " 
            "e, and dt options and use units of frames instead of units of " 
            "time. Frame numbers are counted from one." 
        },
        { "-ef", FALSE, etINT, {&user_ef},
            "Compare range of frames from bf to ef to all other " 
            "frames. The bf and ef options are applied after the b, " 
            "e, and dt options and use units of frames instead of units of " 
            "time. Frame numbers are counted from one." 
        },
        { "-td", FALSE, etINT, {&user_td},
            "Number of frames used for the time difference of -tdo output. " },
        { "-tdsens", FALSE, etINT, {&user_tdsens},
            "Number of frames corresponding to the maximum time difference "
            "used to calculate sensitivity of -sens output option." },
        { "-setmax", FALSE, etREAL, {&setmax},
            "Set maximum value to threshold the xpm file. Must be greater "
            "than the average inter-structure distance." },
    };
    
    
    
    FILE       *out, *out1, *out2;
    t_trxstatus *status;
    t_topology top;
    t_rgb      rlo, rmd, rhi;
    rvec       *x, *iframe, *jframe, *rframe, *cframe;
    rvec       **frames, *topframe, *trjframe;
    real       *nweights, *iweights, *xticks, *yticks, xpm_max;
    real       ISD, **ISDMat, minDCR, maxDCR;
    real       t, t1, t2, dt, dm, dn;
    double     dISD, maxISD, avgISD, varISD, msqISD, avgSCL, maxSCL, dSNR;
    double     *maxISDF, *avgISDF, *varISDF, *DCR, *VDCR, *SNR, *sens;
    matrix     box;
    int        *maxframe, maxframei, *rnum, *nsum, ref_bf, ref_ef, oobc;
    int        i, j, k, m, n, ij, mn, bf, ef, iatoms, natoms, nframes, nf2;
    int        ePBC, pcalcs, noptions;
    static int nlevels = 81;
    gmx_bool   bDFLT, bAvg, bVar, bMax, bPair, bRef, bFit, bSens, bSNR, bTD;
    gmx_bool   bMap, bISD, bISDMat, bDCR, bVDCR, bCalcDCR, bMinDCR, bMaxDCR;
    gmx_bool   bAvgSCL, bMaxSCL;
    atom_id    *index;
    output_env_t oenv;
    gmx_rmpbc_t  gpbc=NULL;
    char       *ISDM, *grpname, title[256], title2[256], buf[256], *rname;
    const char *leg[]  = { "D" }; 
#define NLEG asize(leg) 
    t_filenm fnm[] = {
        { efTRX, "-f",      NULL,     ffREAD }, 
        { efTPS, NULL,      NULL,     ffREAD },
        { efNDX, NULL,      NULL,     ffOPTRD },
        { efXVG, "-avg",    "avg",    ffOPTWR },
        { efXVG, "-var",    "var",    ffOPTWR },
        { efXVG, "-max",    "max",    ffOPTWR },
        { efXVG, "-pair",   "pair",   ffOPTWR },
        { efXVG, "-ref",    "ref",    ffOPTWR },
        { efXPM, "-map",    "map",    ffOPTWR }, 
        { efDAT, "-isd",    "isd",    ffOPTWR },
        { efXVG, "-decorr", "decorr", ffOPTWR },
        { efXVG, "-mindcr", "mindcr", ffOPTWR },
        { efXVG, "-maxdcr", "maxdcr", ffOPTWR },
        { efXVG, "-avgscl", "avgscl", ffOPTWR },
        { efXVG, "-maxscl", "maxscl", ffOPTWR },
        { efXVG, "-vdcr",   "vdcr",   ffOPTWR },
        { efXVG, "-snr",    "snr",    ffOPTWR },
        { efXVG, "-tdo",    "tdo",    ffOPTWR },
        { efXVG, "-sens",   "sens",   ffOPTWR },
    }; 
#define NFILE asize(fnm)
    int npargs;
    
    CopyRight(stderr,argv[0]);
    npargs = asize(pa);
    
    // Lots of black magic with this one. The oenv is used by many things.
    parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
                      NFILE,fnm,npargs,pa,asize(desc),desc,0,NULL,&oenv);
    
    // Output which files?
    bAvg    = opt2bSet("-avg",    NFILE, fnm);
    bVar    = opt2bSet("-var",    NFILE, fnm);
    bMax    = opt2bSet("-max",    NFILE, fnm);
    bPair   = opt2bSet("-pair",   NFILE, fnm);
    bRef    = opt2bSet("-ref",    NFILE, fnm);
    bMap    = opt2bSet("-map",    NFILE, fnm);
    bISD    = opt2bSet("-isd",    NFILE, fnm);
    bDCR    = opt2bSet("-decorr", NFILE, fnm);
    bMinDCR = opt2bSet("-mindcr", NFILE, fnm);
    bMaxDCR = opt2bSet("-maxdcr", NFILE, fnm);
    bAvgSCL = opt2bSet("-avgscl", NFILE, fnm);
    bMaxSCL = opt2bSet("-maxscl", NFILE, fnm);
    bVDCR   = opt2bSet("-vdcr",   NFILE, fnm);
    bSNR    = opt2bSet("-snr",    NFILE, fnm);
    bTD     = opt2bSet("-tdo",    NFILE, fnm);
    bSens   = opt2bSet("-sens",   NFILE, fnm);
    
    // Check if decorrelation needs to be calculated.
    if (bDCR || bAvgSCL || bMaxSCL || bVDCR || bSNR || bSens)
    {
        bCalcDCR = TRUE;
    }
    else
    {
        bCalcDCR = FALSE;
    }
    // Check if full ISD matrix needs to be calculated.
    if (bMap || bISD || bCalcDCR || bMinDCR || bMaxDCR)
    {
        bISDMat = TRUE;
    }
    else
    {
        bISDMat = FALSE;
    }
    
    
    /* Reads the tpr file. Outputs a ton of info.
     * 
     * I think this is the line that forces you to have a -s at prompt.
     */
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&x,NULL,box,TRUE);
    
    // Asks you to choose a selection of atoms at prompt.
    get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&iatoms,&index,&grpname);
    
    // Save the reference structure if -ref option chosen.
    if (bRef)
    {
        // Set aside new memory to store this frame.
        snew(topframe, iatoms);
        // Set aside new memory for the trajectory frame.
        snew(trjframe, iatoms);
        // Saves the current frame into frames.
        for (j = 0; j < iatoms; j++)
        {
            copy_rvec(x[(int)index[j]], topframe[j]);
        }
    }
    
    
    /* Since ISDM options were being added and removed frequently, I chose 
     * to keep track of them as booleans rather than an enum list.
     * 
     * 
     * If there are no options at command line, do default behavior.
     */
    bDFLT = !(bANG || bDIH || bANGDIH || bPHIPSI || bDRMS || bSRMS || bRMSD || 
              bPCOR || bACOR || bMAMMOTH || bESA || bRG || bSRG || bE2E || 
              bSE2E || bMIR || bRROT || bSDRMS || bANG2 || bDIH2 || 
              bPHIPSI2 || bANGDIH2 || bANGDIH2G || bGMRG || bMRMS);
    
    bFit  =  (bDFLT || bRMSD || bMIR || bSRMS || bMRMS || bPCOR);
    
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
    
    if (bGMRG)
    {
        fprintf(stderr,"\nCalculating geometric mean of Rg.\n");
        ISDM = "GMRG";
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
        fprintf(stderr,"\nUsing size-independent RMSD as ISDM.\n");
        ISDM = "SRMS";
        noptions++;
    }
    
    if (bMRMS)
    {
        fprintf(stderr,"\nUsing mirror-scaled RMSD as ISDM.\n");
        ISDM = "MRMS";
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
    
    if (bANGDIH)
    {
        fprintf(stderr,"\nUsing geometric mean of angles and dihedrals as ISDM.\n");
        ISDM = "ANGDIH";
        noptions++;
    }
    
    if (bANGDIH2)
    {
        fprintf(stderr,"\nUsing root mean square of angles and dihedrals as ISDM.\n");
        ISDM = "ANGDIH2";
        noptions++;
    }
    
    if (bANGDIH2G)
    {
        fprintf(stderr,"\nUsing geometric mean of RMS angles and RMS dihedrals as ISDM.\n");
        ISDM = "ANGDIH2G";
        noptions++;
    }
    
    if (bRROT)
    {
        fprintf(stderr,"\nUsing RMSD with random rotation as ISDM.\n");
        ISDM = "RROT";
        noptions++;
        
        // Additional stuff for option.
        srand(time(NULL));
        // Use up the first few random numbers that usually aren't random.
        rand(); rand(); rand();
    }
    
    if (bMAMMOTH)
    {
        fprintf(stderr,"\nUsing MAMMOTH comparison as ISDM.\n");
        noptions++;
        
        // Additional stuff for option.
        snew(rnum, iatoms);
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
    
    
    // Opens trj. Reads first frame. Returns status. Allocates mem for x.
    fprintf(stderr,"\nCounting the number of frames.\n");
    natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
    
    // Now that we have iatoms, allocate memory for other arrays.
    if (bRROT)
    {
        snew(iframe,iatoms);
    }
    if (bFit)
    {
        snew(jframe,iatoms);
    }
    snew(nweights,natoms);
    snew(iweights,iatoms);
    
    // Makes an array of weights. Necessary for reset_x.
    for (i=0; i<iatoms; i++)
    {
        // Give a value for the weights.
        nweights[(int)index[i]] = 1.0;
        iweights[i] = 1.0;
    }
    
    
    nframes = 0, t2 = 0;
    do
    {
        /* This loop doesn't do anything.
         * 
         * It's just the most reliable way to find the number of frames.
         */
        t1 = t2;
        t2 = t;
        nframes++;
    } while(read_next_x(oenv,status,&t,natoms,x,box));
    // Close the trajectory.
    close_trj(status);
    // Throw an error if there aren't enough frames.
    if (nframes < 2)
    {
        gmx_fatal(FARGS,"\nThe trajectory must have at least 2 frames.\n");
    }
    // Calculating variance requires more than 2 frames.
    if (bVar)
    {
        if (nframes < 3)
        {
            gmx_fatal(FARGS, "\nCalculating variance requires at least 3 "
                      "frames.\n");
        }
    }
    // Find trajectory time steps. Assumes even spacing. Find nframes / 2.
    dt  = t2 - t1;
    nf2 = nframes / 2;
    
    
    
    // Check to see if the user gives a range of frames.
    if (user_bf == -1)
    {
        bf     = 1;
        ref_bf = 1;
    }
    else
    {
        // Error check for first frame being within range 1 to nframes.
        if ((user_bf < 1) || (user_bf > nframes))
        {
            gmx_fatal(FARGS,"\nArgument to bf must be between 1 and last frame.\n");
        }
        if (bISDMat)
        {
            bf = 1;
        }
        else
        {
            bf = user_bf;
        }
        ref_bf = user_bf;
    }
    if (user_ef == -1)
    {
        ef     = nframes;
        ref_ef = nframes;
    }
    else
    {
        // Error check for last frame being within range 1 to nframes.
        if ((user_ef < 1) || (user_ef > nframes))
        {
            gmx_fatal(FARGS,"\nArgument to ef must be between 1 and last frame.\n");
        }
        if (bISDMat)
        {
            ef = nframes;
        }
        else
        {
            ef = user_ef;
        }
        ref_ef = user_ef;
    }
    // Error check for last frame being greater than first frame.
    if (ref_bf > ref_ef)
    {
        gmx_fatal(FARGS,"\nValue of -ef must be greater than or equal to -bf. "
                        "\n");
    }
    // Check for errors for the -td option before calculations begin.
    if (bTD)
    {
        if (user_td == -1)
        {
            gmx_fatal(FARGS,"\nThe -tdo option requires -td to be set. \n");
        }
        if ((user_td >= nframes) || (user_td < 1))
        {
            gmx_fatal(FARGS,"\nValue of -td must be between 1 and "
                            "(nframes - 1). \n");
        }
    }
    // Check for errors for the -tdsens option before calculations begin.
    if (bSens)
    {
        if (user_tdsens == -1)
        {
            gmx_fatal(FARGS,"\nThe -sens option requires -tdsens to be set. \n");
        }
        if ((user_tdsens > (nf2 / 2)) || (user_tdsens < 1))
        {
            gmx_fatal(FARGS,"\nValue of -tdsens must be between 1 and "
                            "(nframes / 4). \n");
        }
    }
    
    
    
    // For -ref option, only do this loop.
    if (bRef)
    {
        // Load first frame to x variable and initialize PBC removal.
        natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
        gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms, box);
        // Inform user and open output file.
        fprintf(stderr,"\nCalculating ISD from reference structure.\n");
        out = xvgropen(opt2fn("-ref", NFILE, fnm), 
                       "ISD From Reference Structure", 
                       "Time", 
                       "ISD", 
                       oenv);
        
        /* This section seems hacky since it is not based on all atoms.
         * It might be better to close the trajectory, reload the tpr
         * structure to x, and save the preprocessed frame in this section.
         * However, the trajectory would have to be reopened. This could be
         * slower and would leave strange output at the command line.
         * 
         * Remove PBC and centers topframe.
         */
        gmx_rmpbc(gpbc, iatoms, box, topframe);
        reset_x(iatoms, NULL, iatoms, NULL, topframe, iweights);
        // Calculation loop for -ref option.
        do
        {
            // Remove PBC and center x.
            gmx_rmpbc(gpbc, natoms, box, x);
            reset_x(natoms, NULL, natoms, NULL, x, nweights);
            
            // Copy topframe for rrot ISDM.
            if (bRROT)
            {
                // Make a copy of the topology frame.
                copy_rvecn(topframe, iframe, 0, iatoms);
                rframe = iframe;
            }
            else
            {
                rframe = topframe;
            }
            // Save current frame into trjframe.
            for (j = 0; j < iatoms; j++)
            {
                copy_rvec(x[(int)index[j]], trjframe[j]);
            }
            
            // Copy trjframe if a fit is required.
            if (bFit)
            {
                // Need to make a copy of the fit frame or bad stuff will happen.
                copy_rvecn(trjframe, jframe, 0, iatoms);
                // Aligns jframe to current reference frame.
                do_fit(iatoms, iweights, topframe, jframe);
                cframe = jframe;
            }
            else
            {
                cframe = trjframe;
            }
            
            /* In this section, we'll put calls to all of the ISDMs.
             * 
             * This function doesn't use the output stored in diff.
             */
            
            // Calls most ISDM options.
            if (bDFLT || bRMSD || bSRMS || bRG || bSRG || bE2E || bSE2E || 
                bMIR || bANG || bDIH || bANGDIH || bPHIPSI || bDRMS || 
                bSDRMS || bPCOR || bACOR || bANG2 || bDIH2 || bANGDIH2 || 
                bPHIPSI2 || bANGDIH2G || bGMRG || bMRMS || bRROT)
            {
                ISD = call_ISDM(iatoms, cframe, rframe, ISDM);
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
            
            
            // Output result to file.
            fprintf(out,"%10f %12.8f \n", t, ISD);
            
        } while(read_next_x(oenv, status, &t, natoms, x, box));
        // Close trajectory and output files.
        close_trj(status);
        ffclose(out);
        // Closes the thing that removes periodic boundary conditions.
        gmx_rmpbc_done(gpbc);
    }
    
    
    
    // Check to see if any other output options were chosen.
    if (!(bAvg || bVar || bMax || bPair || bTD || bISDMat))
    {
        // No other options chosen, so end the program.
        thanx(stderr);
        return 0;
    }
    
    /* Opens trj. Reads first frame. Returns status. Allocates mem for x.
     */
    fprintf(stderr,"\nStoring trajectory to memory.\n");
    // Load first frame to x variable and initialize PBC removal.
    natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
    gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms, box);
    // Initialize index to keep track of current frame.
    i = 0;
    // Create an array to hold all frames.
    snew(frames,nframes);
    do
    {
        // Set aside new memory to store this frame.
        snew(frames[i], iatoms);
        // Removes periodic boundary conditions from x.
        gmx_rmpbc(gpbc, natoms, box, x);
        // Centers x. The NULL arguments are necessary to fit based on subset.
        reset_x(natoms, NULL, natoms, NULL, x, nweights);
        // Saves the current frame into frames.
        for (j = 0; j < iatoms; j++)
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
    
    
    
    // For -tdo option, only do this loop.
    if (bTD)
    {
        // Inform user and open output file.
        fprintf(stderr, "\nCalculating ISD at time difference set by -td. "
                        "\n\n");
        char hdString[32] = "ISD From Frame + td = ";
        char tdString[8];
        sprintf(tdString, "%6i", user_td);
        out = xvgropen(opt2fn("-tdo", NFILE, fnm), 
                       strcat(hdString, tdString), 
                       "Time", 
                       "ISD", 
                       oenv);
        
        /* Main calculation loop.
         */
        // Percentage of calculations complete.
        pcalcs = 1;
        // Loop through time steps.
        for (i = 0; i < (nframes - user_td); i++)
        {
            j = i + user_td;
            /* In this section, we'll put calls to all of the ISDMs.
             */
            
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
                bPHIPSI2 || bANGDIH2G || bGMRG || bMRMS || bRROT)
            {
                ISD = call_ISDM(iatoms, cframe, rframe, ISDM);
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
            
            // Output result to file.
            fprintf(out,"%10f %12.8f \n", (real)i * dt / 1000.0, ISD);
            
            // Update progress output.
            while ((double)i / (nframes - user_td) >= (double)pcalcs / 100.0)
            {
                fprintf(stderr, "\rApproximately %i percent complete.", pcalcs);
                fflush(stderr);
                pcalcs++;
            }
        }
        ffclose(out);
    }
    // Check to see if any other output options were chosen.
    if (!(bAvg || bVar || bMax || bPair || bISDMat))
    {
        // No other options chosen, so end the program.
        thanx(stderr);
        return 0;
    }
    
    
    
    /* Create arrays to hold output per frame.
     */
    snew(avgISDF,nframes);
    avgISD = 0.0;
    snew(maxISDF,nframes);
    maxISD = 0.0;
    if (bPair)
    {
        snew(maxframe,nframes);
    }
    if (bVar)
    {
        snew(varISDF, nframes);
        varISD = 0.0;
    }
    // Needed to store ISD matrix.
    if (bISDMat)
    {
        snew(ISDMat, nframes);
        for (i = 0; i < nframes; i++)
        {
            snew(ISDMat[i], nframes);
        }
    }
    // Needed for the map figure.
    if (bMap)
    {
        snew(xticks, nframes);
        snew(yticks, nframes);
        for (i = 0; i < nframes; i++)
        {
            xticks[i] = (real)i;
            yticks[i] = (real)i;
        }
    }
    
    
    
    /* Main calculation loop.
     */
    fprintf(stderr,"\nCalculating results of means. \n");
    pcalcs = 1;
    
    /* Originally this was designed to only loop through each pair of i and j
     * one time to save half of the calculations. Eventually it became 
     * impractical to make sure that each ISDM was symmetric, so now the
     * algorithm takes the performance hit in favor of accuracy and simplicity.
     */
    
    // Loop through reference frames.
    for (i = (bf - 1); i < ef; i++)
    {
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
                bPHIPSI2 || bANGDIH2G || bGMRG || bMRMS || bRROT)
            {
                ISD = call_ISDM(iatoms, cframe, rframe, ISDM);
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
            
            
            
            // Store ISD to the ISD matrix.
            if (bISDMat)
            {
                ISDMat[i][j] = ISD;
            }
            // Calculate summations with doubles for improved accuracy.
            dISD = (double)ISD;
            // Update mean and max.
            avgISDF[i] += dISD;
            if (dISD > maxISDF[i])
            {
                maxISDF[i] = dISD;
                if (bPair)
                {
                    maxframe[i] = j;
                }
            }
            // If calculating variance, calculate the sum of squares.
            if (bVar)
            {
                varISDF[i] += (dISD * dISD);
            }
        }
        
        // Divide by N-1.
        avgISDF[i] /= (nframes - 1);
        if (bVar)
        {
            // Unbiased sample variance.
            varISDF[i] /= (nframes - 2);
        }
        
        // Update progress output.
        while ((double)(i + 2 - bf) / (ef - bf + 1) >= (double)pcalcs/100)
        {
            fprintf(stderr,"\rApproximately %i percent complete.", pcalcs);
            fflush(stderr);
            pcalcs++;
        }
    }
    
    // Print output to specified files.
    if (bAvg)
    {
        out=xvgropen(opt2fn("-avg", NFILE, fnm), 
                     "ISD Analysis", 
                     "Reference Frame", 
                     "Average Difference", 
                     oenv);
        
        for (i = (ref_bf - 1); i < ref_ef; i++)
        {
            fprintf(out,"%-6i %10f \n", (i+1), avgISDF[i]);
        }
        ffclose(out);
    }
    
    if (bVar)
    {
        out = xvgropen(opt2fn("-var", NFILE, fnm), 
                       "ISD Analysis", 
                       "Reference Frame", 
                       "Variance of Difference", 
                       oenv);
        
        double vCorrect = (double)(nframes - 1.0) / (nframes - 2.0);
        for (i = (ref_bf - 1); i < ref_ef; i++)
        {
            // Unbiased sample variance requires a correction.
            msqISD = avgISDF[i] * avgISDF[i] * vCorrect;
            fprintf(out,"%-6i %12.8f \n", (i+1), (varISDF[i] - msqISD));
        }
        ffclose(out);
    }
    
    if (bMax)
    {
        out=xvgropen(opt2fn("-max", NFILE, fnm), 
                     "ISD Analysis", 
                     "Reference Frame", 
                     "Maximum Difference", 
                     oenv);
        
        for (i = (ref_bf - 1); i < ref_ef; i++)
        {
            fprintf(out,"%-6i %12.8f \n", (i+1), maxISDF[i]);
        }
        ffclose(out);
    }
    
    if (bPair)
    {
        out=xvgropen(opt2fn("-pair", NFILE, fnm), 
                     "Paired Maximally Distant Structures", 
                     "Reference Frame", 
                     "Maximally Different Frame", 
                     oenv);
        
        for (i = (ref_bf - 1); i < ref_ef; i++)
        {
            fprintf(out,"%-6i %-6i \n", (i+1), (maxframe[i] + 1));
        }
        ffclose(out);
    }
    
    if (bISD)
    {
        // Opens the output file.
        out = opt2FILE("-isd", NFILE, fnm, "w");
        
        // Write output.
        for (i = 0; i < nframes; i++)
        {
            fprintf(out, "%12.8f", ISDMat[i][0]);
            for (j = 1; j < nframes; j++)
            {
                fprintf(out, ",%12.8f", ISDMat[i][j]);
            }
            fprintf(out, "\n");
        }
        
        // Close the output file.
        ffclose(out);
    }
    
    if (bMap)
    {
        // Opens the output file.
        out = opt2FILE("-map", NFILE, fnm, "w");
        // Needed for write_xpm.
        rlo.r = 0; rlo.g = 0; rlo.b = 1;
        rmd.r = 1; rmd.g = 1; rmd.b = 0;
        rhi.r = 1; rhi.g = 0; rhi.b = 0;
        
        // Need overall average and max ISD for ALL frames.
        for (i = 0; i < nframes; i++)
        {
            avgISD += avgISDF[i];
            if (maxISD < maxISDF[i])
            {
                maxISD = maxISDF[i];
            }
        }
        avgISD /= nframes;
        
        // Choose the maximum value for the xpm.
        if (setmax == -1.0)
        {
            xpm_max = (real)maxISD;
        }
        else
        {
            if (setmax > avgISD)
            {
                xpm_max = setmax;
            }
            else
            {
                fprintf(stderr, "\nWarning: the argument for -setmax must be "
                                "greater than the average ISD. The maximum "
                                "map value will be reset to the default.\n");
            }
        }
        
        // Write to the output file.
        sprintf(buf,"Frame vs Frame");
        unsigned int tmp = 0;
        write_xpm_header(out, buf, "ISD", "Frame", "Frame", FALSE);
        write_mat_levels3(out, nframes, nframes, &nlevels, 0, (real)avgISD, xpm_max, rlo, rmd, rhi);
        write_xpm_axis(out, "x", tmp & MAT_SPATIAL_X, nframes, xticks);
        write_xpm_axis(out, "y", tmp & MAT_SPATIAL_Y, nframes, yticks);
        write_mat_pixels3(out, nframes, nframes, ISDMat, 0, (real)avgISD, xpm_max, nlevels);
        
        // Close the output file.
        ffclose(out);
        
        // Reset some values.
        avgISD = 0.0;
        maxISD = 0.0;
    }
    
    if (bCalcDCR)
    {
        snew(DCR, nf2);
        for (ij = 1; ij < nf2; ij++)
        {
            for (i = 0; i < (nframes - ij); i++)
            {
                j = i + ij;
                DCR[ij] += (double)ISDMat[i][j];
            }
            DCR[ij] /= (nframes - ij);
        }
    }
    
    if (bDCR)
    {
        out = xvgropen(opt2fn("-decorr", NFILE, fnm), 
                       "ISD Decorrelation", 
                       "Time From Reference Frame (ns)", 
                       "Mean ISD", 
                       oenv);
        
        for (i = 0; i < nf2; i++)
        {
            fprintf(out, "%10f %12.8f \n", (real)i * dt / 1000.0, DCR[i]);
        }
        ffclose(out);
    }
    
    if (bVDCR || bSNR)
    {
        // Solve for the variance of the decorrelation.
        snew(VDCR, nf2);
        double vCorrect;
        for (ij = 1; ij < nf2; ij++)
        {
            for (i = 0; i < (nframes - ij); i++)
            {
                j = i + ij;
                dISD = (double)ISDMat[i][j];
                VDCR[ij] += dISD * dISD;
            }
            // Unbiased sample variance
            VDCR[ij] /= (nframes - ij - 1);
            vCorrect  = (double)(nframes - ij) / (nframes - ij - 1.0);
            VDCR[ij] -= DCR[ij] * DCR[ij] * vCorrect;
        }
    }
    
    if (bVDCR)
    {
        out = xvgropen(opt2fn("-vdcr", NFILE, fnm), 
                       "Variance of ISD Decorrelation", 
                       "Time From Reference Frame (ns)", 
                       "Variance of ISD", 
                       oenv);
        
        for (i = 0; i < nf2; i++)
        {
            fprintf(out, "%10f %12.8f \n", (real)i * dt / 1000.0, VDCR[i]);
        }
        ffclose(out);
    }
    
    if (bSNR)
    {
        out = xvgropen(opt2fn("-snr", NFILE, fnm), 
                       "ISD Decorrelation SNR", 
                       "Time From Reference Frame (ns)", 
                       "SNR: Sample Mean / Corrected Sample Standard Deviation", 
                       oenv);
        
        snew(SNR, nf2);
        for (i = 1; i < nf2; i++)
        {
            SNR[i] = DCR[i] / sqrt(VDCR[i]);
            fprintf(out, "%10f %12.8f \n", (real)i * dt / 1000.0, SNR[i]);
        }
        ffclose(out);
    }
    
    if (bAvgSCL)
    {
        out = xvgropen(opt2fn("-avgscl", NFILE, fnm), 
                       "Decorrelation of ISD", 
                       "Time From Reference Frame (ns)", 
                       "Mean ISD / Mean Decorrelated ISD", 
                       oenv);
        
        for (i = 0; i < nf2; i++)
        {
            avgSCL += DCR[i];
        }
        avgSCL /= nf2;
        
        for (i = 0; i < nf2; i++)
        {
            fprintf(out, "%10f %12.8f \n", (real)i * dt / 1000.0, DCR[i] / avgSCL);
        }
        ffclose(out);
    }
    
    if (bMaxSCL)
    {
        out = xvgropen(opt2fn("-maxscl", NFILE, fnm), 
                       "Decorrelation of ISD", 
                       "Time From Reference Frame (ns)", 
                       "Mean ISD / Maximum Decorrelated ISD", 
                       oenv);
        
        maxSCL = -1.0;
        for (i = 0; i < nf2; i++)
        {
            if (DCR[i] > (double)maxSCL)
            {
                maxSCL = (double)DCR[i];
            }
        }
        
        for (i = 0; i < nf2; i++)
        {
            fprintf(out, "%10f %12.8f \n", (real)i * dt / 1000.0, DCR[i] / maxSCL);
        }
        ffclose(out);
    }
    
    if (bMinDCR)
    {
        out = xvgropen(opt2fn("-mindcr", NFILE, fnm), 
                       "Decorrelation of Minimum ISD", 
                       "Time From Reference Frame (ns)", 
                       "Minimum ISD", 
                       oenv);
        
        for (ij = 1; ij < nf2; ij++)
        {
            minDCR = 1000000.0;
            for (i = 0; i < (nframes - ij); i++)
            {
                j = i + ij;
                if (ISDMat[i][j] < minDCR)
                {
                    minDCR = ISDMat[i][j];
                }
            }
            fprintf(out, "%10f %12.8f \n", (real)ij * dt / 1000.0, minDCR);
        }
        ffclose(out);
    }
    
    if (bMaxDCR)
    {
        out = xvgropen(opt2fn("-maxdcr", NFILE, fnm), 
                       "Decorrelation of Minimum ISD", 
                       "Time From Reference Frame (ns)", 
                       "Minimum ISD", 
                       oenv);
        
        for (ij = 1; ij < nf2; ij++)
        {
            maxDCR = -1.0;
            for (i = 0; i < (nframes - ij); i++)
            {
                j = i + ij;
                if (ISDMat[i][j] > maxDCR)
                {
                    maxDCR = ISDMat[i][j];
                }
            }
            fprintf(out, "%10f %12.8f \n", (real)ij * dt / 1000.0, maxDCR);
        }
        ffclose(out);
    }
    
    if (bSens)
    {
        oobc = 0;
        double DCRi, DCRj, DCRij, DCRdt, ijSens;
        snew(sens, (nf2 - user_tdsens));
        sens[0] = 1.0;
        for (i = 1; i < (nf2 - user_tdsens); i++)
        {
            for (ij = 1; ij < user_tdsens; ij++)
            {
                j = i + ij;
                DCRi  = DCR[i]; DCRj = DCR[j]; DCRij = DCR[ij];
                DCRdt = DCRj - DCRi; ijSens = DCRdt / DCRij;
                if (DCRij == 0.0)
                {
                    if (DCRdt > 0.0)
                    {
                        sens[i] += 1.0;
                        oobc++;
                    }
                    else if (DCRdt < 0.0)
                    {
                        sens[i] -= 1.0;
                        oobc++;
                    }
                }
                else
                {
                    if (ijSens > 1.0)
                    {
                        sens[i] += 1.0;
                        oobc++;
                    }
                    else if (ijSens < -1.0)
                    {
                        sens[i] -= 1.0;
                        oobc++;
                    }
                    else
                    {
                        sens[i] += ijSens;
                    }
                }
            }
            sens[i] /= user_tdsens - 1.0;
        }
        
        out = xvgropen(opt2fn("-sens", NFILE, fnm), 
                       "Sensitivity Analysis", 
                       "Time From Reference Frame (ns)", 
                       "Sensitivity", 
                       oenv);
        
        for (i = 0; i < (nf2 - user_tdsens); i++)
        {
            fprintf(out, "%10f %12.8f \n", (real)i * dt / 1000.0, sens[i]);
        }
        ffclose(out);
    }
    
    
    
    // Sum all mean diff.
    for (i = (ref_bf - 1); i < ref_ef; i++)
    {
        // Sum up the per frame totals to an overall total.
        avgISD += avgISDF[i];
        if (maxISDF[i] > maxISD)
        {
            maxISD = maxISDF[i];
            if (bPair)
            {
                maxframei = i;
            }
        }
    }
    avgISD /= (ref_ef - ref_bf + 1);
    
    if (bVar)
    {
        double vCorrect = (double)(ref_ef - ref_bf + 1.0) * (nframes - 1.0);
        for (i = (ref_bf - 1); i < ref_ef; i++)
        {
            varISD += (nframes - 2.0) * varISDF[i];
        }
        varISD /= vCorrect - 1.0;
        varISD -= (avgISD * avgISD) * vCorrect / (vCorrect - 1.0);
    }
    
    // Output the final information.
    printf("\n\nAverage ISD: %12.8f \n", avgISD);
    if (bVar)
    {
        printf("Variance of ISD: %12.8f \n", varISD);
        printf("SNR of ISD: %12.8f \n", avgISD / sqrt(varISD));
    }
    printf("Maximum ISD: %12.8f \n", maxISD);
    if (bPair)
    {
        printf("Maximally Different Pair of Frames: %-6i %-6i \n", 
               (maxframei + 1), 
               (maxframe[maxframei] + 1));
    }
    if (bSNR)
    {
        for (i = 1; i < nf2; i++)
        {
            dSNR += SNR[i];
        }
        dSNR /= nf2 - 1.0;
        printf("\nSNR of ISD Decorrelation: %12.8f \n", dSNR);
    }
    if (bSens)
    {
        int nSens = (nf2 - user_tdsens - 1) * (user_tdsens - 1);
        printf("\nNumber of out of bounds corrections is %-10i out of %-10i "
               "possible. \n", oobc, nSens);
    }
    
    // Not sure what this one does. Sends to xmgrace?
    do_view(oenv,ftp2fn(efXVG,NFILE,fnm),"-nxy");

    // Closing.
    thanx(stderr);
    return 0;
}



/* The below section is from gromacs code which was released separately under 
 * the LGPL version 2 license.
 */

void write_mat_levels(FILE *out,int n_x, int n_y,int *nlevels,real lo,real hi,t_rgb rlo,t_rgb rhi)
{
    static const char cmap[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrst"
                               "uvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/?";
    int    lenMap = 88;
    int    i,nlo;
    real   invlevel,r,g,b;
    
    if (*nlevels > lenMap*lenMap) {
        fprintf(stderr,"Warning, too many levels (%d) in matrix, using %d only\n",
                *nlevels,(int)(lenMap*lenMap));
        *nlevels=lenMap*lenMap;
    }
    else if (*nlevels < 2) {
        fprintf(stderr,"Warning, too few levels (%d) in matrix, using 2 instead\n",*nlevels);
        *nlevels=2;
    }
    
    fprintf(out,"static char *gromacs_xpm[] = {\n");
    fprintf(out,"\"%d %d   %d %d\",\n",
            n_x,n_y,*nlevels,(*nlevels <= lenMap) ? 1 : 2);
    
    invlevel=1.0/(*nlevels-1);
    for(i=0; (i<*nlevels); i++) {
        nlo=*nlevels-1-i;
        r=(nlo*rlo.r+i*rhi.r)*invlevel;
        g=(nlo*rlo.g+i*rhi.g)*invlevel;
        b=(nlo*rlo.b+i*rhi.b)*invlevel;
        fprintf(out,"\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
                cmap[i % lenMap],(*nlevels <= lenMap) ? ' ' : cmap[i/lenMap],
                                  (unsigned int)((int)((255*r) + 0.5)),
                                  (unsigned int)((int)((255*g) + 0.5)),
                                  (unsigned int)((int)((255*b) + 0.5)),
                                  (nlo*lo+i*hi)*invlevel);
    }
}



void write_mat_levels3(FILE *out,int n_x,int n_y,int *nlevels,real lo,real mid,real hi,t_rgb rlo,t_rgb rmid,t_rgb rhi)
{
    static const char cmap[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrst"
                               "uvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/?";
    int    lenMap = 88;
    int    i,nmid;
    real   r,g,b,clev_lo,clev_hi;
    
    if (*nlevels > lenMap*lenMap) {
        fprintf(stderr,"Warning, too many levels (%d) in matrix, using %d only\n",
                *nlevels,(int)(lenMap*lenMap));
        *nlevels=lenMap*lenMap;
    }
    else if (*nlevels < 2) {
        fprintf(stderr,"Warning, too few levels (%d) in matrix, using 2 instead\n",
                *nlevels);
        *nlevels=2;
    }   
    if (!((mid >= lo) && (mid < hi)))
        gmx_fatal(FARGS,"Lo: %f, Mid: %f, Hi: %f\n",lo,mid,hi);
    
    fprintf(out,"static char *gromacs_xpm[] = {\n");
    fprintf(out,"\"%d %d   %d %d\",\n",
            n_x,n_y,*nlevels,(*nlevels <= lenMap) ? 1 : 2);
    
    nmid    = min(max(0,((mid-lo)/(hi-lo))*((*nlevels)-1)),(*nlevels)-1);
    clev_lo = nmid;
    clev_hi = (*nlevels - 1 - nmid);
    for(i=0; (i<nmid); i++) {
        r   = rlo.r+(i*(rmid.r-rlo.r)/clev_lo);
        g   = rlo.g+(i*(rmid.g-rlo.g)/clev_lo);
        b   = rlo.b+(i*(rmid.b-rlo.b)/clev_lo);
        fprintf(out,"\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
                cmap[i % lenMap],
                (*nlevels <= lenMap) ? ' ' : cmap[i/lenMap],
                 (unsigned int)((int)((255*r) + 0.5)),
                 (unsigned int)((int)((255*g) + 0.5)),
                 (unsigned int)((int)((255*b) + 0.5)),
                 ((nmid - i)*lo + i*mid)/clev_lo);
    }
    for(i=0; (i<(*nlevels-nmid)); i++) {
        r   = rmid.r+(i*(rhi.r-rmid.r)/clev_hi);
        g   = rmid.g+(i*(rhi.g-rmid.g)/clev_hi);
        b   = rmid.b+(i*(rhi.b-rmid.b)/clev_hi);
        fprintf(out,"\"%c%c c #%02X%02X%02X \" /* \"%.3g\" */,\n",
                cmap[(i+nmid) % lenMap],
                (*nlevels <= lenMap) ? ' ' : cmap[(i+nmid)/lenMap],
                 (unsigned int)((int)((255*r) + 0.5)),
                 (unsigned int)((int)((255*g) + 0.5)),
                 (unsigned int)((int)((255*b) + 0.5)),
                 ((*nlevels - 1 - nmid - i)*mid + i*hi)/clev_hi);
    }
}



void write_mat_pixels(FILE *out,int n_x,int n_y,real **matrix,real lo,real hi,int nlevels)
{
    static const char cmap[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrst"
                               "uvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/?";
    int  lenMap = 88;
    int  i,j,c;
    real invlevel;
    
    invlevel=(nlevels-1)/(hi-lo);
    //printf("\n");
    for(j=n_y-1; (j>=0); j--) {
        if(j%(1+n_y/100)==0) 
            fprintf(stderr,"%3d%%\b\b\b\b",(100*(n_y-j))/n_y);
        fprintf(out,"\"");
        for(i=0; (i<n_x); i++) {
            c=gmx_nint((matrix[i][j]-lo)*invlevel);
            if (c<0) c=0;
            if (c>=nlevels) c=nlevels-1;
            if (nlevels <= lenMap)
            {
                fprintf(out,"%c",cmap[c]);
                //printf("%c",matmap[c]);
            }
            else
            {
                fprintf(out,"%c%c",cmap[c % lenMap],cmap[c / lenMap]);
            }
        }
        if (j > 0)
            fprintf(out,"\",\n");
        else
            fprintf(out,"\"\n");
        //printf("\n");
    }
}



void write_mat_pixels3(FILE *out,int n_x,int n_y,real **matrix,real lo,real mid,real hi,int nlevels)
{
    static const char cmap[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrst"
                               "uvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/?";
    int  lenMap = 88;
    int  i,j,c=0,nmid;
    real invlev_lo,invlev_hi;
    
    nmid = min(max(0,((mid-lo)/(hi-lo))*((nlevels)-1)),(nlevels)-1);
    invlev_hi=(nlevels-1-nmid)/(hi-mid);
    invlev_lo=(nmid)/(mid-lo);
    
    for(j=n_y-1; (j>=0); j--) {
        if(j%(1+n_y/100)==0) 
            fprintf(stderr,"%3d%%\b\b\b\b",(100*(n_y-j))/n_y);
        fprintf(out,"\"");
        for(i=0; (i<n_x); i++) {
            if (matrix[i][j] >= mid)
                c=nmid+gmx_nint((matrix[i][j]-mid)*invlev_hi);
            else if (matrix[i][j] >= lo)
                c=gmx_nint((matrix[i][j]-lo)*invlev_lo);
            else
                c = 0;
            
            if (c<0) 
                c=0;
            if (c>=nlevels) 
                c=nlevels-1;
            if (nlevels <= lenMap)
                fprintf(out,"%c",cmap[c]);
            else
                fprintf(out,"%c%c",cmap[c % lenMap],cmap[c / lenMap]);
        }
        if (j > 0)
            fprintf(out,"\",\n");
        else
            fprintf(out,"\"\n");
    }
}

/* The above section is from gromacs code which was released separately under 
 * the LGPL version 2 license.
 */