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




int gmx_isdcalc(int argc,char *argv[])
{
    const char *desc[] = {
        "[TT]g_isdcalc[tt] implements a list of measures designed to ",
        "compare two frames from a trajectory and quantify their difference. ",
        "This tool analyzes the input by comparing each pair of frames in ",
        "the trajectory and recording some statistics. For each reference ",
        "frame, the mean and maximum interstructure distance (ISD) with ",
        "all other frames is recorded. These statistics are optional ",
        "outputs. A subset of frames can be used as reference structures ",
        "(still compared to all other frames) with the bframe and eframe ",
        "options. Variance calculations are optional. An overall mean and ",
        "max difference is sent to stdout. The default measure if one is not ",
        "chosen by the user is RMSD. Only one measure can be chosen at a ",
        "time. The measures of inter-structure distance are called ISDMs."
    };
    
    
    
    static gmx_bool bANG=FALSE, bDIH=FALSE, bANGDIH=FALSE, bDRMS=FALSE;
    static gmx_bool bPHIPSI=FALSE, bSRMS=FALSE, bPCOR=FALSE, bMAMMOTH=FALSE;
    static gmx_bool bACOR=FALSE, bESA=FALSE, bRMSD=FALSE, bMIR=FALSE;
    static gmx_bool bRG=FALSE, bSRG=FALSE, bE2E=FALSE, bSE2E=FALSE;
    static gmx_bool bRROT=FALSE;
    static int      bFrame = -1, eFrame = -1;
    t_pargs pa[] = {
        { "-ang", FALSE, etBOOL, {&bANG},
            "ISDM: Mean cosine of difference of backbone angles for each "
            "set of three atoms. Assumes only CA atoms." },
        { "-dih", FALSE, etBOOL, {&bDIH},
            "ISDM: Mean cosine of difference of backbone dihedrals for "
            "each set of four atoms. Assumes only CA atoms." },
        { "-angdih", FALSE, etBOOL, {&bANGDIH},
            "ISDM: Geometric mean of ang and dih measures." },
        { "-phipsi", FALSE, etBOOL, {&bPHIPSI},
            "ISDM: Mean cosine of difference of phi and psi angles. "
            "Assumes only backbone atoms." },
        { "-drms", FALSE, etBOOL, {&bDRMS},
            "ISDM: Mean difference of the paired distances matrix for "
            "all atoms. Distance RMS(D)" },
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
            "code ported from Matlab to C. For this measure, please cite: \n\n"
            "Liu W, Srivastava A, Zhang J (2011) A Mathematical Framework "
            "for Protein Structure Comparison. PLoS Comput Biol 7(2): "
            "e1001075.\n" },
        { "-bframe", FALSE, etINT, {&bFrame},
            "Compare range of frames from bframe to eframe to all other " 
            "frames. The bframe and eframe options are applied after the b, " 
            "e, and dt options and use units of frames instead of units of " 
            "time. Frame numbers are counted from one." 
        },
        { "-eframe", FALSE, etINT, {&eFrame},
            "Compare range of frames from bframe to eframe to all other " 
            "frames. The bframe and eframe options are applied after the b, " 
            "e, and dt options and use units of frames instead of units of " 
            "time. Frame numbers are counted from one." 
        },
    };
    
    
    
    FILE       *out;
    t_trxstatus *status;
    t_topology top;
    int        ePBC;
    rvec       *x, *iframe, *jframe, *rframe, *cframe, rrot_xyz, xold;
    rvec       **frames;
    real       *nweights, *iweights;
    real       ISD, maxISD, avgISD, varISD, msqi;
    real       *diff, *maxdiff, *avgdiff, *vardiff;
    matrix     box, rrot, rrotx, rroty, rrotz;
    real       t, pi = 3.14159265358979;
    int        *maxframe, *rnum;
    int        i, j, k, m, n, bf, ef, iatoms, natoms, nframes, maxframei;
    int        pcalcs, noptions;
    gmx_bool   bDFLT, bMean, bVar, bMax, bPair, bFit;
    char       *ISDM, *grpname, title[256], title2[256], *rname;
    atom_id    *index;
    output_env_t oenv;
    gmx_rmpbc_t  gpbc=NULL;
    const char *leg[]  = { "D" }; 
#define NLEG asize(leg) 
    t_filenm fnm[] = {
        { efTRX, "-f",    NULL,   ffREAD }, 
        { efTPS, NULL,    NULL,   ffREAD },
        { efNDX, NULL,    NULL,   ffOPTRD },
        { efXVG, "-mean", "mean", ffOPTWR },
        { efXVG, "-var",  "var",  ffOPTWR },
        { efXVG, "-max",  "max",  ffOPTWR },
        { efXVG, "-pair", "pair", ffOPTWR },
    }; 
#define NFILE asize(fnm)
    int npargs;
    
    CopyRight(stderr,argv[0]);
    npargs = asize(pa);
    
    // Lots of black magic with this one. The oenv is used by many things.
    parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
                      NFILE,fnm,npargs,pa,asize(desc),desc,0,NULL,&oenv);
    
    
    /* Reads the tpr file. Outputs a ton of info.
     * 
     * I think this is the line that forces you to have a -s at prompt.
     */
    read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&x,NULL,box,TRUE);
    
    // Asks you to choose a selection of atoms at prompt.
    get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&iatoms,&index,&grpname);
    
    // If there are no options at command line, do default behavior.
    bDFLT = !(bANG || bDIH || bANGDIH || bPHIPSI || bDRMS || bSRMS || bRMSD || 
              bPCOR || bACOR || bMAMMOTH || bESA || bRG || bSRG || bE2E || 
              bSE2E || bMIR || bRROT);
    
    bFit  =  (bDFLT || bRMSD || bMIR || bSRMS || bPCOR);
    
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
    
    if (bPHIPSI)
    {
        fprintf(stderr,"\nUsing phi and psi angles as ISDM.\n");
        ISDM = "PHIPSI";
        noptions++;
    }
    
    if (bDRMS)
    {
        fprintf(stderr,"\nUsing distance RMS as ISDM.\n");
        ISDM = "DRMS";
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
    
    if (bANGDIH)
    {
        fprintf(stderr,"\nUsing geometric mean of angles and dihedrals as ISDM.\n");
        ISDM = "ANGDIH";
        noptions++;
    }
    
    if (bRROT)
    {
        fprintf(stderr,"\nUsing RMSD with random rotation as ISDM.\n");
        noptions++;
        
        // Additional stuff for option.
        srand(time(NULL));
        // Use up the first few random numbers that usually aren't random.
        rrot_xyz[0] = (real)rand();
        rrot_xyz[1] = (real)rand();
        rrot_xyz[2] = (real)rand();
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
    
    // Output which files?
    bMean = opt2bSet("-mean", NFILE, fnm);
    bVar  = opt2bSet("-var",  NFILE, fnm);
    bMax  = opt2bSet("-max",  NFILE, fnm);
    bPair = opt2bSet("-pair", NFILE, fnm);
    
    
    
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
    snew(diff,iatoms);
    
    // Initialize nweights to zeros.
    for (i=0; i<natoms; i++)
    {
        nweights[i] = 0;
    }
    // Makes an array of weights. Necessary for reset_x.
    for (i=0; i<iatoms; i++)
    {
        // Give a value for the weights.
        nweights[(int)index[i]] = 1;
        iweights[i] = 1;
        // While we're at it, initialize diff to zeros.
        diff[i] = 0;
    }
    
    
    nframes = 0;
    do
    {
        /* This loop doesn't do anything.
         * 
         * It's just the most reliable way to find the number of frames.
         */
        nframes++;
    } while(read_next_x(oenv,status,&t,natoms,x,box));
    // Close the trajectory.
    close_trj(status);
    // Throw an error if there aren't enough frames.
    if (nframes < 2)
    {
        gmx_fatal(FARGS,"\nThe trajectory must have at least 2 frames.\n");
    }
    
    // Check to see if the user gives a range of frames.
    if ((bFrame == -1) && (eFrame == -1))
    {
        bf = 1;
        ef = nframes;
    }
    else
    {
        if (bFrame == -1)
        {
            // Set bf to default.
            bf = 1;
        }
        else
        {
            // Error check for first frame being within range 1 to nframes.
            if ((bFrame < 1) || (bFrame > nframes))
            {
                gmx_fatal(FARGS,"\nArgument to bframe must be between 1 and last frame.\n");
            }
            
            // Assign user input to bf.
            bf = bFrame;
        }
        
        if (eFrame == -1)
        {
            // Set ef to default.
            ef = nframes;
        }
        else
        {
            // Error check for last frame being within range 1 to nframes.
            if ((eFrame < 1) || (eFrame > nframes))
            {
                gmx_fatal(FARGS,"\nArgument to eframe must be between 1 and last frame.\n");
            }
            
            // Assign user input to ef.
            ef = eFrame;
        }
        
        if ((bFrame != -1) && (eFrame != -1))
        {
            // Error check for last frame being greater than first frame.
            if (bFrame > eFrame)
            {
                gmx_fatal(FARGS,"\nValue of eframe must be greater than or equal to bframe\n");
            }
        }
    }   
    
    // Create an array to hold all frames.
    snew(frames,nframes);
    
    /* Create arrays to hold output per frame.
     */
    snew(avgdiff,nframes);
    snew(maxdiff,nframes);
    snew(maxframe,nframes);
    // Only allocate if variance is being calculated.
    if (bVar)
    {
        snew(vardiff, nframes);
    }
    
    
    
    /* Opens trj. Reads first frame. Returns status. Allocates mem for x.
     * 
     * Not sure which argument determines which atoms to pull info for.
     */
    fprintf(stderr,"\nStoring trajectory to memory.\n");
    natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
    
    // Initialize index to keep track of current frame.
    i = 0;
    // This is for removing periodic boundary conditions.
    gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms, box);
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
    
    // Initialize output to 0.
    maxISD = 0;
    avgISD = 0;
    
    
    /* Main calculation loop.
     */
    fprintf(stderr,"\nCalculating results. \n");
    
    /* Originally this was designed to only loop through each pair of i and j
     * one time to save half of the calculations. Eventually it became 
     * impractical to make sure that each ISDM was symmetric, so now the
     * algorithm takes the performance hit in favor of accuracy and simplicity.
     */
    pcalcs = 1;
    
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
                bPCOR || bACOR)
            {
                ISD = call_ISDM(iatoms, cframe, rframe, diff, ISDM);
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
                ISD = sqrt(calc_msd(iatoms, cframe, rframe, diff));
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
            
            
            // Update sum for the means of ith frame.
            avgdiff[i] += ISD;
            // Update the max for ith frame.
            if (ISD > maxdiff[i])
            {
                maxdiff[i]  = ISD;
                maxframe[i] = j;
            }
            
            // If calculating variance, calculate the sum of squares.
            if (bVar)
            {
                vardiff[i] += (ISD * ISD);
            }
        }
        
        // Divide by N-1.
        avgdiff[i] /= (nframes - 1);
        if (bVar)
        {
            vardiff[i] /= (nframes - 1);
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
    if (bMean)
    {
        out=xvgropen(opt2fn("-mean", NFILE, fnm), 
                     "ISD Analysis", 
                     "Reference Frame", 
                     "Mean Difference", 
                     oenv);
        
        for (i = (bf - 1); i < ef; i++)
        {
            fprintf(out,"%-6i %10f \n", (i+1), avgdiff[i]);
        }
        ffclose(out);
    }
    
    if (bVar)
    {
        out=xvgropen(opt2fn("-var", NFILE, fnm), 
                     "ISD Analysis", 
                     "Reference Frame", 
                     "Variance of Difference", 
                     oenv);
        
        for (i = (bf - 1); i < ef; i++)
        {
            msqi = avgdiff[i] * avgdiff[i];
            fprintf(out,"%-6i %10f \n", (i+1), (vardiff[i] - msqi));
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
        
        for (i = (bf - 1); i < ef; i++)
        {
            fprintf(out,"%-6i %10f \n", (i+1), maxdiff[i]);
        }
        ffclose(out);
    }
    
    if (bPair)
    {
        out=xvgropen(opt2fn("-pair", NFILE, fnm), 
                     "ISD Analysis", 
                     "Reference Frame", 
                     "Maximally Different Frame", 
                     oenv);
        
        for (i = (bf - 1); i < ef; i++)
        {
            fprintf(out,"%-6i %-6i \n", (i+1), (maxframe[i] + 1));
        }
        ffclose(out);
    }
    
    // Sum all mean diff.
    for (i = (bf - 1); i < ef; i++)
    {
        // Sum up the per frame totals to an overall total.
        avgISD += avgdiff[i];
        if (maxdiff[i] > maxISD)
        {
            maxISD  = maxdiff[i];
            maxframei = i;
        }
    }
    avgISD /= (ef - bf + 1);
    
    if (bVar)
    {
        for (i = (bf - 1); i < ef; i++)
        {
            varISD += vardiff[i];
        }
        varISD /= (ef - bf + 1);
        varISD -= (avgISD * avgISD);
    }
    
    // Output the final information.
    printf("\n\nMean Difference: %10f \n", avgISD);
    if (bVar)
    {
        printf("Variance of Difference: %10f \n", varISD);
    }
    printf("Maximum Difference: %10f \n", maxISD);
    printf("Maximally Different Pair of Frames: %-6i %-6i \n", 
           (maxframei + 1), 
           (maxframe[maxframei] + 1));
    
    // Not sure what this one does. Sends to xmgrace?
    do_view(oenv,ftp2fn(efXVG,NFILE,fnm),"-nxy");

    // Closing.
    thanx(stderr);
    return 0;
}
