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



int gmx_isddecorr(int argc,char *argv[])
{
    const char *desc[] = {
        "[TT]g_isddecorr[tt] measures the mean interstructure ",
        "distance (ISD) for all structures separated by each dt. Options ",
        "scaled and mscaled divide the decorrelation results by the mean ",
        "and max ISD respectively. Option ddecorr measures the mean ratio of ",
        "derivatives for all structures separated by each dt. The double ",
        "is similar in purpose to ddecorr, but the calculation is more ",
        "to performing a second decorrelation on the original output data. ",
        "Option tdo calculates the difference between paired structures ",
        "separated by td time. The method used to compare each pair of "
        "structures is chosen by a variety of ISDMs chosen by flags at the "
        "command line. Only one ISDM can be chosen at a time."
    };
    
    
    
    static gmx_bool bANG=FALSE, bDIH=FALSE, bANGDIH=FALSE, bDRMS=FALSE;
    static gmx_bool bPHIPSI=FALSE, bSRMS=FALSE, bPCOR=FALSE, bMAMMOTH=FALSE;
    static gmx_bool bACOR=FALSE, bESA=FALSE, bRMSD=FALSE, bMIR=FALSE;
    static gmx_bool bRG=FALSE, bSRG=FALSE, bE2E=FALSE, bSE2E=FALSE;
    static gmx_bool bRROT=FALSE, bSDRMS=FALSE;
    static real user_td = -1;
    t_pargs pa[] = {
        { "-ang", FALSE, etBOOL, {&bANG},
        "ISDM: Mean cosine of difference of backbone angles for each "
        "set of three atoms. Assumes only CA atoms." },
        { "-dih", FALSE, etBOOL, {&bDIH},
        "ISDM: Mean cosine of difference of backbone dihedrals for "
        "each set of four atoms. Assumes only CA atoms." },
        { "-angdih", FALSE, etBOOL, {&bANGDIH},
        "ISDM: Geometric mean of ang and dih ISDMs." },
        { "-phipsi", FALSE, etBOOL, {&bPHIPSI},
        "ISDM: Mean cosine of difference of phi and psi angles. "
        "Assumes only backbone atoms." },
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
//         { "-mammoth", FALSE, etBOOL, {&bMAMMOTH},
//         "ISDM: MAMMOTH (MAtching Molecular Models Obtained from "
//         "Theory). Compares segments of residues chosen by sequence "
//         "alignment. Attempts to focus on correct secondary structure "
//         "moreso than tertiary structure. Source code modified for "
//         "compatibility. For this ISDM, please cite: \n\n"
//         "Ortiz, AR, Strauss, CE, Olmea, O (2002). MAMMOTH "
//         "(Matching molecular models obtained from theory): An automated "
//         "method for model comparison. Protein Sci. 11 (11), 2606–2621.\n"},
        { "-esa", FALSE, etBOOL, {&bESA},
        "ISDM: Elastic shape analysis. Based on image analysis. "
        "Warps structure onto the reference structure. Original source "
        "code ported from Matlab to C. For this ISDM, please cite: \n\n"
        "Liu W, Srivastava A, Zhang J (2011) A Mathematical Framework "
        "for Protein Structure Comparison. PLoS Comput Biol 7(2): "
        "e1001075.\n\nAssumes only CA atoms." },
        { "-td", FALSE, etREAL, {&user_td},
        "Time difference used for -tdo output." },
    };
    
    
    
    FILE       *out;
    t_trxstatus *status;
    t_topology top;
    int        ePBC;
    rvec       *x, *iframe, *jframe, rrot_xyz, xold;
    rvec       **frames;
    real       *nweights, *iweights;
    real       ISD, sumISD, avgISD, maxISD, decorri, sumdn;
    real       *diff, *decorr, *decorr2, *snr, *td;
    matrix     box, rrot, rrotx, rroty, rrotz;
    real       t, t1, t2, dt, dm, dn, rgi, rgj, pi = 3.14159265358979;
    int        *rnum, *nsum, noptions, pcalcs, nsums;
    int        i, j, k, m, n, ij, mn, df, iatoms, natoms, nframes, nf2, nf4;
    gmx_bool   bDFLT, bFit;
    gmx_bool   bDecorr, bdDecorr, bScaled, bMScaled, bDouble, bSNR, bTD;
    char       *ISDM, *grpname, title[256], title2[256], *rname;
    atom_id    *index;
    output_env_t oenv;
    gmx_rmpbc_t  gpbc=NULL;
    const char *leg[]  = { "Decorrelation", "Overall Mean" }; 
    #define NLEG asize(leg) 
    t_filenm fnm[] = {
        { efTRX, "-f",       NULL,      ffREAD },
        { efTPS, NULL,       NULL,      ffREAD },
        { efNDX, NULL,       NULL,      ffOPTRD },
        { efXVG, "-decorr",  "decorr",  ffOPTWR },
        { efXVG, "-ddecorr", "ddecorr", ffOPTWR },
        { efXVG, "-scaled",  "scaled",  ffOPTWR },
        { efXVG, "-mscaled", "mscaled", ffOPTWR },
        { efXVG, "-snr",     "snr",     ffOPTWR },
        { efXVG, "-tdo",     "td",      ffOPTWR },
        { efXVG, "-double",  "double",  ffOPTWR }, 
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
              bSE2E || bMIR || bRROT || bSDRMS);
    
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
    
    // Output which graphs?
    bDecorr  = opt2bSet("-decorr",  NFILE, fnm);
    bdDecorr = opt2bSet("-ddecorr", NFILE, fnm);
    bScaled  = opt2bSet("-scaled",  NFILE, fnm);
    bMScaled = opt2bSet("-mscaled", NFILE, fnm);
    bSNR     = opt2bSet("-snr",     NFILE, fnm);
    bTD      = opt2bSet("-tdo",     NFILE, fnm);
    bDouble  = opt2bSet("-double",  NFILE, fnm);
    
    
    /* Opens trj. Reads first frame. Returns status. Allocates mem for x.
     * 
     * Not sure which argument determines which atoms to pull info for.
     */
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
    
    
    // Makes an array of weights. Necessary for reset_x.
    for (i = 0; i < iatoms; i++)
    {
        // Give a value for the weights.
        nweights[(int)index[i]] = 1;
        iweights[i] = 1;
    }
    
    
    nframes = 0;
    t  = 0;
    t2 = 0;
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
    // Find time interval between frames, assume evenly spaced.
    dt = t2 - t1;
    // Half of nframes.
    nf2 = nframes / 2;
    // Arrays to store decorrelation information.
    snew(decorr, (nf2 + 1));
    // Create an array to hold all frames.
    snew(frames, nframes);
    
    
    
    /* Opens trj. Reads first frame. Returns status. Allocates mem for x.
     * 
     * Not sure which argument determines which atoms to pull info for.
     */
    fprintf(stderr,"\nStoring trajectory to memory.\n");
    natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
    
    // Initialize index to keep track of current frame.
    i = 0;
    // This is for removing periodic boundary conditions.
    gpbc = gmx_rmpbc_init(&top.idef,ePBC,natoms,box);
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
    
    
    if (bSNR)
    {
        /* Main calculation loop.
         */
        fprintf(stderr,"\nCalculating results for SNR. \n");
        
        // Allocate memory for calculations.
        snew(snr,nframes);
        
        // Percentage of calculations complete.
        pcalcs = 1;
        // Number of calculations complete so far.
        int ncalcs = 0;
        
        // Different type of loop.
        for (ij = 1; ij <= nf2; ij++)
        {
            sumISD = 0.0;
            // Loop through time steps.
            for (i=0; i < (nframes - ij); i++)
            {
                j = i + ij;
                
                /* In this section, we'll put calls to all of the ISDMs.
                 * 
                 * Each should have its own if statement, so it is only executed
                 * if that option is specified at the command line.
                 */
                // Copy ith frame.
                if (bRROT)
                {
                    // Make a copy of the fit frame.
                    copy_rvecn(frames[i], iframe, 0, iatoms);
                }
                else
                {
                    iframe = frames[i];
                }
                
                // Fit the jth frame.
                if (bFit)
                {
                    // Need to make a copy of the fit frame or bad stuff will happen.
                    copy_rvecn(frames[j], jframe, 0, iatoms);
                    // Aligns jframe to current reference frame.
                    do_fit(iatoms, iweights, frames[i], jframe);
                }
                else
                {
                    jframe = frames[j];
                }
                
                // Calls most ISDM options.
                if (bDFLT || bRMSD || bSRMS || bRG || bSRG || bE2E || bSE2E || 
                    bMIR || bANG || bDIH || bANGDIH || bPHIPSI || bDRMS || 
                    bSDRMS || bPCOR || bACOR)
                {
                    ISD = call_ISDM(iatoms, jframe, iframe, diff, ISDM);
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
                            xold[m] = iframe[k][m];
                        }
                        for (m = 0; m < 3; m++)
                        {
                            iframe[k][m] = 0;
                            for (n = 0; n < 3; n++)
                            {
                                iframe[k][m] += rrot[m][n] * xold[n];
                            }
                        }
                    }
                    // Calculate RMSD after rotation.
                    ISD = sqrt(calc_msd(iatoms, jframe, iframe, diff));
                }
                
                // MAMMOTH. User gives -mammoth option.
                if (bMAMMOTH)
                {
                    // Calculate MAMMOTH comparison.
                    ISD = calc_mammoth(iatoms, jframe, iframe, rnum);
                }
                
                // ESA.
                if (bESA)
                {
                    // Calculate ESA comparison.
                    ISD = calc_esa(iatoms, jframe, iframe);
                }
                
                // Store the difference.
                snr[i] = ISD;
                
                // Using this to get the average difference.
                sumISD += ISD;
            }
            
            // Calculate the average of ISD for this loop.
            avgISD = sumISD / (nframes - ij);
            
            // Now loop again to find the variance.
            sumISD = 0.0;  // This variable can be reused.
            for (i = 0; i < (nframes - ij); i++)
            {
                decorri = snr[i] - avgISD;
                decorri = decorri * decorri;
                sumISD += decorri;
            }
            sumISD /= (nframes - ij);
            decorr[ij] = avgISD / ((real)sqrt((double)sumISD));
            
            // Update progress output.
            ncalcs += (nframes - ij);
            while ((double)ncalcs / (nf2 * (nframes + nframes - nf2 - 1) / 2) >= (double)pcalcs / 100)
            {
                fprintf(stderr, "Approximately %i percent complete. \r", pcalcs);
                fflush(stderr);
                pcalcs++;
            }
        }
        
        // Necessary to output to xvg format. Sets up the header.
        out=xvgropen(opt2fn("-snr", NFILE, fnm), 
                     "Signal to Noise Ratio of Decorrelation", 
                     "Time Difference, Delta T (ns)", 
                        "Signal to Noise Ratio at Delta T", 
                        oenv);
        
        // Process decorrelation information.
        for (i = 1; i <= nf2; i++)
        {
            // Print output.
            fprintf(out,"%10f %10f \n", (dt * i) / 1000.0, decorr[i]);
            // Reset decorr in case it is needed for another output option.
            decorr[i] = 0;
        }
        
        // Close the output file.
        ffclose(out);
    }
    
    
    /* The output for the decorr and scaled options is calculated here.
     * 
     * The output for the double option does not occur here, but the data in 
     * the decorr array must be calculated. The same for the ddecorr option.
     */
    if (bDecorr || bScaled || bMScaled || bDouble || bdDecorr)
    {
        /* Main calculation loop.
         */
        fprintf(stderr,"\nCalculating results for decorrelation. \n");
        
        // Percentage of calculations complete.
        pcalcs = 1;
        
        // Initiate the average and maximum interstructure difference.
        avgISD = 0;
        maxISD = 0;
        
        // Loop through time steps.
        for (i=0; i < nframes; i++)
        {
            sumISD = 0.0;
            // Loop through windows.
            for (j=0; j < nframes; j++)
            {
                // Diffs should be zero when i == j.
                if (i == j)
                {
                    continue;
                }
                
                ij = j - i;
                /* In this section, we'll put calls to all of the ISDMs.
                 * 
                 * Each should have its own if statement, so it is only executed
                 * if that option is specified at the command line.
                 */
                // Skip for i == j (comparing structure with self).
                if (i == j)
                {
                    continue;
                }
                
                // Copy ith frame.
                if (bRROT)
                {
                    // Make a copy of the fit frame.
                    copy_rvecn(frames[i], iframe, 0, iatoms);
                }
                else
                {
                    iframe = frames[i];
                }
                
                // Fit the jth frame.
                if (bFit)
                {
                    // Need to make a copy of the fit frame or bad stuff will happen.
                    copy_rvecn(frames[j], jframe, 0, iatoms);
                    // Aligns jframe to current reference frame.
                    do_fit(iatoms, iweights, frames[i], jframe);
                }
                else
                {
                    jframe = frames[j];
                }
                
                // Calls most ISDM options.
                if (bDFLT || bRMSD || bSRMS || bRG || bSRG || bE2E || bSE2E || 
                    bMIR || bANG || bDIH || bANGDIH || bPHIPSI || bDRMS || 
                    bSDRMS || bPCOR || bACOR)
                {
                    ISD = call_ISDM(iatoms, jframe, iframe, diff, ISDM);
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
                            xold[m] = iframe[k][m];
                        }
                        for (m = 0; m < 3; m++)
                        {
                            iframe[k][m] = 0;
                            for (n = 0; n < 3; n++)
                            {
                                iframe[k][m] += rrot[m][n] * xold[n];
                            }
                        }
                    }
                    // Calculate RMSD after rotation.
                    ISD = sqrt(calc_msd(iatoms, jframe, iframe, diff));
                }
                
                // MAMMOTH. User gives -mammoth option.
                if (bMAMMOTH)
                {
                    // Calculate MAMMOTH comparison.
                    ISD = calc_mammoth(iatoms, jframe, iframe, rnum);
                }
                
                // ESA.
                if (bESA)
                {
                    // Calculate ESA comparison.
                    ISD = calc_esa(iatoms, jframe, iframe);
                }
                
                // Put the differences into the algorithm.
                if ((ij > 0) && (ij <= nf2))
                {
                    // This stores the decorrelation information.
                    decorr[ij] += ISD;
                }
                // Using this to get the average difference.
                sumISD += ISD;
                // Using this to get the maximum difference.
                if (ISD > maxISD)
                {
                    maxISD = ISD;
                }
            }
            
            // Update the average.
            sumISD /= (nframes - 1);
            avgISD += sumISD;
            
            // Update progress output.
            while ((double)(i+1)/nframes >= (double)pcalcs/100)
            {
                fprintf(stderr, "\rApproximately %i percent complete.", pcalcs);
                fflush(stderr);
                pcalcs++;
            }
        }
        // Solve for the average difference.
        avgISD /= nframes;
        
        // Divide by the number of differences summed.
        for (i = 1; i <= nf2; i++)
        {
            decorr[i] = decorr[i] / (nframes - i);
        }
        
        // Special case for rrot.
        if (bRROT)
        {
            for (i = 0; i < nframes; i++)
            {
                j = i;
                // Make a copy of the fit frame.
                copy_rvecn(frames[i], iframe, 0, iatoms);
                jframe = frames[j];
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
                        xold[m] = iframe[k][m];
                    }
                    for (m = 0; m < 3; m++)
                    {
                        iframe[k][m] = 0;
                        for (n = 0; n < 3; n++)
                        {
                            iframe[k][m] += rrot[m][n] * xold[n];
                        }
                    }
                }
                // Calculate RMSD after rotation.
                ISD = sqrt(calc_msd(iatoms, jframe, iframe, diff));
                decorr[0] = ISD;
            }
        }
        
        // Output decorrelation.
        if (bDecorr)
        {
            // Necessary to output to xvg format. Sets up the header.
            out=xvgropen(opt2fn("-decorr", NFILE, fnm), 
                         "Simple Decorrelation", 
                         "Time Difference, Delta T (ns)", 
                         "Mean Difference at Delta T", 
                         oenv);
            
            // Sets up the legend for the xvg file.
            //xvgr_legend(out,NLEG,leg,oenv);
            
            // Process decorrelation information.
            for (i = 0; i <= nf2; i++)
            {
                // Print output.
                fprintf(out,"%10f %10f \n", 
                        (dt * i) / 1000.0, 
                        decorr[i]);
            }
            // Close the output file.
            ffclose(out);
        }
        
        // Output decorrelation scaled by overall mean difference.
        if (bScaled)
        {
            // Necessary to output to xvg format. Sets up the header.
            out=xvgropen(opt2fn("-scaled", NFILE, fnm), 
                         "Simple Decorrelation Scaled by Mean Difference", 
                         "Time Difference, Delta T (ns)", 
                         "Difference at Delta T / Mean Difference", 
                         oenv);
            
            // Sets up the legend for the xvg file.
            xvgr_legend(out,NLEG,leg,oenv);
            
            // Print scaling.
            printf("\n\nAverage difference is %10f ",avgISD);
            
            // Process decorrelation information.
            for (i = 0; i <= nf2; i++)
            {
                // Print output.
                fprintf(out, "%10f %10f %10f \n", 
                        (dt * i) / 1000.0, 
                        (decorr[i] / avgISD), 1.0);
            }
            // Close the output file.
            ffclose(out);
        }
        
        // Output decorrelation scaled by overall maximum difference.
        if (bMScaled)
        {
            // Necessary to output to xvg format. Sets up the header.
            out=xvgropen(opt2fn("-mscaled", NFILE, fnm), 
                         "Simple Decorrelation Scaled by Max Difference", 
                         "Time Difference, Delta T (ns)", 
                         "Difference at Delta T / Maximum Difference", 
                         oenv);
            
            // Sets up the legend for the xvg file.
            //xvgr_legend(out,NLEG,leg,oenv);
            
            // Print scaling.
            printf("\n\nMaximum difference is %10f ",maxISD);
            
            // Process decorrelation information.
            for (i = 0; i <= nf2; i++)
            {
                // Print output.
                fprintf(out, "%10f %10f \n", 
                        (dt * i) / 1000.0, 
                        (decorr[i] / maxISD));
            }
            // Close the output file.
            ffclose(out);
        }
    }
    
    
    if (bTD)
    {
        // Error checking.
        if (user_td == -1)
        {
            gmx_fatal(FARGS,"\nThe -tdo output options requires time -td to be set.\n");
        }
        if (user_td < dt)
        {
            gmx_fatal(FARGS,"\nThe time -td must be greater than trx frame intervals.\n");
        }
        if (user_td > (dt * nframes))
        {
            gmx_fatal(FARGS,"\nThe time -td is greater than the length of the trajectory.\n");
        }
        
        // Calculate the ij to be used.
        ij = (int)(user_td / dt);
        fprintf(stderr,"\nUsing %10f ps as -td to align with trx frames.\n",(dt * ij));
        // Allocate memory to store results.
        snew(td,(nframes - ij));
        
        /* Main calculation loop.
         */
        fprintf(stderr,"\n\nCalculating results for TD. \n");
        
        // Percentage of calculations complete.
        pcalcs = 1;
        
        // Loop through time steps.
        for (i=0; i < (nframes - ij); i++)
        {
            j = i + ij;
            
            /* In this section, we'll put calls to all of the ISDMs.
             * 
             * Each should have its own if statement, so it is only executed
             * if that option is specified at the command line.
             */
            // Skip for i == j (comparing structure with self).
            if (i == j)
            {
                continue;
            }
            
            // Copy ith frame.
            if (bRROT)
            {
                // Make a copy of the fit frame.
                copy_rvecn(frames[i], iframe, 0, iatoms);
            }
            else
            {
                iframe = frames[i];
            }
                
            // Fit the jth frame.
            if (bFit)
            {
                // Need to make a copy of the fit frame or bad stuff will happen.
                copy_rvecn(frames[j], jframe, 0, iatoms);
                // Aligns jframe to current reference frame.
                do_fit(iatoms, iweights, frames[i], jframe);
            }
            else
            {
                jframe = frames[j];
            }
            
            // Calls most ISDM options.
            if (bDFLT || bRMSD || bSRMS || bRG || bSRG || bE2E || bSE2E || 
                bMIR || bANG || bDIH || bANGDIH || bPHIPSI || bDRMS || 
                bSDRMS || bPCOR || bACOR)
            {
                ISD = call_ISDM(iatoms, jframe, iframe, diff, ISDM);
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
                        xold[m] = iframe[k][m];
                    }
                    for (m = 0; m < 3; m++)
                    {
                        iframe[k][m] = 0;
                        for (n = 0; n < 3; n++)
                        {
                            iframe[k][m] += rrot[m][n] * xold[n];
                        }
                    }
                }
                // Calculate RMSD after rotation.
                ISD = sqrt(calc_msd(iatoms, jframe, iframe, diff));
            }
            
            // MAMMOTH. User gives -mammoth option.
            if (bMAMMOTH)
            {
                // Calculate MAMMOTH comparison.
                ISD = calc_mammoth(iatoms, jframe, iframe, rnum);
            }
            
            // ESA.
            if (bESA)
            {
                // Calculate ESA comparison.
                ISD = calc_esa(iatoms, jframe, iframe);
            }
            
            // Store the difference.
            td[i] = ISD;
            
            // Update progress output.
            while ((double)i / (nframes - ij) >= (double)pcalcs/100)
            {
                fprintf(stderr, "\rApproximately %i percent complete.", pcalcs);
                fflush(stderr);
                pcalcs++;
            }
        }
        
        
        
        // Necessary to output to xvg format. Sets up the header.
        out=xvgropen(opt2fn("-tdo", NFILE, fnm), 
                     "Difference at Time -td", 
                     "Time (ns)", 
                     "Difference", 
                     oenv);
        
        // Output information to file.
        for (i = 0; i < (nframes - ij); i++)
        {
            // Print output.
            fprintf(out,"%10f %10f \n", (dt * i) / 1000.0, td[i]);
        }
        // Close the output file.
        ffclose(out);
    }
    
    
    // Need to allocate memory for both output types.
    if (bDouble || bdDecorr)
    {
        nf4 = nf2 / 2;
        snew(decorr2, (nf4 + 1));
    }
    
    /* The purpose of this section is to compare the difference between two 
     * frames separated by ij frames near the origin of comparison (frame 0) 
     * with the difference between two frames separated by ij frames near the 
     * end of the decorrelation data (frame nf2 = nframes / 2). I limit the 
     * value of ij to nf4 = nf2 / 2 (might not be the same as nframes / 4).
     * 
     * The m and n variables here refers to positions in the decorrelation 
     * data. They do not refer to the same thing as i and j in previous loops. 
     * The sensitivity should decrease as the distance from the reference 
     * increases. This attempts to quantify the amount of decrease.
     */
    if (bDouble)
    {
        fprintf(stderr, "\n\nSensitivity calculations...\n");
        sumdn = 0;
        snew(nsum, (nf4 + 1));
        for (df = 1; df <= nf4; df++)
        {
            for (m = 0; m < (nf2 - df - 1); m++)
            {
                for (n = (m + 1); n < (nf2 - df); n++)
                {
                    mn = n - m;
                    
                    if (mn <= nf4)
                    {
                        // Test for discontinuities.
                        dm = decorr[m + df] - decorr[m];
                        if (dm == 0)
                        {
                            continue;
                        }
                        // Average of dn is a correction factor.
                        dn = decorr[n + df] - decorr[n];
                        sumdn += dn;
                        decorri = dn / dm;
                        // Increment the number of sums.
                        nsum[mn]++;
                        
                        /* Test for outliers. This works on the assumption 
                         * that the slope of the decorrelation is always 
                         * decreasing. It may be too strict or not strict 
                         * enough.
                         */
                        if (decorri > 1)
                        {
                            decorr2[mn] +=  1;
                        }
                        else if (decorri < -1)
                        {
                            decorr2[mn] += -1;
                        }
                        else
                        {
                            decorr2[mn] += decorri;
                        }
                    }
                }
            }
        }
        
        
        // Necessary to output to xvg format. Sets up the header.
        out=xvgropen(opt2fn("-double", NFILE, fnm), 
                     "Second Decorrelation of Decorrelation Data", 
                     "Time From Reference Frame (ns)", 
                "Reduction of Sensitivity", 
                oenv);
        
        // Print decorrelation information at 0 time interval.
        fprintf(out,"%10f %10f \n", 0.0, 1.0);
        
        // Print output for the -double option.
        nsums = 0;
        for (i = 1; i <= nf4; i++)
        {
            // Divide by the number of differences summed.
            decorr2[i] = decorr2[i] / nsum[i];
            nsums += nsum[i];
            
            // Print output.
            fprintf(out,"%10f %10f \n", 
                    (dt * i) / 1000.0, 
                    decorr2[i]);
            
            // Reset decorr2 in case it is used again.
            decorr2[i] = 0;
        }
        printf("\n\nDouble decorrelation correction factors: %10f %10f \n", 
               (sumdn / nsums), 
               (sumdn / nsums / avgISD));
    }
    
    
    /* The derivative of the decorrelation data as it decreases.
     */
    if (bdDecorr)
    {
        fprintf(stderr, "\n\nSensitivity ddecorr calculations...\n");
        for (m = 1; m < (nf2 - 2); m++)
        {
            for (n = (m + 1); n < (nf2 - 1); n++)
            {
                mn = n - m;
                
                if (mn < nf4)
                {
                    // Check for discontinuities.
                    dm = decorr[m + 1] - decorr[m - 1];
                    if (dm == 0)
                    {
                        continue;
                    }
                    /* Full equation to take the ratio of derivatives is:
                     * 
                     * (decorr[n+1] - decorr[n-1] / 2*dt) / 
                     * (decorr[m+1] - decorr[m-1] / 2*dt)
                     * 
                     * The 2*dt portions cancel.
                     */
                    decorri = (decorr[n + 1] - decorr[n - 1]) / dm;
                    
                    /* Test for outliers. This works under the assumption that 
                     * the slope should always be decreasing. It may be either 
                     * too strict or not strict enough.
                     */
                    if (decorri > 1)
                    {
                        decorr2[mn] +=  1;
                    }
                    else if (decorri < -1)
                    {
                        decorr2[mn] += -1;
                    }
                    else
                    {
                        decorr2[mn] += decorri;
                    }
                }
            }
        }
        
        // Necessary to output to xvg format. Sets up the header.
        out=xvgropen(opt2fn("-ddecorr", NFILE, fnm), 
                     "Ratio of Derivatives of Decorrelation Data", 
                     "Time From Reference Frame (ns)", 
                     "Reduction of Sensitivity", 
                     oenv);
        
        // Print decorrelation information at 0 time interval.
        fprintf(out,"%10f %10f \n", 0.0, 1.0);
        
        // Print output for the -double option.
        for (i = 1; i < nf4; i++)
        {
            // Divide by the number of differences summed.
            decorr2[i] = decorr2[i] / (nf2 - i - 2);
            
            // Print output.
            fprintf(out,"%10f %10f \n", 
                    (dt * i) / 1000.0, 
                    decorr2[i]);
            
            // Reset decorr2 in case it is used again.
            decorr2[i] = 0;
        }
    }
    
    // Not sure what this one does. Sends to xmgrace?
    //do_view(oenv,ftp2fn(efXVG,NFILE,fnm),"-nxy");
    
    thanx(stderr);
    return 0;
}