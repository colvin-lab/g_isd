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



void mat_mult_mat(real* mat1, real* mat2, int m, int n, int o, real* out)
{
    /* Assume a is an array of doubles m by n in dimensions.
     * Assume b is an array of doubles n by o in dimensions.
     * Out should point to enough memory to store m by o doubles.
     */
    int i, j, k;
    for (i = 0; i < m; i++) {
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
    /* Transpose m by n array of doubles a into n by m array of doubles out.
     */
    int i, j;
    
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            out[(j * n) + i] = mat[(i * m) + j];
        }
    }
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
    static gmx_bool bANG2=FALSE, bDIH2=FALSE, bANGDIH2=FALSE;
    static gmx_bool bRROT=FALSE, bSDRMS=FALSE;
    static real setmax = -1.0;
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
        { "-setmax", FALSE, etREAL, {&setmax},
            "Set maximum value to threshold the xpm file. Must be greater "
            "than the average inter-structure distance." },
    };
    
    
    
    FILE       *out;
    t_trxstatus *status;
    t_topology top;
    int        ePBC;
    rvec       *x, *iframe, *jframe, *rframe, *cframe, rrot_xyz, xold;
    rvec       **frames;
    real       *nweights, *iweights, abscoor, maxcoor;
    real       *diff, ISD, **ISDmat, *P2, *J, *P2J, *B, *BT, *E, *V, *MDS;
    real       **Va, **MDSa;
    double     dISD, *avgdiff, avgISD, maxISD;
    matrix     box, rrot, rrotx, rroty, rrotz;
    real       t, xpm_max, pi = 3.14159265358979;
    int        *maxframe, *rnum, maxcoori;
    int        i, j, k, m, n, p, np, iatoms, natoms, nframes;
    int        percentcalcs, noptions;
    gmx_bool   bDFLT, bFit, bISD, bMDS, bEig, bVec;
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
        { efDAT, "-vec", "eigvecs",  ffOPTWR },
        { efDAT, "-isd", "isdcsv",   ffOPTWR },
        { efDAT, "-mds", "mdscsv",   ffOPTWR },
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
    read_tps_conf(ftp2fn(efTPS, NFILE, fnm), title, &top, &ePBC, &x, NULL, box, TRUE);
    
    // Asks you to choose a selection of atoms at prompt.
    get_index(&top.atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &iatoms, &index, &grpname);
    
    
    // If there are no options at command line, do default behavior.
    bDFLT = !(bANG || bDIH || bANGDIH || bPHIPSI || bDRMS || bSRMS || bRMSD || 
              bPCOR || bACOR || bMAMMOTH || bESA || bRG || bSRG || bE2E || 
              bSE2E || bMIR || bRROT || bSDRMS || bANG2 || bDIH2 || bANGDIH2);
    
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
    
    if (bANGDIH2)
    {
        fprintf(stderr,"\nUsing geometric mean of angles and dihedrals as ISDM.\n");
        ISDM = "ANGDIH2";
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
    if (bRROT)
    {
        snew(iframe,iatoms);
    }
    if (bFit)
    {
        snew(jframe,iatoms);
    }
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
    bVec = opt2bSet("-vec", NFILE, fnm);
    bISD = opt2bSet("-isd", NFILE, fnm);
    bMDS = opt2bSet("-mds", NFILE, fnm);
    
    
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
    snew(avgdiff, nframes);
    snew(ISDmat,  nframes);
    for (i=0; i<nframes; i++)
    {
        avgdiff[i] = 0;
        snew(ISDmat[i],nframes);
    }
    
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
    percentcalcs = 1;
    
    // Loop through reference frames.
    for (i = 0; i < nframes; i++)
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
                bSDRMS || bPCOR || bACOR || bANG2 || bDIH2 || bANGDIH2)
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
            if (dISD > maxISD)
            {
                maxISD = dISD;
            }
            avgdiff[i] += dISD;
            
            // Debugging.
            //printf("On the %i th loop. \n",j);
        }
        
        // Average difference for each frame.
        avgdiff[i] /= (nframes - 1);
        
        // Update progress output.
        while ((double)(i+1)/nframes >= (double)percentcalcs/100)
        {
            fprintf(stderr, "Approximately %i percent complete. \r", percentcalcs);
            fflush(stderr);
            percentcalcs++;
        }
    }
    fprintf(stderr, "\n\n");
    
    // Find the final average of differences.
    for (i=0; i<nframes; i++)
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
    mat_mult_mat(P2,   J, nframes, nframes, nframes, P2J);
    B = P2; // Finished with the memory in P2. Reuse it to store B.
    scl_mult_mat(-0.5, J, nframes, nframes, J);
    mat_mult_mat(J,  P2J, nframes, nframes, nframes, B);
    
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
    snew(MDS,  nframes * np);
    snew(MDSa, nframes);
    for (i = 0; i < nframes; i++)
    {
        MDSa[i] = &MDS[np * i];
    }
    for (i = 0; i < nframes; i++)
    {
        for (j = 0; j < np; j++)
        {
            MDSa[i][j] = sqrt(E[nframes - j - 1]) * Va[nframes - j - 1][i];
        }
    }
    
    // Step 5.
    fprintf(stderr, "MDS step 5 of 5. \r");
    for (j = 0; j < np; j++)
    {
        maxcoor = -1.0;
        for (i = 0; i < nframes; i++)
        {
            abscoor = abs(MDSa[i][j]);
            if (abscoor > maxcoor)
            {
                maxcoor  = abscoor;
                maxcoori = i;
            }
        }
        
        if (MDSa[maxcoori][j] < 0.0)
        {
            for (i = 0; i < nframes; i++)
            {
                MDSa[i][j] *= -1.0;
            }
        }
    }
    fprintf(stderr, "\nMDS Complete. \n\n");
    
    // Output dimensionally reduced coordinates.
    if (bMDS)
    {
        // Opens the output file.
        out = opt2FILE("-mds", NFILE, fnm, "w");
        
        // Write output.
        for (i = 0; i < nframes; i++)
        {
            fprintf(out, "%15.6e", MDSa[i][0]);
            for (j = 1; j < np; j++)
            {
                fprintf(out, ",%15.6e", MDSa[i][j]);
            }
            fprintf(out, "\n");
        }
        
        // Close the output file.
        ffclose(out);
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
    
    // Closing.
    thanx(stderr);
    return 0;
}
