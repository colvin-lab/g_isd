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




int gmx_isdorder(int argc,char *argv[])
{
    const char *desc[] = {
        "[TT]g_isdorder[tt] implements a list of measures designed to ",
        "compare two frames from a trajectory and quantify their difference. ",
        "This tool analyzes the input by comparing each pair of frames in ",
        "the trajectory. The quantified ",
        "difference between two structures is referred to as the inter-",
        "structure distance (ISD), and the methods used to measure the ",
        "ISD are referred to as ISDMs. ",
        "Every frame of the trajectory is compared to every other frame, and ",
        "the resulting differences are summed in an equation which scales ",
        "output between 0 (ordered) and 1 (disordered). Depending on the ",
        "chosen ISDM, this tool can approximate the disorder ",
        "of each atom. If the chosen ",
        "ISDM does not have this functionality, the xvg file will show -1 ",
        "for every atom. The default metric if one is not chosen by the ",
        "user is RMSD. Only one metric can be chosen at a time."
    };
    
    
    
    static gmx_bool bANG=FALSE, bDIH=FALSE, bANGDIH=FALSE, bDRMS=FALSE;
    static gmx_bool bPHIPSI=FALSE, bSRMS=FALSE, bPCOR=FALSE, bMAMMOTH=FALSE;
    static gmx_bool bACOR=FALSE, bESA=FALSE, bRMSD=FALSE, bMIR=FALSE;
    static gmx_bool bRG=FALSE, bSRG=FALSE, bE2E=FALSE, bANGDIH2G=FALSE;
    static gmx_bool bANG2=FALSE, bDIH2=FALSE, bANGDIH2=FALSE, bSE2E=FALSE;
    static gmx_bool bRROT=FALSE, bSDRMS=FALSE, bPHIPSI2=FALSE;
    static gmx_bool user_linear=FALSE, user_fisherstultz=FALSE;
    static gmx_bool bLinear=FALSE, bFisherStultz=FALSE, bOrder=FALSE;
    static real     user_sf = -1.0, user_zero = -1.0, user_one = -1.0;
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
            "e1001075.\n\nAssumes only CA atoms." },
        { "-linear", FALSE, etBOOL, {&user_linear},
            "Overrides the default order parameter equation to use linear "
            "rescaling. The rescaling limits can be set manually with the "
            "-zero and -one options respectively. No change in behavior if "
            "-linear is already the default rescaling equation for the chosen "
            "ISDM. The -linear behavior is the default for ISDM options which "
            "do not have an infinite upper bound:\n\n"
            "-ang, -dih, -angdih, -ang2, -dih2, -angdih2, -angdih2g, -phipsi, "
            "-phipsi2, -sdrms, -srms, -pcor, -acor, -esa\n" },
        { "-fisherstultz", FALSE, etBOOL, {&user_fisherstultz},
            "Overrides the default order parameter equation to use the "
            "Fisher-Stultz equation rescaling. The scaling factor can be "
            "set manually with the -sf option. No change in behavior if "
            "-fisherstultz is already the default rescaling equation for the "
            "chosen ISDM. The -fisherstultz behavior is the default for ISDM "
            "options which have an infinite upper bound:\n\n"
            "-drms, -rmsd\n" },
        { "-sf", FALSE, etREAL, {&user_sf},
            "Overrides the default scaling factor used by the Fisher-Stultz "
            "order equation. Ignored unless the -fisherstultz option is "
            "the default for the chosen ISDM or the -fisherstultz flag is "
            "given manually." },
        { "-zero", FALSE, etREAL, {&user_zero},
            "Overrides the default rescaling limits used by the linear "
            "rescaling option. Sets the value which is rescaled to zero. "
            "Ignored unless the -linear option is the default for the chosen "
            "ISDM or the -linear flag is given manually." },
        { "-one", FALSE, etREAL, {&user_one},
            "Overrides the default rescaling limits used by the linear "
            "rescaling option. Sets the value which is rescaled to one. "
            "Ignored unless the -linear option is the default for the chosen "
            "ISDM or the -linear flag is given manually." },
        { "-order", FALSE, etBOOL, {&bOrder},
            "By default, the order parameter for the entire molecule is "
            "calculated if the -resid option is not specified. By default, "
            "the order parameter for the entire molecule is NOT calculated "
            "if the -resid option is specified. Using the -order option "
            "without the -resid option does not change behavior. Using the "
            "-order option with the -resid option causes both to be "
            "calculated with a corresponding increase of analysis time." },
    };
    
    
    
    FILE       *out;
    t_trxstatus *status;
    t_topology top;
    int        ePBC;
    rvec       *x, *iframe, *jframe, *rframe, *cframe, xold;
    rvec       **frames, *topframe, *trjframe;
    real       *nweights, *iweights;
    real       ISD, *rISD;
    double     dISD, *drISD, ISDorder, ISDsum, msqi;
    double     *rISDsum, *rISDorder;
    double     order_sf, order_zero, order_one;
    matrix     box;
    real       t, pi = 3.14159265358979;
    int        *maxframe, *rnum;
    int        i, j, k, m, n, iatoms, natoms, nframes, nres;
    int        pcalcs, noptions;
    gmx_bool   bDFLT, bResid, bFit;
    char       *ISDM, *grpname, title[256], title2[256], *rname;
    atom_id    *index;
    output_env_t oenv;
    gmx_rmpbc_t  gpbc=NULL;
    const char *leg[]  = { "D" }; 
#define NLEG asize(leg) 
    t_filenm fnm[] = {
        { efTRX,     "-f",    NULL,    ffREAD }, 
        { efTPS,     NULL,    NULL,    ffREAD },
        { efNDX,     NULL,    NULL,   ffOPTRD },
        { efXVG, "-resid", "resid",   ffOPTWR },
    }; 
#define NFILE asize(fnm)
    int npargs;
    
    CopyRight(stderr,argv[0]);
    npargs = asize(pa);
    
    // Lots of black magic with this one. The oenv is used by many things.
    parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
                      NFILE,fnm,npargs,pa,asize(desc),desc,0,NULL,&oenv);
    
    // Output which files?
    bResid = opt2bSet("-resid", NFILE, fnm);
    
    
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
              bSE2E || bMIR || bRROT || bSDRMS || bANG2 || bDIH2 || 
              bPHIPSI2 || bANGDIH2 || bANGDIH2G);
    
    bFit  =  (bDFLT || bRMSD || bMIR || bSRMS || bPCOR);
    
    if (bResid)
    {
        if (bPCOR || bACOR || bESA)
        {
            gmx_fatal(FARGS, "The -resid option does not work with this ISDM.");
        }
    }
    
    // For error checking.
    noptions = 0;
    // Check which ISDM will be used. Default is RMSD.
    if (bDFLT || bRMSD)
    {
        fprintf(stderr,"\nUsing RMSD as ISDM.\n");
        ISDM = "RMSD";
        noptions++;
        order_sf   = 2 * 0.275;
        order_zero = 0.0;
        order_one  = 1.0; // Nonsense.
        bFisherStultz = TRUE;
    }
    
    if (bANG)
    {
        fprintf(stderr,"\nUsing backbone angles as ISDM.\n");
        ISDM = "ANG";
        noptions++;
        order_sf   = 2 * 0.015;
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    if (bDIH)
    {
        fprintf(stderr,"\nUsing backbone dihedrals as ISDM.\n");
        ISDM = "DIH";
        noptions++;
        order_sf   = 2 * 0.15;
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    if (bANG2)
    {
        fprintf(stderr,"\nUsing backbone angles as ISDM.\n");
        ISDM = "ANG2";
        noptions++;
        order_sf   = 2 * 0.015; // Nonsense.
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    if (bDIH2)
    {
        fprintf(stderr,"\nUsing backbone dihedrals as ISDM.\n");
        ISDM = "DIH2";
        noptions++;
        order_sf   = 2 * 0.15; // Nonsense.
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    if (bPHIPSI)
    {
        fprintf(stderr,"\nUsing phi and psi angles as ISDM.\n");
        ISDM = "PHIPSI";
        noptions++;
        order_sf   = 2 * 0.15;
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    if (bPHIPSI2)
    {
        fprintf(stderr,"\nUsing phi and psi angles as ISDM.\n");
        ISDM = "PHIPSI2";
        noptions++;
        order_sf   = 2 * 0.15; // Nonsense.
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    if (bDRMS)
    {
        fprintf(stderr,"\nUsing distance RMS as ISDM.\n");
        ISDM = "DRMS";
        noptions++;
        order_sf   = 2 * 0.13; // Nonsense.
        order_zero = 0.0;
        order_one  = 1.0; // Nonsense.
        bFisherStultz = TRUE;
    }
    
    if (bSDRMS)
    {
        fprintf(stderr,"\nUsing scaled distance RMS as ISDM.\n");
        ISDM = "SDRMS";
        noptions++;
        order_sf   = 2 * 0.13;
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    if (bSRMS)
    {
        fprintf(stderr,"\nUsing scaled RMSD as ISDM.\n");
        ISDM = "SRMS";
        noptions++;
        order_sf   = 2 * 0.2; // Nonsense.
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    if (bPCOR)
    {
        fprintf(stderr,"\nUsing position correlation as ISDM.\n");
        ISDM = "PCOR";
        noptions++;
        order_sf   = 2 * 0.2; // Nonsense.
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    if (bACOR)
    {
        fprintf(stderr,"\nUsing backbone angle correlation as ISDM.\n");
        ISDM = "ACOR";
        noptions++;
        order_sf   = 2 * 0.2; // Nonsense.
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    if (bANGDIH)
    {
        fprintf(stderr,"\nUsing geometric mean of angles and dihedrals as ISDM.\n");
        ISDM = "ANGDIH";
        noptions++;
        order_sf   = 2 * 0.05;
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    if (bANGDIH2)
    {
        fprintf(stderr,"\nUsing root mean square of angles and dihedrals as ISDM.\n");
        ISDM = "ANGDIH2";
        noptions++;
        order_sf   = 2 * 0.05; // Nonsense.
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    if (bANGDIH2G)
    {
        fprintf(stderr,"\nUsing geometric mean of RMS angles and RMS dihedrals as ISDM.\n");
        ISDM = "ANGDIH2G";
        noptions++;
        order_sf   = 2 * 0.05; // Nonsense.
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
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
        order_sf   = 2 * 0.2; // Nonsense.
        order_zero = 0.0;
        order_one  = 1.0;
        bLinear = TRUE;
    }
    
    // Throw an error if multiple -ISDM options were given by the user.
    if (noptions > 1)
    {
        gmx_fatal(FARGS,"\nThis tool only supports using one optional ISDM at a time.\n");
    }
    
    // Throw an error if both -linear and -fisherstultz were given by user.
    if (user_linear && user_fisherstultz)
    {
        gmx_fatal(FARGS,"\nOnly one of -linear or -fisherstultz should be used.\n");
    }
    
    // Check if scaling equation is specified manually.
    if (user_linear)
    {
        bLinear       = TRUE;
        bFisherStultz = FALSE;
    }
    if (user_fisherstultz)
    {
        bLinear       = FALSE;
        bFisherStultz = TRUE;
    }
    
    // Check if scaling factor was specified manually.
    if (bFisherStultz)
    {
        if (user_sf != -1.0)
        {
            if (user_sf <= 0.0)
            {
                gmx_fatal(FARGS,"Scaling factor must be > 0, or we're doomed!!!");
            }
            
            fprintf(stderr, "\nScaling factor set to %f by user.\n", user_sf);
            order_sf = (double)(2 * user_sf);
        }
    }
    
    // Check if rescaling limits were specified manually.
    if (bLinear)
    {
        if ((user_zero != -1.0) || (user_one != -1.0))
        {
            if (user_zero == user_one)
            {
                gmx_fatal(FARGS,"Rescaling limits must not be equal, or we're doomed!!!");
            }
            
            if (user_zero != -1.0)
            {
                order_zero = (double)user_zero;
            }
            if (user_one  != -1.0)
            {
                order_one  = (double)user_one;
            }
            
            fprintf(stderr, "Rescaling limits set to zero: %f and one: %f by user.\n", order_zero, order_one);
        }
    }
    
    // Opens trj. Reads first frame. Returns status. Allocates mem for x.
    fprintf(stderr,"\nCounting the number of frames.\n");
    natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
    
    // Now that we have iatoms, allocate memory for other arrays.
    if (bFit)
    {
        snew(jframe,iatoms);
    }
    if (bResid)
    {
        if (bPHIPSI || bPHIPSI2)
        {
            nres = iatoms / 3;
        }
        else
        {
            nres = iatoms;
        }
        snew(rISD,      nres);
        snew(drISD,     nres);
        snew(rISDsum,   nres);
        snew(rISDorder, nres);
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
    
      
    
    // Create an array to hold all frames.
    snew(frames,nframes);
    
    
    // Load first frame to x variable.
    natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
    
    // This is for removing periodic boundary conditions.
    gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms, box);
    
    /* Opens trj. Reads first frame. Returns status. Allocates mem for x.
     * 
     * Not sure which argument determines which atoms to pull info for.
     */
    fprintf(stderr,"\nStoring trajectory to memory.\n");
    // Initialize index to keep track of current frame.
    i = 0;
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
    ISDorder = 0.0;
    
    
    
    /* Main calculation loop.
     */
    fprintf(stderr,"\nCalculating results. \n");
    pcalcs = 1;
    
    /* Originally this was designed to only loop through each pair of i and j
     * one time to save half of the calculations. Eventually it became 
     * impractical to make sure that each ISDM was symmetric, so now the
     * algorithm takes the performance hit in favor of accuracy and simplicity.
     */
    
    // Loop through reference frames.
    for (i = 0; i < nframes; i++)
    {
        ISDsum = 0.0;
        if (bResid)
        {
            for (k = 0; k < nres; k++)
            {
                rISDsum[k] = 0.0;
            }
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
                if (bResid)
                {
                    ISD = call_ISDM_n(iatoms, cframe, rframe, rISD, ISDM);
                    if (bOrder)
                    {
                        ISD = call_ISDM(iatoms, cframe, rframe, ISDM);
                    }
                }
                else
                {
                    ISD = call_ISDM(iatoms, cframe, rframe, ISDM);
                }
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
            
            
            // Calculate summations with doubles for improved accuracy.
            dISD = (double)ISD;
            if (bResid)
            {
                for (k = 0; k < nres; k++)
                {
                    drISD[k] = (double)rISD[k];
                }
            }
            
            // Update sums.
            if (bLinear)
            {
                if (bResid)
                {
                    if (bOrder)
                    {
                        ISDsum += dISD;
                    }
                    for (k = 0; k < nres; k++)
                    {
                        rISDsum[k] += drISD[k];
                    }
                }
                else
                {
                    ISDsum += dISD;
                }
            }
            if (bFisherStultz)
            {
                if (bResid)
                {
                    if (bOrder)
                    {
                        ISDsum += exp(-dISD / order_sf);
                    }
                    for (k = 0; k < nres; k++)
                    {
                        rISDsum[k] += exp(-drISD[k] / order_sf);
                    }
                }
                else
                {
                    ISDsum += exp(-dISD / order_sf);
                }
            }
        }
        
        if (bLinear)
        {
            if (bResid)
            {
                if (bOrder)
                {
                    ISDsum /= (nframes - 1);
                    ISDorder     += ISDsum;
                }
                for (k = 0; k < nres; k++)
                {
                    rISDsum[k]   /= (nframes - 1);
                    rISDorder[k] += rISDsum[k];
                }
            }
            else
            {
                ISDsum /= (nframes - 1);
            }
        }
        
        if (bFisherStultz)
        {
            if (bResid)
            {
                if (bOrder)
                {
                    ISDsum   /= (nframes - 1);
                    ISDorder += log(1 + ISDsum) / log(2);
                }
                for (k = 0; k < nres; k++)
                {
                    rISDsum[k]   /= (nframes - 1);
                    rISDorder[k] += log(1 + rISDsum[k]) / log(2);
                }
            }
            else
            {
                ISDsum   /= (nframes - 1);
                ISDorder += log(1 + ISDsum) / log(2);
            }
        }
        
        
        // Update progress output.
        while ((double)i / nframes >= (double)pcalcs/100)
        {
            fprintf(stderr,"\rApproximately %i percent complete.", pcalcs);
            fflush(stderr);
            pcalcs++;
        }
    }
    
    // Final rescaling steps.
    if (bLinear)
    {
        if (bResid)
        {
            if (bOrder)
            {
                ISDorder /= nframes;
                if (ISDorder < order_zero)
                {
                    ISDorder = order_zero;
                }
                if (ISDorder > order_one)
                {
                    ISDorder = order_one;
                }
                ISDorder  = (ISDorder - order_zero) / (order_one - order_zero);
            }
            for (k = 0; k < nres; k++)
            {
                rISDorder[k] /= nframes;
                if (rISDorder[k] < order_zero)
                {
                    rISDorder[k] = order_zero;
                }
                if (rISDorder[k] > order_one)
                {
                    rISDorder[k] = order_one;
                }
                rISDorder[k] = (rISDorder[k] - order_zero) / (order_one - order_zero);
            }
        }
        else
        {
            ISDorder /= nframes;
            ISDorder  = (ISDorder - order_zero) / (order_one - order_zero);
        }
    }
    
    
    if (bFisherStultz)
    {
        if (bResid)
        {
            if (bOrder)
            {
                ISDorder /= nframes;
                if (ISDorder < order_zero)
                {
                    ISDorder = order_zero;
                }
                if (ISDorder > order_one)
                {
                    ISDorder = order_one;
                }
                ISDorder = (ISDorder - order_zero) / (order_one - order_zero);
            }
            for (k = 0; k < nres; k++)
            {
                rISDorder[k] /= nframes;
                if (rISDorder[k] < order_zero)
                {
                    rISDorder[k] = order_zero;
                }
                if (rISDorder[k] > order_one)
                {
                    rISDorder[k] = order_one;
                }
                rISDorder[k] = (rISDorder[k] - order_zero) / (order_one - order_zero);
            }
        }
        else
        {
            ISDorder /= nframes;
            ISDorder  = (ISDorder - order_zero) / (order_one - order_zero);
        }
    }
    
    // Print output to specified files.
    if (bResid)
    {
        out=xvgropen(opt2fn("-resid", NFILE, fnm), 
                     "Disorder Parameter", 
                     "Residue", 
                     "Disorder", 
                     oenv);
        
        for (k = 0; k < nres; k++)
        {
            fprintf(out,"%-6i %10f \n", k + 1, rISDorder[k]);
        }
        ffclose(out);
    }
    
    
    // Output the final information.
    if (bLinear)
    {
        printf("\n\nOrder parameter calculated using linear rescaling.\n");
    }
    if (bFisherStultz)
    {
        printf("\n\nOrder parameter calculated using Fisher-Stultz equation.\n");
    }
    printf("Order Parameter: %10f \n", ISDorder);
    
    
    // Not sure what this one does. Sends to xmgrace?
    do_view(oenv,ftp2fn(efXVG,NFILE,fnm),"-nxy");

    // Closing.
    thanx(stderr);
    return 0;
}
