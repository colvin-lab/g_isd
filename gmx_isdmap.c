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


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// The below section is from gromacs code which was released separately under the LGPL version 2 license.

void write_mat_levels(FILE *out,int n_x, int n_y,int *nlevels,real lo,real hi,t_rgb rlo,t_rgb rhi)
{
    static const char matmap[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/?";
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
                matmap[i % lenMap],(*nlevels <= lenMap) ? ' ' : matmap[i/lenMap],
                                  (unsigned int)((int)((255*r) + 0.5)),
                                  (unsigned int)((int)((255*g) + 0.5)),
                                  (unsigned int)((int)((255*b) + 0.5)),
                                  (nlo*lo+i*hi)*invlevel);
    }
}



void write_mat_levels3(FILE *out,int n_x,int n_y,int *nlevels,real lo,real mid,real hi,t_rgb rlo,t_rgb rmid,t_rgb rhi)
{
    static const char matmap[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/?";
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
                matmap[i % lenMap],
                (*nlevels <= lenMap) ? ' ' : matmap[i/lenMap],
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
                matmap[(i+nmid) % lenMap],
                (*nlevels <= lenMap) ? ' ' : matmap[(i+nmid)/lenMap],
                 (unsigned int)((int)((255*r) + 0.5)),
                 (unsigned int)((int)((255*g) + 0.5)),
                 (unsigned int)((int)((255*b) + 0.5)),
                 ((*nlevels - 1 - nmid - i)*mid + i*hi)/clev_hi);
    }
}



void write_mat_pixels(FILE *out,int n_x,int n_y,real **matrix,real lo,real hi,int nlevels)
{
    static const char matmap[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/?";
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
                fprintf(out,"%c",matmap[c]);
                //printf("%c",matmap[c]);
            }
            else
            {
                fprintf(out,"%c%c",matmap[c % lenMap],matmap[c / lenMap]);
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
    static const char matmap[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()-_=+{}|;:',<.>/?";
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
                fprintf(out,"%c",matmap[c]);
            else
                fprintf(out,"%c%c",matmap[c % lenMap],matmap[c / lenMap]);
        }
        if (j > 0)
            fprintf(out,"\",\n");
        else
            fprintf(out,"\"\n");
    }
}



// The above section is from gromacs code which was released separately under the LGPL version 2 license.
/////////////////////////////////////////////////////////////////////////////////////////////////////////



int gmx_isdmap(int argc,char *argv[])
{
    const char *desc[] = {
        "[TT]g_isdmap[tt] implements a list of ISDMs designed to compare ",
        "two frames from a trajectory and quantify their inter-structure ",
        "distance (ISD). This tool calculates the matrix of all ISDs and ",
        "outputs to an xpm format file. The default ISDM if one is not ",
        "chosen by the user is RMSD. Only one ISDM can be chosen at a ",
        "time. The -xpm option is required. The xpm format can be converted ",
        "to eps with the xpm2ps tool. An upper threshold can be set with the ",
        "setmax option."
    };
    
    
    
    static gmx_bool bANG=FALSE, bDIH=FALSE, bANGDIH=FALSE, bDRMS=FALSE;
    static gmx_bool bPHIPSI=FALSE, bSRMS=FALSE, bPCOR=FALSE, bMAMMOTH=FALSE;
    static gmx_bool bACOR=FALSE, bESA=FALSE, bRMSD=FALSE, bMIR=FALSE;
    static gmx_bool bRG=FALSE, bSRG=FALSE, bE2E=FALSE, bSE2E=FALSE;
    static gmx_bool bRROT=FALSE, bSDRMS=FALSE;
    static real setmax = -1;
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
    t_rgb      rlo, rmd, rhi;
    int        ePBC;
    rvec       *x, *iframe, *jframe, *rframe, *cframe, rrot_xyz, xold;
    rvec       **frames;
    real       *nweights, *iweights, *xticks, *yticks;
    real       *diff, ISD, *avgdiff, avgISD, maxISD;
    real       **outmat;
    matrix     box, rrot, rrotx, rroty, rrotz;
    real       t, xpm_max, pi = 3.14159265358979;
    int        *maxframe, *rnum;
    int        i, j, k, m, n, iatoms, natoms, nframes;
    int        percentcalcs, noptions;
    static int nlevels = 81;
    gmx_bool   bDFLT, bFit, bMap, bCsv;
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
        { efXPM, "-map", "isdmap",   ffOPTWR }, 
        { efDAT, "-csv", "isdcsv",   ffOPTWR },
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
    
    // Check for error on -setmax before doing the calculations.
    if (setmax != -1)
    {
        if (setmax <= 0)
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
    bMap = opt2bSet("-map", NFILE, fnm);
    bCsv = opt2bSet("-csv", NFILE, fnm);
    
    
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
    snew(xticks,  nframes);
    snew(yticks,  nframes);
    snew(avgdiff, nframes);
    snew(outmat,  nframes);
    for (i=0; i<nframes; i++)
    {
        xticks[i] = (real)i;
        yticks[i] = (real)i;
        avgdiff[i] = 0;
        snew(outmat[i],nframes);
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
    maxISD = 0;
    avgISD = 0;
    
    
    
    /* Main calculation loop.
     */
    printf("\nCalculating results. \n");
    
    /* Originally this was designed to only loop through each pair of i and j
     * one time to save half of the calculations. Eventually it became 
     * impractical to make sure that each ISDM was symmetric, so now the
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
                outmat[i][j] = 0;
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
                bSDRMS || bPCOR || bACOR)
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
            
            // Add difference to the difference matrix.
            outmat[i][j] = ISD;
            
            // Update the max and avg difference for scaling.
            if (ISD > maxISD)
            {
                maxISD = ISD;
            }
            avgdiff[i] += ISD;
            
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
    
    // Find the final average of differences.
    for (i=0; i<nframes; i++)
    {
        avgISD += avgdiff[i];
    }
    avgISD /= nframes;
    
    
    if (bMap)
    {
        // Opens the output file.
        out = opt2FILE("-map", NFILE, fnm, "w");
        
        // Necessary for write_xpm.
        rlo.r = 0; rlo.g = 0; rlo.b = 1;
        rmd.r = 1; rmd.g = 1; rmd.b = 0;
        rhi.r = 1; rhi.g = 0; rhi.b = 0;
        
        // Choose the maximum value for the xpm.
        if (setmax == -1)
        {
            xpm_max = maxISD;
        }
        else
        {
            if (setmax > avgISD)
            {
                xpm_max = setmax;
            }
            else
            {
                gmx_fatal(FARGS,"\nThe argument for -setmax must be greater "
                          "than the average.\n");
            }
        }
        
        // Write to the output file.
        
        // Should check if these variables are even necessary.
        sprintf(buf,"Frame vs Frame");
        unsigned int tmp = 0;
        
        write_xpm_header(out, buf, "difference", "Frame", "Frame", FALSE);
        write_mat_levels3(out, nframes, nframes, &nlevels, 0, avgISD, xpm_max, rlo, rmd, rhi);
        write_xpm_axis(out, "x", tmp & MAT_SPATIAL_X, nframes, xticks);
        write_xpm_axis(out, "y", tmp & MAT_SPATIAL_Y, nframes, yticks);
        write_mat_pixels3(out, nframes, nframes, outmat, 0, avgISD, xpm_max, nlevels);
        
        //printf("\nrlo.r %f rlo.g %f rlo.b %f rhi.r %f rhi.g %f rhi.b %f \n",rlo.r,rlo.g,rlo.b,rhi.r,rhi.g,rhi.b);
        //printf("Test rlo.g %02X \n",(unsigned int)((int)((255*((((int)20)*rlo.g)*0.05))+0.5)));
        //printf("nlevels is now %d \n",nlevels);
        
        // Close the output file.
        ffclose(out);
    }
    
    if (bCsv)
    {
        // Opens the output file.
        out = opt2FILE("-csv", NFILE, fnm, "w");
        
        // Write output.
        for (i = 0; i < nframes; i++)
        {
            fprintf(out, "%10f", outmat[i][0]);
            for (j = 1; j < nframes; j++)
            {
                fprintf(out, ", %10f", outmat[i][j]);
            }
            fprintf(out, " \n");
        }
        
        // Close the output file.
        ffclose(out);
    }
    
    // Closing.
    thanx(stderr);
    return 0;
}
