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



//extern void mammothmod_(double xfit[], double xref[], int* nc, int* rnum, double* zscore);
double esa_analysis(int iatoms, rvec frame[], rvec rframe[], int sample);



real calc_rmsd(int iatoms, rvec frame[], rvec rframe[])
{
    /* This is a simplified version of something in the gromacs lib.
     * The accuracy is slightly improved by using doubles instead of reals.
     */
    int i, d;
    double xd, msdt;
    
    msdt = 0.0;
    for(i = 0; i < iatoms; i++)
    {
        // Calculations to find msd.
        for(d = 0; d < 3; d++)
        {
            xd     = (double)(frame[i][d] - rframe[i][d]);
            msdt += xd * xd;
        }
    }
    // Normalize for number of atoms.
    return (real)sqrt(msdt / iatoms);
}



real calc_rmsd_n(int iatoms, rvec frame[], rvec rframe[], real rmsd[])
{
    /* This is a simplified version of something in the gromacs lib.
     * The accuracy is slightly improved by using doubles instead of reals.
     */
    int i, d;
    double xd, rmsdi, rmsdt;
    
    rmsdt = 0.0;
    for(i = 0; i < iatoms; i++)
    {
        // Calculations to find msd.
        rmsdi = 0.0;
        for(d = 0; d < 3; d++)
        {
            xd     = (double)(frame[i][d] - rframe[i][d]);
            rmsdi += xd * xd;
        }
        rmsdi   = sqrt(rmsdi);
        rmsd[i] = (real)rmsdi;
        rmsdt  += rmsdi;
    }
    // Normalize for number of atoms.
    return (real)(rmsdt / iatoms);
}



real calc_srms(int iatoms, rvec frame[], rvec rframe[])
{
    /* Size-independent modification of RMSD. Scaled to a value [0,1].
     * 
     * Maiorov and Crippen, 1995.
     * SRMS = D(A,B) / [2*R(A)^2 + 2*R(B)^2 - D(A,B)^2]^(1/2)
     * 
     * The original equation (16) replaces 2*R(B)^2 with R(B)^2 which appears
     * to be a typo.
     */
    int i, d;
    double xd, rc, rr, msdt;
    
    msdt = 0.0; rc = 0.0; rr = 0.0;
    for(i = 0; i < iatoms; i++)
    {
        // Calculations to find MSD and Rg.
        for(d = 0; d < 3; d++)
        {
            // MSD portion.
            xd    = (double)(frame[i][d] - rframe[i][d]);
            msdt += xd * xd;
            // Rg portion.
            rc   += (double)( frame[i][d] *  frame[i][d]);
            rr   += (double)(rframe[i][d] * rframe[i][d]);
        }
    }
    // Normalize for the number of atoms.
    msdt /= iatoms;
    rc   /= iatoms;
    rr   /= iatoms;
    // Equation 16 appeared to have a typo.
    return (real)sqrt(msdt) / sqrt((2 * rc) + (2 * rr) - msdt);
}



real calc_srms_n(int iatoms, rvec frame[], rvec rframe[], real srms[])
{
    /* Size-independent modification of RMSD. Scaled to a value [0,1].
     * 
     * 
     */
    int i, d;
    double xd, rmsdi, rmsdt;
    
    rmsdt = 0.0;
    for(i = 0; i < iatoms; i++)
    {
        // Calculations to find msd.
        rmsdi = 0.0;
        for(d = 0; d < 3; d++)
        {
            xd     = (double)(frame[i][d] - rframe[i][d]);
            rmsdi += xd * xd;
        }
        rmsdi   = sqrt(rmsdi);
        srms[i] = (real)rmsdi;
        rmsdt  += rmsdi;
    }
    // Normalize for number of atoms.
    return (real)(rmsdt / iatoms);
}



real calc_msd(int iatoms, rvec frame[], rvec rframe[])
{
    /* This is a simplified version of something in the gromacs lib.
     * The accuracy is slightly improved by using doubles instead of reals.
     */
    int i, d;
    double xd, msdt;
    
    msdt = 0.0;
    for(i = 0; i < iatoms; i++)
    {
        // Calculations to find msd.
        for(d = 0; d < 3; d++)
        {
            xd    = (double)(frame[i][d] - rframe[i][d]);
            msdt += xd * xd;
        }
    }
    // Normalize for number of atoms.
    return (real)(msdt / iatoms);
}



real calc_msd_n(int iatoms, rvec frame[], rvec rframe[], real msd[])
{
    /* This is a simplified version of something in the gromacs lib.
     * The accuracy is slightly improved by using doubles instead of reals.
     */
    int i, d;
    double xd, msdi, msdt;
    
    msdt = 0.0;
    for(i = 0; i < iatoms; i++)
    {
        // Calculations to find msd.
        msdi = 0.0;
        for(d = 0; d < 3; d++)
        {
            xd    = (double)(frame[i][d] - rframe[i][d]);
            msdi += xd * xd;
        }
        msd[i] = (real)msdi;
        msdt  += msdi;
    }
    // Normalize for number of atoms.
    return (real)(msdt / iatoms);
}



real calc_mrms(int iatoms, rvec frame[], rvec rframe[])
{
    /*
     * The mmsdt is the average msd between the reference and its
     * mirrored structures mirrored around 3 planes. Used to normalize
     * the rmsd value.
     */
    int i, j, d;
    double xd, xdm, msdt, mmsdt;
    
    // Initialize to zero.
    msdt  = 0.0;
    mmsdt = 0.0;
    // Main loop.
    for(i = 0; i < iatoms; i++)
    {
        // Copy over coordinates from xref and rotate them.
        for(d = 0; d < 3; d++)
        {
            // This is the same as multiplying the reference structure by the 
            // negative identity matrix.
            
            // Distance to mirrored structure.
            xdm    = (double)(frame[i][d] + rframe[i][d]);
            mmsdt += xdm * xdm;
            
            // Calculate ith msd.
            xd     = (double)(frame[i][d] - rframe[i][d]);
            msdt  += xd  * xd;
        }
    }
    
    // Normalize for number of atoms.
    msdt  /= iatoms;
    mmsdt /= iatoms;
    // Scale by the mirrored reference msd and take the root.
    return (real)sqrt(msdt / mmsdt);
}



real calc_mrms_n(int iatoms, rvec frame[], rvec rframe[], real mrms[])
{
    /*
     * The mmsdt is the average msd between the reference and its
     * mirrored structures mirrored around 3 planes. Used to normalize
     * the rmsd value.
     */
    int i, j, d;
    double xd, xdm, msdi, mmsdt, mmrms;
    
    // Initialize to zero.
    mmrms = 0.0;
    mmsdt = 0.0;
    
    // Main loop.
    for(i = 0; i < iatoms; i++)
    {
        // Copy over coordinates from xref and rotate them.
        msdi = 0.0;
        for(d = 0; d < 3; d++)
        {
            // This is the same as multiplying the reference structure by the 
            // negative identity matrix.
            
            // Distance to mirrored structure.
            xdm    = (double)(frame[i][d] + rframe[i][d]);
            mmsdt += xdm * xdm;
            
            // Calculate ith msd.
            xd     = (double)(frame[i][d] - rframe[i][d]);
            msdi  += xd  * xd;
        }
        mrms[i] = (real)msdi;
    }
    
    // Normalize for number of atoms.
    mmsdt /= iatoms;
    // Scale by the mirrored reference msd and take the root.
    for (i = 0; i < iatoms; i++)
    {
        mrms[i] = (real)sqrt(((double)mrms[i]) / mmsdt);
        // The total msd is always less than mmsd, but individual atoms may 
        // be greater.
        if (mrms[i] > 1.0)
        {
            mrms[i] = 1.0;
        }
        
        // To calculate the mean srms.
        mmrms += mrms[i];
    }
    return (real)(mmrms / iatoms);
}



real calc_drms(int iatoms, rvec frame[], rvec rframe[])
{
    int i, j;
    double sum_dist, drmsi, diff_leg;
    real ileg, rleg;
    
    // Final sum.
    sum_dist = 0.0;
    // Loop through each atom.
    for (i = 0; i < iatoms; i++)
    {
        drmsi = 0.0;
        // Compare to every other atom.
        for (j = 0; j < iatoms; j++)
        {
            // It makes sense to skip i == j, but it is zero and is simpler 
            // to leave.
            
            // Find the length of the vector made by every pair of atoms.
            ileg = distance2(frame[i], frame[j]);
            rleg = distance2(rframe[i], rframe[j]);
            // Find difference in the distances between the pair of atoms.
            diff_leg = sqrt((double)ileg) - sqrt((double)rleg);
            // Update the sum by squaring the difference.
            drmsi += diff_leg * diff_leg;
        }
        // Divide by n - 1 and update the main sum.
        sum_dist += drmsi / (iatoms - 1);
    }
    
    // Normalize by the number of summed differences and take the root.
    return (real)sqrt(sum_dist / iatoms);
}



real calc_drms_n(int iatoms, rvec frame[], rvec rframe[], real drms[])
{
    int i, j;
    double sum_dist, drmsi, diff_leg;
    real ileg, rleg;
    
    // Don't assume that drms is zeros already.
    for (i = 0; i < iatoms; i++)
    {
        drms[i] = 0.0;
    }
    
    // Final sum.
    sum_dist = 0.0;
    // Loop through each atom.
    for (i = 0; i < iatoms; i++)
    {
        // Compare to every other atom.
        for (j = 0; j < iatoms; j++)
        {
            // It makes sense to skip i == j, but it is zero and simpler 
            // to leave.
            
            // Find the length of the vector made by every pair of atoms.
            ileg = distance2(frame[i], frame[j]);
            rleg = distance2(rframe[i], rframe[j]);
            // Find difference in the distances between the pair of atoms.
            diff_leg = sqrt((double)ileg) - sqrt((double)rleg);
            // Update the sum by squaring the difference.
            drmsi   += diff_leg * diff_leg;
        }
        // Normalize the sum by the number of differences.
        drmsi    /= sqrt(drmsi / (iatoms - 1));
        // Update the main sum.
        drms[i]   = (real)drmsi;
        sum_dist += drmsi;
    }
    
    // Normalize by the number of summed distances.
    return (real)(sum_dist / iatoms);
}



real calc_sdrms(int iatoms, rvec frame[], rvec rframe[])
{
    int i, j;
    double sum_dist, drmsi, diff_leg, Rgi, Rgr, Rg2;
    real ileg, rleg;
    
    //Normalization
    Rgi = 0.0;
    Rgr = 0.0;
    // Final sum.
    sum_dist = 0.0;
    // Loop through each atom.
    for (i = 0; i < iatoms; i++)
    {
        drmsi = 0.0;
        // Compare to every other atom.
        for (j = 0; j < iatoms; j++)
        {
            // It makes sense to skip i == j, but it is zero and is simpler 
            // to leave.
            
            // Find the length of the vector made by every pair of atoms.
            ileg = distance2(frame[i], frame[j]);
            rleg = distance2(rframe[i], rframe[j]);
            // Normalization
            Rgi += (double)ileg;
            Rgr += (double)rleg;
            // Find difference in the distances between the pair of atoms.
            diff_leg = sqrt((double)ileg) - sqrt((double)rleg);
            // Update the sum by squaring the difference.
            drmsi += diff_leg * diff_leg;
        }
        // Divide by n -1 and update the main sum.
        sum_dist += drmsi / (iatoms - 1);
    }
    
    // Solves for Rg which is used to scale for molecule size.
    Rgi = sqrt(Rgi / (2 * iatoms * iatoms));
    Rgr = sqrt(Rgr / (2 * iatoms * iatoms));
    Rg2  = Rgi + Rgr; // Equal to twice the mean Rg.
    // First, normalize by the number of summed differences.
    sum_dist /= iatoms;
    // Take the root and scale by the molecule size (twice the mean Rg).
    sum_dist = sqrt(sum_dist) / Rg2;
    // Expected output is between 0 and 1, but the scaling allows for a small
    // chance of an output being > 1. This would be undesirable.
    if (sum_dist > 1.0)
    {
        return 1.0;
    }
    else
    {
        return sum_dist;
    }
}



real calc_sdrms_n(int iatoms, rvec frame[], rvec rframe[], real sdrms[])
{
    int i, j;
    real ileg, rleg;
    double diff_leg, drmsi, sum_dist, Rgi, Rgr, Rg2;
    
    //Normalization
    Rgi = 0.0;
    Rgr = 0.0;
    // Final sum.
    sum_dist = 0.0;
    // Loop through each atom.
    for (i = 0; i < iatoms; i++)
    {
        // Compare to every other atom.
        for (j = 0; j < iatoms; j++)
        {
            // It makes sense to skip i == j, but it is zero and is simpler 
            // to leave.
            
            // Find the length of the vector made by every pair of atoms.
            ileg = distance2(frame[i], frame[j]);
            rleg = distance2(rframe[i], rframe[j]);
            // Normalization
            Rgi += (double)ileg;
            Rgr += (double)rleg;
            // Find difference in the distances between the pair of atoms.
            diff_leg = sqrt((double)ileg) - sqrt((double)rleg);
            // Update the sum by squaring the difference.
            drmsi += diff_leg * diff_leg;
        }
        // Normalize the sum by the number of differences and take the root.
        drmsi     = sqrt(drmsi / (iatoms - 1));
        sdrms[i]  = (real)drmsi;
        // Update the main sum.
        sum_dist += drmsi;
    }
    
    // Solves for Rg which is used for normalization.
    Rgi = sqrt(Rgi / (2 * iatoms * iatoms));
    Rgr = sqrt(Rgr / (2 * iatoms * iatoms));
    Rg2 = Rgi + Rgr; // Twice the mean Rg.
    // Normalize by the number of summed differences.
    sum_dist /= iatoms;
    // Scale by the molecule size. Twice the geometric mean Rg.
    sum_dist /= Rg2;
    // Repeat per atom.
    for (i = 0; i < iatoms; i++)
    {
        sdrms[i] /= Rg2;
        if (sdrms[i] > 1.0)
        {
            // After scaling there is a small chance of the result being > 1.
            sdrms[i] = 1.0;
        }
    }
    // Expected output is between 0 and 1, but the scaling allows for a small
    // chance of an output being >1. This would be undesirable.
    if (sum_dist > 1.0)
    {
        return 1.0;
    }
    else
    {
        return (real)sum_dist;
    }
}



real calc_ang(int iatoms, rvec frame[], rvec rframe[])
{
    // Error checking.
    if (iatoms < 4)
    {
        gmx_fatal(FARGS, "Need at least 4 atoms in index group to use calc_ang.");
    }
    
    // Initializing variables.
    int i;
    rvec vec1, vec2;
    double cosx, cosx2, cosy, cosy2, cosxy, sumcosxy;
    
    sumcosxy = 0.0;
    // There are n - 2 backbone angles.
    for (i = 1; i < (iatoms - 1); i++)
    {
        /* Find the value of cosx and cosy.
         * 
         * x is the angle of two vectors made by three atom coords
         * in frame. y is the same for rframe.
         */
        rvec_sub(frame[i - 1],  frame[i], vec1);
        rvec_sub(frame[i + 1],  frame[i], vec2);
        cosx  = cos_angle(vec1, vec2);
        cosx2 = cosx * cosx;
        
        rvec_sub(rframe[i - 1], rframe[i], vec1);
        rvec_sub(rframe[i + 1], rframe[i], vec2);
        cosy  = cos_angle(vec1, vec2);
        cosy2 = cosy * cosy;
        
        /* Find the value of cos(x - y) using trig identities.
         * 
         * cos(x - y) = cos(x) * cos(y) + sin(x) * sin(y)
         * sin(x) = sin(acos(cos(x)))
         * sin(acos(z)) = sqrt(1 - (z ^ 2))
         *
         * Therefore:
         * cos(x - y) = 
         * cos(x) * cos(y) + sqrt(1 - (cos(x) ^ 2)) * sqrt(1 - (cos(y) ^ 2)) 
         */
        
        // Need to check for numerical precision to avoid negative root.
        if (cosx2 > 1.0)
        {
            cosx2 = 1.0;
        }
        if (cosy2 > 1.0)
        {
            cosy2 = 1.0;
        }
        // This results in a number between -1.0 and +1.0.
        cosxy = cosx * cosy + sqrt(1.0 - cosx2) * sqrt(1.0 - cosy2);
        // Update sum.
        sumcosxy += cosxy;
    }
    
    // Divide by n - 2.
    sumcosxy /= (iatoms - 2);
    // Rescale from [+1.0, -1.0] to [0.0, 1.0].
    sumcosxy  = (sumcosxy - 1.0) / (-2.0);
    // Finished.
    return (real)sumcosxy;
}



real calc_ang_n(int iatoms, rvec frame[], rvec rframe[], real ang[])
{
    // Error checking.
    if (iatoms < 4)
    {
        gmx_fatal(FARGS, "Need at least 4 atoms in index group to use calc_ang.");
    }
    
    // Initializing variables.
    int i;
    rvec vec1, vec2;
    double cosx, cosx2, cosy, cosy2, cosxy, sumcosxy;
    
    sumcosxy        = 0.0;
    ang[0]          = 0.0;
    ang[iatoms - 1] = 0.0;
    // There are n - 2 backbone angles.
    for (i = 1; i < (iatoms - 1); i++)
    {
        /* Find the value of cosx and cosy.
         * 
         * x is the angle of two vectors made by three atom coords
         * in frame. y is the same for rframe.
         */
        rvec_sub(frame[i - 1],  frame[i], vec1);
        rvec_sub(frame[i + 1],  frame[i], vec2);
        cosx  = cos_angle(vec1, vec2);
        cosx2 = cosx * cosx;
        
        rvec_sub(rframe[i - 1], rframe[i], vec1);
        rvec_sub(rframe[i + 1], rframe[i], vec2);
        cosy  = cos_angle(vec1, vec2);
        cosy2 = cosy * cosy;
        
        /* Find the value of cos(x - y) using trig identities.
         * 
         * cos(x - y) = cos(x) * cos(y) + sin(x) * sin(y)
         * sin(x) = sin(acos(cos(x)))
         * sin(acos(z)) = sqrt(1 - (z ^ 2))
         *
         * Therefore:
         * cos(x - y) = 
         * cos(x) * cos(y) + sqrt(1 - (cos(x) ^ 2)) * sqrt(1 - (cos(y) ^ 2)) 
         */
        
        // Need to check for numerical precision to avoid negative root.
        if (cosx2 > 1.0)
        {
            cosx2 = 1.0;
        }
        if (cosy2 > 1.0)
        {
            cosy2 = 1.0;
        }
        // This results in a number between -1.0 and +1.0.
        cosxy = cosx * cosy + sqrt(1.0 - cosx2) * sqrt(1.0 - cosy2);
        // Rescale from [+1.0, -1.0] to [0.0, 1.0].
        cosxy = (cosxy - 1.0) / (-2.0);
        // Update sum.
        ang[i]    = (real)cosxy;
        sumcosxy += cosxy;
    }
    
    // Divide by n - 2.
    sumcosxy /= (iatoms - 2);
    // Finished.
    return (real)sumcosxy;
}



real calc_ang2(int iatoms, rvec frame[], rvec rframe[])
{
    // Error checking.
    if (iatoms < 4)
    {
        gmx_fatal(FARGS, "Need at least 4 atoms in index group to use calc_ang.");
    }
    
    // Initializing variables.
    int i;
    rvec vec1, vec2;
    double irang, irang2, sum_angs;
    double pi20 = 3.14159265358979323846;
    real iang, rang;
    
    sum_angs = 0.0;
    // There are n - 2 backbone angles.
    for (i = 1; i < (iatoms - 1); i++)
    {
        /* Find the value of cosx and cosy.
         * 
         * x is the angle of two vectors made by three atom coords
         * in frame. y is the same for rframe.
         */
        rvec_sub(frame[i - 1],  frame[i], vec1);
        rvec_sub(frame[i + 1],  frame[i], vec2);
        iang   = gmx_angle(vec1, vec2);
        
        rvec_sub(rframe[i - 1], rframe[i], vec1);
        rvec_sub(rframe[i + 1], rframe[i], vec2);
        rang   = gmx_angle(vec1, vec2);
        
        // Absolute value of the difference in angles.
        irang  = (double)(iang - rang);
        irang2 = irang * irang;
        
        // Update sum.
        sum_angs += irang2;
    }
    
    // Divide by n - 2 angles and take the root.
    sum_angs  = sqrt(sum_angs / (iatoms - 2));
    // Rescale from [0, pi] to [0, 1].
    sum_angs /= pi20;
    // Finished.
    return (real)sum_angs;
}



real calc_ang2_n(int iatoms, rvec frame[], rvec rframe[], real ang[])
{
    // Error checking.
    if (iatoms < 4)
    {
        gmx_fatal(FARGS, "Need at least 4 atoms in index group to use calc_ang.");
    }
    
    // Initializing variables.
    int i;
    rvec vec1, vec2;
    double irang, irang2, sum_angs;
    double pi20 = 3.14159265358979323846;
    real iang, rang;
    
    sum_angs        = 0.0;
    ang[0]          = 0.0;
    ang[iatoms - 1] = 0.0;
    // There are n - 2 backbone angles.
    for (i = 1; i < (iatoms - 1); i++)
    {
        /* Find the value of cosx and cosy.
         * 
         * x is the angle of two vectors made by three atom coords
         * in frame. y is the same for rframe.
         */
        rvec_sub(frame[i - 1],  frame[i], vec1);
        rvec_sub(frame[i + 1],  frame[i], vec2);
        iang   = gmx_angle(vec1, vec2);
        
        rvec_sub(rframe[i - 1], rframe[i], vec1);
        rvec_sub(rframe[i + 1], rframe[i], vec2);
        rang   = gmx_angle(vec1, vec2);
        
        // Absolute value of the difference in angles.
        irang  = (double)(iang - rang);
        irang2 = irang * irang;
        // Square root and rescale. Reusing variable irang.
        irang     = sqrt(irang2) / pi20;
        ang[i]    = (real)irang;
        // Update sum.
        sum_angs += irang;
    }
    // Divide by n - 2 backbone angles.
    sum_angs /= iatoms - 2;
    // Finished.
    return (real)sum_angs;
}



real calc_dih(int iatoms, rvec frame[], rvec rframe[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS,"Need at least 6 atoms in index group to use calc_dih.");
    }
    
    // Initializing variables.
    int i;
    double cosxy, sumcosxy;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sumcosxy = 0.0;
    // There are n - 3 dihedral angles.
    for (i = 3; i < iatoms; i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i - 3], frame[i - 2], vec1);
        rvec_sub(frame[i - 1], frame[i - 2], vec2);
        rvec_sub(frame[i - 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i - 3], rframe[i - 2], vec1);
        rvec_sub(rframe[i - 1], rframe[i - 2], vec2);
        rvec_sub(rframe[i - 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for cos(x - y).
        cosxy     = cos((double)(iang - rang));
        // Update sum.
        sumcosxy += cosxy;
    }
    
    // Normalize for number of angles summed.
    sumcosxy /= (iatoms - 3);
    // Rescale from [-1, +1] to [0, 1]
    sumcosxy  = (sumcosxy - 1.0) / (-2.0);
    // Finished.
    return (real)sumcosxy;
}



real calc_dih_n(int iatoms, rvec frame[], rvec rframe[], real dih[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS,"Need at least 6 atoms in index group to use calc_dih.");
    }
    
    // Initializing variables.
    int i;
    double cosxy, sumcosxy;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sumcosxy = 0.0;
    for (i = 0; i < iatoms; i++)
    {
        dih[i] = 0.0;
    }
    // There are n - 3 dihedral angles.
    for (i = 3; i < iatoms; i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i - 3], frame[i - 2], vec1);
        rvec_sub(frame[i - 1], frame[i - 2], vec2);
        rvec_sub(frame[i - 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i - 3], rframe[i - 2], vec1);
        rvec_sub(rframe[i - 1], rframe[i - 2], vec2);
        rvec_sub(rframe[i - 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for cos(x - y).
        cosxy       = cos((double)(iang - rang));
        // Rescale from [-1, +1] to [0, 1].
        cosxy      += (cosxy - 1.0) / (-2.0);
        // Array dih stores the means of two dihedral angles.
        dih[i - 2] += (real)cosxy;
        dih[i - 1] += (real)cosxy;
        // Update sum.
        sumcosxy   += cosxy;
    }
    // Array dih stores the means of two dihedral angles.
    // *** First and last indeces = 0.0.
    // *** Second and second last indeces only include one dihedral angle.
    sumcosxy += (double)dih[1];
    sumcosxy += (double)dih[iatoms - 2];
    for (i = 2; i < (iatoms - 2); i++)
    {
        dih[i]   /= 2.0;
        sumcosxy += (double)dih[i];
    }
    // Divide by the number of summed indeces.
    sumcosxy /= (iatoms - 2);
    // Finished.
    return (real)sumcosxy;
}



real calc_dih2(int iatoms, rvec frame[], rvec rframe[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS,"Need at least 6 atoms in index group to use calc_dih.");
    }
    
    // Initializing variables.
    int i;
    double dihi, dihi2, sum_dihs;
    double pi20 = 3.14159265358979323846;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sum_dihs = 0.0;
    // There are n - 3 dihedral angles.
    for (i = 3; i < iatoms; i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i - 3], frame[i - 2], vec1);
        rvec_sub(frame[i - 1], frame[i - 2], vec2);
        rvec_sub(frame[i - 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i - 3], rframe[i - 2], vec1);
        rvec_sub(rframe[i - 1], rframe[i - 2], vec2);
        rvec_sub(rframe[i - 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for (iang - rang) ^ 2.
        dihi  = (double)(iang - rang);
        // Rescale dihi from [-2.0 * pi, +2.0 * pi] to [-pi, +pi].
        if (dihi > pi20)
        {
            dihi -= 2 * pi20;
        }
        else if (dihi < -pi20)
        {
            dihi += 2 * pi20;
        }
        dihi2 = dihi * dihi;
        // Update sum.
        sum_dihs += dihi2;
    }
    
    // Divide by n - 3 and take sqrt.
    sum_dihs  = sqrt(sum_dihs / (iatoms - 3));
    // Rescale from [0, pi] to [0, 1].
    sum_dihs /= pi20;
    // Finished.
    return (real)sum_dihs;
}



real calc_dih2_n(int iatoms, rvec frame[], rvec rframe[], real dih[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS,"Need at least 6 atoms in index group to use calc_dih.");
    }
    
    // Initializing variables.
    int i;
    double dihi, dihi2, sum_dihs;
    double pi20 = 3.14159265358979323846;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sum_dihs = 0.0;
    for (i = 0; i < iatoms; i++)
    {
        dih[i] = 0.0;
    }
    // There are n - 3 dihedral angles.
    for (i = 3; i < iatoms; i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i - 3], frame[i - 2], vec1);
        rvec_sub(frame[i - 1], frame[i - 2], vec2);
        rvec_sub(frame[i - 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i - 3], rframe[i - 2], vec1);
        rvec_sub(rframe[i - 1], rframe[i - 2], vec2);
        rvec_sub(rframe[i - 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for (iang - rang) ^ 2.
        dihi  = (double)(iang - rang);
        // Rescale dihi from [-2.0 * pi, +2.0 * pi] to [-pi, +pi].
        if (dihi > pi20)
        {
            dihi -= 2 * pi20;
        }
        else if (dihi < -pi20)
        {
            dihi += 2 * pi20;
        }
        dihi2 = dihi * dihi;
        // Take root, rescale, and update output. Reusing variable dihi.
        dihi        = sqrt(dihi2) / pi20;
        // Array dih stores the means of two dihedral angles.
        dih[i - 1] += (real)dihi;
        dih[i - 2] += (real)dihi;
    }
    // Array dih stores the means of two dihedral angles.
    // *** First and last indeces = 0.0.
    // *** Second and second last indeces only include one dihedral angle.
    sum_dihs += dih[1];
    sum_dihs += dih[iatoms - 2];
    for (i = 2; i < (iatoms - 2); i++)
    {
        dih[i] /= 2.0;
        sum_dihs += (double)dih[i];
    }
    // Divide by the number of summed indeces.
    return (real)(sum_dihs / (iatoms - 2));
}



real calc_angdih(int iatoms, rvec frame[], rvec rframe[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS, "Need at least 6 atoms in index group to use calc_angdih.");
    }
    // Initializing variables.
    int i;
    double pi20 = 3.14159265358979323846;
    double cosx, cosx2, cosy, cosy2, cosxy, sum_angs, sum_dihs;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sum_angs = 0.0; sum_dihs = 0.0;
    // There are n - 2 backbone angles.
    for (i = 1; i < (iatoms - 1); i++)
    {
        /* Find the value of cosx and cosy.
         * 
         * x is the angle of two vectors made by three atom coords
         * in frame. y is the same for rframe.
         */
        rvec_sub(frame[i - 1],  frame[i], vec1);
        rvec_sub(frame[i + 1],  frame[i], vec2);
        cosx  = cos_angle(vec1, vec2);
        cosx2 = cosx * cosx;
        
        rvec_sub(rframe[i - 1], rframe[i], vec1);
        rvec_sub(rframe[i + 1], rframe[i], vec2);
        cosy  = cos_angle(vec1, vec2);
        cosy2 = cosy * cosy;
        
        /* Find the value of cos(x - y) using trig identities.
         * 
         * cos(x - y) = cos(x) * cos(y) + sin(x) * sin(y)
         * sin(x) = sin(acos(cos(x)))
         * sin(acos(z)) = sqrt(1 - (z ^ 2))
         *
         * Therefore:
         * cos(x - y) = 
         * cos(x) * cos(y) + sqrt(1 - (cos(x) ^ 2)) * sqrt(1 - (cos(y) ^ 2)) 
         */
        
        // Need to check for numerical precision to avoid negative root.
        if (cosx2 > 1.0)
        {
            cosx2 = 1.0;
        }
        if (cosy2 > 1.0)
        {
            cosy2 = 1.0;
        }
        // This results in a number between -1.0 and +1.0.
        cosxy = cosx * cosy + sqrt(1.0 - cosx2) * sqrt(1.0 - cosy2);
        // Update sum.
        sum_angs += cosxy;
    }
    
    // There are n - 3 dihedral angles.
    for (i = 3; i < iatoms; i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i - 3], frame[i - 2], vec1);
        rvec_sub(frame[i - 1], frame[i - 2], vec2);
        rvec_sub(frame[i - 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i - 3], rframe[i - 2], vec1);
        rvec_sub(rframe[i - 1], rframe[i - 2], vec2);
        rvec_sub(rframe[i - 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for cos(x - y).
        cosxy = cos((double)(iang - rang));
        // Update sum.
        sum_dihs += cosxy;
    }
    
    // Divide by number of angles summed.
    sum_angs  /= (iatoms - 2);
    sum_dihs  /= (iatoms - 3);
    // Rescale from [-1, +1] to [0, 1].
    sum_angs   = (sum_angs - 1.0) / (-2.0);
    sum_dihs   = (sum_dihs - 1.0) / (-2.0);
    // Return geometric mean.
    return (real)sqrt(sum_angs * sum_dihs);
}



real calc_angdih_n(int iatoms, rvec frame[], rvec rframe[], real angdih[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS, "Need at least 6 atoms in index group to use calc_angdih.");
    }
    // Initializing variables.
    int i;
    double pi20 = 3.14159265358979323846;
    double cosx, cosx2, cosy, cosy2, cosxy, sum_angs, sum_dihs;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sum_angs = 0.0; sum_dihs = 0.0;
    for (i = 0; i < iatoms; i++)
    {
        angdih[i] = 0.0;
    }
    // There are n - 2 backbone angles.
    for (i = 1; i < (iatoms - 1); i++)
    {
        /* Find the value of cosx and cosy.
         * 
         * x is the angle of two vectors made by three atom coords
         * in frame. y is the same for rframe.
         */
        rvec_sub(frame[i - 1],  frame[i], vec1);
        rvec_sub(frame[i + 1],  frame[i], vec2);
        cosx  = cos_angle(vec1, vec2);
        cosx2 = cosx * cosx;
        
        rvec_sub(rframe[i - 1], rframe[i], vec1);
        rvec_sub(rframe[i + 1], rframe[i], vec2);
        cosy  = cos_angle(vec1, vec2);
        cosy2 = cosy * cosy;
        
        /* Find the value of cos(x - y) using trig identities.
         * 
         * cos(x - y) = cos(x) * cos(y) + sin(x) * sin(y)
         * sin(x) = sin(acos(cos(x)))
         * sin(acos(z)) = sqrt(1 - (z ^ 2))
         *
         * Therefore:
         * cos(x - y) = 
         * cos(x) * cos(y) + sqrt(1 - (cos(x) ^ 2)) * sqrt(1 - (cos(y) ^ 2)) 
         */
        
        // Need to check for numerical precision to avoid negative root.
        if (cosx2 > 1.0)
        {
            cosx2 = 1.0;
        }
        if (cosy2 > 1.0)
        {
            cosy2 = 1.0;
        }
        // This results in a number between -1.0 and +1.0.
        cosxy = cosx * cosy + sqrt(1.0 - cosx2) * sqrt(1.0 - cosy2);
        // Update sum.
        sum_angs += cosxy;
    }
    
    // There are n - 3 dihedral angles.
    for (i = 3; i < iatoms; i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i - 3], frame[i - 2], vec1);
        rvec_sub(frame[i - 1], frame[i - 2], vec2);
        rvec_sub(frame[i - 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i - 3], rframe[i - 2], vec1);
        rvec_sub(rframe[i - 1], rframe[i - 2], vec2);
        rvec_sub(rframe[i - 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for cos(x - y).
        cosxy = cos((double)(iang - rang));
        // Update sum.
        sum_dihs += cosxy;
    }
    
    // Divide by number of angles summed.
    sum_angs  /= (iatoms - 2);
    sum_dihs  /= (iatoms - 3);
    // Rescale from [-1, +1] to [0, 1].
    sum_angs   = (sum_angs - 1.0) / (-2.0);
    sum_dihs   = (sum_dihs - 1.0) / (-2.0);
    // Return geometric mean.
    return (real)sqrt(sum_angs * sum_dihs);
}



real calc_angdih2(int iatoms, rvec frame[], rvec rframe[], double ang_multi)
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS, "Need at least 6 atoms in index group to use calc_angdih.");
    }
    // Initializing variables.
    int i;
    double pi20 = 3.14159265358979323846;
    double irang, irang2, dihi, dihi2, sum_angs, sum_dihs, sum_angdih;
    double sum_denom;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sum_angs = 0.0; sum_dihs = 0.0;
    // There are n - 2 backbone angles.
    for (i = 1; i < (iatoms - 1); i++)
    {
        /* Find the value of cosx and cosy.
         * 
         * x is the angle of two vectors made by three atom coords
         * in frame. y is the same for rframe.
         */
        rvec_sub(frame[i - 1],  frame[i], vec1);
        rvec_sub(frame[i + 1],  frame[i], vec2);
        iang   = gmx_angle(vec1, vec2);
        
        rvec_sub(rframe[i - 1], rframe[i], vec1);
        rvec_sub(rframe[i + 1], rframe[i], vec2);
        rang   = gmx_angle(vec1, vec2);
        
        // Squared difference of angles.
        irang  = iang - rang;
        irang2 = irang * irang;
        
        // Update sum.
        sum_angs += irang2;
    }
    sum_angs *= ang_multi;
    
    // There are n - 3 dihedral angles.
    for (i = 3; i < iatoms; i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i - 3], frame[i - 2], vec1);
        rvec_sub(frame[i - 1], frame[i - 2], vec2);
        rvec_sub(frame[i - 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i - 3], rframe[i - 2], vec1);
        rvec_sub(rframe[i - 1], rframe[i - 2], vec2);
        rvec_sub(rframe[i - 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for the difference between the two angles.
        dihi  = (double)(iang - rang);
        // Rescale dihi from [-2.0 * pi, +2.0 * pi] to [-pi, +pi].
        if (dihi > pi20)
        {
            dihi -= 2 * pi20;
        }
        else if (dihi < -pi20)
        {
            dihi += 2 * pi20;
        }
        dihi2 = dihi * dihi;
        // Update sum.
        sum_dihs += dihi2;
    }
    
    // Calculate denominator based on the number of angles summed.
    sum_denom   = (ang_multi + 1.0) * iatoms - (2.0 * ang_multi - 3.0); 
    // Divide by angles summed and take sqrt.
    sum_angdih  = sqrt((sum_angs + sum_dihs) / sum_denom);
    // Rescale from [0, pi] to [0, 1].
    sum_angdih /= pi20;
    // Output.
    return (real)sum_angdih;
}



real calc_angdih2_n(int iatoms, rvec frame[], rvec rframe[], real angdih[], double ang_multi)
{
    // Returns the root mean of the squared difference of angles.
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS, "Need at least 6 atoms in index group to use calc_angdih.");
    }
    // Initializing variables.
    int i;
    double pi20 = 3.14159265358979323846;
    double irang, irang2, dihi, dihi2, sum_angdih;
    double sum_denom;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sum_angdih = 0.0;
    for (i = 0; i < iatoms; i++)
    {
        angdih[i] = 0.0;
    }
    // There are n - 3 dihedral angles.
    for (i = 3; i < iatoms; i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i - 3], frame[i - 2], vec1);
        rvec_sub(frame[i - 1], frame[i - 2], vec2);
        rvec_sub(frame[i - 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i - 3], rframe[i - 2], vec1);
        rvec_sub(rframe[i - 1], rframe[i - 2], vec2);
        rvec_sub(rframe[i - 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for the difference between the two angles.
        dihi  = (double)(iang - rang);
        // Rescale dihi from [-2.0 * pi, +2.0 * pi] to [-pi, +pi].
        if (dihi > pi20)
        {
            dihi -= 2 * pi20;
        }
        else if (dihi < -pi20)
        {
            dihi += 2 * pi20;
        }
        dihi2 = dihi * dihi;
        // Update output array.
        angdih[i - 2] += (real)dihi2;
        angdih[i - 1] += (real)dihi2;
    }
    
    // The dih portion of angdih stores the means of two dihedral angles.
    // *** First and last indeces = 0.0.
    // *** Second and second last indeces only include one dihedral angle.
    for (i = 2; i < (iatoms - 2); i++)
    {
        angdih[i] /= 2.0;
    }
    
    // There are n - 2 backbone angles.
    for (i = 1; i < (iatoms - 1); i++)
    {
        /* Find the angle between consecutive pairs of backbone atoms.
         */
        rvec_sub(frame[i - 1],  frame[i], vec1);
        rvec_sub(frame[i + 1],  frame[i], vec2);
        iang   = gmx_angle(vec1, vec2);
        
        rvec_sub(rframe[i - 1], rframe[i], vec1);
        rvec_sub(rframe[i + 1], rframe[i], vec2);
        rang   = gmx_angle(vec1, vec2);
        
        // Squared difference of angles.
        irang      = iang - rang;
        irang2     = irang * irang;
        // Apply multiplier. Update output array.
        irang2    *= ang_multi;
        angdih[i] += (real)irang2;
    }
    
    // Finalize output array.
    for (i = 1; i < (iatoms - 1); i++)
    {
        // Take the mean of ang and dih. Adjust for ang_multi.
        angdih[i]  /= ang_multi + 1.0;
        // Take sqrt.
        angdih[i]   = sqrt(angdih[i]);
        // Rescale from [0, pi] to [0, 1].
        angdih[i]  /= pi20;
        // To calculate the mean of the indeces.
        sum_angdih += (double)angdih[i];
    }
    // Divide by the number of summed indeces.
    return (real)(sum_angdih / (iatoms / 2));
}



real calc_angdih2g(int iatoms, rvec frame[], rvec rframe[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS, "Need at least 6 atoms in index group to use calc_angdih.");
    }
    // Initializing variables.
    int i;
    double pi20 = 3.14159265358979323846;
    double irang, irang2, dihi, dihi2, sum_angs, sum_dihs, sum_angdih;
    double sum_denom;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sum_angs = 0.0; sum_dihs = 0.0;
    // There are n - 2 backbone angles.
    for (i = 1; i < (iatoms - 1); i++)
    {
        /* Find the value of cosx and cosy.
         * 
         * x is the angle of two vectors made by three atom coords
         * in frame. y is the same for rframe.
         */
        rvec_sub(frame[i - 1],  frame[i], vec1);
        rvec_sub(frame[i + 1],  frame[i], vec2);
        iang   = gmx_angle(vec1, vec2);
        
        rvec_sub(rframe[i - 1], rframe[i], vec1);
        rvec_sub(rframe[i + 1], rframe[i], vec2);
        rang   = gmx_angle(vec1, vec2);
        
        // Squared difference of angles.
        irang  = iang - rang;
        irang2 = irang * irang;
        
        // Update sum.
        sum_angs += irang2;
    }
    
    // There are n - 3 dihedral angles.
    for (i = 3; i < iatoms; i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i - 3], frame[i - 2], vec1);
        rvec_sub(frame[i - 1], frame[i - 2], vec2);
        rvec_sub(frame[i - 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i - 3], rframe[i - 2], vec1);
        rvec_sub(rframe[i - 1], rframe[i - 2], vec2);
        rvec_sub(rframe[i - 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for the difference between the two angles.
        dihi  = (double)(iang - rang);
        // Rescale dihi from [-2.0 * pi, +2.0 * pi] to [-pi, +pi].
        if (dihi > pi20)
        {
            dihi -= 2 * pi20;
        }
        else if (dihi < -pi20)
        {
            dihi += 2 * pi20;
        }
        dihi2 = dihi * dihi;
        // Update sum.
        sum_dihs += dihi2;
    }
    
    // Divide by angles summed.
    sum_angs   /= iatoms - 2;
    sum_dihs   /= iatoms - 3;
    // Geometric mean.
    sum_angdih  = sqrt(sqrt(sum_angs) * sqrt(sum_dihs));
    // Rescale from [0, pi] to [0, 1].
    sum_angdih /= pi20;
    // Output.
    return (real)sum_angdih;
}



real calc_angdih2g_n(int iatoms, rvec frame[], rvec rframe[], real angdih[])
{
    // Returns the root mean of the squared difference of angles.
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS, "Need at least 6 atoms in index group to use calc_angdih.");
    }
    // Initializing variables.
    int i;
    double pi20 = 3.14159265358979323846;
    double irang, irang2, dihi, dihi2, sum_angdih;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sum_angdih = 0.0;
    for (i = 0; i < iatoms; i++)
    {
        angdih[i] = 0.0;
    }
    // There are n - 3 dihedral angles.
    for (i = 3; i < iatoms; i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i - 3], frame[i - 2], vec1);
        rvec_sub(frame[i - 1], frame[i - 2], vec2);
        rvec_sub(frame[i - 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i - 3], rframe[i - 2], vec1);
        rvec_sub(rframe[i - 1], rframe[i - 2], vec2);
        rvec_sub(rframe[i - 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for the difference between the two angles.
        dihi  = (double)(iang - rang);
        // Rescale dihi from [-2.0 * pi, +2.0 * pi] to [-pi, +pi].
        if (dihi > pi20)
        {
            dihi -= 2 * pi20;
        }
        else if (dihi < -pi20)
        {
            dihi += 2 * pi20;
        }
        dihi2 = dihi * dihi;
        // Update output array.
        angdih[i - 2] += (real)dihi2;
        angdih[i - 1] += (real)dihi2;
    }
    
    // The dih portion of angdih stores the means of two dihedral angles.
    // *** First and last indeces = 0.0.
    // *** Second and second last indeces only include one dihedral angle.
    angdih[1]          = sqrt(angdih[1]);
    angdih[iatoms - 2] = sqrt(angdih[iatoms - 2]);
    for (i = 2; i < (iatoms - 2); i++)
    {
        angdih[i] = sqrt(angdih[i] / 2.0);
    }
    
    // There are n - 2 backbone angles.
    for (i = 1; i < (iatoms - 1); i++)
    {
        /* Find the angle between consecutive pairs of backbone atoms.
         */
        rvec_sub(frame[i - 1],  frame[i], vec1);
        rvec_sub(frame[i + 1],  frame[i], vec2);
        iang   = gmx_angle(vec1, vec2);
        
        rvec_sub(rframe[i - 1], rframe[i], vec1);
        rvec_sub(rframe[i + 1], rframe[i], vec2);
        rang   = gmx_angle(vec1, vec2);
        
        // Squared difference of angles.
        irang      = iang - rang;
        irang2     = irang * irang;
        angdih[i] *= (real)sqrt(irang2);
    }
    
    // Finalize output array.
    for (i = 1; i < (iatoms - 1); i++)
    {
        // Geometric mean of ang and dih.
        // Rescale from [0, pi] to [0, 1].
        angdih[i]   = sqrt(angdih[i]) / pi20;
        // To calculate the mean of the indeces.
        sum_angdih += (double)angdih[i];
    }
    // Divide by the number of summed indeces.
    return (real)(sum_angdih / (iatoms - 2));
}



real calc_rmsdih(int iatoms, rvec frame[], rvec rframe[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS, "Need at least 6 atoms in index group to use calc_rmsdih.");
    }
    // Initializing variables.
    int i;
    double pi20 = 3.14159265358979323846;
    double ir2, rr2, irrr, itheta, rtheta, dphi, thetaphi;
    rvec vec1, vec2, vec3, pvec1, pvec2;
//    double cosiphi, cosrphi, cosdphi;
    real iphi, rphi;
    
    double sum_rmsdih = 0.0;
    // There are n - 3 backbone dihedral angles.
    for (i = 3; i < iatoms; i++)
    {
        /* Find the value of r2.
         */
        ir2 = (double)distance2(frame[i - 1],  frame[i]);
        rr2 = (double)distance2(rframe[i - 1], rframe[i]);
        
        /* Find the value of theta.
         * 
         * itheta is the angle of two vectors made by three atom coords 
         * in frame. rtheta is the same for rframe.
         */
        rvec_sub(frame[i - 2],  frame[i - 1], vec1);
        rvec_sub(frame[i],      frame[i - 1], vec2);
        itheta = pi20 - (double)gmx_angle(vec1, vec2);
        
        rvec_sub(rframe[i - 2], rframe[i - 1], vec1);
        rvec_sub(rframe[i],     rframe[i - 1], vec2);
        rtheta = pi20 - (double)gmx_angle(vec1, vec2);
        
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i - 3], frame[i - 2], vec1);
        rvec_sub(frame[i - 1], frame[i - 2], vec2);
        rvec_sub(frame[i - 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iphi = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iphi *= -1.0;
        }
//        cosiphi = (double)cos_angle(pvec1, pvec2);
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i - 3], rframe[i - 2], vec1);
        rvec_sub(rframe[i - 1], rframe[i - 2], vec2);
        rvec_sub(rframe[i - 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rphi = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rphi *= -1.0;
        }
//        cosrphi = (double)cos_angle(pvec1, pvec2);
        
        // Solve for the difference between the two angles.
        dphi = (double)(iphi - rphi);
        // Rescale dihi from [-2.0 * pi, +2.0 * pi] to [-pi, +pi].
        if (dphi > pi20)
        {
            dphi -= 2 * pi20;
        }
        else if (dphi < -pi20)
        {
            dphi += 2 * pi20;
        }
//        cosdphi  = cosiphi * cosrphi;
//        cosdphi += sqrt(1 - cosiphi * cosiphi) * sqrt(1 - cosrphi * cosrphi);
        
        
        // Calculate the distance between two points in spherical coords.
        irrr        = 2 * sqrt(ir2) * sqrt(rr2);
        thetaphi    = sin(itheta) * sin(rtheta) * cos(dphi);
        thetaphi   += cos(itheta) * cos(rtheta);
        // Update the sum.
        sum_rmsdih += ir2 + rr2 - irrr * thetaphi;
    }
    
    
        // There are n - 3 backbone dihedral angles.
    for (i = 0; i < (iatoms - 3); i++)
    {
        /* Find the value of r2.
         */
        ir2 = (double)distance2(frame[i + 1],  frame[i]);
        rr2 = (double)distance2(rframe[i + 1], rframe[i]);
        
        /* Find the value of theta.
         * 
         * itheta is the angle of two vectors made by three atom coords 
         * in frame. rtheta is the same for rframe.
         */
        rvec_sub(frame[i + 2],  frame[i + 1], vec1);
        rvec_sub(frame[i],      frame[i + 1], vec2);
        itheta = pi20 - (double)gmx_angle(vec1, vec2);
        
        rvec_sub(rframe[i + 2], rframe[i + 1], vec1);
        rvec_sub(rframe[i],     rframe[i + 1], vec2);
        rtheta = pi20 - (double)gmx_angle(vec1, vec2);
        
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i + 3], frame[i + 2], vec1);
        rvec_sub(frame[i + 1], frame[i + 2], vec2);
        rvec_sub(frame[i + 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iphi = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iphi *= -1.0;
        }
//        cosiphi = (double)cos_angle(pvec1, pvec2);
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i + 3], rframe[i + 2], vec1);
        rvec_sub(rframe[i + 1], rframe[i + 2], vec2);
        rvec_sub(rframe[i + 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rphi = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rphi *= -1.0;
        }
//         cosrphi = (double)cos_angle(pvec1, pvec2);
        
        // Solve for the difference between the two angles.
        dphi = (double)(iphi - rphi);
        // Rescale dihi from [-2.0 * pi, +2.0 * pi] to [-pi, +pi].
        if (dphi > pi20)
        {
            dphi -= 2 * pi20;
        }
        else if (dphi < -pi20)
        {
            dphi += 2 * pi20;
        }
//         cosdphi  = cosiphi * cosrphi;
//         cosdphi += sqrt(1 - cosiphi * cosiphi) * sqrt(1 - cosrphi * cosrphi);
        
        
        // Calculate the distance between two points in spherical coords.
        irrr        = 2 * sqrt(ir2) * sqrt(rr2);
        thetaphi    = sin(itheta) * sin(rtheta) * cos(dphi);
        thetaphi   += cos(itheta) * cos(rtheta);
        // Update the sum.
        sum_rmsdih += ir2 + rr2 - irrr * thetaphi;
    }
    
    // Calculate the final value.
    sum_rmsdih /= 2 * (double)iatoms - 6;
    sum_rmsdih  = sqrt(sum_rmsdih);
    
    // Output.
    return (real)sum_rmsdih;
}



real calc_rmsdih_n(int iatoms, rvec frame[], rvec rframe[], real rmsdih[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS, "Need at least 6 atoms in index group to use calc_rmsdih.");
    }
    // Initializing variables.
    int i;
    double pi20 = 3.14159265358979323846;
    double ir2, rr2, irrr, itheta, rtheta, dphi, thetaphi;
    rvec vec1, vec2, vec3, pvec1, pvec2;
//    double cosiphi, cosrphi, cosdphi;
    real iphi, rphi;
    
    double sum_rmsdih = 0.0;
    for (i = 0; i < iatoms; i++)
    {
        rmsdih[i] = 0.0;
    }
    // There are n - 3 backbone dihedral angles.
    for (i = 3; i < iatoms; i++)
    {
        /* Find the value of r2.
         */
        ir2 = (double)distance2(frame[i - 1],  frame[i]);
        rr2 = (double)distance2(rframe[i - 1], rframe[i]);
        
        /* Find the value of theta.
         * 
         * itheta is the angle of two vectors made by three atom coords 
         * in frame. rtheta is the same for rframe.
         */
        rvec_sub(frame[i - 2],  frame[i - 1], vec1);
        rvec_sub(frame[i],      frame[i - 1], vec2);
        itheta = pi20 - (double)gmx_angle(vec1, vec2);
        
        rvec_sub(rframe[i - 2], rframe[i - 1], vec1);
        rvec_sub(rframe[i],     rframe[i - 1], vec2);
        rtheta = pi20 - (double)gmx_angle(vec1, vec2);
        
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i - 3], frame[i - 2], vec1);
        rvec_sub(frame[i - 1], frame[i - 2], vec2);
        rvec_sub(frame[i - 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iphi = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iphi *= -1.0;
        }
//        cosiphi = (double)cos_angle(pvec1, pvec2);
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i - 3], rframe[i - 2], vec1);
        rvec_sub(rframe[i - 1], rframe[i - 2], vec2);
        rvec_sub(rframe[i - 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rphi = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rphi *= -1.0;
        }
//        cosrphi = (double)cos_angle(pvec1, pvec2);
        
        // Solve for the difference between the two angles.
        dphi = (double)(iphi - rphi);
        // Rescale dihi from [-2.0 * pi, +2.0 * pi] to [-pi, +pi].
        if (dphi > pi20)
        {
            dphi -= 2 * pi20;
        }
        else if (dphi < -pi20)
        {
            dphi += 2 * pi20;
        }
//        cosdphi  = cosiphi * cosrphi;
//        cosdphi += sqrt(1 - cosiphi * cosiphi) * sqrt(1 - cosrphi * cosrphi);
        
        
        // Calculate the distance between two points in spherical coords.
        irrr        = 2 * sqrt(ir2) * sqrt(rr2);
        thetaphi    = sin(itheta) * sin(rtheta) * cos(dphi);
        thetaphi   += cos(itheta) * cos(rtheta);
        // Update the sum.
        rmsdih[i]  += ir2 + rr2 - irrr * thetaphi;
        sum_rmsdih += ir2 + rr2 - irrr * thetaphi;
    }
    
    
        // There are n - 3 backbone dihedral angles.
    for (i = 0; i < (iatoms - 3); i++)
    {
        /* Find the value of r2.
         */
        ir2 = (double)distance2(frame[i + 1],  frame[i]);
        rr2 = (double)distance2(rframe[i + 1], rframe[i]);
        
        /* Find the value of theta.
         * 
         * itheta is the angle of two vectors made by three atom coords 
         * in frame. rtheta is the same for rframe.
         */
        rvec_sub(frame[i + 2],  frame[i + 1], vec1);
        rvec_sub(frame[i],      frame[i + 1], vec2);
        itheta = pi20 - (double)gmx_angle(vec1, vec2);
        
        rvec_sub(rframe[i + 2], rframe[i + 1], vec1);
        rvec_sub(rframe[i],     rframe[i + 1], vec2);
        rtheta = pi20 - (double)gmx_angle(vec1, vec2);
        
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on each plane.
         * 
         * Find the angle between the two planes by calculating two
         * vectors normal to each plane.
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i + 3], frame[i + 2], vec1);
        rvec_sub(frame[i + 1], frame[i + 2], vec2);
        rvec_sub(frame[i + 1], frame[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iphi = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iphi *= -1.0;
        }
//         cosiphi = (double)cos_angle(pvec1, pvec2);
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i + 3], rframe[i + 2], vec1);
        rvec_sub(rframe[i + 1], rframe[i + 2], vec2);
        rvec_sub(rframe[i + 1], rframe[i],     vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rphi = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign. Result lies in range [-pi, pi].
        if (iprod(vec1, pvec2) < 0.0)
        {
            rphi *= -1.0;
        }
//         cosrphi = (double)cos_angle(pvec1, pvec2);
        
        // Solve for the difference between the two angles.
        dphi = (double)(iphi - rphi);
        // Rescale dihi from [-2.0 * pi, +2.0 * pi] to [-pi, +pi].
        if (dphi > pi20)
        {
            dphi -= 2 * pi20;
        }
        else if (dphi < -pi20)
        {
            dphi += 2 * pi20;
        }
//         cosdphi  = cosiphi * cosrphi;
//         cosdphi += sqrt(1 - cosiphi * cosiphi) * sqrt(1 - cosrphi * cosrphi);
        
        
        // Calculate the distance between two points in spherical coords.
        irrr        = 2 * sqrt(ir2) * sqrt(rr2);
        thetaphi    = sin(itheta) * sin(rtheta) * cos(dphi);
        thetaphi   += cos(itheta) * cos(rtheta);
        // Update the sum.
        rmsdih[i]  += ir2 + rr2 - irrr * thetaphi;
        sum_rmsdih += ir2 + rr2 - irrr * thetaphi;
    }
    
    rmsdih[0] = sqrt(rmsdih[0]);
    for (i = 1; i < (iatoms - 1); i++)
    {
        rmsdih[i] /= 2;
        rmsdih[i]  = sqrt(rmsdih[i]);
    }
    rmsdih[iatoms - 1] = sqrt(rmsdih[iatoms - 1]);
    
    // Calculate the final value.
    sum_rmsdih /= 2 * (double)iatoms - 6;
    sum_rmsdih  = sqrt(sum_rmsdih);
    
    // Output.
    return (real)sum_rmsdih;
}



real calc_phipsi(int iatoms, rvec frame[], rvec rframe[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS,"Need at least 6 atoms in index group to use calc_phipsi.");
    }
    
    // Initializing variables.
    int i, iC1 ,iC2, iN1, iN2, iCa;
    double cosxy, sum_phipsi;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sum_phipsi = 0.0;
    // Phi angles.
    for (i = 1; i < (iatoms / 3); i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on those planes.
         * 
         * Find the angle between the two planes by finding two
         * vectors normal to each plane. 
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // The atom indeces.
        iC1 = 3 * i - 1;
        iN1 = 3 * i;
        iCa = 3 * i + 1;
        iC2 = 3 * i + 2;
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[iC1], frame[iN1], vec1);
        rvec_sub(frame[iCa], frame[iN1], vec2);
        rvec_sub(frame[iCa], frame[iC2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[iC1], rframe[iN1], vec1);
        rvec_sub(rframe[iCa], rframe[iN1], vec2);
        rvec_sub(rframe[iCa], rframe[iC2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for cos(x - y).
        cosxy = cos((double)(iang - rang));
        // Update phipsi sum.
        sum_phipsi += cosxy;
    }
    
    // Psi angles.
    for (i = 0; i < ((iatoms / 3) - 1); i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on those planes.
         * 
         * Find the angle between the two planes by finding two
         * vectors normal to each plane. 
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // The atom indeces.
        iN1 = 3 * i;
        iCa = 3 * i + 1;
        iC1 = 3 * i + 2;
        iN2 = 3 * i + 3;
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[iN1], frame[iCa], vec1);
        rvec_sub(frame[iC1], frame[iCa], vec2);
        rvec_sub(frame[iC1], frame[iN2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[iN1], rframe[iCa], vec1);
        rvec_sub(rframe[iC1], rframe[iCa], vec2);
        rvec_sub(rframe[iC1], rframe[iN2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for cos(x - y).
        cosxy = cos((double)(iang - rang));
        // Update phipsi sum.
        sum_phipsi += cosxy;
    }
    
    // Divide by the number of angles summed.
    sum_phipsi /= 2 * ((iatoms / 3) - 1);
    // Rescale from [+1, -1] to [0, 1].
    sum_phipsi  = (sum_phipsi - 1.0) / (-2.0);
    return (real)sum_phipsi;
}



real calc_phipsi_n(int iatoms, rvec frame[], rvec rframe[], real phipsi[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS,"Need at least 6 atoms in index group to use calc_phipsi.");
    }
    
    // Initializing variables.
    int i, iC1 ,iC2, iN1, iN2, iCa;
    double cosxy, sum_phipsi;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sum_phipsi = 0.0;
    for (i = 0; i < (iatoms / 3) - 1; i++)
    {
        phipsi[i] = 0.0;
    }
    // Phi angles.
    for (i = 1; i < (iatoms / 3); i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on those planes.
         * 
         * Find the angle between the two planes by finding two
         * vectors normal to each plane. 
         * 
         * Most tricks only return an angle between 0 and
         * pi rather than the full 2*pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // The atom indeces.
        iC1 = 3 * i - 1;
        iN1 = 3 * i;
        iCa = 3 * i + 1;
        iC2 = 3 * i + 2;
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[iC1], frame[iN1], vec1);
        rvec_sub(frame[iCa], frame[iN1], vec2);
        rvec_sub(frame[iCa], frame[iC2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[iC1], rframe[iN1], vec1);
        rvec_sub(rframe[iCa], rframe[iN1], vec2);
        rvec_sub(rframe[iCa], rframe[iC2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for the dot product and normalize.
        cosxy = (real)(cos(((double)iang) - ((double)rang)) - 1.0) / (-2.0);
        // Add sum to phipsi.
        if (cosxy < 0.0)
        {
            phipsi[i] = 0.0;
        }
        else
        {
            phipsi[i] = cosxy;
        }
    }
    
    // Psi angles.
    for (i = 0; i < ((iatoms / 3) - 1); i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on those planes.
         * 
         * Find the angle between the two planes by finding two
         * vectors normal to each plane. 
         * 
         * Most tricks only return an angle between 0 and
         * pi rather than the full 2*pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // The atom indeces.
        iN1 = 3 * i;
        iCa = 3 * i + 1;
        iC1 = 3 * i + 2;
        iN2 = 3 * i + 3;
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[iN1], frame[iCa], vec1);
        rvec_sub(frame[iC1], frame[iCa], vec2);
        rvec_sub(frame[iC1], frame[iN2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[iN1], rframe[iCa], vec1);
        rvec_sub(rframe[iC1], rframe[iCa], vec2);
        rvec_sub(rframe[iC1], rframe[iN2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for the dot product and normalize range: 0 to 1.
        cosxy = (real)(cos(((double)iang) - ((double)rang)) - 1.0)/(-2.0);
        // Update phipsi.
        if (i == 0) // First residue has no phi angle.
        {
            if (cosxy < 0.0)
            {
                phipsi[0] = 0.0;
            }
            else
            {
                phipsi[0] = cosxy;
            }
        }
        else
        {
            if (cosxy < 0.0)
            {
                phipsi[i] = phipsi[i] / 2.0;
            }
            else
            {
                phipsi[i] = (phipsi[i] + ((real)cosxy)) / 2.0;
            }
        }
    }
    for (i = 0; i < iatoms; i++)
    {
        sum_phipsi += (double)phipsi[i];
    }
    return (real)(sum_phipsi / iatoms);
}



real calc_phipsi2(int iatoms, rvec frame[], rvec rframe[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS,"Need at least 6 atoms in index group to use calc_phipsi.");
    }
    
    // Initializing variables.
    int i, iC1 ,iC2, iN1, iN2, iCa;
    double irang, irang2, sum_phipsi;
    double pi20 = 3.14159265358979323846;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sum_phipsi = 0.0;
    // Phi angles.
    for (i = 1; i < (iatoms / 3); i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on those planes.
         * 
         * Find the angle between the two planes by finding two
         * vectors normal to each plane. 
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // The atom indeces.
        iC1 = 3 * i - 1;
        iN1 = 3 * i;
        iCa = 3 * i + 1;
        iC2 = 3 * i + 2;
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[iC1], frame[iN1], vec1);
        rvec_sub(frame[iCa], frame[iN1], vec2);
        rvec_sub(frame[iCa], frame[iC2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[iC1], rframe[iN1], vec1);
        rvec_sub(rframe[iCa], rframe[iN1], vec2);
        rvec_sub(rframe[iCa], rframe[iC2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Square the difference.
        irang  = (double)(iang - rang);
        if (irang >  pi20)
        {
            irang -= 2 * pi20;
        }
        if (irang < -pi20)
        {
            irang += 2 * pi20;
        }
        irang2 = irang * irang;
        // Update phipsi sum.
        sum_phipsi += irang2;
    }
    
    // Psi angles.
    for (i = 0; i < ((iatoms / 3) - 1); i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on those planes.
         * 
         * Find the angle between the two planes by finding two
         * vectors normal to each plane. 
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // The atom indeces.
        iN1 = 3 * i;
        iCa = 3 * i + 1;
        iC1 = 3 * i + 2;
        iN2 = 3 * i + 3;
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[iN1], frame[iCa], vec1);
        rvec_sub(frame[iC1], frame[iCa], vec2);
        rvec_sub(frame[iC1], frame[iN2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[iN1], rframe[iCa], vec1);
        rvec_sub(rframe[iC1], rframe[iCa], vec2);
        rvec_sub(rframe[iC1], rframe[iN2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Square of difference of angles limited to range [-pi, +pi].
        irang  = (double)(iang - rang);
        if (irang >  pi20)
        {
            irang -= 2 * pi20;
        }
        if (irang < -pi20)
        {
            irang += 2 * pi20;
        }
        irang2 = irang * irang;
        // Update phipsi sum.
        sum_phipsi += irang2;
    }
    
    // Divide by the number of angles summed and sqrt.
    sum_phipsi  = sqrt(sum_phipsi / (2 * ((iatoms / 3) - 1)));
    // Rescale from [0, pi] to [0, 1].
    sum_phipsi /= pi20;
    return (real)sum_phipsi;
}



real calc_phipsi2_n(int iatoms, rvec frame[], rvec rframe[], real phipsi[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS,"Need at least 6 atoms in index group to use calc_phipsi.");
    }
    
    // Initializing variables.
    int i, iC1 ,iC2, iN1, iN2, iCa;
    double irang, irang2, sum_phipsi;
    double pi20 = 3.14159265358979323846;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    sum_phipsi = 0.0;
    for (i = 0; i < (iatoms / 3); i++)
    {
        phipsi[i] = 0.0;
    }
    // Phi angles.
    for (i = 1; i < (iatoms / 3); i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on those planes.
         * 
         * Find the angle between the two planes by finding two
         * vectors normal to each plane. 
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // The atom indeces.
        iC1 = 3 * i - 1;
        iN1 = 3 * i;
        iCa = 3 * i + 1;
        iC2 = 3 * i + 2;
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[iC1], frame[iN1], vec1);
        rvec_sub(frame[iCa], frame[iN1], vec2);
        rvec_sub(frame[iCa], frame[iC2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[iC1], rframe[iN1], vec1);
        rvec_sub(rframe[iCa], rframe[iN1], vec2);
        rvec_sub(rframe[iCa], rframe[iC2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Square the difference.
        irang  = (double)(iang - rang);
        if (irang >  pi20)
        {
            irang -= 2 * pi20;
        }
        if (irang < -pi20)
        {
            irang += 2 * pi20;
        }
        irang2    = irang * irang;
        phipsi[i] = irang2;
    }
    
    // Psi angles.
    for (i = 0; i < ((iatoms / 3) - 1); i++)
    {
        /* Use four atom coordinates to make two planes defined
         * by three atom coordinates on those planes.
         * 
         * Find the angle between the two planes by finding two
         * vectors normal to each plane. 
         * 
         * Most algorithms only return an angle between 0 and
         * pi rather than the full 2 * pi possible range. The
         * method of using the sign of the inner product of
         * the vector normal to one plane with a vector on the
         * other plane comes from bondfree.c in the gromacs
         * source.
         * 
         * Note that since the result is the difference of two
         * angles, the directionality of the sign is unimportant
         * as long as it is consistent.
         */
        // The atom indeces.
        iN1 = 3 * i;
        iCa = 3 * i + 1;
        iC1 = 3 * i + 2;
        iN2 = 3 * i + 3;
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[iN1], frame[iCa], vec1);
        rvec_sub(frame[iC1], frame[iCa], vec2);
        rvec_sub(frame[iC1], frame[iN2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[iN1], rframe[iCa], vec1);
        rvec_sub(rframe[iC1], rframe[iCa], vec2);
        rvec_sub(rframe[iC1], rframe[iN2], vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1, pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1, pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Square of difference of angles limited to range [-pi, +pi].
        irang  = (double)(iang - rang);
        if (irang >  pi20)
        {
            irang -= 2 * pi20;
        }
        if (irang < -pi20)
        {
            irang += 2 * pi20;
        }
        irang2      = irang * irang;
        phipsi[i]  += irang2;
    }
    for (i = 1; i < ((iatoms / 3) - 1); i++)
    {
        phipsi[i] /= 2.0;
    }
    for (i = 0; i < (iatoms / 3); i++)
    {
        phipsi[i]   = sqrt(phipsi[i]) / pi20;
        sum_phipsi += phipsi[i];
    }
    // Divide by number of summed indeces.
    return (real)(sum_phipsi / (iatoms / 3));
}



real calc_mmsd(int iatoms, rvec frame[], rvec rframe[])
{
    /*
     * MSD between fit structure and mirror of reference structure.
     * 
     * Assumes they have already been aligned.
     */
    int i, d;
    real xd, xdm, mmsdi;
    double mmsdt;
    
    // Initialize to zero.
    mmsdt = 0;
    
    // Main loop.
    for(i = 0; i < iatoms; i++)
    {
        // Copy over coordinates from xref and rotate them.
        for(d = 0; d < 3; d++)
        {
            /* Same as multiplication by (-1,-1,-1) diag matrix.
             * 
             * d = 0: x coordinate, etc
             * 
             * xdm = distances from fit coord to mirrored ref coord.
             */
            xdm    = frame[i][d] + rframe[i][d];
            // Add to the sum.
            mmsdt += xdm * xdm;
        }
    }
    
    // Normalize for number of atoms.
    return (real)sqrt(mmsdt / iatoms);
}



real calc_corr_atoms(int iatoms,rvec frame[],rvec rframe[])
{   
    int i,d;
    rvec corr_numer,corr_denom,corr_denmr;
    real corr_atoms;
    
    // Initialize corr_numer and corr_denom to zero.
    corr_atoms = 0;
    for (d=0;d<3;d++)
    {
        corr_numer[d] = 0;
        corr_denom[d] = 0;
        corr_denmr[d] = 0;
    }
    
    // Main loop.
    for (i=0;i<iatoms;i++)
    {
        for (d=0;d<3;d++)
        {
            corr_numer[d] += frame[i][d]  * rframe[i][d];
            corr_denom[d] += frame[i][d]  * frame[i][d];
            corr_denmr[d] += rframe[i][d] * rframe[i][d];
        }
    }
    
    // Solve output correlation.
    for (d=0;d<3;d++)
    {
        corr_atoms += corr_numer[d]/(sqrt(corr_denom[d])*sqrt(corr_denmr[d]));
    }
    corr_atoms /= 3;
    corr_atoms = 1 - corr_atoms;
    
    if (corr_atoms > 1)
    {
        corr_atoms = 1;
    }
    
    return corr_atoms;
}



real calc_corr_angs(int iatoms,rvec frame[],rvec rframe[])
{
    // Error checking.
    if (iatoms < 4)
    {
        // Not sure if this works. Depends on how FARGS works.
        gmx_fatal(FARGS,"Need at least 4 atoms in index group to use calc_corr_angs.");
    }
    
    // Initializing variables.
    int i;
    real corr_numer,corr_denom,corr_denmr,corr_angs;
    rvec vec1,vec2;
    real meanf,meanr;
    
    real *cosx,*cosy;
    snew(cosx,iatoms);
    snew(cosy,iatoms);
    
    // Initialize variables to zero.
    corr_angs  = 0;
    corr_numer = 0;
    corr_denom = 0;
    corr_denmr = 0;
    meanf = 0;
    meanr = 0;
    for (i=0; i<iatoms; i++)
    {
        cosx[i] = 0;
        cosy[i] = 0;
    }
    
    
    // Loop from the third atom to the last atom. Forward direction.
    for (i=2; i < iatoms; i++)
    {
        /* Find the value of cosx and cosy.
         * 
         * x is the angle of two vectors made by three atom coords
         * in frame[]. y is the same for rframe[].
         */
        rvec_sub(frame[i-2],frame[i-1],vec1);
        rvec_sub(frame[i],  frame[i-1],vec2);
        cosx[i] += iprod(vec1,vec2)/(norm(vec1)*norm(vec2));
        
        rvec_sub(rframe[i-2],rframe[i-1],vec1);
        rvec_sub(rframe[i],  rframe[i-1],vec2);
        cosy[i] += iprod(vec1,vec2)/(norm(vec1)*norm(vec2));
    }
    
/* This may eventually be used for per atom results, but it would require
 * finding the correlation over all frames, so it requires some changes to
 * other source code files.
 */ 
    
//     // Loop from third last atom to the first atom. Backward direction.
//     for (i=(iatoms-3); i >= 0; i--)
//     {
//         /* Find the value of cosx and cosy.
//          * 
//          * x is the angle of two vectors made by three atom coords
//          * in frame[]. y is the same for rframe[].
//          */
//         rvec_sub(frame[i+2],frame[i+1],vec1);
//         rvec_sub(frame[i],  frame[i+1],vec2);
//         cosx[i] += iprod(vec1,vec2)/(norm(vec1)*norm(vec2));
//         
//         rvec_sub(rframe[i+2],rframe[i+1],vec1);
//         rvec_sub(rframe[i],  rframe[i+1],vec2);
//         cosy[i] += iprod(vec1,vec2)/(norm(vec1)*norm(vec2));
//     }
//     
//     // The first and last two angles do not have double measurements.
//     // Divide all other angles by two.
//     for (i=2; i<(iatoms-2); i++)
//     {
//         cosx[i] /= 2.0;
//         cosy[i] /= 2.0;
//     }
    
    // Now find the means.
    for (i=0; i<iatoms; i++)
    {
        meanf += cosx[i];
        meanr += cosy[i];
    }
    // Cast iatoms to (real) here?
    meanf /= iatoms;
    meanr /= iatoms;
    
    // Solve for correlation coefficient.
    for (i=0; i<iatoms; i++)
    {
        corr_numer += (cosx[i] - meanf) * (cosy[i] - meanr);
        corr_denom += (cosx[i] - meanf) * (cosx[i] - meanf);
        corr_denmr += (cosy[i] - meanr) * (cosy[i] - meanr);
    }
    // Can't remember if there is a reason why I did it this way. Really odd.
    // Also, should check for divide by zero errors just in case.
    corr_angs  = corr_numer;
    corr_angs /= (real)(sqrt(corr_denom) * sqrt(corr_denmr));
    
    // Free memory.
    sfree(cosx);
    sfree(cosy);
    
    // Only keeping positive correlation information.
    // Also checking for rounding errors.
    if ((1.0 - corr_angs) < 0.0)
    {
        return 0.0;
    }
    else if ((1.0 - corr_angs) > 1.0)
    {
        return 1.0;
    }
    else
    {
        return (1.0 - corr_angs);
    }
}



real calc_rg(int iatoms, rvec frame[])
{
    // Assumes frame[] to be centered.
    int i, d;
    double rg = 0;
    for (i = 0; i < iatoms; i++)
    {
        for (d = 0; d < 3; d++)
        {
            rg += frame[i][d] * frame[i][d];
        }
    }
    return (real)sqrt(rg / iatoms);
}



real calc_rrot(int iatoms, rvec frame[], rvec rframe[])
{
    int    i, j, k;
    double rx, ry, rz;
    double pi20 = 3.14159265358979323846;
    rvec   xold;
    matrix rotx, roty, rotz, rotm;
    
    // Solve for three random numbers.
    rx = 2.0 * pi20 * rand() / RAND_MAX - pi20;
    ry = 2.0 * pi20 * rand() / RAND_MAX - pi20;
    rz = 2.0 * pi20 * rand() / RAND_MAX - pi20;
    // Reset matrices to zero.
    clear_mat(rotx);
    clear_mat(roty);
    clear_mat(rotz);
    
    /*      Rx = rotx[rows][cols]
     * 
     *      |   1.0   |   0.0   |   0.0   |
     * Rx = |   0.0   |  cos(x) |  sin(x) |
     *      |   0.0   | -sin(x) |  cos(x) |
     */
    rotx[0][0] = 1.0;
    rotx[1][1] = (real)cos(rx);
    rotx[2][2] = rotx[1][1];
    rotx[1][2] = (real)sin(rx);
    rotx[2][1] = -1.0 * rotx[1][2];
    
    /*      Ry = roty[rows][cols]
     * 
     *      |  cos(x) |   0.0   | -sin(x) |
     * Ry = |   0.0   |   1.0   |   0.0   |
     *      |  sin(x) |   0.0   |  cos(x) |
     */
    roty[1][1] = 1.0;
    roty[0][0] = (real)cos(ry);
    roty[2][2] = roty[0][0];
    roty[2][0] = (real)sin(ry);
    roty[0][2] = -1.0 * roty[2][0];
    
    /*      Rz = rotz[rows][cols]
     * 
     *      |  cos(x) |  sin(x) |   0.0   |
     * Rz = | -sin(x) |  cos(x) |   0.0   |
     *      |   0.0   |   0.0   |   1.0   |
     */
    rotz[2][2] = 1.0;
    rotz[0][0] = (real)cos(rz);
    rotz[1][1] = rotz[0][0];
    rotz[0][1] = (real)sin(rz);
    rotz[1][0] = -1.0 * rotz[0][1];
    
    // Multiply rotation matrices.
    mmul(rotx, roty, rotm);
    copy_mat(rotm, rotx);
    mmul(rotx, rotz, rotm);
    // Apply random rotation.
    for (k = 0; k < iatoms; k++)
    {
        for (i = 0; i < 3; i++)
        {
            xold[i] = rframe[k][i];
        }
        for (i = 0; i < 3; i++)
        {
            rframe[k][i] = 0.0;
            for (j = 0; j < 3; j++)
            {
                rframe[k][i] += rotm[i][j] * xold[j];
            }
        }
    }
    // Calculate RMSD after rotation.
    return calc_rmsd(iatoms, frame, rframe);
}



real calc_mammoth(int iatoms,rvec frame[],rvec rframe[],int* rnum)
{
    return 0.0; //temporary, should be removed by patch
//     /* This is intended to be a wrapper function to call mammoth.
//      * 
//      * The input data needs to be reformatted then fed to the fortran
//      * function for mammoth.
//      */
//     double zscore = 0.0;
//     double norm = (double)iatoms;
//     double zmax;
//     double *xfit,*xref;
//     int i,j;
//     
//     snew(xfit,3*iatoms);
//     snew(xref,3*iatoms);
//     
//     for (i=0; i<iatoms; i++)
//     {
//         for (j=0; j<3; j++)
//         {
//             // Convert nm to A and float to double.
//             xfit[i*3 + j] = (double)(10.0 * frame[i][j]);
//             xref[i*3 + j] = (double)(10.0 * rframe[i][j]);
//         }
//     }
//     
//     mammothmod_(xfit,xref,&iatoms,rnum,&zscore);
//     
//     // I could save a few cycles by doing this outside of the function. Meh.
//     zmax = (0.80006)*(pow(norm,0.6882)) - (5.9788)*(pow(norm,-0.1089));
//     
//     //We want to return (1 - zscore/zmax).
//     zscore /= zmax;
//     zscore = 1.0 - zscore;
//     
//     sfree(xfit);
//     sfree(xref);
//     
//     // Check for rounding errors.
//     if (zscore > 1.0)
//     {
//         return 1.0;
//     }
//     else if (zscore < 0.0)
//     {
//         return 0.0;
//     }
//     else
//     {
//         return (real)zscore;
//     }
}



real calc_esa(int iatoms,rvec frame[],rvec rframe[])
{
    // For larger chains of atoms, we need to change the 200.
    // It should become a command line option at some point.
    // Everything is done in external library, libesa.
    return (real)esa_analysis(iatoms,frame,rframe,200);
}



real call_ISDM(int iatoms, rvec cframe[], rvec rframe[], const char *ISDM)
{
    // Default behavior (no -ISDM option) is RMSD.
    if (strcmp(ISDM, "RMSD") == 0)
    {
        // Calculate RMSD.
        return calc_rmsd(iatoms, cframe, rframe);
    }
    
    // Scaled RMSD. User gives -srms option.
    if (strcmp(ISDM, "SRMS") == 0)
    {
        // Calculate scaled RMSD.
        return calc_srms(iatoms, cframe, rframe);
    }
    
    // Mirror scaled RMSD. User gives -mrms option.
    if (strcmp(ISDM, "MRMS") == 0)
    {
        // Calculate scaled RMSD.
        return calc_srms(iatoms, cframe, rframe);
    }
    
    // Difference of Rg. User gives -rg option.
    if (strcmp(ISDM, "RG") == 0)
    {
        real rg, rgr, ISD;
        // Rg of jth frame.
        rg  = calc_rg(iatoms, cframe);
        // Rg of ith frame.
        rgr = calc_rg(iatoms, rframe);
        // Find difference.
        ISD = rg - rgr;
        // Find absolute value.
        return sqrt(ISD * ISD);
    }
    
    // Scaled difference of Rg. User gives -srg option.
    if (strcmp(ISDM, "SRG") == 0)
    {
        real rg, rgr, ISD;
        // Rg of jth frame.
        rg  = calc_rg(iatoms, cframe);
        // Subtract Rg of ith frame.
        rgr = calc_rg(iatoms, rframe);
        // Find difference.
        ISD = rg - rgr;
        // Find scaled absolute value.
        ISD  = 2.0 * sqrt(ISD * ISD) / (rg + rgr);
        return ISD;
    }
    
    //  Mean Rg of the two input frames. Used internally.
    if (strcmp(ISDM, "MRG") == 0)
    {
        real rg, rgr, ISD;
        // Rg of jth frame.
        rg  = calc_rg(iatoms, cframe);
        // Rg of ith frame.
        rgr = calc_rg(iatoms, rframe);
        // Find mean.
        ISD = rg + rgr / 2.0;
        return ISD;
    }
    
    //  Geometric mean Rg of the two input frames.
    if (strcmp(ISDM, "GMRG") == 0)
    {
        real rg, rgr, ISD;
        // Rg of jth frame.
        rg  = calc_rg(iatoms, cframe);
        // Rg of ith frame.
        rgr = calc_rg(iatoms, rframe);
        // Geometric mean.
        ISD = sqrt(rg * rg + rgr * rgr);
        return ISD;
    }
    
    // Difference of end-to-end distance. User gives -e2e option.
    if (strcmp(ISDM, "E2E") == 0)
    {
        real ISD;
        // Find difference of end-to-end distances.
        ISD  = (real)sqrt(distance2(cframe[(iatoms - 1)], cframe[0]));
        ISD -= (real)sqrt(distance2(rframe[(iatoms - 1)], rframe[0]));
        // Find absolute value.
        return sqrt(ISD * ISD);
    }
    
    // Scaled difference of end-to-end distance. User gives -se2e option.
    if (strcmp(ISDM, "SE2E") == 0)
    {
        real rg, rgr, ISD;
        // Rg of jth frame.
        rg  = calc_rg(iatoms, cframe);
        // Rg of ith frame.
        rgr = calc_rg(iatoms, rframe);
        // Find difference of end-to-end distances.
        ISD  = (real)sqrt(distance2(cframe[(iatoms - 1)], cframe[0]));
        ISD -= (real)sqrt(distance2(rframe[(iatoms - 1)], rframe[0]));
        // Find scaled absolute value.
        ISD  = sqrt(ISD * ISD) / (rg + rgr);
        if (ISD > 1.0)
        {
            return 1.0;
        }
        return ISD;
    }
    
    // Mirrored RMSD. User gives -mir option.
    if (strcmp(ISDM, "MIR") == 0)
    {
        // Calculate mirrored RMSD.
        return (real)calc_mmsd(iatoms, cframe, rframe);
    }
    
    // Angles. User gives -ang option.
    if (strcmp(ISDM, "ANG") == 0)
    {
        // Calculate dot product of differences in angles.
        return calc_ang(iatoms, cframe, rframe);
    }
    
    // Angles. User gives -ang option.
    if (strcmp(ISDM, "ANG2") == 0)
    {
        // Calculate differences of backbone angles.
        return calc_ang2(iatoms, cframe, rframe);
    }
    
    // Dihedrals. User gives -dih option.
    if (strcmp(ISDM, "DIH") == 0)
    {
        // Calculate dot product of difference in dihedrals.
        return calc_dih(iatoms, cframe, rframe);
    }
    
    // Dihedrals. User gives -dih option.
    if (strcmp(ISDM, "DIH2") == 0)
    {
        // Calculate differences of dihedrals.
        return calc_dih2(iatoms, cframe, rframe);
    }
    
    // Angles and dihedrals. User gives -angdih option.
    if (strcmp(ISDM, "ANGDIH") == 0)
    {
        // Combination of ang and dih.
        return calc_angdih(iatoms, cframe, rframe);
    }
    
    // Combination of ang and dih. User gives -angdih2 option.
    if (strcmp(ISDM, "ANGDIH2") == 0)
    {
        // Combination of ang and dih.
        return calc_angdih2(iatoms, cframe, rframe, 2.0);
    }
    
    // Combination of ang and dih. User gives -angdih2 option.
    if (strcmp(ISDM, "ANGDIH2G") == 0)
    {
        // Combination of ang and dih.
        return calc_angdih2g(iatoms, cframe, rframe);
    }
    
    // Combination of ang and dih. User gives -angdih2m option.
    if (strcmp(ISDM, "ANGDIH2M") == 0)
    {
        // Combination of ang and dih.
        return calc_angdih2(iatoms, cframe, rframe, 10.0);
    }
    
    // Angles and dihedrals. User gives -angdih option.
    if (strcmp(ISDM, "RMSDIH") == 0)
    {
        // Combination of ang and dih.
        return calc_rmsdih(iatoms, cframe, rframe);
    }
    
    // Phi psi angles. User gives -phipsi option.
    if (strcmp(ISDM, "PHIPSI") == 0)
    {
        // Calculate dot product of difference in dihedrals.
        return calc_phipsi(iatoms, cframe, rframe);
    }
    
    // Phi psi angles. User gives -phipsi option.
    if (strcmp(ISDM, "PHIPSI2") == 0)
    {
        // Calculate dot product of difference in dihedrals.
        return calc_phipsi2(iatoms, cframe, rframe);
    }
    
    // Atom to atom distances. User gives -drms option.
    if (strcmp(ISDM, "DRMS") == 0)
    {
        // Calculate differences in all atom to atom distances.
        return calc_drms(iatoms, cframe, rframe);
    }
    
    // Atom to atom distances. User gives -sdrms option.
    if (strcmp(ISDM, "SDRMS") == 0)
    {
        /* Calculate differences in all atom to atom distances.
         * 
         * Scaled by Rgi + Rgj.
         * 
         * A few cycles could be saved by storing Rgi, but most
         * necessary info to calculate Rgi needs to be calculated
         * for every comparison anyway.
         */
        return calc_sdrms(iatoms, cframe, rframe);
    }
    
    // Position correlation. User gives -pcor option.
    if (strcmp(ISDM, "PCOR") == 0)
    {
        // Calculate position correlation.
        return calc_corr_atoms(iatoms,cframe,rframe);
    }
    
    // Backbone angle correlation. User gives -acor option.
    if (strcmp(ISDM, "ACOR") == 0)
    {
        // Calculate backbone angle correlation.
        return calc_corr_angs(iatoms,cframe,rframe);
    }
    
    if (strcmp(ISDM, "RROT") == 0)
    {
        // Randomly rotated RMSD.
        return calc_rrot(iatoms, cframe, rframe);
    }
}


real call_ISDM_n(int iatoms, rvec cframe[], rvec rframe[], real rISD[], const char *ISDM)
{
    // Default behavior (no -ISDM option) is RMSD.
    if (strcmp(ISDM, "RMSD") == 0)
    {
        // Calculate RMSD.
        return calc_rmsd_n(iatoms, cframe, rframe, rISD);
    }
    
    // Scaled RMSD. User gives -srms option.
    if (strcmp(ISDM, "SRMS") == 0)
    {
        // Calculate scaled RMSD.
        return calc_srms_n(iatoms, cframe, rframe, rISD);
    }
    
    // Mirror scaled RMSD. User gives -mrms option.
    if (strcmp(ISDM, "MRMS") == 0)
    {
        // Calculate scaled RMSD.
        return calc_srms_n(iatoms, cframe, rframe, rISD);
    }
    
    // Angles. User gives -ang option.
    if (strcmp(ISDM, "ANG") == 0)
    {
        // Calculate dot product of differences in angles.
        return calc_ang_n(iatoms, cframe, rframe, rISD);
    }
    
    // Angles. User gives -ang option.
    if (strcmp(ISDM, "ANG2") == 0)
    {
        // Calculate differences of backbone angles.
        return calc_ang2_n(iatoms, cframe, rframe, rISD);
    }
    
    // Dihedrals. User gives -dih option.
    if (strcmp(ISDM, "DIH") == 0)
    {
        // Calculate dot product of difference in dihedrals.
        return calc_dih_n(iatoms, cframe, rframe, rISD);
    }
    
    // Dihedrals. User gives -dih option.
    if (strcmp(ISDM, "DIH2") == 0)
    {
        // Calculate differences of dihedrals.
        return calc_dih2_n(iatoms, cframe, rframe, rISD);
    }
    
    // Angles and dihedrals. User gives -angdih option.
    if (strcmp(ISDM, "ANGDIH") == 0)
    {
        // Combination of ang and dih.
        return calc_angdih_n(iatoms, cframe, rframe, rISD);
    }
    
    // Combination of ang and dih. User gives -angdih2 option.
    if (strcmp(ISDM, "ANGDIH2") == 0)
    {
        // Combination of ang and dih.
        return calc_angdih2_n(iatoms, cframe, rframe, rISD, 2.0);
    }
    
    // Combination of ang and dih. User gives -angdih2 option.
    if (strcmp(ISDM, "ANGDIH2G") == 0)
    {
        // Combination of ang and dih.
        return calc_angdih2g_n(iatoms, cframe, rframe, rISD);
    }
    
    // Combination of ang and dih. User gives -angdih2m option.
    if (strcmp(ISDM, "ANGDIH2M") == 0)
    {
        // Combination of ang and dih.
        return calc_angdih2_n(iatoms, cframe, rframe, rISD, 10.0);
    }
    
    // Combination of ang and dih. User gives -angdih2 option.
    if (strcmp(ISDM, "RMSDIH") == 0)
    {
        // Combination of ang and dih.
        return calc_rmsdih_n(iatoms, cframe, rframe, rISD);
    }
    
    // Phi psi angles. User gives -phipsi option.
    if (strcmp(ISDM, "PHIPSI") == 0)
    {
        // Calculate dot product of difference in dihedrals.
        return calc_phipsi_n(iatoms, cframe, rframe, rISD);
    }
    
    // Phi psi angles. User gives -phipsi option.
    if (strcmp(ISDM, "PHIPSI2") == 0)
    {
        // Calculate dot product of difference in dihedrals.
        return calc_phipsi2_n(iatoms, cframe, rframe, rISD);
    }
    
    // Atom to atom distances. User gives -drms option.
    if (strcmp(ISDM, "DRMS") == 0)
    {
        // Calculate differences in all atom to atom distances.
        return calc_drms_n(iatoms, cframe, rframe, rISD);
    }
    
    // Atom to atom distances. User gives -sdrms option.
    if (strcmp(ISDM, "SDRMS") == 0)
    {
        /* Calculate differences in all atom to atom distances.
         * 
         * Scaled by Rgi + Rgj.
         * 
         * A few cycles could be saved by storing Rgi, but most
         * necessary info to calculate Rgi needs to be calculated
         * for every comparison anyway.
         */
        return calc_sdrms_n(iatoms, cframe, rframe, rISD);
    }
}
