/* Modified 04 Oct 2013
 * 
 * Tim Connolly - tconnolly@ucmerced.edu
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



//extern void mammothmod_(double xfit[],double xref[],int* nc,int* rnum,double* zscore);
double esa_analysis(int iatoms,rvec frame[],rvec rframe[],int sample);



real calc_drms(int iatoms,rvec frame[],rvec rframe[],real drms[])
{
    int i,j;
    real ileg,rleg,diff_leg,sum_dist,Rgi,Rgr,Rg2;
    
    // Don't assume that drms is zeros already.
    for (i=0; i < iatoms; i++)
    {
        drms[i] = 0.0;
    }
    
    //Normalization
    Rgi = 0.0;
    Rgr = 0.0;
    // Final sum.
    sum_dist = 0.0;
    // Loop through each atom.
    for (i=0; i < iatoms; i++)
    {
        // Compare to every other atom.
        for (j=0; j < iatoms; j++)
        {
            // It makes sense to skip i==j, but it is zero and is simpler to leave.
            
            // Find the length of the vector made by every pair of atoms.
            ileg = distance2(frame[i],frame[j]);
            rleg = distance2(rframe[i],rframe[j]);
            // Normalization
            Rgi += ileg;
            Rgr += rleg;
            // Find difference in the distances between the pair of atoms.
            diff_leg = sqrt(ileg) - sqrt(rleg);
            // Update the sum by squaring the difference.
            drms[i] += diff_leg * diff_leg;
        }
        // Normalize the sum by the number of differences.
        drms[i] /= (iatoms-1);
        // Takes the square root.
        drms[i] = sqrt(drms[i]);
        // Update the main sum.
        sum_dist += drms[i];
    }
    
    // Solves for Rg which is used for normalization.
    Rgi = sqrt(Rgi / (2 * iatoms * iatoms));
    Rgr = sqrt(Rgr / (2 * iatoms * iatoms));
    Rg2 = 2 * sqrt(Rgi*Rgr);
    // First, normalize by the number of summed differences.
    sum_dist /= iatoms;
    // Second, scale by the molecule size. Uses twice the geometric mean of Rg's.
    sum_dist /= Rg2;
    // Expected output is between 0 and 1, but the scaling allows for a small
    // chance of an output being >1. This would be undesirable.
    if (sum_dist > 1.0)
    {
        sum_dist = 1.0;
    }
    // Repeat per atom.
    for (i=0; i<iatoms; i++)
    {
        drms[i] /= Rg2;
        if (drms[i] > 1.0)
        {
            drms[i] = 1.0;
        }
    }
    // Output.
    return sum_dist;
}



real calc_rg(int iatoms, rvec frame[])
{
    int i, d;
    real rg = 0;
    for (i = 0; i < iatoms; i++)
    {
        for (d = 0; d < 3; d++)
        {
            rg += frame[i][d] * frame[i][d];
        }
    }
    return (real)sqrt(rg / iatoms);
}



real calc_ang(int iatoms, rvec frame[], rvec rframe[], real ang[])
{
    // Error checking.
    if (iatoms < 4)
    {
        // Not sure if this works. Depends on how FARGS works.
        gmx_fatal(FARGS, "Need at least 4 atoms in index group to use calc_ang.");
    }
    
    // Initializing variables.
    int i;
    real sumcosxy,cosxy;
    rvec vec1,vec2;
    real cosx,cosx2,cosy,cosy2;
    
    sumcosxy = 0.0;
    
    // Loop from the third atom to the last atom. Forward direction.
    for (i = 1; i < (iatoms - 1); i++)
    {
        /* Find the value of cosx and cosy.
         * 
         * x is the angle of two vectors made by three atom coords
         * in frame[]. y is the same for rframe[].
         */
        rvec_sub(frame[i-1],  frame[i], vec1);
        rvec_sub(frame[i+1],  frame[i], vec2);
        cosx = iprod(vec1, vec2) / (norm(vec1)*norm(vec2));
        
        rvec_sub(rframe[i-1], rframe[i], vec1);
        rvec_sub(rframe[i+1], rframe[i], vec2);
        cosy = iprod(vec1, vec2) / (norm(vec1)*norm(vec2));
        
        /* Find the value of cos(x-y) using trig identities.
         * 
         * cos(x-y) = cos(x)*cos(y) + sin(x)*sin(y)
         * sin(x) = sin(acos(cos(x)))
         * sin(acos(z)) = sqrt(1 - z^2)
         *
         * Therefore:
         * cos(x-y) = 
         * cos(x)*cos(y) + sqrt(1-(cos(x))^2)*sqrt(1-(cos(y))^2) 
         */
        // First check for boundary conditions.
        cosx2 = cosx*cosx;
        if (cosx2 >  1)
        {
            cosx2 =  1.0;
        }
        
        cosy2 = cosy*cosy;
        if (cosy2 >  1)
        {
            cosy2 =  1.0;
        }
        
        // This results in a number between -1 and +1.
        cosxy = cosx*cosy + sqrt(1.0-cosx2)*sqrt(1.0-cosy2);
        // Normalize to between 0 and 1.
        cosxy = (cosxy - 1.0)/(-2.0);
        // For average sum.
        sumcosxy += cosxy;
    }
    
    // Normalize for number of angles summed.
    sumcosxy /= (iatoms-2);
    // Boundary check.
    if (sumcosxy > 1.0)
    {
        return 1.0;
    }
    if (sumcosxy < 0.0)
    {
        return 0.0;
    }
    // Finished.
    return sumcosxy;
}



real calc_dih(int iatoms,rvec frame[],rvec rframe[],real dih[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS,"Need at least 6 atoms in index group to use calc_dih.");
    }
    
    // Initializing variables.
    int i;
    real cosxy,sumcosxy;
    rvec vec1,vec2,vec3,pvec1,pvec2;
    real iang,rang;
    
    sumcosxy = 0.0;
    
    // Loop from the third atom to the last atom. Forward direction.
    for (i=3; i < iatoms; i++)
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
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i-3],frame[i-2],vec1);
        rvec_sub(frame[i-1],frame[i-2],vec2);
        rvec_sub(frame[i-1],frame[i],  vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1,vec2,pvec1);
        cprod(vec2,vec3,pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1,pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1,pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i-3],rframe[i-2],vec1);
        rvec_sub(rframe[i-1],rframe[i-2],vec2);
        rvec_sub(rframe[i-1],rframe[i],  vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1,vec2,pvec1);
        cprod(vec2,vec3,pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1,pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1,pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for the dot product and normalize.
        cosxy = (real)(cos(((double)rang) - ((double)iang)) - 1.0)/(-2.0);
        // Overall summation.
        sumcosxy += cosxy;
    }
    
    // Normalize for number of angles summed.
    sumcosxy /= (iatoms-3);
    // Boundary check.
    if (sumcosxy > 1)
    {
        return 1.0;
    }
    if (sumcosxy < 0)
    {
        return 0.0;
    }
    
    // Finished.
    return sumcosxy;
}



real calc_angdih(int iatoms, rvec frame[], rvec rframe[], real angdih[])
{
    // Error checking.
    if (iatoms < 6)
    {
        // Not sure if this works. Depends on how FARGS works.
        gmx_fatal(FARGS, "Need at least 6 atoms in index group to use calc_angdih.");
    }
    // Initializing variables.
    int i;
    real angs, dihs, cosxy;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real cosx, cosx2, cosy, cosy2, iang, rang;
    
    angs = 0.0; dihs = 0.0;
    
    // Loop from the third atom to the last atom. Calculates backbone angles.
    for (i = 1; i < (iatoms - 1); i++)
    {
        /* Find the value of cosx and cosy.
         * 
         * x is the angle of two vectors made by three atom coords
         * in frame[]. y is the same for rframe[].
         */
        rvec_sub(frame[i-1],  frame[i], vec1);
        rvec_sub(frame[i+1],  frame[i], vec2);
        cosx = iprod(vec1, vec2) / (norm(vec1)*norm(vec2));
        
        rvec_sub(rframe[i-1], rframe[i], vec1);
        rvec_sub(rframe[i+1], rframe[i], vec2);
        cosy = iprod(vec1, vec2) / (norm(vec1)*norm(vec2));
        
        /* Find the value of cos(x-y) using trig identities.
         * 
         * cos(x-y) = cos(x)*cos(y) + sin(x)*sin(y)
         * sin(x) = sin(acos(cos(x)))
         * sin(acos(z)) = sqrt(1 - z^2)
         *
         * Therefore:
         * cos(x-y) = 
         * cos(x)*cos(y) + sqrt(1-(cos(x))^2)*sqrt(1-(cos(y))^2) 
         */
        // First check for boundary conditions.
        cosx2 = cosx*cosx;
        if (cosx2 >  1)
        {
            cosx2 =  1.0;
        }
        
        cosy2 = cosy*cosy;
        if (cosy2 >  1)
        {
            cosy2 =  1.0;
        }
        
        // This results in a number between -1 and +1.
        cosxy = cosx*cosy + sqrt(1.0-cosx2)*sqrt(1.0-cosy2);
        // Normalize to between 0 and 1.
        cosxy = (cosxy - 1.0)/(-2.0);
        // For average sum.
        angs += cosxy;
    }
    
    // Loop from the third atom to the last atom. Calculate backbone dihedrals.
    for (i=3; i < iatoms; i++)
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
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[i-3],frame[i-2],vec1);
        rvec_sub(frame[i-1],frame[i-2],vec2);
        rvec_sub(frame[i-1],frame[i],  vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1,vec2,pvec1);
        cprod(vec2,vec3,pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1,pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1,pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[i-3],rframe[i-2],vec1);
        rvec_sub(rframe[i-1],rframe[i-2],vec2);
        rvec_sub(rframe[i-1],rframe[i],  vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1,vec2,pvec1);
        cprod(vec2,vec3,pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1,pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1,pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for the dot product and normalize.
        cosxy = (real)(cos(((double)rang) - ((double)iang)) - 1.0)/(-2.0);
        // Overall summation.
        dihs += cosxy;
    }
    
    // Divide by number of angles summed.
    angs /= (iatoms - 2);
    dihs /= (iatoms - 3);
    
    // Boundary check.
    if (angs > 1)
    {
        angs = 1.0;
    }
    if (angs < 0)
    {
        angs = 0.0;
    }
    if (dihs > 1)
    {
        dihs = 1.0;
    }
    if (dihs < 0)
    {
        dihs = 0.0;
    }
    
    // Finished. Return geometric mean of the two.
    return sqrt(angs * dihs);
}



real calc_phipsi(int iatoms, rvec frame[], rvec rframe[], real phipsi[])
{
    // Error checking.
    if (iatoms < 6)
    {
        gmx_fatal(FARGS,"Need at least 6 atoms in index group to use calc_phipsi.");
    }
    
    // Initializing variables.
    int i, iC1 ,iC2, iN1, iN2, iCa;
    real cosxy, phi, psi;
    rvec vec1, vec2, vec3, pvec1, pvec2;
    real iang, rang;
    
    phi = 0.0; psi = 0.0;
    
    // Phi angles.
    for (i=1; i < iatoms/3; i++)
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
        iC1 = 3*i - 1;
        iN1 = 3*i;
        iCa = 3*i + 1;
        iC2 = 3*i + 2;
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[iC1],frame[iN1],vec1);
        rvec_sub(frame[iCa],frame[iN1],vec2);
        rvec_sub(frame[iCa],frame[iC2],vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1,vec2,pvec1);
        cprod(vec2,vec3,pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1,pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1,pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[iC1],rframe[iN1],vec1);
        rvec_sub(rframe[iCa],rframe[iN1],vec2);
        rvec_sub(rframe[iCa],rframe[iC2],vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1, vec2, pvec1);
        cprod(vec2, vec3, pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1,pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1,pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for the dot product and normalize.
        cosxy = (real)(cos(((double)iang) - ((double)rang)) - 1.0)/(-2.0);
        // Final sum.
        phi += cosxy;
    }
    
    
    
    // Psi angles.
    for (i=1; i < iatoms/3; i++)
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
        iN1 = 3*i - 3;
        iCa = 3*i - 2;
        iC1 = 3*i - 1;
        iN2 = 3*i;
        // Convert four atom coordinates into three vectors.
        rvec_sub(frame[iN1],frame[iCa],vec1);
        rvec_sub(frame[iC1],frame[iCa],vec2);
        rvec_sub(frame[iC1],frame[iN2],vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1,vec2,pvec1);
        cprod(vec2,vec3,pvec2);
        // Find the angle between plane vectors.
        iang = gmx_angle(pvec1,pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1,pvec2) < 0.0)
        {
            iang *= -1.0;
        }
        
        /* Repeat for the reference frame.
         */
        // Make three vectors.
        rvec_sub(rframe[iN1],rframe[iCa],vec1);
        rvec_sub(rframe[iC1],rframe[iCa],vec2);
        rvec_sub(rframe[iC1],rframe[iN2],vec3);
        // Convert three vectors into two normal to each plane.
        cprod(vec1,vec2,pvec1);
        cprod(vec2,vec3,pvec2);
        // Find the angle between plane vectors.
        rang = gmx_angle(pvec1,pvec2);
        // Calculate and apply the sign.
        if (iprod(vec1,pvec2) < 0.0)
        {
            rang *= -1.0;
        }
        // Solve for the dot product and normalize range: 0 to 1.
        cosxy = (real)(cos(((double)iang) - ((double)rang)) - 1.0)/(-2.0);
        // Final sum.
        psi += cosxy;
    }
    
    // Divide by the number of angles summed.
    phi /= (iatoms / 3) - 1;
    psi /= (iatoms / 3) - 1;
    
    // Boundary check.
    if (phi > 1)
    {
        phi = 1.0;
    }
    if (phi < 0)
    {
        phi = 0.0;
    }
    if (psi > 1)
    {
        psi = 1.0;
    }
    if (psi < 0)
    {
        psi = 0.0;
    }
    
    // Finished. Return geometric mean of phi and psi.
    return sqrt(phi * psi);
}



real calc_msd(int iatoms,rvec x[],rvec xref[],real msd[])
{
    /* This is a simplified version of something in the gromacs lib.
     * 
     * While I was at it, I changed the return from rmsd to msd.
     * 
     * Inputs:
     * Compare the coordinates in x to the reference coordinates in xref.
     * Both should contain iatoms of coordinates in rvec format.
     * 
     * Outputs:
     * Array msd stores the msd per atom.
     */
    int i,d;
    real xd,msdt,msdi;
    
    msdt=0;
    for(i=0;i<iatoms;i++)
    {
        // Calculations to find msd.
        msdi=0;
        for(d=0;d<3;d++)
        {
            xd = x[i][d] - xref[i][d];
            msdi += xd * xd;
        }
        msd[i] = msdi;
        msdt += msdi;
    }
    // Normalize for number of atoms.
    msdt /= iatoms;
    
    return msdt;
}




real calc_scaled_msd(int iatoms,rvec x[],rvec xref[],real msd[])
{
    /*
     * The mirror_msd is the average msd between the reference and its
     * mirrored structures mirrored around 3 planes. Used to normalize
     * the msd value.
     */
    int i,j,d;
    real xd,xd1,xd2,msdt,msdi,mirror_msd;
    
    // Initialize to zero.
    msdt = 0;
    mirror_msd = 0;
    
    // Main loop.
    for(i=0;i<iatoms;i++)
    {
        // Initialize ith sum to zero.
        msdi = 0;
        // Copy over coordinates from xref and rotate them.
        for(d=0;d<3;d++)
        {
            /* Same as multiplication by (-1,-1,-1) diag matrix.
             * 
             * d = 0: x coordinate, etc
             * 
             * xd1 = distances from fit coord to mirrored ref coord.
             * xd2 = distances from ref coord to mirrored fit coord.
             * 
             * Originally this did something slightly different. Now it's
             * pointless to separate xd1 and xd2.
             * 
             * Changing this. Since xd1^2 = xd2^2, most of the following
             * is unnecessary. I'll leave the lines of code in case we
             * want to change the scaling to something else.
             */
            xd1 = x[i][d] + xref[i][d];
            //xd2 = -1 * xd1;
            // Add to the sum.
            mirror_msd += xd1 * xd1; // + xd2 * xd2
            
            // Calculate ith msd.
            xd = x[i][d] - xref[i][d];
            msdi += xd * xd;
        }
        msdt  += msdi;
    }
    
    // Normalize for number of atoms.
    msdt /= iatoms;
    mirror_msd /= iatoms;
    /* Divide by two since both ref and fit structure are used as a mirrored
     * reference. (No longer necessary due to changes above.)
     */
    //mirror_msd /= 2;
    
    // Scale by the mirrored reference msd.
    msdt /= mirror_msd;
    
    return msdt;
}



real calc_mirror_msd(int iatoms,rvec frame[],rvec rframe[],real mirror_msd[])
{
    /*
     * MSD between fit structure and mirror of reference structure.
     * 
     * Assumes they have already been aligned.
     */
    int i, d;
    real xd, xdm, mirror_msdt, mirror_msdi;
    
    // Initialize to zero.
    mirror_msdt = 0;
    
    // Main loop.
    for(i=0;i<iatoms;i++)
    {
        // Initialize ith sum to zero.
        mirror_msdi = 0;
        // Copy over coordinates from xref and rotate them.
        for(d=0;d<3;d++)
        {
            /* Same as multiplication by (-1,-1,-1) diag matrix.
             * 
             * d = 0: x coordinate, etc
             * 
             * xdm = distances from fit coord to mirrored ref coord.
             */
            xdm = frame[i][d] + rframe[i][d];
            // Add to the sum.
            mirror_msdi += xdm * xdm;
        }
        mirror_msdt  += mirror_msdi;
    }
    
    // Normalize for number of atoms.
    mirror_msdt /= iatoms;
    
    return mirror_msdt;
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

real call_ISDM(int iatoms, rvec cframe[], rvec rframe[], real diff[], 
                 const char *ISDM)
{
    // Default behavior (no -ISDM option) is RMSD.
    if (strcmp(ISDM, "RMSD") == 0)
    {
        // Calculate RMSD.
        return sqrt(calc_msd(iatoms,cframe,rframe,diff));
    }
    
    // Scaled RMSD. User gives -srms option.
    if (strcmp(ISDM, "SRMS") == 0)
    {
        // Calculate scaled RMSD.
        return sqrt(calc_scaled_msd(iatoms,cframe,rframe,diff));
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
    }
    
    // Difference of end-to-end distance. User gives -e2e option.
    if (strcmp(ISDM, "E2E") == 0)
    {
        real ISD;
        // Find difference of end-to-end distances.
        ISD  = (real)sqrt(distance2(cframe[(iatoms - 1)],cframe[0]));
        ISD -= (real)sqrt(distance2(rframe[(iatoms - 1)],rframe[0]));
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
        ISD  = (real)sqrt(distance2(cframe[(iatoms - 1)],cframe[0]));
        ISD -= (real)sqrt(distance2(rframe[(iatoms - 1)],rframe[0]));
        // Find scaled absolute value.
        return sqrt(ISD * ISD) / (rg + rgr);
    }
    
    // Mirrored RMSD. User gives -mir option.
    if (strcmp(ISDM, "MIR") == 0)
    {
        // Calculate mirrored RMSD.
        return (real)sqrt(calc_mirror_msd(iatoms,cframe,rframe,diff));
    }
    
    // Angles. User gives -ang option.
    if (strcmp(ISDM, "ANG") == 0)
    {
        // Calculate dot product of differences in angles.
        return calc_ang(iatoms,cframe,rframe,diff);
    }
    
    // Dihedrals. User gives -dih option.
    if (strcmp(ISDM, "DIH") == 0)
    {
        // Calculate dot product of difference in dihedrals.
        return calc_dih(iatoms,cframe,rframe,diff);
    }
    
    // Dihedrals. User gives -dih option.
    if (strcmp(ISDM, "ANGDIH") == 0)
    {
        // Calculate dot product of difference in dihedrals.
        return calc_angdih(iatoms,cframe,rframe,diff);
    }
    
    // Phi psi angles. User gives -phipsi option.
    if (strcmp(ISDM, "PHIPSI") == 0)
    {
        // Calculate dot product of difference in dihedrals.
        return calc_phipsi(iatoms,cframe,rframe,diff);
    }
    
    // Atom to atom distances. User gives -drms option.
    if (strcmp(ISDM, "DRMS") == 0)
    {
        /* Calculate differences in all atom to atom distances.
         * 
         * Scaled by 2 * sqrt(Rgi * Rgj).
         * 
         * A few cycles could be saved by storing Rgi, but most
         * necessary info to calculate Rgi needs to be calculated
         * for every comparison anyway.
         */
        return calc_drms(iatoms,cframe,rframe,diff);
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
}
