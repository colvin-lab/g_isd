/* Modified Aug 07 2013
 *
 * Robert L. Wang - rwang6@ucmerced.edu
 * 
 * This code was originally a MATLAB code that was translated
 * into C code in order to work with the GROMACS libraries to perform
 * elastic shape analysis as described in: 
 * 
 * Liu W, Srivastava A, Zhang J (2011) A Mathematical Framework for 
 * Protein Structure Comparison. PLoS Comput Biol 7(2): e1001075.
 * doi:10.1371/journal.pcbi.1001075
 * 
 * Original MATLAB code was provided by the authors upon request.
 * All functions aside from the svd were rewritten in C to perform their
 * intended duties.
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
//#include "xdrfile_xtc.h"
//#include "compframes.h"

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

//Used in the svd function
#define NNBRS   23

const int Nbrs[NNBRS][2] = { 
        { 1, 1 }, 
        { 1, 2 },
        { 2, 1 },
        { 2, 3 },
        { 3, 2 },
        { 1, 3 },
        { 3, 1 },
        { 1, 4 },
        { 3, 4 },
        { 4, 3 },
        { 4, 1 },
        { 1, 5 },
        { 2, 5 },
        { 3, 5 },
        { 4, 5 },
        { 5, 4 },
        { 5, 3 },
        { 5, 2 },
        { 5, 1 },
        { 1, 6 },
        { 5, 6 },
        { 6, 5 },
        { 6, 1 }
};



double esa_analysis(int iatoms, rvec frame[], rvec rframe[], int sample);
void resample_curve(rvec*,rvec*, int, int);
void curve_to_q(rvec*, double**, int);
void q_rot(double**, double**, int, double**);
void invert_gamma(double*, int, double*);
void Group_Action_by_Gamma_Coord(rvec*, double*, rvec*, int);
void transpmatmn(double**, int, int, double**);
void matmultmno(double**, double**, int, int, int, double**);
double det3(double** a);
void DynamicProgrammingQ(double* q1, double* q2, double* yy, int m, int z);
double CostFn2(double *q1, double *q2, double *q2L, int k, int l, int i, int j, int n, int N, int M, double lam);
void thomas(double *x, double *a, double *b, double *c, int n);
void spline(double *D, double *y, int n);
void lookupspline(double *t, int *k, double dist, double len, int n);
double evalspline(double t, double D[2], double y[2]);
static double PYTHAG(double a, double b);
int dsvd(double **a, int m, int n, double *w, double **v);

/*
 * Note: This esa_analysis function was originally a piece of MATLAB code called
 * mygeod.m. It is the main structure of the program that calls all of the other functions
 * in order to perform the elastic shape analysis. It has been ported to C and rewritten to
 * fit with all of the other functions.
 */

/******************************************************************************************************************/
double esa_analysis(int iatoms, rvec frame[], rvec rframe[], int sample)
{
    // Form the q function for representing curves and find best rotation.
    double **q1, **q2, **q2n, **Ot, **Ot_t, **q2r;
    double *q1vec, *q2vec, *gam, *gam_I;
    double inner_prod_q1, inner_prod_q2, step, sum, dist, tmp, tmp1, tmp2;
    rvec   *x1, *x2, *x2n;
    int    i,j;
    
    // Memory allocation.
    snew(x1,    sample);
    snew(x2,    sample);
    //snew(x2_n,  sample);
    snew(gam,   sample);
    snew(gam_I, sample);
    snew(q1vec, 3*sample);
    snew(q2vec, 3*sample);
    
    snew(q1,  sample);
    snew(q2,  sample);
    snew(q2r, sample);
    for (i = 0; i < sample; i++)
    {
        snew(q1[i],  3);
        snew(q2[i],  3);
        snew(q2r[i], 3);
    }
    snew(Ot,   3);
    snew(Ot_t, 3);
    for (i = 0; i < 3; i++)
    {
        snew(Ot[i],   3);
        snew(Ot_t[i], 3);
    }
    
    // Resample frames to curves.
    resample_curve(frame,  x1, sample, iatoms);
    resample_curve(rframe, x2, sample, iatoms);
    
    // Curve (x1 and x2) to q.
    curve_to_q(x1, q1, sample);
    curve_to_q(x2, q2, sample);
    
    // Determine orthogonal rotation matrix Ot and transpose to Ot_t.
    q_rot(q1, q2, sample, Ot);
    transpmatmn(Ot, 3, 3, Ot_t);
    
    //Multiply q2 by rotation matrix Ot
    matmultmno(q2, Ot_t, sample, 3, 3, q2r);
    for (i = 0; i < sample; i++)
    {
        for (j = 0; j < 3; j++)
        {
            q2[i][j] = q2r[i][j];
        }
    }
    
    // Determine inner product by integrating q over 0 -> 1 using the trapezoidal rule.
    inner_prod_q1 = 0;
    inner_prod_q2 = 0;
    step = 1.0 / (double)(sample-1);
    tmp1 = q1[0][0]*q1[0][0] + q1[0][1]*q1[0][1] + q1[0][2]*q1[0][2];
    tmp2 = q2[0][0]*q2[0][0] + q2[0][1]*q2[0][1] + q2[0][2]*q2[0][2];
    for (i = 1; i < sample; i++)
    {
        tmp  = tmp1;
        tmp1 = q1[i][0]*q1[i][0] + q1[i][1]*q1[i][1] + q1[i][2]*q1[i][2];
        tmp += tmp1;
        inner_prod_q1 += (tmp * step / 2);
        
        tmp  = tmp2;
        tmp2 = q2[i][0]*q2[i][0] + q2[i][1]*q2[i][1] + q2[i][2]*q2[i][2];
        tmp += tmp2;
        inner_prod_q2 += (tmp * step / 2);
    }
    inner_prod_q1 = sqrt(inner_prod_q1);
    inner_prod_q2 = sqrt(inner_prod_q2);
    
    // Reformat q to vectors and run DynamicProgrammingQ.
    for (i = 0; i < sample; i++)
    {
        for (j = 0; j < 3; j++)
        {
            q1vec[3*i+j] = q1[i][j]/inner_prod_q1;
            q2vec[3*i+j] = q2[i][j]/inner_prod_q2;
        }
    }
    DynamicProgrammingQ(q1vec, q2vec, gam, 3, sample);
    
    // x1 and q2 are no longer used. Assign the memory for reuse.
    x2n = x1;
    q2n = q2;
    
    // Invert gam. Convert x2 to x2n using function. Convert curve to q2n.
    invert_gamma(gam, sample, gam_I);
    Group_Action_by_Gamma_Coord(x2, gam_I, x2n, sample);
    curve_to_q(x2n,q2n,sample); 
    
    // Determining orthogonal rotation matrix Ot and transpose to Ot_t.
    q_rot(q1, q2n, sample, Ot);
    transpmatmn(Ot, 3, 3, Ot_t);
    
    //Multiply q2n by rotation matrix Ot.
    matmultmno(q2n, Ot_t, sample, 3, 3, q2r);
    for (i = 0; i < sample; i++)
    {
        for (j = 0; j < 3; j++)
        {
            q2n[i][j] = q2r[i][j];
        }
    }
    
    // Determine distance between q1 and q2.
    sum = 0;
    for (i = 0; i < sample; i++)
    {
        tmp = 0;
        for (j = 0; j < 3; j++)
        {
            tmp += q1[i][j] * q2[i][j];
        }
        sum += tmp / (double)sample;
    }
    dist = acos(sum);
    
    // Free memory.
    sfree(x1);
    sfree(x2);
    sfree(gam);
    sfree(gam_I);
    sfree(q1vec);
    sfree(q2vec);
    
    for (i = 0; i < sample; i++)
    {
        sfree(q1[i]);
        sfree(q2[i]);
        sfree(q2r[i]);
    }
    sfree(q1);
    sfree(q2);
    sfree(q2r);
    for (i = 0; i < 3; i++)
    {
        sfree(Ot[i]);
        sfree(Ot_t[i]);
    }
    sfree(Ot);
    sfree(Ot_t);
    
    // Output.
    return dist;
}

/*
 * Note: This resample_curve function was originally a piece of
 * MATLAB code called ReSampleCurve.m that was ported into C and modified
 * to work with the rest of this program.
 */

/******************************************************************************************************************/
void resample_curve(rvec* frame, rvec* resampled, int sample, int natoms)
{
    int i,j,k;
    double del,newdel,sum;
    double arr[3];
    double *cumdel;
    
    snew(cumdel,natoms);
    
    cumdel[0] = 0;
    sum = 0.0;
    
    /* At some point, we should make this a command line option.
     * The original code does this, but it doesn't seem like a good idea.
     * 
     * // Before doing anything, preprocess the frames to smooth the curve.
     * for(i = 0; i < 2; i++)
     * {
     *     // Do this twice.
     *     for(j = 1; j < natoms-1; j++)
     *     {
     *         for(k = 0; k < 3; k++)
     *         {
     *             frame[j][k] = (frame[j-1][k] + 2*frame[j][k] + frame[j+1][k])/4.0;
     *         }
     *     }
     * }
     */
    
    //Note: The following segment of code calculates the L2 norm of the
    //coordinates. This was originally a MATLAB function called norm which has
    //been rewritten in C.

    //Calculate L2 norm     
    for (i = 1; i < natoms; i++)
    {
        arr[0] = frame[i][0] - frame[i-1][0];
        arr[1] = frame[i][1] - frame[i-1][1];
        arr[2] = frame[i][2] - frame[i-1][2];
        del = sqrt(arr[0]*arr[0] + arr[1]*arr[1] + arr[2]*arr[2]);
        sum += del;
        cumdel[i] = del + cumdel[i-1];
    }
    for (i = 0; i < natoms; i++)
    {
        cumdel[i] /= sum;
    }
    
    //First and last points of the resampled curve are the same as the first and last points of the frame
    for (k = 0; k < 3; k++)
    {
        resampled[0][k] = frame[0][k];
        resampled[sample-1][k] = frame[natoms-1][k];
    }

    //Note: This interpolating step is a MATLAB function called interp1 that has been
    //rewritten to do the same thing in C
    
    //Interpolate
    //y = y_0 + (y_1 - y_0)*[(x - x_0)/(x_1-x_0)]
    j = 0;
    for (i = 1; i < (sample-1); i++)
    {
        //The new points should be evenly spaced from point 0 to point 199
        newdel = (double)i / (sample - 1);
        while ((newdel > cumdel[j+1]) && (j < natoms))
        {
            j++;
        }
        
        // Solve interpolation.
        for(k = 0; k < 3; k++)
        {
            resampled[i][k] = frame[j][k] + (frame[j+1][k] - frame[j][k]) * ((newdel - cumdel[j]) / (cumdel[j+1]-cumdel[j]));
        }
    }
    
    // Free memory.
    sfree(cumdel);
}

/*
 * Note: This curve_to_q function was originally a piece of
 * MATLAB code called curve_to_q.m that was ported into C and modified
 * to work with the rest of this program.
 */

/******************************************************************************************************************/
void curve_to_q(rvec* curve, double** q, int sample)
{
    double step, len, tmp;
    double **grad;
    int i;
    
    snew(grad, sample);
    for (i = 0; i < sample; i++)
    {
        snew(grad[i], 3);
    }
    
    len = 0;
    step = 1.0 / (double)sample;

    //Note: The following calculates the gradient of the curve.
    //This was originally a MATLAB function called gradient and was
    //rewritten in C to perform the same task.
    
    //Compute the gradient of the curve
    grad[0][0] = (curve[1][0] - curve[0][0])/step;
    grad[0][1] = (curve[1][1] - curve[0][1])/step;
    grad[0][2] = (curve[1][2] - curve[0][2])/step;
    
    for (i = 1; i < sample-1; i++)
    {
        grad[i][0] = (double)(curve[i+1][0] - curve[i-1][0])/(2*step);
        grad[i][1] = (double)(curve[i+1][1] - curve[i-1][1])/(2*step);
        grad[i][2] = (double)(curve[i+1][2] - curve[i-1][2])/(2*step);
    }
    
    grad[sample-1][0] = (curve[sample-1][0] - curve[sample-2][0])/step;
    grad[sample-1][1] = (curve[sample-1][1] - curve[sample-2][1])/step;
    grad[sample-1][2] = (curve[sample-1][2] - curve[sample-2][2])/step;
    
    //Determine len
    for (i = 0; i < sample; i++)
    {
        tmp = sqrt(grad[i][0]*grad[i][0] + grad[i][1]*grad[i][1] + grad[i][2]*grad[i][2]);
        len += tmp;
    }
    len /= sample;
    
    for (i = 0; i < sample; i++) 
    {
        grad[i][0] /= len;
        grad[i][1] /= len;
        grad[i][2] /= len;
    }
    
    // Solve for q.
    for (i = 0; i < sample; i++)
    {
        tmp = sqrt(grad[i][0]*grad[i][0] + grad[i][1]*grad[i][1] + grad[i][2]*grad[i][2]);
        tmp = sqrt(tmp);
        
        if (tmp > 0.0001)
        {
            q[i][0] = grad[i][0]/tmp;
            q[i][1] = grad[i][1]/tmp;
            q[i][2] = grad[i][2]/tmp;
        }
        else
        {
            // Not the same as the original code. Originally this was a multiply.
            // We think it was probably a bug.
            q[i][0] = grad[i][0] / 0.0001;
            q[i][1] = grad[i][1] / 0.0001;
            q[i][2] = grad[i][2] / 0.0001;
        }
    }
    
    // Free memory.
    for (i = 0; i < sample; i++)
    {
        sfree(grad[i]);
    }
    sfree(grad);
}

/*
 * Note: This q_rot function was originally a segment of
 * MATLAB code in the mygeod.m function that performed the task of
 * finding the optimal rotation of the curve. It was extracted as a separate
 * function since it was called several times, and was modified
 * to work with the rest of this program.
 */

/******************************************************************************************************************/
void q_rot(double** q1, double** q2, int sample, double** Ot)
{
    double **q1_t, **A, **V, **V_t;
    double det_A;
    double S[3];
    int i;
    
    // Memory allocation.
    snew(q1_t, 3);
    snew(A, 3);
    snew(V, 3);
    snew(V_t, 3);
    for(i = 0; i < 3; i++)
    {
        snew(q1_t[i], sample);
        snew(A[i], 3);
        snew(V[i], 3);
        snew(V_t[i],3);
    }
    
    // Transpose q1 and store it in q1_t.
    transpmatmn(q1, sample, 3, q1_t);
    
    // Multiply q1_t * q2, store to A, then take the determinant of A.
    matmultmno(q1_t, q2, 3, sample, 3, A);
    det_A = det3(A);
    
    //Perform svd, A becomes U
    dsvd(A, 3, 3, S, V);
    
    // Solve for Ot.
    if (det_A > 0)
    {
        transpmatmn(V, 3, 3, V_t);
        matmultmno(A, V_t, 3, 3, 3, Ot);
    }
    else
    {
        V[0][2] = -V[0][2];
        V[1][2] = -V[1][2];
        V[2][2] = -V[2][2];
        transpmatmn(V, 3, 3, V_t);
        matmultmno(A, V_t, 3, 3, 3, Ot);
    }
    
    // Free memory.
    for (i = 0; i < 3; i++)
    {
        sfree(q1_t[i]);
        sfree(A[i]);
        sfree(V[i]);
        sfree(V_t[i]);
    }
    sfree(q1_t);
    sfree(A);
    sfree(V);
    sfree(V_t);
}

/*
 * Note: This invert_gamma function was originally a piece of
 * MATLAB code called invertgamma.m that was ported into C and modified
 * to work with the rest of this program.
 */

/******************************************************************************************************************/
void invert_gamma(double* gam, int sample, double* gam_I)
{
    double xi, xj1, xj2;
    double n1, n2;
    int i,j;
    
    //for(i = 1; i <= sample; i++) {
    //    x[(i-1)] = i / (double)sample;
    //}
    
    gam_I[0] = 1.0 / (double)sample;
    gam_I[sample-1] = 1.0;
    
    //Note: This interpolating step is a MATLAB function called interp1 that has been
    //rewritten to do the same thing in C

    //Interpolate
    //y = y_0 + (y_1 - y_0)*[(x - x_0)/(x_1-x_0)]
    j = 0;
    for(i = 1; i < (sample-1); i++)
    {
        // Solve for x, y_0, y_1.
        xi = (i + 1) / (double)sample;
        while ((xi > gam[j+1]) && (j < sample))
        {
            j++;
        }
        xj1 = (j + 1) / (double)sample;
        xj2 = (j + 2) / (double)sample;
        
        // Solve interpolation.
        gam_I[i] = xj1 + (xj2 - xj1) * ((xi - gam[j]) / (gam[j+1]-gam[j]));
    }
    
    // Normalize gam_I to between 0 and 1.
    n1 = gam_I[0];
    n2 = gam_I[sample-1];
    for(i = 0; i < sample; i++) {
        gam_I[i] = (gam_I[i]-n1)/(n2-n1);
    }
}

/*
 * Note: This Group_Action_by_Gamma_Coord function was originally a
 * piece of MATLAB code called Group_Action_by_Gamma_Coord.m
 * that has been ported into C code and modified to work with the
 * rest of this program
 */

/******************************************************************************************************************/
void Group_Action_by_Gamma_Coord(rvec* x2, double* gam_I, rvec* x2_n, int sample)
{
    double x_0, x_1;
    int i,j,k;
    
    //First and last elements of x2 and x2_n are the same
    for (i = 0; i < 3; i++) {
        x2_n[0][i] = x2[0][i];
        x2_n[sample-1][i] = x2[sample-1][i];
    }

    //Note: This interpolating step is a MATLAB function called interp1 that has been
    //rewritten to do the same thing in C
    
    //Interpolate   
    //      //y = y_0 + (y_1 - y_0)*[(x - x_0)/(x_1-x_0)]
    j = 0;
    x_1 = 1.0 / (double)(sample-1);
    for (i = 1; i < (sample-1); i++)
    {
        while ((gam_I[i] > x_1) && (j < sample))
        {
            j++;
            x_1 = (j+1) / (double)(sample-1);
        }
        x_0 = j / (double)(sample-1);
        for (k = 0; k < 3; k++)
        {
            x2_n[i][k] = x2[j][k] + (x2[j+1][k] - x2[j][k]) * ((gam_I[i] - x_0) / (x_1 - x_0));
        }
    }
}

/*
 * Note: The following functions det3, matmultmno, transpmatmn
 * are basic functions for matrix operations that are called several
 * times throughout the program.
 */

/******************************************************************************************************************/
double det3(double** a)
{
    /* Assume a is a 3 by 3 matrix
     * Will determine the determinant of a
     */
    double sum = 0;
    sum = a[0][0] * (a[1][1]*a[2][2] - a[1][2]*a[2][1])
    -a[0][1] * (a[1][0]*a[2][2] - a[1][2]*a[2][0])
    +a[0][2] * (a[1][0]*a[2][1] - a[1][1]*a[2][0]);
    
    return sum;
}
/******************************************************************************************************************/
void matmultmno(double** a, double** b, int m, int n, int o, double** out)
{
    /* Assume a is an array of reals m by n in dimensions.
     * Assume b is an array of reals n by o in dimensions.
     * Out should point to enough memory to store m by o reals.
     */
    int i, j, k;
    for (i = 0; i < m; i++) {
        for (k = 0; k < o; k++) {
            out[i][k] = 0;
            for (j = 0; j < n; j++) {
                out[i][k] += a[i][j] * b[j][k];
            }
        }
    }
}
/******************************************************************************************************************/
void transpmatmn(double** a, int m,int n, double** out)
{
    /* Transpose m by n array of reals a into n by m array of reals out.
     */
    int i, j;
    
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            out[j][i] = a[i][j];
        }
    }
}

/*
 * Note: This DynamicProgrammingQ function as well as the following several functions
 * were originally a separate piece of C code called dynamicprogrammingq.c. 
 * It was modified to fit with this program and moved into the same source file
 */

void DynamicProgrammingQ(double* q1, double* q2, double* yy, int m, int z)
{
    //m x z matrix
    int i, j, k, l, n, M, N, Eidx, Fidx, Ftmp, Fmin, Num, *Path, *x, *y, cnt;
    double *q2L, *D, *tmp, *E, Etmp, Emin, t, a, b, lam = 0;
    
    n = m;
    N = z;
    
    M = 5*N;
    q2L = malloc(n*M*sizeof(double));
    
    D = malloc(2*N*sizeof(double));
    tmp = D + N;
    
    a = 1.0/N;
    b = 1.0;
    
    for (i = 0; i < n; ++i) {
        for (j = 0; j < N; ++j) {
            tmp[j] = q2[n*j + i];
        }
        
        spline(D, tmp, N);
        
        for (j = 0; j < M; ++j) {
            /* XXX: Extrapolation 1/M < 1/N */
            lookupspline(&t, &k, (j+1.0)/M - a, b - a, N);
            q2L[n*j + i] = evalspline(t, D+k, tmp+k);
        }
    }
    
    free(D);
    
    E = calloc(N*N, sizeof(double));
    Path = malloc(2*N*N*sizeof(int));
    
    for (i = 0; i < N; ++i) {
        E[N*i + 0] = 1;
        E[N*0 + i] = 1;
        Path[N*(N*0 + i) + 0] = -1;
        Path[N*(N*0 + 0) + i] = -1;
        Path[N*(N*1 + i) + 0] = -1;
        Path[N*(N*1 + 0) + i] = -1;
    }
    E[N*0 + 0] = 0;
    
    for (j = 1; j < N; ++j) {
        for (i = 1; i < N; ++i) {
            
            Emin = 100000;
            Eidx = 0;
            
            for (Num = 0; Num < NNBRS; ++Num) {
                k = i - Nbrs[Num][0];
                l = j - Nbrs[Num][1];
                
                if (k >= 0 && l >= 0) {
                    Etmp = E[N*l + k] + CostFn2(q1,q2,q2L,k,l,i,j,n,N,M,lam);
                    if (Num == 0 || Etmp < Emin) {
                        Emin = Etmp;
                        Eidx = Num;
                    }
                }
            }
            
            E[N*j + i] = Emin;
            Path[N*(N*0 + j) + i] = i - Nbrs[Eidx][0];
            Path[N*(N*1 + j) + i] = j - Nbrs[Eidx][1];
        }
    }
    
    free(E);
    free(q2L);
    
    /* XXX: x, y assumed to be at most length N */
    x = malloc(2*N*sizeof(int));
    y = x + N;
    
    x[0] = N-1;
    y[0] = N-1;
    
    cnt = 1;
    while (x[cnt-1] > 0) {
        y[cnt] = Path[N*(N*0 + x[cnt-1]) + y[cnt-1]];
        x[cnt] = Path[N*(N*1 + x[cnt-1]) + y[cnt-1]];
        ++cnt;
    }
    
    free(Path);
    for (i = 0, j = cnt-1; i < j; ++i, --j) {
        k = x[i];
        x[i] = x[j];
        x[j] = k;
        
        k = y[i];
        y[i] = y[j];
        y[j] = k;
    }
    
    for (i = 0; i < N; ++i) {
        
        Fmin = 100000;
        Fidx = 0;
        
        for (j = 0; j < cnt; ++j) {
            Ftmp = (int)fabs(i - x[j]);
            if (j == 0 || Ftmp < Fmin) {
                Fmin = Ftmp;
                Fidx = j;
            }
        }
        
        if (x[Fidx] == i) {
            yy[i] = (y[Fidx]+1);
        }
        else {
            if (x[Fidx] > i) {
                a = x[Fidx] - i;
                b = i - x[Fidx-1];
                yy[i] = (a*(y[Fidx-1]+1) + b*(y[Fidx]+1))/(a+b);
            }
            else {
                a = i - x[Fidx];
                b = x[Fidx+1] - i;
                yy[i] = (a*(y[Fidx+1]+1) + b*(y[Fidx]+1))/(a+b);
            }
        }
        
        yy[i] /= N;
    }
    
    free(x);
}

double CostFn2(double *q1, double *q2, double *q2L, int k, int l, int i, int j, int n, int N, int M, double lam)
{
    double m = (j-l)/(double)(i-k), sqrtm = sqrt(m), E = 0, y, tmp, ip, fp;
    int x, idx, d;
    
    for (x = k; x <= i; ++x) {
        y = (x-k)*m + l + 1;
        fp = modf(y*M/N, &ip);
        idx = (int)(ip + (fp >= 0.5)) - 1;
        
        for (d = 0; d < n; ++d) {
            tmp = q1[n*x + d] - sqrtm*q2L[n*idx + d];
            E += tmp*tmp;
        }
    }
    
    return E/N;
}

void thomas(double *x, double *a, double *b, double *c, int n)
{
    double tmp;
    int i;
    
    c[0] /= b[0];
    x[0] /= b[0];
    
    for (i = 1; i < n; ++i) {
        tmp = 1/(b[i] - c[i-1] * a[i]);
        c[i] *= tmp;
        x[i] = (x[i] - x[i-1] * a[i])*tmp;
    }
    
    for (i = n-2; i >= 0; --i) {
        x[i] -= c[i]*x[i+1];
    }
}

void spline(double *D, double *y, int n)
{
    int i;
    double *a, *b, *c;
    
    a = malloc(3*n*sizeof(double));
    b = a + n;
    c = b + n;
    
    if (n < 4) {
        a[0] = 0;
        b[0] = 2;
        c[0] = 1;
        D[0] = 3*(y[1]-y[0]);
        
        a[n-1] = 1;
        b[n-1] = 2;
        c[n-1] = 0;
        D[n-1] = 3*(y[n-1]-y[n-2]);
    }
    else {
        a[0] = 0;
        b[0] = 2;
        c[0] = 4;
        D[0] = -5*y[0] + 4*y[1] + y[2];
        
        a[n-1] = 4;
        b[n-1] = 2;
        c[n-1] = 0;
        D[n-1] = 5*y[n-1] - 4*y[n-2] - y[n-3];
    }
    
    for (i = 1; i < n-1; ++i) {
        a[i] = 1;
        b[i] = 4;
        c[i] = 1;
        D[i] = 3*(y[i+1]-y[i-1]);
    }
    
    thomas(D, a, b, c, n);
    
    free(a);
}

void lookupspline(double *t, int *k, double dist, double len, int n)
{
    *t = (n-1)*dist/len;
    *k = (int)floor(*t);
    
    *k = (*k > 0)*(*k);
    *k += (*k > n-2)*(n-2-*k);
    
    *t -= *k;
}

double evalspline(double t, double D[2], double y[2])
{
    double c[4];
    
    c[0] = y[0];
    c[1] = D[0];
    c[2] = 3*(y[1]-y[0])-2*D[0]-D[1];
    c[3] = 2*(y[0]-y[1])+D[0]+D[1];
    
    return t*(t*(t*c[3] + c[2]) + c[1]) + c[0];
}


/* 
 * Note: This svd code was modified from the original source provided by 
 * Diane Cook. The svd code was accessed on Apr 29 2013 from:
 * http://svn.lirec.eu/libs/magicsquares/src/SVD.cpp
 * and tailored to fit this program
 */

/************************************************************
 *                                                          *
 *  Permission is hereby granted  to  any  individual   or  *
 *  institution   for  use,  copying, or redistribution of  *
 *  this code and associated documentation,  provided       *
 *  that   such  code  and documentation are not sold  for  *
 *  profit and the  following copyright notice is retained  *
 *  in the code and documentation:                          *
 *     Copyright (c) held by Dianne Cook                    *
 *  All Rights Reserved.                                    *
 *                                                          *
 *  Questions and comments are welcome, and I request       *
 *  that you share any modifications with me.               *
 *                                                          *
 *                Dianne Cook                               *
 *             dicook@iastate.edu                           *
 *                                                          *
 ************************************************************/

/* 
 * svdcomp - SVD decomposition routine. 
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a 
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is 
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to dsvd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a
 *   n = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
*/


// Changed these names from MAX and SIGN to make them less likely to overlap.
#define MAXXY(x,y) ((x)>(y)?(x):(y))
#define SIGNAB(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) 


static double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}

int dsvd(double **a, int m, int n, double *w, double **v)
{
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
//    double *rv1;
  
    if (m < n) 
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }
  
    //rv1 = (double *)malloc((unsigned int) n*sizeof(double));
                double rv1[n];
/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) 
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
            for (k = i; k < m; k++) 
                scale += fabs(a[k][i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k][i] = (a[k][i]/scale);
                    s += (a[k][i] * a[k][i]);
                }
                f = a[i][i];
                g = -SIGNAB(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (f - g);
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += (a[k][i] * a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k][j] += (f * a[k][i]);
                    }
                }
                for (k = i; k < m; k++) 
                    a[k][i] = (a[k][i]*scale);
            }
        }
        w[i] = (float)(scale * g);
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs(a[i][k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i][k] = (a[i][k]/scale);
                    s += (a[i][k] * a[i][k]);
                }
                f = a[i][l];
                g = -SIGNAB(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (f - g);
                for (k = l; k < n; k++) 
                    rv1[k] = a[i][k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += (a[j][k] * a[i][k]);
                        for (k = l; k < n; k++) 
                            a[j][k] += (s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++) 
                    a[i][k] = (a[i][k]*scale);
            }
        }
        anorm = MAXXY(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++)
                    v[j][i] = ((a[i][j] / a[i][l]) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += (a[i][k] * v[k][j]);
                    for (k = l; k < n; k++) 
                        v[k][j] += (s * v[k][i]);
                }
            }
            for (j = l; j < n; j++) 
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = w[i];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i][j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += (a[k][i] * a[k][j]);
                    f = (s / a[i][i]) * g;
                    for (k = i; k < m; k++) 
                        a[k][j] += (f * a[k][i]);
                }
            }
            for (j = i; j < m; j++) 
                a[j][i] = (a[j][i]*g);
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }
    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                        {
                            y = a[j][nm];
                            z = a[j][i];
                            a[j][nm] = (y * c + z * s);
                            a[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < n; j++) 
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
//                free((void*) rv1);
                fprintf(stderr, "No convergence after 30,000! iterations \n");
                return(0);
            }
    
            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGNAB(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) 
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = (y * c + z * s);
                    a[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    //free( rv1);
    return(1);
}
