#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
 int i, j, k, m, n, r;
 double *datar1, *datar2, *datal1, * datal2;
 double temp;
 if (nrhs != 2)
 mexErrMsgTxt("need and only need 2 inputs.");



    /* Find the dimensions of the data */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    datar1 = mxGetPr(prhs[0]);
    datar2 = mxGetPr(prhs[1]);
    r = (int)(*datar2);    
//     r = int (&prhs[1]);

    /* Create an mxArray for the output data */
    plhs[0] = mxCreateDoubleMatrix(r, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(r, n, mxREAL);

//     /* Retrieve the input data */

//     /* Create a pointer to the output data */
    datal1 = mxGetPr(plhs[0]);
    datal2 = mxGetPr(plhs[1]);
//     
    for (i = 0; i < r; i++)   
    for (k = 0; k < n; k++)
    {
        datal1[i +  k * r] = 1000;
        
    }

     /* Put data in the output array */
       
    for (j = 0; j < n; j++)
    for (i = 0; i < m; i++) 
    {
//         
        if  (datal1[(r-1) +  j * r ] > datar1[i + j * m])
         {
             
            temp = datar1[i + j * m];
            datal1[(r-1) + j * r] = temp;
            datal2[(r-1) + j * r] = (double)i;
            for (k = r-2; k >-1; k--)
            {
                if (temp < datal1[k + j * r] && k >0)
                {
                   datal1[k+1 +  j * r ] = datal1[k + j * r];
                   datal2[k+1 + j * r] = datal2[k + j * r];
                }
                else if (temp < datal1[k + j * r] && k ==0)
                {
                    datal1[k+1 +  j * r ] = datal1[k + j * r];
                    datal1[k +  j * r ] = temp;
                    datal2[k+1 +  j * r ] = datal2[k + j * r];
                    datal2[k +  j * r ] = (double)(i+ 1);
                    break;
                }
                else
                {
                   datal1[k+1 +  j * r ] = temp;
                   datal2[k+1 +  j * r ] = (double)(i+ 1);
                   break;
                }
            }
//             
        }
    }
}
   
   