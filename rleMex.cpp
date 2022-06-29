
#include "mex.h"
#include <string.h>
#include <iostream>
#include <vector>
using namespace std;
//#include "MatlabDataArray.hpp"



void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    const short unsigned int *A  = (uint16_t*) mxGetData(prhs[0]);

    const int n = mxGetN(prhs[0]);
    const int m = mxGetM(prhs[0]);
    short unsigned int i, j, k, count;
    short unsigned int last;
    vector<short unsigned int> vChngCnts, vChngInds, vStart;

    for (i=0;i<n;i++) {
        last=A[i];
        vChngInds.push_back(i+1);
        vStart.push_back(last);
        count=1;
        for (j=1;j<m;j++){
            k=i + j*n;

            if ( A[k] != last ) {
                vChngCnts.push_back(count); // previous
                vChngInds.push_back(k+1);   // next
                last=(int) A[k];
                count=1;
                vStart.push_back(last);
            } else {
                count++;
                if (j==m) {
                    vChngCnts.push_back(count); // previous
                }
            }
        }
        vChngCnts.push_back(count); // previous
    }

    // use populate mx
    mwSize dims0[2] = {vStart.size(),1};
    plhs[0] = mxCreateNumericArray(2, dims0, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(plhs[0]), &(vStart[0]), vStart.size()*sizeof(unsigned short int));


    mwSize dims1[2] = {vChngCnts.size(),1};
    plhs[1] = mxCreateNumericArray(2, dims1, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(plhs[1]), &(vChngCnts[0]), vChngCnts.size()*sizeof(unsigned short int));

    mwSize dims2[2] = {vChngInds.size(),1};
    plhs[2] = mxCreateNumericArray(2, dims2, mxUINT16_CLASS, mxREAL);
    memcpy(mxGetData(plhs[2]), &(vChngInds[0]), vChngInds.size()*sizeof(unsigned short int));

    return;

}
