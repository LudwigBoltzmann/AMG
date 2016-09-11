#ifndef ___MATRIXOPERATION_H___
#define ___MATRIXOPERATION_H___

#include "Common.h"

namespace AlgebraicMultigrid
{
    void spmv(integer dim, intVec& ptr, intVec& ind, reaVec& val, reaVec& x, reaVec& res)
    {
#pragma omp parallel for
        for (integer i = 0; i < dim; i++) {
            integer start = ptr[i];
            integer end = ptr[i + 1];
            real    sum = 0.0;
            for (integer j = start; j < end; j++) {
                sum += x[ind[j]] * val[j];
            }
            res[i] = sum;
        }
    }

    void spmv(integer dim, real alpha, intVec& ptr, intVec& ind, reaVec& val, reaVec& x, real beta, reaVec& res)
    {
#pragma omp parallel for
        for (integer i = 0; i < dim; i++) {
            integer start = ptr[i];
            integer end = ptr[i + 1];
            real sum = 0.0;
            for (integer j = start; j < end; j++) {
                sum += x[ind[j]] * val[j];
            }
            res[i] = alpha * sum + beta * res[i];
        }
    }

    void residual(integer dim, intVec& ptr, intVec& ind, reaVec& val, reaVec& x, reaVec& b, reaVec& res)
    {
#pragma omp parallel for
        for (integer i = 0; i < dim; i++) {
            integer start = ptr[i];
            integer end = ptr[i + 1];
            real sum = 0.0;
            for (integer j = start; j < end; j++) {
                sum += x[ind[j]] * val[j];
            }
            res[i] = b[i] - sum;
        }
    }

}




#endif