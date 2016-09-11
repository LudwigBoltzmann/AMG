#ifndef ___MATRIXOPERATION_H___
#define ___MATRIXOPERATION_H___

#include "Common.h"

namespace AlgebraicMultigrid
{
    inline void spmv(integer dim, intVec& ptr, intVec& ind, reaVec& val, reaVec& x, reaVec& res)
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

    inline void spmv(integer dim, real alpha, intVec& ptr, intVec& ind, reaVec& val, reaVec& x, real beta, reaVec& res)
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

    inline void residual(integer dim, intVec& ptr, intVec& ind, reaVec& val, reaVec& x, reaVec& b, reaVec& res)
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

    inline real norm(reaVec& x)
    {
        double sum = 0.0;
        for (auto it = x.begin(); it != x.end(); it++) {
            sum += *it * *it;
        }
        return sqrt(sum);
    }

    inline real innerProduct(reaVec& x, reaVec& y)
    {
        integer dim = x.size();
        real sum = 0.0;
#pragma omp parallel
        {
            real partial_sum = 0.0;
#pragma omp for
            for (int i = 0; i < dim; i++) {
                partial_sum += x[i] * y[i];
            }
#pragma omp atomic 
            sum += partial_sum;
#pragma omp barrier
        }
        return sum;
        

    }

}




#endif