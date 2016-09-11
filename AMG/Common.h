#ifndef ___COMMON_H___
#define ___COMMON_H___

#include <vector>
#include <list>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#pragma warning(disable: 4267)

namespace AlgebraicMultigrid
{
    typedef int                         integer;
    typedef double                      real;
    typedef std::vector<integer>        intVec;
    typedef std::vector<real>           reaVec;
    typedef std::pair<integer, real>    pair_i_r;

    inline integer poisson(
        integer n,
        intVec    &ptr, intVec &ind, reaVec    &val, reaVec    &rhs
        )
    {
        integer    vecleng = n * n;
        real       h = 1.0 / (n - 1);

        ptr.clear(); ptr.reserve(vecleng + 1); ptr.push_back(0);
        ind.clear(); ind.reserve(vecleng * 5); 
        val.clear(); val.reserve(vecleng * 5);

        rhs.resize(vecleng);

        for (integer j = 0, k = 0; j < n; ++j) {
            for (integer i = 0; i < n; ++i, ++k) {
                if (i == 0 || i == n - 1 || j == 0 || j == n - 1) {
                    ind.push_back(k);
                    val.push_back(1.0);

                    rhs[k] = 0.0;
                }
                else {
                    ind.push_back(k - n);
                    val.push_back(-1.0 / (h * h));

                    ind.push_back(k - 1);
                    val.push_back(-1.0 / (h * h));

                    ind.push_back(k);
                    val.push_back(4.0 / (h * h));

                    ind.push_back(k + 1);
                    val.push_back(-1.0 / (h * h));

                    ind.push_back(k + n);
                    val.push_back(-1.0 / (h * h));

                    rhs[k] = 1.0;
                }

                ptr.push_back(ind.size());
            }
        }
        return vecleng;
    }
}






#endif