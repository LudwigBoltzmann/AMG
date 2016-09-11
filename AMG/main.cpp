#include "Common.h"
#include "ConjugateGradient.h"
#include "Matrix.h"

int main(int argc, char* argv[])
{
    using namespace AlgebraicMultigrid;

    if (argc == 1) {
        std::cerr << "needs arguments (example : \" 100 100 \")" << std::endl;
        return -1;
    }

    integer dim = 0, vecleng;
    dim = atoi(argv[1]);

    intVec ptr;
    intVec ind;
    reaVec val;
    reaVec x;
    reaVec b;

    vecleng = poisson(dim, ptr, ind, val, b);
    
    



    return 0;
}