#include "Multislice.cuh"
namespace PRISM{
    void buildMultisliceOutput_gpu(Parameters <PRISM_FLOAT_PRECISION> &pars,
                                   Array3D <complex<PRISM_FLOAT_PRECISION>> &trans,
                                   Array2D <complex<PRISM_FLOAT_PRECISION>> &PsiProbeInit,
                                   Array2D <PRISM_FLOAT_PRECISION> &alphaInd) {
        using namespace std;
        cout << "test gpu " << endl;

    };
}