#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include "vec3.h"
#include "mcintegratorparameters.h"
#include <functional>
#include <vector>
using std::function;
using std::vector;

class WaveFunction
{
public:
    WaveFunction();
    ~WaveFunction();
    function<double(MCIntegratorParameters *parameters, vector<vec3> &positions)> evaluate;
    function<double(MCIntegratorParameters *parameters, vector<vec3> &positions)> localEnergy;
    function<void(MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &gradient)> gradient;
    void setNumericalGradient();
    void setNumericalLocalEnergy();
};

struct WaveFunctions
{
public:
    static WaveFunction Helium();
};

class HydrogenParameters : public MCIntegratorParameters
{
public:
    double alpha = 0.7;
};

class HeliumParameters : public MCIntegratorParameters
{
public:
    double alpha = 0.7;
    double beta = 0.7;
    HeliumParameters() {
        charge = 2.0;
    }
};
#endif // WAVEFUNCTION_H
