#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include "vec3.h"
#include <functional>
#include <vector>
using std::function;
using std::vector;
class MCIntegratorParameters;

class WaveFunction
{
private:

public:
    WaveFunction();
    ~WaveFunction();
    function<double(MCIntegratorParameters *parameters, vector<vec3> &positions)> evaluate;
    function<double(MCIntegratorParameters *parameters, vector<vec3> &positions)> localEnergy;
    function<void(MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &gradient)> gradient;
    void setNumericalGradient();
    void setNumericalLocalEnergy();
};

#endif // WAVEFUNCTION_H
