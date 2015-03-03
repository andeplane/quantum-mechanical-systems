#ifndef MCINTEGRATORPARAMETERS_H
#define MCINTEGRATORPARAMETERS_H
#include <cmath>

class MCIntegratorParameters
{
public:
    double timestep = 1.0;
    double diffusionConstant = 0.5;
    double charge = 1.0;
    double variance = 2*diffusionConstant*timestep;
    double standardDeviation = sqrt(variance);
    bool   importanceSampling = false;
};

#endif // MCINTEGRATORPARAMETERS_H
