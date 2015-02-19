#include <iostream>
#include <cmath>
#include <vector>

#include "random.h"
#include "vec3.h"
#include "mcintegrator.h"
using namespace std;

class HydrogenParameters : public MCIntegratorParameters
{
public:
    double alpha = 0.7;
};

class HeliumParameters : public MCIntegratorParameters
{
public:
    double alpha = 0.7;
};

void hydrogen() {
    MCIntegrator integrator;

    // Parameters
    // unsigned int numberOfCycles = 5e5;
    unsigned int numberOfCycles = 1e6;
    double stepLength = 0.1;
    HydrogenParameters parameters;

    vector<vec3> initialPositions(1);
    initialPositions[0].randomUniform(0, stepLength);

    integrator.setParameters(&parameters);
    integrator.setInitialPositions(initialPositions);

    integrator.setTrialFunction([&](MCIntegratorParameters *parameters, vector<vec3> &positions) {
        HydrogenParameters *params = (HydrogenParameters *)parameters;
        double rho = positions[0][0];

        return params->alpha*rho*exp(-params->alpha*rho);
    });

//    integrator.setLocalEnergy([&](MCIntegratorParameters *parameters, vector<vec3> &positions) {
//        HydrogenParameters *params = (HydrogenParameters *)parameters;
//        double rho = positions[0][0];
//        return 1.0/rho*(params->alpha - 1.0) - 0.5*params->alpha*params->alpha;;
//    });

    integrator.setWalkerFunction([&](vec3 &position) {
        double step = stepLength*(Random::nextDouble() - 0.5);

        position[0] += step;
        // Negative positions aren't allowed for radial coordinates
        while(position[0] < 0) {
            step = stepLength*(Random::nextDouble() - 0.5);
            position[0] += step;
        }
    });

    vector<double> alphas = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};

    for(auto alpha : alphas) {
        parameters.alpha = alpha;
        MCResult result = integrator.integrate(numberOfCycles);
        cout << "Result for alpha=" << alpha << ": " << result.mean() << " (acceptance: " << integrator.numberOfAcceptedMoves() / double(numberOfCycles) << ", variance: " << result.variance() << ", sigma/N: " << result.standardDeviation() / sqrt(result.numberOfSamples()) << ") " << endl;
    }
}

//void helium() {
//    MCIntegrator integrator;

//    // Parameters
//    unsigned int numberOfCycles = 1e8;
//    double stepLength = 0.1;
//    HeliumParameters parameters;

//    vector<vec3> initialPositions(2);
//    initialPositions[0].randomUniform(0, stepLength);

//    integrator.setParameters(&parameters);
//    integrator.setInitialPositions(initialPositions);

//    integrator.setTrialFunction([&](MCIntegratorParameters *parameters, vector<vec3> &positions) {
//        HydrogenParameters *params = (HydrogenParameters *)parameters;
//        double rho = positions[0][0];

//        return params->alpha*rho*exp(-params->alpha*rho);
//    });

//    integrator.setLocalEnergy([&](MCIntegratorParameters *parameters, vector<vec3> &positions) {
//        HydrogenParameters *params = (HydrogenParameters *)parameters;
//        double rho = positions[0][0];
//        return 1.0/rho*(params->alpha - 1.0) - 0.5*params->alpha*params->alpha;;
//    });

//    integrator.setWalkerFunction([&](vec3 &position) {
//        double step = stepLength*(Random::nextDouble() - 0.5);

//        position[0] += step;
//        // Negative positions aren't allowed for radial coordinates
//        while(position[0] < 0) {
//            step = stepLength*(Random::nextDouble() - 0.5);
//            position[0] += step;
//        }
//    });

//    vector<double> alphas = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};

//    for(auto alpha : alphas) {
//        parameters.alpha = alpha;
//        MCResult result = integrator.integrate(numberOfCycles);
//        cout << "Result for alpha=" << alpha << ": " << result.mean() << " (acceptance: " << integrator.numberOfAcceptedMoves() / double(numberOfCycles) << ", variance: " << result.variance() << ", sigma/N: " << result.standardDeviation() / sqrt(result.numberOfSamples()) << ") " << endl;
//    }
//}

int main()
{
    hydrogen();

    return 0;
}

