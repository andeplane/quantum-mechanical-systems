#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

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
    double beta = 0.7;
    HeliumParameters() {
        charge = 2.0;
    }
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
    initialPositions[0][1] = 0;
    initialPositions[0][2] = 0;

    integrator.setParameters(&parameters);
    integrator.setPositions(initialPositions);

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

    integrator.setWalkerFunction([&](MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces, int particleIndex) {
        double step = stepLength*(Random::nextDouble() - 0.5);

        positions[particleIndex][0] += step;
        // Negative positions aren't allowed for radial coordinates
        while(positions[particleIndex][0] < 0) {
            step = stepLength*(Random::nextDouble() - 0.5);
            positions[particleIndex][0] += step;
        }
    });

    vector<double> alphas = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};

    for(auto alpha : alphas) {
        parameters.alpha = alpha;
        MCResult result = integrator.integrate(numberOfCycles);
        cout << "Result for alpha=" << alpha << ": " << result.mean() << " (acceptance: " << integrator.numberOfAcceptedMoves() / double(numberOfCycles) << ", variance: " << result.variance() << ", sigma/sqrt(N): " << result.standardDeviation() / sqrt(result.numberOfSamples()) << ") " << endl;
    }
}

void helium() {
    MCIntegrator integrator;

    auto trial2 = [&](MCIntegratorParameters *parameters, vector<vec3> &positions) {
        HeliumParameters *params = (HeliumParameters *)parameters;
        vec3 deltaR = positions[0];
        deltaR.subtract(positions[1]);
        double r1 = positions[0].length();
        double r2 = positions[1].length();
        double r12 = deltaR.length();

        return exp(-params->alpha*(r1 + r2))*exp(r12/(2*(1 + params->beta*r12)));
    };


    integrator.setLocalEnergy([&](MCIntegratorParameters *parameters, vector<vec3> &positions) {
        HeliumParameters *params = (HeliumParameters *)parameters;
        double r1DotR2 = positions[0].dot(positions[1]);
        double r1 = positions[0].length();
        double r2 = positions[1].length();
        vec3 deltaR = positions[0];
        deltaR.subtract(positions[1]);
        double r12 = deltaR.length();
        double oneOverR12 = 1.0 / r12;
        double oneOverR1 = 1.0 / r1;
        double oneOverR2 = 1.0 / r2;

        double EL1 = (params->alpha - params->charge)*(oneOverR1 + oneOverR2) + oneOverR12 - params->alpha*params->alpha;
        double oneOverTwoTimeOnePlusBetaR12Sq = 0.5/( (1.0 + params->beta*r12)*(1.0 + params->beta*r12));

        double EL2 = EL1 + oneOverTwoTimeOnePlusBetaR12Sq*(params->alpha*(r1 + r2)*oneOverR12*(1.0 - r1DotR2*oneOverR1*oneOverR2) - oneOverTwoTimeOnePlusBetaR12Sq - 2.0*oneOverR12 + 2.0*params->beta/(1.0 + params->beta*r12));
        return EL2;
    });

    // This is the drift term, F = 2*1/Psi * grad(Psi)
    auto driftF = [&](MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces) {
        HeliumParameters *params = (HeliumParameters *)parameters;

        vec3 deltaR = positions[0];
        deltaR.subtract(positions[1]);
        double r12 = deltaR.length();
        deltaR.normalize();
        quantumForces[0] = deltaR/((1.0 + params->beta*r12)*(1.0 + params->beta*r12)) - 2.0*params->alpha*positions[0].normalized();
        quantumForces[1] = -deltaR/((1.0 + params->beta*r12)*(1.0 + params->beta*r12)) - 2.0*params->alpha*positions[1].normalized();
    };

    integrator.setTrialFunction(trial2);
    integrator.setQuantumForce(driftF);

    integrator.setWalkerFunction([&](MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces, int particleIndex) {
        // stepLength = variance = 2*D*dt
        // sqrt(stepLength) = standardDeviation = sqrt(2*D*dt)
        // In importance sampling, we want to add D*F*dt = 0.5*variance*F
        vec3 deltaR;
        deltaR.randomGaussian(0.0, parameters->standardDeviation);
        deltaR += 0.5*parameters->variance*quantumForces[particleIndex];

        // This gives rNew = rOld + gaussian*sqrt(2*Ddt) + D*dt*F (eq 14.24 in MHJ Computational Physics, 2014)
        positions[particleIndex].add(deltaR);
    });

    // Parameters
    HeliumParameters parameters;
    unsigned int numberOfCycles = 1e6;

    parameters.timestep = 0.1;
    parameters.diffusionConstant = 0.5;
    parameters.variance = 2*parameters.diffusionConstant*parameters.timestep;
    parameters.standardDeviation = sqrt(parameters.variance);
    parameters.importanceSampling = true;

    vector<vec3> initialPositions(2);
    initialPositions[0].randomGaussian(0, parameters.standardDeviation);
    initialPositions[1].randomGaussian(0, parameters.standardDeviation);

    integrator.setParameters(&parameters);
    integrator.setPositions(initialPositions);
    integrator.setTrialFunction(trial2);
    parameters.alpha = 1.843;
    parameters.beta = 0.347;

    ofstream file("output.txt");
    unsigned int iMax = 5;
    unsigned int jMax = 5;
    for(unsigned int i=0; i<iMax; i++) {
        double min = 1.83;
        double max = 1.85;
        double delta = (max - min)/iMax;

        double alpha = min + delta*i;
        for(unsigned int j=0; j<jMax; j++) {
            min = 0.34;
            max = 0.35;
            delta = (max - min)/jMax;
            double beta = min + delta*j;
            parameters.alpha = alpha;
            parameters.beta = beta;
            MCResult result = integrator.integrate(numberOfCycles);
            cout << "Result for (alpha, beta)=(" << alpha << ", " << beta << "): " << result.mean() << " (acceptance: " << integrator.numberOfAcceptedMoves() / double(numberOfCycles) << ", variance: " << result.variance() << ", sigma/sqrt(N): " << result.standardDeviation() / sqrt(result.numberOfSamples()) << ") " << endl;
            file << alpha << " " << beta << " " << result.mean() << " " << result.variance() << " " << result.standardDeviation() / sqrt(result.numberOfSamples()) << endl;
        }
    }

    MCResult result = integrator.integrate(numberOfCycles);
    cout << "Result for (alpha, beta)=(" << parameters.alpha << ", " << parameters.beta << "): " << result.mean() << " (acceptance: " << integrator.numberOfAcceptedMoves() / double(numberOfCycles) << ", variance: " << result.variance() << ", sigma/sqrt(N): " << result.standardDeviation() / sqrt(result.numberOfSamples()) << ") " << endl;
}

int main()
{
    // hydrogen();
    helium();

    return 0;
}

