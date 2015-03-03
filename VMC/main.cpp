#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include "random.h"
#include "vec3.h"
#include "mcintegrator.h"
#include "wavefunction.h"
using namespace std;

void helium() {
    MCIntegrator integrator;
    integrator.setWaveFunction(WaveFunctions::Helium());

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
    unsigned int numberOfCycles = 1e7;

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

