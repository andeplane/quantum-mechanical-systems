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

    integrator.setWalkerFunction([&](MCIntegratorParameters *parameters, vector<vec3> &positions, int particleIndex) {
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

    // Parameters
    unsigned int numberOfCycles = 1e7;
    double timestep = 0.05;
    double D = 0.5;
    double variance = 2*D*timestep;

    double standardDeviation = sqrt(variance);
    HeliumParameters parameters;

    vector<vec3> initialPositions(2);
    initialPositions[0].randomGaussian(0, standardDeviation);
    initialPositions[1].randomGaussian(0, standardDeviation);

    integrator.setParameters(&parameters);
    integrator.setInitialPositions(initialPositions);

    auto trial1 = [&](MCIntegratorParameters *parameters, vector<vec3> &positions) {
        HeliumParameters *params = (HeliumParameters *)parameters;
        double r1 = positions[0].length();
        double r2 = positions[1].length();

        return exp(-params->alpha*(r1 + r2));
    };

    auto trial2 = [&](MCIntegratorParameters *parameters, vector<vec3> &positions) {
        HeliumParameters *params = (HeliumParameters *)parameters;
        vec3 deltaR = positions[0];
        deltaR.subtract(positions[1]);
        double r1 = positions[0].length();
        double r2 = positions[1].length();
        double r12 = deltaR.length();

        return exp(-params->alpha*(r1 + r2))*exp(0.5*r12/((1 + params->beta*r12)));
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

    auto driftF = [&](MCIntegratorParameters *parameters, vector<vec3> &positions, int particleIndex) {
        HeliumParameters *params = (HeliumParameters *)parameters;

        vec3 deltaR = positions[particleIndex];
        int particle2Index = (particleIndex + 1) & 1; // Fastmod exploiting that if b=2^n, a % b = a & (b-1). In this case b = 2, so a%2 = a&1.
        deltaR.subtract(positions[particle2Index]);
        double r12 = deltaR.length();
        deltaR.normalize();
        return deltaR/((1.0 + params->beta*r12)*(1.0 + params->beta*r12)) - 2.0*params->alpha*positions[particleIndex].normalized();
    };

    auto driftFNum = [&](MCIntegratorParameters *parameters, vector<vec3> &positions, int particleIndex) {
        vec3 F;
        vector<vec3> positionsPlus = positions;
        vector<vec3> positionsMinus = positions;
        double waveFunction = trial2(parameters, positions);
        double oneOverWaveFunction = 1.0 / waveFunction;

        double h = 0.000001;
        for(int i=0; i<3; i++) {
            positionsPlus[particleIndex][i] += h;
            positionsMinus[particleIndex][i] -= h;
            double derivative = (trial2(parameters, positionsPlus) - trial2(parameters, positionsMinus)) / (2*h);
            F[i] = 2*oneOverWaveFunction*derivative;

            positionsPlus[particleIndex][i] -= h;
            positionsMinus[particleIndex][i] += h;
        }
        return F;
    };

    integrator.setTrialFunction(trial2);

    integrator.setWalkerFunction([&](MCIntegratorParameters *parameters, vector<vec3> &positions, int particleIndex) {
        // stepLength = variance = 2*D*dt
        // sqrt(stepLength) = standardDeviation = sqrt(2*D*dt)
        // In importance sampling, we want to add D*F*dt = 0.5*stepLength*F
        vec3 deltaR;
        deltaR.randomGaussian(0.0, standardDeviation);
        vec3 F = driftF(parameters, positions, particleIndex);
//        vec3 FNum = driftFNum(parameters, positions, particleIndex);
//        vec3 delta = F - FNum;
//        cout << "deltaR=" << deltaR << " whereas F=" << F << " and FNum=" << FNum << " which gives deltaNorm=" << delta.length() << endl;
        // deltaR += 0.5*variance*F;

        positions[particleIndex].add(deltaR);
        // cout << "Now at " << positions[particleIndex] << " (r=" << positions[particleIndex].length() << ")" << endl;
    });

    ofstream file("output.txt");

    integrator.setTrialFunction(trial2);
    unsigned int iMax = 100;
    unsigned int jMax = 100;
    for(unsigned int i=0; i<iMax; i++) {
        // double alpha = 1.4 + 0.025*i;
        double min = 1.80;
        double max = 1.90;
        double delta = (max - min)/iMax;

        double alpha = min + delta*i;
        for(unsigned int j=0; j<jMax; j++) {
            min = 0.2;
            max = 0.5;
            delta = (max - min)/jMax;
            double beta = min + delta*j;
            parameters.alpha = alpha;
            parameters.beta = beta;

            parameters.alpha = 1.843;
            parameters.beta = 0.347;

            MCResult result = integrator.integrate(numberOfCycles);
            cout << "Result for (alpha, beta)=(" << alpha << ", " << beta << "): " << result.mean() << " (acceptance: " << integrator.numberOfAcceptedMoves() / double(numberOfCycles) << ", variance: " << result.variance() << ", sigma/sqrt(N): " << result.standardDeviation() / sqrt(result.numberOfSamples()) << ") " << endl;
            file << alpha << " " << beta << " " << result.mean() << " " << result.variance() << " " << result.standardDeviation() / sqrt(result.numberOfSamples()) << endl;
            return;
        }
    }
}

void heliumAlphaBetaKnown() {
    MCIntegrator integrator;

    // Parameters
    unsigned int numberOfCycles = 1e7;
    double timestep = 0.05;
    double D = 0.5;
    double variance = 2*D*timestep;

    double standardDeviation = sqrt(variance);
    HeliumParameters parameters;

    vector<vec3> initialPositions(2);
    initialPositions[0].randomGaussian(0, standardDeviation);
    initialPositions[1].randomGaussian(0, standardDeviation);

    integrator.setParameters(&parameters);
    integrator.setInitialPositions(initialPositions);

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
    auto driftF = [&](MCIntegratorParameters *parameters, vector<vec3> &positions, int particleIndex) {
        HeliumParameters *params = (HeliumParameters *)parameters;

        vec3 deltaR = positions[particleIndex];
        int particle2Index = (particleIndex + 1) & 1; // Fastmod exploiting that if b=2^n, a % b = a & (b-1). In this case b = 2, so a%2 = a&1.
        deltaR.subtract(positions[particle2Index]);
        double r12 = deltaR.length();
        deltaR.normalize();
        return deltaR/((1.0 + params->beta*r12)*(1.0 + params->beta*r12)) - 2.0*params->alpha*positions[particleIndex].normalized();
    };

    integrator.setTrialFunction(trial2);

    integrator.setWalkerFunction([&](MCIntegratorParameters *parameters, vector<vec3> &positions, int particleIndex) {
        // stepLength = variance = 2*D*dt
        // sqrt(stepLength) = standardDeviation = sqrt(2*D*dt)
        // In importance sampling, we want to add D*F*dt = 0.5*stepLength*F
        vec3 deltaR;
        deltaR.randomGaussian(0.0, standardDeviation);
        vec3 F = driftF(parameters, positions, particleIndex);
        deltaR += 0.5*variance*F;

        // This gives rNew = rOld + gaussian*2*Ddt + D*dt*F (eq 14.24 in MHJ Computational Physics, 2014)
        positions[particleIndex].add(deltaR);
    });

    integrator.setTrialFunction(trial2);
    parameters.alpha = 1.843;
    parameters.beta = 0.347;
    MCResult result = integrator.integrate(numberOfCycles);
    cout << "Result for (alpha, beta)=(" << parameters.alpha << ", " << parameters.beta << "): " << result.mean() << " (acceptance: " << integrator.numberOfAcceptedMoves() / double(numberOfCycles) << ", variance: " << result.variance() << ", sigma/sqrt(N): " << result.standardDeviation() / sqrt(result.numberOfSamples()) << ") " << endl;
}

int main()
{
    // hydrogen();
    // helium();

    heliumAlphaBetaKnown();

    return 0;
}

