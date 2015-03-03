#include "wavefunction.h"
#include "mcintegrator.h"
#include <iostream>

WaveFunction::WaveFunction()
{
    setNumericalLocalEnergy();
    setNumericalGradient();
    evaluate = [&](MCIntegratorParameters *parameters, vector<vec3> &positions) {
        std::cout << "Evaluate function not implemented in WaveFunction object. Aborting!" << std::endl;
        exit(1);
        return 0;
    };
}

WaveFunction::~WaveFunction()
{

}

void WaveFunction::setNumericalGradient()
{
    gradient = [&](MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &gradient) {
        vector<vec3> positionsPlus = positions;
        vector<vec3> positionsMinus = positions;

        double h = 0.000001;
        double oneOverTwoTimesH = 1.0 / (2.0*h);
        const unsigned int numberOfParticles = positions.size();
        for(unsigned int particle=0; particle<numberOfParticles; particle++) {
            for(int dimension=0; dimension<3; dimension++) {
                positionsPlus[particle][dimension] += h;
                positionsMinus[particle][dimension] -= h;

                // Numerical derivative: y' ≈ (y(x+h) - y(x-h)) / (2*h)
                double derivative = (evaluate(parameters, positionsPlus) - evaluate(parameters, positionsMinus)) * oneOverTwoTimesH;

                gradient[particle][dimension] = derivative;
                positionsPlus[particle][dimension] -= h;
                positionsMinus[particle][dimension] += h;
            }
        }
    };
}

void WaveFunction::setNumericalLocalEnergy()
{
    localEnergy = [&](MCIntegratorParameters *parameters, vector<vec3> &positions) {
        const double waveFunctionCurrent = evaluate(parameters, positions);

        const double h = 0.001;
        const double oneOverHSquared = 1e6;

        // Kinetic energy
        const unsigned int numberOfParticles = positions.size();
        double kineticEnergy = 0;
        for(unsigned int i = 0; i < numberOfParticles; i++) {
            for(unsigned int j = 0; j < 3; j++) {
                // Increase a small change in j'th direction
                positions[i][j] += h;
                double waveFunctionPlus = evaluate(parameters, positions);
                positions[i][j] -= 2.0*h;
                double waveFunctionMinus = evaluate(parameters, positions);

                // This is the standard u'' ≈ ( f(x+h) - 2*f(x) + f(x-h) ) / h^2
                kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);

                positions[i][j] += h;
            }
        }
        kineticEnergy = 0.5 * oneOverHSquared * kineticEnergy / waveFunctionCurrent;

        // Potential energy
        double potentialEnergy = 0;

        for(unsigned int i = 0; i < numberOfParticles; i++) {
            // Coulomb repulsion
            potentialEnergy -= parameters->charge / positions[i].length();

            // Contribution from electron-electron potential
            for(unsigned int j = i + 1; j < numberOfParticles; j++) {
                vec3 deltaR = positions[i];
                deltaR.subtract(positions[j]);
                potentialEnergy += 1.0 / deltaR.length();
            }
        }

        return kineticEnergy + potentialEnergy;
    };
}

