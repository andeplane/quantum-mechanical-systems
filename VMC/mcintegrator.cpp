#include "mcintegrator.h"
#include "random.h"
#include <cmath>


WaveFunction MCIntegrator::waveFunction() const
{
    return m_waveFunction;
}

void MCIntegrator::setWaveFunction(const WaveFunction &waveFunction)
{
    m_waveFunction = waveFunction;
}
MCIntegrator::MCIntegrator()
{
    setDefaultLocalEnergy();
    setDefaultQuantumForce();
}

MCIntegrator::~MCIntegrator()
{
    m_positions.clear();
}

MCResult MCIntegrator::integrate(unsigned int numberOfMCCycles)
{
    if(!m_parameters || !m_trialFunction || !m_localEnergy || !m_positions.size()) {
        std::cout << "Integrator not initialized. Missing at least one of: parameters, trialFunction, localEnergy or empty positions array, aborting." << std::endl;
        exit(1);
    }

    MCResult result;
    m_numberOfAcceptedMoves = 0;
    double energy = m_localEnergy(m_parameters, m_positions);

    const unsigned int numberOfParticles = this->numberOfParticles();

    vector<vec3> quantumForcesOld(numberOfParticles);
    vector<vec3> quantumForcesNew(numberOfParticles);
    m_quantumForce(m_parameters, m_positions, quantumForcesOld);

    for(unsigned int cycle=0; cycle<numberOfMCCycles; cycle++) {
        double waveFunctionOld = m_trialFunction(m_parameters, m_positions);

        for(unsigned int particle1=0; particle1 < numberOfParticles; particle1++) {
            vec3 oldPos = m_positions[particle1]; // Keep old position in case step is rejected
            m_walkerFunction(m_parameters, m_positions, quantumForcesOld, particle1); // Choose new position
            double waveFunctionNew = m_trialFunction(m_parameters, m_positions); // Compute wave function value in this position
            double greensFunctionAcceptanceRatioFactor = 1.0; // Will be used if Fokker-Planck importance sampling is enabled

            if(m_parameters->importanceSampling) {
                // Compute quantum force from the Fokker-Planck Langevin
                m_quantumForce(m_parameters, m_positions, quantumForcesNew);
                vec3 deltaR = m_positions[particle1];
                deltaR.subtract(oldPos);
                greensFunctionAcceptanceRatioFactor = 0;
                greensFunctionAcceptanceRatioFactor = m_parameters->diffusionConstant*m_parameters->timestep*0.25*(quantumForcesOld[particle1].lengthSquared() - quantumForcesNew[particle1].lengthSquared()) - 0.5*deltaR.dot(quantumForcesOld[particle1] + quantumForcesNew[particle1]);
//                for(int j=0; j<3; j++) {
//                    int i = particle1;
//                    greensFunctionAcceptanceRatioFactor += 0.5*(quantumForcesOld[i][j]+quantumForcesNew[i][j])*(m_parameters->diffusionConstant*m_parameters->timestep*0.5*(quantumForcesOld[i][j]-quantumForcesNew[i][j])-m_positions[i][j]+oldPos[j]);
//                }
                greensFunctionAcceptanceRatioFactor = exp(greensFunctionAcceptanceRatioFactor);
            }

            if(Random::nextDouble() < greensFunctionAcceptanceRatioFactor*waveFunctionNew*waveFunctionNew / (waveFunctionOld*waveFunctionOld)) {
                quantumForcesOld = quantumForcesNew;
                m_numberOfAcceptedMoves++;
                waveFunctionOld = waveFunctionNew;
                energy = m_localEnergy(m_parameters, m_positions);
            } else {
                m_positions[particle1] = oldPos;
            }

            result.addDataPoint(energy);
        }
    }

    m_numberOfAcceptedMoves /= m_positions.size();

    return result;
}

void MCIntegrator::setDefaultLocalEnergy() {
    m_localEnergyNumerical = [&](MCIntegratorParameters *parameters, vector<vec3> &positions) {
        const double waveFunctionCurrent = m_trialFunction(parameters, positions);
        const unsigned int numberOfParticles = positions.size();

        const double h = 0.001;
        const double oneOverHSquared = 1e6;
        // Kinetic energy
        double kineticEnergy = 0;
        for(unsigned int i = 0; i < numberOfParticles; i++) {
            for(unsigned int j = 0; j < 3; j++) {
                // Increase a small change in j'th direction
                positions[i][j] += h;
                double waveFunctionPlus = m_trialFunction(parameters, positions);
                positions[i][j] -= 2.0*h;
                double waveFunctionMinus = m_trialFunction(parameters, positions);

                // This is the standard u'' ≈ (f+ - 2f + f-)/h^2
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

    m_localEnergy = m_localEnergyNumerical;
}

void MCIntegrator::setDefaultQuantumForce()
{
    m_quantumForceNumerical = [&](MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces) {
        vector<vec3> positionsPlus = positions;
        vector<vec3> positionsMinus = positions;
        double waveFunction = m_trialFunction(parameters, positions);
        double oneOverWaveFunction = 1.0 / waveFunction;

        double h = 0.000001;
        double oneOverTwoTimesH = 1.0 / (2.0*h);
        const unsigned int numberOfParticles = positions.size();
        for(unsigned int particle=0; particle<numberOfParticles; particle++) {
            for(int dimension=0; dimension<3; dimension++) {
                positionsPlus[particle][dimension] += h;
                positionsMinus[particle][dimension] -= h;

                // Numerical derivative: y' ≈ (y(x+h) - y(x-h)) / (2*h)
                double derivative = (m_trialFunction(parameters, positionsPlus) - m_trialFunction(parameters, positionsMinus)) * oneOverTwoTimesH;

                // Quantum force is 2*grad(psi) / psi
                quantumForces[particle][dimension] = 2*oneOverWaveFunction*derivative;

                positionsPlus[particle][dimension] -= h;
                positionsMinus[particle][dimension] += h;
            }
        }
    };

    m_quantumForce = m_quantumForceNumerical;
}

function<double (MCIntegratorParameters *parameters, vector<vec3> &positions)> MCIntegrator::trialFunction() const
{
    return m_trialFunction;
}

void MCIntegrator::setTrialFunction(const function<double(MCIntegratorParameters *parameters, vector<vec3> &positions)> &trialFunction)
{
    m_trialFunction = trialFunction;
}

function<double (MCIntegratorParameters *parameters, vector<vec3> &positions)> MCIntegrator::localEnergy() const
{
    return m_localEnergy;
}

void MCIntegrator::setLocalEnergy(const function<double(MCIntegratorParameters *parameters, vector<vec3> &positions)> &localEnergy)
{
    m_localEnergy = localEnergy;
}

unsigned int MCIntegrator::numberOfParticles() const
{
    return m_positions.size();
}

unsigned int MCIntegrator::numberOfAcceptedMoves() const
{
    return m_numberOfAcceptedMoves;
}


vector<vec3> MCIntegrator::positions() const
{
    return m_positions;
}

void MCIntegrator::setPositions(const vector<vec3> &positions)
{
    m_positions = positions;
}

MCIntegratorParameters *MCIntegrator::parameters() const
{
    return m_parameters;
}

void MCIntegrator::setParameters(MCIntegratorParameters *parameters)
{
    m_parameters = parameters;
}

function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces)> MCIntegrator::quantumForce() const
{
    return m_quantumForce;
}

void MCIntegrator::setQuantumForce(const function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces)> &quantumForce)
{
    m_quantumForce = quantumForce;
}

function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces, int particleIndex)> MCIntegrator::walkerFunction() const
{
    return m_walkerFunction;
}

void MCIntegrator::setWalkerFunction(const function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces, int particleIndex)> &walkerFunction)
{
    m_walkerFunction = walkerFunction;
}
