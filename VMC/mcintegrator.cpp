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
void MCIntegrator::computeQuantumForces(vector<vec3> &quantumForces)
{
    double waveFunctionValue = m_waveFunction.evaluate(m_parameters, m_positions);
    double twoOverWaveFunctionValue = 2.0/waveFunctionValue;
    m_waveFunction.gradient(m_parameters, m_positions, quantumForces);

    const unsigned int numberOfParticles = m_positions.size();
    for(unsigned int particle=0; particle<numberOfParticles; particle++) {
        quantumForces[particle] *= twoOverWaveFunctionValue;
    }
}

MCIntegrator::MCIntegrator()
{

}

MCIntegrator::~MCIntegrator()
{
    m_positions.clear();
}

MCResult MCIntegrator::integrate(unsigned int numberOfMCCycles)
{
    if(!m_parameters || !m_positions.size()) {
        std::cout << "Integrator not initialized. Missing parameters or empty positions array, aborting." << std::endl;
        exit(1);
    }

    MCResult result;
    m_numberOfAcceptedMoves = 0;
    double energy = m_waveFunction.localEnergy(m_parameters, m_positions);

    const unsigned int numberOfParticles = this->numberOfParticles();

    m_quantumForcesOld.resize(numberOfParticles);
    m_quantumForcesNew.resize(numberOfParticles);
    if(m_parameters->importanceSampling) computeQuantumForces(m_quantumForcesOld);

    for(unsigned int cycle=0; cycle<numberOfMCCycles; cycle++) {
        double waveFunctionOld = m_waveFunction.evaluate(m_parameters, m_positions);

        for(unsigned int particle1=0; particle1 < numberOfParticles; particle1++) {
            vec3 oldPos = m_positions[particle1]; // Keep old position in case step is rejected
            m_walkerFunction(m_parameters, m_positions, m_quantumForcesOld, particle1); // Choose new position
            double waveFunctionNew = m_waveFunction.evaluate(m_parameters, m_positions); // Compute wave function value in this position
            double greensFunctionAcceptanceRatioFactor = 1.0; // Will be used if Fokker-Planck importance sampling is enabled

            if(m_parameters->importanceSampling) {
                // Compute quantum force from the Fokker-Planck Langevin
                computeQuantumForces(m_quantumForcesNew);
                vec3 deltaR = m_positions[particle1];
                deltaR.subtract(oldPos);
                greensFunctionAcceptanceRatioFactor = 0;
                greensFunctionAcceptanceRatioFactor = m_parameters->diffusionConstant*m_parameters->timestep*0.25*(m_quantumForcesOld[particle1].lengthSquared() - m_quantumForcesNew[particle1].lengthSquared()) - 0.5*deltaR.dot(m_quantumForcesOld[particle1] + m_quantumForcesNew[particle1]);
                greensFunctionAcceptanceRatioFactor = exp(greensFunctionAcceptanceRatioFactor);
            }

            if(Random::nextDouble() < greensFunctionAcceptanceRatioFactor*waveFunctionNew*waveFunctionNew / (waveFunctionOld*waveFunctionOld)) {
                // Random walk move is accepted
                m_quantumForcesOld = m_quantumForcesNew;
                m_numberOfAcceptedMoves++;
                waveFunctionOld = waveFunctionNew;
                energy = m_waveFunction.localEnergy(m_parameters, m_positions);
            } else {
                // Random walk move is rejected
                m_positions[particle1] = oldPos;
            }

            result.addDataPoint(energy);
        }
    }

    m_numberOfAcceptedMoves /= m_positions.size();

    return result;
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

function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces, int particleIndex)> MCIntegrator::walkerFunction() const
{
    return m_walkerFunction;
}

void MCIntegrator::setWalkerFunction(const function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces, int particleIndex)> &walkerFunction)
{
    m_walkerFunction = walkerFunction;
}
