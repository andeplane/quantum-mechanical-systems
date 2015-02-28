#include "mcintegrator.h"
#include "random.h"
#include <cmath>

MCResult MCIntegrator::integrate(unsigned int numberOfMCCycles)
{
    if(!m_parameters || !m_trialFunction || !m_localEnergy || !m_positions.size()) {
        std::cout << "Integrator not initialized. Missing at least one of: parameters, trialFunction, localEnergy or numberOfParticles, aborting." << std::endl;
        exit(1);
    }
    MCResult result;
    m_numberOfAcceptedMoves = 0;
    double energy = m_localEnergy(m_parameters, m_positions);

    for(unsigned int cycle=0; cycle<numberOfMCCycles; cycle++) {
        double waveFunctionOld = m_trialFunction(m_parameters, m_positions);

        unsigned int numberOfParticles = m_positions.size();
        for(unsigned int particle=0; particle < numberOfParticles; particle++) {
            vec3 oldPos = m_positions[particle];
            m_walkerFunction(m_parameters, m_positions, particle);

            double waveFunctionNew = m_trialFunction(m_parameters, m_positions);

            if(true || Random::nextDouble() < waveFunctionNew*waveFunctionNew / (waveFunctionOld*waveFunctionOld)) {
                m_numberOfAcceptedMoves++;
                waveFunctionOld = waveFunctionNew;
                energy = m_localEnergy(m_parameters, m_positions);
            } else {
                m_positions[particle] = oldPos;
            }

            result.addDataPoint(energy);
        }
    }

    m_numberOfAcceptedMoves /= m_positions.size();

    return result;
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

function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, int particleIndex)> MCIntegrator::walkerFunction() const
{
    return m_walkerFunction;
}

void MCIntegrator::setWalkerFunction(const function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, int particleIndex)> &walkerFunction)
{
    m_walkerFunction = walkerFunction;
}
MCIntegrator::MCIntegrator() :
    m_parameters(NULL),
    m_numberOfAcceptedMoves(0)
{
    setDefaultLocalEnergy();
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

MCIntegrator::~MCIntegrator()
{

}



