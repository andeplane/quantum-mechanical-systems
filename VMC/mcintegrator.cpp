#include "mcintegrator.h"
#include "random.h"

#include <cmath>

// Calculating variance by http://www.johndcook.com/blog/standard_deviation/
MCResult::MCResult() :
    m_numSamples(0),
    m_sum(0),
    m_M(0),
    m_S(0)
{

}

void MCResult::addDataPoint(double value)
{
    m_sum += value;
    m_numSamples++;

    double lastM = m_M;
    m_M += (value - lastM)/m_numSamples;
    m_S += (value - lastM)*(value - m_M);
}

unsigned int MCResult::numberOfSamples()
{
    return m_numSamples;
}

double MCResult::mean()
{
    return m_sum / std::max((unsigned int)1, m_numSamples);
}

double MCResult::variance()
{
    // var(x) = <x^2> - <x>^2
    // Calculating variance by http://www.johndcook.com/blog/standard_deviation/
    return m_S/(m_numSamples - 1);
}

double MCResult::standardDeviation()
{
    return sqrt(variance());
}

MCResult MCIntegrator::integrate(unsigned int numberOfMCCycles)
{
    if(!m_parameters || !m_trialFunction || !m_localEnergy || !m_initialPositions.size()) {
        std::cout << "Integrator not initialized. Missing at least one of: parameters, trialFunction, localEnergy or numberOfParticles, aborting." << std::endl;
        exit(1);
    }
    MCResult result;
    m_numberOfAcceptedMoves = 0;
    vector<vec3> positions = m_initialPositions;
    double energy = m_localEnergy(m_parameters, positions);

    for(unsigned int cycle=0; cycle<numberOfMCCycles; cycle++) {
        double waveFunctionOld = m_trialFunction(m_parameters, positions);

        unsigned int numberOfParticles = positions.size();
        for(unsigned int particle=0; particle < numberOfParticles; particle++) {
            vec3 oldPos = positions[particle];
            m_walkerFunction(positions[particle]);

            double waveFunctionNew = m_trialFunction(m_parameters, positions);

            if(Random::nextDouble() < waveFunctionNew*waveFunctionNew / (waveFunctionOld*waveFunctionOld)) {
                m_numberOfAcceptedMoves++;
                waveFunctionOld = waveFunctionNew;
                energy = m_localEnergy(m_parameters, positions);
            } else {
                positions[particle] = oldPos;
            }

            result.addDataPoint(energy);
        }
    }

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
    return m_initialPositions.size();
}

unsigned int MCIntegrator::numberOfAcceptedMoves() const
{
    return m_numberOfAcceptedMoves;
}


vector<vec3> MCIntegrator::initialPositions() const
{
    return m_initialPositions;
}

void MCIntegrator::setInitialPositions(const vector<vec3> &initialPositions)
{
    m_initialPositions = initialPositions;
}

MCIntegratorParameters *MCIntegrator::parameters() const
{
    return m_parameters;
}

void MCIntegrator::setParameters(MCIntegratorParameters *parameters)
{
    m_parameters = parameters;
}

function<void (vec3 &position)> MCIntegrator::walkerFunction() const
{
    return m_walkerFunction;
}

void MCIntegrator::setWalkerFunction(const function<void (vec3 &position)> &walkerFunction)
{
    m_walkerFunction = walkerFunction;
}
MCIntegrator::MCIntegrator() :
    m_parameters(NULL),
    m_numberOfAcceptedMoves(0)
{

}

MCIntegrator::~MCIntegrator()
{

}


