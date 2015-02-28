#pragma once
#include "vec3.h"

#include <functional>
#include <vector>
using std::function;
using std::vector;

class MCResult
{
private:
    unsigned int m_numSamples;
    double m_sum;
    double m_M;
    double m_S;
public:
    MCResult();
    void addDataPoint(double value);
    unsigned int numberOfSamples();
    double mean();
    double variance();
    double standardDeviation();
};

class MCIntegratorParameters
{
public:
    double charge = 1.0;
};

class MCIntegrator
{
private:
    MCIntegratorParameters *m_parameters;
    unsigned int m_numberOfAcceptedMoves;
    vector<vec3> m_initialPositions;
    function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, int particleIndex)> m_walkerFunction;
    function<double(MCIntegratorParameters *parameters, vector<vec3> &positions)> m_trialFunction;
    function<double(MCIntegratorParameters *parameters, vector<vec3> &positions)> m_localEnergy;
    function<double(MCIntegratorParameters *parameters, vector<vec3> &positions)> m_localEnergyNumerical;
public:
    MCIntegrator();
    ~MCIntegrator();


    MCResult integrate(unsigned int numberOfMCCycles);
    function<double (MCIntegratorParameters *parameters, vector<vec3> &positions)> trialFunction() const;
    void setTrialFunction(const function<double (MCIntegratorParameters *parameters, vector<vec3> &positions)> &trialFunction);
    function<double (MCIntegratorParameters *parameters, vector<vec3> &positions)> localEnergy() const;
    void setLocalEnergy(const function<double (MCIntegratorParameters *parameters, vector<vec3> &positions)> &localEnergy);
    function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, int particleIndex)> walkerFunction() const;
    void setWalkerFunction(const function<void(MCIntegratorParameters *parameters, vector<vec3> &positions, int particleIndex)> &walkerFunction);
    unsigned int numberOfParticles() const;
    unsigned int numberOfAcceptedMoves() const;
    void setNumberOfAcceptedMoves(unsigned int numberOfAcceptedMoves);
    vector<vec3> initialPositions() const;
    void setInitialPositions(const vector<vec3> &initialPositions);
    MCIntegratorParameters *parameters() const;
    void setParameters(MCIntegratorParameters *parameters);
    void setDefaultLocalEnergy();
};
