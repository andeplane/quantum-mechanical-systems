#pragma once
#include "vec3.h"
#include "mcresult.h"
#include "wavefunction.h"
#include "mcintegratorparameters.h"
#include <functional>
#include <vector>
#include <cmath>
using std::function;
using std::vector;

class MCIntegrator
{
private:
    MCIntegratorParameters *m_parameters = NULL;
    unsigned int m_numberOfAcceptedMoves = 0;
    vector<vec3> m_positions;
    vector<vec3> m_quantumForcesOld;
    vector<vec3> m_quantumForcesNew;
    WaveFunction m_waveFunction;

    function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces, int particleIndex)> m_walkerFunction;
    void computeQuantumForces(vector<vec3> &quantumForces);
public:
    MCIntegrator();
    ~MCIntegrator();
    MCResult integrate(unsigned int numberOfMCCycles);
    function<double (MCIntegratorParameters *parameters, vector<vec3> &positions)> trialFunction() const;
    void setTrialFunction(const function<double (MCIntegratorParameters *parameters, vector<vec3> &positions)> &trialFunction);
    function<double (MCIntegratorParameters *parameters, vector<vec3> &positions)> localEnergy() const;
    void setLocalEnergy(const function<double (MCIntegratorParameters *parameters, vector<vec3> &positions)> &localEnergy);
    unsigned int numberOfParticles() const;
    unsigned int numberOfAcceptedMoves() const;
    void setNumberOfAcceptedMoves(unsigned int numberOfAcceptedMoves);
    vector<vec3> positions() const;
    void setPositions(const vector<vec3> &positions);
    MCIntegratorParameters *parameters() const;
    void setParameters(MCIntegratorParameters *parameters);
    function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces, int particleIndex)> walkerFunction() const;
    void setWalkerFunction(const function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces, int particleIndex)> &walkerFunction);
    WaveFunction waveFunction() const;
    void setWaveFunction(const WaveFunction &waveFunction);
};
