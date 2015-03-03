#pragma once
#include "vec3.h"
#include "mcresult.h"
#include "wavefunction.h"
#include <functional>
#include <vector>
#include <cmath>
using std::function;
using std::vector;

class MCIntegratorParameters
{
public:
    double timestep = 1.0;
    double diffusionConstant = 0.5;
    double charge = 1.0;
    double variance = 2*diffusionConstant*timestep;
    double standardDeviation = sqrt(variance);
    bool   importanceSampling = false;
};

class MCIntegrator
{
private:
    MCIntegratorParameters *m_parameters = NULL;
    unsigned int m_numberOfAcceptedMoves = 0;
    vector<vec3> m_positions;
    WaveFunction m_waveFunction;

    function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces)> m_quantumForce;
    function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces)> m_quantumForceNumerical;
    function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces, int particleIndex)> m_walkerFunction;
    function<double(MCIntegratorParameters *parameters, vector<vec3> &positions)> m_trialFunction;
    function<double(MCIntegratorParameters *parameters, vector<vec3> &positions)> m_localEnergy;
    function<double(MCIntegratorParameters *parameters, vector<vec3> &positions)> m_localEnergyNumerical;
    vector<vec3> computeQuantumForces(MCIntegratorParameters *parameters, vector<vec3> &positions);
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
    void setDefaultLocalEnergy();
    void setDefaultQuantumForce();
    function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces)> quantumForce() const;
    void setQuantumForce(const function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces)> &quantumForce);
    function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces, int particleIndex)> walkerFunction() const;
    void setWalkerFunction(const function<void (MCIntegratorParameters *parameters, vector<vec3> &positions, vector<vec3> &quantumForces, int particleIndex)> &walkerFunction);
    WaveFunction waveFunction() const;
    void setWaveFunction(const WaveFunction &waveFunction);
};
