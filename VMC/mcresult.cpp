#include "mcresult.h"
#include <cmath>     // sqrt
#include <algorithm> // std::max
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
