#ifndef MCRESULT_H
#define MCRESULT_H


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

#endif // MCRESULT_H
