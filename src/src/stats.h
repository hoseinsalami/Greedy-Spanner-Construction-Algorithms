#ifndef STATS_H
#define STATS_H

#include <ostream>

class stats;
class stats {
public:
    stats();
    virtual ~stats();
    virtual void add(double value, double weight=1);
    virtual double mean() const ;
    virtual double variance() const;
    virtual double radius95p() const;

protected:
    double sum, sum2, sumweight;
};

std::ostream & operator << (std::ostream &o, const stats &s);

#endif
