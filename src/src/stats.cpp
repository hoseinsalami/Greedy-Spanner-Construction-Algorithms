#include <iostream>
#include <cstdio>
#include <cmath>

#include "stats.h"

/**
 * \class stats
 * Used for keeping track of statistics.
 * Calculates the mean and variance of all elements add()ed.
 */

/**
 * Creates a new and empty stats object with the given name.
 */
stats::stats() :
    sum(0),
    sum2(0),
    sumweight(0) {
}

/**
 * Destroys the stats object.
 */
stats::~stats() {
}

/**
 * Adds a value to this stats object. 
 * 
 * \param value the value that is added.
 * \param weight the weight for this \c value used to calculate the mean.
 */
void stats::add(double value, double weight) {
    sum+=value * weight;
    sum2+=value * value * weight;
    sumweight+=weight;
}

/**
 * Return the current mean.
 */
double stats::mean() const {
    if (sumweight<1e-9) return 0;
    return sum/sumweight;
}

/**
 * Return the current variance. 
 */
double stats::variance() const {
    if (sumweight<1+1e-9) return 0;
    return (sum2 - sum * sum / sumweight) / sumweight;
}

double stats::radius95p() const {
    if (sumweight<1+1e-9) return 0;
    double sd_mu=std::sqrt(variance()/sumweight);
    return sd_mu*1.959964; // R: qnorm(.025) = qnorm(.975) = 1.959964
}

/**
 * Output the state of this stats object.
 */
std::ostream & operator << (std::ostream &o, const stats &s) {
    char tmp[160];
    double m=s.mean();
    double r=s.radius95p();
    //std::sprintf();
	sprintf(tmp, "Mean=%lf, SD=%lf, 95%% confidence interval=[%lf,%lf]", m, sqrt(s.variance()), m - r, m + r);
    o << tmp;
    return o;
}
