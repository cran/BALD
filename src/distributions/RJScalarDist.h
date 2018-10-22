/*
##################################################################################################
##                                                                                              ##
##    BALD is an R-package.                                                                     ##
##    It is a Bayesian time series model of loss development.                                   ##
##    Features include skewed Student-t distribution with time-varying scale parameters,        ##
##    an expert prior for the calendar year effect,                                             ##
##    and accommodation for structural breaks in the consumption path of development years.     ##
##    It is an update for the older package lossDev as it has been stopped supported.           ##
##                                                                                              ##
##    Copyright © 2018 Frank A. Schmid,                                                         ##
##                                                                                              ##
##    This file is part of BALD.                                                                ##
##                                                                                              ##
##    lossDev is free software: you can redistribute it and/or modify                           ##
##    it under the terms of the GNU General Public License as published by                      ##
##    the Free Software Foundation, either version 3 of the License, or                         ##
##    (at your option) any later version.                                                       ##
##                                                                                              ##
##    This program is distributed in the hope that it will be useful,                           ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of                            ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                             ##
##    GNU General Public License for more details.                                              ##
##                                                                                              ##
##    You should have received a copy of the GNU General Public License                         ##
##    along with this program.  If not, see <https://www.gnu.org/licenses/>.                    ##
##                                                                                              ##
##################################################################################################
*/


#ifndef RJ_SCALAR_DIST_H_
#define RJ_SCALAR_DIST_H_

#include <distribution/ScalarDist.h>

namespace jags {

struct RNG;

/**
 * @short Scalar Distribution using R math library infrastructure.
 *
 * A subclass of RScalarDist has to implement the d,p,q, and r virtual
 * member functions. These are based on the d-p-q-r functions provided
 * by libRmath.
 *
 * The JAGS versions of most (but not all) scalar distributions extend
 * the distribution families in libRmath by allowing the distribution
 * to be bounded.
 */
class RJScalarDist : public ScalarDist
{
    const Support _support;
    const bool _discrete;
    unsigned int _npar;
    double calPlower(double, std::vector<double const *> const &) const;
    double calPupper(double, std::vector<double const *> const &) const;
public:
    /**
     * Constructor
     *
     * @param name BUGS language name of distribution
     *
     * @param npar Number of parameters, excluding upper and lower bound
     *
     * @param support Support of distribution
     *
     * @param discrete Boolean flag indicating whether the distribution is
     *        discrete valued.
     */
    RJScalarDist(std::string const &name, unsigned int npar, Support support,
		bool discrete=false);
    double logDensity(double x, PDFType type,
			 std::vector<double const *> const &parameters,
			 double const *lower, double const *upper) const;
    double randomSample(std::vector<double const *> const &parameters,
			double const *lower, double const *upper,
			RNG *rng) const;
    /**
     * Returns the median. Note that this function can be overloaded
     * by a subclass if necessary.
     */
    double typicalValue(std::vector<double const *> const &parameters,
			double const *lower, double const *upper) const;
    /**
     * Density function, ignoring bounds
     * @param x value at which to evaluate the density
     * @param parameters Array of parameters
     * @param give_log Indicates whether to return log density.
     */
    virtual double d(double x, std::vector<double const *> const &parameters,
		     bool give_log) const = 0;
    /**
     * Distribution function, ignoring bounds
     * @param x quantile at which to evaluate the distribution function
     * @param parameters Array of parameters
     * @param lower If true, return value is P[X <= x]. Otherwise
     * P[X > x]
     * @param give_log Indicates whether to return log probabability
     */
    virtual double p(double x, std::vector<double const *> const &parameters,
		     bool lower, bool give_log) const = 0;
    /**
     * Quantile function, ignoring bounds
     * @param p probability for which to evaluate quantile
     * @param parameters Array of parameters
     * @param log_p Indicates whether p is given as log(p).
     */
    virtual double q(double p, std::vector<double const *> const &parameters,
		     bool lower, bool log_p) const = 0;
    /**
     * Random number generation, ignoring bounds
     * @param parameters Array of parameters
     */
    virtual double
	r(std::vector<double const *> const &parameters, RNG *rng) const = 0;
    /**
     * All RScalarDist distributions can be bounded
     */
    bool canBound() const;
    /**
     * RScalarDist distributions are defined to have support on the integers
     * or on the real line by the constructor
     */
    bool isDiscreteValued(std::vector<bool> const &mask) const;
    /**
     * Alternative function for determining whether the distribution is
     * discrete-valued.
     */
    bool discrete() const;
    /**
     * Returns the number of parameters of the distribution
     */
    unsigned int npar() const;
};
}
#endif /* SCALAR_DIST_RMATH_H_ */
