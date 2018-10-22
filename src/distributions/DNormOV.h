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


#ifndef DNORMOV_H_
#define DNORMOV_H_

#include "RJScalarDist.h"

/**
 * <pre>
 * x ~ dnorm(mu, tau)
 * f(x | mu, tau) = sqrt(tau) * exp(-1/2 * tau * (x - mu)^2)
 * </pre>
 * @short Normal distribution
 */
namespace jags {
class DNormOV : public RJScalarDist {
 public:
  DNormOV();

  double d(double x, std::vector<double const *> const &parameters,
	   bool give_log) const;
  double p(double q, std::vector<double const *> const &parameters, bool lower,
	   bool give_log) const;
  double q(double p, std::vector<double const *> const &parameters, bool lower,
	   bool log_p) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  /**
   * Checks that tau > 0
   */
  bool checkParameterValue(std::vector<double const *> const &parameters) const;
  /**
   * Exploits the capacity to sample truncted normal distributions
   * that is built into the JAGS library, overloading the generic
   * functionality of RScalarDist.
   */
  double randomSample(std::vector<double const *> const &par,
		      double const *lower, double const *upper,
		      RNG *rng) const;
};
}
#endif /* DNORMOV_H_ */
