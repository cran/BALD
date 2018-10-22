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


#ifndef DTOV_H_
#define DTOV_H_

#include "RJScalarDist.h"

/**
 * t-distribution on k degrees of freedom, with median mu and
 * scale parameter tau.
 * <pre>
 * f(x|mu, tau, k)
 * f(x|0,1,k) = Gamma((k+1)/2) / (sqrt(k*pi) Gamma(k/2)) (1 + x^2/k)^-((k+1)/2)
 * </pre>
 * @short t distribution
 */

namespace jags {
class DTOV : public RJScalarDist {
 public:
  DTOV();

  double d(double x, std::vector<double const *> const &parameters,
	   bool log) const;
  double p(double x, std::vector<double const *> const &parameters, bool lower,
	   bool log) const;
  double q(double x, std::vector<double const *> const &parameters, bool lower,
	   bool log) const;
  double r(std::vector<double const *> const &parameters, RNG *rng) const;
  /**
   * Check that tau > 0 and k > 0
   */
  bool checkParameterValue(std::vector<double const *> const &parameters) const;

};
}
#endif /* DTOV_H_ */
