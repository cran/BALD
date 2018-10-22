##################################################################################################
##                                                                                              ##
##    BALD is an R-package.                                                                     ##
##    It is a Bayesian time series model of loss development.                                   ##
##    Features include skewed Student-t distribution with time-varying scale parameters,        ##
##    an expert prior for the calendar year effect,                                             ##
##    and accommodation for structural breaks in the consumption path of development years.     ##
##    It is an update for the older package lossDev as it has been stopped supported.           ##
##                                                                                              ##
##    Copyright (c) 2018 Frank A. Schmid,                                                       ##
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


##' A package for robust stochastic loss development.
##'
##' \pkg{BALD} makes available a Bayesian time series model of loss development, estimated using MCMC.
##' This package is intended used for Property and Casualty (PC) actuaries.
##' PC actuaries aggregate historical loss experience in a triangle two dimensions form.
##' The two dimensions are accident years and development years. Most of the time, the loss experience is paid loss and/or incurred loss.
##' Incurred loss is the sum of paid loss and case reserve put up for the estimate of future payment of the claims.
##' 
##' Accident Years is on the row dimension. It is the year when the incident of the claims happened.
##' 
##' Development Year is on the column dimension. It is the time when the claims payments progresses.
##' 
##' Calendar Year is on the diagonal of the triangle form. 
##' It is by the calendar year view of loss experience and often times the insurance companies reporting balance sheet and income statement views.
##' 
##' The features of this package include skewed Student-t distribution with time-varying scale parameter for the likelihood of the loss distribution with given parameters.
##' User can put expert prior for the calendar year effect, and structural break in the consumption path of services/development years.
##' This is an update for the older package lossDev as it has been stopped supported
##' 
##' Please read the vignette for guidance on usage.  And \dQuote{Robust Loss Development Using MCMC} by Schmid, Frank A. for theory.
##'
##' @references Schmid, Frank A., \dQuote{Robust Loss Development Using MCMC,}, 2009
##' @name BALD
##' @docType package
NULL
