/* $Id: GaussianPeakShape.cpp 2310 2009-09-01 15:29:43Z bkausler $ */

/*
 * GaussianPeakShape.cpp
 *
 * Copyright (c) 2009 Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 *
 * This file is part of ms++.
 *
 * ms++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ms++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ms++. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <cmath>

#include <ms++/Error.h>
#include <ms++/PeakShape.h>

using namespace ms;

double GaussianPeakShape::at(const double xCoordinate) const {
    return std::exp(-(xCoordinate * xCoordinate) / (2 * sigma_ * sigma_));
}

double GaussianPeakShape::getSupportThreshold() const {
    return this->getSigma() * this->getSigmaFactorForSupportThreshold();
}

// construction
GaussianPeakShape::GaussianPeakShape(const double sigma, const double sigmaFactorForSupportThreshold) 
    : sigma_(sigma), sigmaFactorForSupportThreshold_(sigmaFactorForSupportThreshold) {
    mspp_precondition(sigma > 0, "GaussianPeakShape::GaussianPeakShape(): sigma has to be positive.");
    mspp_precondition(sigmaFactorForSupportThreshold > 0, "GaussianPeakShape::GaussianPeakShape(): sigmaFactorForSupportThreshold has to be positive.");
}


// setter/getter
void GaussianPeakShape::setSigma(const double sigma) {
    mspp_precondition(sigma > 0, "GaussianPeakShape::GaussianPeakShape(): Parameter sigma has to be positive."); 
    sigma_ = sigma; 
}


void GaussianPeakShape::setFwhm(const double fwhm) {
    mspp_precondition(fwhm > 0, "GaussianPeakShape::GaussianPeakShape(): Parameter fwhm has to be positive.");
    sigma_ = fwhm / sigmaToFwhmConversionFactor();
}
double GaussianPeakShape::getFwhm() const {
    return sigma_ * sigmaToFwhmConversionFactor();
}

void GaussianPeakShape::setSigmaFactorForSupportThreshold(const double factor) {
    mspp_precondition(factor > 0, "GaussianPeakShape::GaussianPeakShape(): sigmaFactorForSupportThreshold has to be positive.");
    sigmaFactorForSupportThreshold_ = factor;
}
double GaussianPeakShape::getSigmaFactorForSupportThreshold() const {
    return sigmaFactorForSupportThreshold_;
}

double GaussianPeakShape::sigmaToFwhmConversionFactor() const {
    return 2 * sqrt(2 * std::log(2.));
}
