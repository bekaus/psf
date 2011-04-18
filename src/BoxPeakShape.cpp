/* $Id: BoxPeakShape.cpp 2310 2009-09-01 15:29:43Z bkausler $ */

/*
 * BoxPeakShape.cpp
 *
 * Copyright (c) 2009 Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * Copyright (c) 2009 Marc Kirchner <marc.kirchner@childrens.harvard.edu>
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

double BoxPeakShape::at(const double xCoordinate) const {
    // this is the only difference between the Box and tha Gaussian
    return 1.0;
}

double BoxPeakShape::getSupportThreshold() const {
    return this->getSigma() * this->getSigmaFactorForSupportThreshold();
}

// construction
BoxPeakShape::BoxPeakShape(const double sigma, const double sigmaFactorForSupportThreshold) 
    : sigma_(sigma), sigmaFactorForSupportThreshold_(sigmaFactorForSupportThreshold) {
    mspp_precondition(sigma > 0, "BoxPeakShape::BoxPeakShape(): sigma has to be positive.");
    mspp_precondition(sigmaFactorForSupportThreshold > 0, "BoxPeakShape::BoxPeakShape(): sigmaFactorForSupportThreshold has to be positive.");
}


// setter/getter
void BoxPeakShape::setSigma(const double sigma) {
    mspp_precondition(sigma > 0, "BoxPeakShape::BoxPeakShape(): Parameter sigma has to be positive."); 
    sigma_ = sigma; 
}


void BoxPeakShape::setFwhm(const double fwhm) {
    mspp_precondition(fwhm > 0, "BoxPeakShape::BoxPeakShape(): Parameter fwhm has to be positive.");
    sigma_ = fwhm / sigmaToFwhmConversionFactor();
}
double BoxPeakShape::getFwhm() const {
    return sigma_ * sigmaToFwhmConversionFactor();
}

void BoxPeakShape::setSigmaFactorForSupportThreshold(const double factor) {
    mspp_precondition(factor > 0, "BoxPeakShape::BoxPeakShape(): sigmaFactorForSupportThreshold has to be positive.");
    sigmaFactorForSupportThreshold_ = factor;
}
double BoxPeakShape::getSigmaFactorForSupportThreshold() const {
    return sigmaFactorForSupportThreshold_;
}

double BoxPeakShape::sigmaToFwhmConversionFactor() const {
    return 2 * sqrt(2 * std::log(2.));
}
