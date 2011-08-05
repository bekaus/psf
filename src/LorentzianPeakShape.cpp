/* $Id$ */

/*
 * LorentzianPeakShape.cpp
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
#include "psf/PeakShape.h"

using namespace psf;

double LorentzianPeakShape::at(const double xCoordinate) const {
    return fwhm_ / ((xCoordinate * xCoordinate) + (fwhm_*fwhm_));
}

double LorentzianPeakShape::getSupportThreshold() const {
    return this->getFwhm() * this->getFwhmFactorForSupportThreshold();
}

// construction
LorentzianPeakShape::LorentzianPeakShape(const double fwhm, const double fwhmFactorForSupportThreshold) 
    : fwhm_(fwhm), fwhmFactorForSupportThreshold_(fwhmFactorForSupportThreshold) {
    mspp_precondition(fwhm > 0, "LorentzianPeakShape::LorentzianPeakShape(): Parameter fwhm has to be positive.");
    mspp_precondition(fwhmFactorForSupportThreshold > 0, "LorentzianPeakShape::LorentzianPeakShape(): fwhmFactorForSupportThreshold has to be positive.");
}


// setter/getter
void LorentzianPeakShape::setFwhm(const double fwhm) {
    mspp_precondition(fwhm > 0, "LorentzianPeakShape::LorentzianPeakShape(): Parameter fwhm has to be positive.");
    fwhm_ = fwhm;
}
double LorentzianPeakShape::getFwhm() const {
    return fwhm_;
}

void LorentzianPeakShape::setFwhmFactorForSupportThreshold(const double factor) {
    mspp_precondition(factor > 0, "LorentzianPeakShape::LorentzianPeakShape(): Parameter fwhmFactorForSupportThreshold has to be positive.");
    fwhmFactorForSupportThreshold_ = factor;
}
double LorentzianPeakShape::getFwhmFactorForSupportThreshold() const {
    return fwhmFactorForSupportThreshold_;
}
