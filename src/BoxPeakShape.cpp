#include <cmath>

#include <ms++/Error.h>
#include "psf/PeakShape.h"

using namespace psf;

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
