#include <cmath>

#include <ms++/Error.h>
#include <psf/PeakShape.h>

using namespace psf;

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
