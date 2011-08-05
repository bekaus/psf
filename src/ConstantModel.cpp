#include <psf/Error.h>
#include "psf/PeakParameter.h"

namespace psf
{
unsigned int ConstantModel::numberOfParameters() {
    return numberOfParameters_;
}

void ConstantModel::setParameter(unsigned index, double value) {
    psf_precondition(index < numberOfParameters(), "ConstantModel::setParameter(): Parameter index out-of-range.");
    a_ = value;
}
double ConstantModel::getParameter(unsigned index) {
    psf_precondition(index < numberOfParameters(), "ConstantModel::getParameter(): Parameter index out-of-range.");
    return a_;
}

double ConstantModel::at(const double x) const {
    return a_;
}

GeneralizedSlope ConstantModel::slopeInParameterSpaceFor(double x) const {
    double slope[] = {1., 0.};
    return GeneralizedSlope(slope, slope + 2);
}

// setter / getter
void ConstantModel::setA(const double a) {
    a_ = a;
}
double ConstantModel::getA() const {
    return a_;
}

} /* namespace psf */
