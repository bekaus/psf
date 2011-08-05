#include <cmath>

#include <psf/PeakParameter.h>

using namespace psf;

unsigned int QuadraticModel::numberOfParameters() {
    return numberOfParameters_;
}

void QuadraticModel::setParameter(unsigned index, double value) {
    psf_precondition(index < numberOfParameters(), "QuadraticModel::setParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        a_ = value;
    }
    else {
        b_ = value;    
    }
}
double QuadraticModel::getParameter(unsigned index) {
    psf_precondition(index < numberOfParameters(), "QuadraticModel::getParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        return a_;
    }
    else {
        return b_;    
    }
}

double QuadraticModel::at(const double x) const {
    return a_ * x*x + b_;
}

GeneralizedSlope QuadraticModel::slopeInParameterSpaceFor(double x) const {
    double slope[] = {x * x, 1., 0.};
    return GeneralizedSlope(slope, slope + 3);
}

// setter / getter
void QuadraticModel::setA(const double a) {
    a_ = a;
}
double QuadraticModel::getA() const {
    return a_;
}

void QuadraticModel::setB(const double b) {
    b_ = b;
}
double QuadraticModel::getB() const {
    return b_;
}
