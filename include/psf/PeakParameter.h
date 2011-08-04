/* $Id: PeakParameter.h 2671 2009-11-09 21:46:52Z bkausler $ */

/*
 * PeakParameter.h
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

#ifndef __PEAKPARAMETER_H__
#define __PEAKPARAMETER_H__

#include <utility>
#include <vector>

#include <ms++/Error.h>
#include <ms++/Log.h>
#include <ms++/SparseSpectrum.h>
#include <psf/SpectrumAlgorithm.h>

#include <vigra/windows.h>
#include <vigra/matrix.hxx>
#include <vigra/regression.hxx>

namespace ms
{
/**
 * The slope (including the bias) of a multidimensional linear function.
 *
 * Example: For the function @f$ f(x_1,x_2) = ax_1 + bx_2 + c + 3.0 @f$, the generalized slope is
 * @f$ (a,b,c,3.0) @f$ with a bias of 3.0 .
 */
typedef std::vector<double> GeneralizedSlope;


// class ConstantModel
/**
 * @f$ f(x) = a @f$
 *
 * ConstantModel is a policy class. The model depends on one parameter a.
 * @see http://en.wikipedia.org/wiki/Policy_based_design
 * @see ms::PeakParameterModel
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class MSPP_EXPORT ConstantModel
{
public:
    /**
     * This model has one parameter.
     */
    unsigned int numberOfParameters();

    void setParameter(unsigned index, double value);
    double getParameter(unsigned index);

protected:
    /**
     * Value of the model at position x.
     */
    double at(const double x) const;

    /**
     * Protected non-virtual destructor.
     *
     * This is a so called 'policy class'. The host class is inheriting from it to add its
     * functionality. But it is not intended to be used as a base class. In case the user
     * utilizes the policy as a baseclass anyway, we have to prevent calling 'delete' on it
     * (deleting a baseclass with a non-virtual destructor leads to undefined behaviour).
     * So we declare the destructor 'protected'.
     * @see http://en.wikipedia.org/wiki/Policy_based_design
     */
    ~ConstantModel() {}

    /**
     * Equal to @f$ (1.,0.) @f$.
     */
    GeneralizedSlope slopeInParameterSpaceFor(double x) const;

public:
    ConstantModel() : a_(0.1) {}

    void setA(const double a);
    double getA() const;

private:
    static const unsigned numberOfParameters_ = 1;
    double a_;
};


// class LinearSqrtModel()
/**
 * @f$ f(x) = a\cdot x\sqrt{x} + b @f$
 *
 * LinearSqrtModel is a policy class. The model depends on two parameter a and b.
 * @see http://en.wikipedia.org/wiki/Policy_based_design
 * @see ms::PeakParameterModel
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class MSPP_EXPORT LinearSqrtModel
{
public:
    /**
     * This model has two parameters.
     */
    unsigned int numberOfParameters();

    void setParameter(unsigned index, double value);
    double getParameter(unsigned index);

protected:
    /**
     * Value of the model at position x.
     */
    double at(const double x) const;

    /**
     * Protected non-virtual destructor.
     *
     * This is a so called 'policy class'. The host class is inheriting from it to add its
     * functionality. But it is not intended to be used as a base class. In case the user
     * utilizes the policy as a baseclass anyway, we have to prevent calling 'delete' on it
     * (deleting a baseclass with a non-virtual destructor leads to undefined behaviour).
     * So we declare the destructor 'protected'.
     * @see http://en.wikipedia.org/wiki/Policy_based_design
     */
    ~LinearSqrtModel() {}

    /**
     * Equal to @f$ (x\sqrt{x},1.,0.) @f$.
     */
    GeneralizedSlope slopeInParameterSpaceFor(double x) const;

public:
    LinearSqrtModel() : a_(0.1), b_(0.1) {}

    void setA(const double a);
    double getA() const;

    void setB(const double b);
    double getB() const;

private:
    static const unsigned numberOfParameters_ = 2;
    double a_, b_;
};


// class LinearSqrtOriginModel()
/**
 * @f$ f(x) = a\cdot x\sqrt{x}@f$
 *
 * LinearSqrtOriginModel is a policy class. The model depends on one parameter a.
 * @see http://en.wikipedia.org/wiki/Policy_based_design
 * @see ms::PeakParameterModel
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class MSPP_EXPORT LinearSqrtOriginModel
{
public:
    /**
     * This model has two parameters.
     */
    unsigned int numberOfParameters();

    void setParameter(unsigned index, double value);
    double getParameter(unsigned index);

protected:
    /**
     * Value of the model at position x.
     */
    double at(const double x) const;

    /**
     * Protected non-virtual destructor.
     *
     * This is a so called 'policy class'. The host class is inheriting from it to add its
     * functionality. But it is not intended to be used as a base class. In case the user
     * utilizes the policy as a baseclass anyway, we have to prevent calling 'delete' on it
     * (deleting a baseclass with a non-virtual destructor leads to undefined behaviour).
     * So we declare the destructor 'protected'.
     * @see http://en.wikipedia.org/wiki/Policy_based_design
     */
    ~LinearSqrtOriginModel() {}

    /**
     * Equal to @f$ (x\sqrt{x},0.) @f$.
     */
    GeneralizedSlope slopeInParameterSpaceFor(double x) const;

public:
    LinearSqrtOriginModel() : a_(0.1) {}

    void setA(const double a);
    double getA() const;

private:
    static const unsigned numberOfParameters_ = 1;
    double a_;
};


// class SqrtModel
/**
 * @f$ f(x) = a\cdot \sqrt{x} + b @f$
 *
 * SqrtModel is a policy class. The model depends on two parameter a and b.
 * @see http://en.wikipedia.org/wiki/Policy_based_design
 * @see ms::PeakParameterModel
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class MSPP_EXPORT SqrtModel
{
public:
    /**
     * This model has two parameters.
     */
    unsigned int numberOfParameters();

    void setParameter(unsigned index, double value);
    double getParameter(unsigned index);

protected:
    /**
     * Value of the model at position x.
     */
    double at(const double x) const;

    /**
     * Protected non-virtual destructor.
     *
     * This is a so called 'policy class'. The host class is inheriting from it to add its
     * functionality. But it is not intended to be used as a base class. In case the user
     * utilizes the policy as a baseclass anyway, we have to prevent calling 'delete' on it
     * (deleting a baseclass with a non-virtual destructor leads to undefined behaviour).
     * So we declare the destructor 'protected'.
     * @see http://en.wikipedia.org/wiki/Policy_based_design
     */
    ~SqrtModel() {}

    /**
     * Equal to @f$ (\sqrt{x},1.,0.) @f$.
     */    
    GeneralizedSlope slopeInParameterSpaceFor(double x) const;

public:
    SqrtModel() : a_(0.1), b_(0.1) {}

    void setA(const double a);
    double getA() const;

    void setB(const double b);
    double getB() const;

private:
    static const unsigned numberOfParameters_ = 2;
    double a_, b_;
};


// class QuadraticModel
/**
 * @f$ f(x) = a\cdot x^2 + b @f$
 *
 * QuadraticModel is a policy class. The model depends on two parameter a and b.
 * @see http://en.wikipedia.org/wiki/Policy_based_design
 * @see ms::PeakParameterModel
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class MSPP_EXPORT QuadraticModel
{
public:
    /**
     * This model has two parameters.
     */
    unsigned int numberOfParameters();

    void setParameter(unsigned index, double value);
    double getParameter(unsigned index);

protected:
    /**
     * Value of the model at position x.
     */
    double at(const double x) const;

    /**
     * Protected non-virtual destructor.
     *
     * This is a so called 'policy class'. The host class is inheriting from it to add its
     * functionality. But it is not intended to be used as a base class. In case the user
     * utilizes the policy as a baseclass anyway, we have to prevent calling 'delete' on it
     * (deleting a baseclass with a non-virtual destructor leads to undefined behaviour).
     * So we declare the destructor 'protected'.
     * @see http://en.wikipedia.org/wiki/Policy_based_design
     */
    ~QuadraticModel() {}

    /**
     * Equal to @f$ (x^2,1.,0.) @f$.
     */ 
    GeneralizedSlope slopeInParameterSpaceFor(double x) const;

public:
    QuadraticModel() : a_(0.1), b_(0.1) {}

    void setA(const double a);
    double getA() const;

    void setB(const double b);
    double getB() const;

private:
    static const unsigned numberOfParameters_ = 2;
    double a_, b_;
};



// class ParameterModel
/**
 * The interface for the peak parameter policy to be used in ms::PeakParameterFwhm.
 *
 * @attention This interface exists only for the purpose of documentation. Don't inherit
 *      from it, but simply implement it.
 * @see ms::PeakParameterFwhm
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class ParameterModel
{
// mandatory interface
protected:
     /**
     * Value of the model at position x.
     */
     virtual double at(double x) const = 0;

    /**
     * Protected non-virtual destructor.
     *
     * This is a so called 'policy class'. The host class is inheriting from it to add its
     * functionality. But it is not intended to be used as a base class. In case the user
     * utilizes the policy as a baseclass anyway, we have to prevent calling 'delete' on it
     * (deleting a baseclass with a non-virtual destructor leads to undefined behaviour).
     * So we declare the destructor 'protected'.
     * @see http://en.wikipedia.org/wiki/Policy_based_design
     */
     virtual ~ParameterModel() {}

// optional interface
public:
     /**
      * The number of parameters in the model.
      *
      * May be any positive number including zero.
      * This is equivalent to the dimension of parameter space.
      * @attention: Don't implement, if the model has no parameters. This will cause
      * desirable compile time errors, if one tries to learn a model without parameters.
      */
     virtual unsigned int numberOfParameters() = 0;

     /**
      * 0 <= index < this->numberOfParameters() .
      * 
      * @attention: Don't implement, if the model has no parameters. This will cause
      * desirable compile time errors, if one tries to learn a model without parameters.
      * @throw ms::PreconditionViolation Parameter index is out of range or the number of
      *                                  parameters is zero.
      */
     virtual void setParameter(unsigned index, double value) = 0;
     /**
      * 0 <= index < this->numberOfParameters() .
      * 
      * @attention: Don't implement, if the model has no parameters. This will cause
      * desirable compile time errors, if one tries to learn a model without parameters.
      * @throw ms::PreconditionViolation Parameter index is out of range or the number of
      *                                  parameters is zero.
      */
     virtual double getParameter(unsigned index);
     

protected:
    /**
     * The slope of the linear function representing the model in parameter space.
     *
     * In parameter space the coordinates @f$ \vec{x} @f$ and
     * regular parameters @f$ \vec{p} @f$ of the model are switching their roles. Some models may have
     * a linear representation in this space: @f$ model(\vec{p}) = \vec{m}\cdot\vec{p} @f$.
     * @f$ \vec{m} @f$ is the generalized slope including the bias. In a n-dimensional
     * parameter space,  the  generalized slope is (n+1)-dimensional. So, the parameter vector 
     * @f$ \vec{p} @f$  has n+1 elements with the last element always set to unity.
     * Example: For the model @f$ f(x_1,x_2) = ax_1 + bx_2 + c @f$ with parameters a,b and c,
     * the generalized slope in parameter space is
     * @f$ (x_1,x_2,1,0) @f$ and @f$ \vec{p}= (a,b,c,1)^T @f$.
     * @attention This function is part of the optional interface, since not every model may have
     * a linear representation in parameter space.
     *
     * @see ms::PeakParameterModel::GeneralizedSlope
     * 
     * @param x The coordinate x now playing the role of a parameter.
     * @return The generalized slope, a mulitdimensional vector including the bias.
     */
    virtual GeneralizedSlope slopeInParameterSpaceFor(double x) const = 0;
};



// class PeakParameterFwhm
/**
 * 'Full width at half maximum' peak shape parameter.
 * 
 * The FWHM can be used to paramterize a peak shape.
 * 
 * Usually, the FWHM depends on the mass over charge ratio in a mass spectrum. This 
 * dependency can be linear, quadratic or any other function. Therefore, the
 * PeakParameterFwhm is a template, which has to be instantiated with a concrete
 * ParameterModel. Note, that the model should yield non-negative values for non-negative
 * parameters and mz values (There is no such thing as a negative FWHM). Else, you will
 * probably get very unexpected behaviour. 
 *
 * @see ms::ConstantModel, 
 *      ms::LinearSqrtModel,
        ms::LinearSqrtOriginModel
 *      ms::SqrtModel,
 *      ms::QuadraticModel
 * @see ms::PeakShape
 * @see ms::GaussianPeakShape
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
template <typename ParameterModel>
class MSPP_EXPORT PeakParameterFwhm : public ParameterModel 
{
public:
    PeakParameterFwhm() : minimalPeakHeightToLearnFrom_(0) {}

    /**
     * The FWHM at a specific mass channel.
     *
     * @param mz Mass channel; has to be positive.
     * @throw ms::PreconditionViolation Parameter mz is not positive.
     * @throw ms::PostconditionViolation The computed fwhm is negative or zero. This may
     *      be caused by an invalid ParameterModel.
     */
    double at(const double mz) const {
        mspp_precondition(mz > 0, "PeakParameterFwhm::at(): Parameter mz has to be positive.");
        double fwhm = this->ParameterModel::at(mz);
        mspp_postcondition(fwhm > 0, "PeakParameterFwhm::at(): Model returned negative or zero fwhm.");
        return fwhm;
    }

    /**
     * Calibrates the internal model for a specific mass spectrum.
     *
     * Use this function to learn from a sequence of SparseSpectrum::Elements.
     * To learn from a whole SparseSpectrum just set first to spectrum.begin() and last to
     * spectrum.end().
     *
     * Note, that there is no internal error threshold or similar for the quality of the
     * calibration. The calibration is performed as long as it is possible in any way.
     * To achieve a good result, one should filter out the noise of the input spectrum and/
     * or use high quality data in the first place.
     *
     * The distance (last - first) may not be negative, else the behaviour is undefined.
     * It doesn't violate the preconditions, if (last - first) is zero or small. Nevertheless,
     * it increases the chance of a Starvation exception to happen.
     *
     * @param first Points to the first Element of the sequence.
     * @param last Points to one past the last Element of the sequence.
     *
     * @throw ms::Starvation To few or bad data extracted from the input sequence to make a
     *                       calibration possible.
     */
    void learnFrom(SparseSpectrum::const_iterator first, SparseSpectrum::const_iterator last);

    // setMinimalPeakHeightToLearnFrom()
    /**
     * Only use peaks with a minimum absolute intensity to learn from.
     *
     * @param minimalHeight May be negative, albeit that's not meaningful.
     */
    void setMinimalPeakHeightToLearnFrom(double minimalHeight);
    
    // getMinimalPeakHeightToLearnFrom()
    /**
     * Only peaks with this absolute minimal intensity are used for learning.
     */
    double getMinimalPeakHeightToLearnFrom();

private:
    static const double fractionOfMaximum_;

    double minimalPeakHeightToLearnFrom_; 
    
    typedef std::pair<SparseSpectrum::Element::first_type, SparseSpectrum::Element::first_type> MzWidthPair_;
    typedef std::vector<MzWidthPair_> MzWidthPairs_;

    /**
     * Fit the parameter model to measured mz-width pairs.
     *
     * Internally, a linear regression using least squares is done. The PeakParameterModel
     * has to support the optional slopeInParameterSpaceFor() function for this regression
     * to work. Else, calling this function will not compile.
     *
     * If new regression methods should be added, move the implementation of the learn_
     * function to policy classes.
     *
     * @throw ms::PreconditionViolation Parameter pairs is an empty vector.
     * @throw ms::InvariantViolation Numerical regression algorithm failed.
     */
    void learn_(const MzWidthPairs_& pairs);
};

/**
 * Fwhm as it occurs in an Orbitrap mass spectrum.
 */
typedef PeakParameterFwhm<LinearSqrtModel> OrbitrapFwhm;
/**
 * Fwhm as it occurs in an Orbitrap mass spectrum. Is zero at zero Dalton.
 */
typedef PeakParameterFwhm<LinearSqrtOriginModel> OrbitrapWithOriginFwhm;

/**
 * Fwhm as it occurs in a FT-ICR mass spectrum.
 */
typedef PeakParameterFwhm<QuadraticModel> FtIcrFwhm;

/**
 * Fwhm as it occurs in a TOF mass spectrum.
 *
 * The specific Time-of-Flight mass analyzer should measure time internally (not velocity
 * or energy), for this peak parameter to be applicable.
 */
typedef PeakParameterFwhm<SqrtModel> TofFwhm; 

/**
 * A FWHM independent of the mass channel.
 */
typedef PeakParameterFwhm<ConstantModel> ConstantFwhm;



////////////////////
/* implementation */
////////////////////

template <typename ParameterModel>
const double PeakParameterFwhm<ParameterModel>::fractionOfMaximum_ = 0.5;  

// learnFrom()
template <typename ParameterModel>
void PeakParameterFwhm<ParameterModel>::learnFrom(SparseSpectrum::const_iterator first, SparseSpectrum::const_iterator last) {
    using namespace std;
 
    // sample some FWHMs from the spectrum
    MzWidthPairs_ pairs = measureFullWidths(first, last, fractionOfMaximum_, getMinimalPeakHeightToLearnFrom());        

    if(pairs.empty()) {
        throw ms::Starvation("PeakParameterFwhm::learnFrom(): No (Mz | FWHM) could be measured in input spectrum to learn from.");
    }

    // fit the PeakParameterModel to the measured data
    try {
        learn_(pairs);
    } catch(const ms::InvariantViolation& e) {
		MSPP_UNUSED(e);
        MSPP_LOG(logWARNING) << "PeakParameterFwhm::learnFrom(): Numerical regression failed.";
        throw ms::Starvation("PeakParameterFwhm::learnFrom(): Regression of the parameter model for the measured (Mz | FWHM) pairs failed.");
    }
    
    MSPP_LOG(logINFO) << "Learned peak parameter FWHM from spectrum. FWHM at 400 Th is now " << at(400)  << " Th. This corresponds to a resolution of " << 400./at(400) << ".";
}

template <typename ParameterModel>
void PeakParameterFwhm<ParameterModel>::setMinimalPeakHeightToLearnFrom(const double minimalHeight) {
    minimalPeakHeightToLearnFrom_ = minimalHeight;
}

template <typename ParameterModel>
double PeakParameterFwhm<ParameterModel>::getMinimalPeakHeightToLearnFrom() {
    return minimalPeakHeightToLearnFrom_;
}

// learn_()
template <typename ParameterModel>
void PeakParameterFwhm<ParameterModel>::learn_(const MzWidthPairs_& pairs) {
    using namespace vigra;
    mspp_precondition(pairs.empty() == false, "PeakParameterFwhm::learn_(): Called with empty input vector. This is not supposed to happen. A bug in the code preceding the call of learn_ probably caused it.");

    // We now fit the PeakParameterModel to the spectrum using linear regression.    
    // We minimize the residue |A*x - b|^2.
    // b is a column vector with all measured widths.
    // x is a column vector with the model parameters, for example (a, b)^T.
    // A is a matrix such that A*x resembles the original model in every row for the
    // corresponding mz value.
    // x gets optimized.

    /* Construct A and b */

    // Models with no parameter should not implement the numberOfParameters() function.
    // Calling it here while trying to learn this non-existing parameters would cause a
    // desirable compile time error. Nevertheless, we have to guard against wrongly implemented
    // ParameterModels.
    mspp_invariant(this->ParameterModel::numberOfParameters() > 0, "PeakParameterFwhm::learn_(): Number of model parameters is not greater than zero.");
    
    // A: #rows is number of measured pairs; #columns is dimension of parameter space
    linalg::Matrix<double> A(pairs.size(), this->ParameterModel::numberOfParameters());
    // b: column vector with as many elements as measured pairs
    linalg::Matrix<double> b(pairs.size(), 1);
    GeneralizedSlope slope;

    // Calc GeneralizedSlope for every measured pair and store it as rows of A
    // Store the measured width in b
    for(MzWidthPairs_::size_type pairIndex = 0; pairIndex < pairs.size(); ++pairIndex) {
        slope = this->ParameterModel::slopeInParameterSpaceFor(pairs[pairIndex].first);

        // copy the slope into a row of A
        // (we ignore the bias, because it can't be optimized.)
        mspp_invariant((slope.size()-1) == static_cast<GeneralizedSlope::size_type>(A.columnCount()), "PeakParameterFwhm::learn_(): Generalized slope has different dimension than the space, it is living in.");        
        for(linalg::Matrix<double>::difference_type_1 column = 0; column < A.columnCount(); ++column) {      
            A(static_cast<linalg::Matrix<double>::difference_type_1>(pairIndex), column) = slope.at(static_cast<GeneralizedSlope::size_type>(column));   
        }       

        b(static_cast<linalg::Matrix<double>::difference_type_1>(pairIndex), 0) = pairs[pairIndex].second;    
    }

    /* do least squares */
    // result: the optimized parameters
    // Note, that we don't include the bias.
    linalg::Matrix<double> x(this->ParameterModel::numberOfParameters(),1);
    
    // We have to enforce a positive FWHM for positive mz values, so we use a non-negative
    // least squares with x guaranteed to be non-negative. The model then has to yield
    // positive values too, of course. But that is in the responsibility of the caller. 
    linalg::nonnegativeLeastSquares(A, b, x);

    /* set the fitted parameters */
    for(linalg::Matrix<double>::difference_type_1 index = 0; index < x.rowCount(); ++index) {
        MSPP_LOG(logDEBUG2) << "PeakParameterFwhm::learn_(): Parameter " << index << " found: " <<  x(index, 0);
        
        mspp_invariant(index >= 0, "PeakParameterFwhm::learn_(): index may not be negative before calling setParameter(index,0).");
        this->ParameterModel::setParameter(index, x(static_cast<unsigned>(index), 0));    
    }
}

} /* namespace ms */

#endif /*__PEAKPARAMETER_H__*/
