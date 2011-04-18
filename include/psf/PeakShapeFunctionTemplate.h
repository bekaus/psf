/* $Id: PeakShapeFunctionTemplate.h 2613 2009-10-26 23:34:00Z bkausler $ */

/*
 * PeakShapeFunctionTemplate.h
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

#ifndef __PEAKSHAPEFUNCTIONTEMPLATE_H__
#define __PEAKSHAPEFUNCTIONTEMPLATE_H__

#include "psf/PeakShapeFunction.h"

namespace ms
{
/**
 * The peak shape functions, which are implementing the abstract PeakShapeFunction interface.
 *
 * The 'box' type is only used for unit testing. The other types are included in the core
 * library.
 *
 * @author Bernhard X. Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * @date 2009-04-14
 */
enum MSPP_EXPORT PeakShapeFunctionTypes {box, gaussian, orbi, orbiBox, tof};

/**
 * Encapsulates the PeakShapeFunctionTypes and provides conversion functions.
 *
 * @author Bernhard X. Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * @date 2009-04-14
 */
class MSPP_EXPORT PeakShapeFunctionType
{
    public:
        /**
         * You may use the constructor in an implict way: 'PeakShapeFunctionType type = gaussian;'.
         */
        PeakShapeFunctionType(PeakShapeFunctionTypes type) : type_(type) {};

        PeakShapeFunctionTypes toEnum();

        /**
         * @return Depending on the type: 'box', 'gaussian', 'orbi', 'orbiBox' or 'time-of-flight'.
         *         If the type is not known, 'unknown'.
         */
        std::string toString();

    protected:
        PeakShapeFunctionTypes type_;
};




/**
 * Implements a generic ms::PeakShapeFunction.
 *
 * This generic implementation of the ms::PeakShapeFunction interface hast to be parameterized
 * by three template parameters.
 *
 * Use 'typedef PeakShapeFunctionTamplate<...> MyPeakShapeFunction' to implement new peak shape
 * functions quickly. Several peak shape functions defined in this way, can be found in the headerfile 'PeakShapeFunctions.h'.
 *
 * At the moment, this template supports up to two general paramters a,b to parameterize the
 * peak shape function. The actual number depends on the PeakParameterT used. If the PeakParameterT doesn't
 * support a or b, you will get a compile time error trying to use them.
 *
 * Furthermore, some PeakParameterT support autocalibration and/or setting a desired value at a specific
 * coordinate in the spectrum (for example the FWHM at 400 Th). In case of no support, your code will
 * not compile after calling such a function. (This is a desirable behaviour.)
 *
 * @param PeakShapeT Describes the general form of a peak. Can be found in the headerfile 'PeakShape.h'
 * @param PeakParameterT The PeakShapeT has to be parameterized according to peak parameters depending
 *                       on the type of mass spectrometer to describe. Can be found in the
 *                       headerfile 'PeakParameter.h'
 * @param PeakShapeFunctionTypeT The proper name of the peak shape function. Can be found
 *                               in the headerfile 'PeakShapeFunction.h'.
 *
 * @see ms::OrbitrapPeakShapeFunction
 * @see ms::GaussianPeakShapeFunction
 */
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
class MSPP_EXPORT PeakShapeFunctionTemplate
{
public:
    PeakShapeFunctionTemplate();
    explicit PeakShapeFunctionTemplate(const double a);
    PeakShapeFunctionTemplate(const double a, const double b);

    /**
     * operator()
     * @param referenceMass the m/z value at the center of the PSF
     * @param observedMass the m/z value of the mass for which the value of the PSF is desired
     * @return the value of the PSF at (observedMass-referenceMass)
     */
    double operator()(const double referenceMass, const double observedMass) const;

    /**
     * Return the width of the PSF support at a specific m/z value.
     *
     * The threshold is a relative distance measured from the center of the peak shape
     * and is symmetrical around the center. After the threshold, the peak shape function
     * is set to zero.
     * @param mz the m/z value of the PSF center
     * @return the width of the PSF support at the given m/z position
     */
    double getSupportThreshold(const double mz) const;

    /**
     * Returns the actual implementation type of the abstract PeakShapeFunction interface.
     *
     * @author Bernhard X. Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
     * @date 2009-04-14
     */
    PeakShapeFunctionType getType();

    void setA(const double a);
    double getA() const;

    void setB(const double b);
    double getB() const;

    /**
     * Autocalibrate peak shape function parameters using regression.
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
     * The elements in the spectrum have to be in ascending order of their mz value. Furthermore,
     * no elements with duplicate mz values are allowed. Else, the behaviour is undefined.
     *
     * @param first Points to the first Element of the spectrum.
     * @param last Points to one past the last Element of the spectrum.
     *
     * @throw ms::Starvation To few or bad data extracted from the input spectrum to make a
     *                       calibration possible.
     */
    void calibrateFor(SparseSpectrum::const_iterator first, SparseSpectrum::const_iterator last);
    
    // setMinimalPeakHeightForCalibration()
    /**
     * Use only peaks with a minimal intensity for autocalibration.
     *
     * The internally used peak parameters have to support this feature. Else, calling the
     * function would result in a compile time error.
     *
     * @param minimalHeight The minimal absolute intensity. Negative values are possible,
     *                      albeit not meaningful.
     */
    void setMinimalPeakHeightForCalibration(double minimalHeight);

    // getMinimalPeakHeightForCalibration()
    /**
     * Uses only peaks with a minimal intensity for autocalibration.
     *
     * @see PeakShapeFunctionTemplate::setMinimalPeakHeightForCalibration
     */
    double getMinimalPeakHeightForCalibration();

private:
    mutable PeakShapeT peakshape_;
    PeakParameterT peakparameter_;
};



////////////////////
/* Implementation */
////////////////////

// operator()
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
double 
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
operator()(const double referenceMass, const double observedMass) const {
    double supportThreshold = this->getSupportThreshold( referenceMass );        
    double massDifference = observedMass - referenceMass;

    if((-supportThreshold <= massDifference) && (massDifference <= supportThreshold)) {
        peakshape_.setFwhm(peakparameter_.at(referenceMass));
        return peakshape_.at(massDifference);
    }
    else {
        return 0.0;
    }
}

// getSupportThreshold()
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
double 
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
getSupportThreshold(const double mz) const {
    peakshape_.setFwhm(peakparameter_.at(mz));
    return peakshape_.getSupportThreshold();        
}

// getType()
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
inline
PeakShapeFunctionType
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
getType() {
    return PeakShapeFunctionTypeT;
}

// PeakShapeFunctionTemplate()
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
PeakShapeFunctionTemplate() {
}

template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
PeakShapeFunctionTemplate(const double a) {
    this->setA(a);
}

template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
PeakShapeFunctionTemplate(const double a, const double b) {
    this->setA(a);
    this->setB(b);
}

// setA()
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
inline
void
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
setA(const double a) {
    peakparameter_.setA(a);
}
// getA()
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
inline
double
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
getA() const {
    return peakparameter_.getA();
}

// setB()
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
inline
void
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
setB(const double b) {
    peakparameter_.setB(b);
}
// getB()
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
inline 
double
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
getB() const {
    return peakparameter_.getB();
}

// calibrateFor()
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
void
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
calibrateFor(SparseSpectrum::const_iterator first, SparseSpectrum::const_iterator last) {
    peakparameter_.learnFrom(first, last);
}

// setMinimalPeakHeightForCalibration()
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
void 
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
setMinimalPeakHeightForCalibration(const double minimalHeight) {
    peakparameter_.setMinimalPeakHeightToLearnFrom(minimalHeight);
}

// getMinimalPeakHeightForCalibration()
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
double
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
getMinimalPeakHeightForCalibration() {
    return peakparameter_.getMinimalPeakHeightToLearnFrom();
}

} /* namespace ms */

#endif /*__PEAKSHAPEFUNCTIONTEMPLATE_H__*/

