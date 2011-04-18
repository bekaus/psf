/* $Id: PeakShape.h 2613 2009-10-26 23:34:00Z bkausler $ */

/*
 * PeakShape.h
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

#ifndef __PEAKSHAPE_H__
#define __PEAKSHAPE_H__

#include <ms++/config.h>

// for friend declaration further below
struct peakshapeTestSuite;

namespace ms
{

// class PeakShape
/**
 * The shape of a spectral peak.
 *
 * @attention
 * This abstract class exists only for documentation reasons. You shouldn't inherit from
 * it, but just
 * implement the interface in your own peak shapes. That is, because the peak shapes are
 * used as template parameters and directly compiled into the classes, which are
 * using the peak shapes.
 * Therefore, there is no need for a virtual inheritance tree at runtime.
 *
 * In a perfect world ions of the same mass to charge ratio would appear as a sharp stick
 * in a mass spectrum with its intensity proportional to the number of ions.
 * Unfortunately this stick gets blurred out under physical conditions due to an imprecise
 * measurement process. The concrete shape of this blurred peak depends on the type of mass
 * spectrometer used. 
 *
 * Many different theoretical descriptions for peak shapes exist. PeakShape provides an
 * abstract interface for implementing these theories.
 *
 * We assume no special normalization of the peak shape's area. This can speed up calculations.
 *
 * @see ms::GaussianPeakShape
 * 
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * @date 2009-07-09
 */
class PeakShape
{
public:
    // base class needs virtual destructor
    virtual ~PeakShape() {}
    // at()
    /**
     * Returns the height of the peak shape at a x coordinate.
     *
     * The position of the true mass is at x coordinate zero.
     *
     * The absolute value of the peak height is arbitrary. Ony the height relative to
     * other x coordinates is important, i.e. the peak shape is not normalized.
     */
    virtual double at(const double xCoordinate) const = 0;

    // getSupportThreshold()
    /**
     * Return the peak shape support.
     *
     * This threshold is a positive distance measured from the true mass at x coordinate zero
     * and is symmetrical around the center. After the threshold, the height of the peak
     * shape is assumed to be so low, that it could be set to zero for every practical purposes.
     *
     * In the case of an asymmetrical peakshape, the bigger value is chosen for the
     * threshold.
     *
     * You may use this information to speed up calculations depending on the peak shape at
     * at specific coordinate.
     */
    virtual double getSupportThreshold() const = 0;
};



// class BoxPeakShape
/**
 * A Gaussian-based box peak shape.
 *
 * The idea is that centroided data can more efficiently and probably
 * also more accurately be fit using a box peak shape. This is 
 * because the integration over the true signal PSF is carried out
 * by the instrument, and depending on the vendor's propietray 
 * centroiding algorithms, the true m_0 may shift slightly back and
 * forth. In this case, as the PSF intensities have already been 
 * integrated into the centroid intensity, it is not suitable to
 * weight the observation by its devaition from the expected zero:
 * we either take all or nothig.
 * Nonetheless, the width of the function needs to be adapted properly
 * and that is why we rely on the Gaussian that would be used
 * for profile data. 
 *
 * @see ms::PeakShape
 * 
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * @author Marc Kirchner <marc.kirchner@childrens.harvard.edu>
 * @date 2009-10-02
 */
class MSPP_EXPORT BoxPeakShape
{
public:
    double at(const double xCoordinate) const;

    // getSupportThreshold()
    /**
     * The support threshold for the box is calculated based on
     * a Gaussian according to 'sigma x sigmaFactorForSupportThreshold'.
     */
    double getSupportThreshold() const;

public:
    explicit BoxPeakShape(const double sigma = 0.1, const double sigmaFactorForSupportThreshold = 3.0);
    
    // setSigma()
    /**
     * Sets the sigma parameter of the underlying Gaussian.
     *
     * @param sigma Has to be positive.
     * @throw ms::PreconditionViolation Sigma is not positive.
     */
    void setSigma(const double sigma);

    // getSigma()
    /**
     * Gets the sigma parameter of the Gaussian the underlies the box.
     */
    double getSigma() const { return sigma_; }

    // setFwhm()
    /**
     * Sets the full width at half maximum.
     *
     * This changes the sigma parameter according to @f$ \mathrm{FWHM}=2\sqrt{2\ln2}\cdot\sigma @f$.
     *
     * @param fwhm Has to be positive.
     * @throw ms::PreconditionViolation fwhm is not positive.
     */
    void setFwhm(const double fwhm);
    
    // getFwhm()
    /**
     * Gets the full width at half maximum.
     */
    double getFwhm() const;

    // setSigmaFactorForSupportThreshold()
    /**
     * Sets the factor for the threshold calculation.
     *
     * @param factor Has to be positive.
     * @throw ms::PreconditionViolation factor is not positive.
     */
    void setSigmaFactorForSupportThreshold(const double factor);

    // getSigmaFactorForSupportThreshold()
    /**
     * Returns the factor used in the support threshold calculation.
     */
    double getSigmaFactorForSupportThreshold() const;

protected:
    // fwhmToSigmaConversionFactor()
    /**
     * @f$ 2\sqrt{2\ln2} @f$
     */
    double sigmaToFwhmConversionFactor() const;

private:
    double sigma_;
    double sigmaFactorForSupportThreshold_;

friend struct ::peakshapeTestSuite;
};

// class GaussianPeakShape
/**
 * A gaussian peak shape.
 *
 * The gaussian is @f$ e^{-\frac{x^2}{2\cdot\sigma^{2}}} @f$.
 *
 * The gaussian is centered around zero.
 *
 * @see ms::PeakShape
 * 
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * @date 2009-07-09
 */
class MSPP_EXPORT GaussianPeakShape
{
public:
    double at(const double xCoordinate) const;

    // getSupportThreshold()
    /**
     * The support threshold for the gaussian is calculated according to
     * 'sigma x sigmaFactorForSupportThreshold'.
     */
    double getSupportThreshold() const;

public:
    explicit GaussianPeakShape(const double sigma = 0.1, const double sigmaFactorForSupportThreshold = 3.0);
    
    // setSigma()
    /**
     * Sets the sigma parameter of the gaussian.
     *
     * @param sigma Has to be positive.
     * @throw ms::PreconditionViolation Sigma is not positive.
     */
    void setSigma(const double sigma);

    // getSigma()
    /**
     * Gets the sigma parameter of the gaussian.
     */
    double getSigma() const { return sigma_; }

    // setFwhm()
    /**
     * Sets the full width at half maximum.
     *
     * This changes the sigma parameter according to @f$ \mathrm{FWHM}=2\sqrt{2\ln2}\cdot\sigma @f$.
     *
     * @param fwhm Has to be positive.
     * @throw ms::PreconditionViolation fwhm is not positive.
     */
    void setFwhm(const double fwhm);
    
    // getFwhm()
    /**
     * Gets the full width at half maximum.
     */
    double getFwhm() const;

    // setSigmaFactorForSupportThreshold()
    /**
     * Sets the factor for the threshold calculation.
     *
     * @param factor Has to be positive.
     * @throw ms::PreconditionViolation factor is not positive.
     */
    void setSigmaFactorForSupportThreshold(const double factor);

    // getSigmaFactorForSupportThreshold()
    /**
     * Returns the factor used in the support threshold calculation.
     */
    double getSigmaFactorForSupportThreshold() const;

protected:
    // fwhmToSigmaConversionFactor()
    /**
     * @f$ 2\sqrt{2\ln2} @f$
     */
    double sigmaToFwhmConversionFactor() const;

private:
    double sigma_;
    double sigmaFactorForSupportThreshold_;

friend struct ::peakshapeTestSuite;
};



// class LorentzianPeakShape
/**
 * A lorentzian peak shape.
 *
 * The Lorentzian is @f$ \frac{\mathrm{fwhm}}{x^2+\mathrm{fwhm}^2} @f$.
 *
 * The Lorentzian is centered around zero.
 *
 * @see ms::PeakShape
 * 
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class MSPP_EXPORT LorentzianPeakShape
{
public:
    double at(const double xCoordinate) const;

    // getSupportThreshold()
    /**
     * The support threshold for the gaussian is calculated according to
     * 'sigma x sigmaFactorForSupportThreshold'.
     */
    double getSupportThreshold() const;

public:
    explicit LorentzianPeakShape(double fwhm = 0.1, const double fwhmFactorForSupportThreshold = 5.0);
    
    // setFwhm()
    /**
     * Sets the full width at half maximum.
     *
     * @param fwhm Has to be positive.
     * @throw ms::PreconditionViolation fwhm is not positive.
     */
    void setFwhm(double fwhm);
    
    // getFwhm()
    /**
     * Gets the full width at half maximum.
     */
    double getFwhm() const;

    // setFwhmFactorForSupportThreshold()
    /**
     * Sets the factor for the threshold calculation.
     *
     * @param factor Has to be positive.
     * @throw ms::PreconditionViolation factor is not positive.
     */
    void setFwhmFactorForSupportThreshold(const double factor);

    // getFwhmFactorForSupportThreshold()
    /**
     * Returns the factor used in the support threshold calculation.
     */
    double getFwhmFactorForSupportThreshold() const;

private:
    double fwhm_;
    double fwhmFactorForSupportThreshold_;

friend struct ::peakshapeTestSuite;
};

} /* namespace ms */

#endif /*__PEAKSHAPE_H__*/

