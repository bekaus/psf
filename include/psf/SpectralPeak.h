/* $Id: SpectralPeak.h 2616 2009-10-28 08:54:01Z bkausler $ */

/*
 * SpectralPeak.h
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

#ifndef __SPECTRALPEAK_H__
#define __SPECTRALPEAK_H__

#include <algorithm>
#include <iterator>

#include <ms++/config.h>
#include <ms++/SparseSpectrum.h>
#include <psf/SpectrumAlgorithm.h>
#include <psf/Predicates.h>

namespace ms
{

/**
 * Aggregates algorithms working on a spectral peak.
 *
 * A spectral peak is a single peak in a mass spectrum; in contrary to a monoisotopic peak,
 * which represents a whole isotope pattern.
 *
 * A spectral peak is represented as a sequence of ms::SparseSpectrum::Element. The
 * algorithms work with iterators on this sequence. The Elements have to be in ascending
 * order of their m/z values.
 * There are no further requirements to constitute a sequence a spectral peak. For 
 * example, even a set of equiabundant elements can be seen as a spectral peak.
 *
 * @see ms::Peak
 */
namespace SpectralPeak
{
    // SpectralPeak::height()
    /**
     * The height of a spectral peak.
     *
     * The heighest abundance in a sequence of spectral elements is detected.
     * This abundance is interpreted as the peak height.
     *
     * The minimal requirements of the sequence are:
     * @li At least one element.
     *
     * @param firstElement Points to the first element of the sequences of spectral elements
     *                     to be considered a peak.
     * @param lastElement Points to the last element of the sequences of spectral elements
     *                    to be considered a peak.
     * @return The peak height, in units of the element abundance.
     *
     * @throw ms::PreconditionViolation Minimal sequence requirements aren't met.
     *
     * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
     */
    template< typename FwdIter, typename IntensityExtractor >
    MSPP_EXPORT
    typename IntensityExtractor::result_type
    height(const IntensityExtractor&, FwdIter firstElement, FwdIter lastElement);     

    // SpectralPeak::lowness()
    /**
      * The lowness of a spectral peak.
      *
      * The highest element in a sequence of spectral elements is detected and the two
      * elements with the lowest abundance are searched; one on the left and one on the
      * right of the maximum (the maximum itself can be selected). The more abundant peak of these two is chosen.
      * One minus the ratio of this abundance to the maximum is called  the 'peak lowness'.
      *
      * An equiabundant sequence of spectral elements has a lowness of 0.0. In contrast,
      * a maximum flanked by two elements with almost zero abundance has a lowness of
      * almost 1.00.
      *
      * The minimal requirements of the sequence are:
      * @li At least one element. The lowness is then 0.0.
      *       
      * @param firstElement Points to the first element of the sequences of spectral elements
      *                     to be considered a peak.
      * @param lastElement Points to the last element of the sequences of spectral elements
      *                    to be considered a peak.
      * @return The lowness is between 0.0 and 1.0, borders included.
      *
      * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
      */
    template< typename FwdIter, typename IntensityExtractor >
    MSPP_EXPORT 
    double lowness(const IntensityExtractor&, FwdIter firstElement, FwdIter lastElement);

    // SpectralPeak::fullWidthAtFractionOfMaximum()
    /**
      * The full width at a fraction of the maximum of a spectral peak.
      *
      * The most abundant element in a sequence of spectral elements is found. Coming from 
      * the left and right, the two elements which are nearest to the fraction of the
      * abundance maximum are search: abundance fraction <= element abundance. The two
      * elements on both sides are then linearly interpolated with their neighboring elements
      * just below the fraction of the maximum. (In case of fraction == element abundance, this is avoided.)
      * The distance in m/z dimension between the two interpolated elements
      * flanking the maximum is then returned as the full width at the fraction of the maximum.
      *
      * For example, setting the fraction to 0.5 just returns the full width at half maximum.
      *
      * The minimal requirements of the sequence are:
      * @li At least one abundance maximum. If there are multiple maxima, the first such element
      *     is chosen.
      * @li At least one element below the fraction of the maximum on both flanks of the maximum.
      *     The maximum itself may be chosen as the element above the fraction.
      *
      * @param firstElement Points to the first element of the sequences of spectral elements
      *                     to be considered a peak.
      * @param lastElement Points to the last elmement of the sequences of spectral elements
      *                    to be considered a peak.
      * @param fraction Has to be between 0.0 and 1.0, borders included.
      * @return The full width at the fraction of the maximum, in units of m/z.
      *
      * @throw ms::PreconditionViolation Parameter fraction not in the required range.
      * @throw ms::Starvation The sequence doesn't satisfy the minimal requirements.
      *
      * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
      */
    template< typename FwdIter, typename MzExtractor, typename IntensityExtractor >
    MSPP_EXPORT
    typename MzExtractor::result_type 
    fullWidthAtFractionOfMaximum(const MzExtractor&, const IntensityExtractor&, FwdIter firstElement, FwdIter lastElement, const double fraction);

} /* namespace SpectralPeak */



/* implementation */

// height()
template< typename FwdIter, typename IntensityExtractor >
typename IntensityExtractor::result_type
SpectralPeak::height(const IntensityExtractor& get_int, FwdIter firstElement, FwdIter lastElement) {
    mspp_precondition(distance(firstElement, lastElement) >=0, "SpectralPeak::height(): Distance between first and last input element is not nonnegative.");
    
    // Compare elements by intensity
    LessByExtractor<typename IntensityExtractor::element_type, IntensityExtractor> comp(get_int);
    // find maximum intensity
    FwdIter maximum = max_element(firstElement, ++lastElement, comp);

    return get_int(*maximum);
}



// lowness()
template< typename FwdIter, typename IntensityExtractor >
double ms::SpectralPeak::lowness(const IntensityExtractor& get_int, FwdIter firstElement, FwdIter lastElement) {
    // comply to STL standards.
    FwdIter last = lastElement;
    ++last;    

    // Compare elements by intensity
    LessByExtractor<typename IntensityExtractor::element_type, IntensityExtractor> comp(get_int);

    // find maximum intensity
    FwdIter maximum = max_element(firstElement, last, comp);

    // find least abundant element right of the maximum
    FwdIter rightMinimum = min_element(maximum, last, comp);
    // and to the left (both times with the maximum included as possible minimum)
    FwdIter leftMinimum = min_element(firstElement, ++maximum, comp);
    --maximum; // STL required [first, last)

    // more abundant element of the two
    typename IntensityExtractor::element_type moreAbundantOne = std::max(*leftMinimum, *rightMinimum, comp);

    return 1. - (get_int(moreAbundantOne)/get_int(*maximum));
}



// fullWidthAtFractionOfMaximum(): private implementation details
namespace 
{
// findElementBelowTargetAbundance()
/**
 * Find Element below the target abundance given the Element above or on the target abundance.
 *
 * We are working on the sequence [firstElement, above]. Every Element besides the 'above'
 * element has to be lower than the 'target' abundance. The above Element may actually
 * be exactly on or above the target abundance. If these preconditions are not met,
 * the behaviour is undefined.
 *
 * @param firstElement Input const_iterator to the first SparseSpectrum::Element of the sequence to be searched.
 * @param above Input const_iterator to the SparseSpectrum::Element on or above the target abundance. The distance 
 *              above - firstElement may not be negative. Else, the behaviour is undefined.
 * @param target The target abundance
 *
 * @return The SparseSpectrum::Element below the target abundance nearest (in mz dimension) to the above Element.
 *         Is the above Element exactly on the target, below is the same as the above Element.
 *
 * @throw ms::Starvation No Element below found.
 */      
template <typename const_InIter, typename IntensityExtractor>
const_InIter findElementBelowTargetAbundance(const IntensityExtractor&, const const_InIter firstElement, const const_InIter above, const typename IntensityExtractor::result_type target);

// interpolateElements()
/**
 * Takes two elements and blends them together with a specific target intensity.
 *
 * The interpolation is done linearly in the mz dimension so that the result has the 
 * target abundance. The order of element1 and element2 is not important. If element1 and element2 have the same mz value,
 * then the mz value is directly return without interpolation.
 * If element1 and element2 differ in mz, they have to differ in abundance, too.
 *
 * @return The mz value of the interpolated element with an abundance of 'target'.
 *
 * @throw ms::InvariantViolation Element1 and element2 differ in mz but not in abundance.
 *                               No interpolation to a general target is then possible.
 */
template< typename MzExtractor, typename IntensityExtractor, typename element_type>
typename MzExtractor::result_type interpolateElements(const MzExtractor&, const IntensityExtractor&, const element_type& element1, const element_type& element2, const typename IntensityExtractor::result_type target);

} /* anonymous namespace */

// fullWidthAtFractionOfMaximum()
template< typename FwdIter, typename MzExtractor, typename IntensityExtractor >
typename MzExtractor::result_type 
SpectralPeak::
fullWidthAtFractionOfMaximum(const MzExtractor& get_mz, const IntensityExtractor& get_int, FwdIter firstElement, FwdIter lastElement, const double fraction) {
    mspp_precondition(0. <= fraction && fraction <= 1., "fullWidthAtFractionOfMaximum(): Fraction parameter out of range.");

    /* Determine target intensity */
    // Compare elements by intensity
    LessByExtractor<typename IntensityExtractor::element_type, IntensityExtractor> comp(get_int);
    // find maximum intensity
    FwdIter maximum = max_element(firstElement, lastElement + 1, comp);
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): Spectral peak maximum detected at (mz, intensity): " << get_mz(*maximum) << " ," << get_int(*maximum); 
    // calc target intensity
    const typename IntensityExtractor::result_type target = get_int(*maximum) * fraction;
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): Fraction of maximal intensity is: " << target;
    
    // we need that further below for finding elements
    MoreThanValue<typename IntensityExtractor::element_type, IntensityExtractor> compScalar(get_int, target);

    /* find utter left element nearest above or on target */
    // target <= above == !(above < target) 
    FwdIter aboveOnLeft = find_if(firstElement, maximum + 1, compScalar);
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): aboveOnLeft detected at (mz, intensity): " << get_mz(*aboveOnLeft) << " ," << get_int(*aboveOnLeft); 
    // determine belowOnLeft
    FwdIter belowOnLeft = findElementBelowTargetAbundance(get_int, firstElement, aboveOnLeft, target);
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): belowOnLeft detected at (mz, intensity): " << get_mz(*belowOnLeft) << " ," << get_int(*belowOnLeft);

    /* find utter right element */
    // we now start searching from the right
    std::reverse_iterator<FwdIter> rlast(lastElement + 1);
    std::reverse_iterator<FwdIter> rmaximum(maximum);
    // target <= above == !(above < target) 
    std::reverse_iterator<FwdIter> aboveOnRight = find_if(rlast, rmaximum, compScalar);
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): aboveOnRight detected at (mz, intensity): " << get_mz(*aboveOnRight) << " ," << get_int(*aboveOnRight);     
    // determine belowOnRight
    std::reverse_iterator<FwdIter> belowOnRight = findElementBelowTargetAbundance(get_int, rlast, aboveOnRight, target);
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): belowOnRight detected at (mz, intensity): " << get_mz(*belowOnRight) << " ," << get_int(*belowOnRight); 
    
    /* interpolate below and above elements */
    typename MzExtractor::result_type leftInterpolated = interpolateElements(get_mz, get_int, *belowOnLeft, *aboveOnLeft, target);
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): leftInterpolated is: " << leftInterpolated;
    typename MzExtractor::result_type rightInterpolated = interpolateElements(get_mz, get_int, *belowOnRight, *aboveOnRight, target);
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): rightInterpolated is: " << rightInterpolated;

    return rightInterpolated - leftInterpolated;
}



namespace
{
    // findElementBelowTargetAbundance() 
    template <typename const_InIter, typename IntensityExtractor>
    const_InIter findElementBelowTargetAbundance(const IntensityExtractor& get_int, const const_InIter firstElement, const const_InIter above, const typename IntensityExtractor::result_type target) {
        const_InIter below;

        // Check if the above Element coincides with the firstElement
        if (firstElement == above) {
            // Rule out special case, where abundance of above Element equals target abundance
            if(target < get_int(*above)) { 
                throw Starvation("fullWidthAtFractionOfMaximum(): No elements on the left below target abundance.");
            }
            // special case: target == above
            else {
                below = above;
                MSPP_LOG(logDEBUG2) << "findElementBelowTargetAbundance(): Target abundance equals abundance of element above. Setting below element equal to above element.";
            }        
        } else {
            // Choose nearest neighbor in mz dimension (guaranteed to be below target abundance by our imposed preconditions).
            below = above - 1;
        }

        return below;
    }

    // interpolateElements()
    template< typename MzExtractor, typename IntensityExtractor, typename element_type>
    typename MzExtractor::result_type interpolateElements(const MzExtractor& get_mz, const IntensityExtractor& get_int, const element_type& element1, const element_type& element2, const typename IntensityExtractor::result_type target) {
        if (get_mz(element1) == get_mz(element2)) {
            return get_mz(element2);
        } else {
            mspp_invariant(get_int(element1) != get_int(element2), "interpolateElements(): Illegal abundance state: below < target && target <= above && above == below."); 

            // abundance = slope * mz + shift      
            double slope = 1.0 * (get_int(element2) - get_int(element1)) / (get_mz(element2) - get_mz(element1));
            MSPP_LOG(logDEBUG2) << "interpolateElements(): slope of linear interpolation: " << slope;              

            // just take one of the two Elements to determine the shift
            double shift = get_int(element1) - slope * get_mz(element1);
            MSPP_LOG(logDEBUG2) << "interpolateElements(): shift of linear interpolation: " << shift; 

            // => mz = (abundance - shift)/slope
            return static_cast<typename MzExtractor::result_type>( (target - shift) / slope );
        }
    }
} /* anonymous namespace */

} /* namespace ms */

#endif /*__SPECTRALPEAK_H__*/
