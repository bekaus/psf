/* $Id: SpectralPeak.cpp 2418 2009-09-21 16:02:56Z bkausler $ */

/*
 * SpectralPeak.cpp
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

#include <algorithm>
#include <functional>
#include <iterator>
#include <vector>

#include <ms++/Error.h>
#include <ms++/Log.h>
#include <ms++/SpectralPeak.h>

using namespace std;
using namespace ms;

// height()
SparseSpectrum::Element::second_type
ms::SpectralPeak::height(SparseSpectrum::const_iterator firstElement, SparseSpectrum::const_iterator lastElement) {
    mspp_precondition(distance(firstElement, lastElement) >=0, "SpectralPeak::height(): Distance between first and last input element is not nonnegative.");
    
    // Compare elements by abundance
    SparseSpectrum::LessThanAbundance<SparseSpectrum::Element, SparseSpectrum::Element> comp; 

    // find maximum abundance
    SparseSpectrum::const_iterator maximum = max_element(firstElement, ++lastElement, comp);

    return maximum->abundance;
}



// lowness()
double ms::SpectralPeak::lowness(SparseSpectrum::const_iterator firstElement, SparseSpectrum::const_iterator lastElement) {
    // comply to STL standards.
    SparseSpectrum::const_iterator last = lastElement;
    ++last;    

    // Compare elements by abundance
    SparseSpectrum::LessThanAbundance<SparseSpectrum::Element, SparseSpectrum::Element> comp; 

    // find maximum abundance
    SparseSpectrum::const_iterator maximum = max_element(firstElement, last, comp);

    // find least abundant element right of the maximum
    SparseSpectrum::const_iterator rightMinimum = min_element(maximum, last, comp);
    // and to the left (both times with the maximum included as possible minimum)
    SparseSpectrum::const_iterator leftMinimum = min_element(firstElement, ++maximum, comp);
    --maximum; // STL required [first, last)

    // more abundant element of the two
    SparseSpectrum::Element moreAbundantOne = max(*leftMinimum, *rightMinimum, comp);

    return 1. - (moreAbundantOne.abundance/maximum->abundance);
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
template <typename const_InIter>
const_InIter findElementBelowTargetAbundance(const const_InIter firstElement, const const_InIter above, const ms::SparseSpectrum::Element::second_type target);

// interpolateElements()
/**
 * Takes two elements and blends them together with a specific target abundance.
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
SparseSpectrum::Element::first_type interpolateElements(const SparseSpectrum::Element& element1, const SparseSpectrum::Element& element2, const ms::SparseSpectrum::Element::second_type target);

} /* anonymous namespace */
     
// fullWidthAtFractionOfMaximum()
ms::SparseSpectrum::Element::first_type 
ms::SpectralPeak::
fullWidthAtFractionOfMaximum(SparseSpectrum::const_iterator firstElement, SparseSpectrum::const_iterator lastElement, const double fraction) {
    mspp_precondition(0. <= fraction && fraction <= 1., "fullWidthAtFractionOfMaximum(): Fraction parameter out of range.");

    /* Determine target abundance */
    // Compare elements by abundance
    SparseSpectrum::LessThanAbundance<SparseSpectrum::Element, SparseSpectrum::Element> comp;
    // find maximum abundance
    SparseSpectrum::const_iterator maximum = max_element(firstElement, lastElement + 1, comp);
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): Spectral peak maximum detected at (mz, abundance): " << maximum->mz << " ," << maximum->abundance; 
    // calc target abundance
    ms::SparseSpectrum::Element::second_type target = maximum->abundance * fraction;
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): Fraction of maximal abundance is: " << target;
    
    // we need that further below for finding elements
    SparseSpectrum::LessThanAbundance<SparseSpectrum::Element,double> compScalar;

    /* find utter left element nearest above or on target */
    // target <= above == !(above < target) 
    SparseSpectrum::const_iterator aboveOnLeft = find_if(firstElement, maximum + 1, not1(bind2nd(compScalar, target)));
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): aboveOnLeft detected at (mz, abundance): " << aboveOnLeft->mz << " ," << aboveOnLeft->abundance; 
    // determine belowOnLeft
    SparseSpectrum::const_iterator belowOnLeft = findElementBelowTargetAbundance(firstElement, aboveOnLeft, target);
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): belowOnLeft detected at (mz, abundance): " << belowOnLeft->mz << " ," << belowOnLeft->abundance;

    /* find utter right element */
    // we now start searching from the right
    reverse_iterator<SparseSpectrum::const_iterator> rlast(lastElement + 1);
    reverse_iterator<SparseSpectrum::const_iterator> rmaximum(maximum);
    // target <= above == !(above < target) 
    reverse_iterator<SparseSpectrum::const_iterator> aboveOnRight = find_if(rlast, rmaximum, not1(bind2nd(compScalar, target)));
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): aboveOnRight detected at (mz, abundance): " << aboveOnRight->mz << " ," << aboveOnRight->abundance;     
    // determine belowOnRight
    reverse_iterator<SparseSpectrum::const_iterator> belowOnRight = findElementBelowTargetAbundance(rlast, aboveOnRight, target);
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): belowOnRight detected at (mz, abundance): " << belowOnRight->mz << " ," << belowOnRight->abundance; 
    
    /* interpolate below and above elements */
    SparseSpectrum::Element::first_type leftInterpolated = interpolateElements(*belowOnLeft, *aboveOnLeft, target);
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): leftInterpolated is: " << leftInterpolated;
    SparseSpectrum::Element::first_type rightInterpolated = interpolateElements(*belowOnRight, *aboveOnRight, target);
    MSPP_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): rightInterpolated is: " << rightInterpolated;

    return rightInterpolated - leftInterpolated;
}



/* private implementation details */

namespace
{
    // findElementBelowTargetAbundance() 
    template <typename const_InIter>
    const_InIter findElementBelowTargetAbundance(const const_InIter firstElement, const const_InIter above, const ms::SparseSpectrum::Element::second_type target) {
        const_InIter below;

        // Check if the above Element coincides with the firstElement
        if (firstElement == above) {
            // Rule out special case, where abundance of above Element equals target abundance
            if(target < above->abundance) { 
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
    SparseSpectrum::Element::first_type interpolateElements(const SparseSpectrum::Element& element1, const SparseSpectrum::Element& element2, const ms::SparseSpectrum::Element::second_type target) {
        if (element1.mz == element2.mz) {
            return element2.mz;
        } else {
            mspp_invariant(element1.abundance != element2.abundance, "interpolateElements(): Illegal abundance state: below < target && target <= above && above == below."); 

            // abundance = slope * mz + shift      
            double slope = (element2.abundance - element1.abundance) / (element2.mz - element1.mz);
            MSPP_LOG(logDEBUG2) << "interpolateElements(): slope of linear interpolation: " << slope;              

            // just take one of the two Elements to determine the shift
            double shift = element1.abundance - slope * element1.mz;
            MSPP_LOG(logDEBUG2) << "interpolateElements(): shift of linear interpolation: " << shift; 

            // => mz = (abundance - shift)/slope
            return static_cast<SparseSpectrum::Element::first_type>( (target - shift) / slope );
        }
    }
} /* anonymous namespace */
