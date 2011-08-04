/*
 * SpectrumAlgorithm.h
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

#ifndef __SPECTRUMALGORITHM_H__
#define __SPECTRUMALGORITHM_H__
#include <ms++/config.h>

#include <iterator>
#include <utility>
#include <vector>

#include <ms++/Log.h>
#include <ms++/Error.h>

#include <psf/SpectralPeak.h>
#include <psf/LessByExtractor.h>

namespace ms
{

// findBump()
/**
 * Finds the first 'bump' in a sequence.
 *
 * A bump is a range in a sequence containing a (local) maximum and stricly decreasing values
 * to the left and to the right of the maximum.
 *
 * The smallest possible bump consists of only three elements: .'.
 * 
 * @param first Iterator pointing to the first element of the input sequence. Has to point to
 *              an existing element. Else the behaviour is undefined.
 * @param last Iterator pointing to one past the last element of the input sequence. Else the
 *             behaviour is undefined.
 * @param comp comp(*iter, *(iter+1)) is called to compare two elements. Think of it as a 'less than' operator <.
 * 
 * @return The first bump found in the input sequence: [pair.first, pair.second]. Returns a pair of 'last', if no
 *         bump is found.
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * @date 2009-07-16
 */
template<typename FwdIter, typename Compare>
MSPP_EXPORT std::pair<FwdIter, FwdIter> findBump(FwdIter first, FwdIter last, Compare comp);



// measureFullWidths()
 /**
 * Sample the full width at a fraction of the maximum.
 *
 * Goes through a (sub-)spectrum and measures the full width at a fraction of the maximum for
 * every spectral peak considered pure.
 *
 * A pure peak fulfills the requirements of a 'bump' and is at least as low as the fraction
 * of its maximum. The true peak maximum is estimated as the most abundant element of the 
 * bump.
 *
 * Note: You should see this really as a measurement in the physical sense, meaning even
 * in the case of an exactly calculatable width, this function may return a slightly 
 * different value due to rounding errors and similar effects.
 *
 * The distance (last - first) may not be negative, else the behaviour is undefined.
 * If (last - first) is zero, an empty vector is returned. 
 *
 * @param first The first Element of the sequence to sample from.
 * @param last One past the last Element of the sequence to sample from.
 * @param fraction Has to be between 0.0 and 1.0, borders included.
 * @param minimalPeakHeight Even negative values are allowed, albeit in general not meaningful.
 * @return Pairs of (mz | width at mz) in ascending order of mz. If no pairs found, the
 *         vector is empty.
 *
 * @throw ms::PreconditionViolation Parameter fraction is out of the required range.
 *
 * @see ms::findBump
 * @see ms::SpectralPeak::lowness
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
MSPP_EXPORT
std::vector<std::pair<SparseSpectrum::Element::first_type, SparseSpectrum::Element::first_type> > 
measureFullWidths(SparseSpectrum::const_iterator first, SparseSpectrum::const_iterator last, double fraction, double minimalPeakHeight = 0);

template<typename FwdIter, typename MzExtractor, typename IntensityExtractor> 
MSPP_EXPORT
std::vector<std::pair<MzExtractor::result_type, MzExtractor::result_type> > 
measureFullWidths2(MzExtrator get_mz, IntensityExtractor get_int, FwdIter first, FwdIter last, double fraction, IntensityExtractor::result_type minimalPeakHeight = 0);



/* implementation */

// findBump()
template<typename FwdIter, typename Compare>
std::pair<FwdIter, FwdIter> findBump(FwdIter first, FwdIter last, Compare comp) {
    // First, let's assume the bump starts right on the first element.
    FwdIter leftEdge = first;
    // Furthermore, we are starting out on the first element, too.
    FwdIter currentElement = first;

    // Since we are working with ForwardIterators, we cannot rewind and have to store the
    // next position separately
    const FwdIter onePastLastElement = last;
    FwdIter onePastCurrentElement = currentElement; std::advance(onePastCurrentElement, 1);   

    // state of our search
    bool onIncreasingSlope = false;
    bool foundBumpTop = false;

    // Now let's go through the sequence element by element.
    while(onePastCurrentElement != onePastLastElement) {
        mspp_invariant(std::distance(currentElement, onePastCurrentElement) == 1, "findBump(): Distance between current and next elment is not 1.");
        
        /* Current element is smaller than its right neighbor. */      
        if(comp(*currentElement, *(onePastCurrentElement))) {
            // We already passed the top of a bump and are finished!
            if(foundBumpTop) {
                // currentElement points to the right edge of the bump
                break;
            }
            else {
                // We are on the bottom of an increasing slope.
                if(!onIncreasingSlope) {
                    onIncreasingSlope = true;
                    // That's were bumps are starting!
                    leftEdge = currentElement;
                }
            }
        }
        
        /* Current element is bigger than its right neighbor. */
        else if(comp(*(onePastCurrentElement), *currentElement)) { 
            // Were we on an increasing slope?
            if(onIncreasingSlope) {
                // So we are on the top of the bump!
                foundBumpTop = true;
            }
            // We are on a decreasing slope.
            else {
                onIncreasingSlope = false;
            }
        }

        /* Current element and its neighbor are equal. */
        else {
            // Great, we are finished!
            if(foundBumpTop) {
                // currentElement points to the right edge of our bump.
                break;
            }
            
            // That's bad. We have to start our search again.
            else {            
                leftEdge = onePastCurrentElement;
                onIncreasingSlope = false;
            }
        }

        ++currentElement;
        ++onePastCurrentElement;
    }
    
    if(foundBumpTop) { 
        // we return [leftEdge, rightEdge]
        return std::make_pair(leftEdge, currentElement);
    }
    else {
        return std::make_pair(last, last);
    }
}



// measureFullWidths()
template<typename FwdIter, typename MzExtractor, typename IntensityExtractor> 
MSPP_EXPORT
std::vector<std::pair<MzExtractor::result_type, MzExtractor::result_type> > 
measureFullWidths2(MzExtrator get_mz, IntensityExtractor get_int, FwdIter first, FwdIter last, double fraction, IntensityExtractor::result_type minimalPeakHeight = 0) {
    typedef MzExtractor::result_type Mz;
    typedef IntensityExtractor::result_type Intensity;

    mspp_precondition(0. <= fraction && fraction <= 1., 
        "measureFullWidths(): Parameter fraction out of required range.");

    // The result
    std::vector<std::pair<Mz, Mz> > widths;

    // Check for empty spectrum or only one element
    if((last - first) < 1) {
        return widths;
    }

    // prerequisites    
    const double requiredLowness = 1. - fraction;
    std::pair<FwdIter, FwdIter> bump;
    Mz positionOfMaximum = 0;
    Intensity bumpHeight = 0;    
    Mz width = 0;  
    SparseSpectrum::LessThanAbundance<SparseSpectrum::Element, SparseSpectrum::Element> comp; // FIXME
    
    // go through all bumps in the spectrum */       
    while(first < last) {
        bump = findBump(first, last, comp);
        // No new bump found
        if (bump.first == last) {
            break;
        }
        
        // calc full width if bump is low enough and has a minimal height
        mspp_invariant(bump.first <= bump.second && bump.second < last, "Bump in illegal state.");
        bumpHeight = get_int(std::max_element(bump.first, bump.second + 1, comp));       
        if(SpectralPeak::lowness(bump.first, bump.second) >= requiredLowness && bumpHeight >= minimalPeakHeight) {
            // we don't have to try for exceptions here, because a bump fulfills the preconditions
            width = SpectralPeak::fullWidthAtFractionOfMaximum(bump.first, bump.second, fraction);
            positionOfMaximum = get_mz(std::max_element(bump.first, bump.second + 1, comp));
            MSPP_LOG(logDEBUG) << "measureFullWidths(): Measured peak (mz | width): (" << positionOfMaximum << " | " << width << ")";
            widths.push_back(std::make_pair(positionOfMaximum, width));
        }
        
        // last element of the bump may be the first of the next one
        first = bump.second;
    }
    
    return widths;
}

} /* namespace ms */

#endif /*__SPECTRUMALGORITHM_H__*/

