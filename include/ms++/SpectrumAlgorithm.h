/*$Id: SpectrumAlgorithm.h 2613 2009-10-26 23:34:00Z bkausler $*/

/*
 * SpectrumAlgorithm.h
 *
 * Copyright (c) 2009 Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * Copyright (c) 2008 Marc Kirchner <marc.kirchner@iwr.uni-heidelberg.de>
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

#include <ms++/Error.h>
#include <ms++/SparseSpectrum.h>

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



// class SpectrumSplitter
/**
 * Breaks up a sparse spectrum at gaps.
 */
class MSPP_EXPORT SpectrumSplitter
{
public:
    std::vector<SparseSpectrum> operator()(SparseSpectrum& s, const double deltamz);
private:
    struct Diff {
        double operator()(const SparseSpectrum::value_type& lhs, const SparseSpectrum::value_type& rhs) {
            return(lhs.mz - rhs.mz);
        }
    };
};



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

} /* namespace ms */

#endif /*__SPECTRUMALGORITHM_H__*/

