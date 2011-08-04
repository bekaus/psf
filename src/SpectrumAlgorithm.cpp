/* $Id: SpectrumAlgorithm.cpp 2307 2009-09-01 13:14:17Z bkausler $ */

/*
 * SpectrumAlgorithm.cpp
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

#include "psf/SpectrumAlgorithm.h"

#include <algorithm>
#include <utility>
#include <vector>

#include <ms++/Error.h>
#include <ms++/Log.h>
#include <ms++/SparseSpectrum.h>
#include <psf/SpectralPeak.h>

using namespace ms;

// measureFullWidths()
std::vector<std::pair<SparseSpectrum::Element::first_type, SparseSpectrum::Element::first_type> >
ms::measureFullWidths(SparseSpectrum::const_iterator first, SparseSpectrum::const_iterator last, const double fraction, const double minimalPeakHeight) {
    typedef SparseSpectrum::Element::first_type Mz;
    typedef SparseSpectrum::Element::second_type Abundance;

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
    std::pair<SparseSpectrum::const_iterator, SparseSpectrum::const_iterator> bump;
    Mz positionOfMaximum = 0;
    Abundance bumpHeight = 0;    
    Mz width = 0;  
    SparseSpectrum::LessThanAbundance<SparseSpectrum::Element, SparseSpectrum::Element> comp;
    
    // go through all bumps in the spectrum */       
    while(first < last) {
        bump = findBump(first, last, comp);
        // No new bump found
        if (bump.first == last) {
            break;
        }
        
        // calc full width if bump is low enough and has a minimal height
        mspp_invariant(bump.first <= bump.second && bump.second < last, "Bump in illegal state.");
        bumpHeight = std::max_element(bump.first, bump.second + 1, comp)->abundance;       
        if(SpectralPeak::lowness(bump.first, bump.second) >= requiredLowness && bumpHeight >= minimalPeakHeight) {
            // we don't have to try for exceptions here, because a bump fulfills the preconditions
            width = SpectralPeak::fullWidthAtFractionOfMaximum(bump.first, bump.second, fraction);
            positionOfMaximum = std::max_element(bump.first, bump.second + 1, comp)->mz;
            MSPP_LOG(logDEBUG) << "measureFullWidths(): Measured peak (mz | width): (" << positionOfMaximum << " | " << width << ")";
            widths.push_back(std::make_pair(positionOfMaximum, width));
        }
        
        // last element of the bump may be the first of the next one
        first = bump.second;
    }
    
    return widths;
}
