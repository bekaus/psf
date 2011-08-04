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

#include <ms++/config.h>
#include <ms++/SparseSpectrum.h>
#include <psf/SpectrumAlgorithm.h>
#include <psf/LessByExtractor.h>

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
    MSPP_EXPORT double lowness(SparseSpectrum::const_iterator firstElement, SparseSpectrum::const_iterator lastElement);

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
	MSPP_EXPORT
    SparseSpectrum::Element::first_type 
    fullWidthAtFractionOfMaximum(SparseSpectrum::const_iterator firstElement, SparseSpectrum::const_iterator lastElement, const double fraction);

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

} /* namespace ms */

#endif /*__SPECTRALPEAK_H__*/
