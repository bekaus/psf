/*
 * Spectrum.h
 *
 * Copyright (c) 2011 Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
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

#ifndef __SPECTRUM_H__
#define __SPECTRUM_H__

/**
 * @page spectrum Spectrum Sample Implementation
 *
 *
 *
 * @section spectrumdatatye Spectrum data type
 *
 * libpsf can work with any datatype representing a mass spectrum as long as you provide
 * mz and intensity extractors, that adhere to the expected interface.
 * We provide a sample implementation of a spectrum. Feel free to use your own.
 *
 *
 * @section extractorinterface Extractor interface
 *
 * struct MyExtractor {
 *  typedef ... result_type; // numeric type used for the value (either mz or intensity)
 *  typedef ... element_type; // type of entries in a spectrum
 *
 *  result_type operator()( const element_type& e ) const;
 * };
 *
 *
 *
 * @author Bernhard X. Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 * @date 2011
 */


namespace ms
{


/**
* a single entry in a mass spectrum, characterized by \f$mz\f$ and intensity
* values
*/
struct SpectrumElement {
    typedef double mz_type;
    typedef double intensity_type;

    SpectrumElement(const double m, const double i) : mz(m), intensity(i) {}

    double mz;
    double intensity;
};

struct MzExtractor {
  typedef SpectrumElement::mz_type result_type;
  typedef SpectrumElement element_type;

  result_type operator()( const SpectrumElement& e ) const {
    return e.mz;
  };
};

struct IntensityExtractor {
  typedef SpectrumElement::intensity_type result_type;
  typedef SpectrumElement element_type;

  result_type operator()( const SpectrumElement& e ) const {
    return e.intensity;
  };
};



/**
 * A mass spectrum is a sequence of SpectrumElements ordered by mz.
 */ 
typedef std::vector<SpectrumElement> Spectrum;

} /* namespace ms */

#endif /*__SPECTRUM_H__*/
