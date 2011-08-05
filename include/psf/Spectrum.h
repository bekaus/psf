#ifndef __SPECTRUM_H__
#define __SPECTRUM_H__

#include <fstream>
#include <istream>
#include <string>

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


namespace psf
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

std::istream& operator>>(std::istream& is, Spectrum& s) {
    double mz, intensity;
    if (is.good()) {
        while (is >> mz >> intensity) {
	    if (intensity > 0) {
		s.push_back(SpectrumElement(mz, intensity));
	    }
        }
    }
    return is;
}

void loadSpectrumElements(Spectrum& s, const std::string& filename){
    std::ifstream ifs(filename.c_str());
    if (ifs.good()) {
        ifs >> s;
    }
}


} /* namespace psf */

#endif /*__SPECTRUM_H__*/
