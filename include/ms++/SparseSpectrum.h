/*$Id: SparseSpectrum.h 2618 2009-10-28 09:06:40Z bkausler $*/

/*
* SparseSpectrum.h
*
* Copyright (c) 2009 Marc Kirchner <marc.kirchner@childrens.harvard.edu>
* Copyright (c) 2007, 2008 Marc Kirchner <marc.kirchner@iwr.uni-heidelberg.de>
* Copyright (c) 2007 Bjoern Voss <bjoern.voss@iwr.uni-heidelberg.de>
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

// Workaround: "Safe STL" in MSVC causes compile time errors.
#ifdef _MSC_VER
# undef _HAS_ITERATOR_DEBUGGING
# define _HAS_ITERATOR_DEBUGGING 0
#
# undef _SECURE_SCL
# define _SECURE_SCL 0
#endif 

#ifndef __SPARSESPECTRUM_H__
#define __SPARSESPECTRUM_H__
#include <ms++/config.h>

#include <algorithm>
#include <vector>
#include <functional>
#include <string>
#ifdef _MSC_VER
#include <cassert>
#endif

#include <ms++/Collection.h>
#include <ms++/ModelMatrix.h>

class SparseSpectrumTest;

namespace ms
{
 
/**
* a single entry in a mass spectrum, characterized by \f$mz\f$ and abundance
* values
*/
struct SparseSpectrumElement {
    typedef double first_type;
    typedef double second_type;

    SparseSpectrumElement(const double m, const double a) : mz(m), abundance(a) {}
    SparseSpectrumElement(const std::pair<double, double>& p) : mz(p.first), abundance(p.second) {}

    bool operator==(const SparseSpectrumElement& p) const {
        return mz == p.mz && abundance == p.abundance;
    }

    double mz;
    double abundance;
};

/**
* rAa sparse representation of a mass spectrum.
*/
class MSPP_EXPORT SparseSpectrum : public Collection<SparseSpectrumElement>
{
public:
	friend class ::SparseSpectrumTest;
    
    typedef SparseSpectrumElement Element;

    /** @return the SparseSpectrum containing the m/z range [beginMz, endMz)
    *  @param beginMz The lower bound of mz range
    *  @param endMz The upper bound of mz range
    *  @return A new SparseSpectrum object
    */
    SparseSpectrum subset(const double beginMz, const double endMz) const;

    /** orders two ms::SparseSpectrum::Element based on their \f$mz\f$ values
    */
    template<class A, class B>
    struct LessThanMz : std::binary_function<A, B, bool> {
        bool operator()(const A& lhs, const B& rhs);
    };
    /** orders two ms::SparseSpectrum::Element based on their abundance values
    */
    template<class A, class B>
    struct LessThanAbundance : std::binary_function<A, B, bool> {
        bool operator()(const A& lhs, const B& rhs);
    };

    /** Compare two ms::SparseSpectrum objects based on their retention time
    */
    template<class A, class B>
    struct LessThanRt : std::binary_function<A, B, bool> {
        bool operator()(const A& lhs, const B& rhs);
    };

    /** Compare two ms::SparseSpectrum objects based on their MS level
    */
    template<class A, class B>
    struct LessThanMsLevel : std::binary_function<A, B, bool> {
        bool operator()(const A& lhs, const B& rhs);
    };

    /** Check if two ms::SparseSpectrum objects have the same MS level
    */
    struct EqualMsLevel : std::unary_function<SparseSpectrum, bool> {
        EqualMsLevel(unsigned int msLevel) : msLevel_(msLevel) {}
        bool operator()(const SparseSpectrum& s) const {
            return msLevel_ == s.msLevel_;
        }
        unsigned int msLevel_;
    };

    /** add the abundance of element @a e to value @a val
    */
    struct PlusAbundance : std::binary_function<double, ms::SparseSpectrum::value_type, double> {
        double operator()(const double& val, const ms::SparseSpectrum::value_type& e) {
            return val + e.abundance;
        }
    };

    SparseSpectrum();
    SparseSpectrum(const std::vector<double>& mz, const std::vector<double>& abundances);

    // SparseSepctrum(string)
    /**
     * Construct SparseSpectrum from a file.
     *
     * The format of the file consists of entries of two whitespace separated
     * floating-point numbers.
     * The first number corresponds to the m/z value in Thomson and the second to an
     * absolute abundance. No special order of the entries is assumed.
     *
     * For example:<br />
     * <CODE>
     * 600.043986 0.000000<br />
     * 600.045374 0.000013<br />
     * 600.046762 1359.585662<br />
     * 600.048150 1689.438339<br />
     * 600.049539 1848.110212<br />
     * 600.050927 1823.459819<br />
     * 600.052315 1636.566124<br />
     * 600.053703 1294.969852<br />
     * 600.055091 0.000000<br />
     * 600.056479 0.000000<br />
     * </CODE>
     *
     * @param filename The filename optionally preceded by an absolute or relative path.
     */
    SparseSpectrum(const std::string& filename);

    SparseSpectrum(iterator first, iterator last) : Collection<SparseSpectrumElement>(first, last) {}
    SparseSpectrum(const_iterator first, const_iterator last) : Collection<SparseSpectrumElement>(first, last) {}
    ~SparseSpectrum() {}
    
    bool operator==(const SparseSpectrum& s) const;
    
    void clear() {
        c_.clear();
        rt_ = 0.0;
        msLevel_ = 0;
    }
    void shiftTo(const double to);
    void shiftBy(const double diff);
    void shiftMaxToMonoisotopicMass();
    double mz0() const {
        return (*(c_.begin())).mz;
    }

    void setRetentionTime(const double rt) {
        rt_ = rt;
    }
    double getRetentionTime() const {
        return rt_;
    }
    void setMsLevel(const unsigned int l) {
        msLevel_ = l;
    }
    unsigned int getMsLevel() const {
        return msLevel_;
    }

    ModelMatrix& getAbundances(ModelMatrix& m) const;
    ModelMatrix& getMz(ModelMatrix& m) const;

    /** Get the maximum abundance peak
    *   @return Iterator to the maximum abundance element
    */
    iterator getMaxAbundancePeak();
    const_iterator getMaxAbundancePeak() const;

    /** Get the sum of abundance
    *   @return The sum of abundance
    */
    double getTotalAbundance();

    /** Get the mean mz weighted by abundance
    *   @return The mean mz weighted by abundance
    */
    double getMeanMz();

    /** Merges two SparseSpectrum objects.
    *   @param SparseSpectrum instance to merge with.
    *   @pre Requires the current instance and the merge instance to be sorted.
    *
    * This merges two SparseSpectra, combining their abundances if
    * they exhibit identical masses.
    */
    void merge(const SparseSpectrum& ss);

    /** Merges two SparseSpectrum objects.
    *   @param SparseSpectrum instance to merge with.
    *   @pre Requires the current instance and the merge instance to be sorted.
    *   @deprecated Please use SparseSpectrum::merge() instead.
    *   @see SparseSpectrum::merge
    *
    * This merges two SparseSpectra, combining their abundances if
    * they exhibit identical masses.
    */
    void mergeWith(const SparseSpectrum& other) { merge(other); }

    /** Splice a SparseSpectrum into two spectra.
    *   @param first Iterator to first SparseSpectrum::Element that is 
    *          supposed to be moved
    *   @param last Iterator to the first SparseSepctrum::Element that
    *          is NOT to be moved (or end(), following STL semantics)
    *   @note The elements are copied (the underlying data structure 
    *         is a vector, not a list)
    */
    template <typename Out>
    void splice(const iterator first, const iterator last, Out out) {
        std::copy(first, last, out);
        c_.erase(first, last);
    }

    class MzAccessor {
        public:
        typedef double value_type;
        MzAccessor(std::vector<Element> *spectrum) { s = spectrum; }
        double& operator[](size_type pos) const {
            return (*s)[pos].mz;
        }
        size_type size() const { return s->size(); }
        protected:
            std::vector<Element>* s;
    };
    MzAccessor mzAccessor() const { return MzAccessor(const_cast<std::vector<Element>*>(&c_)); }
    
    /**
    * functor to add a constant \f$mz\f$ value to a ms::SparseSpectrum::Element
    */
    struct AddConst : std::unary_function<value_type, value_type> {
        AddConst(double val):val_(val) {}
        value_type operator()(const value_type& p) {
			return(Element(p.mz + val_, p.abundance));
		}
        double val_;
    };

    private:
        double rt_;
        unsigned int msLevel_;

};

MSPP_EXPORT std::istream& operator>>(std::istream& is, ms::SparseSpectrum& s);
MSPP_EXPORT std::ostream& operator<<(std::ostream& os, ms::SparseSpectrum& s);
MSPP_EXPORT std::ostream& operator<<(std::ostream& os, ms::SparseSpectrum::Element& e);

template<>
struct SparseSpectrum::LessThanAbundance<SparseSpectrum::Element, SparseSpectrum::Element> :
std::binary_function<Element, Element, bool> {
    bool operator()(const Element& lhs, const SparseSpectrum::Element& rhs) const {
        return lhs.abundance < rhs.abundance;
    }
};
template<>
struct SparseSpectrum::LessThanAbundance<SparseSpectrum::Element, double> :
std::binary_function<Element, double, bool> {
    bool operator()(const Element& lhs, const double& rhs) const {
        return lhs.abundance < rhs;
    }
};
template<>
struct SparseSpectrum::LessThanAbundance<double, SparseSpectrum::Element> :
std::binary_function<double, Element, bool> {
    bool operator()(const double& lhs, const Element& rhs) const {
        return lhs < rhs.abundance;
    }
};

template<>
struct SparseSpectrum::LessThanMz<SparseSpectrum::Element, SparseSpectrum::Element> :
std::binary_function<Element, Element, bool> {
    bool operator()(const Element& lhs, const Element& rhs) const {
        return lhs.mz < rhs.mz;
    }
};
template<>
struct SparseSpectrum::LessThanMz<SparseSpectrum::Element, double> :
std::binary_function<Element, double, bool> {
    bool operator()(const Element& lhs, const double& rhs) const {
        return lhs.mz < rhs;
    }
};
template<>
struct SparseSpectrum::LessThanMz<double, SparseSpectrum::Element> :
std::binary_function<double, Element, bool> {
    bool operator()(const double& lhs, const Element& rhs) const {
        return lhs < rhs.mz;
    }
};

template<>
struct SparseSpectrum::LessThanMsLevel<SparseSpectrum, SparseSpectrum> :
std::binary_function<SparseSpectrum, SparseSpectrum, bool> {
    bool operator()(const SparseSpectrum& lhs, const SparseSpectrum& rhs) const {
        return lhs.getMsLevel() < rhs.getMsLevel();
    }
};
template<>
struct SparseSpectrum::LessThanMsLevel<SparseSpectrum, double> :
std::binary_function<SparseSpectrum, double, bool> {
    bool operator()(const SparseSpectrum& lhs, const double& rhs) const {
        return lhs.getMsLevel() < rhs;
    }
};
template<>
struct SparseSpectrum::LessThanMsLevel<double, SparseSpectrum> :
std::binary_function<double, SparseSpectrum, bool> {
    bool operator()(const double& lhs, const SparseSpectrum& rhs) const {
        return lhs < rhs.getMsLevel();
    }
};

template<>
struct SparseSpectrum::LessThanRt<SparseSpectrum, SparseSpectrum> :
std::binary_function<SparseSpectrum, SparseSpectrum, bool> {
    bool operator()(const SparseSpectrum& lhs, const SparseSpectrum& rhs) const {
        return lhs.getRetentionTime() < rhs.getRetentionTime();
    }
};
template<>
struct SparseSpectrum::LessThanRt<SparseSpectrum, double> :
std::binary_function<SparseSpectrum, double, bool> {
    bool operator()(const SparseSpectrum& lhs, const double& rhs) const {
        return lhs.getRetentionTime() < rhs;
    }
};
template<>
struct SparseSpectrum::LessThanRt<double, SparseSpectrum> :
std::binary_function<double, SparseSpectrum, bool> {
    bool operator()(const double& lhs, const SparseSpectrum& rhs) const {
        return lhs < rhs.getRetentionTime();
    }
};

} /* namespace ms */


#endif /*__SPARSESPECTRUM_H__*/
