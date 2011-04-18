/*$Id: SparseSpectrum.cpp 2681 2009-11-15 22:07:24Z mkirchner $*/

/*
 * SparseSpectrum.cpp
 *
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

// Workaround: Predicate checks cause compile errors in debug mode using MSVC
// See: http://connect.microsoft.com/VisualStudio/feedback/ViewFeedback.aspx?FeedbackID=98738
// Note: This problem should be gone with the new C++0x standard
#ifdef _MSC_VER
# undef _HAS_ITERATOR_DEBUGGING
# define _HAS_ITERATOR_DEBUGGING 0
#endif

#include "ms++/SparseSpectrum.h"

#include <numeric>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#ifndef _MSC_VER
#include <ext/algorithm>
#else
template<typename _ForwardIterator, typename _StrictWeakOrdering>
bool
is_sorted(_ForwardIterator __first, _ForwardIterator __last, _StrictWeakOrdering __comp)
{
    if (__first == __last)
        return true;
    
    _ForwardIterator __next = __first;
    for (++__next; __next != __last; __first = __next, ++__next)
        if (__comp(*__next, *__first))
            return false;
        return true;
}
#endif

#include <ms++/ModelMatrix.h>
#include <ms++/Error.h>

ms::SparseSpectrum::SparseSpectrum()
{
    rt_ = 0;
    msLevel_ = 0;
}

ms::SparseSpectrum::SparseSpectrum(const std::string& filename)
{
    // call default constructor
    SparseSpectrum();
    // read from file
    std::ifstream ifs(filename.c_str());
    if (ifs.good()) {
        ifs >> *this;
    }
}

ms::SparseSpectrum::SparseSpectrum(const std::vector<double>& mz, const std::vector<double>& abundances)
{
    // call default constructor
    SparseSpectrum();
    // read from vectors
    if (mz.size() == abundances.size()) {
        for (unsigned int i = 0; i < mz.size(); i++) {
            c_.push_back(Element(mz[i], abundances[i]));
        }
    }
}

bool ms::SparseSpectrum::operator==(const ms::SparseSpectrum& s) const {
    return c_ == s.c_ && msLevel_ == s.msLevel_ && rt_ == s.rt_;
}

ms::SparseSpectrum ms::SparseSpectrum::subset(const double beginMz, const double endMz) const
{
    mspp_precondition(is_sorted(begin(), end(), LessThanMz<SparseSpectrum::Element, SparseSpectrum::Element>()), "SparseSpectrum must be sorted by mz before trying to get a subset");
    const_iterator first = std::lower_bound(begin(), end(), beginMz, LessThanMz<SparseSpectrum::Element, double>());
    const_iterator last  = std::upper_bound(first, end(), endMz, LessThanMz<double, SparseSpectrum::Element>());
    SparseSpectrum s(first, last);
    s.setRetentionTime(rt_);
    s.setMsLevel(msLevel_);
    return s;
}

void ms::SparseSpectrum::shiftTo(const double to)
{
    double diff = to - mz0();
    transform(c_.begin(), c_.end(), c_.begin(), ms::SparseSpectrum::AddConst(diff));
}

void ms::SparseSpectrum::shiftBy(const double by)
{
    transform(c_.begin(), c_.end(), c_.begin(), ms::SparseSpectrum::AddConst(by));
}

void ms::SparseSpectrum::shiftMaxToMonoisotopicMass(void)
{
    iterator maxIdx = std::max_element(c_.begin(), c_.end(), LessThanAbundance<SparseSpectrum::Element, SparseSpectrum::Element>());
    if (maxIdx != c_.begin()) {
        shiftBy(mz0() - (maxIdx->mz));
    }
}

std::ostream& ms::operator<<(std::ostream& os, ms::SparseSpectrum& p)
{
    for (ms::SparseSpectrum::iterator i = p.begin(); i != p.end(); ++i) {
        os << i->mz << " " << i->abundance << std::endl;
    }
    return(os);
}

std::istream& ms::operator>>(std::istream& is, ms::SparseSpectrum& s)
{
    double mz, ab;
    if (is.good()) {
        while (is >> mz >> ab) {
            // only push_back if abundance is > 0
            if (ab > 0.0) {
                s.push_back(ms::SparseSpectrum::Element(mz, ab));
            }
        }
    }

    return is;
}

ms::SparseSpectrum::const_iterator ms::SparseSpectrum::getMaxAbundancePeak() const
{
    return std::max_element(begin(), end(), ms::SparseSpectrum::LessThanAbundance<SparseSpectrum::Element, SparseSpectrum::Element>());
}

ms::SparseSpectrum::iterator ms::SparseSpectrum::getMaxAbundancePeak()
{
    return std::max_element(begin(), end(), ms::SparseSpectrum::LessThanAbundance<SparseSpectrum::Element, SparseSpectrum::Element>());
}

double ms::SparseSpectrum::getTotalAbundance()
{
    return std::accumulate(begin(), end(), 0.0, ms::SparseSpectrum::PlusAbundance());
}

double ms::SparseSpectrum::getMeanMz()
{
    ModelMatrix matMz, matAbundances;
    getMz(matMz);
    getAbundances(matAbundances);

    return (matMz.transpose()*matAbundances / getTotalAbundance())(0, 0);
}

void ms::SparseSpectrum::merge(const SparseSpectrum& other)
{
    iterator s1 = begin();
    const_iterator s2 = other.begin();

    while (s1 != end() && s2 != other.end()) {
        if (s1->mz < s2->mz) {
            //s1 is behind
            ++s1;
        }
        else if (s1->mz > s2->mz) {
            //s2 is behind
            s1 = insert(s1, Element(s2->mz, s2->abundance));
            ++s2;
        }
        else { /*s1->mz == s2->mz */
            *s1 = Element(s1->mz, s1->abundance + s2->abundance);
            ++s1;
            ++s2;
        }
    }
    // if we terminated because we reached the end of
    // this->sepctrum_, add the remains
    std::copy(s2, other.end(), std::back_inserter(c_));
}

ms::ModelMatrix& ms::SparseSpectrum::getAbundances(ms::ModelMatrix& m) const {
    m.reshape(c_.size(), 1, 0);
    vigra::MultiArrayIndex j = 0;
    for (const_iterator i = c_.begin(); i != c_.end(); ++i, ++j) {
        m(j, 0) = i->abundance;
    }
    return m;
}

ms::ModelMatrix& ms::SparseSpectrum::getMz(ms::ModelMatrix& m) const {
    m.reshape(c_.size(), 1, 0);
    vigra::MultiArrayIndex j = 0;
    for (const_iterator i = c_.begin(); i != c_.end(); ++i, ++j) {
        m(j, 0) = i->mz;
    }
    return m;
}

std::ostream& operator<<(std::ostream& os, ms::SparseSpectrum::Element& e) {
    os << e.mz << " " << e.abundance;
    return os;
}
