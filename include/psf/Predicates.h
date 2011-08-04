/*
 * Predicates.h
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

#ifndef __PREDICATES_H__
#define __PREDICATES_H__

namespace ms
{

// class LessByExtractor
/**
 * Compare two elements with regard to a certain aspect.
 *
 * Typically, an element in a spectrum represents more than one value 
 * (such as m/z, intensity, time etc.) Use this functor to compare two
 * elements with regard to one of these values.
 */
template< typename Element, typename Extractor >
class MSPP_EXPORT LessByExtractor {
  public:
  LessByExtractor( const Extractor& e ) : extract_(e) {};
  bool operator()( const Element& lhs, const Element& rhs ) const {
    return extract_(lhs) < extract_(rhs);    
  };

  private:
  Extractor extract_;
};



template< typename Element, typename Extractor >
class MSPP_EXPORT MoreThanValue {
   public:
   MoreThanValue( const Extractor& e, typename Extractor::result_type val ) : extract_(e), val_(val) {};
   bool operator()( const Element& e ) const {
     return val_ < extract_(e);    
  };

  private:
  Extractor extract_;
  typename Extractor::result_type val_;
};

} /* namespace ms */

#endif /*__PREDICATES_H__*/
