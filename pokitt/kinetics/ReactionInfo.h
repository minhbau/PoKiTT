/*
 * The MIT License
 *
 * Copyright (c) 2016-2017 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/**
 *  \file   ReactionInfo.h
 *  \date   Oct 15, 2016
 *  \author mike
 */

#ifndef REACTIONINFO_H_
#define REACTIONINFO_H_

#include <pokitt/CanteraObjects.h>

namespace pokitt
{

  // the functional form for rate constant k
  // if the activation energy is 0 and "b" is an integer, we can avoid evaluating an exponential
  enum RateConstantForm{
    CONSTANT,
    LINEAR,
    QUADRATIC,
    RECIPROCAL,
    ARRHENIUS
  };

  // which terms are non-negligible in evaluating Troe falloff
  enum TroeTerms{
    NONE,
    T1,
    T2,
    T12,
    T3,
    T13,
    T23,
    T123
  };

  // similar to the sum macros, this multiplies concentrations in a single kernel
  enum ReactionOrder{
    ONE,
    TWO,
    ONE_ONE,
    TWO_ONE,
    ONE_TWO,
    ONE_ONE_ONE,
    OTHER
  };

  typedef std::vector< RxnData::SpeciesRxnData > SpecDataVecT; // from CanteraObjects.h

  struct ReactionInfo{
    ReactionInfo( const RxnData& data );
    RateConstantForm kForm;
    TroeTerms troeForm;
    ReactionOrder forwardOrder;
    ReactionOrder reverseOrder;
    int sumReactantStoich;
    int sumProductStoich;
  };

}

#endif /* REACTIONINFO_H_ */
