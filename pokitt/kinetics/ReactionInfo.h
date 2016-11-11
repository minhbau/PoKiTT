/*
 * The MIT License
 *
 * Copyright (c) 2016 The University of Utah
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


namespace pokitt
{

  // the functional form for rate constant k
  // if the activation energy is 0 and "b" is an integer, we can avoid evaluating an exponentional
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

  ReactionInfo::ReactionInfo( const RxnData& dat )
  {
    const SpecDataVecT& rxts = dat.reactants;
    switch( rxts.size() ){
      case 1:
        if(      rxts[0].stoich == 1 ) forwardOrder = ONE;
        else if( rxts[0].stoich == 2 ) forwardOrder = TWO;
        else                           forwardOrder = OTHER;
        break;
      case 2:
        if(      rxts[0].stoich == 1 && rxts[1].stoich == 1 ) forwardOrder = ONE_ONE;
        else if( rxts[0].stoich == 1 && rxts[1].stoich == 2 ) forwardOrder = ONE_TWO;
        else if( rxts[0].stoich == 2 && rxts[1].stoich == 1 ) forwardOrder = TWO_ONE;
        else                                                  forwardOrder = OTHER;
        break;
      case 3:
        if( rxts[0].stoich == 1 && rxts[1].stoich == 1 && rxts[2].stoich == 1 ) forwardOrder = ONE_ONE_ONE;
        else                                                                    forwardOrder = OTHER;
        break;
      default:
        forwardOrder = OTHER;
    }
    sumReactantStoich = 0;
    for( size_t s=0; s<rxts.size(); ++s ){
      sumReactantStoich += std::abs( rxts[s].stoich );
    }

    const SpecDataVecT& prods = dat.products;
    switch( prods.size() ){
      case 1:
        if(      prods[0].stoich == -1 ) reverseOrder = ONE;
        else if( prods[0].stoich == -2 ) reverseOrder = TWO;
        else                             reverseOrder = OTHER;
        break;
      case 2:
        if(      prods[0].stoich == -1 && prods[1].stoich == -1 ) reverseOrder = ONE_ONE;
        else if( prods[0].stoich == -1 && prods[1].stoich == -2 ) reverseOrder = ONE_TWO;
        else if( prods[0].stoich == -2 && prods[1].stoich == -1 ) reverseOrder = TWO_ONE;
        else                                                      reverseOrder = OTHER;
        break;
      case 3:
        if( prods[0].stoich == -1 && prods[1].stoich == -1 && prods[2].stoich == -1 ) reverseOrder = ONE_ONE_ONE;
        else                                                                          reverseOrder = OTHER;
        break;
      default:
        reverseOrder = OTHER;
    }
    sumProductStoich = 0;
    for( size_t s=0; s<prods.size(); ++s ){
      sumProductStoich += std::abs( prods[s].stoich );
    }

    /*
     * Here we are checking if the rate constant is Arrhenius in form
     * If it is, we need to perform an exponential evaluation (expensive)
     * If not, then we determine which power of temperature to use
     */
    if( fabs( dat.kFwdCoefs[2] ) < 1e-6 ){ // i.e. 0 activation energy
      if(      fabs( dat.kFwdCoefs[1]     ) < 1e-6 ) kForm = CONSTANT;   // i.e. b=0
      else if( fabs( dat.kFwdCoefs[1] - 1 ) < 1e-6 ) kForm = LINEAR;     // i.e. b=1
      else if( fabs( dat.kFwdCoefs[1] - 2 ) < 1e-6 ) kForm = QUADRATIC;  // i.e. b=2
      else if( fabs( dat.kFwdCoefs[1] + 1 ) < 1e-6 ) kForm = RECIPROCAL; // i.e. b=-1
      else kForm = ARRHENIUS; // activation energy but non-integer value for b
    }
    else kForm = ARRHENIUS;

    // now we check which terms of the troe function are negligible
    const double* troe = dat.troeParams;
    switch( dat.type ){
      case ELEMENTARY: troeForm = NONE; break;
      case THIRD_BODY: troeForm = NONE; break;
      case LINDEMANN:  troeForm = NONE; break;
      case TROE:
        troeForm = NONE;
        if( std::abs(troe[1]) > 1e-8 ){  // 1 is present
          if( std::abs(troe[2]) > 1e-8 ){ // 1 and 2 are present
            if( std::abs(troe[3]) > 1e-8 )    troeForm = T123;
            else                              troeForm = T12;
          }
          else if( std::abs(troe[3]) > 1e-8 ) troeForm = T13;
          else                                troeForm = T1;
        }
        else{ // 1 is not present
          if( std::abs(troe[2]) > 1e-8 ){ // 1 not present, 2 is present
            if( std::abs(troe[3]) > 1e-8 )    troeForm = T23;
            else                              troeForm = T2;
          }
          else if( std::abs(troe[3]) > 1e-8 ) troeForm = T3;
        }
        break;
      default: {
        std::ostringstream msg;
        msg << __FILE__ << " : " << __LINE__
            << "\n Unknown reaction type somehow evaded detection\n"
            << "This should have been caught during construction of CanteraObjects\n"
            << "Check your xml input file to ensure it is not corrupted\n";
        throw std::runtime_error( msg.str() );
      }
    } // switch( dat.type )
  }

}



#endif /* REACTIONINFO_H_ */
