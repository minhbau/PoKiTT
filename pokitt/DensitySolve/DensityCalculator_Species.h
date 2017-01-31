/*
 * The MIT License
 *
 * Copyright (c) 2012-2016 The University of Utah
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

#ifndef DensityCalculator_Species_h
#define DensityCalculator_Species_h

//#include <sci_defs/wasatch_defs.h>
//
//#ifndef HAVE_POKITT
//// kill compilation if we don't have pokitt.
//#error Solving density from species mass fractions requires PoKitt.
//#endif

#include <expression/ExprLib.h>
#include <expression/ExprPatch.h>
#include <expression/Expression.h>
#include <expression/ExpressionFactory.h>
#include <expression/FieldManagerList.h>
#include <spatialops/structured/MatVecFields.h>


namespace WasatchCore{

/**
 *  \class DensityFromSpecies
 *  \author Josh McConnell
 *  \date January 2017
 *
 *  \brief When transporting \f$\rho Y_i\f$ and \f$\rho h\f$, this expression will
 *  calculate \f$\rho\f$,  \f$\Y_i\f$, and \f$\h\f$.  Note that \f$f\f$ and \f$h\f$ are
 *  calculated by other expressions once \f$\rho\f$ is known.
 *
 *  Given the density-weighted species mass fractions, \f$\rho Y_i\f$ and density-weighted
 *  enthalpy, \f$\rho h\f$, this expression finds \f$\rho\f$, \f$\h\f$, and each \f$\Y_i\f$
 *  heat loss \f$\gamma\f$. A guess for density is given in the following form:
 *  \f[
 *    \rho = \mathcal{G}_\rho (Y_j, h) = \frac{pM}{RT}
 *  \f]
 *  In residual form, the equations to be solved are
 *  \f[
 *   r_{\phi_i} = (\rho phi_i) \mathcal{G}_\rho - phi_i
 *  \f]
 *  Where \f$\phi = \left\{ Y_i,h \}\right\f$
 *  and elements of the Jacobian matrix are
 *  \f[
 *  J_{ij}=\frac{\partial r_i}{\partial \phi_j}
 *  \f]
 */
template< typename FieldT >
class DensityFromSpecies : public Expr::Expression<FieldT>
{
  DECLARE_FIELDS(FieldT, rhoOld_, hOld_, rhoH_, temp_, pressure_)
  DECLARE_VECTOR_OF_FIELDS(FieldT, yiOld_)
  DECLARE_VECTOR_OF_FIELDS(FieldT, rhoYi_)

  //
  Expr::ExpressionFactory  localFactory_;
  Expr::ExprPatch*         patchPtr_;
  Expr::TimeStepper*       integratorPtr_;

  bool setupHasRun_;

  // tags to fields that exist on a local patch
  const Expr::Tag hGuessTag_, rhoGuessTag_, mmwTag_, tGuessTag_, pTag_;
  Expr::TagList yiGuessTags_;

  void setup();

  const int nSpec_, nEq_;
  const std::vector<double> mw_;
  std::vector<double> mwInv_;

  DensityFromSpecies( const Expr::Tag&     rhoOldTag,
                      const Expr::Tag&     hOldTag,
                      const Expr::Tag&     rhoHTag,
                      const Expr::Tag&     tOldTag,
                      const Expr::Tag&     pTag,
                      const Expr::TagList& yiOldTags,
                      const Expr::TagList& rhoYiTags );

public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a DensityFromSpecies expression
     *  @param resultTag the tag for the value that this expression computes
     */
    Builder( const Expr::Tag&     resultTag,
             const Expr::Tag&     rhoOldTag,
             const Expr::Tag&     hOldTag,
             const Expr::Tag&     rhoHTag,
             const Expr::Tag&     tOldTag,
             const Expr::Tag&     pTag,
             const Expr::TagList& yiOldTags,
             const Expr::TagList& rhoYiTags,
             const double         rtol,
             const unsigned       maxiter );

    Expr::ExpressionBase* build() const{
      std::cout<<"Expr::build()\n";
      return new DensityFromSpecies<FieldT>( rhoOldTag_, hOldTag_, rhoHTag_, tOldTag_, pTag_, yiOldTags_, rhoYiTags_ );
    }

    private:
      const Expr::Tag rhoOldTag_, hOldTag_, rhoHTag_, tOldTag_, pTag_;
      const Expr::TagList yiOldTags_, rhoYiTags_;
      const double rTol_;       // relative error tolerance
      const unsigned maxIter_; // maximum number of iterations
    };

  void evaluate();
  ~DensityFromSpecies();
  };
}//namespace WasatchCore

#endif /*DensityCalculator_Species_h */
