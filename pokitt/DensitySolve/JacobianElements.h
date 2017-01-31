/*
 * JacobianElements.h
 *
 *  Created on: Jan 30, 2017
 *      Author: josh
 */

#ifndef POKITT_DENSITYSOLVE_JACOBIANELEMENTS_H_
#define POKITT_DENSITYSOLVE_JACOBIANELEMENTS_H_

#include <expression/Expression.h>

/**
 *  \class PartialJacobian_Species
 *  \brief calculate values of each element of the Jacobian
 */
template<typename FieldT>
class PartialJacobian_Species : public Expr::Expression<FieldT>
{
  DECLARE_FIELDS(FieldT, rhoGuess_, phi_, mmw_ )
  const double mwTerm_;

  PartialJacobian_Species( const Expr::Tag& rhoGuessTag,
                           const Expr::Tag& phiTag,
                           const Expr::Tag& mmwTag,
                           const double     mw_i,
                           const double     mw_n );

public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a PartialRhoPartialYi expression
     *  @param resultTag the tag for the value that this expression computes
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& rhoOldTag,
             const Expr::Tag& phiTag,
             const Expr::Tag& mmwTag,
             const double     mw_i,
             const double     mw_n );

    Expr::ExpressionBase* build() const{
      return new PartialJacobian_Species( rhoOldTag_, phiTag_, mmwTag_, mw_i_, mw_n_);
    }

  private:
    const Expr::Tag rhoOldTag_, phiTag_, mmwTag_;
    const double mw_i_, mw_n_;
  };

  void evaluate();

};

/**
 *  \class PartialJacobian_Enthalpy
 *  \brief calculate values of each element of the Jacobian
 */
template<typename FieldT>
class PartialJacobian_Enthalpy : public Expr::Expression<FieldT>
{
  DECLARE_FIELDS(FieldT, rhoGuess_, phi_, cp_, t_ )

    PartialJacobian_Enthalpy( const Expr::Tag& rhoGuessTag,
                              const Expr::Tag& phiTag,
                              const Expr::Tag& cpTag,
                              const Expr::Tag& temperatureTag );

public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a PartialJacobian_Enthalpy expression
     *  @param resultTag the tag for the value that this expression computes
     */
    Builder( const Expr::Tag& resultTag,
             const Expr::Tag& rhoGuessTag,
             const Expr::Tag& phiTag,
             const Expr::Tag& cpTag,
             const Expr::Tag& temperatureTag );

    Expr::ExpressionBase* build() const{
      return new PartialJacobian_Enthalpy( rhoOldTag_, phiTag_, cpTag_, temperatureTag_ );
    }

  private:
    const Expr::Tag rhoOldTag_, phiTag_, cpTag_, temperatureTag_;
  };

  void evaluate();

};




#endif /* POKITT_DENSITYSOLVE_JACOBIANELEMENTS_H_ */
