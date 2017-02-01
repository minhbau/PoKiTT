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
 *  \author Josh McConnell
 *  \date   January 2017
 *
 *  \brief Computes elements of a matrix used to obtain the the Jacobian used in the
 *         density solve when species mass fractions and enthalpy are transported in
 *         conservation form. This expression calculates
 *  \f[
 *   A_{i,j} = - \phi_i\frac{\partial \rho }{partial \Y_j}
 *  \f]
 *
 *  where \f$\left[A\right] = \left[J\right] + \rho \left[I\right]\f$, and
 *  \f$ \phi=\left\{Y_j,h\right\}\f$
 */
template<typename FieldT>
class PartialJacobian_Species : public Expr::Expression<FieldT>
{
  DECLARE_FIELDS(FieldT, rho_, mmw_ )
  DECLARE_VECTOR_OF_FIELDS(FieldT, phi_ )
  const double mwTerm_;

  PartialJacobian_Species( const Expr::TagList& phiTags,
                           const Expr::Tag&     rhoTag,
                           const Expr::Tag&     mmwTag,
                           const double         mw_i,
                           const double         mw_n );

public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a PartialJacobian_Species expression
     *  @param resultTag the tag for the value that this expression computes
     */
    Builder( const Expr::TagList& resultTags,
             const Expr::TagList& phiTags,
             const Expr::Tag&     rhoTag,
             const Expr::Tag&     mmwTag,
             const double         mw_i,
             const double         mw_n );

    Expr::ExpressionBase* build() const{
      return new PartialJacobian_Species( phiTags_, rhoTag_, mmwTag_, mw_i_, mw_n_);
    }

  private:
    const Expr::TagList phiTags_;
    const Expr::Tag rhoTag_, mmwTag_;
    const double mw_i_, mw_n_;
  };

  void evaluate();

};

//=============================================================================
//=============================================================================
/**
 *  \class PartialJacobian_Enthalpy
 *  \author Josh McConnell
 *  \date   January 2017
 *
 *  \brief Computes elements of a matrix used to obtain the the Jacobian used in the
 *         density solve when species mass fractions and enthalpy are transported in
 *         conservation form. This expression calculates
 *  \f[
 *   A_{i,j} = - \phi_i\frac{\partial \rho }{partial \h}
 *  \f]
 *
 *  where \f$\left[A\right] = \left[J\right] + \rho \left[I\right]\f$, and
 *  \f$ \phi=\left\{Y_j,h\right\}\f$. This used when we solve for h (rather than
 *  temperature).
 */
template<typename FieldT>
class PartialJacobian_Enthalpy : public Expr::Expression<FieldT>
{
  DECLARE_FIELDS(FieldT, rho_, cp_, t_ )
  DECLARE_VECTOR_OF_FIELDS(FieldT, phi_ )

    PartialJacobian_Enthalpy( const Expr::TagList& phiTags,
                              const Expr::Tag&     rhoTag,
                              const Expr::Tag&     cpTag,
                              const Expr::Tag&     temperatureTag );

public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a PartialJacobian_Enthalpy expression
     *  @param resultTag the tag for the value that this expression computes
     */
    Builder( const Expr::TagList& resultTag,
             const Expr::TagList& phiTags,
             const Expr::Tag&     rhoTag,
             const Expr::Tag&     cpTag,
             const Expr::Tag&     temperatureTag );

    Expr::ExpressionBase* build() const{
      return new PartialJacobian_Enthalpy( phiTags_, rhoTag_, cpTag_, temperatureTag_ );
    }

  private:
    const Expr::TagList phiTags_;
    const Expr::Tag rhoTag_, cpTag_, temperatureTag_;
  };

  void evaluate();

};

//=============================================================================
//=============================================================================
/**
 *  \class PartialJacobian_Temperature
 *  \author Josh McConnell
 *  \date   January 2017
 *
 *  \brief Computes elements of a matrix used to obtain the the Jacobian used in the
 *         density solve when species mass fractions and enthalpy are transported in
 *         conservation form. This expression calculates
 *  \f[
 *   A_{i,j} = - \phi_i\frac{\partial \rho }{partial \T}
 *  \f]
 *
 *  where \f$\left[A\right] = \left[J\right] + \rho \left[I\right]\f$, and
 *  \f$ \phi=\left\{Y_j,h\right\}\f$. This used when we solve for temperature (rather
 *  than enthalpy).
 */
template<typename FieldT>
class PartialJacobian_Temperature : public Expr::Expression<FieldT>
{
  DECLARE_FIELDS(FieldT, rho_, t_ )
  DECLARE_VECTOR_OF_FIELDS(FieldT, phi_ )

    PartialJacobian_Temperature( const Expr::TagList& phiTags,
                              const Expr::Tag&     rhoTag,
                              const Expr::Tag&     temperatureTag );

public:
  class Builder : public Expr::ExpressionBuilder
  {
  public:
    /**
     *  @brief Build a PartialJacobian_Temperature expression
     *  @param resultTag the tag for the value that this expression computes
     */
    Builder( const Expr::TagList& resultTag,
             const Expr::TagList& phiTags,
             const Expr::Tag&     rhoTag,
             const Expr::Tag&     temperatureTag );

    Expr::ExpressionBase* build() const{
      return new PartialJacobian_Temperature( phiTags_, rhoTag_, temperatureTag_ );
    }

  private:
    const Expr::TagList phiTags_;
    const Expr::Tag rhoTag_, temperatureTag_;
  };

  void evaluate();

};

#endif /* POKITT_DENSITYSOLVE_JACOBIANELEMENTS_H_ */
