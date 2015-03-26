/*
 * TagNames.h
 *
 *  Created on: Feb 24, 2015
 *      Author: Nathan
 */

#ifndef EXAMPLES_REACTIONDIFFUSION_TAGMANAGER_H_
#define EXAMPLES_REACTIONDIFFUSION_TAGMANAGER_H_

#include <boost/lexical_cast.hpp>
#include <expression/Tag.h>

namespace pokitt{

enum TagName{
  RHO, YI, RHOYI, MMW,

  T, P, H,

  R, D, LAM,

  J, Q, QX, QY,

  XCOORD, YCOORD
};

enum TagListName{
  YI_N, RHOYI_N,

  H_N,

  R_N, D_N,

  JX_N, JY_N
};

class TagManager{

  const Expr::Tag rho_, yi_, rhoYi_, mmw_;
  const Expr::Tag t_,   p_,  h_;
  const Expr::Tag r_,   d_,  lam_;
  const Expr::Tag j_,   q_,  qx_,    qy_;
  const Expr::Tag x_,   y_;

  Expr::TagList yiTags_, rhoYiTags_;
  Expr::TagList hiTags_;
  Expr::TagList rTags_,  dTags_;
  Expr::TagList jxTags_, jyTags_;

  const Expr::Tag empty_;
  const Expr::TagList emptyList_;

public:

  const Expr::Tag& operator[]( TagName name ) const{
    switch( name ){
    case RHO:    return rho_;
    case YI:     return yi_;
    case RHOYI:  return rhoYi_;
    case MMW:    return mmw_;
    case T:      return t_;
    case P:      return p_;
    case H:      return h_;
    case R:      return r_;
    case D:      return d_;
    case LAM:    return lam_;
    case J:      return j_;
    case Q:      return q_;
    case QX:     return qx_;
    case QY:     return qy_;
    case XCOORD: return x_;
    case YCOORD: return y_;
    default:     return empty_;
    }
  }

  const Expr::TagList& operator[]( TagListName name ) const{
    switch( name ){
    case YI_N:    return yiTags_;
    case RHOYI_N: return rhoYiTags_;
    case H_N:     return hiTags_;
    case R_N:     return rTags_;
    case D_N:     return dTags_;
    case JX_N:    return jxTags_;
    case JY_N:    return jyTags_;
    default:      return emptyList_;
    }
  }

  TagManager( const int nSpec,
              const Expr::Tag& rhoYi, const Expr::Tag& h, const Expr::Tag& yi, const Expr::Tag& t,
              const Expr::Tag& p,
              const Expr::Tag& rho,
              const Expr::Tag& mmw,
              const Expr::Tag& r,
              const Expr::Tag& lam,
              const Expr::Tag& d,
              const Expr::Tag& j,
              const Expr::Tag& q,
              const Expr::Tag& x,
              const Expr::Tag& y ):
                rhoYi_(rhoYi),
                h_(h),
                yi_(yi),
                t_(t),
                p_(p),
                rho_(rho),
                mmw_(mmw),
                r_(r),
                lam_(lam),
                d_(d),
                j_(j),
                q_(q),
                x_(x),
                y_(y),
                qx_( Expr::Tag( q, "x" ) ),
                qy_( Expr::Tag( q, "y" ) ),
                empty_( Expr::Tag() ),
                emptyList_( Expr::tag_list( Expr::Tag() ) )

  {
    using Expr::Tag;
    for( size_t n=0; n<nSpec; ++n ){
      std::string spec = boost::lexical_cast<std::string>(n);
      hiTags_.push_back( Tag( h_.name() + spec, Expr::STATE_NONE ) );
      yiTags_.push_back( Tag( yi_, spec         ) );
      jxTags_.push_back( Tag( j_,  "x" + spec   ) );
      jyTags_.push_back( Tag( j_,  "y" + spec   ) );
      rTags_.push_back(  Tag( r_,  spec         ) );
      dTags_.push_back(  Tag( d_,  spec         ) );
      rhoYiTags_.push_back( Tag( rhoYi_, spec ) );
    }
    rhoYiTags_.pop_back();
  }

};

} // namespace pokitt
#endif /* EXAMPLES_REACTIONDIFFUSION_TAGMANAGER_H_ */
