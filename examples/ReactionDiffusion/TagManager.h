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
using Expr::Tag;
using Expr::TagList;

namespace pokitt{

  enum TagName{
    RHOYI,
    H,
    YI,
    T,
    P,
    RHO,
    MMW,
    R,
    LAM,
    D,
    J,
    Q,
    QX,
    QY,
    XCOORD,
    YCOORD,
  };

class TagManager{
  const Tag rhoYi_;
  const Tag h_;

  const Tag yi_;
  const Tag t_;
  const Tag p_;
  const Tag rho_;
  const Tag mmw_;
  const Tag r_;
  const Tag lam_;
  const Tag d_;
  const Tag j_;
  const Tag q_;
  const Tag qx_;
  const Tag qy_;

  const Tag x_;
  const Tag y_;

  TagList rTags_;
  TagList dTags_;
  TagList yiTags_;
  TagList rhoYiTags_;
  TagList jxTags_;
  TagList jyTags_;
  TagList hiTags_;

public:

  const Tag& operator[](TagName name) const{
    switch( name ){
    case RHOYI:  return rhoYi_; break;
    case H:      return h_;     break;
    case YI:     return yi_;    break;
    case T:      return t_;     break;
    case P:      return p_;     break;
    case RHO:    return rho_;   break;
    case MMW:    return mmw_;   break;
    case R:      return r_;     break;
    case LAM:    return lam_;   break;
    case D:      return d_;     break;
    case J:      return j_;     break;
    case Q:      return q_;     break;
    case QX:     return qx_;    break;
    case QY:     return qy_;    break;
    case XCOORD: return x_;     break;
    case YCOORD: return y_;     break;
    }
  }

  const TagList& rN()              const { return rTags_;     }
  const TagList& dN()              const { return dTags_;     }
  const TagList& yiN()             const { return yiTags_;    }
  const TagList& rhoYiN()          const { return rhoYiTags_; }
  const TagList& jxN()             const { return jxTags_;    }
  const TagList& hiN()             const { return hiTags_;    }
  const TagList& jyN()             const { return jyTags_;    }

  const Tag& yiN   ( const int n ) const { return yiTags_[n];    }
  const Tag& rhoYiN( const int n ) const { return rhoYiTags_[n]; }
  const Tag& rN    ( const int n ) const { return rTags_[n];     }
  const Tag& jxN   ( const int n ) const { return jxTags_[n];    }
  const Tag& jyN   ( const int n ) const { return jyTags_[n];    }
  const Tag& hiN   ( const int n ) const { return hiTags_[n];    }

  TagManager( const int nSpec,
              const Tag& rhoYi,
              const Tag& h,
              const Tag& yi,
              const Tag& t,
              const Tag& p,
              const Tag& rho,
              const Tag& mmw,
              const Tag& r,
              const Tag& lam,
              const Tag& d,
              const Tag& j,
              const Tag& q,
              const Tag& x,
              const Tag& y ):
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
                qx_( Tag( q, "x" ) ),
                qy_( Tag( q, "y" ) )

{
    for( size_t n=0; n<nSpec; ++n ){
      std::string spec = boost::lexical_cast<std::string>(n);
      hiTags_.push_back( Tag( h_.name() + spec, Expr::STATE_NONE ) );
      yiTags_.push_back( Tag( yi_, spec         ) );
      jxTags_.push_back( Tag( j_,  "x" + spec   ) );
      jyTags_.push_back( Tag( j_,  "y" + spec   ) );
      rTags_.push_back(  Tag( r_,  spec         ) );
      dTags_.push_back(  Tag( d_,  spec         ) );
      if( n != (nSpec - 1) ){
        rhoYiTags_.push_back( Tag( rhoYi_, spec ) );
      }
    }
}
};

}
#endif /* EXAMPLES_REACTIONDIFFUSION_TAGMANAGER_H_ */
