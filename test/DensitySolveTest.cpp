
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include <expression/ExprLib.h>
#include <expression/ExprPatch.h>
#include <expression/Expression.h>
#include <expression/ExpressionFactory.h>
#include <expression/FieldManagerList.h>
typedef Expr::ExprPatch  PatchT;

#include <pokitt/CanteraObjects.h>
#include <pokitt/SpeciesN.h>
#include <pokitt/thermo/Enthalpy.h>
#include <pokitt/thermo/Density.h>
#include <pokitt/MixtureMolWeight.h>

#include <spatialops/Nebo.h>
using SpatialOps::operator<<=;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <pokitt/DensitySolve/DensityCalculator_species.h>

typedef SpatialOps::SVolField FieldT;

using ConstantT       = Expr::ConstantExpr<FieldT>::Builder;
using LinearT         = Expr::LinearFunction<FieldT>::Builder;
using DensSolveT      = WasatchCore::DensityFromSpecies<FieldT>::Builder;
using EnthalpyT       = pokitt::Enthalpy<FieldT>::Builder;
using DensityT        = pokitt::Density<FieldT>::Builder;
using MixMWT          = pokitt::MixtureMolWeight<FieldT>::Builder;
using SpeciesNT       = pokitt::SpeciesN<FieldT>::Builder;
using PlaceHolderT    = Expr::PlaceHolder<FieldT>::Builder;

//==============================================================================
bool
setup_test()
{
  std::cout <<"calling setup_test...\n";

  std::string inputFileName = "h2o2.xml";
  const CanteraObjects::Setup setup( "Mix", inputFileName );
  CanteraObjects::setup_cantera( setup );

  const size_t nSpec = CanteraObjects::number_species(); std::cout<<"\nnSpec: " <<nSpec<<std::endl;
  const std::vector<std::string>& spNames = CanteraObjects::species_names();
  const std::vector<std::string>  spNamesNM1(spNames.begin(),spNames.end()-1);

  PatchT patch(1);
  Expr::ExpressionFactory factory;
//  Expr::ExpressionTree tree( 0, factory, patch.id(), "outer_tree" );
  Expr::ExpressionTree tree( factory, 0, "tree" );
  Expr::FieldManagerList& fml = patch.field_manager_list();

  const double oldDensity = 0.3;
  const double newDensity = 0.32;


  std::vector<double> yiOld( nSpec, 0.0 );
  std::vector<double> rhoYi( nSpec-1, 0.0 );
  for( int i=0; i<nSpec; ++i ){
    yiOld[i] = 1.         / static_cast<double>( nSpec );
    std::cout<<"\nsetting yi: " <<i<<std::endl;
    if( i < nSpec-1){
      std::cout<<"\n setting rhoYi: " <<i<<std::endl;
      rhoYi[i] = newDensity / static_cast<double>( nSpec );
    }
  }
  std::cout<<"\ncheckpoint 1\n";
  const double offset = 1e-6;
  const double oldTemp = 1300;
  const double newTemp = 1301;
  const double pressure = 101325;

  const double rTol = 1e-5;
  const int iter = 5;


  const Expr::TagList yiOldTags = Expr::tag_list( spNames, Expr::STATE_NONE, "Y_" );
  const Expr::TagList rhoYiTags = Expr::tag_list( spNamesNM1, Expr::STATE_NONE, "rhoY_" );

  const Expr::Tag pressTag  ("pressure",  Expr::STATE_NONE);
  const Expr::Tag oldDensTag("old_density",  Expr::STATE_NONE);
  const Expr::Tag newDensTag("new_density",  Expr::STATE_NONE);
  const Expr::Tag oldTempTag("old_temp",  Expr::STATE_NONE);
  const Expr::Tag newTempTag("new_temp",  Expr::STATE_NONE);
  const Expr::Tag oldEnthTag("old_enthalpy",  Expr::STATE_NONE);
  const Expr::Tag rhoEnthTag("rhoEnthalpy",  Expr::STATE_NONE);
  const Expr::Tag mixMWTag  ("mixMolWeight",  Expr::STATE_NONE);
  std::cout<<"\ncheckpoint 2\n";

  for( int i=0; i<nSpec-1; ++i ){
    std::cout<<"\nregistering yi and rhoYi: " <<i<<std::endl;
    tree.insert_tree( factory.register_expression( new ConstantT( yiOldTags[i], yiOld[i] ) ) );
    tree.insert_tree( factory.register_expression( new ConstantT( rhoYiTags[i], rhoYi[i] ) ) );
  }
  tree.insert_tree( factory.register_expression( new SpeciesNT( yiOldTags[nSpec-1], yiOldTags ) ) );

  tree.insert_tree( factory.register_expression( new ConstantT( pressTag, pressure ) ) );
  tree.insert_tree( factory.register_expression( new ConstantT( oldTempTag, oldTemp ) ) );
  tree.insert_tree( factory.register_expression( new EnthalpyT( oldEnthTag, oldTempTag, yiOldTags ) ) );
//  tree.insert_tree( factory.register_expression( new ConstantT( newDensTag, newDensity ), true ) );
  tree.insert_tree( factory.register_expression( new LinearT( rhoEnthTag, oldEnthTag, newDensity, 0 ) ) );
  tree.insert_tree( factory.register_expression( new MixMWT( mixMWTag, yiOldTags, pokitt::FractionType::MASS ) ) );

  tree.insert_tree( factory.register_expression( new DensityT( oldDensTag, oldTempTag, pressTag, mixMWTag ) ) );

  tree.insert_tree( factory.register_expression( new DensSolveT( newDensTag, oldDensTag, oldEnthTag, rhoEnthTag, oldTempTag, pressTag,  yiOldTags, rhoYiTags, rTol, iter ) ) );
//  tree.insert_tree( factory.register_expression( new DensityT( newDensTag, oldTempTag, pressTag, mixMWTag ) ) );
  std::cout<<"\ncheckpoint 3\n";

  tree.register_fields( fml );
  tree.lock_fields(fml);
  fml.allocate_fields( patch.field_info() );
  tree.bind_fields( fml );
  std::cout<<"\ncheckpoint 4\n";

//  const FieldT& output    = fml.field_manager<FieldT>().field_ref( outputTag );

  {
    std::ofstream ofile("outer_tree.dot");
    tree.write_tree( ofile );
  }

  const int nExec = 1;
  std::cout<<"\ncheckpoint 5\n";
  std::cout<<"\ncheckpoint 5.1\n";

  std::cout<<"\ncheckpoint 6\n";
    tree.execute_tree();
    std::cout<<"done \n";


  return true;
}

//==============================================================================

int main(){

  setup_test();

  return 0;
}
