// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Select.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/moleculegroup.h"

#include "SireMol/molecules.h"

#include "SireMol/parser.h"

#include "SireMol/select.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "select.h"

#include "select.h"

SireMol::Select __copy__(const SireMol::Select &other){ return SireMol::Select(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Select_class(){

    { //::SireMol::Select
        typedef bp::class_< SireMol::Select, bp::bases< SireBase::Property > > Select_exposer_t;
        Select_exposer_t Select_exposer = Select_exposer_t( "Select", "This is the only publicly visible selector class. This provides a\nfront-end interface to selecting atoms and molecules\n\nAuthor: Christopher Woods\n", bp::init< >("Construct an empty selection (will select nothing)") );
        bp::scope Select_scope( Select_exposer );
        Select_exposer.def( bp::init< QString const & >(( bp::arg("str") ), "Construct a selection based on the passed string") );
        Select_exposer.def( bp::init< SireMol::Select const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::Select::objectType
        
            typedef ::QString ( ::SireMol::Select::*objectType_function_type)(  ) const;
            objectType_function_type objectType_function_value( &::SireMol::Select::objectType );
            
            Select_exposer.def( 
                "objectType"
                , objectType_function_value
                , "" );
        
        }
        Select_exposer.def( bp::self != bp::self );
        { //::SireMol::Select::operator()
        
            typedef ::SireMol::SelectResult ( ::SireMol::Select::*__call___function_type)( ::SireMol::MolGroupsBase const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::Select::operator() );
            
            Select_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("molgroups"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Select::operator()
        
            typedef ::SireMol::SelectResult ( ::SireMol::Select::*__call___function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::Select::operator() );
            
            Select_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Select::operator()
        
            typedef ::SireMol::SelectResult ( ::SireMol::Select::*__call___function_type)( ::SireMol::Molecules const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::Select::operator() );
            
            Select_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("molecules"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Select::operator()
        
            typedef ::SireMol::SelectResult ( ::SireMol::Select::*__call___function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::Select::operator() );
            
            Select_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Select::operator()
        
            typedef ::SireMol::SelectResult ( ::SireMol::Select::*__call___function_type)( ::SireMol::SelectResult const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::Select::operator() );
            
            Select_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("result"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::Select::operator=
        
            typedef ::SireMol::Select & ( ::SireMol::Select::*assign_function_type)( ::SireMol::Select const & ) ;
            assign_function_type assign_function_value( &::SireMol::Select::operator= );
            
            Select_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Select_exposer.def( bp::self == bp::self );
        { //::SireMol::Select::resetTokens
        
            typedef void ( *resetTokens_function_type )(  );
            resetTokens_function_type resetTokens_function_value( &::SireMol::Select::resetTokens );
            
            Select_exposer.def( 
                "resetTokens"
                , resetTokens_function_value
                , "Clear all user-set tokens" );
        
        }
        { //::SireMol::Select::setToken
        
            typedef void ( *setToken_function_type )( ::QString const &,::QString const & );
            setToken_function_type setToken_function_value( &::SireMol::Select::setToken );
            
            Select_exposer.def( 
                "setToken"
                , setToken_function_value
                , ( bp::arg("token"), bp::arg("selection") )
                , "Set a user token that will be substituted for the passed selection, e.g.\nsetToken(protein, molecules with resname alai)\nwould allow you to use protein to refer to any molecules that contain\nresidues called alai\nNote that the token is set globally for all searches\n" );
        
        }
        { //::SireMol::Select::toString
        
            typedef ::QString ( ::SireMol::Select::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Select::toString );
            
            Select_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMol::Select::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Select::typeName );
            
            Select_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMol::Select::what
        
            typedef char const * ( ::SireMol::Select::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::Select::what );
            
            Select_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        Select_exposer.staticmethod( "resetTokens" );
        Select_exposer.staticmethod( "setToken" );
        Select_exposer.staticmethod( "typeName" );
        Select_exposer.def( "__copy__", &__copy__);
        Select_exposer.def( "__deepcopy__", &__copy__);
        Select_exposer.def( "clone", &__copy__);
        Select_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::Select >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Select_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::Select >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Select_exposer.def( "__str__", &__str__< ::SireMol::Select > );
        Select_exposer.def( "__repr__", &__str__< ::SireMol::Select > );
    }

}
