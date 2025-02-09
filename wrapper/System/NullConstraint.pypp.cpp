// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "NullConstraint.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/numberproperty.h"

#include "SireBase/propertylist.h"

#include "SireError/errors.h"

#include "SireMaths/maths.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireStream/streamdata.hpp"

#include "SireSystem/errors.h"

#include "constraint.h"

#include "delta.h"

#include "system.h"

#include <QDebug>

#include "constraint.h"

SireSystem::NullConstraint __copy__(const SireSystem::NullConstraint &other){ return SireSystem::NullConstraint(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_NullConstraint_class(){

    { //::SireSystem::NullConstraint
        typedef bp::class_< SireSystem::NullConstraint, bp::bases< SireSystem::Constraint, SireBase::Property > > NullConstraint_exposer_t;
        NullConstraint_exposer_t NullConstraint_exposer = NullConstraint_exposer_t( "NullConstraint", "The null constraint", bp::init< >("Null constructor") );
        bp::scope NullConstraint_scope( NullConstraint_exposer );
        NullConstraint_exposer.def( bp::init< SireSystem::NullConstraint const & >(( bp::arg("other") ), "Copy constructor") );
        NullConstraint_exposer.def( bp::self != bp::self );
        { //::SireSystem::NullConstraint::operator=
        
            typedef ::SireSystem::NullConstraint & ( ::SireSystem::NullConstraint::*assign_function_type)( ::SireSystem::NullConstraint const & ) ;
            assign_function_type assign_function_value( &::SireSystem::NullConstraint::operator= );
            
            NullConstraint_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        NullConstraint_exposer.def( bp::self == bp::self );
        { //::SireSystem::NullConstraint::toString
        
            typedef ::QString ( ::SireSystem::NullConstraint::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireSystem::NullConstraint::toString );
            
            NullConstraint_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation" );
        
        }
        { //::SireSystem::NullConstraint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::NullConstraint::typeName );
            
            NullConstraint_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        NullConstraint_exposer.staticmethod( "typeName" );
        NullConstraint_exposer.def( "__copy__", &__copy__);
        NullConstraint_exposer.def( "__deepcopy__", &__copy__);
        NullConstraint_exposer.def( "clone", &__copy__);
        NullConstraint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::NullConstraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullConstraint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::NullConstraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NullConstraint_exposer.def( "__str__", &__str__< ::SireSystem::NullConstraint > );
        NullConstraint_exposer.def( "__repr__", &__str__< ::SireSystem::NullConstraint > );
    }

}
