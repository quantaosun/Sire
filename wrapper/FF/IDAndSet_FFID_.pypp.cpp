// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "IDAndSet_FFID_.pypp.hpp"

namespace bp = boost::python;

#include "SireFF/errors.h"

#include "SireStream/datastream.h"

#include "ffid.h"

#include "ffidx.h"

#include "ffname.h"

#include "forcefields.h"

#include "ffid.h"

SireID::IDAndSet<SireFF::FFID> __copy__(const SireID::IDAndSet<SireFF::FFID> &other){ return SireID::IDAndSet<SireFF::FFID>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_IDAndSet_FFID__class(){

    { //::SireID::IDAndSet< SireFF::FFID >
        typedef bp::class_< SireID::IDAndSet< SireFF::FFID >, bp::bases< SireFF::FFID, SireID::ID > > IDAndSet_FFID__exposer_t;
        IDAndSet_FFID__exposer_t IDAndSet_FFID__exposer = IDAndSet_FFID__exposer_t( "IDAndSet_FFID_", bp::init< >() );
        bp::scope IDAndSet_FFID__scope( IDAndSet_FFID__exposer );
        IDAndSet_FFID__exposer.def( bp::init< SireFF::FFID const & >(( bp::arg("id") )) );
        IDAndSet_FFID__exposer.def( bp::init< SireFF::FFID const &, SireFF::FFID const & >(( bp::arg("id0"), bp::arg("id1") )) );
        IDAndSet_FFID__exposer.def( bp::init< QList< SireFF::FFIdentifier > const & >(( bp::arg("ids") )) );
        IDAndSet_FFID__exposer.def( bp::init< SireID::IDAndSet< SireFF::FFID > const & >(( bp::arg("ids") )) );
        IDAndSet_FFID__exposer.def( bp::init< SireID::IDAndSet< SireFF::FFID > const & >(( bp::arg("other") )) );
        { //::SireID::IDAndSet< SireFF::FFID >::IDs
        
            typedef SireID::IDAndSet< SireFF::FFID > exported_class_t;
            typedef ::QSet< SireFF::FFIdentifier > const & ( ::SireID::IDAndSet< SireFF::FFID >::*IDs_function_type)(  ) const;
            IDs_function_type IDs_function_value( &::SireID::IDAndSet< SireFF::FFID >::IDs );
            
            IDAndSet_FFID__exposer.def( 
                "IDs"
                , IDs_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireID::IDAndSet< SireFF::FFID >::hash
        
            typedef SireID::IDAndSet< SireFF::FFID > exported_class_t;
            typedef ::uint ( ::SireID::IDAndSet< SireFF::FFID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireID::IDAndSet< SireFF::FFID >::hash );
            
            IDAndSet_FFID__exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireID::IDAndSet< SireFF::FFID >::isNull
        
            typedef SireID::IDAndSet< SireFF::FFID > exported_class_t;
            typedef bool ( ::SireID::IDAndSet< SireFF::FFID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireID::IDAndSet< SireFF::FFID >::isNull );
            
            IDAndSet_FFID__exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireID::IDAndSet< SireFF::FFID >::map
        
            typedef SireID::IDAndSet< SireFF::FFID > exported_class_t;
            typedef ::QList< SireFF::FFIdx > ( ::SireID::IDAndSet< SireFF::FFID >::*map_function_type)( ::SireID::IDAndSet< SireFF::FFID >::SearchObject const & ) const;
            map_function_type map_function_value( &::SireID::IDAndSet< SireFF::FFID >::map );
            
            IDAndSet_FFID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("obj") ) );
        
        }
        IDAndSet_FFID__exposer.def( bp::self != bp::other< SireID::ID >() );
        IDAndSet_FFID__exposer.def( bp::self != bp::self );
        IDAndSet_FFID__exposer.def( bp::self != bp::other< SireFF::FFID >() );
        { //::SireID::IDAndSet< SireFF::FFID >::operator=
        
            typedef SireID::IDAndSet< SireFF::FFID > exported_class_t;
            typedef ::SireID::IDAndSet< SireFF::FFID > & ( ::SireID::IDAndSet< SireFF::FFID >::*assign_function_type)( ::SireID::IDAndSet< SireFF::FFID > const & ) ;
            assign_function_type assign_function_value( &::SireID::IDAndSet< SireFF::FFID >::operator= );
            
            IDAndSet_FFID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireID::IDAndSet< SireFF::FFID >::operator=
        
            typedef SireID::IDAndSet< SireFF::FFID > exported_class_t;
            typedef ::SireID::IDAndSet< SireFF::FFID > & ( ::SireID::IDAndSet< SireFF::FFID >::*assign_function_type)( ::SireFF::FFID const & ) ;
            assign_function_type assign_function_value( &::SireID::IDAndSet< SireFF::FFID >::operator= );
            
            IDAndSet_FFID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        IDAndSet_FFID__exposer.def( bp::self == bp::other< SireID::ID >() );
        IDAndSet_FFID__exposer.def( bp::self == bp::self );
        IDAndSet_FFID__exposer.def( bp::self == bp::other< SireFF::FFID >() );
        { //::SireID::IDAndSet< SireFF::FFID >::toString
        
            typedef SireID::IDAndSet< SireFF::FFID > exported_class_t;
            typedef ::QString ( ::SireID::IDAndSet< SireFF::FFID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireID::IDAndSet< SireFF::FFID >::toString );
            
            IDAndSet_FFID__exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireID::IDAndSet< SireFF::FFID >::typeName
        
            typedef SireID::IDAndSet< SireFF::FFID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireID::IDAndSet< SireFF::FFID >::typeName );
            
            IDAndSet_FFID__exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireID::IDAndSet< SireFF::FFID >::what
        
            typedef SireID::IDAndSet< SireFF::FFID > exported_class_t;
            typedef char const * ( ::SireID::IDAndSet< SireFF::FFID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireID::IDAndSet< SireFF::FFID >::what );
            
            IDAndSet_FFID__exposer.def( 
                "what"
                , what_function_value );
        
        }
        IDAndSet_FFID__exposer.staticmethod( "typeName" );
        IDAndSet_FFID__exposer.def( "__copy__", &__copy__);
        IDAndSet_FFID__exposer.def( "__deepcopy__", &__copy__);
        IDAndSet_FFID__exposer.def( "clone", &__copy__);
        IDAndSet_FFID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireID::IDAndSet<SireFF::FFID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDAndSet_FFID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireID::IDAndSet<SireFF::FFID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDAndSet_FFID__exposer.def( "__str__", &__str__< ::SireID::IDAndSet<SireFF::FFID> > );
        IDAndSet_FFID__exposer.def( "__repr__", &__str__< ::SireID::IDAndSet<SireFF::FFID> > );
        IDAndSet_FFID__exposer.def( "__hash__", &::SireID::IDAndSet<SireFF::FFID>::hash );
    }

}
