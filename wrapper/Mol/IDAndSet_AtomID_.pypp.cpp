// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "IDAndSet_AtomID_.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "atomid.h"

#include "atomidentifier.h"

#include "chain.h"

#include "cutgroup.h"

#include "editor.hpp"

#include "groupatomids.h"

#include "molatomid.h"

#include "moleculegroup.h"

#include "moleculegroups.h"

#include "molecules.h"

#include "molinfo.h"

#include "mover.hpp"

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include "withatoms.h"

#include <QDebug>

#include "atomid.h"

SireID::IDAndSet<SireMol::AtomID> __copy__(const SireID::IDAndSet<SireMol::AtomID> &other){ return SireID::IDAndSet<SireMol::AtomID>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_IDAndSet_AtomID__class(){

    { //::SireID::IDAndSet< SireMol::AtomID >
        typedef bp::class_< SireID::IDAndSet< SireMol::AtomID >, bp::bases< SireMol::AtomID, SireID::ID > > IDAndSet_AtomID__exposer_t;
        IDAndSet_AtomID__exposer_t IDAndSet_AtomID__exposer = IDAndSet_AtomID__exposer_t( "IDAndSet_AtomID_", bp::init< >() );
        bp::scope IDAndSet_AtomID__scope( IDAndSet_AtomID__exposer );
        IDAndSet_AtomID__exposer.def( bp::init< SireMol::AtomID const & >(( bp::arg("id") )) );
        IDAndSet_AtomID__exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("id0"), bp::arg("id1") )) );
        IDAndSet_AtomID__exposer.def( bp::init< QList< SireMol::AtomIdentifier > const & >(( bp::arg("ids") )) );
        IDAndSet_AtomID__exposer.def( bp::init< SireID::IDAndSet< SireMol::AtomID > const & >(( bp::arg("ids") )) );
        IDAndSet_AtomID__exposer.def( bp::init< SireID::IDAndSet< SireMol::AtomID > const & >(( bp::arg("other") )) );
        { //::SireID::IDAndSet< SireMol::AtomID >::IDs
        
            typedef SireID::IDAndSet< SireMol::AtomID > exported_class_t;
            typedef ::QSet< SireMol::AtomIdentifier > const & ( ::SireID::IDAndSet< SireMol::AtomID >::*IDs_function_type)(  ) const;
            IDs_function_type IDs_function_value( &::SireID::IDAndSet< SireMol::AtomID >::IDs );
            
            IDAndSet_AtomID__exposer.def( 
                "IDs"
                , IDs_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireID::IDAndSet< SireMol::AtomID >::hash
        
            typedef SireID::IDAndSet< SireMol::AtomID > exported_class_t;
            typedef ::uint ( ::SireID::IDAndSet< SireMol::AtomID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireID::IDAndSet< SireMol::AtomID >::hash );
            
            IDAndSet_AtomID__exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireID::IDAndSet< SireMol::AtomID >::isNull
        
            typedef SireID::IDAndSet< SireMol::AtomID > exported_class_t;
            typedef bool ( ::SireID::IDAndSet< SireMol::AtomID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireID::IDAndSet< SireMol::AtomID >::isNull );
            
            IDAndSet_AtomID__exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireID::IDAndSet< SireMol::AtomID >::map
        
            typedef SireID::IDAndSet< SireMol::AtomID > exported_class_t;
            typedef ::QList< SireMol::AtomIdx > ( ::SireID::IDAndSet< SireMol::AtomID >::*map_function_type)( ::SireID::IDAndSet< SireMol::AtomID >::SearchObject const & ) const;
            map_function_type map_function_value( &::SireID::IDAndSet< SireMol::AtomID >::map );
            
            IDAndSet_AtomID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("obj") ) );
        
        }
        IDAndSet_AtomID__exposer.def( bp::self != bp::other< SireID::ID >() );
        IDAndSet_AtomID__exposer.def( bp::self != bp::self );
        IDAndSet_AtomID__exposer.def( bp::self != bp::other< SireMol::AtomID >() );
        { //::SireID::IDAndSet< SireMol::AtomID >::operator=
        
            typedef SireID::IDAndSet< SireMol::AtomID > exported_class_t;
            typedef ::SireID::IDAndSet< SireMol::AtomID > & ( ::SireID::IDAndSet< SireMol::AtomID >::*assign_function_type)( ::SireID::IDAndSet< SireMol::AtomID > const & ) ;
            assign_function_type assign_function_value( &::SireID::IDAndSet< SireMol::AtomID >::operator= );
            
            IDAndSet_AtomID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireID::IDAndSet< SireMol::AtomID >::operator=
        
            typedef SireID::IDAndSet< SireMol::AtomID > exported_class_t;
            typedef ::SireID::IDAndSet< SireMol::AtomID > & ( ::SireID::IDAndSet< SireMol::AtomID >::*assign_function_type)( ::SireMol::AtomID const & ) ;
            assign_function_type assign_function_value( &::SireID::IDAndSet< SireMol::AtomID >::operator= );
            
            IDAndSet_AtomID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        IDAndSet_AtomID__exposer.def( bp::self == bp::other< SireID::ID >() );
        IDAndSet_AtomID__exposer.def( bp::self == bp::self );
        IDAndSet_AtomID__exposer.def( bp::self == bp::other< SireMol::AtomID >() );
        { //::SireID::IDAndSet< SireMol::AtomID >::toString
        
            typedef SireID::IDAndSet< SireMol::AtomID > exported_class_t;
            typedef ::QString ( ::SireID::IDAndSet< SireMol::AtomID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireID::IDAndSet< SireMol::AtomID >::toString );
            
            IDAndSet_AtomID__exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireID::IDAndSet< SireMol::AtomID >::typeName
        
            typedef SireID::IDAndSet< SireMol::AtomID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireID::IDAndSet< SireMol::AtomID >::typeName );
            
            IDAndSet_AtomID__exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireID::IDAndSet< SireMol::AtomID >::what
        
            typedef SireID::IDAndSet< SireMol::AtomID > exported_class_t;
            typedef char const * ( ::SireID::IDAndSet< SireMol::AtomID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireID::IDAndSet< SireMol::AtomID >::what );
            
            IDAndSet_AtomID__exposer.def( 
                "what"
                , what_function_value );
        
        }
        IDAndSet_AtomID__exposer.staticmethod( "typeName" );
        IDAndSet_AtomID__exposer.def( "__copy__", &__copy__);
        IDAndSet_AtomID__exposer.def( "__deepcopy__", &__copy__);
        IDAndSet_AtomID__exposer.def( "clone", &__copy__);
        IDAndSet_AtomID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireID::IDAndSet<SireMol::AtomID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDAndSet_AtomID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireID::IDAndSet<SireMol::AtomID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDAndSet_AtomID__exposer.def( "__str__", &__str__< ::SireID::IDAndSet<SireMol::AtomID> > );
        IDAndSet_AtomID__exposer.def( "__repr__", &__str__< ::SireID::IDAndSet<SireMol::AtomID> > );
        IDAndSet_AtomID__exposer.def( "__hash__", &::SireID::IDAndSet<SireMol::AtomID>::hash );
    }

}
