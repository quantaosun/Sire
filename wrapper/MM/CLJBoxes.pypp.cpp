// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CLJBoxes.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMM/cljboxes.h"

#include "SireMM/cljdelta.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireVol/aabox.h"

#include "SireVol/periodicbox.h"

#include "SireVol/space.h"

#include "cljboxes.h"

#include "tostring.h"

#include <QDebug>

#include <QElapsedTimer>

#include "cljboxes.h"

SireMM::CLJBoxes __copy__(const SireMM::CLJBoxes &other){ return SireMM::CLJBoxes(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_CLJBoxes_class(){

    { //::SireMM::CLJBoxes
        typedef bp::class_< SireMM::CLJBoxes > CLJBoxes_exposer_t;
        CLJBoxes_exposer_t CLJBoxes_exposer = CLJBoxes_exposer_t( "CLJBoxes", bp::init< >() );
        bp::scope CLJBoxes_scope( CLJBoxes_exposer );
        CLJBoxes_exposer.def( bp::init< SireUnits::Dimension::Length >(( bp::arg("box_size") )) );
        CLJBoxes_exposer.def( bp::init< SireMM::CLJAtoms const & >(( bp::arg("atoms") )) );
        CLJBoxes_exposer.def( bp::init< SireMM::CLJAtoms const &, SireMM::CLJAtoms const & >(( bp::arg("atoms0"), bp::arg("atoms1") )) );
        CLJBoxes_exposer.def( bp::init< SireMM::CLJAtoms const &, SireUnits::Dimension::Length >(( bp::arg("atoms"), bp::arg("box_size") )) );
        CLJBoxes_exposer.def( bp::init< SireMM::CLJAtoms const &, SireMM::CLJAtoms const &, SireUnits::Dimension::Length >(( bp::arg("atoms0"), bp::arg("atoms1"), bp::arg("box_size") )) );
        CLJBoxes_exposer.def( bp::init< SireMM::CLJBoxes const & >(( bp::arg("other") )) );
        { //::SireMM::CLJBoxes::add
        
            typedef ::QVector< SireMM::CLJBoxIndex > ( ::SireMM::CLJBoxes::*add_function_type)( ::SireMM::CLJAtoms const & ) ;
            add_function_type add_function_value( &::SireMM::CLJBoxes::add );
            
            CLJBoxes_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("atoms") ) );
        
        }
        { //::SireMM::CLJBoxes::at
        
            typedef ::SireMM::CLJAtom ( ::SireMM::CLJBoxes::*at_function_type)( ::SireMM::CLJBoxIndex const & ) const;
            at_function_type at_function_value( &::SireMM::CLJBoxes::at );
            
            CLJBoxes_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("idx") ) );
        
        }
        { //::SireMM::CLJBoxes::atoms
        
            typedef ::SireMM::CLJAtoms ( ::SireMM::CLJBoxes::*atoms_function_type)(  ) const;
            atoms_function_type atoms_function_value( &::SireMM::CLJBoxes::atoms );
            
            CLJBoxes_exposer.def( 
                "atoms"
                , atoms_function_value );
        
        }
        { //::SireMM::CLJBoxes::atoms
        
            typedef ::SireMM::CLJAtoms ( ::SireMM::CLJBoxes::*atoms_function_type)( ::QVector< SireMM::CLJBoxIndex > const & ) const;
            atoms_function_type atoms_function_value( &::SireMM::CLJBoxes::atoms );
            
            CLJBoxes_exposer.def( 
                "atoms"
                , atoms_function_value
                , ( bp::arg("atoms") ) );
        
        }
        { //::SireMM::CLJBoxes::boxAt
        
            typedef ::SireMM::CLJBox ( ::SireMM::CLJBoxes::*boxAt_function_type)( int ) const;
            boxAt_function_type boxAt_function_value( &::SireMM::CLJBoxes::boxAt );
            
            CLJBoxes_exposer.def( 
                "boxAt"
                , boxAt_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMM::CLJBoxes::boxAt
        
            typedef ::SireMM::CLJBox ( ::SireMM::CLJBoxes::*boxAt_function_type)( ::SireMM::CLJBoxIndex const & ) const;
            boxAt_function_type boxAt_function_value( &::SireMM::CLJBoxes::boxAt );
            
            CLJBoxes_exposer.def( 
                "boxAt"
                , boxAt_function_value
                , ( bp::arg("index") ) );
        
        }
        { //::SireMM::CLJBoxes::boxAt
        
            typedef ::SireMM::CLJBox ( ::SireMM::CLJBoxes::*boxAt_function_type)( ::SireMaths::Vector const & ) const;
            boxAt_function_type boxAt_function_value( &::SireMM::CLJBoxes::boxAt );
            
            CLJBoxes_exposer.def( 
                "boxAt"
                , boxAt_function_value
                , ( bp::arg("coords") ) );
        
        }
        { //::SireMM::CLJBoxes::boxDimensions
        
            typedef ::QVector< SireVol::AABox > ( ::SireMM::CLJBoxes::*boxDimensions_function_type)(  ) const;
            boxDimensions_function_type boxDimensions_function_value( &::SireMM::CLJBoxes::boxDimensions );
            
            CLJBoxes_exposer.def( 
                "boxDimensions"
                , boxDimensions_function_value );
        
        }
        { //::SireMM::CLJBoxes::boxDimensionsAt
        
            typedef ::SireVol::AABox ( ::SireMM::CLJBoxes::*boxDimensionsAt_function_type)( int ) const;
            boxDimensionsAt_function_type boxDimensionsAt_function_value( &::SireMM::CLJBoxes::boxDimensionsAt );
            
            CLJBoxes_exposer.def( 
                "boxDimensionsAt"
                , boxDimensionsAt_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMM::CLJBoxes::boxDimensionsAt
        
            typedef ::SireVol::AABox ( ::SireMM::CLJBoxes::*boxDimensionsAt_function_type)( ::SireMM::CLJBoxIndex const & ) const;
            boxDimensionsAt_function_type boxDimensionsAt_function_value( &::SireMM::CLJBoxes::boxDimensionsAt );
            
            CLJBoxes_exposer.def( 
                "boxDimensionsAt"
                , boxDimensionsAt_function_value
                , ( bp::arg("index") ) );
        
        }
        { //::SireMM::CLJBoxes::boxDimensionsAt
        
            typedef ::SireVol::AABox ( ::SireMM::CLJBoxes::*boxDimensionsAt_function_type)( ::SireMaths::Vector const & ) const;
            boxDimensionsAt_function_type boxDimensionsAt_function_value( &::SireMM::CLJBoxes::boxDimensionsAt );
            
            CLJBoxes_exposer.def( 
                "boxDimensionsAt"
                , boxDimensionsAt_function_value
                , ( bp::arg("coords") ) );
        
        }
        { //::SireMM::CLJBoxes::boxes
        
            typedef ::QVector< SireMM::CLJBox > ( ::SireMM::CLJBoxes::*boxes_function_type)(  ) const;
            boxes_function_type boxes_function_value( &::SireMM::CLJBoxes::boxes );
            
            CLJBoxes_exposer.def( 
                "boxes"
                , boxes_function_value );
        
        }
        { //::SireMM::CLJBoxes::get
        
            typedef ::SireMM::CLJAtoms ( ::SireMM::CLJBoxes::*get_function_type)( ::QVector< SireMM::CLJBoxIndex > const & ) const;
            get_function_type get_function_value( &::SireMM::CLJBoxes::get );
            
            CLJBoxes_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("atoms") ) );
        
        }
        { //::SireMM::CLJBoxes::getDistance
        
            typedef float ( ::SireMM::CLJBoxes::*getDistance_function_type)( ::SireMM::CLJBoxIndex const &,::SireMM::CLJBoxIndex const & ) const;
            getDistance_function_type getDistance_function_value( &::SireMM::CLJBoxes::getDistance );
            
            CLJBoxes_exposer.def( 
                "getDistance"
                , getDistance_function_value
                , ( bp::arg("box0"), bp::arg("box1") ) );
        
        }
        { //::SireMM::CLJBoxes::getDistance
        
            typedef float ( ::SireMM::CLJBoxes::*getDistance_function_type)( ::SireVol::Space const &,::SireMM::CLJBoxIndex const &,::SireMM::CLJBoxIndex const & ) const;
            getDistance_function_type getDistance_function_value( &::SireMM::CLJBoxes::getDistance );
            
            CLJBoxes_exposer.def( 
                "getDistance"
                , getDistance_function_value
                , ( bp::arg("space"), bp::arg("box0"), bp::arg("box1") ) );
        
        }
        { //::SireMM::CLJBoxes::getDistance
        
            typedef float ( ::SireMM::CLJBoxes::*getDistance_function_type)( ::SireVol::Space const &,::SireMM::CLJBoxIndex const &,::SireMM::CLJBoxIndex const &,::quint32,::quint32,::quint32 ) const;
            getDistance_function_type getDistance_function_value( &::SireMM::CLJBoxes::getDistance );
            
            CLJBoxes_exposer.def( 
                "getDistance"
                , getDistance_function_value
                , ( bp::arg("space"), bp::arg("box0"), bp::arg("box1"), bp::arg("nx"), bp::arg("ny"), bp::arg("nz") ) );
        
        }
        { //::SireMM::CLJBoxes::getDistances
        
            typedef ::QVector< SireMM::CLJBoxDistance > ( *getDistances_function_type )( ::SireVol::Space const &,::SireMM::CLJBoxes const & );
            getDistances_function_type getDistances_function_value( &::SireMM::CLJBoxes::getDistances );
            
            CLJBoxes_exposer.def( 
                "getDistances"
                , getDistances_function_value
                , ( bp::arg("space"), bp::arg("boxes") ) );
        
        }
        { //::SireMM::CLJBoxes::getDistances
        
            typedef ::QVector< SireMM::CLJBoxDistance > ( *getDistances_function_type )( ::SireVol::Space const &,::SireMM::CLJBoxes const &,::SireUnits::Dimension::Length );
            getDistances_function_type getDistances_function_value( &::SireMM::CLJBoxes::getDistances );
            
            CLJBoxes_exposer.def( 
                "getDistances"
                , getDistances_function_value
                , ( bp::arg("space"), bp::arg("boxes"), bp::arg("cutoff") ) );
        
        }
        { //::SireMM::CLJBoxes::getDistances
        
            typedef ::QVector< SireMM::CLJBoxDistance > ( *getDistances_function_type )( ::SireVol::Space const &,::SireMM::CLJBoxes const &,::SireMM::CLJBoxes const & );
            getDistances_function_type getDistances_function_value( &::SireMM::CLJBoxes::getDistances );
            
            CLJBoxes_exposer.def( 
                "getDistances"
                , getDistances_function_value
                , ( bp::arg("space"), bp::arg("boxes0"), bp::arg("boxes1") ) );
        
        }
        { //::SireMM::CLJBoxes::getDistances
        
            typedef ::QVector< SireMM::CLJBoxDistance > ( *getDistances_function_type )( ::SireVol::Space const &,::SireMM::CLJBoxes const &,::SireMM::CLJBoxes const &,::SireUnits::Dimension::Length );
            getDistances_function_type getDistances_function_value( &::SireMM::CLJBoxes::getDistances );
            
            CLJBoxes_exposer.def( 
                "getDistances"
                , getDistances_function_value
                , ( bp::arg("space"), bp::arg("boxes0"), bp::arg("boxes1"), bp::arg("cutoff") ) );
        
        }
        { //::SireMM::CLJBoxes::getDistances
        
            typedef ::QVector< SireMM::CLJBoxDistance > ( *getDistances_function_type )( ::SireVol::Space const &,::SireMM::CLJAtoms const &,::SireMM::CLJBoxes const & );
            getDistances_function_type getDistances_function_value( &::SireMM::CLJBoxes::getDistances );
            
            CLJBoxes_exposer.def( 
                "getDistances"
                , getDistances_function_value
                , ( bp::arg("space"), bp::arg("atoms0"), bp::arg("boxes1") ) );
        
        }
        { //::SireMM::CLJBoxes::getDistances
        
            typedef ::QVector< SireMM::CLJBoxDistance > ( *getDistances_function_type )( ::SireVol::Space const &,::SireMM::CLJAtoms const &,::SireMM::CLJBoxes const &,::SireUnits::Dimension::Length );
            getDistances_function_type getDistances_function_value( &::SireMM::CLJBoxes::getDistances );
            
            CLJBoxes_exposer.def( 
                "getDistances"
                , getDistances_function_value
                , ( bp::arg("space"), bp::arg("atoms0"), bp::arg("boxes1"), bp::arg("cutoff") ) );
        
        }
        { //::SireMM::CLJBoxes::getitem
        
            typedef ::SireMM::CLJAtom ( ::SireMM::CLJBoxes::*getitem_function_type)( ::SireMM::CLJBoxIndex const & ) const;
            getitem_function_type getitem_function_value( &::SireMM::CLJBoxes::getitem );
            
            CLJBoxes_exposer.def( 
                "getitem"
                , getitem_function_value
                , ( bp::arg("idx") ) );
        
        }
        { //::SireMM::CLJBoxes::isEmpty
        
            typedef bool ( ::SireMM::CLJBoxes::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::CLJBoxes::isEmpty );
            
            CLJBoxes_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMM::CLJBoxes::length
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::CLJBoxes::*length_function_type)(  ) const;
            length_function_type length_function_value( &::SireMM::CLJBoxes::length );
            
            CLJBoxes_exposer.def( 
                "length"
                , length_function_value );
        
        }
        { //::SireMM::CLJBoxes::nAtoms
        
            typedef int ( ::SireMM::CLJBoxes::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMM::CLJBoxes::nAtoms );
            
            CLJBoxes_exposer.def( 
                "nAtoms"
                , nAtoms_function_value );
        
        }
        { //::SireMM::CLJBoxes::nOccupiedBoxes
        
            typedef int ( ::SireMM::CLJBoxes::*nOccupiedBoxes_function_type)(  ) const;
            nOccupiedBoxes_function_type nOccupiedBoxes_function_value( &::SireMM::CLJBoxes::nOccupiedBoxes );
            
            CLJBoxes_exposer.def( 
                "nOccupiedBoxes"
                , nOccupiedBoxes_function_value );
        
        }
        { //::SireMM::CLJBoxes::occupiedBoxIndicies
        
            typedef ::QVector< SireMM::CLJBoxIndex > ( ::SireMM::CLJBoxes::*occupiedBoxIndicies_function_type)(  ) const;
            occupiedBoxIndicies_function_type occupiedBoxIndicies_function_value( &::SireMM::CLJBoxes::occupiedBoxIndicies );
            
            CLJBoxes_exposer.def( 
                "occupiedBoxIndicies"
                , occupiedBoxIndicies_function_value );
        
        }
        { //::SireMM::CLJBoxes::occupiedBoxes
        
            typedef ::SireMM::CLJBoxes::Container const & ( ::SireMM::CLJBoxes::*occupiedBoxes_function_type)(  ) const;
            occupiedBoxes_function_type occupiedBoxes_function_value( &::SireMM::CLJBoxes::occupiedBoxes );
            
            CLJBoxes_exposer.def( 
                "occupiedBoxes"
                , occupiedBoxes_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        CLJBoxes_exposer.def( bp::self != bp::self );
        CLJBoxes_exposer.def( bp::self + bp::self );
        { //::SireMM::CLJBoxes::operator=
        
            typedef ::SireMM::CLJBoxes & ( ::SireMM::CLJBoxes::*assign_function_type)( ::SireMM::CLJBoxes const & ) ;
            assign_function_type assign_function_value( &::SireMM::CLJBoxes::operator= );
            
            CLJBoxes_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        CLJBoxes_exposer.def( bp::self == bp::self );
        { //::SireMM::CLJBoxes::operator[]
        
            typedef ::SireMM::CLJAtom ( ::SireMM::CLJBoxes::*__getitem___function_type)( ::SireMM::CLJBoxIndex const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::CLJBoxes::operator[] );
            
            CLJBoxes_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idx") ) );
        
        }
        { //::SireMM::CLJBoxes::remove
        
            typedef void ( ::SireMM::CLJBoxes::*remove_function_type)( ::QVector< SireMM::CLJBoxIndex > const & ) ;
            remove_function_type remove_function_value( &::SireMM::CLJBoxes::remove );
            
            CLJBoxes_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("atoms") ) );
        
        }
        { //::SireMM::CLJBoxes::squeeze
        
            typedef ::SireMM::CLJBoxes ( ::SireMM::CLJBoxes::*squeeze_function_type)(  ) const;
            squeeze_function_type squeeze_function_value( &::SireMM::CLJBoxes::squeeze );
            
            CLJBoxes_exposer.def( 
                "squeeze"
                , squeeze_function_value );
        
        }
        { //::SireMM::CLJBoxes::take
        
            typedef ::SireMM::CLJAtoms ( ::SireMM::CLJBoxes::*take_function_type)( ::QVector< SireMM::CLJBoxIndex > const & ) ;
            take_function_type take_function_value( &::SireMM::CLJBoxes::take );
            
            CLJBoxes_exposer.def( 
                "take"
                , take_function_value
                , ( bp::arg("atoms") ) );
        
        }
        { //::SireMM::CLJBoxes::toString
        
            typedef ::QString ( ::SireMM::CLJBoxes::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::CLJBoxes::toString );
            
            CLJBoxes_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMM::CLJBoxes::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::CLJBoxes::typeName );
            
            CLJBoxes_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMM::CLJBoxes::what
        
            typedef char const * ( ::SireMM::CLJBoxes::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::CLJBoxes::what );
            
            CLJBoxes_exposer.def( 
                "what"
                , what_function_value );
        
        }
        CLJBoxes_exposer.staticmethod( "getDistances" );
        CLJBoxes_exposer.staticmethod( "typeName" );
        CLJBoxes_exposer.def( "__copy__", &__copy__);
        CLJBoxes_exposer.def( "__deepcopy__", &__copy__);
        CLJBoxes_exposer.def( "clone", &__copy__);
        CLJBoxes_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CLJBoxes >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJBoxes_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CLJBoxes >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJBoxes_exposer.def( "__str__", &__str__< ::SireMM::CLJBoxes > );
        CLJBoxes_exposer.def( "__repr__", &__str__< ::SireMM::CLJBoxes > );
        CLJBoxes_exposer.def( "__getitem__", &::SireMM::CLJBoxes::getitem );
    }

}
