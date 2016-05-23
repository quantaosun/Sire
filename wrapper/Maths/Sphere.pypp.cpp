// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Sphere.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "rangenerator.h"

#include "sphere.h"

#include <QDebug>

#include <QElapsedTimer>

#include "sphere.h"

SireMaths::Sphere __copy__(const SireMaths::Sphere &other){ return SireMaths::Sphere(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Sphere_class(){

    { //::SireMaths::Sphere
        typedef bp::class_< SireMaths::Sphere > Sphere_exposer_t;
        Sphere_exposer_t Sphere_exposer = Sphere_exposer_t( "Sphere", bp::init< >() );
        bp::scope Sphere_scope( Sphere_exposer );
        Sphere_exposer.def( bp::init< double const & >(( bp::arg("radius") )) );
        Sphere_exposer.def( bp::init< SireMaths::Vector const &, double const & >(( bp::arg("position"), bp::arg("radius") )) );
        Sphere_exposer.def( bp::init< SireMaths::Sphere const & >(( bp::arg("other") )) );
        { //::SireMaths::Sphere::center
        
            typedef ::SireMaths::Vector const & ( ::SireMaths::Sphere::*center_function_type)(  ) const;
            center_function_type center_function_value( &::SireMaths::Sphere::center );
            
            Sphere_exposer.def( 
                "center"
                , center_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMaths::Sphere::combinedVolume
        
            typedef double ( *combinedVolume_function_type )( ::QVector< SireMaths::Sphere > const & );
            combinedVolume_function_type combinedVolume_function_value( &::SireMaths::Sphere::combinedVolume );
            
            Sphere_exposer.def( 
                "combinedVolume"
                , combinedVolume_function_value
                , ( bp::arg("spheres") ) );
        
        }
        { //::SireMaths::Sphere::combinedVolumeMC
        
            typedef double ( *combinedVolumeMC_function_type )( ::QVector< SireMaths::Sphere > const &,::qint64 );
            combinedVolumeMC_function_type combinedVolumeMC_function_value( &::SireMaths::Sphere::combinedVolumeMC );
            
            Sphere_exposer.def( 
                "combinedVolumeMC"
                , combinedVolumeMC_function_value
                , ( bp::arg("spheres"), bp::arg("nsamples")=(::qint64)(-1) ) );
        
        }
        { //::SireMaths::Sphere::contains
        
            typedef bool ( ::SireMaths::Sphere::*contains_function_type)( ::SireMaths::Vector const & ) const;
            contains_function_type contains_function_value( &::SireMaths::Sphere::contains );
            
            Sphere_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("point") ) );
        
        }
        { //::SireMaths::Sphere::contains
        
            typedef bool ( ::SireMaths::Sphere::*contains_function_type)( ::SireMaths::Sphere const & ) const;
            contains_function_type contains_function_value( &::SireMaths::Sphere::contains );
            
            Sphere_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::Sphere::intersectionVolume
        
            typedef double ( ::SireMaths::Sphere::*intersectionVolume_function_type)( ::SireMaths::Sphere const & ) const;
            intersectionVolume_function_type intersectionVolume_function_value( &::SireMaths::Sphere::intersectionVolume );
            
            Sphere_exposer.def( 
                "intersectionVolume"
                , intersectionVolume_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::Sphere::intersectionVolume
        
            typedef double ( ::SireMaths::Sphere::*intersectionVolume_function_type)( ::SireMaths::Sphere const &,::SireMaths::Sphere const & ) const;
            intersectionVolume_function_type intersectionVolume_function_value( &::SireMaths::Sphere::intersectionVolume );
            
            Sphere_exposer.def( 
                "intersectionVolume"
                , intersectionVolume_function_value
                , ( bp::arg("other0"), bp::arg("other1") ) );
        
        }
        { //::SireMaths::Sphere::intersects
        
            typedef bool ( ::SireMaths::Sphere::*intersects_function_type)( ::SireMaths::Sphere const & ) const;
            intersects_function_type intersects_function_value( &::SireMaths::Sphere::intersects );
            
            Sphere_exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("other") ) );
        
        }
        Sphere_exposer.def( bp::self != bp::self );
        Sphere_exposer.def( bp::self == bp::self );
        { //::SireMaths::Sphere::position
        
            typedef ::SireMaths::Vector const & ( ::SireMaths::Sphere::*position_function_type)(  ) const;
            position_function_type position_function_value( &::SireMaths::Sphere::position );
            
            Sphere_exposer.def( 
                "position"
                , position_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMaths::Sphere::radius
        
            typedef double ( ::SireMaths::Sphere::*radius_function_type)(  ) const;
            radius_function_type radius_function_value( &::SireMaths::Sphere::radius );
            
            Sphere_exposer.def( 
                "radius"
                , radius_function_value );
        
        }
        { //::SireMaths::Sphere::setCenter
        
            typedef void ( ::SireMaths::Sphere::*setCenter_function_type)( ::SireMaths::Vector const & ) ;
            setCenter_function_type setCenter_function_value( &::SireMaths::Sphere::setCenter );
            
            Sphere_exposer.def( 
                "setCenter"
                , setCenter_function_value
                , ( bp::arg("center") ) );
        
        }
        { //::SireMaths::Sphere::setPosition
        
            typedef void ( ::SireMaths::Sphere::*setPosition_function_type)( ::SireMaths::Vector const & ) ;
            setPosition_function_type setPosition_function_value( &::SireMaths::Sphere::setPosition );
            
            Sphere_exposer.def( 
                "setPosition"
                , setPosition_function_value
                , ( bp::arg("position") ) );
        
        }
        { //::SireMaths::Sphere::setRadius
        
            typedef void ( ::SireMaths::Sphere::*setRadius_function_type)( double ) ;
            setRadius_function_type setRadius_function_value( &::SireMaths::Sphere::setRadius );
            
            Sphere_exposer.def( 
                "setRadius"
                , setRadius_function_value
                , ( bp::arg("radius") ) );
        
        }
        { //::SireMaths::Sphere::surfaceArea
        
            typedef double ( ::SireMaths::Sphere::*surfaceArea_function_type)(  ) const;
            surfaceArea_function_type surfaceArea_function_value( &::SireMaths::Sphere::surfaceArea );
            
            Sphere_exposer.def( 
                "surfaceArea"
                , surfaceArea_function_value );
        
        }
        { //::SireMaths::Sphere::toString
        
            typedef ::QString ( ::SireMaths::Sphere::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMaths::Sphere::toString );
            
            Sphere_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMaths::Sphere::translate
        
            typedef ::SireMaths::Sphere ( ::SireMaths::Sphere::*translate_function_type)( ::SireMaths::Vector const & ) const;
            translate_function_type translate_function_value( &::SireMaths::Sphere::translate );
            
            Sphere_exposer.def( 
                "translate"
                , translate_function_value
                , ( bp::arg("delta") ) );
        
        }
        { //::SireMaths::Sphere::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::Sphere::typeName );
            
            Sphere_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMaths::Sphere::volume
        
            typedef double ( ::SireMaths::Sphere::*volume_function_type)(  ) const;
            volume_function_type volume_function_value( &::SireMaths::Sphere::volume );
            
            Sphere_exposer.def( 
                "volume"
                , volume_function_value );
        
        }
        { //::SireMaths::Sphere::what
        
            typedef char const * ( ::SireMaths::Sphere::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMaths::Sphere::what );
            
            Sphere_exposer.def( 
                "what"
                , what_function_value );
        
        }
        Sphere_exposer.staticmethod( "combinedVolume" );
        Sphere_exposer.staticmethod( "combinedVolumeMC" );
        Sphere_exposer.staticmethod( "typeName" );
        Sphere_exposer.def( "__copy__", &__copy__);
        Sphere_exposer.def( "__deepcopy__", &__copy__);
        Sphere_exposer.def( "clone", &__copy__);
        Sphere_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::Sphere >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Sphere_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::Sphere >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Sphere_exposer.def( "__str__", &__str__< ::SireMaths::Sphere > );
        Sphere_exposer.def( "__repr__", &__str__< ::SireMaths::Sphere > );
    }

}
