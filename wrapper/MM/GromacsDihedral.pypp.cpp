// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "GromacsDihedral.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/conditional.h"

#include "SireCAS/exp.h"

#include "SireCAS/sum.h"

#include "SireCAS/trigfuncs.h"

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "amberparams.h"

#include "gromacsparams.h"

#include "gromacsparams.h"

SireMM::GromacsDihedral __copy__(const SireMM::GromacsDihedral &other){ return SireMM::GromacsDihedral(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_GromacsDihedral_class(){

    { //::SireMM::GromacsDihedral
        typedef bp::class_< SireMM::GromacsDihedral > GromacsDihedral_exposer_t;
        GromacsDihedral_exposer_t GromacsDihedral_exposer = GromacsDihedral_exposer_t( "GromacsDihedral", "This class holds all of the information about a Gromacs Dihedral\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope GromacsDihedral_scope( GromacsDihedral_exposer );
        GromacsDihedral_exposer.def( bp::init< int >(( bp::arg("function_type") ), "Construct a dihedral that is of the specified type, but the parameters have yet\nto be resolved. This is because Gromacs can indicate the required type of\nfunction in the molecule specification, without providing the parameters") );
        GromacsDihedral_exposer.def( bp::init< int, double, bp::optional< double, double, double, double, double > >(( bp::arg("function_type"), bp::arg("k0"), bp::arg("k1")=0, bp::arg("k2")=0, bp::arg("k3")=0, bp::arg("k4")=0, bp::arg("k5")=0 ), "Construct an dihedral of the specified function type with specified parameters\n(the order should be the same as in the Gromacs Manual, table 5.5)") );
        GromacsDihedral_exposer.def( bp::init< int, QList< double > const & >(( bp::arg("function_type"), bp::arg("params") ), "Construct a dihedral of the specified function type by interpreting the parameter\ndata from the passed list of parameter values. These should be in the\nsame order as in the Gromacs Manual, table 5.5") );
        GromacsDihedral_exposer.def( bp::init< SireCAS::Expression const &, SireCAS::Symbol const & >(( bp::arg("dihedral"), bp::arg("phi") ), "Construct from the passed dihedral, using phi as the symbol for the phi value") );
        GromacsDihedral_exposer.def( bp::init< SireMM::GromacsDihedral const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::GromacsDihedral::assertResolved
        
            typedef void ( ::SireMM::GromacsDihedral::*assertResolved_function_type)(  ) const;
            assertResolved_function_type assertResolved_function_value( &::SireMM::GromacsDihedral::assertResolved );
            
            GromacsDihedral_exposer.def( 
                "assertResolved"
                , assertResolved_function_value
                , "Assert that the parameters for this dihedral have been resolved" );
        
        }
        { //::SireMM::GromacsDihedral::at
        
            typedef double ( ::SireMM::GromacsDihedral::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMM::GromacsDihedral::at );
            
            GromacsDihedral_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , "Return the ith parameter for this dihedral" );
        
        }
        { //::SireMM::GromacsDihedral::construct
        
            typedef ::QList< SireMM::GromacsDihedral > ( *construct_function_type )( ::SireCAS::Expression const &,::SireCAS::Symbol const & );
            construct_function_type construct_function_value( &::SireMM::GromacsDihedral::construct );
            
            GromacsDihedral_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("dihedral"), bp::arg("phi") )
                , "Construct from the passed dihedral, using phi as the symbol for the phi value" );
        
        }
        { //::SireMM::GromacsDihedral::constructImproper
        
            typedef ::QList< SireMM::GromacsDihedral > ( *constructImproper_function_type )( ::SireCAS::Expression const &,::SireCAS::Symbol const & );
            constructImproper_function_type constructImproper_function_value( &::SireMM::GromacsDihedral::constructImproper );
            
            GromacsDihedral_exposer.def( 
                "constructImproper"
                , constructImproper_function_value
                , ( bp::arg("dihedral"), bp::arg("phi") )
                , "Construct from the passed improper, using phi as the symbol for the phi value" );
        
        }
        { //::SireMM::GromacsDihedral::count
        
            typedef int ( ::SireMM::GromacsDihedral::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMM::GromacsDihedral::count );
            
            GromacsDihedral_exposer.def( 
                "count"
                , count_function_value
                , "Return the number of parameters associated with this dihedral type" );
        
        }
        { //::SireMM::GromacsDihedral::functionType
        
            typedef int ( ::SireMM::GromacsDihedral::*functionType_function_type)(  ) const;
            functionType_function_type functionType_function_value( &::SireMM::GromacsDihedral::functionType );
            
            GromacsDihedral_exposer.def( 
                "functionType"
                , functionType_function_value
                , "Return the Gromacs ID number for the function type for this dihedral. See table\n5.5 in the Gromacs manual for information" );
        
        }
        { //::SireMM::GromacsDihedral::functionTypeString
        
            typedef ::QString ( ::SireMM::GromacsDihedral::*functionTypeString_function_type)(  ) const;
            functionTypeString_function_type functionTypeString_function_value( &::SireMM::GromacsDihedral::functionTypeString );
            
            GromacsDihedral_exposer.def( 
                "functionTypeString"
                , functionTypeString_function_value
                , "Return the string description of the function type for this dihedral" );
        
        }
        { //::SireMM::GromacsDihedral::hash
        
            typedef ::uint ( ::SireMM::GromacsDihedral::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMM::GromacsDihedral::hash );
            
            GromacsDihedral_exposer.def( 
                "hash"
                , hash_function_value
                , "Return a hash for this bond" );
        
        }
        { //::SireMM::GromacsDihedral::isAngleTorsionCrossTerm
        
            typedef bool ( ::SireMM::GromacsDihedral::*isAngleTorsionCrossTerm_function_type)(  ) const;
            isAngleTorsionCrossTerm_function_type isAngleTorsionCrossTerm_function_value( &::SireMM::GromacsDihedral::isAngleTorsionCrossTerm );
            
            GromacsDihedral_exposer.def( 
                "isAngleTorsionCrossTerm"
                , isAngleTorsionCrossTerm_function_value
                , "Return whether or not this dihedral is a angletorsion cross term" );
        
        }
        { //::SireMM::GromacsDihedral::isCosine
        
            typedef bool ( ::SireMM::GromacsDihedral::*isCosine_function_type)(  ) const;
            isCosine_function_type isCosine_function_value( &::SireMM::GromacsDihedral::isCosine );
            
            GromacsDihedral_exposer.def( 
                "isCosine"
                , isCosine_function_value
                , "Return whether or not this is a cosine-series dihedral" );
        
        }
        { //::SireMM::GromacsDihedral::isImproperAngleTerm
        
            typedef bool ( ::SireMM::GromacsDihedral::*isImproperAngleTerm_function_type)(  ) const;
            isImproperAngleTerm_function_type isImproperAngleTerm_function_value( &::SireMM::GromacsDihedral::isImproperAngleTerm );
            
            GromacsDihedral_exposer.def( 
                "isImproperAngleTerm"
                , isImproperAngleTerm_function_value
                , "Return whether or not this dihedral is really an improper angle term" );
        
        }
        { //::SireMM::GromacsDihedral::isResolved
        
            typedef bool ( ::SireMM::GromacsDihedral::*isResolved_function_type)(  ) const;
            isResolved_function_type isResolved_function_value( &::SireMM::GromacsDihedral::isResolved );
            
            GromacsDihedral_exposer.def( 
                "isResolved"
                , isResolved_function_value
                , "Return whether or not the parameters for this dihedral are resolved" );
        
        }
        { //::SireMM::GromacsDihedral::isSimple
        
            typedef bool ( ::SireMM::GromacsDihedral::*isSimple_function_type)(  ) const;
            isSimple_function_type isSimple_function_value( &::SireMM::GromacsDihedral::isSimple );
            
            GromacsDihedral_exposer.def( 
                "isSimple"
                , isSimple_function_value
                , "Return whether or not this is a simple dihedral function, based only on the\nsize of the torsion" );
        
        }
        { //::SireMM::GromacsDihedral::needsResolving
        
            typedef bool ( ::SireMM::GromacsDihedral::*needsResolving_function_type)(  ) const;
            needsResolving_function_type needsResolving_function_value( &::SireMM::GromacsDihedral::needsResolving );
            
            GromacsDihedral_exposer.def( 
                "needsResolving"
                , needsResolving_function_value
                , "Return whether or not this parameter needs resolving" );
        
        }
        GromacsDihedral_exposer.def( bp::self != bp::self );
        GromacsDihedral_exposer.def( bp::self < bp::self );
        GromacsDihedral_exposer.def( bp::self <= bp::self );
        { //::SireMM::GromacsDihedral::operator=
        
            typedef ::SireMM::GromacsDihedral & ( ::SireMM::GromacsDihedral::*assign_function_type)( ::SireMM::GromacsDihedral const & ) ;
            assign_function_type assign_function_value( &::SireMM::GromacsDihedral::operator= );
            
            GromacsDihedral_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        GromacsDihedral_exposer.def( bp::self == bp::self );
        GromacsDihedral_exposer.def( bp::self > bp::self );
        GromacsDihedral_exposer.def( bp::self >= bp::self );
        { //::SireMM::GromacsDihedral::operator[]
        
            typedef double ( ::SireMM::GromacsDihedral::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::GromacsDihedral::operator[] );
            
            GromacsDihedral_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMM::GromacsDihedral::parameters
        
            typedef ::QList< double > ( ::SireMM::GromacsDihedral::*parameters_function_type)(  ) const;
            parameters_function_type parameters_function_value( &::SireMM::GromacsDihedral::parameters );
            
            GromacsDihedral_exposer.def( 
                "parameters"
                , parameters_function_value
                , "Return all of the parameters for this dihedral" );
        
        }
        { //::SireMM::GromacsDihedral::size
        
            typedef int ( ::SireMM::GromacsDihedral::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMM::GromacsDihedral::size );
            
            GromacsDihedral_exposer.def( 
                "size"
                , size_function_value
                , "Return the number of parameters associated with this dihedral type" );
        
        }
        { //::SireMM::GromacsDihedral::toAngleTorsionExpression
        
            typedef ::SireCAS::Expression ( ::SireMM::GromacsDihedral::*toAngleTorsionExpression_function_type)( ::SireCAS::Symbol const &,::SireCAS::Symbol const &,::SireCAS::Symbol const & ) const;
            toAngleTorsionExpression_function_type toAngleTorsionExpression_function_value( &::SireMM::GromacsDihedral::toAngleTorsionExpression );
            
            GromacsDihedral_exposer.def( 
                "toAngleTorsionExpression"
                , toAngleTorsionExpression_function_value
                , ( bp::arg("theta0"), bp::arg("theta1"), bp::arg("phi") )
                , "Return this function converted to a SireCAS::Expression using the passed\nsymbol to represent the torsion (phi) and the angles either side of the\ntorsion (theta0 and theta1)" );
        
        }
        { //::SireMM::GromacsDihedral::toExpression
        
            typedef ::SireCAS::Expression ( ::SireMM::GromacsDihedral::*toExpression_function_type)( ::SireCAS::Symbol const & ) const;
            toExpression_function_type toExpression_function_value( &::SireMM::GromacsDihedral::toExpression );
            
            GromacsDihedral_exposer.def( 
                "toExpression"
                , toExpression_function_value
                , ( bp::arg("phi") )
                , "Return this function converted to a SireCAS::Expression using the passed symbol\nto represent the torsion size" );
        
        }
        { //::SireMM::GromacsDihedral::toImproperExpression
        
            typedef ::SireCAS::Expression ( ::SireMM::GromacsDihedral::*toImproperExpression_function_type)( ::SireCAS::Symbol const & ) const;
            toImproperExpression_function_type toImproperExpression_function_value( &::SireMM::GromacsDihedral::toImproperExpression );
            
            GromacsDihedral_exposer.def( 
                "toImproperExpression"
                , toImproperExpression_function_value
                , ( bp::arg("eta") )
                , "Return this function converted to a SireCAS::Expression using the passed\nsymbol to represent the improper angle eta" );
        
        }
        { //::SireMM::GromacsDihedral::toString
        
            typedef ::QString ( ::SireMM::GromacsDihedral::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::GromacsDihedral::toString );
            
            GromacsDihedral_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this dihedral" );
        
        }
        { //::SireMM::GromacsDihedral::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::GromacsDihedral::typeName );
            
            GromacsDihedral_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMM::GromacsDihedral::what
        
            typedef char const * ( ::SireMM::GromacsDihedral::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::GromacsDihedral::what );
            
            GromacsDihedral_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        GromacsDihedral_exposer.staticmethod( "construct" );
        GromacsDihedral_exposer.staticmethod( "constructImproper" );
        GromacsDihedral_exposer.staticmethod( "typeName" );
        GromacsDihedral_exposer.def( "__copy__", &__copy__);
        GromacsDihedral_exposer.def( "__deepcopy__", &__copy__);
        GromacsDihedral_exposer.def( "clone", &__copy__);
        GromacsDihedral_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::GromacsDihedral >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GromacsDihedral_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::GromacsDihedral >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GromacsDihedral_exposer.def( "__str__", &__str__< ::SireMM::GromacsDihedral > );
        GromacsDihedral_exposer.def( "__repr__", &__str__< ::SireMM::GromacsDihedral > );
        GromacsDihedral_exposer.def( "__len__", &__len_size< ::SireMM::GromacsDihedral > );
        GromacsDihedral_exposer.def( "__hash__", &::SireMM::GromacsDihedral::hash );
    }

}
