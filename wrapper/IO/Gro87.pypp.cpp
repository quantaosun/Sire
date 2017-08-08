// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Gro87.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireIO/errors.h"

#include "SireIO/gro87.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "gro87.h"

#include "gro87.h"

SireIO::Gro87 __copy__(const SireIO::Gro87 &other){ return SireIO::Gro87(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Gro87_class(){

    { //::SireIO::Gro87
        typedef bp::class_< SireIO::Gro87, bp::bases< SireIO::MoleculeParser, SireBase::Property > > Gro87_exposer_t;
        Gro87_exposer_t Gro87_exposer = Gro87_exposer_t( "Gro87", "This class holds a parser for reading and writing Gromacs Gro87 files.\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope Gro87_scope( Gro87_exposer );
        Gro87_exposer.def( bp::init< QString const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("filename"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to read in the data from the file called filename. The\npassed property map can be used to pass extra parameters to control\nthe parsing") );
        Gro87_exposer.def( bp::init< QStringList const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("lines"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to read in the data from the passed text lines. The\npassed property map can be used to pass extra parameters to control\nthe parsing") );
        Gro87_exposer.def( bp::init< SireSystem::System const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("system"), bp::arg("map")=SireBase::PropertyMap() ), "Construct this parser by extracting all necessary information from the\npassed SireSystem::System, looking for the properties that are specified\nin the passed property map") );
        Gro87_exposer.def( bp::init< SireIO::Gro87 const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireIO::Gro87::atomNames
        
            typedef ::QVector< QString > ( ::SireIO::Gro87::*atomNames_function_type)(  ) const;
            atomNames_function_type atomNames_function_value( &::SireIO::Gro87::atomNames );
            
            Gro87_exposer.def( 
                "atomNames"
                , atomNames_function_value
                , "Return the names of all of the atoms, in the same order as the coordinates" );
        
        }
        { //::SireIO::Gro87::atomNumbers
        
            typedef ::QVector< long long > ( ::SireIO::Gro87::*atomNumbers_function_type)(  ) const;
            atomNumbers_function_type atomNumbers_function_value( &::SireIO::Gro87::atomNumbers );
            
            Gro87_exposer.def( 
                "atomNumbers"
                , atomNumbers_function_value
                , "Return the numbers of all of the atoms. These are in the same order\nas the coordinates" );
        
        }
        { //::SireIO::Gro87::boxV1
        
            typedef ::SireMaths::Vector ( ::SireIO::Gro87::*boxV1_function_type)(  ) const;
            boxV1_function_type boxV1_function_value( &::SireIO::Gro87::boxV1 );
            
            Gro87_exposer.def( 
                "boxV1"
                , boxV1_function_value
                , "" );
        
        }
        { //::SireIO::Gro87::boxV1
        
            typedef ::SireMaths::Vector ( ::SireIO::Gro87::*boxV1_function_type)( int ) const;
            boxV1_function_type boxV1_function_value( &::SireIO::Gro87::boxV1 );
            
            Gro87_exposer.def( 
                "boxV1"
                , boxV1_function_value
                , ( bp::arg("frame") )
                , "" );
        
        }
        { //::SireIO::Gro87::boxV2
        
            typedef ::SireMaths::Vector ( ::SireIO::Gro87::*boxV2_function_type)(  ) const;
            boxV2_function_type boxV2_function_value( &::SireIO::Gro87::boxV2 );
            
            Gro87_exposer.def( 
                "boxV2"
                , boxV2_function_value
                , "" );
        
        }
        { //::SireIO::Gro87::boxV2
        
            typedef ::SireMaths::Vector ( ::SireIO::Gro87::*boxV2_function_type)( int ) const;
            boxV2_function_type boxV2_function_value( &::SireIO::Gro87::boxV2 );
            
            Gro87_exposer.def( 
                "boxV2"
                , boxV2_function_value
                , ( bp::arg("frame") )
                , "" );
        
        }
        { //::SireIO::Gro87::boxV3
        
            typedef ::SireMaths::Vector ( ::SireIO::Gro87::*boxV3_function_type)(  ) const;
            boxV3_function_type boxV3_function_value( &::SireIO::Gro87::boxV3 );
            
            Gro87_exposer.def( 
                "boxV3"
                , boxV3_function_value
                , "" );
        
        }
        { //::SireIO::Gro87::boxV3
        
            typedef ::SireMaths::Vector ( ::SireIO::Gro87::*boxV3_function_type)( int ) const;
            boxV3_function_type boxV3_function_value( &::SireIO::Gro87::boxV3 );
            
            Gro87_exposer.def( 
                "boxV3"
                , boxV3_function_value
                , ( bp::arg("frame") )
                , "" );
        
        }
        { //::SireIO::Gro87::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::Gro87::*construct_function_type)( ::QString const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::Gro87::construct );
            
            Gro87_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("filename"), bp::arg("map") )
                , "Return the parser that has been constructed by reading in the passed\nfile using the passed properties" );
        
        }
        { //::SireIO::Gro87::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::Gro87::*construct_function_type)( ::QStringList const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::Gro87::construct );
            
            Gro87_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("lines"), bp::arg("map") )
                , "Return the parser that has been constructed by reading in the passed\ntext lines using the passed properties" );
        
        }
        { //::SireIO::Gro87::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::Gro87::*construct_function_type)( ::SireSystem::System const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::Gro87::construct );
            
            Gro87_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("system"), bp::arg("map") )
                , "Return the parser that has been constructed by extract all necessary\ndata from the passed SireSystem::System using the specified properties" );
        
        }
        { //::SireIO::Gro87::coordinates
        
            typedef ::QVector< SireMaths::Vector > ( ::SireIO::Gro87::*coordinates_function_type)(  ) const;
            coordinates_function_type coordinates_function_value( &::SireIO::Gro87::coordinates );
            
            Gro87_exposer.def( 
                "coordinates"
                , coordinates_function_value
                , "Return the coordinates of the atoms for the first frame of the trajectory" );
        
        }
        { //::SireIO::Gro87::coordinates
        
            typedef ::QVector< SireMaths::Vector > ( ::SireIO::Gro87::*coordinates_function_type)( int ) const;
            coordinates_function_type coordinates_function_value( &::SireIO::Gro87::coordinates );
            
            Gro87_exposer.def( 
                "coordinates"
                , coordinates_function_value
                , ( bp::arg("frame") )
                , "Return the coordinates of the atoms at frame frame" );
        
        }
        { //::SireIO::Gro87::formatDescription
        
            typedef ::QString ( ::SireIO::Gro87::*formatDescription_function_type)(  ) const;
            formatDescription_function_type formatDescription_function_value( &::SireIO::Gro87::formatDescription );
            
            Gro87_exposer.def( 
                "formatDescription"
                , formatDescription_function_value
                , "Return a description of the file format" );
        
        }
        { //::SireIO::Gro87::formatName
        
            typedef ::QString ( ::SireIO::Gro87::*formatName_function_type)(  ) const;
            formatName_function_type formatName_function_value( &::SireIO::Gro87::formatName );
            
            Gro87_exposer.def( 
                "formatName"
                , formatName_function_value
                , "Return the format name that is used to identify this file format within Sire" );
        
        }
        { //::SireIO::Gro87::formatSuffix
        
            typedef ::QStringList ( ::SireIO::Gro87::*formatSuffix_function_type)(  ) const;
            formatSuffix_function_type formatSuffix_function_value( &::SireIO::Gro87::formatSuffix );
            
            Gro87_exposer.def( 
                "formatSuffix"
                , formatSuffix_function_value
                , "Return the suffixes that these files are normally associated with" );
        
        }
        { //::SireIO::Gro87::hasCoordinates
        
            typedef bool ( ::SireIO::Gro87::*hasCoordinates_function_type)(  ) const;
            hasCoordinates_function_type hasCoordinates_function_value( &::SireIO::Gro87::hasCoordinates );
            
            Gro87_exposer.def( 
                "hasCoordinates"
                , hasCoordinates_function_value
                , "Return whether or not this file contained coordinate data" );
        
        }
        { //::SireIO::Gro87::hasVelocities
        
            typedef bool ( ::SireIO::Gro87::*hasVelocities_function_type)(  ) const;
            hasVelocities_function_type hasVelocities_function_value( &::SireIO::Gro87::hasVelocities );
            
            Gro87_exposer.def( 
                "hasVelocities"
                , hasVelocities_function_value
                , "Return whether or not this file contained velocity data" );
        
        }
        { //::SireIO::Gro87::nAtoms
        
            typedef int ( ::SireIO::Gro87::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::Gro87::nAtoms );
            
            Gro87_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , "Return the number of atoms whose data is contained in this file" );
        
        }
        { //::SireIO::Gro87::nFrames
        
            typedef int ( ::SireIO::Gro87::*nFrames_function_type)(  ) const;
            nFrames_function_type nFrames_function_value( &::SireIO::Gro87::nFrames );
            
            Gro87_exposer.def( 
                "nFrames"
                , nFrames_function_value
                , "Return the number of frames of the trajectory loaded from the file" );
        
        }
        { //::SireIO::Gro87::nResidues
        
            typedef int ( ::SireIO::Gro87::*nResidues_function_type)(  ) const;
            nResidues_function_type nResidues_function_value( &::SireIO::Gro87::nResidues );
            
            Gro87_exposer.def( 
                "nResidues"
                , nResidues_function_value
                , "" );
        
        }
        Gro87_exposer.def( bp::self != bp::self );
        { //::SireIO::Gro87::operator=
        
            typedef ::SireIO::Gro87 & ( ::SireIO::Gro87::*assign_function_type)( ::SireIO::Gro87 const & ) ;
            assign_function_type assign_function_value( &::SireIO::Gro87::operator= );
            
            Gro87_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Gro87_exposer.def( bp::self == bp::self );
        { //::SireIO::Gro87::residueNames
        
            typedef ::QVector< QString > ( ::SireIO::Gro87::*residueNames_function_type)(  ) const;
            residueNames_function_type residueNames_function_value( &::SireIO::Gro87::residueNames );
            
            Gro87_exposer.def( 
                "residueNames"
                , residueNames_function_value
                , "Return the residue name for each atom (one per atom), in the same\norder as the coordinates" );
        
        }
        { //::SireIO::Gro87::residueNumbers
        
            typedef ::QVector< long long > ( ::SireIO::Gro87::*residueNumbers_function_type)(  ) const;
            residueNumbers_function_type residueNumbers_function_value( &::SireIO::Gro87::residueNumbers );
            
            Gro87_exposer.def( 
                "residueNumbers"
                , residueNumbers_function_value
                , "Return the residue number for each atom (one per atom), in the same\norder as the coordinates" );
        
        }
        { //::SireIO::Gro87::time
        
            typedef double ( ::SireIO::Gro87::*time_function_type)(  ) const;
            time_function_type time_function_value( &::SireIO::Gro87::time );
            
            Gro87_exposer.def( 
                "time"
                , time_function_value
                , "Return the current time of the simulation from which this coordinate\nfile was written. Returns 0 if there is no time set. If there are\nmultiple frames, then the time of the first frame is returned" );
        
        }
        { //::SireIO::Gro87::time
        
            typedef double ( ::SireIO::Gro87::*time_function_type)( int ) const;
            time_function_type time_function_value( &::SireIO::Gro87::time );
            
            Gro87_exposer.def( 
                "time"
                , time_function_value
                , ( bp::arg("frame") )
                , "Return the time for the structure at the specified frame" );
        
        }
        { //::SireIO::Gro87::title
        
            typedef ::QString ( ::SireIO::Gro87::*title_function_type)(  ) const;
            title_function_type title_function_value( &::SireIO::Gro87::title );
            
            Gro87_exposer.def( 
                "title"
                , title_function_value
                , "Return the title of the file" );
        
        }
        { //::SireIO::Gro87::toString
        
            typedef ::QString ( ::SireIO::Gro87::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireIO::Gro87::toString );
            
            Gro87_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this parser" );
        
        }
        { //::SireIO::Gro87::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::Gro87::typeName );
            
            Gro87_exposer.def( 
                "typeName"
                , typeName_function_value
                , "Return the C++ name for this class" );
        
        }
        { //::SireIO::Gro87::velocities
        
            typedef ::QVector< SireMaths::Vector > ( ::SireIO::Gro87::*velocities_function_type)(  ) const;
            velocities_function_type velocities_function_value( &::SireIO::Gro87::velocities );
            
            Gro87_exposer.def( 
                "velocities"
                , velocities_function_value
                , "Return the velocities of the atoms for the first frame of the trajectory" );
        
        }
        { //::SireIO::Gro87::velocities
        
            typedef ::QVector< SireMaths::Vector > ( ::SireIO::Gro87::*velocities_function_type)( int ) const;
            velocities_function_type velocities_function_value( &::SireIO::Gro87::velocities );
            
            Gro87_exposer.def( 
                "velocities"
                , velocities_function_value
                , ( bp::arg("frame") )
                , "Return the velocities of the atoms at frame frame" );
        
        }
        { //::SireIO::Gro87::warnings
        
            typedef ::QStringList ( ::SireIO::Gro87::*warnings_function_type)(  ) const;
            warnings_function_type warnings_function_value( &::SireIO::Gro87::warnings );
            
            Gro87_exposer.def( 
                "warnings"
                , warnings_function_value
                , "Return the warnings encountered when parsing the file. This\nis empty if everything was ok" );
        
        }
        { //::SireIO::Gro87::what
        
            typedef char const * ( ::SireIO::Gro87::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireIO::Gro87::what );
            
            Gro87_exposer.def( 
                "what"
                , what_function_value
                , "Return the C++ name for this class" );
        
        }
        Gro87_exposer.staticmethod( "typeName" );
        Gro87_exposer.def( "__copy__", &__copy__);
        Gro87_exposer.def( "__deepcopy__", &__copy__);
        Gro87_exposer.def( "clone", &__copy__);
        Gro87_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::Gro87 >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Gro87_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::Gro87 >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Gro87_exposer.def( "__str__", &__str__< ::SireIO::Gro87 > );
        Gro87_exposer.def( "__repr__", &__str__< ::SireIO::Gro87 > );
    }

}
