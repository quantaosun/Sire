// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "RanGenerator.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/refcountdata.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "rangenerator.h"

#include "vector.h"

#include <QDebug>

#include <QMutex>

#include <QUuid>

#include <QVector>

#include <boost/noncopyable.hpp>

#include <boost/scoped_array.hpp>

#include <limits>

#include "rangenerator.h"

SireMaths::RanGenerator __copy__(const SireMaths::RanGenerator &other){ return SireMaths::RanGenerator(other); }

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireMaths::RanGenerator&){ return "SireMaths::RanGenerator";}

void register_RanGenerator_class(){

    { //::SireMaths::RanGenerator
        typedef bp::class_< SireMaths::RanGenerator > RanGenerator_exposer_t;
        RanGenerator_exposer_t RanGenerator_exposer = RanGenerator_exposer_t( "RanGenerator", "This class provides a thread-safe, copyable and streamable\nrandom number generator. Copies are guaranteed to produce\ndifferent random number sequences (thus the possibility\nof accidental repeat random numbers is removed).\n\nAuthor: Christopher Woods\n", bp::init< >("Create a randomly seeded generator\n(actually a copy of the global, random generator)") );
        bp::scope RanGenerator_scope( RanGenerator_exposer );
        RanGenerator_exposer.def( bp::init< quint32 >(( bp::arg("seed") ), "Create a generator seeded with seed") );
        RanGenerator_exposer.def( bp::init< QVector< unsigned int > const & >(( bp::arg("seed") ), "Create a generator seeded with seed") );
        RanGenerator_exposer.def( bp::init< SireMaths::RanGenerator const & >(( bp::arg("other") ), "Copy constructor - this takes an explicitly shared\ncopy of other (this is to prevent repeat random numbers\nfrom being generated by implicit copies)") );
        { //::SireMaths::RanGenerator::detach
        
            typedef void ( ::SireMaths::RanGenerator::*detach_function_type)(  ) ;
            detach_function_type detach_function_value( &::SireMaths::RanGenerator::detach );
            
            RanGenerator_exposer.def( 
                "detach"
                , detach_function_value
                , "Detach from shared storage" );
        
        }
        { //::SireMaths::RanGenerator::getState
        
            typedef ::QVector< unsigned int > ( ::SireMaths::RanGenerator::*getState_function_type)(  ) const;
            getState_function_type getState_function_value( &::SireMaths::RanGenerator::getState );
            
            RanGenerator_exposer.def( 
                "getState"
                , getState_function_value
                , "Return the current state of the random number generator.\nUse this if you truly wish to get reproducible sequences\nof random numbers" );
        
        }
        { //::SireMaths::RanGenerator::global
        
            typedef ::SireMaths::RanGenerator const & ( *global_function_type )(  );
            global_function_type global_function_value( &::SireMaths::RanGenerator::global );
            
            RanGenerator_exposer.def( 
                "global"
                , global_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return a reference to the global random number generator\n(shared between all threads)" );
        
        }
        { //::SireMaths::RanGenerator::lock
        
            typedef void ( ::SireMaths::RanGenerator::*lock_function_type)(  ) const;
            lock_function_type lock_function_value( &::SireMaths::RanGenerator::lock );
            
            RanGenerator_exposer.def( 
                "lock"
                , lock_function_value
                , "Take hold of the generator lock. Only you can now generate\nrandom numbers while this lock is held" );
        
        }
        { //::SireMaths::RanGenerator::locked_rand
        
            typedef double ( ::SireMaths::RanGenerator::*locked_rand_function_type)(  ) const;
            locked_rand_function_type locked_rand_function_value( &::SireMaths::RanGenerator::locked_rand );
            
            RanGenerator_exposer.def( 
                "locked_rand"
                , locked_rand_function_value
                , "Return a random real number on [0,1]. Should only be called while\nyou hold the generator lock" );
        
        }
        { //::SireMaths::RanGenerator::locked_rand
        
            typedef double ( ::SireMaths::RanGenerator::*locked_rand_function_type)( double ) const;
            locked_rand_function_type locked_rand_function_value( &::SireMaths::RanGenerator::locked_rand );
            
            RanGenerator_exposer.def( 
                "locked_rand"
                , locked_rand_function_value
                , ( bp::arg("maxval") )
                , "Return a random real number on [0,maxval]. Should only be called while\nyou hold the generator lock" );
        
        }
        { //::SireMaths::RanGenerator::locked_rand
        
            typedef double ( ::SireMaths::RanGenerator::*locked_rand_function_type)( double,double ) const;
            locked_rand_function_type locked_rand_function_value( &::SireMaths::RanGenerator::locked_rand );
            
            RanGenerator_exposer.def( 
                "locked_rand"
                , locked_rand_function_value
                , ( bp::arg("minval"), bp::arg("maxval") )
                , "Return a random real number on [minval,maxval]. Should only be called while\nyou hold the generator lock" );
        
        }
        { //::SireMaths::RanGenerator::locked_randNorm
        
            typedef double ( ::SireMaths::RanGenerator::*locked_randNorm_function_type)(  ) const;
            locked_randNorm_function_type locked_randNorm_function_value( &::SireMaths::RanGenerator::locked_randNorm );
            
            RanGenerator_exposer.def( 
                "locked_randNorm"
                , locked_randNorm_function_value
                , "Return a random number generated from the normal distribution\nwith mean 0 and standard deviation 1. You must hold the generator\nlock when calling this function" );
        
        }
        { //::SireMaths::RanGenerator::locked_randNorm
        
            typedef double ( ::SireMaths::RanGenerator::*locked_randNorm_function_type)( double,double ) const;
            locked_randNorm_function_type locked_randNorm_function_value( &::SireMaths::RanGenerator::locked_randNorm );
            
            RanGenerator_exposer.def( 
                "locked_randNorm"
                , locked_randNorm_function_value
                , ( bp::arg("minval"), bp::arg("maxval") )
                , "Return a random number from the normal distribution\nwith supplied mean and variance. You must hold the generator\nlock when calling this function" );
        
        }
        { //::SireMaths::RanGenerator::locked_vectorOnSphere
        
            typedef ::SireMaths::Vector ( ::SireMaths::RanGenerator::*locked_vectorOnSphere_function_type)(  ) const;
            locked_vectorOnSphere_function_type locked_vectorOnSphere_function_value( &::SireMaths::RanGenerator::locked_vectorOnSphere );
            
            RanGenerator_exposer.def( 
                "locked_vectorOnSphere"
                , locked_vectorOnSphere_function_value
                , "Return a random vector on the unit sphere. You must hold the generator\nlock when calling this function" );
        
        }
        { //::SireMaths::RanGenerator::locked_vectorOnSphere
        
            typedef ::SireMaths::Vector ( ::SireMaths::RanGenerator::*locked_vectorOnSphere_function_type)( double ) const;
            locked_vectorOnSphere_function_type locked_vectorOnSphere_function_value( &::SireMaths::RanGenerator::locked_vectorOnSphere );
            
            RanGenerator_exposer.def( 
                "locked_vectorOnSphere"
                , locked_vectorOnSphere_function_value
                , ( bp::arg("radius") )
                , "Return a random vector on the sphere with radius radius.\nYou must hold the generator lock when calling this function" );
        
        }
        { //::SireMaths::RanGenerator::nrand
        
            typedef ::QVector< double > ( ::SireMaths::RanGenerator::*nrand_function_type)( int ) const;
            nrand_function_type nrand_function_value( &::SireMaths::RanGenerator::nrand );
            
            RanGenerator_exposer.def( 
                "nrand"
                , nrand_function_value
                , ( bp::arg("n") )
                , "Return an array of n random numbers on [0,1]" );
        
        }
        { //::SireMaths::RanGenerator::nrand
        
            typedef ::QVector< double > ( ::SireMaths::RanGenerator::*nrand_function_type)( int,double ) const;
            nrand_function_type nrand_function_value( &::SireMaths::RanGenerator::nrand );
            
            RanGenerator_exposer.def( 
                "nrand"
                , nrand_function_value
                , ( bp::arg("n"), bp::arg("maxval") )
                , "Return an array of n random numbers on [0,maxval]" );
        
        }
        { //::SireMaths::RanGenerator::nrand
        
            typedef ::QVector< double > ( ::SireMaths::RanGenerator::*nrand_function_type)( int,double,double ) const;
            nrand_function_type nrand_function_value( &::SireMaths::RanGenerator::nrand );
            
            RanGenerator_exposer.def( 
                "nrand"
                , nrand_function_value
                , ( bp::arg("n"), bp::arg("minval"), bp::arg("maxval") )
                , "Return an array of n random numbers on [minval,maxval]" );
        
        }
        { //::SireMaths::RanGenerator::nrand
        
            typedef void ( ::SireMaths::RanGenerator::*nrand_function_type)( ::QVector< double > & ) const;
            nrand_function_type nrand_function_value( &::SireMaths::RanGenerator::nrand );
            
            RanGenerator_exposer.def( 
                "nrand"
                , nrand_function_value
                , ( bp::arg("result") )
                , "Fill the passed array of doubles with random numbers. This replaces each\nvalue in the array with a random number on [0,1]" );
        
        }
        { //::SireMaths::RanGenerator::nrand
        
            typedef void ( ::SireMaths::RanGenerator::*nrand_function_type)( ::QVector< double > &,double ) const;
            nrand_function_type nrand_function_value( &::SireMaths::RanGenerator::nrand );
            
            RanGenerator_exposer.def( 
                "nrand"
                , nrand_function_value
                , ( bp::arg("result"), bp::arg("maxval") )
                , "Fill the passed array of doubles with random numbers. This replaces each\nvalue in the array with a random number on [0,maxval]" );
        
        }
        { //::SireMaths::RanGenerator::nrand
        
            typedef void ( ::SireMaths::RanGenerator::*nrand_function_type)( ::QVector< double > &,double,double ) const;
            nrand_function_type nrand_function_value( &::SireMaths::RanGenerator::nrand );
            
            RanGenerator_exposer.def( 
                "nrand"
                , nrand_function_value
                , ( bp::arg("result"), bp::arg("minval"), bp::arg("maxval") )
                , "Fill the passed array of doubles with random numbers. This replaces each\nvalue in the array with a random number on [minval,maxval]" );
        
        }
        { //::SireMaths::RanGenerator::nrandNorm
        
            typedef void ( ::SireMaths::RanGenerator::*nrandNorm_function_type)( ::QVector< double > &,double,double ) const;
            nrandNorm_function_type nrandNorm_function_value( &::SireMaths::RanGenerator::nrandNorm );
            
            RanGenerator_exposer.def( 
                "nrandNorm"
                , nrandNorm_function_value
                , ( bp::arg("result"), bp::arg("mean"), bp::arg("variance") )
                , "Fill the passed array with random numbers drawn from the normal\ndistribution with supplied mean and variance" );
        
        }
        { //::SireMaths::RanGenerator::nrandNorm
        
            typedef ::QVector< double > ( ::SireMaths::RanGenerator::*nrandNorm_function_type)( int,double,double ) const;
            nrandNorm_function_type nrandNorm_function_value( &::SireMaths::RanGenerator::nrandNorm );
            
            RanGenerator_exposer.def( 
                "nrandNorm"
                , nrandNorm_function_value
                , ( bp::arg("n"), bp::arg("mean"), bp::arg("variance") )
                , "Return an array of N random numbers drawn from the normal distribution with\nsupplied mean and variance" );
        
        }
        { //::SireMaths::RanGenerator::nvectorOnSphere
        
            typedef void ( ::SireMaths::RanGenerator::*nvectorOnSphere_function_type)( ::QVector< SireMaths::Vector > & ) const;
            nvectorOnSphere_function_type nvectorOnSphere_function_value( &::SireMaths::RanGenerator::nvectorOnSphere );
            
            RanGenerator_exposer.def( 
                "nvectorOnSphere"
                , nvectorOnSphere_function_value
                , ( bp::arg("result") )
                , "Fill the passed array with random vectors on a unit sphere" );
        
        }
        { //::SireMaths::RanGenerator::nvectorOnSphere
        
            typedef void ( ::SireMaths::RanGenerator::*nvectorOnSphere_function_type)( ::QVector< SireMaths::Vector > &,double ) const;
            nvectorOnSphere_function_type nvectorOnSphere_function_value( &::SireMaths::RanGenerator::nvectorOnSphere );
            
            RanGenerator_exposer.def( 
                "nvectorOnSphere"
                , nvectorOnSphere_function_value
                , ( bp::arg("result"), bp::arg("radius") )
                , "Fill the passed array with random vectors on a sphere with radius radius" );
        
        }
        { //::SireMaths::RanGenerator::nvectorOnSphere
        
            typedef ::QVector< SireMaths::Vector > ( ::SireMaths::RanGenerator::*nvectorOnSphere_function_type)( int ) const;
            nvectorOnSphere_function_type nvectorOnSphere_function_value( &::SireMaths::RanGenerator::nvectorOnSphere );
            
            RanGenerator_exposer.def( 
                "nvectorOnSphere"
                , nvectorOnSphere_function_value
                , ( bp::arg("n") )
                , "Return an array of n random vectors on a unit sphere" );
        
        }
        { //::SireMaths::RanGenerator::nvectorOnSphere
        
            typedef ::QVector< SireMaths::Vector > ( ::SireMaths::RanGenerator::*nvectorOnSphere_function_type)( int,double ) const;
            nvectorOnSphere_function_type nvectorOnSphere_function_value( &::SireMaths::RanGenerator::nvectorOnSphere );
            
            RanGenerator_exposer.def( 
                "nvectorOnSphere"
                , nvectorOnSphere_function_value
                , ( bp::arg("n"), bp::arg("radius") )
                , "Return an array of n random vectors on a sphere of radius radius" );
        
        }
        RanGenerator_exposer.def( bp::self != bp::self );
        { //::SireMaths::RanGenerator::operator=
        
            typedef ::SireMaths::RanGenerator & ( ::SireMaths::RanGenerator::*assign_function_type)( ::SireMaths::RanGenerator const & ) ;
            assign_function_type assign_function_value( &::SireMaths::RanGenerator::operator= );
            
            RanGenerator_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        RanGenerator_exposer.def( bp::self == bp::self );
        { //::SireMaths::RanGenerator::rand
        
            typedef double ( ::SireMaths::RanGenerator::*rand_function_type)(  ) const;
            rand_function_type rand_function_value( &::SireMaths::RanGenerator::rand );
            
            RanGenerator_exposer.def( 
                "rand"
                , rand_function_value
                , "Return a random real number on [0,1]" );
        
        }
        { //::SireMaths::RanGenerator::rand
        
            typedef double ( ::SireMaths::RanGenerator::*rand_function_type)( double ) const;
            rand_function_type rand_function_value( &::SireMaths::RanGenerator::rand );
            
            RanGenerator_exposer.def( 
                "rand"
                , rand_function_value
                , ( bp::arg("maxval") )
                , "Return a random real number on [0,maxval]" );
        
        }
        { //::SireMaths::RanGenerator::rand
        
            typedef double ( ::SireMaths::RanGenerator::*rand_function_type)( double,double ) const;
            rand_function_type rand_function_value( &::SireMaths::RanGenerator::rand );
            
            RanGenerator_exposer.def( 
                "rand"
                , rand_function_value
                , ( bp::arg("minval"), bp::arg("maxval") )
                , "Return a random real number on [minval,maxval]" );
        
        }
        { //::SireMaths::RanGenerator::rand53
        
            typedef double ( ::SireMaths::RanGenerator::*rand53_function_type)(  ) const;
            rand53_function_type rand53_function_value( &::SireMaths::RanGenerator::rand53 );
            
            RanGenerator_exposer.def( 
                "rand53"
                , rand53_function_value
                , "Return a high-precision random real number on [0,1)" );
        
        }
        { //::SireMaths::RanGenerator::rand53
        
            typedef double ( ::SireMaths::RanGenerator::*rand53_function_type)( double ) const;
            rand53_function_type rand53_function_value( &::SireMaths::RanGenerator::rand53 );
            
            RanGenerator_exposer.def( 
                "rand53"
                , rand53_function_value
                , ( bp::arg("maxval") )
                , "Return a high-precision random real number on [0,1)" );
        
        }
        { //::SireMaths::RanGenerator::rand53
        
            typedef double ( ::SireMaths::RanGenerator::*rand53_function_type)( double,double ) const;
            rand53_function_type rand53_function_value( &::SireMaths::RanGenerator::rand53 );
            
            RanGenerator_exposer.def( 
                "rand53"
                , rand53_function_value
                , ( bp::arg("minval"), bp::arg("maxval") )
                , "Return a high-precision random real number on [minval,maxval)" );
        
        }
        { //::SireMaths::RanGenerator::randBool
        
            typedef bool ( ::SireMaths::RanGenerator::*randBool_function_type)(  ) const;
            randBool_function_type randBool_function_value( &::SireMaths::RanGenerator::randBool );
            
            RanGenerator_exposer.def( 
                "randBool"
                , randBool_function_value
                , "Return a random true or false value" );
        
        }
        { //::SireMaths::RanGenerator::randInt
        
            typedef ::quint32 ( ::SireMaths::RanGenerator::*randInt_function_type)(  ) const;
            randInt_function_type randInt_function_value( &::SireMaths::RanGenerator::randInt );
            
            RanGenerator_exposer.def( 
                "randInt"
                , randInt_function_value
                , "Return a random 32bit unsigned integer in [0,2^32 - 1]" );
        
        }
        { //::SireMaths::RanGenerator::randInt
        
            typedef ::quint32 ( ::SireMaths::RanGenerator::*randInt_function_type)( ::quint32 ) const;
            randInt_function_type randInt_function_value( &::SireMaths::RanGenerator::randInt );
            
            RanGenerator_exposer.def( 
                "randInt"
                , randInt_function_value
                , ( bp::arg("maxval") )
                , "Return a random 32bit unsigned integer in [0,maxval]" );
        
        }
        { //::SireMaths::RanGenerator::randInt
        
            typedef ::qint32 ( ::SireMaths::RanGenerator::*randInt_function_type)( ::qint32,::qint32 ) const;
            randInt_function_type randInt_function_value( &::SireMaths::RanGenerator::randInt );
            
            RanGenerator_exposer.def( 
                "randInt"
                , randInt_function_value
                , ( bp::arg("minval"), bp::arg("maxval") )
                , "Return a random 32bit integer in [minval,maxval]" );
        
        }
        { //::SireMaths::RanGenerator::randInt64
        
            typedef ::quint64 ( ::SireMaths::RanGenerator::*randInt64_function_type)(  ) const;
            randInt64_function_type randInt64_function_value( &::SireMaths::RanGenerator::randInt64 );
            
            RanGenerator_exposer.def( 
                "randInt64"
                , randInt64_function_value
                , "Return a random 64bit unsigned integer on [0,2^64 - 1]" );
        
        }
        { //::SireMaths::RanGenerator::randInt64
        
            typedef ::quint64 ( ::SireMaths::RanGenerator::*randInt64_function_type)( ::quint64 ) const;
            randInt64_function_type randInt64_function_value( &::SireMaths::RanGenerator::randInt64 );
            
            RanGenerator_exposer.def( 
                "randInt64"
                , randInt64_function_value
                , ( bp::arg("maxval") )
                , "Return a random 64bit unsigned integer on [0,maxval]" );
        
        }
        { //::SireMaths::RanGenerator::randInt64
        
            typedef ::qint64 ( ::SireMaths::RanGenerator::*randInt64_function_type)( ::qint64,::qint64 ) const;
            randInt64_function_type randInt64_function_value( &::SireMaths::RanGenerator::randInt64 );
            
            RanGenerator_exposer.def( 
                "randInt64"
                , randInt64_function_value
                , ( bp::arg("minval"), bp::arg("maxval") )
                , "Return a random 64bit integer on [minval,maxval]" );
        
        }
        { //::SireMaths::RanGenerator::randNorm
        
            typedef double ( ::SireMaths::RanGenerator::*randNorm_function_type)(  ) const;
            randNorm_function_type randNorm_function_value( &::SireMaths::RanGenerator::randNorm );
            
            RanGenerator_exposer.def( 
                "randNorm"
                , randNorm_function_value
                , "Return a random number generated from the normal distribution\nwith mean 0 and standard deviation 1" );
        
        }
        { //::SireMaths::RanGenerator::randNorm
        
            typedef double ( ::SireMaths::RanGenerator::*randNorm_function_type)( double,double ) const;
            randNorm_function_type randNorm_function_value( &::SireMaths::RanGenerator::randNorm );
            
            RanGenerator_exposer.def( 
                "randNorm"
                , randNorm_function_value
                , ( bp::arg("mean"), bp::arg("variance") )
                , "Return a random number from the normal distribution\nwith supplied mean and variance." );
        
        }
        { //::SireMaths::RanGenerator::seed
        
            typedef void ( ::SireMaths::RanGenerator::*seed_function_type)(  ) ;
            seed_function_type seed_function_value( &::SireMaths::RanGenerator::seed );
            
            RanGenerator_exposer.def( 
                "seed"
                , seed_function_value
                , "See the generator with a new, random seed - this will detach\nthis explicitly shared copy of the generator" );
        
        }
        { //::SireMaths::RanGenerator::seed
        
            typedef void ( ::SireMaths::RanGenerator::*seed_function_type)( ::quint32 ) ;
            seed_function_type seed_function_value( &::SireMaths::RanGenerator::seed );
            
            RanGenerator_exposer.def( 
                "seed"
                , seed_function_value
                , ( bp::arg("seed") )
                , "Seed the generator with s  - this will detach\nthis explicitly shared copy of the generator" );
        
        }
        { //::SireMaths::RanGenerator::seed
        
            typedef void ( ::SireMaths::RanGenerator::*seed_function_type)( ::QVector< unsigned int > const & ) ;
            seed_function_type seed_function_value( &::SireMaths::RanGenerator::seed );
            
            RanGenerator_exposer.def( 
                "seed"
                , seed_function_value
                , ( bp::arg("seed") )
                , "Seed the generator with seed - this will detach\nthis explicitly shared copy of the generator" );
        
        }
        { //::SireMaths::RanGenerator::seed
        
            typedef void ( ::SireMaths::RanGenerator::*seed_function_type)( ::SireMaths::RanGenerator const & ) ;
            seed_function_type seed_function_value( &::SireMaths::RanGenerator::seed );
            
            RanGenerator_exposer.def( 
                "seed"
                , seed_function_value
                , ( bp::arg("other") )
                , "Seed the generator with another generator - this\nreally just copies the generator as they are\nall explicit copies of one another" );
        
        }
        { //::SireMaths::RanGenerator::seedGlobal
        
            typedef void ( *seedGlobal_function_type )(  );
            seedGlobal_function_type seedGlobal_function_value( &::SireMaths::RanGenerator::seedGlobal );
            
            RanGenerator_exposer.def( 
                "seedGlobal"
                , seedGlobal_function_value
                , "Seed the global random number generator" );
        
        }
        { //::SireMaths::RanGenerator::seedGlobal
        
            typedef void ( *seedGlobal_function_type )( ::quint32 );
            seedGlobal_function_type seedGlobal_function_value( &::SireMaths::RanGenerator::seedGlobal );
            
            RanGenerator_exposer.def( 
                "seedGlobal"
                , seedGlobal_function_value
                , ( bp::arg("seed") )
                , "Seed the global random number generator" );
        
        }
        { //::SireMaths::RanGenerator::seedGlobal
        
            typedef void ( *seedGlobal_function_type )( ::QVector< unsigned int > const & );
            seedGlobal_function_type seedGlobal_function_value( &::SireMaths::RanGenerator::seedGlobal );
            
            RanGenerator_exposer.def( 
                "seedGlobal"
                , seedGlobal_function_value
                , ( bp::arg("seed") )
                , "Seed the global random number generator" );
        
        }
        { //::SireMaths::RanGenerator::seedGlobal
        
            typedef void ( *seedGlobal_function_type )( ::SireMaths::RanGenerator const & );
            seedGlobal_function_type seedGlobal_function_value( &::SireMaths::RanGenerator::seedGlobal );
            
            RanGenerator_exposer.def( 
                "seedGlobal"
                , seedGlobal_function_value
                , ( bp::arg("other") )
                , "Seed the global random number generator" );
        
        }
        { //::SireMaths::RanGenerator::setState
        
            typedef void ( ::SireMaths::RanGenerator::*setState_function_type)( ::QVector< unsigned int > const & ) ;
            setState_function_type setState_function_value( &::SireMaths::RanGenerator::setState );
            
            RanGenerator_exposer.def( 
                "setState"
                , setState_function_value
                , ( bp::arg("state") )
                , "Load the state into this generator - the state must have\nbeen produced by the getState() function above.\nThis will detach this copy from shared storage.\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMaths::RanGenerator::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::RanGenerator::typeName );
            
            RanGenerator_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMaths::RanGenerator::unlock
        
            typedef void ( ::SireMaths::RanGenerator::*unlock_function_type)(  ) const;
            unlock_function_type unlock_function_value( &::SireMaths::RanGenerator::unlock );
            
            RanGenerator_exposer.def( 
                "unlock"
                , unlock_function_value
                , "Release the generator lock" );
        
        }
        { //::SireMaths::RanGenerator::vectorOnSphere
        
            typedef ::SireMaths::Vector ( ::SireMaths::RanGenerator::*vectorOnSphere_function_type)(  ) const;
            vectorOnSphere_function_type vectorOnSphere_function_value( &::SireMaths::RanGenerator::vectorOnSphere );
            
            RanGenerator_exposer.def( 
                "vectorOnSphere"
                , vectorOnSphere_function_value
                , "Return a random vector on the unit sphere" );
        
        }
        { //::SireMaths::RanGenerator::vectorOnSphere
        
            typedef ::SireMaths::Vector ( ::SireMaths::RanGenerator::*vectorOnSphere_function_type)( double ) const;
            vectorOnSphere_function_type vectorOnSphere_function_value( &::SireMaths::RanGenerator::vectorOnSphere );
            
            RanGenerator_exposer.def( 
                "vectorOnSphere"
                , vectorOnSphere_function_value
                , ( bp::arg("radius") )
                , "Return a random vector on the sphere with radius radius" );
        
        }
        { //::SireMaths::RanGenerator::what
        
            typedef char const * ( ::SireMaths::RanGenerator::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMaths::RanGenerator::what );
            
            RanGenerator_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        RanGenerator_exposer.staticmethod( "global" );
        RanGenerator_exposer.staticmethod( "seedGlobal" );
        RanGenerator_exposer.staticmethod( "typeName" );
        RanGenerator_exposer.def( "__copy__", &__copy__);
        RanGenerator_exposer.def( "__deepcopy__", &__copy__);
        RanGenerator_exposer.def( "clone", &__copy__);
        RanGenerator_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::RanGenerator >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RanGenerator_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::RanGenerator >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RanGenerator_exposer.def( "__str__", &pvt_get_name);
        RanGenerator_exposer.def( "__repr__", &pvt_get_name);
    }

}
