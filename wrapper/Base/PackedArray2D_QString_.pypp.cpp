// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "PackedArray2D_QString_.pypp.hpp"

namespace bp = boost::python;

#include "packedarrays.h"

#include "packedarrays.h"

SireBase::PackedArray2D<QString> __copy__(const SireBase::PackedArray2D<QString> &other){ return SireBase::PackedArray2D<QString>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_PackedArray2D_QString__class(){

    { //::SireBase::PackedArray2D< QString >
        typedef bp::class_< SireBase::PackedArray2D< QString > > PackedArray2D_QString__exposer_t;
        PackedArray2D_QString__exposer_t PackedArray2D_QString__exposer = PackedArray2D_QString__exposer_t( "PackedArray2D_QString_", bp::init< >() );
        bp::scope PackedArray2D_QString__scope( PackedArray2D_QString__exposer );
        PackedArray2D_QString__exposer.def( bp::init< SireBase::PackedArray2D< QString >::Array const & >(( bp::arg("array") )) );
        PackedArray2D_QString__exposer.def( bp::init< QVector< SireBase::detail::PackedArray2D_Array< QString > > const & >(( bp::arg("arrays") )) );
        PackedArray2D_QString__exposer.def( bp::init< QVector< QString > const & >(( bp::arg("values") )) );
        PackedArray2D_QString__exposer.def( bp::init< QVector< QVector< QString > > const & >(( bp::arg("values") )) );
        PackedArray2D_QString__exposer.def( bp::init< SireBase::PackedArray2D< QString > const &, SireBase::PackedArray2D< QString > const & >(( bp::arg("array0"), bp::arg("array1") )) );
        PackedArray2D_QString__exposer.def( bp::init< SireBase::PackedArray2D< QString > const & >(( bp::arg("other") )) );
        { //::SireBase::PackedArray2D< QString >::append
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef void ( ::SireBase::PackedArray2D< QString >::*append_function_type)( ::SireBase::PackedArray2D< QString >::Array const & ) ;
            append_function_type append_function_value( &::SireBase::PackedArray2D< QString >::append );
            
            PackedArray2D_QString__exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("array") ) );
        
        }
        { //::SireBase::PackedArray2D< QString >::append
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef void ( ::SireBase::PackedArray2D< QString >::*append_function_type)( ::SireBase::PackedArray2D< QString > const & ) ;
            append_function_type append_function_value( &::SireBase::PackedArray2D< QString >::append );
            
            PackedArray2D_QString__exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("arrays") ) );
        
        }
        { //::SireBase::PackedArray2D< QString >::append
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef void ( ::SireBase::PackedArray2D< QString >::*append_function_type)( ::QVector< QString > const & ) ;
            append_function_type append_function_value( &::SireBase::PackedArray2D< QString >::append );
            
            PackedArray2D_QString__exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("array") ) );
        
        }
        { //::SireBase::PackedArray2D< QString >::append
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef void ( ::SireBase::PackedArray2D< QString >::*append_function_type)( ::QVector< QVector< QString > > const & ) ;
            append_function_type append_function_value( &::SireBase::PackedArray2D< QString >::append );
            
            PackedArray2D_QString__exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("arrays") ) );
        
        }
        { //::SireBase::PackedArray2D< QString >::assertValidIndex
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef void ( ::SireBase::PackedArray2D< QString >::*assertValidIndex_function_type)( ::quint32 ) const;
            assertValidIndex_function_type assertValidIndex_function_value( &::SireBase::PackedArray2D< QString >::assertValidIndex );
            
            PackedArray2D_QString__exposer.def( 
                "assertValidIndex"
                , assertValidIndex_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireBase::PackedArray2D< QString >::at
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef ::SireBase::PackedArray2D< QString >::Array const & ( ::SireBase::PackedArray2D< QString >::*at_function_type)( ::quint32 ) const;
            at_function_type at_function_value( &::SireBase::PackedArray2D< QString >::at );
            
            PackedArray2D_QString__exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireBase::PackedArray2D< QString >::at
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef ::QString const & ( ::SireBase::PackedArray2D< QString >::*at_function_type)( ::quint32,::quint32 ) const;
            at_function_type at_function_value( &::SireBase::PackedArray2D< QString >::at );
            
            PackedArray2D_QString__exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireBase::PackedArray2D< QString >::count
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef int ( ::SireBase::PackedArray2D< QString >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireBase::PackedArray2D< QString >::count );
            
            PackedArray2D_QString__exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireBase::PackedArray2D< QString >::fromVariant
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef ::SireBase::PackedArray2D< QString > ( *fromVariant_function_type )( ::SireBase::PackedArray2D< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireBase::PackedArray2D< QString >::fromVariant );
            
            PackedArray2D_QString__exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("variant") ) );
        
        }
        { //::SireBase::PackedArray2D< QString >::isEmpty
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef bool ( ::SireBase::PackedArray2D< QString >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireBase::PackedArray2D< QString >::isEmpty );
            
            PackedArray2D_QString__exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireBase::PackedArray2D< QString >::nArrays
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef int ( ::SireBase::PackedArray2D< QString >::*nArrays_function_type)(  ) const;
            nArrays_function_type nArrays_function_value( &::SireBase::PackedArray2D< QString >::nArrays );
            
            PackedArray2D_QString__exposer.def( 
                "nArrays"
                , nArrays_function_value );
        
        }
        { //::SireBase::PackedArray2D< QString >::nValues
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef int ( ::SireBase::PackedArray2D< QString >::*nValues_function_type)(  ) const;
            nValues_function_type nValues_function_value( &::SireBase::PackedArray2D< QString >::nValues );
            
            PackedArray2D_QString__exposer.def( 
                "nValues"
                , nValues_function_value );
        
        }
        { //::SireBase::PackedArray2D< QString >::nValues
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef int ( ::SireBase::PackedArray2D< QString >::*nValues_function_type)( ::quint32 ) const;
            nValues_function_type nValues_function_value( &::SireBase::PackedArray2D< QString >::nValues );
            
            PackedArray2D_QString__exposer.def( 
                "nValues"
                , nValues_function_value
                , ( bp::arg("i") ) );
        
        }
        PackedArray2D_QString__exposer.def( bp::self != bp::self );
        { //::SireBase::PackedArray2D< QString >::operator()
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef ::QString const & ( ::SireBase::PackedArray2D< QString >::*__call___function_type)( ::quint32,::quint32 ) const;
            __call___function_type __call___function_value( &::SireBase::PackedArray2D< QString >::operator() );
            
            PackedArray2D_QString__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireBase::PackedArray2D< QString >::operator=
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef ::SireBase::PackedArray2D< QString > & ( ::SireBase::PackedArray2D< QString >::*assign_function_type)( ::SireBase::PackedArray2D< QString > const & ) ;
            assign_function_type assign_function_value( &::SireBase::PackedArray2D< QString >::operator= );
            
            PackedArray2D_QString__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        PackedArray2D_QString__exposer.def( bp::self == bp::self );
        { //::SireBase::PackedArray2D< QString >::operator[]
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef ::SireBase::PackedArray2D< QString >::Array const & ( ::SireBase::PackedArray2D< QString >::*__getitem___function_type)( ::quint32 ) const;
            __getitem___function_type __getitem___function_value( &::SireBase::PackedArray2D< QString >::operator[] );
            
            PackedArray2D_QString__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireBase::PackedArray2D< QString >::remove
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef void ( ::SireBase::PackedArray2D< QString >::*remove_function_type)( ::quint32 ) ;
            remove_function_type remove_function_value( &::SireBase::PackedArray2D< QString >::remove );
            
            PackedArray2D_QString__exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireBase::PackedArray2D< QString >::removeAll
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef void ( ::SireBase::PackedArray2D< QString >::*removeAll_function_type)( ::QVarLengthArray< int, 256 > const & ) ;
            removeAll_function_type removeAll_function_value( &::SireBase::PackedArray2D< QString >::removeAll );
            
            PackedArray2D_QString__exposer.def( 
                "removeAll"
                , removeAll_function_value
                , ( bp::arg("idxs") ) );
        
        }
        { //::SireBase::PackedArray2D< QString >::size
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef int ( ::SireBase::PackedArray2D< QString >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireBase::PackedArray2D< QString >::size );
            
            PackedArray2D_QString__exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireBase::PackedArray2D< QString >::toQVector
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef ::QVector< QString > ( ::SireBase::PackedArray2D< QString >::*toQVector_function_type)(  ) const;
            toQVector_function_type toQVector_function_value( &::SireBase::PackedArray2D< QString >::toQVector );
            
            PackedArray2D_QString__exposer.def( 
                "toQVector"
                , toQVector_function_value );
        
        }
        { //::SireBase::PackedArray2D< QString >::toQVectorVector
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef ::QVector< QVector< QString > > ( ::SireBase::PackedArray2D< QString >::*toQVectorVector_function_type)(  ) const;
            toQVectorVector_function_type toQVectorVector_function_value( &::SireBase::PackedArray2D< QString >::toQVectorVector );
            
            PackedArray2D_QString__exposer.def( 
                "toQVectorVector"
                , toQVectorVector_function_value );
        
        }
        { //::SireBase::PackedArray2D< QString >::toString
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef ::QString ( ::SireBase::PackedArray2D< QString >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireBase::PackedArray2D< QString >::toString );
            
            PackedArray2D_QString__exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireBase::PackedArray2D< QString >::toVariant
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef ::SireBase::PackedArray2D< QVariant > ( ::SireBase::PackedArray2D< QString >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireBase::PackedArray2D< QString >::toVariant );
            
            PackedArray2D_QString__exposer.def( 
                "toVariant"
                , toVariant_function_value );
        
        }
        { //::SireBase::PackedArray2D< QString >::update
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef void ( ::SireBase::PackedArray2D< QString >::*update_function_type)( ::quint32,::SireBase::PackedArray2D< QString >::Array const & ) ;
            update_function_type update_function_value( &::SireBase::PackedArray2D< QString >::update );
            
            PackedArray2D_QString__exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("i"), bp::arg("array") ) );
        
        }
        { //::SireBase::PackedArray2D< QString >::update
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef void ( ::SireBase::PackedArray2D< QString >::*update_function_type)( ::quint32,::QVector< QString > const & ) ;
            update_function_type update_function_value( &::SireBase::PackedArray2D< QString >::update );
            
            PackedArray2D_QString__exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("i"), bp::arg("array") ) );
        
        }
        { //::SireBase::PackedArray2D< QString >::updateAll
        
            typedef SireBase::PackedArray2D< QString > exported_class_t;
            typedef void ( ::SireBase::PackedArray2D< QString >::*updateAll_function_type)( ::QVarLengthArray< int, 256 > const &,::SireBase::PackedArray2D< QString > const & ) ;
            updateAll_function_type updateAll_function_value( &::SireBase::PackedArray2D< QString >::updateAll );
            
            PackedArray2D_QString__exposer.def( 
                "updateAll"
                , updateAll_function_value
                , ( bp::arg("idxs"), bp::arg("arrays") ) );
        
        }
        PackedArray2D_QString__exposer.staticmethod( "fromVariant" );
        PackedArray2D_QString__exposer.def( "__copy__", &__copy__);
        PackedArray2D_QString__exposer.def( "__deepcopy__", &__copy__);
        PackedArray2D_QString__exposer.def( "clone", &__copy__);
        PackedArray2D_QString__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::PackedArray2D<QString> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PackedArray2D_QString__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::PackedArray2D<QString> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PackedArray2D_QString__exposer.def( "__str__", &__str__< ::SireBase::PackedArray2D<QString> > );
        PackedArray2D_QString__exposer.def( "__repr__", &__str__< ::SireBase::PackedArray2D<QString> > );
        PackedArray2D_QString__exposer.def( "__len__", &__len_size< ::SireBase::PackedArray2D<QString> > );
    }

}
