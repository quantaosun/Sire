// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Selector_Chain_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "chain.h"

#include "chaineditor.h"

#include "chainresid.h"

#include "evaluator.h"

#include "groupatomids.h"

#include "molecule.h"

#include "mover.hpp"

#include "mover_metaid.h"

#include "residue.h"

#include "selector.hpp"

#include "chain.h"

SireMol::Selector<SireMol::Chain> __copy__(const SireMol::Selector<SireMol::Chain> &other){ return SireMol::Selector<SireMol::Chain>(other); }

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_Selector_Chain__class(){

    { //::SireMol::Selector< SireMol::Chain >
        typedef bp::class_< SireMol::Selector< SireMol::Chain >, bp::bases< SireMol::MoleculeView, SireBase::Property > > Selector_Chain__exposer_t;
        Selector_Chain__exposer_t Selector_Chain__exposer = Selector_Chain__exposer_t( "Selector_Chain_", bp::init< >() );
        bp::scope Selector_Chain__scope( Selector_Chain__exposer );
        Selector_Chain__exposer.def( bp::init< SireMol::Chain const & >(( bp::arg("view") )) );
        Selector_Chain__exposer.def( bp::init< SireMol::MoleculeData const & >(( bp::arg("moldata") )) );
        Selector_Chain__exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::AtomSelection const & >(( bp::arg("moldata"), bp::arg("selected_atoms") )) );
        Selector_Chain__exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::Chain::ID const & >(( bp::arg("moldata"), bp::arg("viewid") )) );
        Selector_Chain__exposer.def( bp::init< SireMol::MoleculeData const &, QList< SireMol::ChainIdx > const & >(( bp::arg("moldata"), bp::arg("idxs") )) );
        Selector_Chain__exposer.def( bp::init< SireMol::Selector< SireMol::Chain > const & >(( bp::arg("other") )) );
        { //::SireMol::Selector< SireMol::Chain >::add
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*add_function_type)( ::SireMol::Selector< SireMol::Chain > const & ) const;
            add_function_type add_function_value( &::SireMol::Selector< SireMol::Chain >::add );
            
            Selector_Chain__exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::add
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*add_function_type)( ::SireMol::Chain const & ) const;
            add_function_type add_function_value( &::SireMol::Selector< SireMol::Chain >::add );
            
            Selector_Chain__exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("view") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::add
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*add_function_type)( ::SireMol::Chain::ID const & ) const;
            add_function_type add_function_value( &::SireMol::Selector< SireMol::Chain >::add );
            
            Selector_Chain__exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("id") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::at
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Chain ( ::SireMol::Selector< SireMol::Chain >::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMol::Selector< SireMol::Chain >::at );
            
            Selector_Chain__exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::at
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*at_function_type)( int,int ) const;
            at_function_type at_function_value( &::SireMol::Selector< SireMol::Chain >::at );
            
            Selector_Chain__exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i"), bp::arg("j") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::contains
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Chain >::*contains_function_type)( ::SireMol::Selector< SireMol::Chain > const & ) const;
            contains_function_type contains_function_value( &::SireMol::Selector< SireMol::Chain >::contains );
            
            Selector_Chain__exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::contains
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Chain >::*contains_function_type)( ::SireMol::Chain const & ) const;
            contains_function_type contains_function_value( &::SireMol::Selector< SireMol::Chain >::contains );
            
            Selector_Chain__exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("view") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::contains
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Chain >::*contains_function_type)( ::SireMol::Chain::ID const & ) const;
            contains_function_type contains_function_value( &::SireMol::Selector< SireMol::Chain >::contains );
            
            Selector_Chain__exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("id") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::count
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef int ( ::SireMol::Selector< SireMol::Chain >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::Selector< SireMol::Chain >::count );
            
            Selector_Chain__exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::evaluate
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Evaluator ( ::SireMol::Selector< SireMol::Chain >::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Selector< SireMol::Chain >::evaluate );
            
            Selector_Chain__exposer.def( 
                "evaluate"
                , evaluate_function_value );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::evaluate
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Evaluator ( ::SireMol::Selector< SireMol::Chain >::*evaluate_function_type)( int ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Selector< SireMol::Chain >::evaluate );
            
            Selector_Chain__exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::evaluate
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Evaluator ( ::SireMol::Selector< SireMol::Chain >::*evaluate_function_type)( int,int ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Selector< SireMol::Chain >::evaluate );
            
            Selector_Chain__exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("i"), bp::arg("j") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::hasMetadata
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Chain >::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Selector< SireMol::Chain >::hasMetadata );
            
            Selector_Chain__exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::hasMetadata
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Chain >::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Selector< SireMol::Chain >::hasMetadata );
            
            Selector_Chain__exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::hasProperty
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Chain >::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::Selector< SireMol::Chain >::hasProperty );
            
            Selector_Chain__exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::index
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Chain::Index ( ::SireMol::Selector< SireMol::Chain >::*index_function_type)( int ) const;
            index_function_type index_function_value( &::SireMol::Selector< SireMol::Chain >::index );
            
            Selector_Chain__exposer.def( 
                "index"
                , index_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::intersection
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*intersection_function_type)( ::SireMol::Selector< SireMol::Chain > const & ) const;
            intersection_function_type intersection_function_value( &::SireMol::Selector< SireMol::Chain >::intersection );
            
            Selector_Chain__exposer.def( 
                "intersection"
                , intersection_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::intersection
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*intersection_function_type)( ::SireMol::Chain const & ) const;
            intersection_function_type intersection_function_value( &::SireMol::Selector< SireMol::Chain >::intersection );
            
            Selector_Chain__exposer.def( 
                "intersection"
                , intersection_function_value
                , ( bp::arg("view") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::intersection
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*intersection_function_type)( ::SireMol::Chain::ID const & ) const;
            intersection_function_type intersection_function_value( &::SireMol::Selector< SireMol::Chain >::intersection );
            
            Selector_Chain__exposer.def( 
                "intersection"
                , intersection_function_value
                , ( bp::arg("id") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::intersects
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Chain >::*intersects_function_type)( ::SireMol::Selector< SireMol::Chain > const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Selector< SireMol::Chain >::intersects );
            
            Selector_Chain__exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::intersects
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Chain >::*intersects_function_type)( ::SireMol::Chain const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Selector< SireMol::Chain >::intersects );
            
            Selector_Chain__exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("view") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::intersects
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Chain >::*intersects_function_type)( ::SireMol::Chain::ID const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Selector< SireMol::Chain >::intersects );
            
            Selector_Chain__exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("id") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::invert
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*invert_function_type)(  ) const;
            invert_function_type invert_function_value( &::SireMol::Selector< SireMol::Chain >::invert );
            
            Selector_Chain__exposer.def( 
                "invert"
                , invert_function_value );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::isEmpty
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Chain >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::Selector< SireMol::Chain >::isEmpty );
            
            Selector_Chain__exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::metadataKeys
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::QStringList ( ::SireMol::Selector< SireMol::Chain >::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Selector< SireMol::Chain >::metadataKeys );
            
            Selector_Chain__exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::metadataKeys
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::QStringList ( ::SireMol::Selector< SireMol::Chain >::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Selector< SireMol::Chain >::metadataKeys );
            
            Selector_Chain__exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::move
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Chain > > ( ::SireMol::Selector< SireMol::Chain >::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMol::Selector< SireMol::Chain >::move );
            
            Selector_Chain__exposer.def( 
                "move"
                , move_function_value );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::move
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Chain > > ( ::SireMol::Selector< SireMol::Chain >::*move_function_type)( int ) const;
            move_function_type move_function_value( &::SireMol::Selector< SireMol::Chain >::move );
            
            Selector_Chain__exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::move
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Chain > > ( ::SireMol::Selector< SireMol::Chain >::*move_function_type)( int,int ) const;
            move_function_type move_function_value( &::SireMol::Selector< SireMol::Chain >::move );
            
            Selector_Chain__exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("i"), bp::arg("j") ) );
        
        }
        Selector_Chain__exposer.def( bp::self != bp::self );
        { //::SireMol::Selector< SireMol::Chain >::operator()
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Chain ( ::SireMol::Selector< SireMol::Chain >::*__call___function_type)( int ) const;
            __call___function_type __call___function_value( &::SireMol::Selector< SireMol::Chain >::operator() );
            
            Selector_Chain__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::operator()
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*__call___function_type)( int,int ) const;
            __call___function_type __call___function_value( &::SireMol::Selector< SireMol::Chain >::operator() );
            
            Selector_Chain__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") ) );
        
        }
        Selector_Chain__exposer.def( bp::self + bp::self );
        Selector_Chain__exposer.def( bp::self + bp::other< SireMol::ChainID >() );
        Selector_Chain__exposer.def( bp::self + bp::other< SireMol::Chain >() );
        Selector_Chain__exposer.def( bp::self - bp::self );
        Selector_Chain__exposer.def( bp::self - bp::other< SireMol::ChainID >() );
        Selector_Chain__exposer.def( bp::self - bp::other< SireMol::Chain >() );
        { //::SireMol::Selector< SireMol::Chain >::operator=
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > & ( ::SireMol::Selector< SireMol::Chain >::*assign_function_type)( ::SireMol::Selector< SireMol::Chain > const & ) ;
            assign_function_type assign_function_value( &::SireMol::Selector< SireMol::Chain >::operator= );
            
            Selector_Chain__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::operator=
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > & ( ::SireMol::Selector< SireMol::Chain >::*assign_function_type)( ::SireMol::Chain const & ) ;
            assign_function_type assign_function_value( &::SireMol::Selector< SireMol::Chain >::operator= );
            
            Selector_Chain__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("view") )
                , bp::return_self< >() );
        
        }
        Selector_Chain__exposer.def( bp::self == bp::self );
        { //::SireMol::Selector< SireMol::Chain >::operator[]
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Chain ( ::SireMol::Selector< SireMol::Chain >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Selector< SireMol::Chain >::operator[] );
            
            Selector_Chain__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::propertyKeys
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::QStringList ( ::SireMol::Selector< SireMol::Chain >::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::Selector< SireMol::Chain >::propertyKeys );
            
            Selector_Chain__exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::selectedAll
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef bool ( ::SireMol::Selector< SireMol::Chain >::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::Selector< SireMol::Chain >::selectedAll );
            
            Selector_Chain__exposer.def( 
                "selectedAll"
                , selectedAll_function_value );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::selection
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::AtomSelection ( ::SireMol::Selector< SireMol::Chain >::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMol::Selector< SireMol::Chain >::selection );
            
            Selector_Chain__exposer.def( 
                "selection"
                , selection_function_value );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::selection
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::AtomSelection ( ::SireMol::Selector< SireMol::Chain >::*selection_function_type)( int ) const;
            selection_function_type selection_function_value( &::SireMol::Selector< SireMol::Chain >::selection );
            
            Selector_Chain__exposer.def( 
                "selection"
                , selection_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::selection
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::AtomSelection ( ::SireMol::Selector< SireMol::Chain >::*selection_function_type)( int,int ) const;
            selection_function_type selection_function_value( &::SireMol::Selector< SireMol::Chain >::selection );
            
            Selector_Chain__exposer.def( 
                "selection"
                , selection_function_value
                , ( bp::arg("i"), bp::arg("j") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::selector
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*selector_function_type)(  ) const;
            selector_function_type selector_function_value( &::SireMol::Selector< SireMol::Chain >::selector );
            
            Selector_Chain__exposer.def( 
                "selector"
                , selector_function_value );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::selector
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*selector_function_type)( int ) const;
            selector_function_type selector_function_value( &::SireMol::Selector< SireMol::Chain >::selector );
            
            Selector_Chain__exposer.def( 
                "selector"
                , selector_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::selector
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*selector_function_type)( int,int ) const;
            selector_function_type selector_function_value( &::SireMol::Selector< SireMol::Chain >::selector );
            
            Selector_Chain__exposer.def( 
                "selector"
                , selector_function_value
                , ( bp::arg("i"), bp::arg("j") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::subtract
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*subtract_function_type)( ::SireMol::Selector< SireMol::Chain > const & ) const;
            subtract_function_type subtract_function_value( &::SireMol::Selector< SireMol::Chain >::subtract );
            
            Selector_Chain__exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::subtract
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*subtract_function_type)( ::SireMol::Chain const & ) const;
            subtract_function_type subtract_function_value( &::SireMol::Selector< SireMol::Chain >::subtract );
            
            Selector_Chain__exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("view") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::subtract
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Chain > ( ::SireMol::Selector< SireMol::Chain >::*subtract_function_type)( ::SireMol::Chain::ID const & ) const;
            subtract_function_type subtract_function_value( &::SireMol::Selector< SireMol::Chain >::subtract );
            
            Selector_Chain__exposer.def( 
                "subtract"
                , subtract_function_value
                , ( bp::arg("id") ) );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::toString
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef ::QString ( ::SireMol::Selector< SireMol::Chain >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Selector< SireMol::Chain >::toString );
            
            Selector_Chain__exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::Selector< SireMol::Chain >::typeName
        
            typedef SireMol::Selector< SireMol::Chain > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Selector< SireMol::Chain >::typeName );
            
            Selector_Chain__exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        Selector_Chain__exposer.staticmethod( "typeName" );
        Selector_Chain__exposer.def( "__copy__", &__copy__);
        Selector_Chain__exposer.def( "__deepcopy__", &__copy__);
        Selector_Chain__exposer.def( "clone", &__copy__);
        Selector_Chain__exposer.def( "__str__", &__str__< ::SireMol::Selector<SireMol::Chain> > );
        Selector_Chain__exposer.def( "__repr__", &__str__< ::SireMol::Selector<SireMol::Chain> > );
        Selector_Chain__exposer.def( "__len__", &__len_count< ::SireMol::Selector<SireMol::Chain> > );
    }

}
