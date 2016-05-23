//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "SireMaths_registrars.h"

#include "freeenergyaverage.h"
#include "accumulator.h"
#include "align.h"
#include "histogram.h"
#include "trigmatrix.h"
#include "nmatrix.h"
#include "plane.h"
#include "vectorproperty.h"
#include "matrix.h"
#include "distvector.h"
#include "vector.h"
#include "axisset.h"
#include "torsion.h"
#include "triangle.h"
#include "complex.h"
#include "nvector.h"
#include "quaternion.h"
#include "rangenerator.h"
#include "n4matrix.h"
#include "line.h"
#include "sphere.h"

#include "Helpers/objectregistry.hpp"

void register_SireMaths_objects()
{

    ObjectRegistry::registerConverterFor< SireMaths::FreeEnergyAverage >();
    ObjectRegistry::registerConverterFor< SireMaths::BennettsFreeEnergyAverage >();
    ObjectRegistry::registerConverterFor< SireMaths::NullAccumulator >();
    ObjectRegistry::registerConverterFor< SireMaths::Average >();
    ObjectRegistry::registerConverterFor< SireMaths::AverageAndStddev >();
    ObjectRegistry::registerConverterFor< SireMaths::ExpAverage >();
    ObjectRegistry::registerConverterFor< SireMaths::Median >();
    ObjectRegistry::registerConverterFor< SireMaths::RecordValues >();
    ObjectRegistry::registerConverterFor< SireMaths::Transform >();
    ObjectRegistry::registerConverterFor< SireMaths::Histogram >();
    ObjectRegistry::registerConverterFor< SireMaths::TrigMatrix >();
    ObjectRegistry::registerConverterFor< SireMaths::NMatrix >();
    ObjectRegistry::registerConverterFor< SireMaths::Plane >();
    ObjectRegistry::registerConverterFor< SireMaths::VectorProperty >();
    ObjectRegistry::registerConverterFor< SireMaths::VectorArrayProperty >();
    ObjectRegistry::registerConverterFor< SireMaths::Matrix >();
    ObjectRegistry::registerConverterFor< SireMaths::DistVector >();
    ObjectRegistry::registerConverterFor< SireMaths::Vector >();
    ObjectRegistry::registerConverterFor< SireMaths::AxisSet >();
    ObjectRegistry::registerConverterFor< SireMaths::Torsion >();
    ObjectRegistry::registerConverterFor< SireMaths::Triangle >();
    ObjectRegistry::registerConverterFor< SireMaths::Complex >();
    ObjectRegistry::registerConverterFor< SireMaths::NVector >();
    ObjectRegistry::registerConverterFor< SireMaths::Quaternion >();
    ObjectRegistry::registerConverterFor< SireMaths::RanGenerator >();
    ObjectRegistry::registerConverterFor< SireMaths::N4Matrix >();
    ObjectRegistry::registerConverterFor< SireMaths::Line >();
    ObjectRegistry::registerConverterFor< SireMaths::Sphere >();

}

