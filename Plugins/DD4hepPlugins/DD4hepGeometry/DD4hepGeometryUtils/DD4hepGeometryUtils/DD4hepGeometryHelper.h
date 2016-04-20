///////////////////////////////////////////////////////////////////
// DD4hepGeometryHelper.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPGEOMETRYHELPER_H
#define ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPGEOMETRYHELPER_H

// Core module
#include "Algebra/AlgebraDefinitions.h"
// Geometry Module
#include "Volumes/VolumeBounds.h"
// DD4hep
#include "DD4hep/Detector.h"
// Gaudi
#include "GaudiKernel/GaudiException.h"

class MsgStream;

namespace Acts {
    
    /** @ class DD4hepGeometryHelper
     
     Provides helper function to translate the DD4hep geometry into the ACTS Tracking Geometry.
     
     @author julia.hrdinka@cern.ch
     */
    
    class DD4hepGeometryHelper {
    
    public:
        /** constructor */
        DD4hepGeometryHelper()
	    {}
        /** destructor */
        ~DD4hepGeometryHelper()
        {}
        /**helper method to extract the transformation matrix from a DD4hep DetElement*/
        static std::shared_ptr<Acts::Transform3D> extractTransform(DD4hep::Geometry::DetElement& detElement);
        /**helper method to extract the volume boundaries of a cylindrical volume*/
        static std::shared_ptr<const Acts::VolumeBounds> extractVolumeBounds(DD4hep::Geometry::DetElement& detElement);
    };
}

#endif //ACTS_DD4HEPGEOMETRYTOOLS_DD4HEPGEOMETRYHELPER_H