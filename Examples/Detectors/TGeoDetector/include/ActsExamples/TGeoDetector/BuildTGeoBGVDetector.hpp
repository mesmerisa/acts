// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/SurfaceBinningMatcher.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsExamples/TGeoDetector/BuildTGeoBGVDetector.hpp"
#include "ActsExamples/TGeoDetector/TGeoDetectorOptions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/CylinderLayer.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
//#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/GeometryStatics.hpp"
//#include "Acts/Geometry/GlueVolumesDescriptor.hpp"
#include "Acts/Geometry/ILayerArrayCreator.hpp"
#include "Acts/Geometry/ITrackingVolumeArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
//#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Geometry/ConeLayer.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Geometry/ConeLayer.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"


#include "TGeoBBox.h"
#include "TGeoCone.h"
#include "TGeoTube.h"

#include "TGeoNode.h"
#include "TGeoVolume.h"

#include <list>
#include <vector>
#include <cmath>

#include "TGeoManager.h"

using MutableTrackingVolumePtr = std::shared_ptr<Acts::TrackingVolume>;
using MutableTrackingVolumeVector = std::vector<MutableTrackingVolumePtr>;


using LayerPtr = std::shared_ptr<const Acts::Layer>;
using LayerArray = Acts::BinnedArray<LayerPtr>;
using LayerVector = std::vector<LayerPtr>;

namespace ActsExamples {
namespace TGeo {

/// @brief global method to build the generic tracking geometry
// from a TGeo object.
///
/// It does *currently* not translate the material, this has
/// to be done with a material mapping stage
///
/// @tparam variable_map_t is the variable map
///
/// @param vm is the variable map from the options
template <typename variable_maps_t>
std::shared_ptr<const Acts::TrackingGeometry> buildTGeoBGVDetector(
    variable_maps_t& vm, const Acts::GeometryContext& context,
    std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>&
        detElementStore,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) {
  Acts::Logging::Level surfaceLogLevel =
      Acts::Logging::Level(vm["geo-surface-loglevel"].template as<size_t>());
  Acts::Logging::Level layerLogLevel =
      Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
  Acts::Logging::Level volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());
      

  // configure surface array creator
  Acts::SurfaceArrayCreator::Config sacConfig;
  auto surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(
      sacConfig,
      Acts::getDefaultLogger("SurfaceArrayCreator", surfaceLogLevel));
  // configure the proto layer helper
  Acts::ProtoLayerHelper::Config plhConfig;
  auto protoLayerHelper = std::make_shared<const Acts::ProtoLayerHelper>(
      plhConfig, Acts::getDefaultLogger("ProtoLayerHelper", layerLogLevel));
  // configure the layer creator that uses the surface array creator
  Acts::LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<const Acts::LayerCreator>(
      lcConfig, Acts::getDefaultLogger("LayerCreator", layerLogLevel));
  // configure the layer array creator
  Acts::LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
      lacConfig, Acts::getDefaultLogger("LayerArrayCreator", layerLogLevel));
  // tracking volume array creator
  Acts::TrackingVolumeArrayCreator::Config tvacConfig;
  auto tVolumeArrayCreator =
      std::make_shared<const Acts::TrackingVolumeArrayCreator>(
          tvacConfig,
          Acts::getDefaultLogger("TrackingVolumeArrayCreator", volumeLogLevel));
  // configure the cylinder volume helper
  Acts::CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper =
      std::make_shared<const Acts::CylinderVolumeHelper>(
          cvhConfig,
          Acts::getDefaultLogger("CylinderVolumeHelper", volumeLogLevel));

  //-------------------------------------------------------------------------------------
  // list the volume builders
  std::list<std::shared_ptr<const Acts::ITrackingVolumeBuilder>> volumeBuilders;


  // logger for later 
  //std::unique_ptr<const Acts::Logger> myLogger
  //    = Acts::getDefaultLogger("MyLogger", Acts::Logging::INFO);


  std::string rootFileName = vm["geo-tgeo-filename"].template as<std::string>();
  
  // import the file from
  TGeoManager::Import(rootFileName.c_str());

   // Bail out if you have no gGeoManager
  if (gGeoManager == nullptr) {
    //ACTS_WARNING("No gGeoManager found - bailing out.");
    std::cout << "No gGeoManager found - bailing out." << std::endl;
    return NULL;
  }

  // get the world
  auto BGVworld = gGeoManager->FindVolumeFast("wl");
  // get the worlds shape and get it's bounding box
  TGeoShape* WShape = BGVworld->GetShape();
  auto WorldShape = dynamic_cast<TGeoBBox*>(WShape);
  
  std::cout << "World dims: " << WorldShape->GetDX() << " " <<  WorldShape->GetDY() << " " << WorldShape->GetDZ() << " " << std::endl;
  
  // set the rMax parameter
  double rMax = WorldShape->GetDY()*10.0*0.25; // conversion to mm
  // set the min z of the interaction region.
  double zMinInteraction = -100.0; //-WorldShape->GetDZ()*10.0;
  
  // get the node and shape of the exit window (transition region)
  TGeoNode* EWNode = BGVworld->GetNode("exitWindConpv"); 
  TGeoShape* EWshape = EWNode->GetVolume()->GetShape();
  // get the exit windows bounding box 
  auto EWShape = dynamic_cast<TGeoBBox*>(EWshape);
  double EWHalfLength = EWShape->GetDZ()*10;
  // get the translation of the exit window
  const Double_t* EWtranslation = EWNode->GetMatrix()->GetTranslation();   
  std::cout << "Bound box shape of EW: "<< EWShape->GetDX()*10 << " " <<  EWShape->GetDY()*10 << " " << EWHalfLength << " " << std::endl;
  std::cout << "Trans of EW: "<< EWtranslation[0]*10 << " " <<  EWtranslation[1]*10 << " " << EWtranslation[2]*10 << " " << std::endl;
    
  double zMaxInteraction = 10.0*(EWtranslation[2] - EWShape->GetDZ() ); // maximum of the interaction volume is the begin of the exit window
  double zMaxTransition = zMaxInteraction + 10*2.0*EWShape->GetDZ(); // end of interaction vol plus length of EW

  auto EWConeShape = dynamic_cast<TGeoCone*>(EWShape);
  std::cout << "Cone dims: "<< EWConeShape->GetDz()*10  << " "  <<  EWConeShape->GetRmax1()*10 << " " << EWConeShape->GetRmin1()*10 << " " <<  EWConeShape->GetRmax2()*10 << " " << EWConeShape->GetRmin2()*10 << std::endl;  
  //::cout << EWtranslation[0] << " " << EWtranslation[1] << " " << EWtranslation[2] << " " <<std::endl;

  // get node and shape beam pipe 
  TGeoNode* LRBPNode = BGVworld->GetNode("lowRBPp"); 
  TGeoShape* LRBPshape = LRBPNode->GetVolume()->GetShape();
  // get the beam pipes bounding box 
  auto LRBPShape = dynamic_cast<TGeoBBox*>(LRBPshape);
  const Double_t* LRBPtranslation = LRBPNode->GetMatrix()->GetTranslation(); 
  
  std::cout << "Bound box shape of LRBP: "<< LRBPShape->GetDX()*10 << " " <<  LRBPShape->GetDY()*10 << " " << LRBPShape->GetDZ()*10 << " " << std::endl;
  
  double zMaxBeamPipe = 10.0*(LRBPtranslation[2] + LRBPShape->GetDZ()); // max z of the beampipe is it's z position plus it's length
  double rBeamPipe = LRBPShape->GetDX()*10.0; // radius of the beam pipe 
  
  Acts::GeometryContext geoCtx; // XXX needed?
  
  // Create the ConeBounds first
  double upper = EWConeShape->GetRmax1()*10 - rBeamPipe;
  //std::cout << "upper: " << upper << std::endl;
  double lower = 2.0*EWHalfLength; // length of cone
  //double helper_ang = atan(lower/upper);
  //double ConeOpeningAng = 2*(M_PI*0.5 - helper_ang);//2*(M_PI*0.5 - atan(lower/upper)); //2.0 * atan(upper/lower);//
  double ConeOpeningAng = atan(upper/lower);
  //std::cout << "lower " << lower << std::endl;
  
  //std::cout << "Cone opening angle " << ConeOpeningAng/M_PI*180 << std::endl; 
  
  //auto coneBounds = std::make_shared<Acts::ConeBounds>(ConeOpeningAng, zMaxInteraction, zMaxTransition, M_PI, 0.0);
  double LenCutOffTip = rBeamPipe*cos(ConeOpeningAng)/(sin(ConeOpeningAng));
  auto coneBounds = std::make_shared<Acts::ConeBounds>(ConeOpeningAng, LenCutOffTip, LenCutOffTip + 2.0*EWHalfLength);
  std::cout << "rMax  "<< rMax << " zMinInteraction " <<  zMinInteraction << " zMaxInteraction " << zMaxInteraction  << " zMaxTransition " << zMaxTransition << " zMaxBeampipe " << zMaxBeamPipe << " rBeamPipe " << rBeamPipe << std::endl;
  //Rmin2 is the inner radius of the larger side

  // create the rotation and translation for the cone, so it has the right orientation and place:
  Acts::Vector3D xAxis(-1, 0, 0);
  Acts::Vector3D yAxis(0, 1, 0);
  Acts::Vector3D zAxis(0, 0, -1);
  Acts::RotationMatrix3D mat; 
  mat.col(0) = xAxis;
  mat.col(1) = yAxis;
  mat.col(2) = zAxis;         
  Acts::Translation3D trans(EWtranslation[0]*10,EWtranslation[1]*10,EWtranslation[2]*10 + LenCutOffTip + EWHalfLength);   
  std::cout << "z-len of cut off cone part " << LenCutOffTip  << std::endl;
  std::cout << "translation from tgeo " << 10*EWtranslation[0] << " " << 10*EWtranslation[1] << " " << 10*EWtranslation[2] << " " <<std::endl;   
  Acts::Transform3D trafo(trans * mat);

  // create the cone layer 
  auto coneLayer = Acts::ConeLayer::create(trafo, coneBounds, nullptr);
  auto exitWindowTransform = Acts::Transform3D( Acts::Translation3D(0.,0., 0.5 * (zMaxInteraction + zMaxTransition)) * Acts::Transform3D::Identity());

  //auto exitWindowBounds = std::make_shared<const Acts::CylinderVolumeBounds>(0.,rMax, 0.5 * (zMaxInteraction + zMaxTransition));
  auto exitWindowBounds = std::make_shared<const Acts::CylinderVolumeBounds>(0.,rMax, EWHalfLength);

  //auto coneLayerArray = std::make_shared<Acts::BinnedArray<std::shared_ptr<Acts::Layer>>>(coneLayer); 
  //auto coneLayerArray = std::unique_ptr<const Acts::BinnedArray<std::shared_ptr<const Acts::Layer>>>(coneLayer);

  std::unique_ptr<const LayerArray> coneLayerArray = std::make_unique<const Acts::BinnedArrayXD<LayerPtr>>(coneLayer);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // Create tracking volumes
  
  // Create the volume which hosts the exit window
  auto exitWindowVolume = Acts::TrackingVolume::create(exitWindowTransform, exitWindowBounds, nullptr,
                            std::move(coneLayerArray), nullptr, {},
                            "ExitWindowVolume");                   

  MutableTrackingVolumeVector mtvVectorgap = {};

  // we should now have the parameters needed for the 4 Volumes
  // interaction volume
  auto interactionVolume = 
    cylinderVolumeHelper->createGapTrackingVolume(geoCtx, mtvVectorgap, nullptr, 0., rMax, zMinInteraction, zMaxInteraction, 0, true, "InteractionVolume");       

  // beampipe volume
  auto bpVolume = 
    cylinderVolumeHelper->createGapTrackingVolume(geoCtx, mtvVectorgap, nullptr, 0., rBeamPipe, zMaxTransition, zMaxBeamPipe, 0, true, "BeamPipe");
   
    /*Acts::CylinderVolumeHelper::createTrackingVolume(
    const GeometryContext& gctx, const LayerVector& layers,
    MutableTrackingVolumeVector mtvVector,
    std::shared_ptr<const IVolumeMaterial> volumeMaterial, double rMin,
    double rMax, double zMin, double zMax, const std::string& volumeName,
    BinningType bType) const {*/
    
    
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // Creating the layers for the BGV detector
  auto layerBuilderConfigs =
    ActsExamples::Options::readTGeoLayerBuilderConfigs<variable_maps_t>(vm);

  // remember the layer builders to collect the detector elements
  std::vector<std::shared_ptr<const Acts::TGeoLayerBuilder>> tgLayerBuilders;
  //std::shared_ptr<const Acts::TGeoLayerBuilder> tgLayerBuilders;
 
  for (auto& lbc : layerBuilderConfigs) {
    std::shared_ptr<const Acts::LayerCreator> layerCreatorLB = nullptr;

    if (lbc.autoSurfaceBinning) {
      std::cout << "autoSurfaceBinning "  << std::endl;
      // Configure surface array creator (optionally) per layer builder
      // (in order to configure them to work appropriately)
      Acts::SurfaceArrayCreator::Config sacConfigLB;
      sacConfigLB.surfaceMatcher = lbc.surfaceBinMatcher;
      auto surfaceArrayCreatorLB =
          std::make_shared<const Acts::SurfaceArrayCreator>(
              sacConfigLB, Acts::getDefaultLogger(
                               lbc.configurationName + "SurfaceArrayCreator",
                               surfaceLogLevel));
      // configure the layer creator that uses the surface array creator
      Acts::LayerCreator::Config lcConfigLB;
      lcConfigLB.surfaceArrayCreator = surfaceArrayCreatorLB;
      layerCreatorLB = std::make_shared<const Acts::LayerCreator>(
          lcConfigLB,
          Acts::getDefaultLogger(lbc.configurationName + "LayerCreator",
                                 layerLogLevel));
    }
    // Configure the proto layer helper
    Acts::ProtoLayerHelper::Config plhConfigLB;
    auto protoLayerHelperLB = std::make_shared<const Acts::ProtoLayerHelper>(
        plhConfigLB,
        Acts::getDefaultLogger(lbc.configurationName + "ProtoLayerHelper",
                               layerLogLevel));

    lbc.layerCreator =
        (layerCreatorLB != nullptr) ? layerCreatorLB : layerCreator;
    lbc.protoLayerHelper =
        (protoLayerHelperLB != nullptr) ? protoLayerHelperLB : protoLayerHelper;

    auto layerBuilder = std::make_shared<const Acts::TGeoLayerBuilder>(
        lbc, Acts::getDefaultLogger(lbc.configurationName + "LayerBuilder",
                                    layerLogLevel));
    // remember the layer builder
    tgLayerBuilders.push_back(layerBuilder);
  }
  
  LayerVector BGVDetLayerVector;
  
  for (auto& lBuilder : tgLayerBuilders) {
    BGVDetLayerVector = lBuilder->positiveLayers(context);  
    //std::cout << "det elements when building TGeo Detector " << std::endl;
    auto detElements = lBuilder->detectorElements();
    detElementStore.insert(detElementStore.begin(), detElements.begin(),
                           detElements.end());
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////
    
  MutableTrackingVolumeVector mtvVector = {};
 
  // BGV volume 
  auto bgvVolume = cylinderVolumeHelper->createTrackingVolume(geoCtx, BGVDetLayerVector, mtvVector, nullptr, rBeamPipe, rMax, zMaxTransition, zMaxBeamPipe, "BGV"); // XXX changed

  // first we need to stuff bpVolume & bgvVolume into a container volume
  auto bpBgvContainer = cylinderVolumeHelper->createContainerTrackingVolume(geoCtx, { bpVolume, bgvVolume} );
  // combine volumes to world container:          
  auto worldContainer = cylinderVolumeHelper->createContainerTrackingVolume(geoCtx, { interactionVolume, exitWindowVolume , bpBgvContainer } );
  // Make a TrackingGeometry out of it:
  return std::make_shared<Acts::TrackingGeometry>(worldContainer, mdecorator.get()); 
   
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

}

}  // namespace TGeo
}  // namespace ActsExamples
