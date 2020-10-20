// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/AdaptiveMultiVertexFinderAlgorithmFromTrajBGV.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include "VertexingHelpers.hpp"

ActsExamples::AdaptiveMultiVertexFinderAlgorithmFromTrajBGV::
    AdaptiveMultiVertexFinderAlgorithmFromTrajBGV(const Config& cfg,
                                       Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("AdaptiveMultiVertexFinder", lvl),
      m_cfg(cfg) {
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input track trajectories collection");
  }
  if (m_cfg.outputProtoVertices.empty()) {
    throw std::invalid_argument("Missing output proto vertices collection");
  }
}

ActsExamples::ProcessCode
ActsExamples::AdaptiveMultiVertexFinderAlgorithmFromTrajBGV::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // retrieve input tracks and convert into the expected format
  // XXX
  //const auto& inputTrackParameters =
  //    ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackParameters);
  
  //////////////////////////////////////////////////////////////////////////////////////
  // retrieve input trajectories and convert into the expected format    
  const auto& trajectoriesContainer =
      ctx.eventStore.get<TrajectoryContainer>(m_cfg.inputTrajectories);
    
  std::vector<Acts::BoundTrackParameters> trackParameters;
      
  int traj_counter = -1;
  
  std::vector<int> empty_traj;
  empty_traj.clear();
  
  for (const auto& traj : trajectoriesContainer ) {
     traj_counter += 1;
     std::cout << "traj " << std::endl;
     //std::cout << traj.hasTrackParameters() << std::endl;
     const auto& [trackTips, mj] = traj.trajectory();

     if (trackTips.empty()) {
       ACTS_WARNING("Empty multiTrajectory.");
       empty_traj.push_back(traj_counter);
       continue;
     }
     // Check the size of the trajectory entry indices. For track fitting, there
     // should be at most one trajectory
     if (trackTips.size() > 1) {
       ACTS_ERROR("Track fitting should not result in multiple trajectories.");
       return ProcessCode::ABORT;
     }
     auto& trackTip = trackTips.front();
     // Select reco track with fitted parameters
     if (not traj.hasTrackParameters(trackTip)) {
       ACTS_WARNING("No fitted track parameters.");
       empty_traj.push_back(traj_counter);
       continue;
     }
     if (not traj.hasTrajectory(trackTip)) {
       ACTS_WARNING("No trajectory.");
       empty_traj.push_back(traj_counter);
       continue;
     }
     if (not traj.trackParameters(trackTip).covariance()) {
       ACTS_WARNING("No covariance matrix.");
       empty_traj.push_back(traj_counter);
       continue;     
     }
     
     //std::cout <<  "cov    " << traj.trackParameters(trackTip).covariance().empty() << std::endl;
     //const auto& fittedParameters = traj.trackParameters(trackTip);
     trackParameters.push_back(traj.trackParameters(trackTip));   
  }
  
  for(auto btr_ : empty_traj) {
     std::cout << "empty traj    " << btr_ << std::endl;
  }
   ///////////////////////////////////////////////////////////////////////////////////////////////////    
      
      
  const auto& inputTrackPointers =
      makeTrackParametersPointerContainer(trackParameters);      
   
  // XXX 
  //const auto& inputTrackPointers =
  //    makeTrackParametersPointerContainer(inputTrackParameters);

  //////////////////////////////////////////////
  /* Full tutorial example code for reference */
  //////////////////////////////////////////////

  // Set up EigenStepper
  Acts::ConstantBField bField(m_cfg.bField);
  Acts::EigenStepper<Acts::ConstantBField> stepper(bField);

  // Set up the propagator
  using Propagator = Acts::Propagator<Acts::EigenStepper<Acts::ConstantBField>>;
  auto propagator = std::make_shared<Propagator>(stepper);

  // Set up ImpactPointEstimator
  using IPEstimator =
      Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;
  IPEstimator::Config ipEstimatorCfg(bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  // Set up the helical track linearizer
  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Set up deterministic annealing with user-defined temperatures
  std::vector<double> temperatures{16.0, 8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  //std::vector<double> temperatures{64., 16., 4., 2., 1.5, 1.};
  Acts::AnnealingUtility::Config annealingConfig(temperatures);
  Acts::AnnealingUtility annealingUtility(annealingConfig);

  // Set up the vertex fitter with user-defined annealing
  using Fitter =
      Acts::AdaptiveMultiVertexFitter<Acts::BoundTrackParameters, Linearizer>;
  Fitter::Config fitterCfg(ipEstimator);
  fitterCfg.annealingTool = annealingUtility;
  Fitter fitter(fitterCfg);

  // Set up the vertex seed finder
  using SeedFinder = Acts::TrackDensityVertexFinder<
      Fitter, Acts::GaussianTrackDensity<Acts::BoundTrackParameters>>;
  SeedFinder seedFinder;

  // The vertex finder type
  using Finder = Acts::AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              linearizer);
                              
  // We do not want to use a beamspot constraint here    // XXX
  finderConfig.useBeamSpotConstraint = false;
  //finderConfig.useBeamSpotConstraint = true;  
  finderConfig.maxIterations = 1000;         
  finderConfig.tracksMaxZinterval = 100.0 * Acts::UnitConstants::mm;

  // Instantiate the finder
  Finder finder(finderConfig);
  // The vertex finder state
  Finder::State state;

  // Default vertexing options, this is where e.g. a constraint could be set
  using VertexingOptions = Acts::VertexingOptions<Acts::BoundTrackParameters>;
  VertexingOptions finderOpts(ctx.geoContext, ctx.magFieldContext);
  
  /*Acts::Vector3D constraintPos = Acts::Vector3D(0, 0, 0);
  Acts::Vertex<Acts::BoundTrackParameters> vertexConstraint(constraintPos);
  Acts::SymMatrix3D constraintCov =
        Acts::Vector3D(100 * Acts::UnitConstants::mm, 100 * Acts::UnitConstants::mm,
                       100 * Acts::UnitConstants::mm)
            .asDiagonal();
  vertexConstraint.setCovariance(constraintCov);
  
  finderOpts.vertexConstraint = vertexConstraint;*/
  
  /*/// Vertex constraint position
  Acts::Vector3D constraintPos = Acts::Vector3D(0, 0, 0);
  /// Vertex constraint covariance matrix
  Acts::SymMatrix3D constraintCov =
        Acts::Vector3D(100 * Acts::UnitConstants::mm, 100 * Acts::UnitConstants::mm,
                       100 * Acts::UnitConstants::mm)
            .asDiagonal();*/

  // find vertices
  auto result = finder.find(inputTrackPointers, finderOpts, state);
  if (not result.ok()) {
    ACTS_ERROR("Error in vertex finder: " << result.error().message());
    return ProcessCode::ABORT;
  }
  auto vertices = *result;

  // show some debug output
  ACTS_INFO("Found " << vertices.size() << " vertices in event");
  for (const auto& vtx : vertices) {
    ACTS_INFO("Found vertex at " << vtx.fullPosition().transpose() << " with "
                                 << vtx.tracks().size() << " tracks.");
  }

  // store proto vertices extracted from the found vertices
  ctx.eventStore.add(m_cfg.outputProtoVertices,
                     makeProtoVertices(trackParameters, vertices));
                     
  ctx.eventStore.add(m_cfg.outputFoundVertices, std::move(vertices));                   

  return ActsExamples::ProcessCode::SUCCESS;
}
