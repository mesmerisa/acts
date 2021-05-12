// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/VertexFitterAlgorithmFromTraj.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/SimHit.hpp"

#include <stdexcept>

ActsExamples::VertexFitterAlgorithmFromTraj::VertexFitterAlgorithmFromTraj(
    const Config& cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("VertexFit", lvl), m_cfg(cfg) {
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input trajectory collection");
  }
  if (m_cfg.inputProtoVertices.empty()) {
    throw std::invalid_argument("Missing input proto vertices collection");
  }
}

ActsExamples::ProcessCode ActsExamples::VertexFitterAlgorithmFromTraj::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
  using PropagatorOptions = Acts::PropagatorOptions<>;
  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  using VertexFitter =
      Acts::FullBilloirVertexFitter<Acts::BoundTrackParameters, Linearizer>;
  using VertexFitterOptions =
      Acts::VertexingOptions<Acts::BoundTrackParameters>;
  using HitSimHitsMap = IndexMultimap<Index>;
       
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);     
  const auto& hitSimHitsMap =
      ctx.eventStore.get<HitSimHitsMap>(m_cfg.inputMeasurementSimHitsMap);
      
  // Set up EigenStepper
  Acts::EigenStepper<> stepper(m_cfg.bField);

  // Setup the propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);
  PropagatorOptions propagatorOpts(ctx.geoContext, ctx.magFieldContext,
                                   Acts::LoggerWrapper{logger()});
  // Setup the vertex fitter
  VertexFitter::Config vertexFitterCfg;
  VertexFitter vertexFitter(vertexFitterCfg);
  VertexFitter::State state(m_cfg.bField->makeCache(ctx.magFieldContext));
  // Setup the linearizer
  Linearizer::Config ltConfig(m_cfg.bField, propagator);
  Linearizer linearizer(ltConfig);

  ////////////////////////////////////////////////////////////////////////////////
  
  //const auto& trackParameters =
  //    ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackParameters);
  
  ///////////////////////////////
  
  const auto& trajectories =
      ctx.eventStore.get<TrajectoriesContainer>(m_cfg.inputTrajectories);

  std::vector<Acts::BoundTrackParameters> inputTrackParameters;
   inputTrackParameters.clear();
  
                        
  int traj_counter = -1;
  
  std::vector<int> empty_traj;
  empty_traj.clear();
  
  //std::cout << "size of trajectories at start: " << trajectories.size() << std::endl;
  
  
  for (size_t itraj = 0; itraj < trajectories.size(); ++itraj) {
    const auto& traj = trajectories[itraj];
  
     traj_counter += 1;
     
     //std::cout << traj.hasTrackParameters() << std::endl;
     //const auto& [trackTips, mj] = traj;
     //m_trajNr = itraj;

     // The trajectory entry indices and the multiTrajectory
     const auto& mj = traj.multiTrajectory();
     const auto& trackTips = traj.tips();
     
     if (trackTips.empty()) {
       ACTS_WARNING("Empty multiTrajectory.");
       empty_traj.push_back(traj_counter);
       continue;
     }
     // Check the size of the trajectory entry indices. For track fitting, there
     // should be at most one trajectory
     //if (trackTips.size() > 1) {
     //  ACTS_ERROR("Track fitting should not result in multiple trajectories."); 
     //  return ProcessCode::ABORT;
     //}
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
       
      // std::cout << "trajectory counter " << traj_counter << std::endl; 
     }
     const auto fitted_params = traj.trackParameters(trackTip);
     inputTrackParameters.push_back(fitted_params);
     
     if(trajectories.size() >= 2) {
     
       /*auto ts = mj.getTrackState(1);
     
       const auto& surface =fitted_params.referenceSurface();

       // get the truth hits corresponding to this trackState
       // Use average truth in the case of multiple contributing sim hits
       const auto hitIdx = ts.uncalibrated().index();
       auto indices = makeRange(hitSimHitsMap.equal_range(hitIdx));
     
       auto [truthLocal, truthPos4, truthUnitDir] =
            averageSimHits(ctx.geoContext, surface, simHits, indices);
     
       const auto trackicovi = *traj.trackParameters(trackTip).covariance();
        
     
        float truthLOC0 = truthLocal[Acts::ePos0];
        float truthLOC1 = truthLocal[Acts::ePos1];
        //float truthTIME = truthPos4[Acts::eTime];
        float truthPHI = Acts::VectorHelpers::phi(truthUnitDir);
        float truthTHETA = Acts::VectorHelpers::theta(truthUnitDir);
        
    
        
       std::cout << truthLOC0 << " ";
       std::cout << truthLOC1 << " ";
       std::cout << truthPHI << " ";
       std::cout << truthTHETA << " ";
     
     
       
       std::cout << fitted_params.parameters()[Acts::eBoundLoc0] << " ";
       std::cout << fitted_params.parameters()[Acts::eBoundLoc1] << " ";
       std::cout << fitted_params.parameters()[Acts::eBoundPhi] << " ";
       std::cout << fitted_params.parameters()[Acts::eBoundTheta] << " ";
       std::cout << fitted_params.parameters()[Acts::eBoundQOverP] << " ";
       
     
       
        std::cout << sqrt(trackicovi(Acts::eBoundLoc0, Acts::eBoundLoc0)) << " ";
        std::cout << sqrt(trackicovi(Acts::eBoundLoc1, Acts::eBoundLoc1)) << " ";
        std::cout << sqrt(trackicovi(Acts::eBoundPhi, Acts::eBoundPhi)) << " ";
        std::cout << sqrt(trackicovi(Acts::eBoundTheta, Acts::eBoundTheta)) << " ";
        std::cout << sqrt(trackicovi(Acts::eBoundQOverP, Acts::eBoundQOverP)) << " ";
       
        std::cout << std::endl;
        //std::cout << " -------------- " << std::endl;*/
     }
        
        
        
     //const auto& measvec = traj.calibrated();
     //trackParameters.push_back(traj.trackParameters(trackTip));
   
  }

  for(auto btr_ : empty_traj) {
     std::cout << "empty trajectory index   " << btr_ << std::endl;
  }
   //std::cout << " ------------------------------------------- " << std::endl;
  
  
   std::vector<Acts::Vertex<Acts::BoundTrackParameters>> outputVertexCollection;
  
    outputVertexCollection.clear();
  
  
  
  //////////////////////////////////////////////////////////////////////////////////
  
  const auto& protoVertices =
      ctx.eventStore.get<ProtoVertexContainer>(m_cfg.inputProtoVertices);
  std::vector<const Acts::BoundTrackParameters*> inputTrackPtrCollection;

  for (const auto& protoVertex : protoVertices) {
    // un-constrained fit requires at least two tracks
    if ((not m_cfg.doConstrainedFit) and (protoVertex.size() < 2)) {
      ACTS_WARNING(
          "Skip un-constrained vertex fit on proto-vertex with less than two "
          "tracks");
      continue;
    }
//////////////////////////////////////////////////////////////////////////////////
    // select input tracks for the input proto vertex
    inputTrackPtrCollection.clear();
    //inputTrackPtrCollection.reserve(protoVertex.size());
      
    int track_index_not_empy = -1;
    for (const auto& trackIdx : protoVertex) {
      track_index_not_empy += 1;
      if (std::count(empty_traj.begin(), empty_traj.end(), trackIdx)) {
        track_index_not_empy -= 1;
        continue;
      }
      inputTrackPtrCollection.push_back(&inputTrackParameters[track_index_not_empy]);     
    }
//////////////////////////////////////////////////////////////////////////////////

    // select input tracks for the input proto vertex
    /*inputTrackPtrCollection.clear();
    inputTrackPtrCollection.reserve(protoVertex.size());
    for (const auto& trackIdx : protoVertex) {
      inputTrackPtrCollection.push_back(&trackParameters[trackIdx]);
    }*/
    
   

    Acts::Vertex<Acts::BoundTrackParameters> fittedVertex;
    if (!m_cfg.doConstrainedFit) {
      VertexFitterOptions vfOptions(ctx.geoContext, ctx.magFieldContext);

      auto fitRes = vertexFitter.fit(inputTrackPtrCollection, linearizer,
                                     vfOptions, state);
      if (fitRes.ok()) {
        fittedVertex = *fitRes;
        outputVertexCollection.push_back(fittedVertex);
      } else {
        ACTS_ERROR("Error in vertex fit.");
        ACTS_ERROR(fitRes.error().message());
      }
    } else {
      // Vertex constraint
      Acts::Vertex<Acts::BoundTrackParameters> theConstraint;

      theConstraint.setCovariance(m_cfg.constraintCov);
      theConstraint.setPosition(m_cfg.constraintPos);

      // Vertex fitter options
      VertexFitterOptions vfOptionsConstr(ctx.geoContext, ctx.magFieldContext,
                                          theConstraint);

      auto fitRes = vertexFitter.fit(inputTrackPtrCollection, linearizer,
                                     vfOptionsConstr, state);
      if (fitRes.ok()) {
        fittedVertex = *fitRes;
        outputVertexCollection.push_back(fittedVertex);
      } else {
        ACTS_ERROR("Error in vertex fit with constraint.");
        ACTS_ERROR(fitRes.error().message());
      }
    }

    ACTS_INFO("Fitted Vertex " << fittedVertex.fullPosition().transpose());
    ACTS_INFO("Tracks at fitted Vertex: " << fittedVertex.tracks().size());
    
    
  }
  
  ctx.eventStore.add(m_cfg.outputFittedVertices, std::move(outputVertexCollection));
  
  
  return ProcessCode::SUCCESS;
}
