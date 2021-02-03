// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/VertexFitterAlgorithmFromTrajBGV.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
//#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

ActsExamples::VertexFitterAlgorithmFromTrajBGV::VertexFitterAlgorithmFromTrajBGV(
    const Config& cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("VertexFit", lvl), m_cfg(cfg) {
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input trajectory collection");
  }
  if (m_cfg.inputProtoVertices.empty()) {
    throw std::invalid_argument("Missing input proto vertices collection");
  }
}

ActsExamples::ProcessCode ActsExamples::VertexFitterAlgorithmFromTrajBGV::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  using MagneticField = Acts::ConstantBField;
  using Stepper = Acts::EigenStepper<MagneticField>;
  using Propagator = Acts::Propagator<Stepper>;
  using PropagatorOptions = Acts::PropagatorOptions<>;
  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  using VertexFitter =
      Acts::FullBilloirVertexFitter<Acts::BoundTrackParameters, Linearizer>;
  using VertexFitterOptions =
      Acts::VertexingOptions<Acts::BoundTrackParameters>;

  // Setup the magnetic field
  MagneticField bField(m_cfg.bField);
  // Setup the propagator with void navigator
  auto propagator = std::make_shared<Propagator>(Stepper(bField));
  PropagatorOptions propagatorOpts(ctx.geoContext, ctx.magFieldContext,
                                   Acts::LoggerWrapper{logger()});
  // Setup the vertex fitter
  VertexFitter::Config vertexFitterCfg;
  VertexFitter vertexFitter(vertexFitterCfg);
  VertexFitter::State state(ctx.magFieldContext);
  // Setup the linearizer
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  const auto& trajectories =
      ctx.eventStore.get<TrajectoriesContainer>(m_cfg.inputTrajectories);
  
  
  std::vector<Acts::BoundTrackParameters> trackParameters;

   /* // Get the majority truth particle for this trajectory
    const auto particleHitCount = traj.identifyMajorityParticle(trackTip);
    if (particleHitCount.empty()) {
      ACTS_WARNING("No truth particle associated with this trajectory.");
      continue;
    }
    // Find the truth particle for the majority barcode
    const auto ip = particles.find(particleHitCount.front().particleId);
    if (ip == particles.end()) {
      ACTS_WARNING("Majority particle not found in the particles collection.");
      continue;
    }

    // Record this majority particle ID of this trajectory
    reconParticleIds.push_back(ip->particleId());
    // Fill the residual plots
    m_resPlotTool.fill(m_resPlotCache, ctx.geoContext, *ip,
                       traj.trackParameters(trackTip));*/
                       
  int traj_counter = -1;
  
  std::vector<int> empty_traj;
  empty_traj.clear();
  
  
  for (size_t itraj = 0; itraj < trajectories.size(); ++itraj) {
    const auto& traj = trajectories[itraj];
  
     traj_counter += 1;
     std::cout << "traj " << std::endl;
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
       
       
     }
     /*std::cout << "------------------------------------" << std::endl;
     std::vector<double> chi2_arr;
     mj.visitBackwards(trackTip, [&](const auto& state1) {
       std::cout << state1.chi2() << std::endl;
       chi2_arr.push_back(state1.chi2());
     });    
     // cut on chi2 ? 
     if (chi2_arr[2] > 20) continue;

     */
     //std::cout <<  "cov    " << traj.trackParameters(trackTip).covariance().empty() << std::endl;
     //const auto& fittedParameters = traj.trackParameters(trackTip);

     trackParameters.push_back(traj.trackParameters(trackTip));
   
  }

  
  for(auto btr_ : empty_traj) {
     std::cout << "empty    " << btr_ << std::endl;
  }
   
   
      
      
  const auto& protoVertices =
      ctx.eventStore.get<ProtoVertexContainer>(m_cfg.inputProtoVertices);
      
  std::vector<const Acts::BoundTrackParameters*> inputTrackPtrCollection;

  std::vector<Acts::Vertex<Acts::BoundTrackParameters>> outputVertexCollection;
  
  outputVertexCollection.clear();

  std::cout << "proto len and empty len   " << protoVertices.size() << " " << empty_traj.size()  << std::endl;   

  for (const auto& protoVertex : protoVertices) {
    // un-constrained fit requires at least two tracks
    if ((not m_cfg.doConstrainedFit) and protoVertex.size() < 2) {
      ACTS_WARNING(
          "Skip un-constrained vertex fit on proto-vertex with less than two "
          "tracks");
      continue;
    }

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
      inputTrackPtrCollection.push_back(&trackParameters[track_index_not_empy]);     
    }
    
    std::cout << "traj after proto loop 1 " << std::endl;
    Acts::Vertex<Acts::BoundTrackParameters> fittedVertex;

    std::cout << "traj coll input length " << inputTrackPtrCollection.size() << std::endl;

    if (!m_cfg.doConstrainedFit and inputTrackPtrCollection.size() > 1) { // XXX added inputTrackPtrCollection.size() >1
      VertexFitterOptions vfOptions(ctx.geoContext, ctx.magFieldContext);
      std::cout << "unconst traj after proto loop 2 " << std::endl;
      auto fitRes = vertexFitter.fit(inputTrackPtrCollection, linearizer,
                                     vfOptions, state);
      std::cout << "unconst traj after proto loop 3 " << std::endl;
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
      std::cout << "const traj after proto  loop 5 " << std::endl;
      theConstraint.setCovariance(m_cfg.constraintCov);
      theConstraint.setPosition(m_cfg.constraintPos);
      std::cout << "const traj after proto  loop 6 " << std::endl;
      // Vertex fitter options
      VertexFitterOptions vfOptionsConstr(ctx.geoContext, ctx.magFieldContext,
                                          theConstraint);
      std::cout << "const traj after proto  loop 7 " << std::endl;
      auto fitRes = vertexFitter.fit(inputTrackPtrCollection, linearizer,
                                     vfOptionsConstr, state);
                                     
      std::cout << "const traj after proto  loop 8 " << std::endl;                               
      if (fitRes.ok()) {
        fittedVertex = *fitRes;
        outputVertexCollection.push_back(fittedVertex);
      } else {
        ACTS_ERROR("Error in vertex fit with constraint.");
        ACTS_ERROR(fitRes.error().message());
      }
    }

    ACTS_INFO("Fitted Vertex " << fittedVertex.fullPosition().transpose());
  }

  // store vertices 
  
  //fittedVertex

  ctx.eventStore.add(m_cfg.outputFittedVertices, std::move(outputVertexCollection));

  return ProcessCode::SUCCESS;
}
