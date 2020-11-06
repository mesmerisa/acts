// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/ParticleSmearingBGV.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimHit.hpp"

#include <cmath>
#include <vector>

ActsExamples::ParticleSmearingBGV::ParticleSmearingBGV(const Config& cfg,
                                                 Acts::Logging::Level lvl)
    : BareAlgorithm("ParticleSmearingBGV", lvl), m_cfg(cfg) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputTrackCandidates.empty()) {
    throw std::invalid_argument("Missing input track candidate collection");
  }
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input hit collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing output track parameters collection");
  }
}

ActsExamples::ProcessCode ActsExamples::ParticleSmearingBGV::execute(
    const AlgorithmContext& ctx) const {
  // setup input and output containers
  //XXX const auto& particles =
  //XXX    ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  
  const auto& trackCandidates =
      ctx.eventStore.get<ProtoTrackContainer>(m_cfg.inputTrackCandidates);

  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputMeasurements);
//auto inserted = simHits.emplace_hint(simHits.end(), simGeometryId,
//                                            simParticleId, simPos4, simMom4,
//                                             simMom4 + simDelta4, simIndex);
  
  
  for (auto cand : trackCandidates) {
     std::cout << cand.front() << " " << cand.back() << std::endl;    
     auto first_hit_ind = cand.front();
     auto last_hit_ind = cand.back();
     int counter = 0;
     for (auto sh: simHits) {
         if(counter == first_hit_ind) {
            std::cout << "first hit, id " << first_hit_ind << " pos: "<< sh.position() << std::endl;
         
         }   
         if(counter == last_hit_ind) {
         
         } 
     } 
  }
  
  
  
  /*for (auto sh : simHits) {
   std::cout <<"index: " << sh.index() << " position: " << sh.position() << std::endl;   
      
      
  }*/
  
  /*std::cout << "test position " << *(simHits.begin()) << std::endl;*/

  TrackParametersContainer parameters;
  // XXX parameters.reserve(particles.size());
  parameters.reserve(trackCandidates.size());

  // setup random number generator and standard gaussian
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);
  std::normal_distribution<double> stdNormal(0.0, 1.0);

  /*for (auto&& [vtxId, vtxParticles] : groupBySecondaryVertex(particles)) {
    // a group contains at least one particle by construction. assume that all
    // particles within the group originate from the same position and use it to
    // as the refernce position for the perigee frame.
    auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(
        vtxParticles.begin()->position());

    for (const auto& particle : vtxParticles) {
      const auto time = particle.time();
      const auto phi = Acts::VectorHelpers::phi(particle.unitDirection());
      const auto theta = Acts::VectorHelpers::theta(particle.unitDirection());
      const auto pt = particle.transverseMomentum();
      const auto p = particle.absMomentum();
      const auto q = particle.charge();
      
      std::cout << "mom " << p << std::endl;
      std::cout << "q " << q << std::endl;

      // compute momentum-dependent resolutions
      const double sigmaD0 =
          m_cfg.sigmaD0 +
          m_cfg.sigmaD0PtA * std::exp(-1.0 * std::abs(m_cfg.sigmaD0PtB) * pt);
      const double sigmaZ0 =
          m_cfg.sigmaZ0 +
          m_cfg.sigmaZ0PtA * std::exp(-1.0 * std::abs(m_cfg.sigmaZ0PtB) * pt);
      const double sigmaP = m_cfg.sigmaPRel * p;
      // var(q/p) = (d(1/p)/dp)² * var(p) = (-1/p²)² * var(p)
      const double sigmaQOverP = sigmaP / (p * p);
      // shortcuts for other resolutions
      const double sigmaT0 = m_cfg.sigmaT0;
      const double sigmaPhi = m_cfg.sigmaPhi;
      const double sigmaTheta = m_cfg.sigmaTheta;

      Acts::BoundVector params = Acts::BoundVector::Zero();
      // smear the position/time
      params[Acts::eBoundLoc0] = sigmaD0 * stdNormal(rng);
      params[Acts::eBoundLoc1] = sigmaZ0 * stdNormal(rng);
      params[Acts::eBoundTime] = time + sigmaT0 * stdNormal(rng);
      // smear direction angles phi,theta ensuring correct bounds
      const double deltaPhi = sigmaPhi * stdNormal(rng);
      const double deltaTheta = sigmaTheta * stdNormal(rng);
      const auto [newPhi, newTheta] =
          Acts::detail::ensureThetaBounds(phi + deltaPhi, theta + deltaTheta);
      params[Acts::eBoundPhi] = newPhi;
      params[Acts::eBoundTheta] = newTheta;
      // compute smeared absolute momentum vector
      const double newP = std::max(0.0, p + sigmaP * stdNormal(rng));
      params[Acts::eBoundQOverP] = (q != 0) ? (q / newP) : (1 / newP);

      // build the track covariance matrix using the smearing sigmas
      Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();
      cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = sigmaD0 * sigmaD0;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = sigmaZ0 * sigmaZ0;
      cov(Acts::eBoundTime, Acts::eBoundTime) = sigmaT0 * sigmaT0;
      cov(Acts::eBoundPhi, Acts::eBoundPhi) = sigmaPhi * sigmaPhi;
      cov(Acts::eBoundTheta, Acts::eBoundTheta) = sigmaTheta * sigmaTheta;
      cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = sigmaQOverP * sigmaQOverP;

      parameters.emplace_back(perigee, params, q, cov);
    }
  }
*/
  ctx.eventStore.add(m_cfg.outputTrackParameters, std::move(parameters));
  return ProcessCode::SUCCESS;
}
