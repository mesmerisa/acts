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

#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>

/*// See https://stackoverflow.com/a/8513803/2706707
template <typename... Containers>
auto zip( Containers&&... containers ) 
  -> boost::iterator_range <boost::zip_iterator <decltype( boost::make_tuple( std::begin( containers )... ) )> >
{
  auto zip_begin = boost::make_zip_iterator( boost::make_tuple( std::begin( containers )... ) );
  auto zip_end   = boost::make_zip_iterator( boost::make_tuple( std::end(   containers )... ) );
  return boost::make_iterator_range( zip_begin, zip_end );
}*/

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
   const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  
  const auto& trackCandidates =
      ctx.eventStore.get<ProtoTrackContainer>(m_cfg.inputTrackCandidates);

  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputMeasurements);
//auto inserted = simHits.emplace_hint(simHits.end(), simGeometryId,
//                                            simParticleId, simPos4, simMom4,
//                                             simMom4 + simDelta4, simIndex);
  
  std::vector<std::vector<double>> first_hit_pos;
  std::vector<std::vector<double>> last_hit_pos;
  
  for (auto cand : trackCandidates) {
     //std::cout << "-------------------------------" << std::endl; 
     //std::cout << cand.front() << " " << cand.back() << std::endl;    
     auto first_hit_ind = cand.front();
     auto last_hit_ind = cand.back();
     int counter = 0;
     for (auto sh: simHits) {
         if(counter == first_hit_ind) {
              //std::cout << "first hit, id " << first_hit_ind << " pos: "<< sh.position() << std::endl;
              first_hit_pos.push_back({sh.position()[0], sh.position()[1], sh.position()[2]});
              break;
          }
         counter ++;
     }
     counter = 0;
     for (auto sh: simHits) {
         if(counter == last_hit_ind) {
              last_hit_pos.push_back({sh.position()[0], sh.position()[1], sh.position()[2]});
              //std::cout << "last hit, id " << last_hit_ind << " pos: "<< sh.position() << std::endl;
              break;
         }
         counter ++;
     }
  }
  
  //std::cout << "len of cands: " << trackCandidates.size() << std::endl;
  //std::cout << "len of position array " << last_hit_pos.size() << " " << first_hit_pos.size() << std::endl;

  TrackParametersContainer parameters;
  // XXX parameters.reserve(particles.size());
  parameters.reserve(trackCandidates.size());

  // setup random number generator and standard gaussian
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);
  std::normal_distribution<double> stdNormal(0.0, 1.0);
  
  
  
  
//  for (auto&& [vtxId, vtxParticles] : groupBySecondaryVertex(particles)) {
    // a group contains at least one particle by construction. assume that all
    // particles within the group originate from the same position and use it to
    // as the refernce position for the perigee frame.
    
    auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3D{0., 0., 0.}); // I put the perigee at the origin
    
    for (std::size_t n = 0; n < std::min( first_hit_pos.size(), last_hit_pos.size() ); n++) {
       ////////////////////////////////////////////////////////////////////////////////////////////////// 
       //auto  first_p = first_hit_pos[ n ];
       std::vector<double> first_p = {0,0,0};
       auto  last_p = last_hit_pos[ n ];
       //for(int i=0; i<first_p.size(); ++i) std::cout << first_p[i] << ' ';
       //std::cout <<  std::endl;
       //for(int i=0; i<last_p.size(); ++i) std::cout <<  last_p[i] << ' ';
       //std::cout <<  std::endl;
 
       std::vector<double> vec_cand;
       vec_cand.push_back(last_p[0] - first_p[0]); // << std::endl;
       vec_cand.push_back(last_p[1] - first_p[1]);
       vec_cand.push_back(last_p[2] - first_p[2]);
       
       //for(int i=0; i<vec_cand.size(); ++i) std::cout << vec_cand[i] << ' ';
       //std::cout <<  std::endl;
       
       // calculate phi and theta of the vector pointing of the first to the last hit:
       double phi = std::atan2(vec_cand[1],vec_cand[0]);
       double r = std::sqrt(vec_cand[0]*vec_cand[0] + vec_cand[1]*vec_cand[1] + vec_cand[2]*vec_cand[2]);
       double theta =std::acos(vec_cand[2]/r);
       //std::cout << "r, phi, theta " << r << " " << phi << " " << theta << std::endl;
       double time = 0.0; 
       double pt = 0.4; //0.4; // approx the mean value of the input distr., in GeV
       double p = 3; //3;  // same as above 
       double q = 1; //1;
       ////////////////////////////////////////////////////////////////////////////////////////////////
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
      const auto [newPhi, newTheta] = Acts::detail::normalizePhiTheta(
          phi + sigmaPhi * stdNormal(rng), theta + sigmaTheta * stdNormal(rng));
      params[Acts::eBoundPhi] = newPhi;
      params[Acts::eBoundTheta] = newTheta;
       /*// smear direction angles phi,theta ensuring correct bounds
       const double deltaPhi = sigmaPhi * stdNormal(rng);
       const double deltaTheta = sigmaTheta * stdNormal(rng);
       const auto [newPhi, newTheta] =
          Acts::detail::ensureThetaBounds(phi + deltaPhi, theta + deltaTheta);
       params[Acts::eBoundPhi] = newPhi;
       params[Acts::eBoundTheta] = newTheta;*/
       
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
    
    
    
    /*for (auto zi : zip( first_hit_pos, last_hit_pos )) {
      auto [first_p, last_p] = zi;
      std::cout << "= " << *first_p << ", " << *last_p << ", " << std::endl;
    }*/
    /*for (auto&& [vtxId, vtxParticles] : groupBySecondaryVertex(particles)) {
      for (const auto& particle : vtxParticles) {
          const auto time2 = particle.time();
          const auto phi2 = Acts::VectorHelpers::phi(particle.unitDirection());
          const auto theta2 = Acts::VectorHelpers::theta(particle.unitDirection());
          const auto pt2 = particle.transverseMomentum();
          const auto p2 = particle.absMomentum();
          const auto q2 = particle.charge();
          
          std::cout << "true --------    phi, theta "  << " " << phi2 << " " << theta2 << std::endl;
      }
  }*/
      /*// compute momentum-dependent resolutions
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
    }*/
//  }


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
