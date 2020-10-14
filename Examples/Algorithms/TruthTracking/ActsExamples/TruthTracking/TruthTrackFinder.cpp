// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

using namespace ActsExamples;

TruthTrackFinder::TruthTrackFinder(const Config& cfg, Acts::Logging::Level lvl)
    : BareAlgorithm("TruthTrackFinder", lvl), m_cfg(cfg) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputHitParticlesMap.empty()) {
    throw std::invalid_argument("Missing input hit-particles map collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing output proto tracks collection");
  }
}

ProcessCode TruthTrackFinder::execute(const AlgorithmContext& ctx) const {
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
  
  std::cout << "TruthTrackFinder::execute " <<  std::endl; 
  
  // prepare input collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputHitParticlesMap);
  // compute particle_id -> {hit_id...} map from the
  // hit_id -> {particle_id...} map on the fly.
  const auto& particleHitsMap = invertIndexMultimap(hitParticlesMap);
std::cout << "TruthTrackFinder::execute 1" <<  std::endl; 
  // prepare output collection
  ProtoTrackContainer tracks;
  tracks.reserve(particles.size());

  // create prototracks for all input particles
  for (const auto& particle : particles) {
    std::cout << "TruthTrackfinder... particles loop " << std::endl;
    // find the corresponding hits for this particle
    const auto& hits =
        makeRange(particleHitsMap.equal_range(particle.particleId()));
    // fill hit indices to create the proto track
    ProtoTrack track;
    track.reserve(hits.size());
    for (const auto& hit : hits) {
      std::cout << "TruthTrackfinder... hits loop " << std::endl; 
      std::cout << hit.first << std::endl; 
      std::cout << hit.second << std::endl; 
      std::cout <<  std::endl; 
      track.emplace_back(hit.second);
    }
    // add proto track to the output collection
    tracks.emplace_back(std::move(track));
  }
  std::cout << "TruthTrackFinder::execute 5" <<  std::endl; 
  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(tracks));
  std::cout << "TruthTrackFinder::execute last " <<  std::endl; 
  return ProcessCode::SUCCESS;
}
