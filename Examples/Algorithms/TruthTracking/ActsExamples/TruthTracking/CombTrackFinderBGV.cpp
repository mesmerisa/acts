// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/CombTrackFinderBGV.hpp"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

using namespace ActsExamples;

void print(const std::vector<std::vector<int>>& v) {
  std::cout << "{ ";
  for (const auto& p : v) {
    std::cout << "(";
    for (const auto& e : p) {
      std::cout << e << " ";
    }
    std::cout << ") ";
  }
  std::cout << "}" << std::endl;
}

auto product(const std::vector<std::vector<int>>& lists) {
  std::vector<std::vector<int>> result;
  if (std::find_if(std::begin(lists), std::end(lists), 
    [](auto e) -> bool { return e.size() == 0; }) != std::end(lists)) {
    return result;
  }
  for (auto& e : lists[0]) {
    result.push_back({ e });
  }
  for (size_t i = 1; i < lists.size(); ++i) {
    std::vector<std::vector<int>> temp;
    for (auto& e : result) {
      for (auto f : lists[i]) {
        auto e_tmp = e;
        e_tmp.push_back(f);
        temp.push_back(e_tmp);
      }
    }
    result = temp;
  }
  return result;
}

CombTrackFinderBGV::CombTrackFinderBGV(const Config& cfg, Acts::Logging::Level lvl)
    : BareAlgorithm("CombTrackFinderBGV", lvl), m_cfg(cfg) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing input hit-particles map collection");
  }
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing input source links collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing output proto tracks collection");
  }
}

ProcessCode CombTrackFinderBGV::execute(const AlgorithmContext& ctx) const {
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
  
  //std::cout << "CombTrackFinderBGV::execute .... hitparticlemap: " <<  std::endl; 
  
  
  const auto& sourceLinks =
      ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);
  
  // prepare input collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
      
  std::vector<int> particles_ids;
  for (auto p_ : particles) {
     //selectedParticles.push_back(p_.particleId().value());
     std::cout << "selected " << p_.particleId().value() << std::endl;
     particles_ids.push_back(p_.particleId().value());
     }    
      
      
      
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);
  // compute particle_id -> {hit_id...} map from the
  // hit_id -> {particle_id...} map on the fly.
  const auto& particleHitsMap = invertIndexMultimap(hitParticlesMap);
  //std::cout << "CombTrackFinderBGV::execute 1" <<  std::endl; 
  // prepare output collection
  ProtoTrackContainer tracks;
  tracks.reserve(particles.size());
  
  auto old_geo_id = (*sourceLinks.begin()).geometryId();
  
  std::vector<std::vector<int>> hit_id_arrays;
  std::vector<int> curr_ids;
  
  for (auto mea_ : sourceLinks) {
      //if (mea_ == *sourceLinks.end()) std::cout << "joopdi doodle deeee" << std::endl;
      if (old_geo_id != mea_.geometryId()) { // || mea_ == *sourceLinks.end()) {
         // if the geometry id changes, create a new list. TODO this we should change later to whenever the layer changes or something.
         old_geo_id = mea_.geometryId();
         hit_id_arrays.push_back(curr_ids);
         curr_ids.clear();
      }
      curr_ids.push_back(mea_.index());
      std::cout << "source " << mea_.geometryId() << " " << mea_.index() << std::endl;
      //for (auto ind_ : curr_ids) std::cout << "index added " << ind_ << std::endl;
  }
  hit_id_arrays.push_back(curr_ids);
  
  std::cout << "hit arrays: -------------------------------------- " << std::endl;
  
  for (auto ind_arr : hit_id_arrays) {
      for (auto ind_ : ind_arr) std::cout <<  ind_ << " ";
      std::cout << std::endl;
  }
  std::cout << "-------------------------------------------------- " << std::endl;
  
  
  
  
  for (auto phm : particleHitsMap) {
    if (std::find(particles_ids.begin(), particles_ids.end(), phm.first.subParticle()) != particles_ids.end() ) {
    //std::cout << " in if ...... hits map " << phm.first  << " " << phm.second << std::endl;  
    //std::cout << "hits map " << phm.first.subParticle()  << " " << phm.second << std::endl; 
    //particleHitsMap.erase(phm);
    
    
    // XXX save index of the entries that should be removed and remove afterwards?
    //particleHitsMap.erase(std::find(particleHitsMap.begin(),particleHitsMap.end(),phm));
    
    }
  }
   
  auto hit_id_arr_combs = product(hit_id_arrays);
  print(product(hit_id_arrays));
  
   ////////////////////////////////////////////////// here are the track candidates created
  for (auto hit_comb : hit_id_arr_combs) {
    ProtoTrack track;    
    for (auto hit : hit_comb) {
      track.emplace_back(hit);
      //old_hit = hit;
    }
     tracks.emplace_back(std::move(track));
  }
  
  /////////////////////////////////////////////////////

 
  ///////////////////////////////////////////////////// only true tracks, use this:
  // create prototracks for all input particles
 /* for (const auto& particle : particles) {
    //std::cout << "CombTrackFinderBGV... particles loop " << std::endl;
    // find the corresponding hits for this particle
    const auto& hits =
        makeRange(particleHitsMap.equal_range(particle.particleId()));
    // fill hit indices to create the proto track
    ProtoTrack track;
    track.reserve(hits.size());
    std::cout << "CombTrackFinderBGV... hits loop " << std::endl; 
    for (const auto& hit : hits) {
      
      
      std::cout << hit.first <<  " " <<  hit.second << std::endl; 
      //std::cout <<  std::endl; 
      track.emplace_back(hit.second);
    }
    // add proto track to the output collection
    tracks.emplace_back(std::move(track));
  }*/
  /////////////////////////////////////////////////////
  /*for(auto tr: tracks) {
	std::cout << "track" << std::endl;  
	for(auto hit: tr) {
		std::cout << "hit " << hit << " ";
    }	  
    std::cout << std::endl;
  }	  */
  
  
  //std::cout << "CombTrackFinderBGV::execute 5" <<  std::endl; 
  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(tracks));
  //std::cout << "CombTrackFinderBGV::execute last " <<  std::endl; 
  return ProcessCode::SUCCESS;
}
