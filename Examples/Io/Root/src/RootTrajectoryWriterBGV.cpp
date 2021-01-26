// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootTrajectoryWriterBGV.hpp"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

#include "detail/AverageSimHits.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::RootTrajectoryWriterBGV::RootTrajectoryWriterBGV(
    const ActsExamples::RootTrajectoryWriterBGV::Config& cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputTrajectories, "RootTrajectoryWriterBGV", lvl),
      m_cfg(cfg),
      m_outputFile(cfg.rootFile) {
  // trajectories collection name is already checked by base ctor
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements input collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-simulated-hits map input collection");
  }
  if (m_cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.outputTreename.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    auto path = joinPaths(m_cfg.outputDir, m_cfg.outputFilename);
    m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + path);
    }
  }
  m_outputFile->cd();
  m_outputTree =
      new TTree(m_cfg.outputTreename.c_str(), m_cfg.outputTreename.c_str());
  m_outputTree_wrong =
      new TTree(m_cfg.outputTreename_wrong.c_str(), m_cfg.outputTreename_wrong.c_str());    
      
  if (m_outputTree == nullptr)
    throw std::bad_alloc();
  else {
    // I/O parameters for "wrong" tree  
    m_outputTree_wrong->Branch("event_nr", &m_eventNr_wrong);
    m_outputTree_wrong->Branch("traj_nr", &m_trajNr_wrong);  
    m_outputTree_wrong->Branch("true_or_not", &true_or_not_wrong);
    m_outputTree_wrong->Branch("chi2", &m_chi2_wrong);
      
    // I/O parameters
    m_outputTree->Branch("event_nr", &m_eventNr);
    m_outputTree->Branch("traj_nr", &m_trajNr);
    m_outputTree->Branch("t_barcode", &m_t_barcode, "t_barcode/l");
    m_outputTree->Branch("t_charge", &m_t_charge);
    m_outputTree->Branch("t_time", &m_t_time);
    m_outputTree->Branch("t_vx", &m_t_vx);
    m_outputTree->Branch("t_vy", &m_t_vy);
    m_outputTree->Branch("t_vz", &m_t_vz);
    m_outputTree->Branch("t_px", &m_t_px);
    m_outputTree->Branch("t_py", &m_t_py);
    m_outputTree->Branch("t_pz", &m_t_pz);
    m_outputTree->Branch("t_theta", &m_t_theta);
    m_outputTree->Branch("t_phi", &m_t_phi);
    m_outputTree->Branch("t_eta", &m_t_eta);
    m_outputTree->Branch("t_pT", &m_t_pT);
    
    m_outputTree->Branch("true_or_not", &true_or_not);

    m_outputTree->Branch("t_x", &m_t_x);
    m_outputTree->Branch("t_y", &m_t_y);
    m_outputTree->Branch("t_z", &m_t_z);
    m_outputTree->Branch("t_r", &m_t_r);
    m_outputTree->Branch("t_dx", &m_t_dx);
    m_outputTree->Branch("t_dy", &m_t_dy);
    m_outputTree->Branch("t_dz", &m_t_dz);
    m_outputTree->Branch("t_eLOC0", &m_t_eLOC0);
    m_outputTree->Branch("t_eLOC1", &m_t_eLOC1);
    m_outputTree->Branch("t_ePHI", &m_t_ePHI);
    m_outputTree->Branch("t_eTHETA", &m_t_eTHETA);
    m_outputTree->Branch("t_eQOP", &m_t_eQOP);
    m_outputTree->Branch("t_eT", &m_t_eT);

    m_outputTree->Branch("nStates", &m_nStates);
    m_outputTree->Branch("nMeasurements", &m_nMeasurements);
    m_outputTree->Branch("volume_id", &m_volumeID);
    m_outputTree->Branch("layer_id", &m_layerID);
    m_outputTree->Branch("module_id", &m_moduleID);
    m_outputTree->Branch("l_x_hit", &m_lx_hit);
    m_outputTree->Branch("l_y_hit", &m_ly_hit);
    m_outputTree->Branch("g_x_hit", &m_x_hit);
    m_outputTree->Branch("g_y_hit", &m_y_hit);
    m_outputTree->Branch("g_z_hit", &m_z_hit);
    m_outputTree->Branch("res_x_hit", &m_res_x_hit);
    m_outputTree->Branch("res_y_hit", &m_res_y_hit);
    m_outputTree->Branch("err_x_hit", &m_err_x_hit);
    m_outputTree->Branch("err_y_hit", &m_err_y_hit);
    m_outputTree->Branch("pull_x_hit", &m_pull_x_hit);
    m_outputTree->Branch("pull_y_hit", &m_pull_y_hit);
    m_outputTree->Branch("dim_hit", &m_dim_hit);

    m_outputTree->Branch("hasFittedParams", &m_hasFittedParams);
    m_outputTree->Branch("eLOC0_fit", &m_eLOC0_fit);
    m_outputTree->Branch("eLOC1_fit", &m_eLOC1_fit);
    m_outputTree->Branch("ePHI_fit", &m_ePHI_fit);
    m_outputTree->Branch("eTHETA_fit", &m_eTHETA_fit);
    m_outputTree->Branch("eQOP_fit", &m_eQOP_fit);
    m_outputTree->Branch("eT_fit", &m_eT_fit);
    m_outputTree->Branch("err_eLOC0_fit", &m_err_eLOC0_fit);
    m_outputTree->Branch("err_eLOC1_fit", &m_err_eLOC1_fit);
    m_outputTree->Branch("err_ePHI_fit", &m_err_ePHI_fit);
    m_outputTree->Branch("err_eTHETA_fit", &m_err_eTHETA_fit);
    m_outputTree->Branch("err_eQOP_fit", &m_err_eQOP_fit);
    m_outputTree->Branch("err_eT_fit", &m_err_eT_fit);

    m_outputTree->Branch("nPredicted", &m_nPredicted);
    m_outputTree->Branch("predicted", &m_prt);
    m_outputTree->Branch("eLOC0_prt", &m_eLOC0_prt);
    m_outputTree->Branch("eLOC1_prt", &m_eLOC1_prt);
    m_outputTree->Branch("ePHI_prt", &m_ePHI_prt);
    m_outputTree->Branch("eTHETA_prt", &m_eTHETA_prt);
    m_outputTree->Branch("eQOP_prt", &m_eQOP_prt);
    m_outputTree->Branch("eT_prt", &m_eT_prt);
    m_outputTree->Branch("res_eLOC0_prt", &m_res_eLOC0_prt);
    m_outputTree->Branch("res_eLOC1_prt", &m_res_eLOC1_prt);
    m_outputTree->Branch("res_ePHI_prt", &m_res_ePHI_prt);
    m_outputTree->Branch("res_eTHETA_prt", &m_res_eTHETA_prt);
    m_outputTree->Branch("res_eQOP_prt", &m_res_eQOP_prt);
    m_outputTree->Branch("res_eT_prt", &m_res_eT_prt);
    m_outputTree->Branch("err_eLOC0_prt", &m_err_eLOC0_prt);
    m_outputTree->Branch("err_eLOC1_prt", &m_err_eLOC1_prt);
    m_outputTree->Branch("err_ePHI_prt", &m_err_ePHI_prt);
    m_outputTree->Branch("err_eTHETA_prt", &m_err_eTHETA_prt);
    m_outputTree->Branch("err_eQOP_prt", &m_err_eQOP_prt);
    m_outputTree->Branch("err_eT_prt", &m_err_eT_prt);
    m_outputTree->Branch("pull_eLOC0_prt", &m_pull_eLOC0_prt);
    m_outputTree->Branch("pull_eLOC1_prt", &m_pull_eLOC1_prt);
    m_outputTree->Branch("pull_ePHI_prt", &m_pull_ePHI_prt);
    m_outputTree->Branch("pull_eTHETA_prt", &m_pull_eTHETA_prt);
    m_outputTree->Branch("pull_eQOP_prt", &m_pull_eQOP_prt);
    m_outputTree->Branch("pull_eT_prt", &m_pull_eT_prt);
    m_outputTree->Branch("g_x_prt", &m_x_prt);
    m_outputTree->Branch("g_y_prt", &m_y_prt);
    m_outputTree->Branch("g_z_prt", &m_z_prt);
    m_outputTree->Branch("px_prt", &m_px_prt);
    m_outputTree->Branch("py_prt", &m_py_prt);
    m_outputTree->Branch("pz_prt", &m_pz_prt);
    m_outputTree->Branch("eta_prt", &m_eta_prt);
    m_outputTree->Branch("pT_prt", &m_pT_prt);

    m_outputTree->Branch("nFiltered", &m_nFiltered);
    m_outputTree->Branch("filtered", &m_flt);
    m_outputTree->Branch("eLOC0_flt", &m_eLOC0_flt);
    m_outputTree->Branch("eLOC1_flt", &m_eLOC1_flt);
    m_outputTree->Branch("ePHI_flt", &m_ePHI_flt);
    m_outputTree->Branch("eTHETA_flt", &m_eTHETA_flt);
    m_outputTree->Branch("eQOP_flt", &m_eQOP_flt);
    m_outputTree->Branch("eT_flt", &m_eT_flt);
    m_outputTree->Branch("res_eLOC0_flt", &m_res_eLOC0_flt);
    m_outputTree->Branch("res_eLOC1_flt", &m_res_eLOC1_flt);
    m_outputTree->Branch("res_ePHI_flt", &m_res_ePHI_flt);
    m_outputTree->Branch("res_eTHETA_flt", &m_res_eTHETA_flt);
    m_outputTree->Branch("res_eQOP_flt", &m_res_eQOP_flt);
    m_outputTree->Branch("res_eT_flt", &m_res_eT_flt);
    m_outputTree->Branch("err_eLOC0_flt", &m_err_eLOC0_flt);
    m_outputTree->Branch("err_eLOC1_flt", &m_err_eLOC1_flt);
    m_outputTree->Branch("err_ePHI_flt", &m_err_ePHI_flt);
    m_outputTree->Branch("err_eTHETA_flt", &m_err_eTHETA_flt);
    m_outputTree->Branch("err_eQOP_flt", &m_err_eQOP_flt);
    m_outputTree->Branch("err_eT_flt", &m_err_eT_flt);
    m_outputTree->Branch("pull_eLOC0_flt", &m_pull_eLOC0_flt);
    m_outputTree->Branch("pull_eLOC1_flt", &m_pull_eLOC1_flt);
    m_outputTree->Branch("pull_ePHI_flt", &m_pull_ePHI_flt);
    m_outputTree->Branch("pull_eTHETA_flt", &m_pull_eTHETA_flt);
    m_outputTree->Branch("pull_eQOP_flt", &m_pull_eQOP_flt);
    m_outputTree->Branch("pull_eT_flt", &m_pull_eT_flt);
    m_outputTree->Branch("g_x_flt", &m_x_flt);
    m_outputTree->Branch("g_y_flt", &m_y_flt);
    m_outputTree->Branch("g_z_flt", &m_z_flt);
    m_outputTree->Branch("px_flt", &m_px_flt);
    m_outputTree->Branch("py_flt", &m_py_flt);
    m_outputTree->Branch("pz_flt", &m_pz_flt);
    m_outputTree->Branch("eta_flt", &m_eta_flt);
    m_outputTree->Branch("pT_flt", &m_pT_flt);
    m_outputTree->Branch("chi2", &m_chi2);

    m_outputTree->Branch("nSmoothed", &m_nSmoothed);
    m_outputTree->Branch("smoothed", &m_smt);
    m_outputTree->Branch("eLOC0_smt", &m_eLOC0_smt);
    m_outputTree->Branch("eLOC1_smt", &m_eLOC1_smt);
    m_outputTree->Branch("ePHI_smt", &m_ePHI_smt);
    m_outputTree->Branch("eTHETA_smt", &m_eTHETA_smt);
    m_outputTree->Branch("eQOP_smt", &m_eQOP_smt);
    m_outputTree->Branch("eT_smt", &m_eT_smt);
    m_outputTree->Branch("res_eLOC0_smt", &m_res_eLOC0_smt);
    m_outputTree->Branch("res_eLOC1_smt", &m_res_eLOC1_smt);
    m_outputTree->Branch("res_ePHI_smt", &m_res_ePHI_smt);
    m_outputTree->Branch("res_eTHETA_smt", &m_res_eTHETA_smt);
    m_outputTree->Branch("res_eQOP_smt", &m_res_eQOP_smt);
    m_outputTree->Branch("res_eT_smt", &m_res_eT_smt);
    m_outputTree->Branch("err_eLOC0_smt", &m_err_eLOC0_smt);
    m_outputTree->Branch("err_eLOC1_smt", &m_err_eLOC1_smt);
    m_outputTree->Branch("err_ePHI_smt", &m_err_ePHI_smt);
    m_outputTree->Branch("err_eTHETA_smt", &m_err_eTHETA_smt);
    m_outputTree->Branch("err_eQOP_smt", &m_err_eQOP_smt);
    m_outputTree->Branch("err_eT_smt", &m_err_eT_smt);
    m_outputTree->Branch("pull_eLOC0_smt", &m_pull_eLOC0_smt);
    m_outputTree->Branch("pull_eLOC1_smt", &m_pull_eLOC1_smt);
    m_outputTree->Branch("pull_ePHI_smt", &m_pull_ePHI_smt);
    m_outputTree->Branch("pull_eTHETA_smt", &m_pull_eTHETA_smt);
    m_outputTree->Branch("pull_eQOP_smt", &m_pull_eQOP_smt);
    m_outputTree->Branch("pull_eT_smt", &m_pull_eT_smt);
    m_outputTree->Branch("g_x_smt", &m_x_smt);
    m_outputTree->Branch("g_y_smt", &m_y_smt);
    m_outputTree->Branch("g_z_smt", &m_z_smt);
    m_outputTree->Branch("px_smt", &m_px_smt);
    m_outputTree->Branch("py_smt", &m_py_smt);
    m_outputTree->Branch("pz_smt", &m_pz_smt);
    m_outputTree->Branch("eta_smt", &m_eta_smt);
    m_outputTree->Branch("pT_smt", &m_pT_smt);
  }
}

ActsExamples::RootTrajectoryWriterBGV::~RootTrajectoryWriterBGV() {
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootTrajectoryWriterBGV::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
    ACTS_INFO("Write trajectories to tree '"
              << m_cfg.outputTreename << "' in '"
              << joinPaths(m_cfg.outputDir, m_cfg.outputFilename) << "'");
    m_outputTree_wrong->Write();
    ACTS_INFO("Write trajectories to tree '"
              << m_cfg.outputTreename_wrong << "' in '"
              << joinPaths(m_cfg.outputDir, m_cfg.outputFilename) << "'");          
  }
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootTrajectoryWriterBGV::writeT(
    const AlgorithmContext& ctx, const TrajectoriesContainer& trajectories) {
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
  using HitSimHitsMap = IndexMultimap<Index>;
  using ConcreteMeasurement =
      Acts::Measurement<ActsExamples::IndexSourceLink, Acts::BoundIndices,
                        Acts::eBoundLoc0>; //, Acts::eBoundLoc1>;

  if (m_outputFile == nullptr)
    return ProcessCode::SUCCESS;

  auto& gctx = ctx.geoContext;
  // Read additional input collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);
  const auto& measurements =
      ctx.eventStore.get<MeasurementContainer>(m_cfg.inputMeasurements);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);
  const auto& hitSimHitsMap =
      ctx.eventStore.get<HitSimHitsMap>(m_cfg.inputMeasurementSimHitsMap);

  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;
  
  
  //for (auto hitpmap : hitParticlesMap)          std::cout << hitpmap.first << std::endl;  
  //std::cout << std::endl;                      

  //std::cout << "#################################################################### number of trajectories " << trajectories.size() << std::endl;

  // Loop over the trajectories
  int counter_else = 0;
  int counter_if = 0;
  for (size_t itraj = 0; itraj < trajectories.size(); ++itraj) {
      
    //std::cout << "START LOOP!" << std::endl;    
      
    const auto& traj = trajectories[itraj];

    if (traj.empty()) {
      ACTS_WARNING("Empty trajectories object " << itraj);
      continue;
    }

    m_trajNr = itraj;

    // The trajectory entry indices and the multiTrajectory
    const auto& mj = traj.multiTrajectory();
    const auto& trackTips = traj.tips();
    // Check the size of the trajectory entry indices. For track fitting, there
    // should be at most one trajectory
    if (trackTips.size() > 1) {
      ACTS_ERROR("Track fitting should not result in multiple trajectories.");
      return ProcessCode::ABORT;
    }

    // Get the entry index for the single trajectory
    auto trackTip = trackTips.front();
    // Collect the trajectory summary info
    auto trajState =
        Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
    m_nMeasurements = trajState.nMeasurements;
    m_nStates = trajState.nStates;

    // Get the majority truth particle to this track
    identifyContributingParticles(hitParticlesMap, traj, trackTip,
                                  particleHitCounts);
                                  
     std::vector<int> selectedParticles;
     selectedParticles.clear();
     
     for (auto p_ : particles) selectedParticles.push_back(p_.particleId().value());
                
     true_or_not = 0;            
                     
     if (particleHitCounts.size() == 1 &&  std::find(selectedParticles.begin(), selectedParticles.end(), particleHitCounts.front().particleId.value()) != selectedParticles.end() ) { //&& inTruthParticleColl) { // if only 1 particle contributes - these are the true tracks.

      true_or_not = 1;
      //std::cout << "true!!!!" << std::endl;   
      // Get the barcode of the majority truth particle
      m_t_barcode = particleHitCounts.front().particleId.value();
      // Find the truth particle via the barcode
      auto ip = particles.find(m_t_barcode);
      if (ip != particles.end()) {
        const auto& particle = *ip;
        ACTS_DEBUG("Find the truth particle with barcode = " << m_t_barcode);
        // Get the truth particle info at vertex
        const auto p = particle.absoluteMomentum();
        m_t_charge = particle.charge();
        m_t_time = particle.time();
        m_t_vx = particle.position().x();
        m_t_vy = particle.position().y();
        m_t_vz = particle.position().z();
        m_t_px = p * particle.unitDirection().x();
        m_t_py = p * particle.unitDirection().y();
        m_t_pz = p * particle.unitDirection().z();
        m_t_theta = theta(particle.unitDirection());
        m_t_phi = phi(particle.unitDirection());
        m_t_eta = eta(particle.unitDirection());
        m_t_pT = p * perp(particle.unitDirection());
        
        counter_if ++;
      } else {
        ACTS_WARNING("Truth particle with barcode = " << m_t_barcode
                                                      << " not found!");
      }
    }
    else {
    true_or_not = 0;
    counter_else ++;      
    //std::cout << "WRONG TRACK! CONTINUE!" << std::endl;
    m_eventNr_wrong = ctx.eventNumber;
    m_trajNr_wrong = itraj;
    true_or_not_wrong = 0;
    mj.visitBackwards(trackTip, [&](const auto& state) {
        if (state.hasFiltered()) m_chi2_wrong.push_back(state.chi2());
    });    
    m_outputTree_wrong->Fill();
    m_chi2_wrong.clear();
    continue;
    }    
    //std::cout << "AFTER CONTINUE!" << std::endl;

    // Get the fitted track parameter
    m_hasFittedParams = false;
    if (traj.hasTrackParameters(trackTip)) {
      m_hasFittedParams = true;
      const auto& boundParam = traj.trackParameters(trackTip);
      const auto& parameter = boundParam.parameters();
      const auto& covariance = *boundParam.covariance();
      m_eLOC0_fit = parameter[Acts::eBoundLoc0];
      m_eLOC1_fit = parameter[Acts::eBoundLoc1];
      m_ePHI_fit = parameter[Acts::eBoundPhi];
      m_eTHETA_fit = parameter[Acts::eBoundTheta];
      m_eQOP_fit = parameter[Acts::eBoundQOverP];
      m_eT_fit = parameter[Acts::eBoundTime];
      m_err_eLOC0_fit = sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0));
      m_err_eLOC1_fit = sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1));
      m_err_ePHI_fit = sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi));
      m_err_eTHETA_fit = sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta));
      m_err_eQOP_fit = sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP));
      m_err_eT_fit = sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime));
    }

    // Get the trackStates on the trajectory
    m_nPredicted = 0;
    m_nFiltered = 0;
    m_nSmoothed = 0;
    mj.visitBackwards(trackTip, [&](const auto& state) {
      // we only fill the track states with non-outlier measurement
      auto typeFlags = state.typeFlags();
      if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
        return true;
      }

      const auto hitIdx = state.uncalibrated().index();
      const auto& fullMeas = measurements[hitIdx];
      const auto& meas = std::get<ConcreteMeasurement>(fullMeas);
      const auto& surface = meas.referenceObject();

      // get the geometry ID
      auto geoID = surface.geometryId();
      m_volumeID.push_back(geoID.volume());
      m_layerID.push_back(geoID.layer());
      m_moduleID.push_back(geoID.sensitive());

      // get local position
      Acts::Vector2 local(meas.parameters()[Acts::eBoundLoc0],
                           meas.parameters()[Acts::eBoundLoc1]);
      // get global position
      Acts::Vector3 mom(1, 1, 1);
      Acts::Vector3 global = surface.localToGlobal(ctx.geoContext, local, mom);

      // get measurement covariance
      auto cov = meas.covariance();
      // float resX = sqrt(cov(Acts::eBoundLoc0, Acts::eBoundLoc0));
      // float resY = sqrt(cov(Acts::eBoundLoc1, Acts::eBoundLoc1));

      // push the measurement info
      m_lx_hit.push_back(local.x());
      m_ly_hit.push_back(local.y());
      m_x_hit.push_back(global.x());
      m_y_hit.push_back(global.y());
      m_z_hit.push_back(global.z());

      // get the truth hits corresponding to this trackState
      // Use average truth in the case of multiple contributing sim hits
      auto indices = makeRange(hitSimHitsMap.equal_range(hitIdx));
      auto [truthLocal, truthPos4, truthUnitDir] =
          detail::averageSimHits(ctx.geoContext, surface, simHits, indices);
      // momemtum averaging makes even less sense than averaging position and
      // direction. use the first momentum or set q/p to zero
      float truthQOP = 0.0f;
      if (not indices.empty()) {
        // we assume that the indices are within valid ranges so we do not need
        // to check their validity again.
        const auto simHitIdx0 = indices.begin()->second;
        const auto& simHit0 = *simHits.nth(simHitIdx0);
        const auto p =
            simHit0.momentum4Before().template segment<3>(Acts::eMom0).norm();
        truthQOP = m_t_charge / p;
      }

      // push the truth hit info
      m_t_x.push_back(truthPos4[Acts::ePos0]);
      m_t_y.push_back(truthPos4[Acts::ePos1]);
      m_t_z.push_back(truthPos4[Acts::ePos2]);
      m_t_r.push_back(perp(truthPos4.template segment<3>(Acts::ePos0)));
      m_t_dx.push_back(truthUnitDir[Acts::eMom0]);
      m_t_dy.push_back(truthUnitDir[Acts::eMom1]);
      m_t_dz.push_back(truthUnitDir[Acts::eMom1]);

      // get the truth track parameter at this track State
      float truthLOC0 = truthLocal[Acts::ePos0];
      float truthLOC1 = truthLocal[Acts::ePos1];
      float truthTIME = truthPos4[Acts::eTime];
      float truthPHI = phi(truthUnitDir);
      float truthTHETA = theta(truthUnitDir);

      // push the truth track parameter at this track State
      m_t_eLOC0.push_back(truthLOC0);
      m_t_eLOC1.push_back(truthLOC1);
      m_t_ePHI.push_back(truthPHI);
      m_t_eTHETA.push_back(truthTHETA);
      m_t_eQOP.push_back(truthQOP);
      m_t_eT.push_back(truthTIME);

      // get the predicted parameter
      bool predicted = false;
      if (state.hasPredicted()) {
        predicted = true;
        m_nPredicted++;
        auto parameters = state.predicted();
        auto covariance = state.predictedCovariance();
        // local hit residual info
        auto H = meas.projector();
        auto resCov = cov + H * covariance * H.transpose();
        auto residual = meas.residual(parameters);
        m_res_x_hit.push_back(residual(Acts::eBoundLoc0));
        m_res_y_hit.push_back(residual(Acts::eBoundLoc1));
        m_err_x_hit.push_back(sqrt(resCov(Acts::eBoundLoc0, Acts::eBoundLoc0)));
        m_err_y_hit.push_back(sqrt(resCov(Acts::eBoundLoc1, Acts::eBoundLoc1)));
        m_pull_x_hit.push_back(
            residual(Acts::eBoundLoc0) /
            sqrt(resCov(Acts::eBoundLoc0, Acts::eBoundLoc0)));
        m_pull_y_hit.push_back(
            residual(Acts::eBoundLoc1) /
            sqrt(resCov(Acts::eBoundLoc1, Acts::eBoundLoc1)));
        m_dim_hit.push_back(state.calibratedSize());

        // predicted parameter
        m_eLOC0_prt.push_back(parameters[Acts::eBoundLoc0]);
        m_eLOC1_prt.push_back(parameters[Acts::eBoundLoc1]);
        m_ePHI_prt.push_back(parameters[Acts::eBoundPhi]);
        m_eTHETA_prt.push_back(parameters[Acts::eBoundTheta]);
        m_eQOP_prt.push_back(parameters[Acts::eBoundQOverP]);
        m_eT_prt.push_back(parameters[Acts::eBoundTime]);

        // predicted residual
        m_res_eLOC0_prt.push_back(parameters[Acts::eBoundLoc0] - truthLOC0);
        m_res_eLOC1_prt.push_back(parameters[Acts::eBoundLoc1] - truthLOC1);
        m_res_ePHI_prt.push_back(parameters[Acts::eBoundPhi] - truthPHI);
        m_res_eTHETA_prt.push_back(parameters[Acts::eBoundTheta] - truthTHETA);
        m_res_eQOP_prt.push_back(parameters[Acts::eBoundQOverP] - truthQOP);
        m_res_eT_prt.push_back(parameters[Acts::eBoundTime] - truthTIME);

        // predicted parameter error
        m_err_eLOC0_prt.push_back(
            sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
        m_err_eLOC1_prt.push_back(
            sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
        m_err_ePHI_prt.push_back(
            sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
        m_err_eTHETA_prt.push_back(
            sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
        m_err_eQOP_prt.push_back(
            sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
        m_err_eT_prt.push_back(
            sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

        // predicted parameter pull
        m_pull_eLOC0_prt.push_back(
            (parameters[Acts::eBoundLoc0] - truthLOC0) /
            sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
        m_pull_eLOC1_prt.push_back(
            (parameters[Acts::eBoundLoc1] - truthLOC1) /
            sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
        m_pull_ePHI_prt.push_back(
            (parameters[Acts::eBoundPhi] - truthPHI) /
            sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
        m_pull_eTHETA_prt.push_back(
            (parameters[Acts::eBoundTheta] - truthTHETA) /
            sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
        m_pull_eQOP_prt.push_back(
            (parameters[Acts::eBoundQOverP] - truthQOP) /
            sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
        m_pull_eT_prt.push_back(
            (parameters[Acts::eBoundTime] - truthTIME) /
            sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

        // further predicted parameter info
        Acts::FreeVector freeParams =
            Acts::detail::transformBoundToFreeParameters(surface, gctx,
                                                         parameters);
        m_x_prt.push_back(freeParams[Acts::eFreePos0]);
        m_y_prt.push_back(freeParams[Acts::eFreePos1]);
        m_z_prt.push_back(freeParams[Acts::eFreePos2]);
        auto p = std::abs(1 / freeParams[Acts::eFreeQOverP]);
        m_px_prt.push_back(p * freeParams[Acts::eFreeDir0]);
        m_py_prt.push_back(p * freeParams[Acts::eFreeDir1]);
        m_pz_prt.push_back(p * freeParams[Acts::eFreeDir2]);
        m_pT_prt.push_back(p * std::hypot(freeParams[Acts::eFreeDir0],
                                          freeParams[Acts::eFreeDir1]));
        m_eta_prt.push_back(
            Acts::VectorHelpers::eta(freeParams.segment<3>(Acts::eFreeDir0)));
      } else {
        // push default values if no predicted parameter
        m_res_x_hit.push_back(-99.);
        m_res_y_hit.push_back(-99.);
        m_err_x_hit.push_back(-99.);
        m_err_y_hit.push_back(-99.);
        m_pull_x_hit.push_back(-99.);
        m_pull_y_hit.push_back(-99.);
        m_dim_hit.push_back(-99.);
        m_eLOC0_prt.push_back(-99.);
        m_eLOC1_prt.push_back(-99.);
        m_ePHI_prt.push_back(-99.);
        m_eTHETA_prt.push_back(-99.);
        m_eQOP_prt.push_back(-99.);
        m_eT_prt.push_back(-99.);
        m_res_eLOC0_prt.push_back(-99.);
        m_res_eLOC1_prt.push_back(-99.);
        m_res_ePHI_prt.push_back(-99.);
        m_res_eTHETA_prt.push_back(-99.);
        m_res_eQOP_prt.push_back(-99.);
        m_res_eT_prt.push_back(-99.);
        m_err_eLOC0_prt.push_back(-99);
        m_err_eLOC1_prt.push_back(-99);
        m_err_ePHI_prt.push_back(-99);
        m_err_eTHETA_prt.push_back(-99);
        m_err_eQOP_prt.push_back(-99);
        m_err_eT_prt.push_back(-99);
        m_pull_eLOC0_prt.push_back(-99.);
        m_pull_eLOC1_prt.push_back(-99.);
        m_pull_ePHI_prt.push_back(-99.);
        m_pull_eTHETA_prt.push_back(-99.);
        m_pull_eQOP_prt.push_back(-99.);
        m_pull_eT_prt.push_back(-99.);
        m_x_prt.push_back(-99.);
        m_y_prt.push_back(-99.);
        m_z_prt.push_back(-99.);
        m_px_prt.push_back(-99.);
        m_py_prt.push_back(-99.);
        m_pz_prt.push_back(-99.);
        m_pT_prt.push_back(-99.);
        m_eta_prt.push_back(-99.);
      }

      // get the filtered parameter
      bool filtered = false;
      if (state.hasFiltered()) {
        filtered = true;
        m_nFiltered++;
        auto parameters = state.filtered();
        auto covariance = state.filteredCovariance();
        // filtered parameter
        m_eLOC0_flt.push_back(parameters[Acts::eBoundLoc0]);
        m_eLOC1_flt.push_back(parameters[Acts::eBoundLoc1]);
        m_ePHI_flt.push_back(parameters[Acts::eBoundPhi]);
        m_eTHETA_flt.push_back(parameters[Acts::eBoundTheta]);
        m_eQOP_flt.push_back(parameters[Acts::eBoundQOverP]);
        m_eT_flt.push_back(parameters[Acts::eBoundTime]);

        // filtered residual
        m_res_eLOC0_flt.push_back(parameters[Acts::eBoundLoc0] - truthLOC0);
        m_res_eLOC1_flt.push_back(parameters[Acts::eBoundLoc1] - truthLOC1);
        m_res_ePHI_flt.push_back(parameters[Acts::eBoundPhi] - truthPHI);
        m_res_eTHETA_flt.push_back(parameters[Acts::eBoundTheta] - truthTHETA);
        m_res_eQOP_flt.push_back(parameters[Acts::eBoundQOverP] - truthQOP);
        m_res_eT_flt.push_back(parameters[Acts::eBoundTime] - truthTIME);

        // filtered parameter error
        m_err_eLOC0_flt.push_back(
            sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
        m_err_eLOC1_flt.push_back(
            sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
        m_err_ePHI_flt.push_back(
            sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
        m_err_eTHETA_flt.push_back(
            sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
        m_err_eQOP_flt.push_back(
            sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
        m_err_eT_flt.push_back(
            sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

        // filtered parameter pull
        m_pull_eLOC0_flt.push_back(
            (parameters[Acts::eBoundLoc0] - truthLOC0) /
            sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
        m_pull_eLOC1_flt.push_back(
            (parameters[Acts::eBoundLoc1] - truthLOC1) /
            sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
        m_pull_ePHI_flt.push_back(
            (parameters[Acts::eBoundPhi] - truthPHI) /
            sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
        m_pull_eTHETA_flt.push_back(
            (parameters[Acts::eBoundTheta] - truthTHETA) /
            sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
        m_pull_eQOP_flt.push_back(
            (parameters[Acts::eBoundQOverP] - truthQOP) /
            sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
        m_pull_eT_flt.push_back(
            (parameters[Acts::eBoundTime] - truthTIME) /
            sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

        // more filtered parameter info
        const Acts::FreeVector freeParams =
            Acts::detail::transformBoundToFreeParameters(surface, gctx,
                                                         parameters);
        m_x_flt.push_back(freeParams[Acts::eFreePos0]);
        m_y_flt.push_back(freeParams[Acts::eFreePos1]);
        m_z_flt.push_back(freeParams[Acts::eFreePos2]);
        const auto p = std::abs(1 / freeParams[Acts::eFreeQOverP]);
        m_px_flt.push_back(p * freeParams[Acts::eFreeDir0]);
        m_py_flt.push_back(p * freeParams[Acts::eFreeDir1]);
        m_pz_flt.push_back(p * freeParams[Acts::eFreeDir2]);
        m_pT_flt.push_back(p * std::hypot(freeParams[Acts::eFreeDir0],
                                          freeParams[Acts::eFreeDir1]));
        m_eta_flt.push_back(
            Acts::VectorHelpers::eta(freeParams.segment<3>(Acts::eFreeDir0)));
        m_chi2.push_back(state.chi2());
        std::cout << "chi2 " << state.chi2() << std::endl;
      } else {
        // push default values if no filtered parameter
        m_eLOC0_flt.push_back(-99.);
        m_eLOC1_flt.push_back(-99.);
        m_ePHI_flt.push_back(-99.);
        m_eTHETA_flt.push_back(-99.);
        m_eQOP_flt.push_back(-99.);
        m_eT_flt.push_back(-99.);
        m_res_eLOC0_flt.push_back(-99.);
        m_res_eLOC1_flt.push_back(-99.);
        m_res_ePHI_flt.push_back(-99.);
        m_res_eTHETA_flt.push_back(-99.);
        m_res_eQOP_flt.push_back(-99.);
        m_res_eT_flt.push_back(-99.);
        m_err_eLOC0_flt.push_back(-99);
        m_err_eLOC1_flt.push_back(-99);
        m_err_ePHI_flt.push_back(-99);
        m_err_eTHETA_flt.push_back(-99);
        m_err_eQOP_flt.push_back(-99);
        m_err_eT_flt.push_back(-99);
        m_pull_eLOC0_flt.push_back(-99.);
        m_pull_eLOC1_flt.push_back(-99.);
        m_pull_ePHI_flt.push_back(-99.);
        m_pull_eTHETA_flt.push_back(-99.);
        m_pull_eQOP_flt.push_back(-99.);
        m_pull_eT_flt.push_back(-99.);
        m_x_flt.push_back(-99.);
        m_y_flt.push_back(-99.);
        m_z_flt.push_back(-99.);
        m_py_flt.push_back(-99.);
        m_pz_flt.push_back(-99.);
        m_pT_flt.push_back(-99.);
        m_eta_flt.push_back(-99.);
        m_chi2.push_back(-99.0);
      }

      // get the smoothed parameter
      bool smoothed = false;
      if (state.hasSmoothed()) {
        smoothed = true;
        m_nSmoothed++;
        auto parameters = state.smoothed();
        auto covariance = state.smoothedCovariance();

        // smoothed parameter
        m_eLOC0_smt.push_back(parameters[Acts::eBoundLoc0]);
        m_eLOC1_smt.push_back(parameters[Acts::eBoundLoc1]);
        m_ePHI_smt.push_back(parameters[Acts::eBoundPhi]);
        m_eTHETA_smt.push_back(parameters[Acts::eBoundTheta]);
        m_eQOP_smt.push_back(parameters[Acts::eBoundQOverP]);
        m_eT_smt.push_back(parameters[Acts::eBoundTime]);

        // smoothed residual
        m_res_eLOC0_smt.push_back(parameters[Acts::eBoundLoc0] - truthLOC0);
        m_res_eLOC1_smt.push_back(parameters[Acts::eBoundLoc1] - truthLOC1);
        m_res_ePHI_smt.push_back(parameters[Acts::eBoundPhi] - truthPHI);
        m_res_eTHETA_smt.push_back(parameters[Acts::eBoundTheta] - truthTHETA);
        m_res_eQOP_smt.push_back(parameters[Acts::eBoundQOverP] - truthQOP);
        m_res_eT_smt.push_back(parameters[Acts::eBoundTime] - truthTIME);

        // smoothed parameter error
        m_err_eLOC0_smt.push_back(
            sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
        m_err_eLOC1_smt.push_back(
            sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
        m_err_ePHI_smt.push_back(
            sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
        m_err_eTHETA_smt.push_back(
            sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
        m_err_eQOP_smt.push_back(
            sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
        m_err_eT_smt.push_back(
            sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

        // smoothed parameter pull
        m_pull_eLOC0_smt.push_back(
            (parameters[Acts::eBoundLoc0] - truthLOC0) /
            sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
        m_pull_eLOC1_smt.push_back(
            (parameters[Acts::eBoundLoc1] - truthLOC1) /
            sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
        m_pull_ePHI_smt.push_back(
            (parameters[Acts::eBoundPhi] - truthPHI) /
            sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
        m_pull_eTHETA_smt.push_back(
            (parameters[Acts::eBoundTheta] - truthTHETA) /
            sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
        m_pull_eQOP_smt.push_back(
            (parameters[Acts::eBoundQOverP] - truthQOP) /
            sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
        m_pull_eT_smt.push_back(
            (parameters[Acts::eBoundTime] - truthTIME) /
            sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

        // further smoothed parameter info
        const Acts::FreeVector freeParams =
            Acts::detail::transformBoundToFreeParameters(surface, gctx,
                                                         parameters);
        m_x_smt.push_back(freeParams[Acts::eFreePos0]);
        m_y_smt.push_back(freeParams[Acts::eFreePos1]);
        m_z_smt.push_back(freeParams[Acts::eFreePos2]);
        const auto p = std::abs(1 / freeParams[Acts::eFreeQOverP]);
        m_px_smt.push_back(p * freeParams[Acts::eFreeDir0]);
        m_py_smt.push_back(p * freeParams[Acts::eFreeDir1]);
        m_pz_smt.push_back(p * freeParams[Acts::eFreeDir2]);
        m_pT_smt.push_back(p * std::hypot(freeParams[Acts::eFreeDir0],
                                          freeParams[Acts::eFreeDir1]));
        m_eta_smt.push_back(
            Acts::VectorHelpers::eta(freeParams.segment<3>(Acts::eFreeDir0)));
      } else {
        // push default values if no smoothed parameter
        m_eLOC0_smt.push_back(-99.);
        m_eLOC1_smt.push_back(-99.);
        m_ePHI_smt.push_back(-99.);
        m_eTHETA_smt.push_back(-99.);
        m_eQOP_smt.push_back(-99.);
        m_eT_smt.push_back(-99.);
        m_res_eLOC0_smt.push_back(-99.);
        m_res_eLOC1_smt.push_back(-99.);
        m_res_ePHI_smt.push_back(-99.);
        m_res_eTHETA_smt.push_back(-99.);
        m_res_eQOP_smt.push_back(-99.);
        m_res_eT_smt.push_back(-99.);
        m_err_eLOC0_smt.push_back(-99);
        m_err_eLOC1_smt.push_back(-99);
        m_err_ePHI_smt.push_back(-99);
        m_err_eTHETA_smt.push_back(-99);
        m_err_eQOP_smt.push_back(-99);
        m_err_eT_smt.push_back(-99);
        m_pull_eLOC0_smt.push_back(-99.);
        m_pull_eLOC1_smt.push_back(-99.);
        m_pull_ePHI_smt.push_back(-99.);
        m_pull_eTHETA_smt.push_back(-99.);
        m_pull_eQOP_smt.push_back(-99.);
        m_pull_eT_smt.push_back(-99.);
        m_x_smt.push_back(-99.);
        m_y_smt.push_back(-99.);
        m_z_smt.push_back(-99.);
        m_px_smt.push_back(-99.);
        m_py_smt.push_back(-99.);
        m_pz_smt.push_back(-99.);
        m_pT_smt.push_back(-99.);
        m_eta_smt.push_back(-99.);
      }

      m_prt.push_back(predicted);
      m_flt.push_back(filtered);
      m_smt.push_back(smoothed);
      return true;
    });  // all states

    // fill the variables for one track to tree
    m_outputTree->Fill();
   
    std::cout << "filled"  << std::endl;
    
    // now reset
    m_t_x.clear();
    m_t_y.clear();
    m_t_z.clear();
    m_t_r.clear();
    m_t_dx.clear();
    m_t_dy.clear();
    m_t_dz.clear();
    m_t_eLOC0.clear();
    m_t_eLOC1.clear();
    m_t_ePHI.clear();
    m_t_eTHETA.clear();
    m_t_eQOP.clear();
    m_t_eT.clear();

    m_volumeID.clear();
    m_layerID.clear();
    m_moduleID.clear();
    m_lx_hit.clear();
    m_ly_hit.clear();
    m_x_hit.clear();
    m_y_hit.clear();
    m_z_hit.clear();
    m_res_x_hit.clear();
    m_res_y_hit.clear();
    m_err_x_hit.clear();
    m_err_y_hit.clear();
    m_pull_x_hit.clear();
    m_pull_y_hit.clear();
    m_dim_hit.clear();

    m_prt.clear();
    m_eLOC0_prt.clear();
    m_eLOC1_prt.clear();
    m_ePHI_prt.clear();
    m_eTHETA_prt.clear();
    m_eQOP_prt.clear();
    m_eT_prt.clear();
    m_res_eLOC0_prt.clear();
    m_res_eLOC1_prt.clear();
    m_res_ePHI_prt.clear();
    m_res_eTHETA_prt.clear();
    m_res_eQOP_prt.clear();
    m_res_eT_prt.clear();
    m_err_eLOC0_prt.clear();
    m_err_eLOC1_prt.clear();
    m_err_ePHI_prt.clear();
    m_err_eTHETA_prt.clear();
    m_err_eQOP_prt.clear();
    m_err_eT_prt.clear();
    m_pull_eLOC0_prt.clear();
    m_pull_eLOC1_prt.clear();
    m_pull_ePHI_prt.clear();
    m_pull_eTHETA_prt.clear();
    m_pull_eQOP_prt.clear();
    m_pull_eT_prt.clear();
    m_x_prt.clear();
    m_y_prt.clear();
    m_z_prt.clear();
    m_px_prt.clear();
    m_py_prt.clear();
    m_pz_prt.clear();
    m_eta_prt.clear();
    m_pT_prt.clear();

    m_flt.clear();
    m_eLOC0_flt.clear();
    m_eLOC1_flt.clear();
    m_ePHI_flt.clear();
    m_eTHETA_flt.clear();
    m_eQOP_flt.clear();
    m_eT_flt.clear();
    m_res_eLOC0_flt.clear();
    m_res_eLOC1_flt.clear();
    m_res_ePHI_flt.clear();
    m_res_eTHETA_flt.clear();
    m_res_eQOP_flt.clear();
    m_res_eT_flt.clear();
    m_err_eLOC0_flt.clear();
    m_err_eLOC1_flt.clear();
    m_err_ePHI_flt.clear();
    m_err_eTHETA_flt.clear();
    m_err_eQOP_flt.clear();
    m_err_eT_flt.clear();
    m_pull_eLOC0_flt.clear();
    m_pull_eLOC1_flt.clear();
    m_pull_ePHI_flt.clear();
    m_pull_eTHETA_flt.clear();
    m_pull_eQOP_flt.clear();
    m_pull_eT_flt.clear();
    m_x_flt.clear();
    m_y_flt.clear();
    m_z_flt.clear();
    m_px_flt.clear();
    m_py_flt.clear();
    m_pz_flt.clear();
    m_eta_flt.clear();
    m_pT_flt.clear();
    m_chi2.clear();

    m_smt.clear();
    m_eLOC0_smt.clear();
    m_eLOC1_smt.clear();
    m_ePHI_smt.clear();
    m_eTHETA_smt.clear();
    m_eQOP_smt.clear();
    m_eT_smt.clear();
    m_res_eLOC0_smt.clear();
    m_res_eLOC1_smt.clear();
    m_res_ePHI_smt.clear();
    m_res_eTHETA_smt.clear();
    m_res_eQOP_smt.clear();
    m_res_eT_smt.clear();
    m_err_eLOC0_smt.clear();
    m_err_eLOC1_smt.clear();
    m_err_ePHI_smt.clear();
    m_err_eTHETA_smt.clear();
    m_err_eQOP_smt.clear();
    m_err_eT_smt.clear();
    m_pull_eLOC0_smt.clear();
    m_pull_eLOC1_smt.clear();
    m_pull_ePHI_smt.clear();
    m_pull_eTHETA_smt.clear();
    m_pull_eQOP_smt.clear();
    m_pull_eT_smt.clear();
    m_x_smt.clear();
    m_y_smt.clear();
    m_z_smt.clear();
    m_px_smt.clear();
    m_py_smt.clear();
    m_pz_smt.clear();
    m_eta_smt.clear();
    m_pT_smt.clear();
  }  // all trajectories


    std::cout << "else counter " << counter_else << std::endl;
    std::cout << "if counter " << counter_if << std::endl;

  return ProcessCode::SUCCESS;
}
