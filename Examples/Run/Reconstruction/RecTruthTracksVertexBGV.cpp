// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/HitSmearing.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFitterPerformanceWriter.hpp"
//#include "ActsExamples/Io/Root/RootTrajectoryWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryParametersWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Definitions/Units.hpp>
#include "ActsExamples/Io/Root/RootVertexAndTrackWriterBGV.hpp"
#include "ActsExamples/Vertexing/VertexFitterAlgorithmFromTrajBGV.hpp"
#include "ActsExamples/TruthTracking/TruthVertexFinder.hpp"

#include "ActsExamples/Vertexing/AdaptiveMultiVertexFinderAlgorithmFromTrajBGV.hpp"
#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  TGeoDetector detector;

 
  // setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc);
  detector.addOptions(desc);
  Options::addBFieldOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel = Options::readLogLevel(vm);
  auto inputDir = vm["input-dir"].as<std::string>();
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  auto rnd = std::make_shared<ActsExamples::RandomNumbers>(
      Options::readRandomNumbersConfig(vm));

  // Setup detector geometry
  auto geometry = Geometry::build(vm, detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (auto cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }
  // Setup the magnetic field
  auto magneticField = Options::readBField(vm);

  // Read particles (initial states) and clusters from CSV files
  auto particleReader = Options::readCsvParticleReaderConfig(vm);
  particleReader.inputStem = "particles_initial";
  particleReader.outputParticles = "particles_initial";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(particleReader, logLevel));
  // Read clusters from CSV files
  auto clusterReaderCfg = Options::readCsvPlanarClusterReaderConfig(vm);
  clusterReaderCfg.trackingGeometry = trackingGeometry;
  clusterReaderCfg.outputClusters = "clusters";
  clusterReaderCfg.outputHitIds = "hit_ids";
  // only simhits are used and the map will be re-created by the digitizer
  clusterReaderCfg.outputMeasurementParticlesMap = "unused-hit_particles_map";
  clusterReaderCfg.outputSimHits = "hits";
  sequencer.addReader(
      std::make_shared<CsvPlanarClusterReader>(clusterReaderCfg, logLevel));

  // Create smeared measurements
  HitSmearing::Config hitSmearingCfg;
  hitSmearingCfg.inputSimHits = clusterReaderCfg.outputSimHits;
  hitSmearingCfg.outputMeasurements = "measurements";
  hitSmearingCfg.outputSourceLinks = "sourcelinks";
  hitSmearingCfg.outputMeasurementParticlesMap = "measurement_particles_map";
  hitSmearingCfg.outputMeasurementSimHitsMap = "measurement_simhits_map";
  hitSmearingCfg.sigmaLoc0 = 50_um;
  hitSmearingCfg.sigmaLoc1 = 50_um;
  hitSmearingCfg.randomNumbers = rnd;
  hitSmearingCfg.trackingGeometry = trackingGeometry;
  sequencer.addAlgorithm(
      std::make_shared<HitSmearing>(hitSmearingCfg, logLevel));

  // Pre-select particles
  // The pre-selection will select truth particles satisfying provided criteria
  // from all particles read in by particle reader for further processing. It
  // has no impact on the truth hits read-in by the cluster reader.
  // @TODO: add options for truth particle selection criteria
  TruthSeedSelector::Config particleSelectorCfg;
  particleSelectorCfg.inputParticles = particleReader.outputParticles;
  particleSelectorCfg.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  particleSelectorCfg.outputParticles = "particles_selected";
  particleSelectorCfg.nHitsMin = 3;
  particleSelectorCfg.etaMin = 1.9;
  particleSelectorCfg.etaMax = 4.4;
  sequencer.addAlgorithm(
      std::make_shared<TruthSeedSelector>(particleSelectorCfg, logLevel));

  // The fitter needs the measurements (proto tracks) and initial
  // track states (proto states). The elements in both collections
  // must match and must be created from the same input particles.
  const auto& inputParticles = particleSelectorCfg.outputParticles;
  // Create truth tracks
  TruthTrackFinder::Config trackFinderCfg;
  trackFinderCfg.inputParticles = inputParticles;
  trackFinderCfg.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  trackFinderCfg.outputProtoTracks = "prototracks";
  sequencer.addAlgorithm(
      std::make_shared<TruthTrackFinder>(trackFinderCfg, logLevel));
  // Create smeared particles states
  ParticleSmearing::Config particleSmearingCfg;
  particleSmearingCfg.inputParticles = inputParticles;
  particleSmearingCfg.outputTrackParameters = "smearedparameters";
  particleSmearingCfg.randomNumbers = rnd;
  // Gaussian sigmas to smear particle parameters
  particleSmearingCfg.sigmaD0 = 20_um;
  particleSmearingCfg.sigmaD0PtA = 30_um;
  particleSmearingCfg.sigmaD0PtB = 0.3 / 1_GeV;
  particleSmearingCfg.sigmaZ0 = 20_um;
  particleSmearingCfg.sigmaZ0PtA = 30_um;
  particleSmearingCfg.sigmaZ0PtB = 0.3 / 1_GeV;
  particleSmearingCfg.sigmaPhi = 0.01_degree;
  particleSmearingCfg.sigmaTheta = 0.001_degree;
  particleSmearingCfg.sigmaPRel = 0.01;
  particleSmearingCfg.sigmaT0 = 1_ns;
  sequencer.addAlgorithm(
      std::make_shared<ParticleSmearing>(particleSmearingCfg, logLevel));

  // setup the fitter
  TrackFittingAlgorithm::Config fitter;
  fitter.inputMeasurements = hitSmearingCfg.outputMeasurements;
  fitter.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  fitter.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  fitter.inputInitialTrackParameters =
      particleSmearingCfg.outputTrackParameters;
  fitter.outputTrajectories = "trajectories";
  fitter.fit = TrackFittingAlgorithm::makeTrackFitterFunction(trackingGeometry,
                                                              magneticField);
  sequencer.addAlgorithm(
      std::make_shared<TrackFittingAlgorithm>(fitter, logLevel));
      
      
      
      
      
  //////////////////////////////////////////////////////////////////////////////////////
  // find true primary vertices w/o secondary particles
  TruthVertexFinder::Config findVertices;
  findVertices.inputParticles = inputParticles; //selectParticles.outputParticles;
  findVertices.outputProtoVertices = "protovertices";
  findVertices.excludeSecondaries = true;
  sequencer.addAlgorithm(
      std::make_shared<TruthVertexFinder>(findVertices, logLevel));
      
  // fit vertices using the Billoir fitter
  VertexFitterAlgorithmFromTrajBGV::Config fitVertices;
  fitVertices.inputTrajectories = fitter.outputTrajectories; // smearParticles.outputTrackParameters;
  fitVertices.inputProtoVertices = findVertices.outputProtoVertices;
  fitVertices.outputFittedVertices = "fitted_vertices";
  fitVertices.doConstrainedFit = false;
  fitVertices.bField = Acts::Vector3D(0_T, 0_T, 0_T);
  sequencer.addAlgorithm(
      std::make_shared<VertexFitterAlgorithmFromTrajBGV>(fitVertices, logLevel));
      
  RootVertexAndTrackWriterBGV::Config writerCfg;
  writerCfg.collection = fitVertices.outputFittedVertices;
  writerCfg.filePath = joinPaths(outputDir, fitVertices.outputFittedVertices + ".root");
  sequencer.addWriter(
      std::make_shared<RootVertexAndTrackWriterBGV>(writerCfg, logLevel));        
  
      
     
   // find vertices
  /*AdaptiveMultiVertexFinderAlgorithmFromTrajBGV::Config findVertices;
  findVertices.inputTrajectories = fitter.outputTrajectories; // XXX smearParticles.outputTrackParameters;
  findVertices.outputProtoVertices = "protovertices";
  findVertices.outputFoundVertices = "found_vertices";
  findVertices.bField = Acts::Vector3D(0_T, 0_T, 0_T);
  sequencer.addAlgorithm(std::make_shared<AdaptiveMultiVertexFinderAlgorithmFromTrajBGV>(
      findVertices, logLevel));    
      
  RootVertexAndTrackWriterBGV::Config writerCfg;
  writerCfg.collection = findVertices.outputFoundVertices; //fitVertices.outputFittedVertices;
  writerCfg.filePath = joinPaths(outputDir,  findVertices.outputFoundVertices + ".root");
  sequencer.addWriter(
      std::make_shared<RootVertexAndTrackWriterBGV>(writerCfg, logLevel));   */   
      
  

  //////////////////////////////////////////////////////////////////////////////////////    
      
      
      
      
      
      
      

 /* // write tracks from fitting
  RootTrajectoryWriter::Config trackWriter;
  trackWriter.inputTrajectories = fitter.outputTrajectories;
  trackWriter.inputParticles = inputParticles;
  trackWriter.inputSimHits = clusterReaderCfg.outputSimHits;
  trackWriter.inputMeasurements = hitSmearingCfg.outputMeasurements;
  trackWriter.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  trackWriter.inputMeasurementSimHitsMap =
      hitSmearingCfg.outputMeasurementSimHitsMap;
  trackWriter.outputDir = outputDir;
  trackWriter.outputFilename = "tracks.root";
  trackWriter.outputTreename = "tracks";
  sequencer.addWriter(
      std::make_shared<RootTrajectoryWriter>(trackWriter, logLevel));*/

  // write reconstruction performance data
  TrackFinderPerformanceWriter::Config perfFinder;
  perfFinder.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  perfFinder.inputParticles = inputParticles;
  perfFinder.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  perfFinder.outputDir = outputDir;
  sequencer.addWriter(
      std::make_shared<TrackFinderPerformanceWriter>(perfFinder, logLevel));
  TrackFitterPerformanceWriter::Config perfFitter;
  perfFitter.inputTrajectories = fitter.outputTrajectories;
  perfFitter.inputParticles = inputParticles;
  perfFitter.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  perfFitter.outputDir = outputDir;
  sequencer.addWriter(
      std::make_shared<TrackFitterPerformanceWriter>(perfFitter, logLevel));

  return sequencer.run();
}
