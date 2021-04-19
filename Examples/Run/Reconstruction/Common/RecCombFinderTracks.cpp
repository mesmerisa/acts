// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Digitization/DigitizationOptions.hpp"
#include "ActsExamples/Digitization/HitSmearing.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFitterPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryParametersWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriterBGV.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingOptions.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/TruthTracking/CombTrackFinderBGV.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Definitions/Units.hpp>
#include "ActsExamples/Io/Root/RootVertexAndTrackWriterBGV.hpp"
#include "ActsExamples/Vertexing/VertexFitterAlgorithmFromTraj.hpp"
#include "ActsExamples/TruthTracking/TruthVertexFinder.hpp"

#include <memory>

#include "RecInput.hpp"

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int runRecCombFinderTracks(int argc, char* argv[],
                      std::shared_ptr<ActsExamples::IBaseDetector> detector) {
  // setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  detector->addOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addFittingOptions(desc);
  Options::addDigitizationOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel = Options::readLogLevel(vm);
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  auto rnd = std::make_shared<const ActsExamples::RandomNumbers>(
      Options::readRandomNumbersConfig(vm));
  auto dirNav = vm["directed-navigation"].as<bool>();

  // Setup detector geometry
  auto geometry = Geometry::build(vm, *detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (auto cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }
  // Setup the magnetic field
  auto magneticField = Options::readMagneticField(vm);

  // Read the sim hits
  auto simHitReaderCfg = setupSimHitReading(vm, sequencer);
  // Read the particles
  auto particleReader = setupParticleReading(vm, sequencer);

  // Run the sim hits smearing
  auto hitSmearingCfg = setupDigitization(vm, sequencer, rnd, trackingGeometry,
                                   simHitReaderCfg.outputSimHits);
                                   
                                
  // Run the sim hits smearing
 // auto hitSmearingCfg = setupSimHitSmearing(
  //    vm, sequencer, rnd, trackingGeometry, simHitReaderCfg.outputSimHits);
      
  // Run the particle selection
  // The pre-selection will select truth particles satisfying provided criteria
  // from all particles read in by particle reader for further processing. It
  // has no impact on the truth hits read-in by the cluster reader.
  TruthSeedSelector::Config particleSelectorCfg;
  particleSelectorCfg.inputParticles = particleReader.outputParticles;
  particleSelectorCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  particleSelectorCfg.inputMeasurements = hitSmearingCfg.outputMeasurements;
  particleSelectorCfg.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  particleSelectorCfg.outputParticles = "particles_selected";
  particleSelectorCfg.nHitsMin = 3;
  
  // 1 ring detector:
  //particleSelectorCfg.etaMin = 2.814708661281755;
  //particleSelectorCfg.etaMax = 3.453280970676263;
  
  // 2 ring detector:  
  particleSelectorCfg.etaMin = 2.428760570711203;
  particleSelectorCfg.etaMax = 3.453280970676263;  
  
  sequencer.addAlgorithm(
      std::make_shared<TruthSeedSelector>(particleSelectorCfg, logLevel));

  // The selected particles
  const auto& inputParticles = particleSelectorCfg.outputParticles;

  // The fitter needs the measurements (proto tracks) and initial
  // track states (proto states). The elements in both collections
  // must match and must be created from the same input particles.
  // Create truth tracks
 /* TruthTrackFinder::Config trackFinderCfg;
  trackFinderCfg.inputParticles = inputParticles;
  trackFinderCfg.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  trackFinderCfg.outputProtoTracks = "prototracks";
  sequencer.addAlgorithm(
      std::make_shared<TruthTrackFinder>(trackFinderCfg, logLevel));*/



  CombTrackFinderBGV::Config trackFinderCfg;
  trackFinderCfg.inputParticles = inputParticles; 
  trackFinderCfg.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  trackFinderCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  trackFinderCfg.outputProtoTracks = "prototracks";
  sequencer.addAlgorithm(
      std::make_shared<CombTrackFinderBGV>(trackFinderCfg, logLevel));

   // Run the particle smearing
  auto particleSmearingCfg =
      setupParticleSmearingCombFinder(vm, sequencer, rnd, inputParticles, trackFinderCfg.outputProtoTracks, simHitReaderCfg.outputSimHits);
  
 /*auto particleSmearingCfg =
      setupParticleSmearing(vm, sequencer, rnd, inputParticles);*/
  

  SurfaceSortingAlgorithm::Config sorterCfg;
  // Setup the surface sorter if running direct navigator
  sorterCfg.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  sorterCfg.inputSimulatedHits = simHitReaderCfg.outputSimHits;
  sorterCfg.inputMeasurementSimHitsMap = hitSmearingCfg.outputMeasurementSimHitsMap;
  sorterCfg.outputProtoTracks = "sortedprototracks";
  if (dirNav) {
    sequencer.addAlgorithm(
        std::make_shared<SurfaceSortingAlgorithm>(sorterCfg, logLevel));
  }

  // setup the fitter
  TrackFittingAlgorithm::Config fitter;
  fitter.inputMeasurements = hitSmearingCfg.outputMeasurements;
  fitter.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  fitter.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  if (dirNav) {
    fitter.inputProtoTracks = sorterCfg.outputProtoTracks;
  }
  fitter.inputInitialTrackParameters =
      particleSmearingCfg.outputTrackParameters;
  fitter.outputTrajectories = "trajectories";
  fitter.directNavigation = dirNav;
  fitter.trackingGeometry = trackingGeometry;
  fitter.dFit = TrackFittingAlgorithm::makeTrackFitterFunction(magneticField);
  fitter.fit = TrackFittingAlgorithm::makeTrackFitterFunction(trackingGeometry,
                                                              magneticField);
  sequencer.addAlgorithm(
      std::make_shared<TrackFittingAlgorithm>(fitter, logLevel));

  // write track states from fitting
  RootTrajectoryStatesWriterBGV::Config trackStatesWriter;
  trackStatesWriter.inputTrajectories = fitter.outputTrajectories;
  trackStatesWriter.inputParticles = inputParticles;
  trackStatesWriter.inputSimHits = simHitReaderCfg.outputSimHits;
  trackStatesWriter.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  trackStatesWriter.inputMeasurementSimHitsMap =
      hitSmearingCfg.outputMeasurementSimHitsMap;
  trackStatesWriter.outputDir = outputDir;
  trackStatesWriter.outputFilename = "trackstates_fitter.root";
  trackStatesWriter.outputTreename = "trackstates_fitter";
  sequencer.addWriter(std::make_shared<RootTrajectoryStatesWriterBGV>(
      trackStatesWriter, logLevel));

 /* // write track parameters from fitting
  RootTrajectoryParametersWriter::Config trackParamsWriter;
  trackParamsWriter.inputTrajectories = fitter.outputTrajectories;
  trackParamsWriter.inputParticles = inputParticles;
  trackParamsWriter.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  trackParamsWriter.outputDir = outputDir;
  trackParamsWriter.outputFilename = "trackparams_fitter.root";
  trackParamsWriter.outputTreename = "trackparams_fitter";
  sequencer.addWriter(std::make_shared<RootTrajectoryParametersWriter>(
      trackParamsWriter, logLevel));

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
      
  //////////////////////////////////////////////////////////////////////////////////////
  // find true primary vertices w/o secondary particles
  TruthVertexFinder::Config findVertices;
  findVertices.inputParticles = inputParticles; //selectParticles.outputParticles;
  findVertices.outputProtoVertices = "protovertices";
  findVertices.excludeSecondaries = true;
  sequencer.addAlgorithm(
      std::make_shared<TruthVertexFinder>(findVertices, logLevel));
      
  // fit vertices using the Billoir fitter
  VertexFitterAlgorithmFromTraj::Config fitVertices(magneticField);
  fitVertices.inputTrajectories = fitter.outputTrajectories; // smearParticles.outputTrackParameters;
  fitVertices.inputProtoVertices = findVertices.outputProtoVertices;
  fitVertices.outputFittedVertices = "fitted_vertices";
  fitVertices.doConstrainedFit = false;
  //fitVertices.bField = Acts::Vector3(0_T, 0_T, 0_T);
  sequencer.addAlgorithm(
      std::make_shared<VertexFitterAlgorithmFromTraj>(fitVertices, logLevel));
      
  RootVertexAndTrackWriterBGV::Config writerCfg;
  writerCfg.collection = fitVertices.outputFittedVertices;
  writerCfg.filePath = joinPaths(outputDir, fitVertices.outputFittedVertices + ".root");
  sequencer.addWriter(
      std::make_shared<RootVertexAndTrackWriterBGV>(writerCfg, logLevel));      
     
   //////////////////////////////////////////////////////////////////////////////////////    
*/
  return sequencer.run();
}
