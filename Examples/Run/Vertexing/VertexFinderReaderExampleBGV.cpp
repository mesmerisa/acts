// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file Find vertices using truth particle information as input
///
/// Reads truth particles from TrackMl files and use the truth information
/// to generate smeared track parameters. Use this pseudo-reconstructed
/// tracks as the input to the vertex finder.

#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Printers/TrackParametersPrinter.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/Vertexing/IterativeVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/AdaptiveMultiVertexFinderAlgorithm.hpp"
#include "ActsExamples/Io/Root/RootVertexAndTrackWriterBGV.hpp"
#include "ActsExamples/Vertexing/VertexFitterAlgorithmBGV.hpp"
#include "ActsExamples/TruthTracking/TruthVertexFinder.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  ParticleSelector::addOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc);
  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  // basic setup
  auto logLevel = Options::readLogLevel(vars);
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));
  Sequencer sequencer(Options::readSequencerConfig(vars));

  // setup particle reader generator
  CsvParticleReader::Config readParticles =
      Options::readCsvParticleReaderConfig(vars);
  readParticles.inputStem = "particles_initial";
  readParticles.outputParticles = "particles";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(readParticles, logLevel));

  // pre-select particles
  ParticleSelector::Config selectParticles = ParticleSelector::readConfig(vars);
  selectParticles.inputParticles = readParticles.outputParticles;
  selectParticles.outputParticles = "particles_selected";
  selectParticles.etaMin = 1.9;
  selectParticles.etaMax = 4.5;
  // smearing only works with charge particles for now
  selectParticles.removeNeutral = true;
  sequencer.addAlgorithm(
      std::make_shared<ParticleSelector>(selectParticles, logLevel));

  // simulate track reconstruction by smearing truth track parameters
  ParticleSmearing::Config smearParticles;
  smearParticles.inputParticles = selectParticles.outputParticles;
  smearParticles.outputTrackParameters = "trackparameters";
  smearParticles.sigmaD0 = 20_um;
  smearParticles.sigmaD0PtA = 30_um;
  smearParticles.sigmaD0PtB = 0.3 / 1_GeV;
  smearParticles.sigmaZ0 = 20_um;
  smearParticles.sigmaZ0PtA = 30_um;
  smearParticles.sigmaZ0PtB = 0.3 / 1_GeV;
  smearParticles.sigmaPhi = 0.001_degree;
  smearParticles.sigmaTheta = 0.0001_degree;
  smearParticles.sigmaPRel = 0.001;
  smearParticles.sigmaT0 = 1_ns;
  smearParticles.randomNumbers = rnd;
  sequencer.addAlgorithm(
      std::make_shared<ParticleSmearing>(smearParticles, logLevel));

  // print input track parameters
  /*TrackParametersPrinter::Config printTracks;
  printTracks.inputTrackParameters = smearParticles.outputTrackParameters;
  sequencer.addAlgorithm(
      std::make_shared<TrackParametersPrinter>(printTracks, logLevel));*/

  // find vertices
  // proto vertices contains indices of tracks that belong to a vertex
  /*IterativeVertexFinderAlgorithm::Config findVertices;
  findVertices.inputTrackParameters = smearParticles.outputTrackParameters;
  findVertices.outputProtoVertices = "protovertices";
  findVertices.bField = Acts::Vector3D(0_T, 0_T, 0.00000001_T);
  sequencer.addAlgorithm(
      std::make_shared<IterativeVertexFinderAlgorithm>(findVertices, logLevel));*/

  // find vertices
  /*AdaptiveMultiVertexFinderAlgorithm::Config findVertices;
  findVertices.inputTrackParameters = smearParticles.outputTrackParameters;
  findVertices.outputProtoVertices = "protovertices";
  findVertices.bField = Acts::Vector3D(0_T, 0_T, 0_T);
  sequencer.addAlgorithm(std::make_shared<AdaptiveMultiVertexFinderAlgorithm>(
      findVertices, logLevel));*/

  
  // find true primary vertices w/o secondary particles
  TruthVertexFinder::Config findVertices;
  findVertices.inputParticles = selectParticles.outputParticles;
  findVertices.outputProtoVertices = "protovertices";
  findVertices.excludeSecondaries = true;
  sequencer.addAlgorithm(
      std::make_shared<TruthVertexFinder>(findVertices, logLevel));

  // fit vertices using the Billoir fitter
  VertexFitterAlgorithmBGV::Config fitVertices;
  fitVertices.inputTrackParameters = smearParticles.outputTrackParameters;
  fitVertices.inputProtoVertices = findVertices.outputProtoVertices;
  fitVertices.outputFittedVertices = "fitted_vertices";
  fitVertices.bField = Acts::Vector3D(0_T, 0_T, 0_T);
  sequencer.addAlgorithm(
      std::make_shared<VertexFitterAlgorithmBGV>(fitVertices, logLevel));
  







  ////////////////////////////////////////////////


  auto outputDir = ensureWritableDirectory(vars["output-dir"].as<std::string>());

  RootVertexAndTrackWriterBGV::Config writerCfg;
  writerCfg.collection = fitVertices.outputFittedVertices;
  writerCfg.filePath = joinPaths(outputDir, findVertices.outputProtoVertices + ".root");
  sequencer.addWriter(
      std::make_shared<RootVertexAndTrackWriterBGV>(writerCfg, logLevel));

  return sequencer.run();
}
