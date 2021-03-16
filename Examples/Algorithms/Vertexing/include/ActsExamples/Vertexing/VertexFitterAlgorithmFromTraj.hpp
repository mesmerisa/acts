
// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"

#include <string>

namespace ActsExamples {

class VertexFitterAlgorithmFromTraj final : public BareAlgorithm {
 public:
  struct Config {
    Config(std::shared_ptr<Acts::MagneticFieldProvider> magneticField)
        : bField(magneticField) {}
    /// Input track parameters collection
    std::string inputTrackParameters;
    /// Input trajectory collection.
    std::string inputTrajectories;
     /// Output fitted vertex collection.
    std::string outputFittedVertices;
    /// Input proto vertex collection
    std::string inputProtoVertices;
    /// The magnetic field
    std::shared_ptr<Acts::MagneticFieldProvider> bField;
    /// Constraint vertex fit bool
    bool doConstrainedFit = false;
    /// Vertex constraint position
    Acts::Vector3 constraintPos = Acts::Vector3(0, 0, 0);
    /// Vertex constraint covariance matrix
    Acts::SymMatrix3 constraintCov =
        Acts::Vector3(3 * Acts::UnitConstants::mm, 3 * Acts::UnitConstants::mm,
                      10 * Acts::UnitConstants::mm)
            .asDiagonal();
  };

  VertexFitterAlgorithmFromTraj(const Config& cfg, Acts::Logging::Level lvl);

  /// Fit the input vertices.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
