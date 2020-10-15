// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/AccumulatedVolumeMaterial.hpp"

#include "Acts/Material/detail/AverageMaterials.hpp"

void Acts::AccumulatedVolumeMaterial::accumulate(const MaterialSlab& mat) {
  m_average = detail::combineSlabs(m_average, mat);
}
