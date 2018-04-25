// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Surfaces/InfiniteBounds.hpp"
#include "ACTS/Utilities/VariantData.hpp"

Acts::variant_data
Acts::InfiniteBounds::toVariantData() const
{
  using namespace std::string_literals;
  return variant_map({{"type", "InfiniteBounds"s}});
}