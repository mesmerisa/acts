// This file is part of the ACTS project.

#ifndef DIGITIZATION_DIGITIZATION_HPP
#define DIGITIZATION_DIGITIZATION_HPP

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "ACTS/Digitization/DigitizationCell.hpp"

namespace Acts {
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
    Graph;

/// @brief merge cells
/// This function recieves digitization cells and merges the cells which are at
/// the same position. In case we have analgue readout the energy is summed up.
/// Furthermore an energycut can be applied. It uses boost
/// connected_components
/// (https://www.boost.org/doc/libs/1_46_1/libs/graph/doc/connected_components.html)
/// @param cells all digitization cells
/// @param anaglogueReadout flag indicating if the module has analgue
/// readout
/// @param energyCut Possible energy cut to be applied
/// @return the merged digitization cells
const std::vector<Acts::DigitizationCell>
mergeCells(std::vector<Acts::DigitizationCell>& cells,
           bool                                 analogueReadout = false,
           double                               energyCut       = 0.);

/// @brief create clusters
/// This function recieves digitization cells and bundles the neighbouring
/// cells. It uses boost connected_components
/// (https://www.boost.org/doc/libs/1_46_1/libs/graph/doc/connected_components.html)
/// @param cells all digitization cells
/// @param commonCorner flag indicating if also cells sharing a common corner
/// should be merged (all cells sharing a common edge are merged per default)
/// @return vector (the different clusters) of vector of digitization cells (the
/// cells which belong to each cluster)
const std::vector<std::vector<Acts::DigitizationCell>>
createClusters(const std::vector<Acts::DigitizationCell>& cells,
               bool                                       commonCorner = false);
}

#endif  // DIGITIZATION_DIGITIZATION_HPP
