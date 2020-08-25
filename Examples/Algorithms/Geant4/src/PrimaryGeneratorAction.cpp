// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "PrimaryGeneratorAction.hpp"

#include <stdexcept>

#include <G4Event.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4RandomDirection.hh>
#include <G4UnitsTable.hh>
#include <Randomize.hh>

#include <CLHEP/Random/RandFlat.h>
#include "G4Exp.hh"

using namespace ActsExamples;

PrimaryGeneratorAction* PrimaryGeneratorAction::s_instance = nullptr;

PrimaryGeneratorAction* PrimaryGeneratorAction::instance() {
  return s_instance;
}

PrimaryGeneratorAction::PrimaryGeneratorAction(const G4String& particleName,
                                               G4double energy,
                                               G4int randomSeed1,
                                               G4int randomSeed2,
                                               std::array<double, 2> etaRange)
    : G4VUserPrimaryGeneratorAction(),
      m_particleGun(std::make_unique<G4ParticleGun>(1)), m_eta(etaRange) {
  if (s_instance) {
    throw std::logic_error(
        "Attempted to duplicate the PrimaryGeneratorAction singleton");
  } else {
    s_instance = this;
  }

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName);
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleEnergy(energy);
  G4UnitDefinition::PrintUnitsTable();

  // set the random seeds
  CLHEP::HepRandom::getTheEngine()->setSeed(randomSeed1, randomSeed2);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  s_instance = nullptr;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  // this function is called at the begining of event
 // G4double phi = -M_PI + G4UniformRand() * 2. * M_PI;
  //G4double theta = G4UniformRand() * M_PI;   
  //G4double eta = CLHEP::RandFlat::shoot(2,4);
  
  G4double phi = CLHEP::RandFlat::shoot(0.0, 2. * M_PI);

  G4double theta_1 = 2.0 * atan(exp(-m_eta[0])); 
  G4double theta_2 = 2.0 * atan(exp(-m_eta[1]));
  G4double theta = 0.0;
  
  if (theta_1 < theta_2) theta = CLHEP::RandFlat::shoot(theta_1, theta_2);
  else theta = CLHEP::RandFlat::shoot(theta_2, theta_1);

  //std::cout << m_eta[0] << " " << m_eta[1] << " " << theta_1 << " " << theta_2 << " " << theta << std::endl;
  //std::cout << G4UniformRand() * M_PI << std::endl;
  
  // build a direction
  m_direction =
      G4ThreeVector(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
  m_position = G4ThreeVector(0., 0., 0.);
  // set to the particle gun and
  m_particleGun->SetParticleMomentumDirection(m_direction);
  m_particleGun->SetParticlePosition(m_position);
  m_particleGun->GeneratePrimaryVertex(anEvent);
}
