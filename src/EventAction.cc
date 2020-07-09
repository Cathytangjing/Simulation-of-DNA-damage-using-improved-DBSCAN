//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class#include <cstdio>
#include <cstdio>
#include <string>
#include <sstream>
#include <chrono>

#include "EventAction.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4Run.hh"

//#include "Analysis.hh"
#include "ClusteringAlgo.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "SizeDistribution.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction() : G4UserEventAction() {
  //default parameter values
  event_energy_ = 0.0;
  event_step_length_ = 0.0;

  // Create clustering algorithm
  // These default values have been tuned for the Physics List G4EmDNAPhysics
  // to reproduce data published by:
  // Francis et al. 2011 Comput. Meth. Programs. Biomed. 2011 101(3)
  clustering_ = new ClusteringAlgo(3.2 * nanometer, 2, 0.2, 0.2, 5 * eV, 37.5 * eV, 0.4);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() {
  delete clustering_;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *) {
  event_energy_ = 0.;
  clustering_->Purge();

  G4cout << "--------Mass: " << (G4LogicalVolumeStore::GetInstance()->
      GetVolume("Target")->GetMass() / kg) << "--------" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *event) {
  std::map<G4int, G4int> sizeDistribution = clustering_->RunClustering();


  G4int ssb_num = clustering_->GetSSB();
  G4int cssbn_num = clustering_->GetCSSBN();
  G4int dsbn_num = clustering_->GetDSBN();
  G4int dsbp_num = clustering_->GetDSBP();
  G4int dsbpp_num = clustering_->GetDSBPP();

  RunAction *run_action =
      (RunAction *) G4RunManager::GetRunManager()->GetUserRunAction();
  run_action->AddEnergyDeposition(event_energy_);
  run_action->AddSingleStrandBreak(ssb_num);
  run_action->AddComplexSingleStrandBreakNew(cssbn_num);
  run_action->AddDoubleStrandBreakNew(dsbn_num);
  run_action->AddDoubleStrandBreakPlus(dsbp_num);
  run_action->AddDoubleStrandBreakPlusPlus(dsbpp_num);

  for (auto it = sizeDistribution.begin(); it != sizeDistribution.end(); it++) {
    for (int i = 0; i < it->second; i++) {
      run_action->AddDistribution(it->first);
    }
  }

  G4double total_energy_deposit = event_energy_ / keV;
  G4double absorbed_dose =
      (event_energy_ / joule) / (G4LogicalVolumeStore::GetInstance()->GetVolume("Target")->GetMass() / kg);

  G4double ssb_yield = static_cast<double>(ssb_num) / absorbed_dose / (6);
  G4double cssbn_yield = static_cast<double>(cssbn_num) / absorbed_dose / (6);
  G4double dsbn_yield = static_cast<double>(dsbn_num) / absorbed_dose / (6);
  G4double dsbp_yield = static_cast<double>(dsbp_num) / absorbed_dose / (6);
  G4double dsbpp_yield = static_cast<double>(dsbpp_num) / absorbed_dose / (6);

  G4cout << "-----------------RESULT-----------------" << G4endl;
  G4cout << "| ssb_yield |" << ssb_yield << " Gbp-1Gy-1" << "| ssb_num |" << ssb_num << G4endl;
  G4cout << "| cssbn_yield |" << cssbn_yield << " Gbp-1Gy-1" << "| cssbn_num |" << cssbn_num << G4endl;
  G4cout << "| dsbn_yield |" << dsbn_yield << " Gbp-1Gy-1" << "| dsbn_num |" << dsbn_num << G4endl;
  G4cout << "| dsbp_yield |" << dsbp_yield << " Gbp-1Gy-1" << "| dsb_num |" << dsbp_num << G4endl;
  G4cout << "| dsbpp_yield |" << dsbpp_yield << " Gbp-1Gy-1" << "| dsb_num |" << dsbpp_num << G4endl;
  if (dsbn_yield != 0) G4cout << "| SSB/DSB |" << ssb_yield / dsbn_yield << G4endl;
  G4cout << "| energy_deposit |" << total_energy_deposit << " keV" << G4endl;
  G4cout << "| absorbed_dose | " << absorbed_dose << " J/Kg" << G4endl;
  G4cout << "----------------------------------------" << G4endl;
}



