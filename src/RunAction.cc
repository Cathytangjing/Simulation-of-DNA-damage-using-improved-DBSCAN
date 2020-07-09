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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
#include <cstdio>

#include <string>
#include <sstream>
#include <chrono>

#include <G4LogicalVolumeStore.hh>
//#include "Analysis.hh"
#include "RunAction.hh"
#include "RunInitObserver.hh"
#include "RunActionMessenger.hh"
#include "G4AccumulableManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"

#include "G4StateManager.hh"
#include "ExceptionHandler.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::string operator*(const std::string &a, const std::string &b) {
  return a + b;
}

RunAction::RunAction()
    : G4UserRunAction(),
      energy_deposition_(0.0),
      ssb_num_(0),
      cssbn_num_(0),
      dsbn_num_(0),
      dsbp_num_(0),
      dsbpp_num_(0),
      distribution_("") {
  fFileName = "clusters_output";
  run_messenger_ = new RunActionMessenger(this);

  G4AccumulableManager *accumulable_manager = G4AccumulableManager::Instance();
  accumulable_manager->RegisterAccumulable(energy_deposition_);
  accumulable_manager->RegisterAccumulable(cssbn_num_);
  accumulable_manager->RegisterAccumulable(dsbn_num_);
  accumulable_manager->RegisterAccumulable(dsbp_num_);
  accumulable_manager->RegisterAccumulable(dsbpp_num_);
  accumulable_manager->RegisterAccumulable(distribution_);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction() {
  delete run_messenger_;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::map<int, int> RunAction::String_Distribution(const std::string &str) {
  std::cout << "String_Distribution: " << str << std::endl;
  std::map<int, int> temp;
  std::stringstream ss(str);
  std::vector<int> sizes;
  int temp_num;
  while (ss >> temp_num) {
    sizes.push_back(temp_num);
  }

  for (auto i : sizes) {
    temp[i] = temp[i] + 1;
  }
  std::cout << "String_Distribution----------" << std::endl;
  for (auto &i : temp) {
    std::cout << i.first << " - " << i.second << std::endl;
  }
  return temp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run *run) {
  G4cout << "### Run " << run->GetRunID() << " starts." << G4endl; // tj modified
  G4cout << "### Events" << G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed() << " to be processed."
         << G4endl;
  //
  RunInitManager::Instance()->Initialize();

  G4AccumulableManager *accumulable_manager = G4AccumulableManager::Instance();
  accumulable_manager->Reset();

  G4StateManager::GetStateManager()->SetExceptionHandler(new ExceptionHandler());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run *run) {
  G4int nofEvents = run->GetNumberOfEvent();
  G4AccumulableManager *accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  G4double total_energy_deposition = energy_deposition_.GetValue();
  G4int total_ssb_num = ssb_num_.GetValue();
  G4int total_cssbn_num = cssbn_num_.GetValue();
  G4int total_dsbn_num = dsbn_num_.GetValue();
  G4int total_dsbp_num = dsbp_num_.GetValue();
  G4int total_dsbpp_num = dsbpp_num_.GetValue();


  std::string total_size_distribution = distribution_.GetValue();
  std::map<int, int> TotalSizeDistribution = String_Distribution(total_size_distribution);

  if (nofEvents != 0) {
    G4double absorbed_dose = (total_energy_deposition / nofEvents / joule)
        / (G4LogicalVolumeStore::GetInstance()->GetVolume("Target")->GetMass() / kg);

    G4double ssb_yield = static_cast<double>(total_ssb_num) / nofEvents / absorbed_dose / 6;
    G4double cssbn_yield = static_cast<double>(total_cssbn_num) / nofEvents / absorbed_dose / 6;
    G4double dsbn_yield = static_cast<double>(total_dsbn_num) / nofEvents / absorbed_dose / 6;
    G4double dsbp_yield = static_cast<double>(total_dsbp_num) / nofEvents / absorbed_dose / 6;
    G4double dsbpp_yield = static_cast<double>(total_dsbpp_num) / nofEvents / absorbed_dose / 6;


    if (IsMaster()) {
      G4cout
          << G4endl
          << "--------------------End of Global Run-----------------------" << G4endl;

      G4cout << "-----------------RESULT-----------------" << G4endl;
      G4cout << "| ssb_yield |" << ssb_yield << " Gbp-1Gy-1"
             << "| ssb_num |" << static_cast<double>(total_ssb_num) / nofEvents << G4endl;
      G4cout << "| cssbn_yield |" << cssbn_yield << " Gbp-1Gy-1"
             << "| cssbn_num |" << static_cast<double>(total_cssbn_num) / nofEvents << G4endl;
      G4cout << "| dsbn_yield |" << dsbn_yield << " Gbp-1Gy-1"
             << "| dsbn_num |" << static_cast<double>(total_dsbn_num) / nofEvents << G4endl;
      G4cout << "| dsbp_yield |" << dsbp_yield << " Gbp-1Gy-1"
             << "| dsbp_num |" << static_cast<double>(total_dsbp_num) / nofEvents << G4endl;
      G4cout << "| dsbpp_yield |" << dsbpp_yield << " Gbp-1Gy-1"
             << "| dsbpp_num |" << static_cast<double>(total_dsbpp_num) / nofEvents << G4endl;
      G4cout << "| SSB/DSB | " << ssb_yield / dsbn_yield << G4endl;
      G4cout << "| energy_deposit |" << total_energy_deposition / nofEvents / keV << " keV" << G4endl;
      G4cout << "absorbed_dose " << absorbed_dose << " J/Kg" << G4endl;
    } else {
      G4cout
          << G4endl
          << "--------------------End of Local Run------------------------";
    }

    G4cout << G4endl << " The run consists of " << nofEvents << " " << "events." << G4endl;

    if (IsMaster()) {
      std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
      std::time_t now_c = std::chrono::system_clock::to_time_t(now);
      std::stringstream string_stream;
      string_stream << std::put_time(std::localtime(&now_c), "%F_%T_");
      std::string path0;
      string_stream >> path0;
    }
  } else {
    G4cout << "None events consisted in the run. " << G4endl;
  }

  G4cout << "Run " << run->GetRunID() << " finished." << G4endl;

}


