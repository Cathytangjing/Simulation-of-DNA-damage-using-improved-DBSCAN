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
//#include <iomanip>

#include "EventAction.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4Run.hh"

#include "Analysis.hh"
#include "ClusteringAlgo.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "SizeDistribution.hh"
#include "csv.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction() : G4UserEventAction() {
  //default parameter values
  event_energy_ = 0.0;
  event_step_length_ = 0.0;

  //yield_factor_ = 0.95;

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
//  for (auto &i : sizeDistribution) {
//    std::cout << i.first << " - " << i.second << std::endl;
//  }

  G4int ssb_num = clustering_->GetSSB();
  G4int cssb_num = clustering_->GetComplexSSB();
  G4int dsb_num = clustering_->GetDSB();
  G4int cssbn_num = clustering_->GetCSSBN();
  G4int dsbn_num = clustering_->GetDSBN();
  G4int dsbp_num = clustering_->GetDSBP();
  G4int dsbpp_num = clustering_->GetDSBPP();

  // before compute SD(Stardard Diviation)
  G4double ssb_num_mean = 875.453;  // change according to different initial energy,now for 0.8MeV,1000,2.5ns,2020.2.28
  G4double dsb_num_mean = 110.539;
  G4double cssb_num_mean = 82.081;
  G4double dsbn_num_mean = 71.483;
  G4double dsbp_num_mean = 32.754;
  G4double dsbpp_num_mean = 6.3;

  G4double pow_ssb_n_diff = pow((ssb_num - ssb_num_mean), 2);
  G4double pow_dsb_n_diff = pow((dsb_num - dsb_num_mean), 2);
  G4double pow_cssb_n_diff = pow((cssb_num - cssb_num_mean), 2);
  G4double pow_dsbn_n_diff = pow((dsbn_num - dsbn_num_mean), 2);
  G4double pow_dsbp_n_diff = pow((dsbp_num - dsbp_num_mean), 2);
  G4double pow_dsbpp_n_diff = pow((dsbpp_num - dsbpp_num_mean), 2);

  RunAction *run_action =
      (RunAction *) G4RunManager::GetRunManager()->GetUserRunAction();//(RunAction*)G4RunManager::GetRunManager()->GetUserRunAction();
  run_action->AddEnergyDeposition(event_energy_);
  run_action->AddSingleStrandBreak(ssb_num);
  run_action->AddComplexSingleStrandBreak(cssb_num);
  run_action->AddDoubleStrandBreak(dsb_num);
  run_action->AddComplexSingleStrandBreakNew(cssbn_num);
  run_action->AddDoubleStrandBreakNew(dsbn_num);
  run_action->AddDoubleStrandBreakPlus(dsbp_num);
  run_action->AddDoubleStrandBreakPlusPlus(dsbpp_num);

  run_action->AddPowSSBNumDiff(pow_ssb_n_diff);
  run_action->AddPowCSSBNumDiff(pow_cssb_n_diff);
  run_action->AddPowDSBNumDiff(pow_dsb_n_diff);
  run_action->AddPowDSBnNumDiff(pow_dsbn_n_diff);
  run_action->AddPowDSBpNumDiff(pow_dsbp_n_diff);
  run_action->AddPowDSBppNumDiff(pow_dsbpp_n_diff);

  for (auto it = sizeDistribution.begin(); it != sizeDistribution.end(); it++) {
    for (int i = 0; i < it->second; i++) {
      run_action->AddDistribution(it->first);
    }
  }

  G4double total_energy_deposit = event_energy_ / keV;
  G4double absorbed_dose =
      (event_energy_ / joule) / (G4LogicalVolumeStore::GetInstance()->GetVolume("Target")->GetMass() / kg);

//  G4double ssb_yield = static_cast<G4double>(ssb_num)*(total_energy_deposit/absorbed_dose)/total_energy_deposit/6.0*yield_factor_;
//  G4double cssb_yield = static_cast<G4double>(cssb_num)*(total_energy_deposit/absorbed_dose)/total_energy_deposit/6.0*yield_factor_;
//  G4double dsb_yield = static_cast<G4double>(dsb_num)*(total_energy_deposit/absorbed_dose)/total_energy_deposit/6.0*yield_factor_;

  G4double ssb_yield = static_cast<double>(ssb_num) / absorbed_dose / (6);//5.985
  G4double cssb_yield = static_cast<double>(cssb_num) / absorbed_dose / (6);
  G4double dsb_yield = static_cast<double>(dsb_num) / absorbed_dose / (6);
  G4double cssbn_yield = static_cast<double>(cssbn_num) / absorbed_dose / (6);
  G4double dsbn_yield = static_cast<double>(dsbn_num) / absorbed_dose / (6);
  G4double dsbp_yield = static_cast<double>(dsbp_num) / absorbed_dose / (6);
  G4double dsbpp_yield = static_cast<double>(dsbpp_num) / absorbed_dose / (6);

  G4double ssb_over_length = static_cast<double>(ssb_num) / (event_step_length_ / um);
  G4double cssb_over_length = static_cast<double>(cssb_num) / (event_step_length_ / um);
  G4double dsb_over_length = static_cast<double>(dsb_num) / (event_step_length_ / um);
  G4double cssbn_over_length = static_cast<double>(cssbn_num) / (event_step_length_ / um);
  G4double dsbn_over_length = static_cast<double>(dsbn_num) / (event_step_length_ / um);
  G4double dsbp_over_length = static_cast<double>(dsbp_num) / (event_step_length_ / um);
  G4double dsbpp_over_length = static_cast<double>(dsbp_num) / (event_step_length_ / um);

  G4cout << "-----------------RESULT-----------------" << G4endl;
  G4cout << "| ssb_yield |" << ssb_yield << " Gbp-1Gy-1" << "| ssb_num |" << ssb_num << G4endl;
  G4cout << "| cssb_yield |" << cssb_yield << " Gbp-1Gy-1" << "| cssb_num |" << cssb_num << G4endl;
  G4cout << "| dsb_yield |" << dsb_yield << " Gbp-1Gy-1" << "| dsb_num |" << dsb_num << G4endl;
  G4cout << "| cssbn_yield |" << cssbn_yield << " Gbp-1Gy-1" << "| cssbn_num |" << cssbn_num << G4endl;
  G4cout << "| dsbn_yield |" << dsbn_yield << " Gbp-1Gy-1" << "| dsbn_num |" << dsbn_num << G4endl;
  G4cout << "| dsbp_yield |" << dsbp_yield << " Gbp-1Gy-1" << "| dsb_num |" << dsbp_num << G4endl;
  G4cout << "| dsbpp_yield |" << dsbpp_yield << " Gbp-1Gy-1" << "| dsb_num |" << dsbpp_num << G4endl;
  if (dsb_yield != 0) G4cout << "| SSB/DSB |" << ssb_yield / dsb_yield << G4endl;
  G4cout << "| energy_deposit |" << total_energy_deposit << " keV" << G4endl;
  G4cout << "| absorbed_dose | " << absorbed_dose << " J/Kg" << G4endl;
  G4cout << "| ssb_over_length |" << ssb_over_length << G4endl;
  G4cout << "| cssb_over_length |" << cssb_over_length << G4endl;
  G4cout << "| dsb_over_length |" << dsb_over_length << G4endl;
  G4cout << "| cssbn_over_length |" << cssbn_over_length << G4endl;
  G4cout << "| dsbn_over_length |" << dsbn_over_length << G4endl;
  G4cout << "| dsbp_over_length |" << dsbp_over_length << G4endl;
  G4cout << "| dsbpp_over_length |" << dsbpp_over_length << G4endl;
  G4cout << "----------------------------------------" << G4endl;

//    G4double ssb_yield = static_cast<double>(ssb_num)/absorbed_dose/(3.6*1e12);
//    G4double cssb_yield = static_cast<double>(cssb_num)/absorbed_dose/(3.6*1e12);
//    G4double dsb_yield = static_cast<double>(dsb_num)/absorbed_dose/(3.6*1e12);
//
//    G4cout << "-----------------RESULT-----------------" << G4endl;
//    G4cout << "| ssb_yield |" << ssb_yield << "Da-1Gy-1" << "| ssb_num |" << ssb_num << G4endl;
//    G4cout << "| cssb_yield |" << cssb_yield << "Da-1Gy-1" << "| cssb_num |" << cssb_num << G4endl;
//    G4cout << "| dsb_yield |" << dsb_yield << "Da-1Gy-1" << "| dsb_num |" << dsb_num << G4endl;
//    G4cout << "| SSB/DSB |" << ssb_yield/dsb_yield << G4endl;
//    G4cout << "| energy_deposit |" << total_energy_deposit << " keV" << G4endl;
//    G4cout << "| absorbed_dose | " << absorbed_dose << " J/Kg" << G4endl;
//    G4cout << "----------------------------------------" << G4endl;

//  print out the result of every event to the file
// {
  std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
  std::time_t now_c = std::chrono::system_clock::to_time_t(now);
  std::stringstream string_stream;
  string_stream << std::put_time(std::localtime(&now_c), "%F_%T_");
  std::string path0;
  string_stream >> path0;

  CSVWriter writer("output/event_" + std::to_string(event->GetEventID()) + "_" + path0 + ".csv", std::vector<std::string>({"Attribute", "Value"}));
  writer.WriteSingle("ID");
  writer.WriteSingle(event->GetEventID());
  writer.WriteSingle("Time");
  writer.WriteSingle(clustering_->GetRuningTime());
  writer.WriteSingle("SSB");
  writer.WriteSingle(ssb_num);
  writer.WriteSingle("CSSB");
  writer.WriteSingle(cssb_num);
  writer.WriteSingle("CSSBN");
  writer.WriteSingle(cssbn_num);
  writer.WriteSingle("DSB");
  writer.WriteSingle(dsb_num);
  writer.WriteSingle("DSBN");
  writer.WriteSingle(dsbn_num);
  writer.WriteSingle("DSBP");
  writer.WriteSingle(dsbp_num);
  writer.WriteSingle("DSBPP");
  writer.WriteSingle(dsbpp_num);
  writer.WriteSingle("SSB_Yield");
  writer.WriteSingle(ssb_yield);
  writer.WriteSingle("CSSB_Yield");
  writer.WriteSingle(cssb_yield);
  writer.WriteSingle("CSSBN_Yield");
  writer.WriteSingle(cssbn_yield);
  writer.WriteSingle("DSB_Yield");
  writer.WriteSingle(dsb_yield);
  writer.WriteSingle("DSBN_Yield");
  writer.WriteSingle(dsbn_yield);
  writer.WriteSingle("DSBP_Yield");
  writer.WriteSingle(dsbp_yield);
  writer.WriteSingle("DSBPP_Yield");
  writer.WriteSingle(dsbpp_yield);
  writer.WriteSingle("SSB_Over_DSB");
  if (dsb_yield != 0) {
    writer.WriteSingle(ssb_yield / dsb_yield);
  } else {
    writer.WriteSingle("NaN");
  }
  writer.WriteSingle("energy_deposit");
  writer.WriteSingle(std::to_string(total_energy_deposit) + " keV");
  writer.WriteSingle("absorbed_dose");
  writer.WriteSingle(std::to_string(absorbed_dose) + " J/Kg");
  writer.WriteSingle("ssb_over_length");
  writer.WriteSingle(ssb_over_length);
  writer.WriteSingle("cssb_over_length");
  writer.WriteSingle(cssb_over_length);
  writer.WriteSingle("dsb_over_length");
  writer.WriteSingle(dsb_over_length);
  writer.WriteSingle("cssbn_over_length");
  writer.WriteSingle(cssbn_over_length);
  writer.WriteSingle("dsbn_over_length");
  writer.WriteSingle(dsbn_over_length);
  writer.WriteSingle("dsbp_over_length");
  writer.WriteSingle(dsbp_over_length);
  writer.WriteSingle("dsbpp_over_length");
  writer.WriteSingle(dsbpp_over_length);

//  std::string path = std::string(
//      "/home/nuccathy/work_dir/clustering_chemical_tj_parallel_4_improved DBSCAN(20190805) with flann - xiaoqinfeng/clustering_chemical_tj_parallel/Result/Yield/")
//      + "res_Event_" + path0 + std::to_string(event->GetEventID()) + ".txt";
//
//  FILE *file = std::fopen(path.c_str(), "w");
//
//  fprintf(file, "ssb_yield %f Gbp-1Gy-1\nssb_num %f\n"
//                "dsb_yield %f Gbp-1Gy-1\ndsb_num %f\n"
//                "cssbn_yield %f Gbp-1Gy-1\ncssbn_num %f\n"
//                "dsbn_yield %f Gbp-1Gy-1\ndsbn_num %f\n"
//                "dsbp_yield %f Gbp-1Gy-1\ndsbp_num %f\n"
//                "dsbpp_yield %f Gbp-1Gy-1\ndsbpp_num %f\n"
//                //"ssb_dst_ratio %f\n"
//                "energy_deposit %f keV\nabsorbed_dose %f J/Kg\n",
//          ssb_yield, ssb_num,
//          dsb_yield, dsb_num,
//          cssbn_yield, cssbn_num,
//          dsbn_yield, dsbn_num,
//          dsbp_yield, dsbp_num,
//          dsbpp_yield, dsbpp_num,
//      //ssb_yield / dsb_yield,
//          total_energy_deposit / keV, absorbed_dose);
//  std::fclose(file);

//   std::string path2 = "sizeDis_Event" + path0 + std::to_string(event->GetEventID()) + ".txt";
//
//   FILE* file2 = std::fopen(path2.c_str(), "w");
//
//   for(auto it = sizeDistribution.begin(); it != sizeDistribution.end(); it++) {
//      fprintf(file2,"%d:%8d\n",it->first, it->second);
//   }
//
//    std::fclose(file2);
//
//   }

}


// analysisManager->FillH1(1,ssb_num);
// analysisManager->FillH1(2,cssb_num);
// analysisManager->FillH1(3,dsb_num);

// while ( !sizeDistribution.empty() )
// {
//   analysisManager->FillH1(4,
//                           sizeDistribution.begin()->first,
//                           sizeDistribution.begin()->second);
//   sizeDistribution.erase(sizeDistribution.begin());
// }

// analysisManager->FillH1(5,
//                         (event_energy_/joule)/
//                         (G4LogicalVolumeStore::GetInstance()->
//                             GetVolume("Target")->GetMass()/kg)
// );

