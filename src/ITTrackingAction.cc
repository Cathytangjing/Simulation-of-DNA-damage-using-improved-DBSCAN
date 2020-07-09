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
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Jing Tang,Qinfeng Xiao et al. improve the DBSCAN algorithm by utilizing
// the KD-Tree to find neighbors of each site to calculate clustered DNA damage.
// To prevent the ·OH radical that produced one DNA damage site from chain
// propagation throughout the cell, here “kill” the ·OH track after one effective
// contact, thereby stopping the track and its secondaries from producing further damage.
// This work is published:
//
/// \file ITTrackingAction.hh
/// \brief Implementation of the ITTrackingAction class

#include <G4LogicalVolumeStore.hh>
#include "G4Track.hh"
#include "G4EventManager.hh"
#include "G4Molecule.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4MoleculeDefinition.hh"
#include "G4Scheduler.hh"
#include "Randomize.hh"

#include "ITTrackingAction.hh"
#include "EventAction.hh"
#include "ClusteringAlgo.hh"

ITTrackingAction::ITTrackingAction() : G4UserTrackingAction() {
  the_event_action_ = (EventAction *) G4EventManager::GetEventManager()->
      GetUserEventAction();
}

void ITTrackingAction::PreUserTrackingAction(const G4Track *track) {
//  G4cout << "Track ID -- " << track->GetTrackID()
//         << "  Molecule Name -- " << GetMolecule(track)->GetName()
//         << "  Track Time -- " << G4BestUnit(track->GetGlobalTime(), "Time")
//         << "  Position -- " << track->GetPosition() << G4endl;

  G4LogicalVolume *target_volume = G4LogicalVolumeStore::GetInstance()->GetVolume("Target");
  G4LogicalVolume *the_volume = track->GetVolume()->GetLogicalVolume();

  G4String molecule_name = GetMolecule(track)->GetName();

  if (molecule_name == "OH^0" && target_volume == the_volume) {
    G4bool result = the_event_action_->clustering_->RegisterChemicalDamage(track->GetPosition());
    if (result) {
      G4Track *track_var = const_cast<G4Track *>(track);
      track_var->SetTrackStatus(G4TrackStatus::fKillTrackAndSecondaries);
    }
  } else if (target_volume != the_volume) {
    G4Track *track_var = const_cast<G4Track *>(track);
    track_var->SetTrackStatus(G4TrackStatus::fKillTrackAndSecondaries);
  }
}

void ITTrackingAction::PostUserTrackingAction(const G4Track *track) {
}
