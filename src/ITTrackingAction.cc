//
// Created by hanzawa on 18-5-3.
//

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
//    G4cout<< molecule_name <<" track ID_origin " << track->GetTrackID() << G4endl;
    G4bool result = the_event_action_->clustering_->RegisterChemicalDamage(track->GetPosition());
    if (result) {
      G4Track *track_var = const_cast<G4Track *>(track);
      track_var->SetTrackStatus(G4TrackStatus::fKillTrackAndSecondaries);
//        G4cout<< molecule_name <<" track ID " << track_var->GetTrackID() << G4endl;
    }
  } else if (target_volume != the_volume) {
    G4Track *track_var = const_cast<G4Track *>(track);
    track_var->SetTrackStatus(G4TrackStatus::fKillTrackAndSecondaries);
//    G4EventManager::GetEventManager()->GetTrackingManager()->GetTrack()->SetTrackStatus(G4TrackStatus::fStopAndKill);

  }
}

void ITTrackingAction::PostUserTrackingAction(const G4Track *track) {
//  G4TrackManyList* allTrackList = G4ITTrackHolder::Instance()->GetMainList();
//  if (allTrackList->size() == 0) {
//    G4cout << allTrackList->size() << G4endl;
//  }
}
