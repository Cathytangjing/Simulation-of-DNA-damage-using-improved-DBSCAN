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
// Jing Tang and Qinfeng Xiao modified on 18-5-3.
///

#ifndef DNA_DAMAGE_CELL_ITTRACKINGINTERACTIVITY_H
#define DNA_DAMAGE_CELL_ITTRACKINGINTERACTIVITY_H

#include <vector>

#include "G4ITTrackingInteractivity.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VTrajectory.hh"

class ITTrackingInteractivity : public G4ITTrackingInteractivity {
 public:
  ITTrackingInteractivity();
  ~ITTrackingInteractivity() override;

  void Initialize() override;
  void StartTracking(G4Track *track) override;
  void AppendStep(G4Track *track, G4Step *step) override;
  void EndTracking(G4Track *track) override;
  void Finalize() override;

  inline void SetUserAction(G4UserTrackingAction *user_tracking_action) {
    user_tracking_action_ = user_tracking_action;
  }

  inline G4UserTrackingAction *GetUserTrackingAction() const {
    return user_tracking_action_;
  }

  inline void SetUserAction(G4UserSteppingAction *user_stepping_action) {
    user_stepping_action_ = user_stepping_action;
  }

  inline G4UserSteppingAction *GetUserSteppingAction() const {
    return user_stepping_action_;
  }

 private:
  G4UserTrackingAction *user_tracking_action_;
  G4UserSteppingAction *user_stepping_action_;

  std::vector<G4VTrajectory *> trajectories_;

  G4int store_trajectory_;
};

#endif //DNA_DAMAGE_CELL_ITTRACKINGINTERACTIVITY_H
