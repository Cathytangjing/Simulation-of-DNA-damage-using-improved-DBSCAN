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
//

#include "G4EventManager.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4TrackingInformation.hh"
#include "G4VTrajectory.hh"
#include "G4Trajectory.hh"
#include "G4SmoothTrajectory.hh"
#include "G4RichTrajectory.hh"
#include "G4IT.hh"
#include "G4Event.hh"
#include "G4VSteppingVerbose.hh"
#include "G4VisManager.hh"
#include "G4ITTrackingInteractivity.hh"

#include "ITTrackingInteractivity.hh"

class G4Trajectory_Lock {
  friend class ITTrackingInteractivity;

  G4Trajectory_Lock() : fpTrajectory(0) { ; }

  ~G4Trajectory_Lock() { ; }

  G4VTrajectory *fpTrajectory;
};

ITTrackingInteractivity::ITTrackingInteractivity() : user_tracking_action_(nullptr), user_stepping_action_(nullptr) {
  store_trajectory_ = 0;

  fVerboseLevel = 0;
}

ITTrackingInteractivity::~ITTrackingInteractivity() {
  G4EventManager *event_manager = G4EventManager::GetEventManager();

  if (event_manager) {
    G4UserTrackingAction *standard_tracking_action = event_manager->GetUserTrackingAction();
    G4UserSteppingAction *standard_stepping_action = event_manager->GetUserSteppingAction();

    if (user_tracking_action_ != standard_tracking_action) delete user_tracking_action_;
    if (user_stepping_action_ != standard_stepping_action) delete user_stepping_action_;
  } else {
    delete user_tracking_action_;
    delete user_stepping_action_;
  }
}

void ITTrackingInteractivity::Initialize() {
  G4TrackingManager *tracking_manager =
      G4EventManager::GetEventManager()->GetTrackingManager();
  store_trajectory_ = tracking_manager->GetStoreTrajectory();
  fVerboseLevel = tracking_manager->GetVerboseLevel();
}

void ITTrackingInteractivity::StartTracking(G4Track *track) {

#ifdef G4VERBOSE
  if (fVerboseLevel) {
    TrackBanner(track, "G4ITTrackingManager::StartTracking : ");
  }
#endif

  if (fVerboseLevel > 0 && (G4VSteppingVerbose::GetSilent() != 1))
    TrackBanner(track);

  // Pre tracking user intervention process.
  if (user_tracking_action_ != 0) {
    user_tracking_action_->PreUserTrackingAction(track);
  }
//#ifdef G4_STORE_TRAJECTORY
  G4TrackingInformation *tracking_info = GetIT(track)->GetTrackingInfo();
  G4Trajectory_Lock *trajectory_lock =
      tracking_info->GetTrajectory_Lock();

  // Construct a trajectory if it is requested
  if (store_trajectory_ && (!trajectory_lock)) {
    trajectory_lock = new G4Trajectory_Lock();
    tracking_info->SetTrajectory_Lock(trajectory_lock);

    G4VTrajectory *trajectory = nullptr;

    // default trajectory concrete class object
    switch (store_trajectory_) {
      default:
      case 1: trajectory = new G4Trajectory(track);
        break;
      case 2: trajectory = new G4SmoothTrajectory(track);
        break;
      case 3: trajectory = new G4RichTrajectory(track);
        break;
    }
    trajectory_lock->fpTrajectory = trajectory;
  }
}

void ITTrackingInteractivity::AppendStep(G4Track *track, G4Step *step) {
  if (user_stepping_action_)
    user_stepping_action_->UserSteppingAction(step);

  if (store_trajectory_) {
    G4TrackingInformation *tracking_info =
        GetIT(track)->GetTrackingInfo();
    G4Trajectory_Lock *trajectory_lock =
        tracking_info->GetTrajectory_Lock();
    trajectory_lock->fpTrajectory->AppendStep(step);
  }
}

void ITTrackingInteractivity::EndTracking(G4Track *track) {
#ifdef G4VERBOSE
  if (fVerboseLevel) {
    TrackBanner(track, "G4ITTrackingManager::EndTracking : ");
  }
#endif
  // Post tracking user intervention process.
  if (user_tracking_action_) {
    user_tracking_action_->PostUserTrackingAction(track);
  }

//#ifdef G4_STORE_TRAJECTORY
  G4TrackingInformation *tracking_info = GetIT(track)->GetTrackingInfo();
  G4Trajectory_Lock *trajectory_lock = tracking_info->GetTrajectory_Lock();

  if (trajectory_lock) {
    G4VTrajectory *&trajectory = trajectory_lock->fpTrajectory;

    if (store_trajectory_ && trajectory) {

#ifdef G4VERBOSE
      if (fVerboseLevel > 10) trajectory->ShowTrajectory();
#endif

      G4TrackStatus is_top = track->GetTrackStatus();

      if (trajectory && (is_top != fStopButAlive) && (is_top != fSuspend)) {
        G4Event *currentEvent = G4EventManager::GetEventManager()
            ->GetNonconstCurrentEvent();

        if (currentEvent) {
          G4TrajectoryContainer *trajectory_ontainer = currentEvent
              ->GetTrajectoryContainer();

          if (!trajectory_ontainer) {
            trajectory_ontainer = new G4TrajectoryContainer;
            currentEvent->SetTrajectoryContainer(trajectory_ontainer);
          }
          trajectory_ontainer->insert(trajectory);
        } else {
          trajectories_.push_back(trajectory);
        }
      }
    }
      // Destruct the trajectory if it was created
    else if ((!store_trajectory_) && trajectory) {
      delete trajectory;
      trajectory = nullptr;
    }
    delete trajectory_lock;
    tracking_info->SetTrajectory_Lock(0);
  }
}

void ITTrackingInteractivity::Finalize() {
  for (std::vector<G4VTrajectory *>::iterator it = trajectories_.begin();
       it != trajectories_.end(); it++) {
    G4VisManager::GetConcreteInstance()->Draw(**it);
  }
}
