//
// Created by hanzawa on 18-5-3.
//

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
