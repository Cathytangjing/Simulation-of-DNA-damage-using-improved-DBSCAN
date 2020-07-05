//
// Created by hanzawa on 18-5-3.
//

#ifndef DNA_DAMAGE_CELL_ITTRACKINGACTION_H
#define DNA_DAMAGE_CELL_ITTRACKINGACTION_H

#include "G4UserTrackingAction.hh"

#include "EventAction.hh"

class ITTrackingAction : public G4UserTrackingAction {
 public:
  ITTrackingAction();
  ~ITTrackingAction() override = default;
  void PreUserTrackingAction(const G4Track *track) override;
  void PostUserTrackingAction(const G4Track *track) override;

 private:
  EventAction *the_event_action_;
};

#endif //DNA_DAMAGE_CELL_ITTRACKINGACTION_H
