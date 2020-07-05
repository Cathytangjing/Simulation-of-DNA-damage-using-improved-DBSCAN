//
// Created by hanzawa on 18-5-3.
//

#ifndef DNA_DAMAGE_CELL_STACKINGACTION_H
#define DNA_DAMAGE_CELL_STACKINGACTION_H

#include "G4UserStackingAction.hh"
#include "globals.hh"

class StackingAction : public G4UserStackingAction {
 public:
  StackingAction();
  ~StackingAction() override = default;
  void NewStage() override;
};

#endif //DNA_DAMAGE_CELL_STACKINGACTION_H
