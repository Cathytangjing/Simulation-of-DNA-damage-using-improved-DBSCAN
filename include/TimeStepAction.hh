//
// Created by hanzawa on 18-5-28.
//

#ifndef DNA_DAMAGE_CELL_TIMESTEPACTION_H_
#define DNA_DAMAGE_CELL_TIMESTEPACTION_H_

#include "G4UserTimeStepAction.hh"

class TimeStepAction : public G4UserTimeStepAction {
 public:
  TimeStepAction();
  virtual ~TimeStepAction();
  TimeStepAction(const TimeStepAction &other);
  TimeStepAction &operator=(const TimeStepAction &other);

  virtual void StartProcessing();

  virtual void UserPreTimeStepAction();
  virtual void UserPostTimeStepAction();

  virtual void UserReactionAction(const G4Track &,
                                  const G4Track &,
                                  const std::vector<G4Track *> *);

  virtual void EndProcessing();
};

#endif // DNA_DAMAGE_CELL_TIMESTEPACTION_H_
