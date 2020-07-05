//
// Created by hanzawa on 18-5-28.
//

#include "TimeStepAction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Scheduler.hh"

TimeStepAction::TimeStepAction() : G4UserTimeStepAction() {
  /**
   * Inform G4ITTimeStepper of the selected minimum time steps
   * eg : from 1 picosecond to 10 picosecond, the minimum time
   * step that the TimeStepper can returned is 0.1 picosecond.
   *
   * Case 1) If the rection model calculates a minimum reaction time
   * bigger than the user defined time step, the reaction model wins
   *
   * Case 2) If an interaction process with the continuous medium
   * calculates a time step less than the selected minimum time step,
   * the interaction process wins
   */

  AddTimeStep(1 * picosecond, 0.1 * picosecond);
  AddTimeStep(10 * picosecond, 1 * picosecond);
  AddTimeStep(100 * picosecond, 10 * picosecond);
  AddTimeStep(1000 * picosecond, 100 * picosecond);
  AddTimeStep(10000 * picosecond, 1000 * picosecond);
}

TimeStepAction::~TimeStepAction() {

}

TimeStepAction::TimeStepAction(const TimeStepAction &other) : G4UserTimeStepAction(other) {

}

TimeStepAction &TimeStepAction::operator=(const TimeStepAction &rhs) {
  if (this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}

void TimeStepAction::StartProcessing() {
// You want to know why the simulation stopped ?
  G4Scheduler::Instance()->WhyDoYouStop();
// At the end of the simulation, information will be printed
// It is better to place this command before the simulation starts
}

void TimeStepAction::UserPreTimeStepAction() {

}

void TimeStepAction::UserPostTimeStepAction() {

}

// Here you can retrieve information related to reactions
void TimeStepAction::UserReactionAction(const G4Track &,
                                        const G4Track &,
                                        const std::vector<G4Track *> *  /*products*/) {
/*
   for (int i = 0 ; i < nbProducts ; i ++)
   {
   G4cout << "Product[" << i << "] : "
   << GetMolecule(products[i])->GetName()
   << G4endl ;
   }
*/
}

void TimeStepAction::EndProcessing() {

}
