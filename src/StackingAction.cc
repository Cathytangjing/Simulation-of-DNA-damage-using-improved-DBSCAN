//
// Created by hanzawa on 18-5-3.
//

#include "G4StackManager.hh"
#include "G4DNAChemistryManager.hh"

#include "StackingAction.hh"

StackingAction::StackingAction() : G4UserStackingAction() {

}

void StackingAction::NewStage() {
  if (stackManager->GetNTotalTrack() == 0) {
    G4cout << "Physics stage ends, chemistry stage starts..." << G4endl;
    G4DNAChemistryManager::Instance()->Run();
  }
}
