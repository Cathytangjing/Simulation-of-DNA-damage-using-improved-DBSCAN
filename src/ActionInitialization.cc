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
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "G4DNAChemistryManager.hh"
#include "G4Scheduler.hh"
#include "G4SystemOfUnits.hh"

#include "ActionInitialization.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "StackingAction.hh"

#include "ITTrackingInteractivity.hh"
#include "ITTrackingAction.hh"
#include "TimeStepAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization()
    : G4VUserActionInitialization() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const {
  SetUserAction(new RunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const {
  SetUserAction(new RunAction());
  SetUserAction(new PrimaryGeneratorAction());
  SetUserAction(new EventAction());
  SetUserAction(new SteppingAction());
  SetUserAction(new StackingAction());

  if (G4DNAChemistryManager::IsActivated()) {
    G4Scheduler::Instance()->SetVerbose(1);
    G4Scheduler::Instance()->SetUserAction(new TimeStepAction());
    // Uncomment and set to stop chemistry stage after:
    // ...given number of time steps
    //G4Scheduler::Instance()->SetMaxNbSteps(1000);
    G4Scheduler::Instance()->SetEndTime(2.5 * nanosecond);

    ITTrackingInteractivity *tracking_interactivity = new ITTrackingInteractivity();
    tracking_interactivity->SetUserAction(new ITTrackingAction());
    G4Scheduler::Instance()->SetInteractivity(tracking_interactivity);
  }

}

