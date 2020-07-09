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
// Jing Tang,Qinfeng Xiao et al. improve the DBSCAN algorithm by utilizing
// the KD-Tree to find neighbors of each site to calculate clustered DNA damage.
// To prevent the ·OH radical that produced one DNA damage site from chain
// propagation throughout the cell, here “kill” the ·OH track after one effective
// contact, thereby stopping the track and its secondaries from producing further damage.
// This work is published:
//
/// \file ITTrackingAction.hh
/// \brief Implementation of the ITTrackingAction class

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
