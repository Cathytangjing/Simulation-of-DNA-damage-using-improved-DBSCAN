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
// Authors: Henri Payno and Yann Perrot
//
// $Id$
//
/// \file ClusteringAlgo.hh
/// \brief Definition of the ClustreringAlgo class

#ifndef ClusteringAlgo_H
#define ClusteringAlgo_H 1

#include "ClusterSBPoints.hh"
#include "SBPoint.hh"

//#include "tree.h"

#include <map>

//typedef KD::Core<3, SBPoint> CORE;

class ClusteringAlgoMessenger;

class ClusteringAlgo {
 public:

  ClusteringAlgo(G4double pEps, G4int pMinPts, G4double pSPointsProb, G4double pChemicalPointsProb,
                 G4double pEMinDamage, G4double pEMaxDamage, G4double reaction_rate);
  ~ClusteringAlgo();

  // Get Set methods
  G4double GetEps() {
    return energy_;
  };
  void SetEps(G4double val) {
    energy_ = val;
  };
  G4int GetMinPts() {
    return min_pts_;
  };
  void SetMinPts(G4int val) {
    min_pts_ = val;
  };
  G4double GetSPointsProb() {
    return damage_prob_;
  };
  void SetSPointsProb(G4double val) {
    damage_prob_ = val;
  };
  G4double GetChemicalPointsProb() {
    return chemical_damage_prob_;
  };
  void SetChemicalPointsProb(G4double val) {
    chemical_damage_prob_ = val;
  };
  G4double GetEMinDamage() {
    return min_damage_energy_;
  };
  void SetEMinDamage(G4double val) {
    min_damage_energy_ = val;
  };
  G4double GetEMaxDamage() {
    return max_damage_energy_;
  };
  void SetEMaxDamage(G4double val) {
    max_damage_energy_ = val;
  };
  G4double GetReactionRate() {
    return reaction_rate_;
  }
  void SetReactionRate(G4double val) {
    reaction_rate_ = val;
  };

  // Register a damage (position, edep)
  void RegisterPhysicalDamage(G4ThreeVector, G4double);
  G4bool RegisterChemicalDamage(G4ThreeVector pPos);

  // Clustering Algorithm
  std::map<G4int, G4int> RunClustering();

  // Clean all data structures
  void Purge();
  G4bool IsCSSBN(std::vector<ClusterSBPoints *>::const_iterator it) const;
  G4bool IsDSBN(std::vector<ClusterSBPoints *>::const_iterator it) const;
  G4bool IsDSBP(std::vector<ClusterSBPoints *>::const_iterator it) const;
  G4bool IsDSBPP(std::vector<ClusterSBPoints *>::const_iterator it) const;
  // Return the number of simple break
  G4int GetSSB() const;
  // Return the number of complex simple break
  G4int GetComplexSSB() const;
  // Return the number of double strand break
  G4int GetDSB() const;
  // Return the number of complex simple break with new method
  G4int GetCSSBN() const;
  // Return the number of double strand break with new method, only size = 2
  G4int GetDSBN() const;
  //Return the number of double strand break plus
  G4int GetDSBP() const;
  //Return the number of double strand break plus plus
  G4int GetDSBPP() const;

  inline G4double GetRuningTime() const {
    return runing_time_;
  }

  // Return a map representing cluster size distribution
  // first G4int : cluster size (1 = SSB)
  // second G4int : counts
  std::map<G4int, G4int> GetClusterSizeDistribution();

 private:

  // Functions to check if SB candidate
  G4bool IsInSensitiveArea();
  G4bool IsInChemicalSensitiveArea();
  G4bool IsEdepSufficient(G4double);
  G4bool IsReactionRate();

  // Check if a SB point can be merged to a cluster, and do it
  bool FindCluster(SBPoint *pPt);
  // Check if two points can be merged
  bool AreOnTheSameCluster(G4ThreeVector, G4ThreeVector, G4double);
  // Merge clusters
  void MergeClusters();
  // Add SSB to clusters
  void IncludeUnassociatedPoints();

  // Parameters to run clustering algorithm
  G4double energy_;         // distance to merge SBPoints
  G4int min_pts_;         // number of SBPoints to create a cluster
  G4double damage_prob_; // probability for a point to be in the sensitive area
  G4double chemical_damage_prob_;
  G4double min_damage_energy_;  // min energy to create a damage
  G4double max_damage_energy_;  // energy to have a probability to create a damage = 1
  G4double reaction_rate_;

  // Data structure containing all SB points
  std::vector<SBPoint *> set_of_points_;
  // Data structure containing all clusters
  std::vector<ClusterSBPoints *> clusters_;
  // ID of the next SB point
  unsigned int next_point_id_;
  G4double runing_time_;

  ClusteringAlgoMessenger *messenger_;
};

#endif

