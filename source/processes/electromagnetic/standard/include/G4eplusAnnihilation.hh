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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eplusAnnihilation
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 02.08.2004
//
// Modifications:
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivanchenko)
// 15-03-05 Update interface according to changings 
//           in G4VEmProcess (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivanchenko)
// 04-05-05 Make class to be default (V.Ivanchenko)
// 04-12-05 SetProposedKineticEnergy(0.) for annihilated positron (mma)
// 09-08-06 add SetModel(G4VEmModel*) (mma)
// 12-09-06, move SetModel(G4VEmModel*) in G4VEmProcess (mma)
//
//
// Class Description:
//
// This class manages the process of e+ annihilation at rest and on fly
// It is possible to enable ApplyCuts and Entanglement options via
// G4EmParameters class using UI commands or C++ interface.
// EM splitting or Russian roulette are allowed is corresponding options
// are enabled.
//
// -------------------------------------------------------------------
//

#ifndef G4eplusAnnihilation_h
#define G4eplusAnnihilation_h 1

#include "G4VEmProcess.hh"

class G4ParticleDefinition;
class G4VPositronAtRestModel;

class G4eplusAnnihilation : public G4VEmProcess
{

public:

  explicit G4eplusAnnihilation(const G4String& name = "annihil");

  ~G4eplusAnnihilation() override;

  G4bool IsApplicable(const G4ParticleDefinition& p) final;

  G4VParticleChange* AtRestDoIt(
                             const G4Track& track,
                             const G4Step& stepData) override;

  G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition
                            ) override;

  // print documentation in html format
  void ProcessDescription(std::ostream&) const override;

  G4eplusAnnihilation & operator=(const G4eplusAnnihilation &right) = delete;
  G4eplusAnnihilation(const G4eplusAnnihilation&) = delete;

protected:

  void InitialiseProcess(const G4ParticleDefinition*) override;

  // Print out of the class parameters
  void StreamProcessInfo(std::ostream& outFile) const override;

private:

  G4VPositronAtRestModel* f2GammaAtRestModel{nullptr};
  G4VPositronAtRestModel* f3GammaAtRestModel{nullptr};
  G4int fEntanglementModelID;
  G4bool isInitialised{false};
  G4bool fEntangled{false};
  G4bool fApplyCuts{false};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
