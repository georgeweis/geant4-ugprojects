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
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

// MY CLASS: added around line 162

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"
#include "G4GeorgeSolid.hh"
#include "G4GeorgeNurbs.hh"
#include "SavedTorusData.hh"
#include "G4Tubs.hh"


namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  G4double env_sizeX = 40 * cm, env_sizeY = 20 * cm, env_sizeZ = 40 * cm;

  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
  // World material
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");



  // =======================================================================
  //   World dimensions (800 mm x 800 mm x 800 mm)
  G4double world_half = 41 * cm;

  auto solidWorld = new G4Box("World", world_half, world_half, world_half);
  auto logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  auto physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false, 0, false);

  // envelope
  auto solidEnv = new G4Box("Envelope",  // its name
                            env_sizeX, env_sizeY, env_sizeZ);  // its size

  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
                                    env_mat,  // its material
                                    "Envelope");  // its name


  // Torus NURBS definition
  SavedTorusData torusData;
  std::vector<std::vector<G4ThreeVector>> torusControlPts = torusData.getTorusControlPts();
  std::vector<std::vector<G4double>> torusWeights = torusData.getTorusWeights();

  std::vector<G4double> torusKnotsU = {0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1};
  std::vector<G4double> torusKnotsV = {0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1};

  int torusDegreeU = 2;
  int torusDegreeV = 2;

  G4GeorgeNurbs* torusNurbs = new G4GeorgeNurbs("NurbsTorus",
                                                 torusControlPts, torusWeights,
                                                 torusKnotsU, torusKnotsV,
                                                 torusDegreeU, torusDegreeV);





  // Logical and physical placement of the torus
  G4Material* torus_mat = nist->FindOrBuildMaterial("G4_Fe");
  auto logicTorus = new G4LogicalVolume(torusNurbs, torus_mat, "NurbsTorusLV");

  new G4PVPlacement(nullptr, G4ThreeVector(), logicTorus, "NurbsTorus", logicWorld, false, 0, true);


  //==========================================================

//  // World dimensions (800 mm x 800 mm x 800 mm)
//  G4double world_half = 4;
//
//  auto solidWorld = new G4Box("World", world_half, world_half, world_half);
//  auto logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
//  auto physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false, 0, false);
//
//  // envelope
//  auto solidEnv = new G4Box("Envelope",  // its name
//                            env_sizeX, env_sizeY, env_sizeZ);  // its size
//
//  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
//                                    env_mat,  // its material
//                                    "Envelope");  // its name
//  // instantiaion of my nurbs cylinder
//  G4double r = 1.0;
//  G4double w = std::sqrt(2) / 2.0;
//
//  std::vector<std::vector<G4ThreeVector>> controlPts = {
//    {
//      G4ThreeVector{0.0, -r,  r},
//      G4ThreeVector{ r,   -r,  r},
//      G4ThreeVector{ r,    0.0, r},
//      G4ThreeVector{ r,    r,  r},
//      G4ThreeVector( 0.0,  r,  r)
//  },
//  {
//    G4ThreeVector{ 0.0, -r, -r},
//    G4ThreeVector{ r,   -r, -r},
//    G4ThreeVector{ r,    0.0, -r},
//    G4ThreeVector{ r,    r, -r},
//    G4ThreeVector{ 0.0,  r, -r}
//  }};
//
//  std::vector<std::vector<G4double>> weights = {
//    {1, w, 1, w, 1},
//    {1, w, 1, w, 1}
//  };
//
//  std::vector<G4double> knotsU = {0, 0, 1, 1};
//  std::vector<G4double> knotsV = {0, 0, 0, 0.5, 0.5, 1, 1, 1};
//
//  int degreeU = 1;
//  int degreeV = 2;
//
//  G4GeorgeNurbs* NurbsCylinder = new G4GeorgeNurbs("Nurbs Half Cylinder",
//                                                   controlPts,
//                                                   weights,
//                                                   knotsU,
//                                                   knotsV,
//                                                   degreeU,
//                                                   degreeV);
//  //analytical cylinder
//  double radius = 1;
//  double halfLength = 10;
//  auto G4Cylinder = new G4Tubs("CylZ", 0, radius, halfLength, 0., 360.*deg);
//
//
//
//
//
//  // Logical and physical placement of the torus
//  G4Material* torus_mat = nist->FindOrBuildMaterial("G4_Fe");
//  auto logicTorus = new G4LogicalVolume(NurbsCylinder, torus_mat, "NurbsTorusLV");
//
//
//
//
//  new G4PVPlacement(nullptr, G4ThreeVector(), logicTorus, "NurbsTorus", logicWorld, false, 0, true);

  //==========================================================

  fScoringVolume = logicTorus;

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
