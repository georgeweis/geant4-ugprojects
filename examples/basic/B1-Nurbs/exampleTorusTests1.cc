//
// The code below is the exact main that i used to perform all the nurbs surface
// tests once the solid was fully implemented. Ive saved it here in case i need
// to revisit any of these tests.


// second copy for torus tests



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
/// \file exampleB1.cc
/// \brief Main program of the B1 example

//#include "ActionInitialization.hh"
//#include "DetectorConstruction.hh"
//#include "QBBC.hh"
//
//#include "G4RunManagerFactory.hh"
//#include "G4SteppingVerbose.hh"
//#include "G4UIExecutive.hh"
//#include "G4UImanager.hh"
//#include "G4VisExecutive.hh"
//// #include "Randomize.hh"
#include "G4GeorgeNurbs.hh"
#include "G4GeorgeSolid.hh"
#include "G4Box.hh"

#include <iostream>
#include <vector>
#include <string>
//#include <nlopt.hpp>



//using namespace B1;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace georgeTorusTests{
std::vector<std::vector<G4ThreeVector>> getTorusControlPts()
{

  std::vector<std::vector<G4ThreeVector>> controlPts =
  {
    {G4ThreeVector{-125, 0, 0},
    G4ThreeVector{-125, 0, 125},
    G4ThreeVector{0, 0, 125},
    G4ThreeVector{125, 0, 125},
    G4ThreeVector{125, 0, 0},
    G4ThreeVector{125, 0, -125},
    G4ThreeVector{0, 0, -125},
    G4ThreeVector{-125, 0, -125},
    G4ThreeVector{-125, 0, 0}},

    {G4ThreeVector{-125, 100, 0},
    G4ThreeVector{-125, 100, 125},
    G4ThreeVector{0, 100, 125},
    G4ThreeVector{125, 100, 125},
    G4ThreeVector{125, 100, 0},
    G4ThreeVector{125, 100, -125},
    G4ThreeVector{0, 100, -125},
    G4ThreeVector{-125, 100, -125},
    G4ThreeVector{-125, 100, 0}},

    {G4ThreeVector{-225, 100, 0},
    G4ThreeVector{-225, 100, 225},
    G4ThreeVector{0, 100, 225},
    G4ThreeVector{225, 100, 225},
    G4ThreeVector{225, 100, 0},
    G4ThreeVector{225, 100, -225},
    G4ThreeVector{0, 100, -225},
    G4ThreeVector{-225, 100, -225},
    G4ThreeVector{-225, 100, 0}},

    {G4ThreeVector{-325, 100, 0},
    G4ThreeVector{-325, 100, 325},
    G4ThreeVector{0, 100, 325},
    G4ThreeVector{325, 100, 325},
    G4ThreeVector{325, 100, 0},
    G4ThreeVector{325, 100, -325},
    G4ThreeVector{0, 100, -325},
    G4ThreeVector{-325, 100, -325},
    G4ThreeVector{-325, 100, 0}},

    {G4ThreeVector{-325, 0, 0},
    G4ThreeVector{-325, 0, 325},
    G4ThreeVector{0, 0, 325},
    G4ThreeVector{325, 0, 325},
    G4ThreeVector{325, 0, 0},
    G4ThreeVector{325, 0, -325},
    G4ThreeVector{0, 0, -325},
    G4ThreeVector{-325, 0, -325},
    G4ThreeVector{-325, 0, 0}},

    {G4ThreeVector{-325, -100, 0},
    G4ThreeVector{-325, -100, 325},
    G4ThreeVector{0, -100, 325},
    G4ThreeVector{325, -100, 325},
    G4ThreeVector{325, -100, 0},
    G4ThreeVector{325, -100, -325},
    G4ThreeVector{0, -100, -325},
    G4ThreeVector{-325, -100, -325},
    G4ThreeVector{-325, -100, 0}},

    {G4ThreeVector{-225, -100, 0},
    G4ThreeVector{-225, -100, 225},
    G4ThreeVector{0, -100, 225},
    G4ThreeVector{225, -100, 225},
    G4ThreeVector{225, -100, 0},
    G4ThreeVector{225, -100, -225},
    G4ThreeVector{0, -100, -225},
    G4ThreeVector{-225, -100, -225},
    G4ThreeVector{-225, -100, 0}},

    {G4ThreeVector{-125, -100, 0},
    G4ThreeVector{-125, -100, 125},
    G4ThreeVector{0, -100, 125},
    G4ThreeVector{125, -100, 125},
    G4ThreeVector{125, -100, 0},
    G4ThreeVector{125, -100, -125},
    G4ThreeVector{0, -100, -125},
    G4ThreeVector{-125, -100, -125},
    G4ThreeVector{-125, -100, 0}},

    {G4ThreeVector{-125, 0, 0},
    G4ThreeVector{-125, 0, 125},
    G4ThreeVector{0, 0, 125},
    G4ThreeVector{125, 0, 125},
    G4ThreeVector{125, 0, 0},
    G4ThreeVector{125, 0, -125},
    G4ThreeVector{0, 0, -125},
    G4ThreeVector{-125, 0, -125},
    G4ThreeVector{-125, 0, 0}},
  };
  return controlPts;
}


std::vector<std::vector<G4double>> getTorusWeights()
{
  std::vector<std::vector<G4double>> weights = {
  {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
  {0.707, 0.707, 0.707, 0.707, 0.707, 0.707, 0.707, 0.707, 0.707},
  {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
  {0.707, 0.707, 0.707, 0.707, 0.707, 0.707, 0.707, 0.707, 0.707},
  {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
  {0.707, 0.707, 0.707, 0.707, 0.707, 0.707, 0.707, 0.707, 0.707},
  {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
  {0.707, 0.707, 0.707, 0.707, 0.707, 0.707, 0.707, 0.707, 0.707},
  {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
};
  return weights;
}

}// namespace georgeTorus




int main(int argc, char** argv)
{


  // useful strings for output messages

  std::cout<<"-------------------------------------------------------------------------"<<std::endl;

  std::string start_section = "----------------------------";
  std::string end_section = "-------------------------------------------------------------------------\n\n";


  // creating a nurbs torus
  std::vector<std::vector<G4ThreeVector>> torusControlPts = georgeTorusTests::getTorusControlPts();
  std::vector<std::vector<G4double>> torusWeights = georgeTorusTests::getTorusWeights();

  std::vector<G4double> torusKnotsU = {0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1};
  std::vector<G4double> torusKnotsV = {0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1};

  int torusDegreeU = 2;
  int torusDegreeV = 2;


  G4GeorgeNurbs* torusNurbs = new G4GeorgeNurbs("Nurbs Torus",
                                                torusControlPts, torusWeights,
                                                torusKnotsU, torusKnotsV,
                                                torusDegreeU, torusDegreeV);



  // checking the bounding box

  std::cout<<start_section<<" bounding box checks"<<start_section<<std::endl;
  G4Box BoundingBox = torusNurbs->GetBoundingBox();
  std::cout << "Bounding box half-lengths: "
          << BoundingBox.GetXHalfLength() << ", "
          << BoundingBox.GetYHalfLength() << ", "
          << BoundingBox.GetZHalfLength() << std::endl;


  G4ThreeVector P_inside_check = G4ThreeVector{10, 10, -400};
  G4ThreeVector v_DTI_check = G4ThreeVector{0, 0, 1};
  std::cout<<"BoundingBox.Inside("<<P_inside_check<<"): "<<BoundingBox.Inside(P_inside_check)<<std::endl;
  std::cout<<"BoundingBox.DistanceToIn("<<P_inside_check<<","<<v_DTI_check<<"): "<<BoundingBox.DistanceToIn(P_inside_check, v_DTI_check)<<std::endl;

  std::cout<<end_section<<std::endl;


  // check distance to in
  std::cout<<start_section<<" DistanceToIn(p,v)"<<start_section<<std::endl;
  torusNurbs->EnableOptVerbose();


  G4ThreeVector p0_out_to_in = G4ThreeVector(230,50,-400);
  G4ThreeVector v_out_to_in = G4ThreeVector{0, 0, 1};

  double d_to_in = torusNurbs->DistanceToIn(p0_out_to_in, v_out_to_in);

  std::cout<<"d_to_in = "<<d_to_in<<std::endl;

  std::cout<<end_section<<std::endl;














  delete torusNurbs;




  std::cout<<"working on BruteForceIntersection branch"<<std::endl;








  return 0;


}















//  // Detect interactive mode (if no arguments) and define UI session
//  //
//  G4UIExecutive* ui = nullptr;
//  if (argc == 1) {
//    ui = new G4UIExecutive(argc, argv);
//  }
//
//  // Optionally: choose a different Random engine...
//  // G4Random::setTheEngine(new CLHEP::MTwistEngine);
//
//  // use G4SteppingVerboseWithUnits
//  G4int precision = 4;
//  G4SteppingVerbose::UseBestUnit(precision);
//
//  // Construct the default run manager
//  //
//  auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
//
//  // Set mandatory initialization classes
//  //
//  // Detector construction
//  runManager->SetUserInitialization(new DetectorConstruction());
//
//  // Physics list
//  auto physicsList = new QBBC;
//  physicsList->SetVerboseLevel(1);
//  runManager->SetUserInitialization(physicsList);
//
//  // User action initialization
//  runManager->SetUserInitialization(new ActionInitialization());
//
//  // Initialize visualization with the default graphics system
//  auto visManager = new G4VisExecutive(argc, argv);
//  // Constructors can also take optional arguments:
//  // - a graphics system of choice, eg. "OGL"
//  // - and a verbosity argument - see /vis/verbose guidance.
//  // auto visManager = new G4VisExecutive(argc, argv, "OGL", "Quiet");
//  // auto visManager = new G4VisExecutive("Quiet");
//  visManager->Initialize();
//
//  // Get the pointer to the User Interface manager
//  auto UImanager = G4UImanager::GetUIpointer();
//
//  // Process macro or start UI session
//  //
//  if (!ui) {
//    // batch mode
//    G4String command = "/control/execute ";
//    G4String fileName = argv[1];
//    UImanager->ApplyCommand(command + fileName);
//  }
//  else {
//    // interactive mode
//    UImanager->ApplyCommand("/control/execute init_vis.mac");
//    ui->SessionStart();
//    delete ui;
//  }
//
//  // Job termination
//  // Free the store: user actions, physics_list and detector_description are
//  // owned and deleted by the run manager, so they should not be deleted
//  // in the main() program !
//
//  delete visManager;
//  delete runManager;
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
