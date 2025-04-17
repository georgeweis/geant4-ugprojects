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

#include <iostream>
#include <vector>
#include <string>
//#include <nlopt.hpp>



//using namespace B1;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  //editing here means I only need to rebuild the example

  //sphere class
  /*
  G4GeorgeSolid* georgeSolid = new G4GeorgeSolid("test", 5.0);
  std::cout<<"georgeSolid->print_info() output: "<<std::endl;
  georgeSolid->print_info();
  delete georgeSolid;
  */

  //nurbs class

  G4double r = 1;
  std::vector<std::vector<G4ThreeVector>> controlPts = {
    {G4ThreeVector{0, 0, r}, G4ThreeVector{0, r, r}, G4ThreeVector{0, r, 0}},
    {G4ThreeVector{r, 0, r}, G4ThreeVector{r, r, r}, G4ThreeVector{r, r, 0}}
    };
  std::vector<std::vector<G4double>> weights = {
    {1, 0.7071, 1},
    {1, 0.7071, 1}
    };

  std::vector<G4double> knotsU = {0,0,1,1};
  std::vector<G4double> knotsV = {0,0,0,1,1,1};

  int degreeU = 1;
  int degreeV = 2;



  G4GeorgeNurbs* georgeNurbs = new G4GeorgeNurbs("test", controlPts, weights, knotsU, knotsV, degreeU, degreeV);


  std::cout<<"\ngeorgeNurbs->PrintVariables() output: "<<std::endl;
  georgeNurbs->PrintVariables();



  G4ThreeVector surfacePt = georgeNurbs->SurfacePoint(0.3,0.3);
  std::cout<<"\ngeorgeNurbs->SurfacePoint(0.3,0.3) output: "<<surfacePt<<std::endl;

  // ✅✅ same output as Python class --> Nurbs_Surface10_arc.py, TEST_NB = 0


  std::cout<<"-------------------------------------------------------------------------"<<std::endl;

  std::string start_section = "----------------------------";
  std::string end_section = "-------------------------------------------------------------------------\n\n";





  // closes knot test
  std::cout<<start_section<<"closest knot test"<<start_section <<"\n"<<std::endl;

  G4ThreeVector const_point = G4ThreeVector(0,0,0.2);

  G4TwoVector uv_closest_knot = georgeNurbs->ClosestKnot(const_point);
  G4ThreeVector p_closest_knot = georgeNurbs->SurfacePoint(uv_closest_knot[0],uv_closest_knot[1]);


  std::cout<<"const_point: "<<const_point<<std::endl;
  std::cout<<"uv_closest_knot: "<<uv_closest_knot<<std::endl;
  std::cout<<"p_closest_knot: "<<p_closest_knot<<std::endl;

  std::cout<<end_section;

  // ✅✅ same output as Python class --> Nurbs_Surface10_arc.py, TEST_NB = 1




  // closes surface point test
  std::cout<<start_section<<"closest surface point test"<<start_section <<"\n"<<std::endl;

  G4ThreeVector const_point1 = G4ThreeVector(0.5,0.2,0.3);
  G4ThreeVector P_surface_opt = georgeNurbs->ClosestPoint(const_point1);
  std::cout<<"const_point: "<<const_point1<<std::endl;
  std::cout<<"P_surface_opt: "<<P_surface_opt<<std::endl;

  std::cout<<end_section;
  // ✅✅ same output as Python class --> Nurbs_Surface10_arc.py, TEST_NB = 2







  // line intersection test
  std::cout<<start_section<<"line intersection test"<<start_section <<"\n"<<std::endl;

  int SUB_TEST_NB = 2;
  G4ThreeVector p0;
  G4ThreeVector nline;
  double lambda_max;

  if (SUB_TEST_NB == 0)
  {
  	// Line which intersects the surface
  	p0 = G4ThreeVector(0.5, 1.0, 1.0);
  	nline = G4ThreeVector(0.0, -1.0, -1.0).unit();
  	lambda_max = 3.0;
  }// ✅✅ same output as Python class

  if (SUB_TEST_NB == 1)
  {
    // Line is too short to reach the surface
  	p0 = G4ThreeVector(0.5, 0.6, 0.6);
  	nline = G4ThreeVector(0.0, -1.0, -1.0).unit();
  	lambda_max = 0.1;
  }// ✅✅ same output as Python class

  if (SUB_TEST_NB == 2)
  {
  	// Infinite line misses the surface completely
  	p0 = G4ThreeVector(0.5, 1.0, 0.6);
  	nline = G4ThreeVector(10.0, -1.0, 1.0).unit();
  	lambda_max = 3.0;
  }// ✅✅ same output as Python class


  G4ThreeVector line_intersec = georgeNurbs->LineIntersection(p0, nline, lambda_max);


  std::cout<<"p0 = "<<p0<<std::endl;
  std::cout<<"nline = "<<nline<<std::endl;
  std::cout<<"lambda_max = "<<lambda_max<<std::endl;

  std::cout<<"line_intersec = "<<line_intersec<<std::endl;
  std::cout<<end_section;

  // ✅✅ same output as Python class --> Nurbs_Surface10_arc.py, TEST_NB = 5



  delete georgeNurbs;

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
