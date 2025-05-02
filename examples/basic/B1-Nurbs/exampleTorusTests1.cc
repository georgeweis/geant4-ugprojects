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
#include "SavedTorusData.hh"

#include <iostream>
#include <vector>
#include <string>
//#include <nlopt.hpp>



using namespace B1;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


int main(int argc, char** argv)
{


  // useful strings for output messages

  std::cout<<"-------------------------------------------------------------------------"<<std::endl;

  std::string start_section = "----------------------------";
  std::string end_section = "-------------------------------------------------------------------------\n\n";


  // creating a nurbs torus
  SavedTorusData torusData;
  std::vector<std::vector<G4ThreeVector>> torusControlPts = torusData.getTorusControlPts();
  std::vector<std::vector<G4double>> torusWeights = torusData.getTorusWeights();

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


  G4ThreeVector P_inside_check = G4ThreeVector{-60,0,-400};
  G4ThreeVector v_DTI_check = G4ThreeVector{0, 0, 1};
  std::cout<<"BoundingBox.Inside("<<P_inside_check<<"): "<<BoundingBox.Inside(P_inside_check)<<std::endl;
  std::cout<<"BoundingBox.DistanceToOut("<<P_inside_check<<","<<v_DTI_check<<"): "<<BoundingBox.DistanceToOut(P_inside_check, v_DTI_check)<<std::endl;

  std::cout<<end_section<<std::endl;















  // check discrete line search algorithm
  std::cout<<start_section<<" check discrete line search algorithm "<<start_section<<std::endl;

  //line to test
  G4ThreeVector P0_start = G4ThreeVector(-60,0,-400);
  G4ThreeVector direction = G4ThreeVector{0, 0.3, 1};
  double line_length = 800;

  // step size
  int max_nb_steps = 20;
  double step_size = line_length/max_nb_steps;

  //declaring variables
  std::vector<double> uvl_opt;
  double residual_opt = 1e-1*10;

  // vectors to be filled for indside/ outside points
  std::vector<G4ThreeVector> inside_pts;
  std::vector<G4ThreeVector> outside_pts;
  std::vector<G4ThreeVector> surface_pts;


  //step process
  G4ThreeVector p_previous = P0_start; // assumed to be inside already
  G4ThreeVector p_current;
  int nb_steps_taken = 0;
  bool intersection_found = false;

  for (int i = 0; i < max_nb_steps; i++)
  {
    p_current = p_previous + step_size*direction;
    nb_steps_taken++;

    if (torusNurbs->Inside(p_current) == kInside){inside_pts.push_back(p_current);}
    if (torusNurbs->Inside(p_current) == kOutside){outside_pts.push_back(p_current);}
    if (torusNurbs->Inside(p_current) == kSurface){surface_pts.push_back(p_current);}
    p_previous = p_current; // update p_current for next iteration
  }


  // printing results for python
  std::cout<<"p0 = np.array(["<<P0_start.x()<<","<<P0_start.y()<<","<<P0_start.z()<<"], dtype = float)"<<"\n"
           <<"n_line = np.array(["<<direction.x()<<","<<direction.y()<<","<<direction.z()<<"], dtype = float)"<<"\n"
           <<"line_length = "<<line_length<<std::endl;


  std::cout<<"inside_pts = np.array([";
  for(int i = 0; i < inside_pts.size(); ++i) {
    std::cout<<"["<<inside_pts[i].x()<<", "<<inside_pts[i].y()<<", "<<inside_pts[i].z()<<"],";
  }
  std::cout<<" ])"<<std::endl;

  std::cout<<"outside_pts = np.array([";
  for(int i = 0; i < outside_pts.size(); ++i) {
    std::cout<<"["<<outside_pts[i].x()<<", "<<outside_pts[i].y()<<", "<<outside_pts[i].z()<<"],";
  }
  std::cout<<" ])"<<std::endl;

  std::cout<<"surface_pts = np.array([";
  for(int i = 0; i < surface_pts.size(); ++i) {
    std::cout<<"["<<surface_pts[i].x()<<", "<<surface_pts[i].y()<<", "<<surface_pts[i].z()<<"],";
  }
  std::cout<<" ])"<<std::endl;

  std::cout<<end_section<<std::endl;

  // check distance to in
  std::cout<<start_section<<" DistanceToIn(p,v)"<<start_section<<std::endl;
  torusNurbs->EnableOptVerbose();


  G4ThreeVector p0_out_to_in = G4ThreeVector(230,50,-400);
  G4ThreeVector v_out_to_in = G4ThreeVector{0, 0, 1};

  double d_to_in = torusNurbs->DistanceToIn(p0_out_to_in, v_out_to_in);

  std::cout<<"d_to_in = "<<d_to_in<<std::endl;

  std::cout<<end_section<<std::endl;








  // linked intersection test ======================================

  std::cout << start_section << "// linking intersection tests" << start_section << "\n" << std::endl;

  std::vector<G4ThreeVector> linked_intersections;

  // P0_cont1 is below torus
  G4ThreeVector P0_cont1 = G4ThreeVector(10,-10,0);
  G4ThreeVector direction_cont = G4ThreeVector(0,0, 1);
  direction_cont = direction_cont.unit();
  linked_intersections.push_back(P0_cont1);


  G4double distance_cont1 = torusNurbs->DistanceToIn(P0_cont1, direction_cont.unit());
  std::cout<<start_section<<"\n distance_cont1: "<<distance_cont1<<"\n"<<std::endl; //works



  //P0_cont2 is bottom of torus
  G4ThreeVector P0_cont2 = P0_cont1 + distance_cont1*direction_cont;
  linked_intersections.push_back(P0_cont2);

  G4double distance_cont2 = torusNurbs->DistanceToOut(P0_cont2, direction_cont.unit());
  std::cout<<start_section<<"\n distance_cont2: "<<distance_cont2<<"\n"<<std::endl;
//
//
  // P0_cont3 is bottom of inner loop
  G4ThreeVector P0_cont3 = P0_cont2 + distance_cont2*direction_cont;
  linked_intersections.push_back(P0_cont3);

//  G4double distance_cont3 = torusNurbs->DistanceToIn(P0_cont3, direction_cont.unit());
//  std::cout<<start_section<<"\n distance_cont3: "<<distance_cont3<<"\n"<<std::endl;
//
//
//  // P0_cont4 is top of innter radius
//  G4ThreeVector P0_cont4 = P0_cont3 + distance_cont3*direction_cont;
//  linked_intersections.push_back(P0_cont4);
//
//  G4double distance_cont4 = torusNurbs->DistanceToOut(P0_cont4, direction_cont.unit());
//  std::cout<<start_section<<"\n distance_cont4: "<<distance_cont3<<"\n"<<std::endl;
//
//  // P0_cont5 is at top of outer radius
//  G4ThreeVector P0_cont5 = P0_cont4 + distance_cont4*direction_cont;
//  linked_intersections.push_back(P0_cont5);
//
//  G4double distance_cont5 = torusNurbs->DistanceToIn(P0_cont5, direction_cont.unit()); // should be large number


  std::cout<<"linked_intersection = np.array([";
  for(int i = 0; i < linked_intersections.size(); ++i) {
    std::cout<<"["<<linked_intersections[i].x()<<", "<<linked_intersections[i].y()<<", "<<linked_intersections[i].z()<<"],";
  }
  std::cout<<" ])"<<std::endl;

















  delete torusNurbs;




  std::cout<<"working on BruteForceIntersection branch"<<std::endl;








//  return 0;


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
