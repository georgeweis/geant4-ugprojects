//
// The code below is the exact main that i used to perform all the nurbs surface
// tests once the solid was fully implemented. Ive saved it here in case i need
// to revisit any of these tests.


















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
  //editing here means I only need to rebuild the example

  //sphere class
  /*
  G4GeorgeSolid* georgeSolid = new G4GeorgeSolid("test", 5.0);
  std::cout<<"georgeSolid->print_info() output: "<<std::endl;
  georgeSolid->print_info();
  delete georgeSolid;
  */

  //nurbs class


  // useful strings for output messages

  std::string start_section = "-------------------------";
  std::string end_section = "----------------------------------------------------------------------\n\n";








  if(false) // wrapped arc nurbs surface to separate from torus tests
  {

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

   std::cout<<"G4ThreeVector const_point = G4ThreeVector(0,0,0.2);"<<std::endl;


  auto [uv_closest_knot, R_closest_knot] = georgeNurbs->ClosestKnot(const_point);

  std::cout<<"std::vector<double> uv_closest_knot = georgeNurbs->ClosestKnot(const_point);"<<std::endl;

  G4ThreeVector p_closest_knot = georgeNurbs->SurfacePoint(uv_closest_knot[0],uv_closest_knot[1]);

  std::cout<<"G4ThreeVector p_closest_knot = georgeNurbs->SurfacePoint(uv_closest_knot[0],uv_closest_knot[1]);"<<std::endl;



  std::cout<<"const_point: "<<const_point<<std::endl;
  std::cout<<"uv_closest_knot: "<<uv_closest_knot.size()<<std::endl;//"("<<uv_closest_knot[0]<<","<<uv_closest_knot[1]<<")"<<std::endl;
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


//  G4ThreeVector line_intersec = georgeNurbs->LineIntersection(p0, nline, lambda_max);
//
//
//  std::cout<<"p0 = "<<p0<<std::endl;
//  std::cout<<"nline = "<<nline<<std::endl;
//  std::cout<<"lambda_max = "<<lambda_max<<std::endl;
//
//  std::cout<<"line_intersec = "<<line_intersec<<std::endl;
//  std::cout<<end_section;

  // ✅✅ same output as Python class --> Nurbs_Surface10_arc.py, TEST_NB = 5


  //surface derivative + normal test
  std::cout << start_section << "surface derivative + normal test" << start_section << "\n" << std::endl;

  G4double u = 0.3;
  G4double v = 0.7;


  G4ThreeVector P_surface = georgeNurbs->SurfacePoint(u, v); // some point on surface P(u,v)

  std::vector<G4ThreeVector> tangents = georgeNurbs->SurfaceDerivatives(u, v); // tangents at that point
  G4ThreeVector tangent_u = tangents[0];
  G4ThreeVector tangent_v = tangents[1];

  // cross product of the tanglent vectors is the normal
  G4ThreeVector normal = tangent_u.cross(tangent_v);

  // using the SurfaceNormal(u,v) function
  G4ThreeVector normal_uv_func = georgeNurbs->SurfaceNormal(u, v);//✅✅ same output as Python class TEST_NB = 6

  // using overridden SurfaceNormal(p) function
  G4ThreeVector normal_p_func = georgeNurbs->SurfaceNormal(P_surface); //✅✅ same output as SurfaceNormal(u, v)



  std::cout << "u: " << u << ", v: " << v << std::endl;
  std::cout << "P_surface: " << P_surface << std::endl;
  std::cout << "tangent_u: " << tangent_u << std::endl;
  std::cout << "tangent_v: " << tangent_v << std::endl;
  std::cout << "normal: " << normal << std::endl;
  std::cout << "normal using SurfaceNormal(u, v): " << normal_uv_func << std::endl;
  std::cout << "normal using SurfaceNormal(p): " << normal_p_func << std::endl;

  std::cout << end_section;
  // ✅✅ same output as Python class --> Nurbs_Surface10_arc.py, TEST_NB = 6

  delete georgeNurbs;

  }// end of arc tests




  // tests of torus shape =======================================================================================

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


//  torusNurbs->PrintVariables();


  //inside/outsidee torus test
  std::cout << start_section << "inside/outsise torus test" << start_section << "\n" << std::endl;

  G4ThreeVector inside_pt = G4ThreeVector(20, 20, 200);
  G4ThreeVector surface_pt = torusNurbs->ClosestPoint(inside_pt);
  G4ThreeVector outside_pt = G4ThreeVector(10, 10, 10);

  std::cout<<"inside_pt = "<<inside_pt<<std::endl;
  std::cout<<"surface_pt = "<<surface_pt<<std::endl;
  std::cout<<"outside_pt = "<<outside_pt<<std::endl;

  std::cout<<"Output from torusNurbs->Inside(inside_pt): "<<torusNurbs->Inside(inside_pt)<<std::endl;
  std::cout<<"Output from torusNurbs->Inside(surface_pt): "<<torusNurbs->Inside(surface_pt)<<std::endl;
  std::cout<<"Output from torusNurbs->Inside(outside_pt): "<<torusNurbs->Inside(outside_pt)<<std::endl;


  std::cout << end_section;
  // ✅✅ seems to work fine --> Nurbs_Surface9_torus.py, TEST_NB = 3


  //inside/outsidee torus test
  std::cout << start_section << "DistanceToIn(p) / DistanceToOut(p)" << start_section << "\n" << std::endl;

  std::cout << start_section<<std::endl;
  std::cout<<"DistanceToIn(inside_pt): "<<torusNurbs->DistanceToIn(inside_pt)<<std::endl;
  std::cout<<"DistanceToIn(surface_pt): "<<torusNurbs->DistanceToIn(surface_pt)<<std::endl;
  std::cout<<"DistanceToIn(outside_pt): "<<torusNurbs->DistanceToIn(outside_pt)<<std::endl;

  std::cout<<"\n";
  std::cout<<"DistanceToOut(inside_pt): "<<torusNurbs->DistanceToOut(inside_pt)<<std::endl;
  std::cout<<"DistanceToOut(surface_pt): "<<torusNurbs->DistanceToOut(surface_pt)<<std::endl;
  std::cout<<"DistanceToOut(outside_pt): "<<torusNurbs->DistanceToOut(outside_pt)<<std::endl;
  std::cout << start_section<<std::endl;


  std::cout << end_section;
  // ✅✅ seems to work fine

  //bounding limits checks
  std::cout << start_section << "bounding limits checks" << start_section << "\n" << std::endl;
  auto bounds = torusNurbs->GetBounds();
  std::cout <<"bmin: " << bounds[0]<<std::endl;
  std::cout <<"bmax: " << bounds[1]<<std::endl;
  std::cout <<"maxExtent: " << torusNurbs->GetMaxExtent()<<std::endl;
  std::cout << end_section;



//  std::cout << "Without verbose" << start_section << "\n" << std::endl;
//  torusNurbs->DisableOptVerbose();
//  G4ThreeVector P0 = G4ThreeVector(10,0,250);
//  G4ThreeVector direction = G4ThreeVector(0, 0, 1);
//  double line_length = 1000;
//
//  auto [uvl_opt, minimised_residual] = torusNurbs->LineIntersectionParams(P0, direction.unit(), line_length, 1);

  std::cout << end_section;



  //final test of important functions
  std::cout << start_section << "final test of important functions" << start_section << "\n" << std::endl;
  /* After getting it to show up in the gui with errors, need to do a test of
    - Inside(p)
    - DistanceToIn(p)
    - DistanceToOut(p)
    - DistanceToIn(p,v)
    - DistanceToOut(p,v)
  */

  // problem:
  // Track stuck, not moving for 25 steps.
  //  Current  phys volume: 'World'
  //   - at position : (0,0,125.0012635833764)
  //     in direction: (0,0,1)
  //    (local position: (0,0,125.0012635833764))
  //    (local direction: (0,0,1)).

  torusNurbs->EnableOptVerbose();



  //Test 1
  std::cout<<start_section<<"Test 1"<<start_section<<std::endl;

  // will just rewite these
  G4ThreeVector P0_test;
  G4ThreeVector direction_test;



  // 1) origin inside?
  P0_test = G4ThreeVector(0,0,0);
  std::cout<<"Check 1"<<start_section<<start_section <<"\n"
           <<"p = "<<P0_test<<"\n"
           <<"Inside(p): "<<torusNurbs->Inside(P0_test)<<"\n"
           <<start_section<<std::endl;

  // 2) distance to out from intersection point
  P0_test = G4ThreeVector(0,0,125.0012635833764);
  direction_test = G4ThreeVector(0,0,1);
  G4double DTO_1 = torusNurbs->DistanceToOut(P0_test, direction_test.unit());

  std::cout<<"Check 2"<<start_section <<start_section<<"\n"
           <<"p = "<<P0_test<<"\n"
           <<"v = "<<direction_test<<"\n"
           <<"DistanceToOut(p,v): "<<DTO_1<<"\n"
           <<start_section<<std::endl;

  std::cout<<start_section<<"\n"<<std::endl;


  //Test 2
  std::cout<<start_section<<"Test 2"<<start_section<<std::endl;
  P0_test = G4ThreeVector(0,0,-400);
  direction_test = G4ThreeVector(0,0,1);

  //check 1 Inside
  P0_test = G4ThreeVector(0,0,0);
  std::cout<<"Check 1"<<start_section<<start_section <<"\n"
           <<"p = "<<P0_test<<"\n"
           <<"Inside(p): "<<torusNurbs->Inside(P0_test)<<"\n"
           <<start_section<<std::endl;


  // check 2 distance to in
  P0_test = G4ThreeVector(0,0,-400);
  G4double DTI_1 = torusNurbs->DistanceToIn(P0_test, direction_test.unit());

  std::cout<<"Check 2"<<start_section <<start_section<<"\n"
           <<"p = "<<P0_test<<"\n"
           <<"v = "<<direction_test<<"\n"
           <<"DistanceToIn(p,v): "<<DTO_1<<"\n"
           <<start_section<<std::endl;


  // check 3 closest knot
  P0_test = G4ThreeVector(0,0,-400);
  std::cout<<"Check 3 - closest knot"<<start_section <<start_section<<"\n";
  auto [closest_knot, R_knot] = torusNurbs->ClosestKnot(P0_test);

  std::cout<<"closest_knot = "<<"("<<closest_knot[0] <<","<<closest_knot[1]<<")"<<"\n"
           <<"R_knot = "<<R_knot<<"\n"
           <<start_section<<std::endl;

  std::cout<<start_section<<"\n"<<std::endl;

  std::cout << end_section;





  //line intersection
 std::cout << start_section << "line intersection" << start_section << "\n" << std::endl;

  std::cout << "With verbose" << start_section << "\n" << std::endl;
  torusNurbs->EnableOptVerbose();
  G4ThreeVector P0_verb = G4ThreeVector(20,20,-300);
  G4ThreeVector direction_verb = G4ThreeVector(0,-1, 1);
  double line_length_verb = 1500;


  auto [uvl_opt_verb, minimised_residual_verb] = torusNurbs->LineIntersectionParams(P0_verb, direction_verb.unit(), line_length_verb, 1);


  // linking intersection tests

//  torusNurbs->DisableOptVerbose();
  if(false)
  {
  std::cout << start_section << "// linking intersection tests" << start_section << "\n" << std::endl;

  std::vector<G4ThreeVector> linked_intersections;

  // P0_cont1 is below torus
  G4ThreeVector P0_cont1 = G4ThreeVector(0,0,-400);
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

  G4double distance_cont3 = torusNurbs->DistanceToIn(P0_cont3, direction_cont.unit());
  std::cout<<start_section<<"\n distance_cont3: "<<distance_cont3<<"\n"<<std::endl;


  // P0_cont4 is top of innter radius
  G4ThreeVector P0_cont4 = P0_cont3 + distance_cont3*direction_cont;
  linked_intersections.push_back(P0_cont4);

  G4double distance_cont4 = torusNurbs->DistanceToOut(P0_cont4, direction_cont.unit());
  std::cout<<start_section<<"\n distance_cont4: "<<distance_cont3<<"\n"<<std::endl;

  // P0_cont5 is at top of outer radius
  G4ThreeVector P0_cont5 = P0_cont4 + distance_cont4*direction_cont;
  linked_intersections.push_back(P0_cont5);

  G4double distance_cont5 = torusNurbs->DistanceToIn(P0_cont5, direction_cont.unit()); // should be large number


  std::cout<<"linked_intersection = np.array([";
  for(int i = 0; i < linked_intersections.size(); ++i) {
    std::cout<<"["<<linked_intersections[i].x()<<", "<<linked_intersections[i].y()<<", "<<linked_intersections[i].z()<<"],";
  }
  std::cout<<" ])"<<std::endl;
  }



  std::cout<<"working on BruteForceIntersection branch"<<std::endl;











  delete torusNurbs;













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
