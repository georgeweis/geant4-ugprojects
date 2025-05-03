//code for cylinder tests


#include "G4GeorgeNurbs.hh"
#include "G4GeorgeSolid.hh"
#include "G4Box.hh"
#include "SavedTorusData.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <iostream>
#include <vector>
#include <string>
//#include <nlopt.hpp>



using namespace B1;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


int main(int argc, char** argv)
{

  // instantiaion of my nurbs cylinder
  G4double r = 1.0;
  G4double w = std::sqrt(2) / 2.0;

  std::vector<std::vector<G4ThreeVector>> controlPts = {
    {
      G4ThreeVector{0.0, -r,  r},
      G4ThreeVector{ r,   -r,  r},
      G4ThreeVector{ r,    0.0, r},
      G4ThreeVector{ r,    r,  r},
      G4ThreeVector( 0.0,  r,  r)
  },
  {
    G4ThreeVector{ 0.0, -r, -r},
    G4ThreeVector{ r,   -r, -r},
    G4ThreeVector{ r,    0.0, -r},
    G4ThreeVector{ r,    r, -r},
    G4ThreeVector{ 0.0,  r, -r}
}};

  std::vector<std::vector<G4double>> weights = {
    {1, w, 1, w, 1},
    {1, w, 1, w, 1}
  };

  std::vector<G4double> knotsU = {0, 0, 1, 1};
  std::vector<G4double> knotsV = {0, 0, 0, 0.5, 0.5, 1, 1, 1};

  int degreeU = 1;
  int degreeV = 2;

  G4GeorgeNurbs* NurbsCylinder = new G4GeorgeNurbs("Nurbs Half Cylinder",
                                                   controlPts,
                                                   weights,
                                                   knotsU,
                                                   knotsV,
                                                   degreeU,
                                                   degreeV);


  // instantiation of a G4 cylinder
  double radius = 1;
  double halfLength = 1;
  auto G4Cylinder = new G4Tubs("CylZ", 0, radius, halfLength, 0., 360.*deg);



  // finding the intersection with a line
  G4ThreeVector p_line = G4ThreeVector{3,0,0};
  G4ThreeVector n_line_unnormalised= G4ThreeVector{-1,0.,0.3};
  G4ThreeVector n_line = n_line_unnormalised.unit();

  NurbsCylinder->EnableOptVerbose();

  double l_nurbs = NurbsCylinder->DistanceToIn(p_line, n_line);
  G4ThreeVector P_nurbs = (p_line + l_nurbs*n_line);
  std::cout<<"l_nurbs: "<<l_nurbs<<std::endl;
  std::cout<<"P_nurbs (at intersection): "<<P_nurbs<<std::endl;

  double  l_ana = G4Cylinder->DistanceToIn(p_line, n_line);
  G4ThreeVector P_ana = (p_line + l_ana*n_line);
  std::cout<<"l_ana: "<<l_nurbs<<std::endl;
  std::cout<<"P_ana (at intersection): "<<P_ana<<std::endl;

  double nurbs_ana_residual = (P_nurbs-P_ana).mag();
  std::cout<<"nurbs_ana_residual: "<<nurbs_ana_residual<< " l_nurbs-l_ana: "<<l_nurbs-l_ana<<std::endl;

  /* Intersection accuracy data
  // creating random lines to check
  double n_y_extent = 0.3;
  double n_z_extent = 0.4;
  std::vector<double> nurbs_ana_residuals;
  double d_to_in_nurbs;
  double d_to_in_ana;
  G4double y_compt_rad;
  G4double z_compt_rad;



  for (int i = 0; i < 4000; ++i) {
    y_compt_rad = -n_y_extent + 2 * n_y_extent * G4UniformRand();
    z_compt_rad = -n_z_extent + 2 * n_z_extent * G4UniformRand();

    G4ThreeVector n_rand(-1.0, y_compt_rad, z_compt_rad);
    n_rand = n_rand.unit();
    d_to_in_nurbs = NurbsCylinder->DistanceToIn(p_line, n_rand);
    d_to_in_ana = G4Cylinder->DistanceToIn(p_line, n_rand);


    nurbs_ana_residuals.push_back(d_to_in_ana - d_to_in_nurbs);
    //std::cout<<d_to_in_nurbs-d_to_in_ana<<std::endl;
  }

  std::cout<<"nurbs_ana_residuals = np.array([";
  for(int i = 0; i < nurbs_ana_residuals.size(); ++i) {
    std::cout<<nurbs_ana_residuals[i]<<",";
  }
  std::cout<<" ])"<<std::endl;

   */













  delete NurbsCylinder;
  delete G4Cylinder;


}
