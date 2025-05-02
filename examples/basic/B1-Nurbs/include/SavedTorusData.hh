//
// Created by George Weis on 29/04/2025.
//




#ifndef GEORGETORUSNURBS_H
#define GEORGETORUSNURBS_H


#include "G4GeorgeNurbs.hh"

namespace B1 {



class SavedTorusData{
public:
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
      {0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107},
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
      {0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107},
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
      {0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107},
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
      {0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107, 0.707107},
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
    };
    return weights;
  }


  std::vector<std::vector<G4ThreeVector>> getOuterTokamakControlPts(){
    std::vector<std::vector<G4ThreeVector>> controlPts = getTorusControlPts();
    for (int i = 0; i < 9; ++i)
    {
      double scale_factor;      // Scale factor
      double end_ctrl_pt_scale_factor = 1.33;
      double adjacent_ctrl_pt_scale_factor = 1.33;



      // Iterate over 4 rows (0, 1, 7, 8)
      for (int j = 0; j < 4; ++j)
      {
        int u;                   // Row index


        if (j == 0)
        {
          u = 0;                // First row
          scale_factor = end_ctrl_pt_scale_factor;   // End knot scaling
        }
        else if (j == 1)
        {
          u = 1;                // Second row
          scale_factor = adjacent_ctrl_pt_scale_factor;  // Interior scaling
        }
        else if (j == 2)
        {
          u = 7;                // Second to last row
          scale_factor = adjacent_ctrl_pt_scale_factor;
        }
        else if (j == 3)
        {
          u = 8;                // Last row
          scale_factor = end_ctrl_pt_scale_factor;
        }

        // Scale x and z components of control point
        controlPts[u][i].setX(controlPts[u][i].x() * scale_factor);
        controlPts[u][i].setZ(controlPts[u][i].z() * scale_factor);
      }
    }
    return controlPts;
    }


  std::vector<std::vector<G4double>> getTokamakWeights(){

    std::vector<std::vector<G4double>> weights = getTorusWeights();

    double weight_scale_factor = 4.0;
    for (int i = 0; i < 9; ++i)
    {
      weights[1][i] = weights[1][i]*weight_scale_factor;
      weights[7][i] = weights[7][i]*weight_scale_factor;
    }
    return weights;
  }



  std::vector<std::vector<G4ThreeVector>> getInnerTokamakControlPts(){

    std::vector<std::vector<G4ThreeVector>> controlPts = getOuterTokamakControlPts();
    std::vector<std::vector<G4ThreeVector>> inner_tokamak_cpts = controlPts;
    double inner_scale_factor = 0.9;

    for (int i = 0; i < 9; ++i)
    {
      // Extract the entire column of control points at index i
      std::vector<G4ThreeVector> c_pts;
      for (int u = 0; u < controlPts.size(); ++u)
      {
        c_pts.push_back(controlPts[u][i]);
      }

      // Compute average radial position (between points 0 and 4)
      double avj_rad = (c_pts[0].mag() + c_pts[4].mag()) / 2.0;

      // Compute average position vector (x and z components)
      double avj_x = (c_pts[0].x() + c_pts[4].x()) / 2.0;
      double avj_y = 0.0; // explicitly 0
      double avj_z = (c_pts[0].z() + c_pts[4].z()) / 2.0;
      G4ThreeVector avj_pos(avj_x, avj_y, avj_z);

      // Loop over all 9 control points along that column
      for (int j = 0; j < 9; ++j)
      {
        // Shift control point relative to average position
        G4ThreeVector transformed_ctrlpt = c_pts[j] - avj_pos;

        // Scale relative position
        G4ThreeVector scaled_transformed_ctrlpt = transformed_ctrlpt * inner_scale_factor;

        // Shift back
        G4ThreeVector scaled_ctrlpt = scaled_transformed_ctrlpt + avj_pos;

        // Store into new control point array
        inner_tokamak_cpts[j][i] = scaled_ctrlpt;
      }
    }
    return inner_tokamak_cpts;
  }

}; // class saved data

} //namespace





#endif //GEORGETORUSNURBS_H
