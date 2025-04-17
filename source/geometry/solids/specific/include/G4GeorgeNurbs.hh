// nurbs class header



#ifndef G4GEORGENURBS_HH
#define G4GEORGENURBS_HH

#include <iostream>
#include <vector>

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"


class G4GeorgeNurbs : public G4VSolid
{
  private:
  std::vector<std::vector<G4ThreeVector>> controlPts;
  std::vector<std::vector<G4double>> weights;
  std::vector<G4double> knotVectorU;
  std::vector<G4double> knotVectorV;
  G4int degreeU;
  G4int degreeV;

  static double LARGE_NUMBER;
  static G4ThreeVector LARGE_THREE_VECTOR;


  public:
  G4GeorgeNurbs(const G4String& name,
                const std::vector<std::vector<G4ThreeVector>>& controlPts_in,
                const std::vector<std::vector<G4double>>& weights_in,
                const std::vector<G4double>& knotVectorU_in,
                const std::vector<G4double>& knotVectorV_in,
                const G4int& degreeU_in,
                const G4int& degreeV_in);

  ~G4GeorgeNurbs() override;

  void PrintVariables() const;


  // Functions required to define a NURBS surface ====================================================

  G4double BasisFunction(G4int i, G4int k, G4double t, const std::vector<G4double>& knotVector) const;
  G4ThreeVector SurfacePoint(G4double u, G4double v) const;



  // Functions for optimisation techniques ============================================================
  G4TwoVector ClosestKnot(const G4ThreeVector& point) const;
  // iterates through the unique knots to find the closest knot to a given point.
  // Used as initial guess for optimisation techniques

  G4ThreeVector ClosestPoint(const G4ThreeVector& point) const;
  // Finds the closest point on the NURBS surface to a given point

  G4ThreeVector LineIntersection(const G4ThreeVector& P0,
                                 const G4ThreeVector& direction,
                                 double lambda_bound) const;
  // finds the point of intersection between a line and a nurbs surface








  // Static functions and structs used in optimisation ======================================

  private:
  // #ifdef GEANT4_USE_NLOPT

  // optimisation for ClosestPoint -------------------------------------
  struct OptimizationContext
  {
    const G4GeorgeNurbs* nurbs;
    G4ThreeVector target_point;
  }; // struct to pass context for closest point optimisation

  static double ResidualToSurfacePoint(const std::vector<double>& uv,
                                       std::vector<double>& grad,
                                       void* data);
  // calculates the magnitude of the residual three-vector between a surface point (defined by uv)
  // and a target point defined within OptimizationContext struct. Static function used because
  // NLopt cannot deal with member functions and arguments grad and data are required for nlopt.




  // optimisation for line intersection -------------------------------------
  struct LineIntersectionContext {
    const G4GeorgeNurbs* nurbs;
    G4ThreeVector P0;
    G4ThreeVector direction; // should be unit vector
    double max_lambda;
  };

  static double ResidualLineDistance(const std::vector<double>& uvl,
                                     std::vector<double>& grad,
                                     void* data);
  // calculates the distance between a point on a line (defined by P0, direction and lambda)
  // and a point on the surface (defined by the prarameters u and v).
  // The first argument is a vector [u,v,l] and defines the parameters to optimise in LineIntersecion




  // #endif














  public:
  // Overidden Base class functions ======================================================================
  EInside Inside(const G4ThreeVector& p) const override;
  G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;


  G4double DistanceToIn(const G4ThreeVector& p) const override;



  G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const override;
  G4double DistanceToOut(const G4ThreeVector& p) const override;
  G4double DistanceToOut(const G4ThreeVector& p,
                         const G4ThreeVector& v,
                         const G4bool calcNorm = false,
                               G4bool* validNorm = nullptr,
                               G4ThreeVector* n = nullptr) const override;

  G4double GetCubicVolume() override;
  G4double GetSurfaceArea() override;

  G4ThreeVector GetPointOnSurface() const override;

  G4VisExtent GetExtent() const override;
  G4bool CalculateExtent( const EAxis pAxis,
                          const G4VoxelLimits& pVoxelLimit,
                          const G4AffineTransform& pTransform,
                                G4double& pMin, G4double& pMax ) const override;

  void DescribeYourselfTo ( G4VGraphicsScene& scene ) const override;
  std::ostream& StreamInfo(std::ostream& os) const override;
  G4GeometryType GetEntityType() const override;



  void ValidateKnotVectors() const;
  // ensure that knot vectors have length nb_controlPts + degree + 1 in both directions
  // called in the constructor and throws an error if invalid









};



#endif // G4GEORGENURBS_HH
