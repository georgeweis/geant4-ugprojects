// nurbs class header



#ifndef G4GEORGENURBS_HH
#define G4GEORGENURBS_HH

#include <iostream>
#include <vector>
#include <tuple>

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"
#include "G4Box.hh"



class G4GeorgeNurbs : public G4VSolid
{
  private:

  // variables required for a nurbs surface
  std::vector<std::vector<G4ThreeVector>> controlPts;
  std::vector<std::vector<G4double>> weights;
  std::vector<G4double> knotVectorU;
  std::vector<G4double> knotVectorV;
  G4int degreeU;
  G4int degreeV;

  // bounding variables
  G4ThreeVector bminCached, bmaxCached;
  G4double maxExtent;
  G4bool boundsCached;

  G4Box boundingBox;
  G4ThreeVector boundingBoxCentre;


  // for output on optimisation process on LineIntersection
  bool optVerbose;

  // static variables
  static double LARGE_NUMBER;
  static G4ThreeVector LARGE_THREE_VECTOR;
  static double CONVERGENCE_TOLERANCE;
  static double SURFACE_TOLERANCE;

  // tracking nb function calls and other data
  static bool trackFunctions;

  static int nbInsideCalls; //
  static int nbSurfaceNormalCalls;


  static int nbDTIpCalls; //
  static int nbDTIpvCalls; //
  static int nbDTIdiscreteCalls; //
  static std::vector<int> nbStepsForDTI; //
  static int nbEstimationsForDTI;

  static int nbDTOpCalls; //
  static int nbDTOpvCalls; //
  static int nbDTOdiscreteCalls; //
  static std::vector<int> nbStepsForDTO; //
  static int nbEstimationsForDTO;


  static int nbClosestPointCalls; //
  static std::vector<int> nbFuncEvalsInClosestPoint; //

  static int nbLineIntersecOptCalls; //
  static std::vector<int> nbFuncEvalsInLineIntersec; //














  public:
  G4GeorgeNurbs(const G4String& name,
                const std::vector<std::vector<G4ThreeVector>>& controlPts_in,
                const std::vector<std::vector<G4double>>& weights_in,
                const std::vector<G4double>& knotVectorU_in,
                const std::vector<G4double>& knotVectorV_in,
                const G4int& degreeU_in,
                const G4int& degreeV_in);

  virtual ~G4GeorgeNurbs();

  void PrintVariables() const;


  // Functions required to define a NURBS surface ====================================================

  G4double BasisFunction(G4int i, G4int k, G4double t, const std::vector<G4double>& knotVector) const;
  // recursive basis function

  G4ThreeVector SurfacePoint(G4double u, G4double v) const;
  // calculates the surface point by applying the formula for a surface point. (uses BasisFunction)

  G4double BasisFunctionDerivative(G4int i, G4int k, G4double t,
                                   const std::vector<G4double>& knotVector) const;
  // calculates the derivative of the basis function. (non-recursive but uses BasisFunction)

  std::vector<G4ThreeVector> SurfaceDerivatives(G4double u, G4double v) const;
  // calculates the surface derivative (tangent vector) using the appropriate formula (uses BasisFunctionDerivative)
  // Returns the derivative in the u and v direction (dS_du and dS_dv) as a vector of G4ThreeVectors




  // Functions for optimisation techniques ============================================================
  std::tuple<std::vector<double>, double> ClosestKnot(const G4ThreeVector& point) const;
  // iterates through the unique knots to find the closest knot to a given point.
  // Generally used as initial guess for optimisation techniques.
  // Returns a tuple of < [u,v] , R > where
  //   - [u,v] ->the surface parameters of the closest knot as a vector<double>
  //   - R is -> distance to the closest knot as a double


  std::tuple<std::vector<double>, double> ClosestPointParams(const G4ThreeVector& point) const;
  // Finds the u and v parameters of the closest point on the NURBS surface to a given point by
  // minimising the magnitude (P_surface(u,v) - point(x,y,z)).
  // Returns a tuple of < [u,v], R> where:
  //   - [u,v] -> vector<double> of the optimised u,v surface parameters
  //   - R -> double representing separation between P_surface(u,v) and point(x,y,z) at optimal value

  G4ThreeVector ClosestPoint(const G4ThreeVector& point) const;
  // Calls ClosestPointParams for u, v and returns the surface point as a position vector.
  // Function adds no logic to but improves readability of code.



  std::tuple<std::vector<double>, double, bool> LineIntersectionOpt(const G4ThreeVector& P0,
                                                                   const G4ThreeVector& direction,
                                                                   double line_length,
                                                                   std::vector<double>& uvl_guess) const;
  // Finds optimal u,v (parametrising surface) and l (parametrising line) of the closest intersection of the line
  // passed though in the argument, and the Nurbs surface. Minimises ResidualLineDistance and returns a tuple of:
  //   - [u,v,l] -> opimised values as vector of doubles
  //   - residual -> minimised P_surface -> P_line distance as a double
  //   - successful_opt -> check whether optimisation produced an error.

  void PrintIntersectionOutcome(std::string output_message, std::vector<double> uvl_guess,
                                std::vector<double> uvl_opt, double residual,
                                const G4ThreeVector& P0, const G4ThreeVector& direction) const;
  // print statement for LineIntersectionOpt outcomes

  std::tuple<std::vector<double>,bool> CheckKnotBounds(std::vector<double>& uvl_opt, std::vector<double>& uvl_guess) const;
  // Checks output of optimised u, v values in line intersection procedures to check for convergence issues due to knot
  // boundaries. Returns a tuple of:
  //   - new_uvl_guess -> suggested new initial guess parameters as a vector of doubles
  //   - reached_knot_boundary -> flag if condition met as a bool

  std::tuple<std::vector<double>, double> LineIntersectionParams(const G4ThreeVector& P0,
                                                                 const G4ThreeVector& direction,
                                                                 double line_length) const;
  // Processes LineIntersectionOpt for DistanceToIn or DistanceToIn calculation and applies logic for boundary cases.
  // Optimsed for torus geometry. The argument next_step_indiator tells the function where to go once logic stream is
  // finished. Input 0 for DistanceToIn calculation or 1 for DistanceToOut calculation.
  // Returns a tuple of:
  //   - [uvl_opt] -> validated optimised parameters as vector of doubles
  //   - R_opt -> residual at optimisation as a double



  // Static functions and structs used in optimisation ======================================

  private:
  // #ifdef GEANT4_USE_NLOPT

  // optimisation for ClosestPoint -------------------------------------
  struct ClostestPointContext
  {
    const G4GeorgeNurbs* nurbs;
    G4ThreeVector target_point;
    mutable int call_count = 0;
  }; // struct to pass context for closest point optimisation

  static double ResidualToSurfacePoint(const std::vector<double>& uv,
                                       std::vector<double>& grad,
                                       void* data);
  // calculates the magnitude of the residual three-vector between a surface point (defined by uv)
  // and a target point defined within ClostestPointContext struct. Static function used because
  // NLopt cannot deal with member functions and arguments grad and data are required for nlopt.


  // optimisation for line intersection -------------------------------------
  struct LineIntersectionContext {
    const G4GeorgeNurbs* nurbs;
    G4ThreeVector P0;
    G4ThreeVector direction; // should be unit vector
    mutable int call_count = 0;
    };

  static double ResidualLineDistance(const std::vector<double>& uvl,
                                     std::vector<double>& grad,
                                     void* data);
  // calculates the distance between a point on a line (defined by P0, direction and lambda)
  // and a point on the surface (defined by the prarameters u and v).
  // The first argument is a vector [u,v,l] and defines the parameters to optimise in LineIntersecion




  // #endif





  // Overridden Base class functions ======================================================================
  public:

  EInside Inside(const G4ThreeVector& p) const override;
  G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;

  G4ThreeVector SurfaceNormal(const double u, const double v) const;
  // finds the surface normal at a surface point defined by u and v


  G4double DistanceToIn(const G4ThreeVector& p) const override;
  G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const override;
  double DiscreteLineSearchToIn(const G4ThreeVector& P0, const G4ThreeVector& direction) const;
  // stepwise process to find intersection if inital optimistion for DistacnceToIn(p,v) fails.


  G4double DistanceToOut(const G4ThreeVector& p) const override;
  G4double DistanceToOut(const G4ThreeVector& p,
                         const G4ThreeVector& v,
                         const G4bool calcNorm = false,
                               G4bool* validNorm = nullptr,
                               G4ThreeVector* n = nullptr) const override;
  double DiscreteLineSearchToOut(const G4ThreeVector& P0_start, const G4ThreeVector& direction) const;
  // stepwise process to find intersection if inital optimistion for DistacnceToOut(p,v) fails.

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


  void SetBoundingLimits();

  std::vector<G4ThreeVector> GetBounds() const;
  G4double GetMaxExtent() const;



  void ValidateKnotVectors() const;
  // ensure that knot vectors have length nb_controlPts + degree + 1 in both directions
  // called in the constructor and throws an error if invalid

  void InitialiseBoundingBox();
  void SetBoundingBox(G4Box boundingBoxIn);
  G4Box GetBoundingBox() const;


  void EnableOptVerbose();
  void DisableOptVerbose();




};



#endif // G4GEORGENURBS_HH
