// nurbs class function definitions

#include "G4GeorgeNurbs.hh"


#include "G4BoundingEnvelope.hh"
#include "G4QuickRand.hh"
#include "G4VisExtent.hh"
#include "G4VGraphicsScene.hh"
#include "G4AffineTransform.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include "G4GeomTools.hh"
#include "G4GeometryTolerance.hh"
#include "G4VPVParameterisation.hh"
#include "G4VoxelLimits.hh"

#include "meshdefs.hh"
#include <G4Box.hh>



#include <algorithm>  // for std::unique_copy
#include <iterator>   // for std::back_inserter
#include <tuple> // object to be returned after optimisation
#include "G4TwoVector.hh"


#include <nlopt.hpp>




double G4GeorgeNurbs::LARGE_NUMBER = 1e6;
G4ThreeVector G4GeorgeNurbs::LARGE_THREE_VECTOR = G4ThreeVector( 1e6,  1e6,  1e6);
double G4GeorgeNurbs::CONVERGENCE_TOLERANCE = 1e-6;





G4GeorgeNurbs::G4GeorgeNurbs(const G4String& name,
               const std::vector<std::vector<G4ThreeVector>>& controlPts_in,
               const std::vector<std::vector<G4double>>& weights_in,
               const std::vector<G4double>& knotVectorU_in,
               const std::vector<G4double>& knotVectorV_in,
               const G4int& degreeU_in,
               const G4int& degreeV_in)
  : G4VSolid(name),
  controlPts(controlPts_in),
  weights(weights_in),
  knotVectorU(knotVectorU_in),
  knotVectorV(knotVectorV_in),
  degreeU(degreeU_in),
  degreeV(degreeV_in)
{
  ValidateKnotVectors();
}

G4GeorgeNurbs::~G4GeorgeNurbs() = default;





void G4GeorgeNurbs::PrintVariables() const
{
  std::cout << "Degree U: " << degreeU << "\n";
  std::cout << "Degree V: " << degreeV << "\n\n";

  // Print control points
  std::cout << "Control Points:\n";
  for (size_t i = 0; i < controlPts.size(); i++) {
    for (size_t j = 0; j < controlPts[i].size(); j++) {
      std::cout << "(" << controlPts[i][j].x() << ", "
                << controlPts[i][j].y() << ", "
                << controlPts[i][j].z() << ")  ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  // Print weights
  std::cout << "Weights:\n";
  for (const auto& row : weights) {
    for (const auto& w : row) {
      std::cout << w << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  // Print knot vectors
  std::cout << "Knot Vector U:\n";
  for (const auto& u : knotVectorU) {
    std::cout << u << " ";
  }
  std::cout << "\n\n";

  std::cout << "Knot Vector V:\n";
  for (const auto& v : knotVectorV) {
    std::cout << v << " ";
  }
  std::cout << "\n";



  //test to check that nlopt has been included
  // try {
  //   // Create a dummy 2D optimizer
  //   nlopt::opt test_opt(nlopt::LN_NELDERMEAD, 2);
  //   std::cout << "NLopt is working! Algorithm: "
  //             << test_opt.get_algorithm() << std::endl;
  // }
  // catch (const std::exception& e) {
  //   std::cerr << "NLopt error: " << e.what() << std::endl;
  // }
}




// Functions required to define a NURBS surface =====================================================


G4double G4GeorgeNurbs::BasisFunction(G4int i, G4int k, G4double t,
                                      const std::vector<G4double>& knotVector) const
{
  // Avoid unnecessary computation if t is out of range
  if (!(knotVector[i] <= t && t <= knotVector[i + k + 1])) {
    return 0.0;
  }

  // Base case k = 0
  if (k == 0) {
    // std::cout<<"k=0 case reached "<<std::endl;

    return (knotVector[i] <= t && t <= knotVector[i + 1]) ? 1.0 : 0.0;
  }

  // Compute denominators
  double denom1 = knotVector[i + k] - knotVector[i];
  double denom2 = knotVector[i + k + 1] - knotVector[i + 1];

  // Compute recursive terms
  double term1 = (denom1 > 0) ? ((t - knotVector[i]) / denom1) * BasisFunction(i, k - 1, t, knotVector) : 0.0;
  double term2 = (denom2 > 0) ? ((knotVector[i + k + 1] - t) / denom2) * BasisFunction(i + 1, k - 1, t, knotVector) : 0.0;


  // std::cout<<"\nterm1 + term2: "<< term1 + term2<<std::endl;
  return term1 + term2;
};

G4ThreeVector G4GeorgeNurbs::SurfacePoint(G4double u, G4double v) const
{
  size_t nu = controlPts.size() - 1;      // Final index in u direction
  size_t nv = controlPts[0].size() - 1;   // Final index in v direction

  G4ThreeVector numerator(0.0, 0.0, 0.0);
  double denominator = 0.0;

  for (size_t i = 0; i <= nu; i++) {
    for (size_t j = 0; j <= nv; j++) {
      double N_i = BasisFunction(i, degreeU, u, knotVectorU);
      double N_j = BasisFunction(j, degreeV, v, knotVectorV);
      double weight = weights[i][j] * N_i * N_j;


      numerator += weight * controlPts[i][j]; // Weighted sum
      denominator += weight;

      // std::cout<<"[i,j] = ["<<i<<","<<j<<"]\n"
      // <<"numerator: ["<<numerator.x()<<","<<numerator.y()<<","<<numerator.z() <<"]\n"
      // <<"denominator: ["<<denominator<<"]\n"<<std::endl;
    }
  }

  if (denominator == 0.0) {
    return LARGE_THREE_VECTOR;
  } else {
    return numerator / denominator;
  }
}

G4double G4GeorgeNurbs::BasisFunctionDerivative(G4int i, G4int k, G4double t,
                                                const std::vector<G4double>& knotVector) const
{
  // Base case: the derivative of a zero-degree (k=0) basis function is always zero
  if (k == 0) {
    return 0.0;
  }

  // Compute denominator terms
  double denom1 = knotVector[i + k] - knotVector[i];
  double denom2 = knotVector[i + k + 1] - knotVector[i + 1];

  // Term 1: uses N_{i,k-1}(t)
  double term1 = 0.0;
  if (denom1 > 0.0) {
    double basis_i_k_minus_1 = BasisFunction(i, k - 1, t, knotVector);
    term1 = (k / denom1) * basis_i_k_minus_1;
  }

  // Term 2: uses N_{i+1,k-1}(t)
  double term2 = 0.0;
  if (denom2 > 0.0) {
    double basis_i_plus1_k_minus_1 = BasisFunction(i + 1, k - 1, t, knotVector);
    term2 = (k / denom2) * basis_i_plus1_k_minus_1;
  }

  // The derivative is the difference between the two weighted terms
  return term1 - term2;
}


std::vector<G4ThreeVector> G4GeorgeNurbs::SurfaceDerivatives(G4double u, G4double v) const
{
  size_t nu = controlPts.size() - 1;
  size_t nv = controlPts[0].size() - 1;

  G4ThreeVector numerator_u(0, 0, 0);
  double denominator_u = 0;

  G4ThreeVector numerator_v(0, 0, 0);
  double denominator_v = 0;

  G4ThreeVector numerator(0, 0, 0);
  double denominator = 0;

  for (size_t i = 0; i <= nu; ++i) {
    for (size_t j = 0; j <= nv; ++j) {
      double N_i = BasisFunction(i, degreeU, u, knotVectorU);
      double N_j = BasisFunction(j, degreeV, v, knotVectorV);

      double dN_i = BasisFunctionDerivative(i, degreeU, u, knotVectorU);
      double dN_j = BasisFunctionDerivative(j, degreeV, v, knotVectorV);

      double weight = weights[i][j];
      const G4ThreeVector& P = controlPts[i][j];

      // Weighted basis functions
      double wN   = weight * N_i * N_j;
      double wN_u = weight * dN_i * N_j;
      double wN_v = weight * N_i * dN_j;

      numerator    += wN   * P;
      denominator  += wN;

      numerator_u  += wN_u * P;
      denominator_u += wN_u;

      numerator_v  += wN_v * P;
      denominator_v += wN_v;
    }
  }

  if (denominator == 0.0) {
    return {LARGE_THREE_VECTOR, LARGE_THREE_VECTOR}; // fallback for invalid evaluation
  }

  // Rational NURBS tangent vectors (quotient rule)
  G4ThreeVector tangent_u = (numerator_u * denominator - numerator * denominator_u) / (denominator * denominator);
  G4ThreeVector tangent_v = (numerator_v * denominator - numerator * denominator_v) / (denominator * denominator);

  return {tangent_u, tangent_v};
}



// Functions using optimisation techniques ================================================================

std::tuple<std::vector<double>, double> G4GeorgeNurbs::ClosestKnot(const G4ThreeVector& point) const
{

  double smallest_R = LARGE_NUMBER;
  std::vector<double> u_v_smallest_R{0.0, 0.0};

  // Unique knot vectors
  std::vector<G4double> unique_u_knots;
  std::vector<G4double> unique_v_knots;

  std::unique_copy(knotVectorU.begin(), knotVectorU.end(), std::back_inserter(unique_u_knots));
  std::unique_copy(knotVectorV.begin(), knotVectorV.end(), std::back_inserter(unique_v_knots));


  for (size_t i = 0; i < unique_u_knots.size(); ++i)
  {
    for (size_t j = 0; j < unique_v_knots.size(); ++j)
    {
      // cycles through each unique knot combination and calculates the distance
      // between that surface point and the given point
      G4ThreeVector P_s_knot = SurfacePoint(unique_u_knots[i], unique_v_knots[j]);
      G4ThreeVector residual_vec = P_s_knot - point;
      double separation = residual_vec.mag();

      // Debugging output
      // G4cout << "\nuknot[" << i << "] = " << unique_u_knots[i] << G4endl;
      // G4cout << "vknot[" << j << "] = " << unique_v_knots[j] << G4endl;
      // G4cout << "separation = " << separation << G4endl;

      if (separation < smallest_R)
      {
        // updates u_v_smallest_R if separation with these parameters is smaller than the current smallest_R

        smallest_R = separation;
        u_v_smallest_R[0] = unique_u_knots[i];
        u_v_smallest_R[1] = unique_v_knots[j];

        // G4cout << "====entered if statement======" << G4endl;
      }
    }
  }
  return std::make_tuple(u_v_smallest_R, smallest_R);
}



std::tuple<std::vector<double>, double> G4GeorgeNurbs::ClosestPointParams(const G4ThreeVector& point) const
{
  // Finds the closest point on the surface to a given point in space

  // 1) set initial guess as closest knot
  auto [uv_closes_knot, R] = ClosestKnot(point);
  std::vector<double> uv_guess = uv_closes_knot;

  // 2) Set up the NLopt optimizer (Nelder-Mead, 2D problem)
  nlopt::opt optimizer(nlopt::LN_NELDERMEAD, 2);
  optimizer.set_xtol_rel(CONVERGENCE_TOLERANCE);

  // 3) Prepare context with reference to this surface and the target point
  OptimizationContext context = {this, point};

  // 4) Set objective function and pass context
  optimizer.set_min_objective(G4GeorgeNurbs::ResidualToSurfacePoint, &context);

  // 5) ensure u and v don't go out of bounds
  optimizer.set_lower_bounds({knotVectorU.front(), knotVectorV.front()});
  optimizer.set_upper_bounds({knotVectorU.back(), knotVectorV.back()});

  // 6) optimising
  std::vector<double> uv = uv_guess;
  double optimised_residual;

  try {
    optimizer.optimize(uv, optimised_residual); // changes value of uv during optimisation


    return std::make_tuple(uv, optimised_residual); // Return closest surface point
  }
  catch (const std::exception& e) {
    // In case of failure, report and return a fallback value
    G4cerr << "NLopt error in ClosestPoint: " << e.what() << G4endl;
    return std::make_tuple(std::vector<double>{LARGE_NUMBER, LARGE_NUMBER}, LARGE_NUMBER);
  }
}

G4ThreeVector G4GeorgeNurbs::ClosestPoint(const G4ThreeVector& point) const
{
  auto [uv_opt, R_opt] = ClosestPointParams(point);
  return SurfacePoint(uv_opt[0], uv_opt[1]);
}


G4ThreeVector G4GeorgeNurbs::LineIntersection(const G4ThreeVector& P0,
                                              const G4ThreeVector& direction,
                                              double lambda_bound) const
{

  std::vector<double> uvl = {0.5, 0.5, 0.0}; // initial guess

  nlopt::opt opt(nlopt::LN_NELDERMEAD, 3);
  opt.set_xtol_rel(1e-6);

  // Set bounds for (u, v, lambda)
  opt.set_lower_bounds({knotVectorU.front(), knotVectorV.front(), -std::abs(lambda_bound)});
  opt.set_upper_bounds({knotVectorU.back(),  knotVectorV.back(),  std::abs(lambda_bound)});

  LineIntersectionContext context = {
    this,
    P0,
    direction.unit(),  // ensure it's a unit vector
    lambda_bound
  };

  opt.set_min_objective(G4GeorgeNurbs::ResidualLineDistance, &context);

  double d_opt;
  try {
    opt.optimize(uvl, d_opt);
    double u_opt = uvl[0];
    double v_opt = uvl[1];
    double lambda_opt = uvl[2];

    if (std::abs(d_opt) > 1e-6 ||
        std::signbit(lambda_opt) != std::signbit(lambda_bound) ||
        std::abs(lambda_opt) > std::abs(lambda_bound)) {
      return LARGE_THREE_VECTOR;
        }

    return SurfacePoint(u_opt, v_opt);
  }
  catch (const std::exception& e) {
    G4cerr << "NLopt error in LineIntersection: " << e.what() << G4endl;
    return LARGE_THREE_VECTOR;
  }
}




// Static functions and structs used in optimisation =============================================================

// #ifdef GEANT4_USE_NLOPT //
double G4GeorgeNurbs::ResidualToSurfacePoint(const std::vector<double>& uv, std::vector<double>& grad, void* data)
{
  // definition of static function used by NLopt to compute the distance (residual)
  // between a surface point and a fixed target point.

  (void)grad; // to suppress warning on build


  // Unpacking context (contains the target point and pointer to the NURBS surface)
  const auto* context = static_cast<G4GeorgeNurbs::OptimizationContext*>(data);
  const G4GeorgeNurbs* nurbs = context->nurbs;
  const G4ThreeVector& target_point = context->target_point;

  // find the surface point at uv (passed as argument)
  G4ThreeVector surf_pt = nurbs->SurfacePoint(uv[0], uv[1]);

  // calculating the distance between the surface point and target point
  G4ThreeVector residual = surf_pt - target_point;
  return residual.mag();
}



double G4GeorgeNurbs::ResidualLineDistance(const std::vector<double>& uvl,
                                           std::vector<double>& grad,
                                           void* data)
{
  // function to optimise in LineIntersection.
  const auto* context = static_cast<LineIntersectionContext*>(data);

  double u = uvl[0]; // u and v parametrises a point on the surface
  double v = uvl[1];
  double l = uvl[2]; // l parametrises a point on the line

  G4ThreeVector P_nurbs = context->nurbs->SurfacePoint(u, v);
  G4ThreeVector P_line = context->P0 + l * context->direction;

  (void)grad; // to supress warning on build

  return (P_line - P_nurbs).mag(); // return distance between those two points
}







// #endif // GEANT4_USE_NLOPT


// Overridden Base class functions ====================================================================

EInside G4GeorgeNurbs::Inside(const G4ThreeVector& p) const
{
  // finding the closest surface point to p
  auto [uv_closest_pt, R_closest_pt] = ClosestPointParams(p); // u,v of closest point (can maybe just use ClosestKnot)
  G4ThreeVector closest_surf_pt = SurfacePoint(uv_closest_pt[0], uv_closest_pt[1]); // [x,y,z] of closest pt


  // finding the residual from p --> closest_surf_pt
  G4ThreeVector residual = closest_surf_pt - p;

  if (residual.mag() < CONVERGENCE_TOLERANCE)
  {
    // if the residual is of order of the CONVERGENCE_TOLERANCE, then p must be on the surface.
    return kSurface;
  }

  // finding normal at the surface point and taking dot product with the residual
  G4ThreeVector normal_at_surf_pt = SurfaceNormal(uv_closest_pt[0], uv_closest_pt[1]);
  G4double residual_dot_normal = residual.dot(normal_at_surf_pt);

  // inside/outside can be determined from the sign of the dot product
  if (residual_dot_normal>0){return kInside;}
  return kOutside;
}


G4ThreeVector G4GeorgeNurbs::SurfaceNormal(const G4ThreeVector& p) const
{
  // finding the u,v parameters of the surface point with ClosestPointParams function.
  // residual_from_p should be 0.
  auto [uv_from_p, residual_from_p] = ClosestPointParams(p);
  G4ThreeVector n = SurfaceNormal(uv_from_p[0], uv_from_p[1]);
  return n;
}



G4ThreeVector G4GeorgeNurbs::SurfaceNormal(const double u, const double v) const
{
  std::vector<G4ThreeVector> tangents = SurfaceDerivatives(u,v);

  G4ThreeVector& tangent_u = tangents[0];
  G4ThreeVector& tangent_v = tangents[1];

  G4ThreeVector normal = tangent_u.cross(tangent_v);
  normal = normal.unit();
  return normal;
}



G4double G4GeorgeNurbs::DistanceToIn(const G4ThreeVector& p) const
{
  // first check if point is already inside or on surface, return 0 if so
  EInside inside_status = Inside(p);
  if (inside_status == kInside || inside_status == kSurface) {
    return 0;
  }

  // find the closest surface point using ClosestPointParams function
  auto [uv, residual_from_p] = ClosestPointParams(p);

  // return the minimised distance
  return residual_from_p;
}

G4double G4GeorgeNurbs::DistanceToIn(const G4ThreeVector& p0, const G4ThreeVector& v) const
{
  return 0.0;
}

G4double G4GeorgeNurbs::DistanceToOut(const G4ThreeVector& p) const
{
  // first check if point is already outside, return 0 if so
  EInside inside_status = Inside(p);
  if (inside_status == kOutside) {
    return 0;
  }

  // find the closest point on the surface using ClosestPointParams function
  auto [uv, residual_from_p] = ClosestPointParams(p);

  // return the minimised distance
  return residual_from_p;


  return 0.0;
}

G4double G4GeorgeNurbs::DistanceToOut( const G4ThreeVector& p,const G4ThreeVector& v,
                                       const G4bool calcNorm,
                                        G4bool* validNorm,
                                        G4ThreeVector* n ) const
{
  return 0.0;
}

G4double G4GeorgeNurbs::GetCubicVolume()
{
  return 5;
}

G4double G4GeorgeNurbs::GetSurfaceArea()
{
  return 10;
}

G4ThreeVector G4GeorgeNurbs::GetPointOnSurface() const
{
  return G4ThreeVector(0, 0, 0);
}

G4VisExtent G4GeorgeNurbs::GetExtent() const
{
  G4double xMin, xMax, yMin, yMax, zMin, zMax;

  xMin = -1;
  xMax = 1;
  yMin = -1;
  yMax = 1;
  zMin = -1;
  zMax = 1;

  return { xMin, xMax, yMin, yMax, zMin, zMax};
}



G4bool G4GeorgeNurbs::CalculateExtent( const EAxis pAxis,
                                  const G4VoxelLimits& pVoxelLimit,
                                  const G4AffineTransform& pTransform,
                                        G4double& pMin, G4double& pMax ) const
{
  G4ThreeVector bmin, bmax;

  // Get bounding box
  BoundingLimits(bmin,bmax);

  // Find extent
  G4BoundingEnvelope bbox(bmin,bmax);
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
}

void G4GeorgeNurbs::DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
  scene.AddSolid (*this);
}

std::ostream& G4GeorgeNurbs::StreamInfo( std::ostream& os ) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4GeorgeNurbs\n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}


G4GeometryType G4GeorgeNurbs::GetEntityType() const
{
  return {"G4GeorgeNurbs"};
}




void G4GeorgeNurbs::ValidateKnotVectors() const
{
  std::size_t expected_knot_length_U = controlPts.size() + degreeU + 1;
  std::size_t expected_knot_length_V = controlPts[0].size() + degreeV + 1;

  if (knotVectorU.size() != expected_knot_length_U || knotVectorV.size() != expected_knot_length_V)
  {
    std::ostringstream oss;
    if (knotVectorU.size() != expected_knot_length_U)
    {
      oss << "Invalid U knot vector length: expected " << expected_knot_length_U
          << ", got " << knotVectorU.size() << ".\n";
    }
    if (knotVectorV.size() != expected_knot_length_V)
    {
      oss << "Invalid V knot vector length: expected " << expected_knot_length_V
          << ", got " << knotVectorV.size() << ".\n";
    }

    G4Exception("G4GeorgeNurbs", "InvalidKnotVectorSize",
                FatalException, oss.str().c_str());
  }
}