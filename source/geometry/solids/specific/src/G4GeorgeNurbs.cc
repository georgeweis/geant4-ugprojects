// nurbs class function definitions

#include "G4GeorgeNurbs.hh"


#include "G4BoundingEnvelope.hh"
#include "G4QuickRand.hh"
#include "Randomize.hh"
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
double G4GeorgeNurbs::SURFACE_TOLERANCE = 1e-1;





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
  SetBoundingLimits();
  optVerbose = false;
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

  /*attempt to use more sample points to stop optimiser failure issues*/
  auto linspace = [](double start, double end, int num) -> std::vector<double> { // lambda funciton to make a quick linspace
        std::vector<double> result;
        if (num == 0) return result;
        if (num == 1) {
            result.push_back(start);
            return result;
        }

        double step = (end - start) / (num - 1);
        for (int i = 0; i < num; ++i) {
            result.push_back(start + i * step);
        }
        return result;
    };

  // create twice as many sample points
  auto sample_u_knots = linspace(knotVectorU.front(), knotVectorU.back(), unique_u_knots.size()*4);
  auto sample_v_knots = linspace(knotVectorV.front(), knotVectorV.back(), unique_u_knots.size()*4);

  // in case of floating point error propogation
  sample_u_knots.back() = unique_u_knots.back();
  sample_v_knots.back() = unique_v_knots.back();

  double separation;
  for (size_t i = 0; i < sample_u_knots.size(); ++i)
  {
    for (size_t j = 0; j < sample_v_knots.size(); ++j)
    {
      // cycles through each unique knot combination and calculates the distance
      // between that surface point and the given point
      G4ThreeVector P_s_knot = SurfacePoint(sample_u_knots[i], sample_v_knots[j]);
      G4ThreeVector residual_vec = P_s_knot - point;
      separation = residual_vec.mag();

      // Debugging output
      // G4cout<<"--------------"<<G4endl;      // G4cout<<"knot: ("<<unique_u_knots[i]<<", "<<unique_v_knots[j]<<"),  r = ("<<separation<<")\n";      // G4cout<<P_s_knot<<G4endl;

      if (separation < smallest_R)
      {
        //G4cout<<"<----------entered if statment"<<G4endl;

        // updates u_v_smallest_R if separation with these parameters is smaller than the current smallest_R
        smallest_R = separation;
        u_v_smallest_R[0] = sample_u_knots[i];
        u_v_smallest_R[1] = sample_v_knots[j];
      }
    }
  }


  // test if knots are at the end point, then do a further check to see which side of the surface is closer.
  // Saves convergence issues due to knot vector boundaries when ClosestKnot is used as a start guess.
  double residual_near_front;
  double residual_near_back;
  if( u_v_smallest_R[0] == sample_u_knots.front() || u_v_smallest_R[0] == sample_u_knots.back())
  {
    // if u is at the end of the knot vector, compare residuals of the knots that are 1 knot away
    int size_u = unique_u_knots.size();
    residual_near_front = (point - SurfacePoint(sample_u_knots[1], u_v_smallest_R[1])).mag();
    residual_near_back = (point - SurfacePoint(sample_u_knots[size_u-2], u_v_smallest_R[1])).mag();

     u_v_smallest_R[0] = (residual_near_front < residual_near_back) // if residual near front of u is closer
                      ? sample_u_knots.front() // return the first knot (usually 0)
                      : sample_u_knots.back(); // else return the last knot (usually 1)
  }

  if( u_v_smallest_R[1] == sample_v_knots.front() || u_v_smallest_R[1] == sample_v_knots.back())
  {
    // if v is at the end of the knot vector, compare residuals of the knots that are 1 knot away
    int size_v = sample_v_knots.size();
    residual_near_front = (point - SurfacePoint(u_v_smallest_R[0], unique_v_knots[1])).mag();
    residual_near_back = (point - SurfacePoint(u_v_smallest_R[0], unique_v_knots[size_v-2])).mag();
    //std::cout<<"BEFORE u_v_smallest_R: ("<<u_v_smallest_R[0] << " "<<u_v_smallest_R[1]<<")"<<std::endl;
    //std::cout<<"residual_near_front: " << residual_near_front << std::endl;    //std::cout<<"residual_near_back: " << residual_near_back << std::endl;

    u_v_smallest_R[1] = (residual_near_front < residual_near_back) // if residual near front of v is closer
                      ? sample_v_knots.front() // return the first knot (usually 0)
                      : sample_v_knots.back(); // else return the last knot (usually 1)
    //std::cout<<"AFTER u_v_smallest_R: ("<<u_v_smallest_R[0] << " "<<u_v_smallest_R[1]<<")"<<std::endl;
  }

  //G4cout<<"closest knot from indide closest knot function: ("<<u_v_smallest_R[0]<<", "<<u_v_smallest_R[1]<<"),  r = ("<<smallest_R<<")\n";
  return std::make_tuple(u_v_smallest_R, smallest_R);
}




std::tuple<std::vector<double>, double> G4GeorgeNurbs::ClosestPointParams(const G4ThreeVector& point) const
{
  // Finds the closest point on the surface to a given point in space

  // 1) set initial guess as closest knot
  auto [uv_closes_knot, R] = ClosestKnot(point);
  std::vector<double> uv_guess = uv_closes_knot;

  // 2) Set up the NLopt optimizer (Nelder-Mead, 2D problem)
  nlopt::opt optimizer(nlopt::LN_COBYLA, 2);
  optimizer.set_xtol_rel(CONVERGENCE_TOLERANCE);

  // 3) Prepare context with reference to this surface and the target point
  ClostestPointContext context = {this, point};

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


std::tuple<std::vector<double>, double, bool> G4GeorgeNurbs::LineIntersectionOpt(const G4ThreeVector& P0,
                                                                                const G4ThreeVector& direction,
                                                                                double line_length,
                                                                                std::vector<double>& uvl_guess) const
{
  // set up optimiser
  nlopt::opt opt(nlopt::LN_BOBYQA, 3);
  opt.set_xtol_rel(CONVERGENCE_TOLERANCE);

  // Set bounds for (u, v, l)
  opt.set_lower_bounds({knotVectorU.front(), knotVectorV.front(), 0}); // lower bounds uv usually 0 and l always 0
  opt.set_upper_bounds({knotVectorU.back(),  knotVectorV.back(),  line_length}); // upper bounds are the end of the
                                                                               // knot vectors and the argument line_length

  // set context of the line.
  LineIntersectionContext context = {
    this,
    P0,
    direction.unit()  // ensure it's a unit vector
  };

  opt.set_min_objective(G4GeorgeNurbs::ResidualLineDistance, &context);



  // starting guess for uvl_opt is passed throgh in the argument
  std::vector<double> uvl_opt  = uvl_guess;

  double residual_opt; // parameter to minimise: distance between P_surface(u,v) and P_line(l)
  bool successful_optimisation; // result of optimisation
  try {

    opt.optimize(uvl_opt, residual_opt);

    if (optVerbose){
    G4cout << "LineIntersectionOpt: objective evaluated "
         << context.call_count << " times." << G4endl;
    }


    // return successfully optimised params with exit_status 0
    successful_optimisation = true;
    return std::make_tuple(uvl_opt, residual_opt, successful_optimisation);

  }
  catch (const std::exception& e) {
    G4cerr << "NLopt error in LineIntersection: " << e.what() << G4endl;
    successful_optimisation = false;
    return std::make_tuple(std::vector<double>{LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER},
                           LARGE_NUMBER,
                           successful_optimisation);
  }
}


void G4GeorgeNurbs::PrintIntersectionOutcome(std::string output_message,
                                             std::vector<double> uvl_guess,
                                             std::vector<double> uvl_opt,
                                             double residual,
                                             const G4ThreeVector& P0,
                                             const G4ThreeVector& direction) const
{
  G4ThreeVector P_surface = SurfacePoint(uvl_opt[0], uvl_opt[1]);
  G4ThreeVector P_line = (P0 + uvl_opt[2]*direction);

  std::cout<<output_message<<"\n"
           <<"Initial uvl guess: ("<<uvl_guess[0]<<", "<<uvl_guess[1]<<", "<<uvl_guess[2]<<")\n"
           <<"P_surface (u,v): ("<<uvl_opt[0]<<","<< uvl_opt[1]<<")"<<"\n"
           <<"P_surface (x,y,z): "<<P_surface<<"\n"
           <<"P_line: "<<P_line<<"\n"
           <<"Convergence l = "<<uvl_opt[2] <<"\n"
           <<"Convergence residual = "<<residual<<"\n"

           <<"\nFor python: \n"
           <<"# line points\n"
           <<"p0 = np.array(["<<P0.x()<<","<<P0.y()<<","<<P0.z()<<"], dtype = float)"<<"\n"
           <<"n_line = np.array(["<<direction.x()<<","<<direction.y()<<","<<direction.z()<<"], dtype = float)"<<"\n"
           <<"# points from convergence\n"
           <<"P_suf = np.array(["<<P_surface.x()<<","<<P_surface.y()<<","<<P_surface.z()<<"], dtype = float)"<<"\n"
           <<"P_line = np.array(["<<P_line.x()<<","<<P_line.y()<<","<<P_line.z()<<"], dtype = float)"<<"\n"
           <<"# l guess\n"
           <<"new_l_guess = "<<uvl_guess[2]<<"\n"
           <<std::endl;

}

std::tuple<std::vector<double>,bool> G4GeorgeNurbs::CheckKnotBounds(std::vector<double>& uvl_opt,
                                                                    std::vector<double>& uvl_guess) const
{
  std::vector<double> new_uvl_guess = uvl_guess; // local copy to be edited with following if statements
  bool reached_knot_boundary = false; // flag to indicate if the optimised u,v params are at the end of their vector.

  // if u,v params are at the end of the knot vector, change initial guess to other end and set reached_knot_boundary to true
  double tolerance = 0.05;
  if(uvl_opt[0] <(knotVectorU.front() + tolerance)){
    new_uvl_guess[0] = knotVectorU.back(); reached_knot_boundary = true;
  }
  else if(uvl_opt[0] >(knotVectorU.back() - tolerance)){
    new_uvl_guess[0] = knotVectorU.front(); reached_knot_boundary = true;
  }

  if(uvl_opt[1] <(knotVectorV.front() + tolerance)){
    new_uvl_guess[1] = knotVectorV.back(); reached_knot_boundary = true;
  }
  else if(uvl_opt[1] >(knotVectorV.back() - tolerance)){
    new_uvl_guess[1] = knotVectorV.front(); reached_knot_boundary = true;
  }

  if (reached_knot_boundary)
  {
    if (optVerbose) {std::cout<<"Knot boundary checks failed."<<std::endl;}

    return std::make_tuple(new_uvl_guess, reached_knot_boundary); // return new inital guess parameters
  }

  if (optVerbose) {std::cout<<"Knot boundary checks passed."<<std::endl;}

  return std::make_tuple(uvl_guess, reached_knot_boundary); // return original inital guess parameters
}




std::tuple<std::vector<double>, double> G4GeorgeNurbs::LineIntersectionParams(const G4ThreeVector& P0,
                                                                              const G4ThreeVector& direction,
                                                                              double line_length,
                                                                              int next_step_indiator) const
{
  // STEP 1: Initial optimisation from start point of line


  // declaring variables to be overwritten

  std::vector<double> uvl_opt;
  double residual_opt;
  bool opt_success;
  std::tuple<std::vector<double>, double> result;
  std::string output_message;

  auto [closest_knot_start_point, R_knot] = ClosestKnot(P0); // guess for uv is the closest knot to P0

  double initial_l_guess = 0;
  if(R_knot > maxExtent/2) // P0 is not close to the surface
  {
    initial_l_guess = R_knot-maxExtent/3; // choose inital l thats closer
  }


  std::vector<double> uvl_guess = {closest_knot_start_point[0], closest_knot_start_point[1] , initial_l_guess};
  double line_length_upper_lim = R_knot;
  std::tie(uvl_opt, residual_opt, opt_success) = LineIntersectionOpt(P0, direction, line_length_upper_lim, uvl_guess);

  // std::cout<<"======================== opt success poutcome: "<<opt_success<<std::endl;


  /* Check 1: is this a intersection point (residual = 0)? */
  if (residual_opt<SURFACE_TOLERANCE) // yes to check 1
  {
    /* Check 2: is this intersection at P0 (meaning the line started ON the surface)?  */
    if (uvl_opt[2] < 2*SURFACE_TOLERANCE) // yes to check 2
    {
      if (optVerbose) {
          output_message = "❌ Convergence to a non-intersecting point at end of STEP 1. ";
          PrintIntersectionOutcome(output_message, uvl_guess, uvl_opt, residual_opt, P0, direction);
      }

      if (next_step_indiator == 0)
      {
        // apply logic for DistanceToIn
        result = PushLGuessForDTI(uvl_opt, maxExtent/4,  P0, direction, line_length);
      }
      else if (next_step_indiator == 1)
      {
        // apply logic for DistanceToOut
        result = PushLGuessForDTO(uvl_opt, P0, direction, line_length);
      }
      return result;
    }

    // yes to check 1 and no to check 2: intersection is likely the first intersection
    if (optVerbose) {
      output_message = "✅ Successful convergence to an intersection at STEP 1, Check 2.";
      PrintIntersectionOutcome(output_message, uvl_guess, uvl_opt, residual_opt, P0, direction);
    }
    return std::make_tuple(uvl_opt, residual_opt);

  }

  /* Check 3: does u_opt or v_opt = the ends of their knot vectors?  */
  auto [uvl_guess_knot, reached_knot_bounds] = CheckKnotBounds(uvl_opt, uvl_guess);
  uvl_guess = uvl_guess_knot;




  if (reached_knot_bounds) // yes to check 3
  {
    // redo optimisation with new initial knot guesses
    std::tie(uvl_opt, residual_opt, opt_success) = LineIntersectionOpt(P0, direction, line_length, uvl_guess);

    /* Check 4: is this now an intersection point (residual = 0)?  */
    if (residual_opt<SURFACE_TOLERANCE) // yes to check 4
    {
      /* Check 5: is this intersection at P0 (meaning the line started ON the surface)?  */
      if (uvl_opt[2] > SURFACE_TOLERANCE) // no to check 5
      {
        if (optVerbose) {
          output_message = "✅ Successful convergence to an intersection at STEP 1, check 3 knot vector change.";
          PrintIntersectionOutcome(output_message, uvl_guess, uvl_opt, residual_opt, P0, direction);
        }
        return std::make_tuple(uvl_opt, residual_opt); // accept point
      }
    }
  }


  // All logic streams now require the intial guess value of l to be push along the line to search for a new intersection.
  // Possible cases: optimised to start of line as a surface point or optimised to a close approach. Neither are intersections.
  // Done in a separate function but Check and STEP numbers follow on from this point.
  if (optVerbose) {
    output_message = "❌ Convergence to a non-intersecting point at end of STEP 1. ";
    PrintIntersectionOutcome(output_message, uvl_guess, uvl_opt, residual_opt, P0, direction);
  }

  if (next_step_indiator == 0)
      {
        // apply logic for DistanceToIn
        result = PushLGuessForDTI(uvl_opt, maxExtent/4,  P0, direction, line_length);
      }
      else if (next_step_indiator == 1)
      {
        // apply logic for DistanceToOut
        result = PushLGuessForDTO(uvl_opt, P0, direction, line_length);
      }
  return result;

}

std::tuple<std::vector<double>, double> G4GeorgeNurbs::PushLGuessForDTI(std::vector<double> uvl_opt_old,
                                                                  double push_length,
                                                                  const G4ThreeVector& P0,
                                                                  const G4ThreeVector& direction,
                                                                  double line_length) const
{
  // STEP 2: push the initial guess along the line and optimise from the 'pushed' point

  (void)push_length; // supress for now, may be used in future development to vary push length based on the nurbs Geometry

  // declaring variables to be overwritten
  std::vector<double> uvl_pushed_opt; // optimised result after pushing
  double residual_pushed_opt; // optimised residual after pushing
  bool opt_success; // successful optimisation flag (unused in checks)

  std::vector<double> closest_knot_pushed;
  double R;
  std::string output_message;

  // push and optimise
  double l_pushed = uvl_opt_old[2] + maxExtent/4;
  G4ThreeVector P_pushed = P0 + l_pushed*direction; // new location from l_pushed
  std::tie(closest_knot_pushed, R) = ClosestKnot(P_pushed);
  std::vector<double> uvl_pushed_guess = {closest_knot_pushed[0], closest_knot_pushed[1] , l_pushed};
  std::tie(uvl_pushed_opt, residual_pushed_opt, opt_success) = LineIntersectionOpt(P0, direction, line_length, uvl_pushed_guess);

  if (optVerbose) {
    output_message = "🌀 Output after first push by "+ std::to_string(l_pushed);
    PrintIntersectionOutcome(output_message, uvl_pushed_guess, uvl_pushed_opt, residual_pushed_opt, P0, direction);
  }

   /* Check 6: has this optimisation with a pushed l landed on the same spot?  */
  if (std::abs(uvl_opt_old[2]-uvl_pushed_opt[2]) < 2*SURFACE_TOLERANCE) // yes to check 6
  {
    /* Check 7: do both of these optimised l values = 0?  */
    if (uvl_opt_old[2] < SURFACE_TOLERANCE) // yes to check 7
    {
      if (optVerbose) {
        std::cout <<"❌ line found to be pointing away from all points on surface at STEP 2, Check 7. \n"
                  <<"Returning LARGE_NUMBERs."<<std::endl;
      output_message = "Last intersection. ";
      PrintIntersectionOutcome(output_message, uvl_pushed_guess, uvl_pushed_opt, residual_pushed_opt, P0, direction);
      }

      return std::make_tuple(std::vector<double>(3, LARGE_NUMBER), LARGE_NUMBER);
    }
    // push again by 2 times the origional push length and optimise again from this new pushed point
    l_pushed = uvl_opt_old[2] + maxExtent/2; // new l is the
    P_pushed = P0 + l_pushed*direction; // new location from l_pushed
    std::tie(closest_knot_pushed, R) = ClosestKnot(P_pushed);
    uvl_pushed_guess = {closest_knot_pushed[0], closest_knot_pushed[1] , l_pushed};
    std::tie(uvl_pushed_opt, residual_pushed_opt, opt_success) = LineIntersectionOpt(P0, direction, line_length, uvl_pushed_guess);

    if (optVerbose) {
      output_message = "🌀 Output after second pushing by "+ std::to_string(l_pushed);
      PrintIntersectionOutcome(output_message, uvl_pushed_guess, uvl_pushed_opt, residual_pushed_opt, P0, direction);
    }
  }




  /* Check 8: is this new pushed minimum an intersction (residual = 0)?  */
  if (residual_pushed_opt<SURFACE_TOLERANCE) // yes to check 8
  {
    if (optVerbose) {
      output_message = "❓ Intersection found at STEP 2, check 8. Requires check to see if its the closest. ";
      PrintIntersectionOutcome(output_message, uvl_pushed_guess, uvl_pushed_opt, residual_pushed_opt, P0, direction);
    }

    auto result = FirstIntersectionTestForDTI(uvl_opt_old, uvl_pushed_opt, P0,direction, line_length);
    return result;
  }

  /* Check 9: do the optimised u or v = the ends of their knot vectors?  */

  auto [uvl_pushed_guess_knot, reached_knot_bounds] = CheckKnotBounds(uvl_pushed_opt, uvl_pushed_guess);
  uvl_pushed_guess = uvl_pushed_guess_knot;
  uvl_pushed_guess[2] = uvl_pushed_opt[2]; // a close point was found, start with same l guess

  if (reached_knot_bounds) // yes to check 9
  {
    // redo optimisation with new initial knot guesses
    std::tie(uvl_pushed_opt, residual_pushed_opt, opt_success) = LineIntersectionOpt(P0, direction, line_length, uvl_pushed_guess);
  }

  /* Check 10: is this new pushed value an intersection (residual = 0)?  */
  if (residual_pushed_opt>SURFACE_TOLERANCE) // no to check 10
  {
    if (optVerbose) {
      std::cout <<"❌ line does not intersect at any points, confirmed at STEP 2, Check 10. \n"
                <<"Returning LARGE_NUMBERs."<<std::endl;
      output_message = "Last intersection. ";
    PrintIntersectionOutcome(output_message, uvl_pushed_guess, uvl_pushed_opt, residual_pushed_opt, P0, direction);
    }
    return std::make_tuple(std::vector<double>(3, LARGE_NUMBER), LARGE_NUMBER);
  }

  // All logic streams now require a test to confirm that the found intersection is actually the first intersection.
  // Done in function FirstIntersectionTestForDTI, with test and check numbers continuing

  if (optVerbose) {
    output_message = "❓ Intersection should have been found at the end of STEP 2. Current values: ";
    PrintIntersectionOutcome(output_message, uvl_pushed_guess, uvl_pushed_opt, residual_pushed_opt, P0, direction);
  }

  auto result = FirstIntersectionTestForDTI(uvl_opt_old, uvl_pushed_opt, P0,direction, line_length);

  return result;
}



std::tuple<std::vector<double>, double> G4GeorgeNurbs::FirstIntersectionTestForDTI(std::vector<double> uvl_opt_old,
                                                                             std::vector<double> uvl_opt_pushed,
                                                                             const G4ThreeVector& P0,
                                                                             const G4ThreeVector& direction,
                                                                             double line_length) const
{
  // STEP 3: check whether the intersection found after pushing is the first intersection. If not, find it.

  // declaring variables
  std::vector<double> uvl_check_opt; // optimised result after pushing
  double residual_check_opt; // optimised residual after pushing
  bool opt_success; // successful optimisation flag (unused in checks)
  std::string output_message;


  // optimise from a new point at some fraction of the way between l_opt_old and l_opt_pushed
  double l_check = uvl_opt_old[2] + 0.6*uvl_opt_pushed[2]; // a little over half the separation
  G4ThreeVector P_check = P0 + l_check*direction; // new location from l_check
  auto [closest_knot_check, R] = ClosestKnot(P_check);
  std::vector<double> uvl_check_guess = {closest_knot_check[0], closest_knot_check[1] , l_check};
  std::tie(uvl_check_opt, residual_check_opt, opt_success) = LineIntersectionOpt(P0, direction, line_length, uvl_check_guess);



  /* Check 11: does this optimise to the same point?  */
  if (std::abs(uvl_opt_pushed[2]-uvl_check_opt[2]) < 2*SURFACE_TOLERANCE) // yes to check 11
  {
    if (optVerbose) {
      output_message = std::string("✅ Successful convergence to a validated intersection at STEP 3, Check 11. \n" )+
                       "(check point same as pushed, accept this point)";
    PrintIntersectionOutcome(output_message, uvl_check_guess, uvl_check_opt, residual_check_opt, P0, direction);
    }

    return std::make_tuple(uvl_check_opt, residual_check_opt); // accept point (uvl_check_opt=uvl_opt_old here)
  }

  /* Check 12: is this point an intersection? */
  if (residual_check_opt<SURFACE_TOLERANCE) // yes to check 12
  {
    if (optVerbose) {
      output_message = std::string("✅ Successful convergence to a validated intersection at STEP 3, Check 12.\n") +
                       "(check point is a different to pushed, accept check point)";
      PrintIntersectionOutcome(output_message, uvl_check_guess, uvl_check_opt, residual_check_opt, P0, direction);
    }
    return std::make_tuple(uvl_check_opt, residual_check_opt); // accept check point
  }

  /* Check 13: do the optimised u or v = the ends of their knot vectors? */

  auto [new_uvl_check_guess, reached_knot_bounds] = CheckKnotBounds(uvl_check_opt, uvl_check_guess);
  uvl_check_guess = new_uvl_check_guess;
  if (!reached_knot_bounds) // no to check 13
  {
    if (optVerbose) {
      output_message = std::string("✅ Successful convergence to a validated intersection at STEP 3, Check 13. \n" )+
                       "(check point is a close approach, accept pushed point)";
      std::vector<double> unknown_initial_guess = {101,101,101};
      PrintIntersectionOutcome(output_message, unknown_initial_guess, uvl_opt_pushed, 0, P0, direction); // inital guess lost
    }
    return std::make_tuple(uvl_opt_pushed, 0); // accept pushed point. Residual value not known here so assumed as 0
  }

  if (reached_knot_bounds) // yes to check 13
  {
    // redo optimisation with new initial knot guesses
    std::tie(uvl_check_opt, residual_check_opt, opt_success) = LineIntersectionOpt(P0, direction, line_length, uvl_check_guess);
  }

  /* Check 14: after changing knot vectors, is this check point an intersection? */
  if (residual_check_opt<SURFACE_TOLERANCE) //yes to check 14
  {
    if (optVerbose) {
      output_message = std::string("✅ Successful convergence to a validated intersection at STEP 3, Check 14. \n")+
                       "(check point is a different to pushed, accept pushed point)";
      PrintIntersectionOutcome(output_message, uvl_check_guess, uvl_check_opt, residual_check_opt, P0, direction);
    }
    return std::make_tuple(uvl_check_opt, residual_check_opt); // accept check point
  }

  if (residual_check_opt>SURFACE_TOLERANCE) // no to check 14
  {
    if (optVerbose) {
      output_message = std::string("✅ Successful convergence to a validated intersection at STEP 3, Check 14. \n" )+
                       "(check point is a close approach, accept pushed point)";
      std::vector<double> unknown_initial_guess = {101,101,101};
      PrintIntersectionOutcome(output_message, unknown_initial_guess, uvl_opt_pushed, 0, P0, direction); // inital guess lost
    }
    return std::make_tuple(uvl_opt_pushed, 0); // accept pushed point. Residual value not known here so assumed as 0
  }

  if (optVerbose) {
    std::cout <<"❌ logic error, case not identified after check 14 \n"
              <<"Returning LARGE_NUMBERs."<<std::endl;
  }
  return std::make_tuple(std::vector<double>(3, LARGE_NUMBER), LARGE_NUMBER);

}


std::tuple<std::vector<double>, double> G4GeorgeNurbs::PushLGuessForDTO(std::vector<double> uvl_opt_old,
                                                                        const G4ThreeVector& P0,
                                                                        const G4ThreeVector& direction,
                                                                        double line_length) const
{
  // STEP 2: push the initial guess by maxExtent and find optimal value from there. Should be the outermost intersection

  // declaring variables to be overwritten
  std::vector<double> uvl_pushed_opt; // optimised result after pushing
  double residual_pushed_opt; // optimised residual after pushing
  bool opt_success; // successful optimisation flag (unused in checks)

  std::vector<double> closest_knot_pushed;
  double R;
  std::string output_message;

  // Push l guess and optimise
  double l_pushed = uvl_opt_old[2] + maxExtent;
  G4ThreeVector P_pushed = P0 + l_pushed*direction; // new location from l_pushed
  std::tie(closest_knot_pushed, R) = ClosestKnot(P_pushed);
  std::vector<double> uvl_pushed_guess = {closest_knot_pushed[0], closest_knot_pushed[1] , l_pushed};
  std::tie(uvl_pushed_opt, residual_pushed_opt, opt_success) = LineIntersectionOpt(P0, direction, line_length, uvl_pushed_guess);


  /* Check 6: is this an intersection (residual = 0)? */
  if (residual_pushed_opt > CONVERGENCE_TOLERANCE) // no to check 6
  {
    /* Check 7: do the optimised u or v = the ends of their knot vectors?  */

    auto [uvl_pushed_guess_knot, reached_knot_bounds] = CheckKnotBounds(uvl_pushed_opt, uvl_pushed_guess);
    uvl_pushed_guess = uvl_pushed_guess_knot;
    uvl_pushed_guess[2] = uvl_pushed_opt[2]; // a close point was found, start with same l guess

    if (reached_knot_bounds) // yes to check 7
    {
      // redo optimisation with new initial knot guesses
      std::tie(uvl_pushed_opt, residual_pushed_opt, opt_success) = LineIntersectionOpt(P0, direction, line_length, uvl_pushed_guess);
    }
  }


  /* Check 8: did the first optimisation land on a point of close approach? */
  if (uvl_opt_old[2]>2*SURFACE_TOLERANCE) // yes to check 8
  {
    if (optVerbose) {
      output_message = std::string("✅ Successful convergence to a validated intersection at STEP 2, Check 8." );
      PrintIntersectionOutcome(output_message, uvl_pushed_guess, uvl_pushed_opt, residual_pushed_opt, P0, direction);
    }
    return std::make_tuple(uvl_pushed_opt, residual_pushed_opt); // accept pushed point
  }

  if (optVerbose) {
      output_message = std::string("🌀 value of pushed optimised value at the end of STEP 2." );
      PrintIntersectionOutcome(output_message, uvl_pushed_guess, uvl_pushed_opt, residual_pushed_opt, P0, direction);
    }

  // all logic streams now require a search for an intersection between initial optimisaion and the pushed optimisation.

  auto result = IntermediateIntersectionSearchForDTO(uvl_opt_old, uvl_pushed_opt, P0,direction, line_length);

  return result;
}

std::tuple<std::vector<double>, double> G4GeorgeNurbs::IntermediateIntersectionSearchForDTO(std::vector<double> uvl_opt_old,
                                                                             std::vector<double> uvl_opt_pushed,
                                                                             const G4ThreeVector& P0,
                                                                             const G4ThreeVector& direction,
                                                                             double line_length) const
{
  // STEP 3: Search for an intersection between the initial optimisation poing and the pushed optimisation point

  // declaring variables
  std::vector<double> uvl_check_opt; // optimised result after pushing
  double residual_check_opt; // optimised residual after pushing
  bool opt_success; // successful optimisation flag (unused in checks)
  std::string output_message;


  // optimise from a new point at some fraction of the way between l_opt_old and l_opt_pushed
  double l_check = uvl_opt_old[2] + 1.0/3.0*uvl_opt_pushed[2]; //
  G4ThreeVector P_check = P0 + l_check*direction; // new location from l_check
  auto [closest_knot_check, R] = ClosestKnot(P_check);
  std::vector<double> uvl_check_guess = {closest_knot_check[0], closest_knot_check[1] , l_check};
  std::tie(uvl_check_opt, residual_check_opt, opt_success) = LineIntersectionOpt(P0, direction, line_length, uvl_check_guess);


  /* Check 9: is this an intersecion? */
  if (residual_check_opt > SURFACE_TOLERANCE) // no to check 9
  {
    /* Check 10: could this non intersection be due to knot boundaries? */
    auto [new_uvl_check_guess, reached_knot_bounds] = CheckKnotBounds(uvl_check_opt, uvl_check_guess);
    uvl_check_guess = new_uvl_check_guess;

    if (reached_knot_bounds) // yes to check 10
    {
      // redo optimisation with new initial knot guesses
      std::tie(uvl_check_opt, residual_check_opt, opt_success) = LineIntersectionOpt(P0, direction, line_length, uvl_check_guess);
    }

   /* Check 11: with resolved knot vector issue, is this now an intersection (residual = 0)? */
    if (residual_check_opt>SURFACE_TOLERANCE) // no to check 11
    {
      if (optVerbose) {
      output_message = std::string("✅ Successful convergence to a validated intersection at STEP 3, Check 11." );
      std::vector<double> unknown_initial_guess = {101,101,101};
      PrintIntersectionOutcome(output_message, unknown_initial_guess, uvl_opt_pushed, 0, P0, direction); // inital guess lost
      }
      return std::make_tuple(uvl_opt_pushed, 0); // accept pushed point. Residual value not known here so assumed as 0
    }
  }
  // only valid intersections remaining after this
  if (optVerbose) {
    output_message = "❓ Intermediate intersection found at STEP 3.";
    PrintIntersectionOutcome(output_message, uvl_check_guess, uvl_check_opt, residual_check_opt, P0, direction);
  }
  /* Check 12: is this a new intersection (different to l_opt_old and l_opt_push? */
  if (std::abs(uvl_opt_old[2]-uvl_check_opt[2])>3*SURFACE_TOLERANCE &&
      std::abs((uvl_opt_pushed[2]-uvl_check_opt[2]))>3*SURFACE_TOLERANCE) // no to check 12
  {
    if (optVerbose) {
      output_message = std::string("✅ Successful convergence to a validated intersection at STEP 3, Check 12 (yes).");
      PrintIntersectionOutcome(output_message, uvl_check_guess, uvl_check_opt, residual_check_opt, P0, direction);
    }
    return std::make_tuple(uvl_check_opt, residual_check_opt); // accept check point



  }

  else // yes to check 12
  {
    if (optVerbose) {
      output_message = std::string("✅ Successful convergence to a validated intersection at STEP 3, Check 12 (no)." );
      std::vector<double> unknown_initial_guess = {101,101,101};
      PrintIntersectionOutcome(output_message, unknown_initial_guess, uvl_opt_pushed, 0, P0, direction); // inital guess lost
    }
    return std::make_tuple(uvl_opt_pushed, 0); // accept pushed point. Residual value not known here so assumed as 0
  }

}





// Static functions  used in optimisation =============================================================


double G4GeorgeNurbs::ResidualToSurfacePoint(const std::vector<double>& uv, std::vector<double>& grad, void* data)
{
  // definition of static function used by NLopt to compute the distance (residual)
  // between a surface point and a fixed target point.

  (void)grad; // Suppresses warning on build, gradient not used here.


  // Unpacking context (contains the target point and pointer to the NURBS surface)
  const auto* context = static_cast<G4GeorgeNurbs::ClostestPointContext*>(data);
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
  context->call_count++;

  double u = uvl[0]; // u and v parametrises a point on the surface
  double v = uvl[1];
  double l = uvl[2]; // l parametrises a point on the line




  G4ThreeVector P_nurbs = context->nurbs->SurfacePoint(u, v);
  G4ThreeVector P_line = context->P0 + l * context->direction;

  (void)grad; // to supress warning on build


  // debug
//  std::cout<<"u,v = "<<u<<" , "<<v<<"\n" <<"l = "<<l<<"\n";
//  std::cout<<"r = "<<(P_line - P_nurbs).mag()<<"\n";
//  std::cout<<"-----------"<<"\n";


  return (P_line - P_nurbs).mag(); // return distance between those two points
}











// Overridden Base class functions ====================================================================

EInside G4GeorgeNurbs::Inside(const G4ThreeVector& p) const
{

  // finding the closest surface point to p
  auto [uv_closest_pt, R_closest_pt] = ClosestPointParams(p); // u,v of closest point (can maybe just use ClosestKnot)
  G4ThreeVector closest_surf_pt = SurfacePoint(uv_closest_pt[0], uv_closest_pt[1]); // [x,y,z] of closest pt


  // finding the residual from p --> closest_surf_pt
  G4ThreeVector residual = closest_surf_pt - p;

  if (residual.mag() < SURFACE_TOLERANCE)
  {
    // G4cout<<" Inside: "<<p<< " kSurface "<<G4endl;
    // if the residual is less than SURFACE_TOLERANCE, then p must be on the surface.
    return kSurface;
  }

  // finding normal at the surface point and taking dot product with the residual
  G4ThreeVector normal_at_surf_pt = SurfaceNormal(uv_closest_pt[0], uv_closest_pt[1]);
  G4double residual_dot_normal = residual.dot(normal_at_surf_pt);

  // inside/outside can be determined from the sign of the dot product
  if (residual_dot_normal>0){
    // G4cout<<" Inside: "<<p<< " kInside "<<G4endl;
    return kInside;
  }
  // G4cout<<" Inside: "<<p<< " kOutside "<<G4endl;
  return kOutside;
}


G4ThreeVector G4GeorgeNurbs::SurfaceNormal(const G4ThreeVector& p) const
{
  // finding the u,v parameters of the surface point with ClosestPointParams function.
  // residual_from_p should be 0.
  auto [uv_from_p, residual_from_p] = ClosestPointParams(p);
  G4ThreeVector n = SurfaceNormal(uv_from_p[0], uv_from_p[1]);

  // G4cout<<" SurfaceNormal: "<<n<<G4endl;
  return n;
}



G4ThreeVector G4GeorgeNurbs::SurfaceNormal(const double u, const double v) const
{
  // overloaded surface mornal function taking parameters u and v
  std::vector<G4ThreeVector> tangents = SurfaceDerivatives(u,v);

  G4ThreeVector& tangent_u = tangents[0];
  G4ThreeVector& tangent_v = tangents[1];

  G4ThreeVector normal = tangent_u.cross(tangent_v);
  normal = normal.unit();
  return normal;
}



G4double G4GeorgeNurbs::DistanceToIn(const G4ThreeVector& p) const
{
  //G4cout<<" DistanceToIn(p): "<<p<<" (Default 10)"<<G4endl;
  //return 10;



  // first check if point is already inside, return 0 if so
  EInside inside_status = Inside(p);
  if (inside_status == kInside ) {
    return 0;
  }

  // find the closest surface point using ClosestPointParams function
  auto [uv, residual_from_p] = ClosestPointParams(p);

  // return the minimised distance
  return residual_from_p;
}

G4double G4GeorgeNurbs::DistanceToIn(const G4ThreeVector& p0, const G4ThreeVector& v) const
{
  EInside inside_status = Inside(p0);
  if (inside_status == kInside) {
    return 0;
  }

  auto [uvl_opt, minimised_residual] = LineIntersectionParams(p0,v,LARGE_NUMBER, 0);

  //G4cout<<" DistanceToIn(p,v) "<<p0<< " " <<v<<" "<< uvl_opt[2]<<G4endl;
  return uvl_opt[2];
}

G4double G4GeorgeNurbs::DistanceToOut(const G4ThreeVector& p) const
{
  //G4cout<<" DistanceToOut(p) "<<p<<" (Default 10)"<<G4endl;
  //return 10;


  // first check if point is already outside, return 0 if so

  EInside inside_status = Inside(p);
  if (inside_status == kOutside) {
    return 0;
  }

  // find the closest point on the surface using ClosestPointParams function
  auto [uv, residual_from_p] = ClosestPointParams(p);

  // return the minimised distance
  return residual_from_p;

}

G4double G4GeorgeNurbs::DistanceToOut( const G4ThreeVector& p,const G4ThreeVector& v,
                                       const G4bool calcNorm,
                                        G4bool* validNorm,
                                        G4ThreeVector* n ) const
{


  // check its not already outside
  EInside inside_status = Inside(p);
  if (inside_status == kOutside) {
    return 0;
  }

  // find intersection of surface point (params u,v) and the line (param l)
  auto [uvl_opt, residual_opt] = LineIntersectionParams(p, v.unit(), maxExtent*2, 1);
  //G4cout<<" DistanceToOut(p,v) "<<p<< " " <<v<< " "<< uvl_opt[2]<<G4endl;


  /* Implementing this optional part completely breaks the programme for some reason, even though the method is velid*/

  // Handle normal computation
  //  if (calcNorm && n && validNorm) {
  //
  //
  //    G4ThreeVector normal = SurfaceNormal(uvl_opt[0], uvl_opt[1]);
  //
  //
  //
  //    if (normal.mag2() > 0) {  // Validate it's not a zero vector
  //      *n = normal.unit();
  //      *validNorm = true;
  //    } else {
  //      *validNorm = false;
  //    }
  //    G4cout<<" (Norm requested in DistanceToOut(p,v), norm: "<<*n<<",  validNorm: "<<*validNorm<< G4endl;
  //
  //  }


  return uvl_opt[2];
}

G4double G4GeorgeNurbs::GetCubicVolume()
{
  return 10;
}

G4double G4GeorgeNurbs::GetSurfaceArea()
{
  return 10;
}

G4ThreeVector G4GeorgeNurbs::GetPointOnSurface() const
{

  G4double rand_u = G4UniformRand();
  G4double rand_v = G4UniformRand();
  G4ThreeVector surf_point = SurfacePoint(rand_u, rand_v);
  //G4cout<<" GetPointOnSurface(): p("<<rand_u<< ", " <<rand_v<< ") = "<< surf_point<<G4endl;

  return surf_point;
}

G4VisExtent G4GeorgeNurbs::GetExtent() const
{
  //G4cout<<" GetExtent() ( bminCached: "<<bminCached<< " , " <<bmaxCached<< " )"<<G4endl;

  return G4VisExtent(
    bminCached.x(), bmaxCached.x(),
    bminCached.y(), bmaxCached.y(),
    bminCached.z(), bmaxCached.z()
  );
}



G4bool G4GeorgeNurbs::CalculateExtent( const EAxis pAxis,
                                  const G4VoxelLimits& pVoxelLimit,
                                  const G4AffineTransform& pTransform,
                                        G4double& pMin, G4double& pMax ) const
{
  G4ThreeVector bmin, bmax;

  bmin = bminCached;
  bmax = bmaxCached;

  // Get bounding box
  //BoundingLimits(bmin,bmax);


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


void G4GeorgeNurbs::SetBoundingLimits()
{
  // Safety check: empty control point grid
  if (controlPts.empty() || controlPts[0].empty()) {
    bminCached = G4ThreeVector(0.0, 0.0, 0.0);
    bmaxCached = G4ThreeVector(0.0, 0.0, 0.0);
    maxExtent = 0.0;
    boundsCached = true;
    return;
  }

  // Start from first control point
  bminCached = controlPts[0][0];
  bmaxCached = controlPts[0][0];


  //iterate through control points to find the bounding box
  for (const auto& row : controlPts) {
    for (const auto& pt : row) {
      if (pt.x() < bminCached.x()) bminCached.setX(pt.x());
      if (pt.x() > bmaxCached.x()) bmaxCached.setX(pt.x());

      if (pt.y() < bminCached.y()) bminCached.setY(pt.y());
      if (pt.y() > bmaxCached.y()) bmaxCached.setY(pt.y());

      if (pt.z() < bminCached.z()) bminCached.setZ(pt.z());
      if (pt.z() > bmaxCached.z()) bmaxCached.setZ(pt.z());
    }
  }

  maxExtent = (bmaxCached - bminCached).mag();
  boundsCached = true;
}

std::vector<G4ThreeVector> G4GeorgeNurbs::GetBounds() const {return {bminCached, bmaxCached};}
G4double G4GeorgeNurbs::GetMaxExtent() const {return maxExtent;};

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


void G4GeorgeNurbs::EnableOptVerbose(){optVerbose = true;}
void G4GeorgeNurbs::DisableOptVerbose(){optVerbose = false;}

