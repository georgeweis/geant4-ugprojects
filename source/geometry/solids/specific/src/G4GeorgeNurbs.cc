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


#ifdef GEANT4_USE_NLOPT
#include <nlopt.hpp> // only includes if GEANT4_USE_NLOPT = ON
#endif



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
  try {
    // Create a dummy 2D optimizer
    nlopt::opt test_opt(nlopt::LN_NELDERMEAD, 2);
    std::cout << "NLopt is working! Algorithm: "
              << test_opt.get_algorithm() << std::endl;
  }
  catch (const std::exception& e) {
    std::cerr << "NLopt error: " << e.what() << std::endl;
  }
}



G4double G4GeorgeNurbs::BasisFunction(G4int i, G4int k, G4double t, const std::vector<G4double>& knotVector) const
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

  const double LARGE_NUMBER = 1e6;

  if (denominator == 0.0) {
    return G4ThreeVector(LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER);
  } else {
    return numerator / denominator;
  }
}














EInside G4GeorgeNurbs::Inside(const G4ThreeVector& p) const
{
  return kInside;
}


G4ThreeVector G4GeorgeNurbs::SurfaceNormal(const G4ThreeVector& p) const
{
  G4ThreeVector n = G4ThreeVector(0, 0, 0);
  return n;
}

G4double G4GeorgeNurbs::DistanceToIn(const G4ThreeVector& p) const
{
  return 0.0;
}

G4double G4GeorgeNurbs::DistanceToIn(const G4ThreeVector& p0, const G4ThreeVector& v) const
{
  return 0.0;
}

G4double G4GeorgeNurbs::DistanceToOut(const G4ThreeVector& p) const
{
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
  int expected_knot_length_U = controlPts.size() + degreeU + 1;
  int expected_knot_length_V = controlPts[0].size() + degreeV + 1;

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