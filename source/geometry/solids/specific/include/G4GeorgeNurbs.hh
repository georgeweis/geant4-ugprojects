// nurbs class header



#ifndef G4GEORGENURBS_HH
#define G4GEORGENURBS_HH

#include <iostream>
#include <vector>

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"


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

  G4double BasisFunction(G4int i, G4int k, G4double t, const std::vector<G4double>& knotVector) const;
  G4ThreeVector SurfacePoint(G4double u, G4double v) const;


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








};



#endif // G4GEORGENURBS_HH
