#ifndef G4GEORGESOLID_H
#define G4GEORGESOLID_H


#include "G4VSolid.hh"

class G4GeorgeSolid : public G4VSolid
{
  private:
  G4ThreeVector centre;
  G4double radius;

  public:
    G4GeorgeSolid(const G4String& name, const G4ThreeVector& centreIn,const G4double& radiusIn);
    G4GeorgeSolid(const G4String& name);
    ~G4GeorgeSolid() override;

    G4ThreeVector getCentre() const;
    G4double getRadius() const;

    void setCentre(const G4ThreeVector& centreIn);
    void setRadius(const G4double& radiusIn);


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

};


#endif