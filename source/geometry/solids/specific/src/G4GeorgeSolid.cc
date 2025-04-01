//class defined a sphere

#include "G4GeorgeSolid.hh"

#include "G4BoundingEnvelope.hh"
#include "G4QuickRand.hh"
#include "G4VisExtent.hh"
#include "G4VGraphicsScene.hh"
#include "G4AffineTransform.hh"

#include "CLHEP/Units/SystemOfUnits.h"

// all includes used in sphere
#include "G4GeomTools.hh"
#include "G4GeometryTolerance.hh"
#include "G4VPVParameterisation.hh"
#include "G4VoxelLimits.hh"

#include "meshdefs.hh"
#include <G4Box.hh>

G4GeorgeSolid::G4GeorgeSolid(const G4String& name, const G4ThreeVector& centreIn,const double& radiusIn):
    G4VSolid(name), centre(centreIn), radius(radiusIn)
{
  if (radius <= 0.0) {
    G4Exception("G4GeorgeSolid", "InvalidRadius",
                FatalException, "Radius must be positive.");
  }
}

G4GeorgeSolid::G4GeorgeSolid(const G4String& name, const double& radiusIn):
    G4VSolid(name), centre(G4ThreeVector(0.0, 0.0,0.0)), radius(radiusIn)
{
  if (radius <= 0.0) {
    G4Exception("G4GeorgeSolid", "InvalidRadius",
                FatalException, "Radius must be positive.");
  }
}

G4GeorgeSolid::G4GeorgeSolid(const G4String& name):
    G4VSolid(name), centre(G4ThreeVector(0.0, 0.0,0.0)), radius(1.0){}

G4GeorgeSolid::~G4GeorgeSolid() = default;



G4ThreeVector G4GeorgeSolid::getCentre() const {return centre;}
G4double G4GeorgeSolid::getRadius() const {return radius;}

void G4GeorgeSolid::setCentre(const G4ThreeVector& centreIn){centre = centreIn;}
void G4GeorgeSolid::setRadius(const G4double& radiusIn)
{
  if (radiusIn <= 0.0) {
    G4Exception("G4GeorgeSolid", "InvalidRadius",
                FatalException, "Radius must be positive.");
  }
   radius= radiusIn;
}



EInside G4GeorgeSolid::Inside(const G4ThreeVector& p) const
{


  if ((p-centre).mag()<radius) {
    return kInside;
  }
  else if ((p-centre).mag()==radius) {
    return kSurface;
  }
  else {
    return kOutside;
  }
}

G4ThreeVector G4GeorgeSolid::SurfaceNormal(const G4ThreeVector& p) const
{
  return (p - centre).unit();
}

G4double G4GeorgeSolid::DistanceToIn(const G4ThreeVector& p) const
{
  // G4Box *g = nullptr;
  // std::cout << g->GetName()<< std::endl;

  G4ThreeVector direction = (p-centre).unit(); // direction to point from centre
  G4ThreeVector p_on_surface = direction*radius;
  G4double distance = (p-p_on_surface).mag();
  if (distance <= 0)
  {
    return 0;
  }
  return distance;
}

G4double G4GeorgeSolid::DistanceToIn(const G4ThreeVector& p0, const G4ThreeVector& v) const
{
  /* found by subbing the equation for a line p = p0 + l*n
   * into the equation for a sphere |p - c|^2 = R^2
   * which has solutions l = -(n.d) ± sqrt((n.d)^2 -(d^2-R^2)) where d = p0 - c  */

  G4ThreeVector n = v.unit();
  G4ThreeVector d = p0-centre;

  G4double B = n.dot(d);                  // A and B defined as in quadratic formula
  G4double A = d.mag2() - radius*radius;
  G4double discriminant = B*B - A;

  if (discriminant < 0) {return kInfinity;} //no intersection since no real solutions

  G4double l1 = -B - std::sqrt(discriminant); //solution closer to start point of line -->--*(  )
  G4double l2 = -B + std::sqrt(discriminant); //solution closer to end point of line -->---(-*)

  if (l1>0) {return l1;} // distance from surface
  else if (l2>0) {return 0;} // point already inside sphere (->-) ie. one neg one pos



  return kInfinity; // particle moving away from sphere - no intersection <---( )
}

G4double G4GeorgeSolid::DistanceToOut(const G4ThreeVector& p) const
{
  G4ThreeVector direction = (p-centre).unit(); // direction to point from centre
  G4ThreeVector p_on_surface = direction*radius;
  G4double distance = (p-p_on_surface).mag();
  return distance;
}

G4double G4GeorgeSolid::DistanceToOut( const G4ThreeVector& p,const G4ThreeVector& v,
                                       const G4bool calcNorm,
                                        G4bool* validNorm,
                                        G4ThreeVector* n ) const
{
  // same maths as for DistanceToIn function
  G4ThreeVector norm = v.unit();
  G4ThreeVector d = p-centre;

  G4double B = norm.dot(d);                  // A and B defined as in quadratic formula
  G4double A = d.mag2() - radius*radius;
  G4double discriminant = B*B - A;

  //this func should only be called if point is already inside
  //it will have a positive and negative solution, need to take the positive solution as this is
  //in the direction of the velocity vector. This should always be the solution with + (ie. l1)

  G4double l1 = -B + std::sqrt(discriminant);
  // G4double l2 = -B - std::sqrt(discriminant);
  if (l1>0) {return l1;}
  //else if (l2>0) {return l2;}

  return kInfinity; // shouldn't occur in a closed surface
}

G4double G4GeorgeSolid::GetCubicVolume()
{
  G4double vol = 4.0/3.0*CLHEP::pi*radius*radius*radius;
  return vol;
}

G4double G4GeorgeSolid::GetSurfaceArea()
{
  G4double SA = 4.0*CLHEP::pi*radius*radius;
  return SA;
}

G4ThreeVector G4GeorgeSolid::GetPointOnSurface() const
{
  //defining random direction
  G4double nx = 2 * G4QuickRand() - 1; // random number between -1 and 1
  G4double ny = 2 * G4QuickRand() - 1;
  G4double nz = 2 * G4QuickRand() - 1;
  G4ThreeVector randomDirection = G4ThreeVector(nx, ny, nz).unit();

  //creating a point on the surface
  G4ThreeVector point_on_surface = randomDirection*radius + centre;
  return point_on_surface;
}


G4VisExtent G4GeorgeSolid::GetExtent() const
{
  return { centre.x() - radius, centre.x() + radius,
         centre.y() - radius, centre.y() + radius,
         centre.z() - radius, centre.z() + radius };
}


G4bool G4GeorgeSolid::CalculateExtent( const EAxis pAxis,
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

void G4GeorgeSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
  scene.AddSolid (*this);
}

std::ostream& G4GeorgeSolid::StreamInfo( std::ostream& os ) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4GeorgeSolid\n"
     << " Parameters: \n"
     << "    radius: " << radius << " cm \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

G4GeometryType G4GeorgeSolid::GetEntityType() const
{
  return {"G4GeorgeSolid"};
}


void G4GeorgeSolid::print_info() const
{
  std::cout << "G4GeorgeSolid: " << GetName() << "\n";
  std::cout << "Radius: " << radius << "\n";
  std::cout << "Centre: " << centre << "\n";
}
