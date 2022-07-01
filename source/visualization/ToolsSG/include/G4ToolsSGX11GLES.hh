//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Guy Barrand 12th March 2021

#ifndef G4TOOLSSGX11GLES_HH
#define G4TOOLSSGX11GLES_HH

#include "G4VGraphicsSystem.hh"

namespace toolx {namespace X11 {class session;}}

class G4ToolsSGX11GLES: public G4VGraphicsSystem {
  typedef G4VGraphicsSystem parent;
public:
  G4ToolsSGX11GLES();
  virtual ~G4ToolsSGX11GLES();
protected:  
  G4ToolsSGX11GLES(const G4ToolsSGX11GLES& a_from):parent(a_from){}
  G4ToolsSGX11GLES& operator=(const G4ToolsSGX11GLES&) {return *this;}
public:
  G4VSceneHandler* CreateSceneHandler(const G4String& name = "");
  G4VViewer* CreateViewer (G4VSceneHandler&, const G4String& name = "");
  G4bool IsUISessionCompatible () const;
protected:  
  void Initialise();
protected:
  toolx::X11::session* fSGSession;
};

#endif
