/random/setSavingFlag 1
/run/verbose 1
/event/verbose 0 
/tracking/verbose 0
#
/run/initialize
/gun/particle pi-
/gun/energy 50 GeV
#
#---
#
/mydet/material_tracker G4_POLYSTYRENE
/mydet/inner_radius_tracker  5.0 cm
/mydet/outer_radius_tracker 25.0 cm
#
/mydet/material_emCalo G4_PbWO4
/mydet/inner_radius_emCalo 30.0 cm
/mydet/outer_radius_emCalo 60.0 cm
#
/mydet/material_hadCalo G4_Cu
/mydet/inner_radius_hadCalo  70.0 cm
/mydet/outer_radius_hadCalo 170.0 cm
#
/mydet/update
/run/beamOn 100
#
#---
#
/mydet/material_tracker G4_lAr
/mydet/inner_radius_tracker  5.0 cm
/mydet/outer_radius_tracker 25.0 cm
#
/mydet/material_emCalo G4_Pb
/mydet/inner_radius_emCalo 30.0 cm
/mydet/outer_radius_emCalo 60.0 cm
#
/mydet/material_hadCalo G4_Fe
/mydet/inner_radius_hadCalo  70.0 cm
/mydet/outer_radius_hadCalo 170.0 cm
#
/mydet/update
/run/beamOn 100
#
#---
#
/mydet/material_tracker G4_Si
/mydet/inner_radius_tracker  5.0 cm
/mydet/outer_radius_tracker 25.0 cm
#
/mydet/material_emCalo G4_W
/mydet/inner_radius_emCalo 30.0 cm
/mydet/outer_radius_emCalo 60.0 cm
#
/mydet/material_hadCalo G4_W
/mydet/inner_radius_hadCalo  70.0 cm
/mydet/outer_radius_hadCalo 170.0 cm
#
/mydet/update
/run/beamOn 100
