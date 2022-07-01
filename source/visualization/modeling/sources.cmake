# - G4modeling module build definition

# Define the Geant4 Module.
geant4_add_module(G4modeling
  PUBLIC_HEADERS
    G4ArrowModel.hh
    G4AttFilterUtils.hh
    G4AttValueFilterT.hh
    G4AttributeFilterT.hh
    G4AxesModel.hh
    G4BoundingExtentScene.hh
    G4CallbackModel.hh
    G4DigiFilterFactories.hh
    G4DigiModel.hh
    G4ElectricFieldModel.hh
    G4GPSModel.hh
    G4HitFilterFactories.hh
    G4HitsModel.hh
    G4LogicalVolumeModel.hh
    G4MagneticFieldModel.hh
    G4PlotterModel.hh
    G4Mesh.hh
    G4ModelApplyCommandsT.hh
    G4ModelColourMap.hh
    G4ModelCommandUtils.hh
    G4ModelCommandsT.hh
    G4ModelCompoundCommandsT.hh
    G4ModelingParameters.hh
    G4ModelingParameters.icc
    G4NullModel.hh
    G4PSHitsModel.hh
    G4PhysicalVolumeMassScene.hh
    G4PhysicalVolumeModel.hh
    G4PhysicalVolumesSearchScene.hh
    G4PseudoScene.hh
    G4TextModel.hh
    G4TouchablePropertiesScene.hh
    G4TouchableUtils.hh
    G4TrajectoriesModel.hh
    G4TrajectoryChargeFilter.hh
    G4TrajectoryDrawByAttribute.hh
    G4TrajectoryDrawByCharge.hh
    G4TrajectoryDrawByOriginVolume.hh
    G4TrajectoryDrawByParticleID.hh
    G4TrajectoryDrawByEncounteredVolume.hh
    G4TrajectoryDrawerUtils.hh
    G4TrajectoryFilterFactories.hh
    G4TrajectoryGenericDrawer.hh
    G4TrajectoryModelFactories.hh
    G4TrajectoryOriginVolumeFilter.hh
    G4TrajectoryParticleFilter.hh
    G4TrajectoryEncounteredVolumeFilter.hh
    G4VAttValueFilter.hh
    G4VFieldModel.hh
    G4VModel.hh
    G4VModel.icc
    G4VModelCommand.hh
    G4VModelFactory.hh
    G4VTrajectoryModel.hh
    G4VisTrajContext.hh
    G4VisTrajContext.icc
  SOURCES
    G4ArrowModel.cc
    G4AttFilterUtils.cc
    G4AxesModel.cc
    G4BoundingExtentScene.cc
    G4DigiFilterFactories.cc
    G4DigiModel.cc
    G4ElectricFieldModel.cc
    G4GPSModel.cc
    G4HitFilterFactories.cc
    G4HitsModel.cc
    G4LogicalVolumeModel.cc
    G4MagneticFieldModel.cc
    G4PlotterModel.cc
    G4Mesh.cc
    G4ModelingParameters.cc
    G4NullModel.cc
    G4PSHitsModel.cc
    G4PhysicalVolumeMassScene.cc
    G4PhysicalVolumeModel.cc
    G4PhysicalVolumesSearchScene.cc
    G4PseudoScene.cc
    G4TextModel.cc
    G4TouchablePropertiesScene.cc
    G4TouchableUtils.cc
    G4TrajectoriesModel.cc
    G4TrajectoryChargeFilter.cc
    G4TrajectoryDrawByAttribute.cc
    G4TrajectoryDrawByCharge.cc
    G4TrajectoryDrawByOriginVolume.cc
    G4TrajectoryDrawByParticleID.cc
    G4TrajectoryDrawByEncounteredVolume.cc
    G4TrajectoryDrawerUtils.cc
    G4TrajectoryFilterFactories.cc
    G4TrajectoryGenericDrawer.cc
    G4TrajectoryModelFactories.cc
    G4TrajectoryOriginVolumeFilter.cc
    G4TrajectoryParticleFilter.cc
    G4TrajectoryEncounteredVolumeFilter.cc
    G4VFieldModel.cc
    G4VModel.cc
    G4VTrajectoryModel.cc
    G4VisTrajContext.cc)

geant4_module_link_libraries(G4modeling
  PUBLIC
    G4graphics_reps
    G4csg
    G4geometrymng
    G4hepgeometry
    G4specsolids
    G4hepnumerics
    G4intercoms
    G4globman
    G4digits
    G4hits
    G4tracking
  PRIVATE
    G4geomBoolean
    G4navigation
    G4event
    G4magneticfield
    G4materials
    G4volumes
    G4run
    G4detutils
    G4detector)
