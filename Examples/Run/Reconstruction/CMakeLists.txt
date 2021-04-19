# shared code 
add_subdirectory(Common)

# Truth fitting for Generic detector 
add_executable(
  ActsExampleTruthTracksGeneric
  GenericRecTruthTracks.cpp)
target_link_libraries(
  ActsExampleTruthTracksGeneric
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorGeneric)

# Truth fitting for Telescope detector
add_executable(
  ActsExampleTruthTracksTelescope
  TelescopeRecTruthTracks.cpp)
target_link_libraries(
  ActsExampleTruthTracksTelescope
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorTelescope)
    
# Truth fitting for TGeo detector
add_executable(
  ActsExampleTruthTracksTGeo
  TGeoRecTruthTracks.cpp)
target_link_libraries(
  ActsExampleTruthTracksTGeo
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorTGeo)   
    
    
# Truth fitting for TGeo detector
add_executable(
  ActsExampleCombFinderTracksTGeo
  TGeoRecCombFinderTracks.cpp)
target_link_libraries(
  ActsExampleCombFinderTracksTGeo
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorTGeo)       
     

# CKF reco for Generic detector 
add_executable(
  ActsExampleCKFTracksGeneric
  GenericRecCKFTracks.cpp)
target_link_libraries(
  ActsExampleCKFTracksGeneric
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorGeneric)

# CKF reco for Telescope detector 
add_executable(
  ActsExampleCKFTracksTelescope
  TelescopeRecCKFTracks.cpp)
target_link_libraries(
  ActsExampleCKFTracksTelescope
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorTelescope)

# Seeding for Generic detector
add_executable(
  ActsExampleSeedingGeneric
  GenericSeedingExample.cpp)
target_link_libraries(
  ActsExampleSeedingGeneric
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorGeneric
)

install(
  TARGETS
    ActsExampleTruthTracksGeneric
    ActsExampleTruthTracksTelescope
    ActsExampleTruthTracksTGeo
    ActsExampleCKFTracksGeneric
    ActsExampleCKFTracksTelescope
    ActsExampleSeedingGeneric 
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
