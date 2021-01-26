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

# Truth fitting for Telescope detector
add_executable(
  ActsExampleTruthTracksTelescopeBGV
  TelescopeBGVRecTruthTracks.cpp)
target_link_libraries(
  ActsExampleTruthTracksTelescopeBGV
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

install(
  TARGETS
    ActsExampleTruthTracksGeneric
    ActsExampleTruthTracksTelescope
    ActsExampleTruthTracksTelescopeBGV
    ActsExampleCKFTracksGeneric
    ActsExampleCKFTracksTelescope
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(
  ActsExampleTruthTracksBGV
  RecTruthTracksBGV.cpp)

target_link_libraries(
  ActsExampleTruthTracksBGV
  PRIVATE
ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesDigitization
    ActsExamplesTrackFitting
    ActsExamplesDetectorTGeo
    ActsExamplesMagneticField
    ActsExamplesTruthTracking
    ActsExamplesIoCsv
    ActsExamplesIoPerformance)

install(
  TARGETS ActsExampleTruthTracksBGV
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  
add_executable(
  ActsExampleCombTrackFinderBGV
  RecCombTrackFinderBGV.cpp)

target_link_libraries(
  ActsExampleCombTrackFinderBGV
  PRIVATE
ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesDigitization
    ActsExamplesTrackFitting
    ActsExamplesDetectorTGeo
    ActsExamplesMagneticField
    ActsExamplesTruthTracking
    ActsExamplesIoCsv
    ActsExamplesIoPerformance)

install(
  TARGETS ActsExampleCombTrackFinderBGV
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})  
  
add_executable(
  ActsExampleTruthTracksVertexBGV
  RecTruthTracksVertexBGV.cpp)

target_link_libraries(
  ActsExampleTruthTracksVertexBGV
  PRIVATE
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesDigitization
    ActsExamplesTrackFitting
    ActsExamplesDetectorTGeo
    ActsExamplesMagneticField
    ActsExamplesTruthTracking
    ActsExamplesIoCsv
    ActsExamplesVertexing
    ActsExamplesIoPerformance)

install(
  TARGETS ActsExampleTruthTracksVertexBGV
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})  
  

add_executable(
    ActsExampleCKFTracksBGV
    RecCKFTracksBGV.cpp)
target_link_libraries(
    ActsExampleCKFTracksBGV
    PRIVATE
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesDigitization
    ActsExamplesTrackFinding
    ActsExamplesDetectorTGeo
    ActsExamplesMagneticField
    ActsExamplesIoCsv
    ActsExamplesIoPerformance)

install(
  TARGETS
  ActsExampleCKFTracksBGV
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)

