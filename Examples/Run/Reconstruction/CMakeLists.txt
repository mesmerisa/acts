add_executable(
  ActsExampleTruthTracks
  RecTruthTracks.cpp)

target_link_libraries(
  ActsExampleTruthTracks
  PRIVATE
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesDigitization
    ActsExamplesTrackFitting
    ActsExamplesDetectorGeneric
    ActsExamplesMagneticField
    ActsExamplesTruthTracking
    ActsExamplesIoCsv
    ActsExamplesIoPerformance)

install(
  TARGETS ActsExampleTruthTracks
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(
    ActsExampleCKFTracks
    RecCKFTracks.cpp)
target_link_libraries(
    ActsExampleCKFTracks
    PRIVATE
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesDigitization
    ActsExamplesTrackFinding
    ActsExamplesDetectorGeneric
    ActsExamplesMagneticField
    ActsExamplesIoCsv
    ActsExamplesIoPerformance)

install(
  TARGETS
  ActsExampleCKFTracks
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

