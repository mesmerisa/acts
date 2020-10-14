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
    ActsExampleCKFTracksBGVnew
    RecCKFTracksBGVnew.cpp)
target_link_libraries(
    ActsExampleCKFTracksBGVnew
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
  ActsExampleCKFTracksBGVnew
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

