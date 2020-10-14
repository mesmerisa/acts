set(_common_libraries
  ActsCore
  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesIoRoot
  ActsExamplesIoJson
  ActsExamplesMagneticField
  Boost::program_options)

add_executable(
  ActsExampleMaterialValidationGeneric
  GenericMaterialValidation.cpp)
target_link_libraries(
  ActsExampleMaterialValidationGeneric
  PRIVATE ${_common_libraries} ActsExamplesPropagation ActsExamplesDetectorGeneric)

add_executable(
  ActsExampleMaterialMappingGeneric
  GenericMaterialMapping.cpp)
target_link_libraries(
  ActsExampleMaterialMappingGeneric
  PRIVATE ${_common_libraries} ActsExamplesMaterialMapping ActsExamplesDetectorGeneric)

add_executable(
  ActsExampleMaterialMappingTGeoBGV
  TGeoMaterialMappingBGV.cpp)
target_link_libraries(
  ActsExampleMaterialMappingTGeoBGV
  PRIVATE ${_common_libraries} ActsExamplesMaterialMapping ActsExamplesDetectorTGeo)

add_executable(
  ActsExampleMaterialValidationTGeoBGV
  TGeoMaterialValidationBGV.cpp)
target_link_libraries(
  ActsExampleMaterialValidationTGeoBGV
  PRIVATE ${_common_libraries} ActsExamplesPropagation ActsExamplesDetectorTGeo)

install(
  TARGETS
    ActsExampleMaterialValidationGeneric
    ActsExampleMaterialMappingGeneric
    ActsExampleMaterialMappingTGeoBGV
    ActsExampleMaterialValidationTGeoBGV    
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
add_subdirectory_if(TGeo ACTS_BUILD_PLUGIN_TGEO)
