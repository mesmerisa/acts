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
  ActsExampleMaterialMappingTGeo
  TGeoMaterialMapping.cpp)
target_link_libraries(
  ActsExampleMaterialMappingTGeo
  PRIVATE ${_common_libraries} ActsExamplesMaterialMapping ActsExamplesDetectorTGeo)

add_executable(
  ActsExampleMaterialValidationTGeo
  TGeoMaterialValidation.cpp)
target_link_libraries(
  ActsExampleMaterialValidationTGeo
  PRIVATE ${_common_libraries} ActsExamplesPropagation ActsExamplesDetectorTGeo)

install(
  TARGETS
    ActsExampleMaterialValidationGeneric
    ActsExampleMaterialMappingGeneric
    ActsExampleMaterialMappingTGeo
    ActsExampleMaterialValidationTGeo    
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
add_subdirectory_if(TGeo ACTS_BUILD_PLUGIN_TGEO)
