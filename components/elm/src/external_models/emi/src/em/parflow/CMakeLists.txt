set(EMI_EM_PARFLOW_SOURCES
  ExternalModelPARFLOW.F90
  elm_parflow.F90
  clm_parflow_interface_data.F90  
  inputread.F90
  mapping.F90
  parflow_model.F90 
)

include_directories(${CMAKE_BINARY_DIR}/elm_stub/shr)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/utils)
include_directories(${CMAKE_BINARY_DIR}/constants)
include_directories(${CMAKE_BINARY_DIR}/emi_data_types)
include_directories(${CMAKE_BINARY_DIR}/emi_data_definition)
include_directories(${CMAKE_BINARY_DIR}/em/base)

include(add_emi_library)
add_emi_library(emi_em_parflow ${EMI_EM_PARFLOW_SOURCES})

set(EMI_LIBRARIES emi_em_parflow;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES emi_em_;mrflow${EMI_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS emi_em_parflow DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
