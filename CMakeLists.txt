cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

project(base)

# Find Geant4 package
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# Find ROOT package
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

# Gather source files
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

file(GLOB MACRO_FILES ${PROJECT_SOURCE_DIR}/macros/*.mac)
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR}/Release) # Necessary copy
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Release/Temportary)

file(GLOB DATA_FILES data/*.dat)
file(COPY ${DATA_FILES} DESTINATION ${PROJECT_BINARY_DIR})

file(GLOB DATA_FILES data/*.root)
file(COPY ${DATA_FILES} DESTINATION ${PROJECT_BINARY_DIR})

# Create executable
add_executable(basic_setup basic_setup.cc ${sources} ${headers})
target_include_directories(basic_setup PRIVATE include)

# Link Geant4 and ROOT libraries
target_link_libraries(basic_setup ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

# Create custom target
add_custom_target(base DEPENDS basic_setup)