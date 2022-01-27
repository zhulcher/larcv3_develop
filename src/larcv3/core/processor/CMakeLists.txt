set(name processor)

# Get all the source files:
file(GLOB SOURCES *.cxx)
file(GLOB HEADERS *.h)

# Add a shared library
add_library(${name} OBJECT ${SOURCES})


find_package(EDepSim REQUIRED) #short term fix
target_link_libraries(${name} EDepSim::edepsim_io)# "/app/edep-sim/io") #short term fix




install (FILES ${HEADERS}
    DESTINATION ${CMAKE_PACKAGE_DIR}/include/larcv3/app/${name})
