#################################################
cmake_minimum_required(VERSION 3.10)
project(GalaxyKinematic)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#################################################
# include project sorce files
set(SRC_DIR ${CMAKE_SOURCE_DIR}/cppCode)
set(SRC_FILES
    ${SRC_DIR}/calccatalog.cpp
    ${SRC_DIR}/calckinematic.cpp
    ${SRC_DIR}/centroid.cpp
    ${SRC_DIR}/configloader.cpp
    ${SRC_DIR}/main.cpp
    ${SRC_DIR}/program.cpp
    ${SRC_DIR}/star.cpp
)
set(HEADER_FILES
    ${SRC_DIR}/calccatalog.h
    ${SRC_DIR}/calckinematic.h
    ${SRC_DIR}/centroid.h
    ${SRC_DIR}/configloader.h
    ${SRC_DIR}/program.h
    ${SRC_DIR}/star.h
)

#################################################
# include helper classes
set(HELPER_DIR ${CMAKE_SOURCE_DIR}/../CppHelperClasses)
set(HELPER_FILES
    ${HELPER_DIR}/TException.cpp
    ${HELPER_DIR}/TFiletable.cpp
    ${HELPER_DIR}/TLogger.cpp
    ${HELPER_DIR}/TMatrix.cpp
    ${HELPER_DIR}/TVersion.cpp
)
set(HELPER_HEADERS
    ${HELPER_DIR}/TException.h
    ${HELPER_DIR}/TFiletable.h
    ${HELPER_DIR}/TLogger.h
    ${HELPER_DIR}/TMatrix.h
    ${HELPER_DIR}/TVersion.h
)

#################################################
# include external libs
include(FetchContent)

# ALGLIB
FetchContent_Declare(alglib SOURCE_DIR ${CMAKE_SOURCE_DIR}/../CppExternalSourceLibs/alglib/3.19.0)
FetchContent_MakeAvailable(alglib)

# PUGIXML
FetchContent_Declare(pugixml SOURCE_DIR ${CMAKE_SOURCE_DIR}/../CppExternalSourceLibs/pugixml/1.14)
FetchContent_MakeAvailable(pugixml)

#################################################
add_executable(GalaxyKinematic
    ${SRC_FILES}
    ${HEADER_FILES}
    ${HELPER_FILES}
    ${HELPER_HEADERS}
)
target_include_directories(GalaxyKinematic PRIVATE
    ${SRC_DIR}
    ${HELPER_DIR}
)
target_link_libraries(GalaxyKinematic PRIVATE
    alglib
    pugixml
)
add_custom_command(TARGET GalaxyKinematic POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:alglib>
        $<TARGET_FILE:pugixml>
        $<TARGET_FILE_DIR:GalaxyKinematic>
)

#################################################
