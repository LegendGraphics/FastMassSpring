cmake_minimum_required( VERSION 2.8 )
project( FastMassSpring )

LIST(APPEND CMAKE_CXX_FLAGS         "-std=c++11")

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

FIND_PACKAGE(3rdParty REQUIRED)

FILE( GLOB Project_SRCS "src/*.*" )
SET( Project_INCLUDE_DIR src/ )

INCLUDE_DIRECTORIES( ${Project_INCLUDE_DIR}
                     ${EIGEN_DIR} )
ADD_EXECUTABLE( ${PROJECT_NAME}
                ${Project_SRCS} )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} )