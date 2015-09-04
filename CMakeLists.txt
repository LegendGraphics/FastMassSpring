cmake_minimum_required( VERSION 2.8 )
project( FastMassSpring )

FIND_PATH( BOOST_DIR "boost" )
FIND_PATH( EIGEN_DIR "eigen" )

FILE( GLOB Project_SRCS "src/*.*" )
SET( Project_INCLUDE_DIR src/ )

INCLUDE_DIRECTORIES( ${Project_INCLUDE_DIR} )
ADD_EXECUTABLE( ${PROJECT_NAME}
                ${Project_SRCS} )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} )