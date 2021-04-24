cmake_minimum_required(VERSION 3.1.0)

project(yeast)

set(INCLUDE_PATHS ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${INCLUDE_PATHS})

# need this to include the auto generated ui_mainwindow.h
# file from the qt5_wrap_ui call below.
set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(CMAKE_VERSION VERSION_LESS "3.7.0")
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif()

find_package(Qt5 COMPONENTS REQUIRED Core Gui Widgets Sql )

# wrap the ui file to a c++ header
#qt5_wrap_ui(WrappedForms ${FORMS})
#qt5_wrap_cpp(WrappedSources ${HEADERS})

#source_group("MOCS" FILES
#	${WrappedForms}
#	${WrappedSources}
#)

add_executable(MapGenes3D "MapGenes3D.cpp")
target_link_libraries(MapGenes3D Qt5::Core	Qt5::Gui Qt5::Widgets Qt5::Sql)

