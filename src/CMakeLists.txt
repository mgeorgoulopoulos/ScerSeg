cmake_minimum_required(VERSION 3.1.0)

project(yeast)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

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

add_executable(01_MapGenes3D "MapGenes3D.cpp")
target_link_libraries(01_MapGenes3D Qt5::Core	Qt5::Gui Qt5::Widgets Qt5::Sql)

add_executable(02_CsvPatcher "CsvPatcher.cpp")
target_link_libraries(02_CsvPatcher Qt5::Core	Qt5::Gui Qt5::Widgets Qt5::Sql)

add_executable(03_CommunityScore "CommunityScore.cpp")
target_link_libraries(03_CommunityScore Qt5::Core	Qt5::Gui Qt5::Widgets Qt5::Sql)

add_executable(04_knn "knn.cpp")
target_link_libraries(04_knn Qt5::Core	Qt5::Gui Qt5::Widgets Qt5::Sql)

add_executable(05_CommunitySvg "CommunitySvg.cpp")
target_link_libraries(05_CommunitySvg Qt5::Core	Qt5::Gui Qt5::Widgets Qt5::Sql)

add_executable(06_CommunityEntropy "CommunityEntropy.cpp")
target_link_libraries(06_CommunityEntropy Qt5::Core	Qt5::Gui Qt5::Widgets Qt5::Sql)

add_executable(07_CommunityEntropyFilter "CommunityEntropyFilter.cpp")
target_link_libraries(07_CommunityEntropyFilter Qt5::Core	Qt5::Gui Qt5::Widgets Qt5::Sql)

add_executable(08_TightCommunityDistances "TightCommunityDistances.cpp")
target_link_libraries(08_TightCommunityDistances Qt5::Core	Qt5::Gui Qt5::Widgets Qt5::Sql)

add_executable(09_PromoterSpheres "PromoterSpheres.cpp")
target_link_libraries(09_PromoterSpheres Qt5::Core	Qt5::Gui Qt5::Widgets Qt5::Sql)




