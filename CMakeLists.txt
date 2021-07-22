cmake_minimum_required(VERSION 3.1...3.15)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(
	Bubble
	VERSION 0.1
	DESCRIPTION "Library for random sampling from arbitrary distributions."
	LANGUAGES CXX)

add_library(bubble INTERFACE)
add_library(Bubble::bubble ALIAS bubble)
target_sources(bubble INTERFACE bubble.hpp)
target_include_directories(bubble INTERFACE ".")

find_package(ROOT 6.16 REQUIRED COMPONENTS Core Foam CONFIG)
find_package(OpenMP 3.0)

add_executable(test test.cpp)
if(${OpenMP_FOUND})
	target_compile_definitions(test PRIVATE BUBBLE_USE_OPENMP)
	if(${CMAKE_VERSION} VERSION_LESS "3.9.0")
		target_link_libraries(test PRIVATE "${OpenMP_CXX_FLAGS}")
		target_compile_options(test PRIVATE "${OpenMP_CXX_FLAGS}")
	else()
		target_link_libraries(test PRIVATE OpenMP::OpenMP_CXX)
	endif()
endif()
target_link_libraries(test PRIVATE bubble ROOT::Core ROOT::Foam)
target_compile_features(test PRIVATE cxx_std_11)
set_target_properties(test PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(
	test
	PRIVATE
	$<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -pedantic>
	$<$<CXX_COMPILER_ID:MSVC>:/W4>)

include(GNUInstallDirs)
install(
	FILES bubble.hpp
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/bubble")
install(
	TARGETS bubble
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

