cmake_minimum_required(VERSION 3.1...3.15)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(
	Bubble
	VERSION 0.1
	DESCRIPTION "Library for random sampling from arbitrary distributions."
	LANGUAGES CXX)

set(Bubble_BUILD_EXAMPLE NO
	CACHE BOOL
	"should the example (using ROOT) be built")

add_library(bubble INTERFACE)
add_library(Bubble::bubble ALIAS bubble)
target_sources(bubble INTERFACE bubble.hpp)
target_include_directories(bubble INTERFACE ".")

if(${Bubble_BUILD_EXAMPLE})
	find_package(ROOT 6.16 REQUIRED COMPONENTS Core Foam CONFIG)
	find_package(OpenMP 3.0)

	add_executable(example example.cpp)
	if(${OpenMP_FOUND})
		target_compile_definitions(example PRIVATE BUBBLE_USE_OPENMP)
		if(${CMAKE_VERSION} VERSION_LESS "3.9.0")
			target_link_libraries(example PRIVATE "${OpenMP_CXX_FLAGS}")
			target_compile_options(example PRIVATE "${OpenMP_CXX_FLAGS}")
		else()
			target_link_libraries(example PRIVATE OpenMP::OpenMP_CXX)
		endif()
	endif()
	target_link_libraries(example PRIVATE bubble ROOT::Core ROOT::Foam)
	target_compile_features(example PRIVATE cxx_std_11)
	set_target_properties(example PROPERTIES CXX_EXTENSIONS OFF)
	target_compile_options(
		example
		PRIVATE
		$<$<OR:$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall -Wextra -pedantic>
		$<$<CXX_COMPILER_ID:MSVC>:/W4>)
endif()

include(GNUInstallDirs)
install(
	FILES bubble.hpp
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/bubble")
install(
	TARGETS bubble
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

