cmake_minimum_required(VERSION 3.16)

# ── Version: single source of truth is _version.py at the repo root ──────────
# Re-configure automatically when _version.py changes
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/_version.py")
file(READ "${CMAKE_SOURCE_DIR}/_version.py" _version_py_contents)
string(REGEX MATCH "__version__ = '([^']+)'" _version_match "${_version_py_contents}")
set(MARLIM3_VERSION "${CMAKE_MATCH_1}")
if(NOT MARLIM3_VERSION)
    message(FATAL_ERROR "Could not parse __version__ from _version.py")
endif()

project(Marlim3 VERSION ${MARLIM3_VERSION} LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(SRC_ROOT ${CMAKE_SOURCE_DIR}/src)

set(SRCDIR  ${SRC_ROOT}/core)
set(FSRCDIR ${SRC_ROOT}/fortran)
set(INCDIR  ${SRC_ROOT}/include)
set(LIBDIR  ${SRC_ROOT}/thirdparty)
set(MODDIR ${CMAKE_BINARY_DIR}/mod)
file(MAKE_DIRECTORY ${MODDIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${MODDIR})

# ── Auto-generate versao.h from _version.py ──────────────────────────────────
string(TIMESTAMP MARLIM3_BUILD_DATE "%d/%m/%Y %H:%M:%S")
if(DEFINED ENV{USER})
    set(MARLIM3_BUILD_USER "$ENV{USER}")
elseif(DEFINED ENV{USERNAME})
    set(MARLIM3_BUILD_USER "$ENV{USERNAME}")
else()
    set(MARLIM3_BUILD_USER "unknown")
endif()
configure_file(
    ${INCDIR}/versao.h.in
    ${INCDIR}/versao.h
    @ONLY
)

file(GLOB CXX_SOURCES ${SRCDIR}/*.cpp)
set(F90_SOURCES
    ${FSRCDIR}/Constants.f90
    ${FSRCDIR}/Newton.f90
    ${FSRCDIR}/PvtSimWAXFileImport.f90
    ${FSRCDIR}/PvtSimCTMFileImport.f90
    ${FSRCDIR}/VLECalculations.f90
    ${FSRCDIR}/PhaseProperties.f90
    ${FSRCDIR}/Obsoletes.f90
    ${FSRCDIR}/FlashExtras.f90
    ${FSRCDIR}/BlackOilFromCompositional.f90
    ${FSRCDIR}/MarlimComposicional.f90
    ${FSRCDIR}/ViscosityCalculationModule.f90
)

add_executable(Marlim3 ${CXX_SOURCES} ${F90_SOURCES})

set_source_files_properties(${CXX_SOURCES} PROPERTIES INCLUDE_DIRECTORIES "${INCDIR};${LIBDIR}")
set_source_files_properties(${F90_SOURCES} PROPERTIES INCLUDE_DIRECTORIES "${MODDIR}")

target_link_libraries(Marlim3 PRIVATE ${MARLIM_FORTRAN_LIBS})

set_source_files_properties(${FSRCDIR}/Newton.f90            PROPERTIES OBJECT_DEPENDS ${FSRCDIR}/Constants.f90)
set_source_files_properties(${FSRCDIR}/VLECalculations.f90   PROPERTIES OBJECT_DEPENDS ${FSRCDIR}/Constants.f90)
set_source_files_properties(${FSRCDIR}/PhaseProperties.f90   PROPERTIES OBJECT_DEPENDS "${FSRCDIR}/VLECalculations.f90;${FSRCDIR}/Constants.f90")

set_target_properties(Marlim3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

install(TARGETS Marlim3 DESTINATION bin)
