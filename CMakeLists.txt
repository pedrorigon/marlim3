cmake_minimum_required(VERSION 3.16)

# Version: single source of truth is _version.py; re-configure when it changes.
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/_version.py")
file(READ "${CMAKE_SOURCE_DIR}/_version.py" _version_py_contents)
string(REGEX MATCH "__version__ = '([^']+)'" _version_match "${_version_py_contents}")
set(MARLIM3_VERSION "${CMAKE_MATCH_1}")
if(NOT MARLIM3_VERSION)
    message(FATAL_ERROR "Could not parse __version__ from _version.py")
endif()

project(Marlim3 VERSION ${MARLIM3_VERSION} LANGUAGES CXX Fortran)

option(MARLIM_STATIC_GNU_RUNTIMES "Prefer static GNU runtimes where supported" ON)
option(MARLIM_MINGW_STATIC "Build a fully static executable with MinGW" ON)
set(MARLIM_FORTRAN_LIBS "" CACHE STRING "Additional Fortran runtime libraries to link explicitly")

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

# Auto-generate versao.h from _version.py.
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

# Auto-generate translations_data.h from translations.json.
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/marlim3/translations.json")
file(READ "${CMAKE_SOURCE_DIR}/marlim3/translations.json" TRANSLATIONS_JSON_CONTENT)
configure_file(
    ${INCDIR}/translations_data.h.in
    ${INCDIR}/translations_data.h
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

# Remove implicit dynamic Fortran runtimes before target creation to preserve static links.
if(MARLIM_STATIC_GNU_RUNTIMES AND UNIX AND NOT APPLE
        AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
        AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(MARLIM_LINUX_STATIC_RUNTIMES ON)
    list(REMOVE_ITEM CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES gfortran quadmath)
else()
    set(MARLIM_LINUX_STATIC_RUNTIMES OFF)
endif()

add_executable(Marlim3 ${CXX_SOURCES} ${F90_SOURCES})

target_include_directories(Marlim3 PRIVATE
    ${INCDIR}
    ${LIBDIR}
    ${MODDIR}
)

target_compile_options(Marlim3 PRIVATE
    $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:-Wall>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-fmessage-length=0>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fmax-errors=3>
    $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>
    $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>>:-ffp-contract=off>
    $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>>:-funroll-loops>
    $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>>:-O0>
    $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>>:-g3>
    $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-O0>
    $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-g3>
)

# MinGW uses compile-only -fopenmp so -static can select libgomp.a instead of libgomp.dll.a.
if(MARLIM_STATIC_GNU_RUNTIMES AND MINGW AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(MARLIM_MINGW_STATIC_RUNTIMES ON)
else()
    set(MARLIM_MINGW_STATIC_RUNTIMES OFF)
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)
if(MARLIM_LINUX_STATIC_RUNTIMES OR MARLIM_MINGW_STATIC_RUNTIMES)
    # Static-runtime paths use compile-only OpenMP; libgomp is linked explicitly below.
    target_compile_options(Marlim3 PRIVATE -fopenmp)
else()
    target_link_libraries(Marlim3 PRIVATE OpenMP::OpenMP_CXX)
endif()

# Static GNU runtimes keep Linux binaries independent from host compiler modules.
if(MARLIM_LINUX_STATIC_RUNTIMES)
    # Link GNU runtimes statically; keep pthread/dl dynamic for static libgomp dependencies.
    target_link_options(Marlim3 PRIVATE -static-libgcc -static-libstdc++)
    target_link_libraries(Marlim3 PRIVATE
        -Wl,-Bstatic gfortran quadmath gomp -Wl,-Bdynamic
        -Wl,--push-state,-Bdynamic pthread dl -Wl,--pop-state)

elseif(MARLIM_STATIC_GNU_RUNTIMES AND MINGW
        AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_options(Marlim3 PRIVATE -static-libgcc -static-libstdc++)
    if(MARLIM_MINGW_STATIC)
        target_link_options(Marlim3 PRIVATE -static)
    endif()
    # Link libgomp/winpthread explicitly so MinGW static builds do not import their DLLs.
    target_link_libraries(Marlim3 PRIVATE gomp pthread)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_link_libraries(Marlim3 PRIVATE gfortran quadmath)
    endif()

elseif(UNIX AND NOT APPLE AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(Marlim3 PRIVATE gfortran quadmath)
endif()

if(MARLIM_FORTRAN_LIBS)
    target_link_libraries(Marlim3 PRIVATE ${MARLIM_FORTRAN_LIBS})
endif()

set_source_files_properties(${FSRCDIR}/Newton.f90            PROPERTIES OBJECT_DEPENDS ${FSRCDIR}/Constants.f90)
set_source_files_properties(${FSRCDIR}/VLECalculations.f90   PROPERTIES OBJECT_DEPENDS ${FSRCDIR}/Constants.f90)
set_source_files_properties(${FSRCDIR}/PhaseProperties.f90   PROPERTIES OBJECT_DEPENDS "${FSRCDIR}/VLECalculations.f90;${FSRCDIR}/Constants.f90")

set_target_properties(Marlim3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Copy executable into the Python package directory after each build.
add_custom_command(TARGET Marlim3 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:Marlim3>
        ${CMAKE_SOURCE_DIR}/marlim3/$<TARGET_FILE_NAME:Marlim3>
    COMMENT "Copying Marlim3 executable to marlim3/ package directory"
)

install(TARGETS Marlim3 DESTINATION bin)
