set( src
    ${CMAKE_CURRENT_SOURCE_DIR}/GLheader.h
    ${CMAKE_CURRENT_SOURCE_DIR}/hair.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/hair.h
    ${CMAKE_CURRENT_SOURCE_DIR}/input.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/input.h
    ${CMAKE_CURRENT_SOURCE_DIR}/pic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pic.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ppm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/show.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/show.h
    ${CMAKE_CURRENT_SOURCE_DIR}/simu.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/simu.h
)


add_executable(
    HairSimu ${src}
)


find_package(OpenGL REQUIRED)

target_link_libraries(HairSimu
    OpenGL::GL
)

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    target_link_libraries(HairSimu "freeglutd")
elseif(CMAKE_BUILD_TYPE MATCHES "Release")
    target_link_libraries(HairSimu "freeglut")
endif()