cmake_minimum_required(VERSION 3.5)

project(lab7 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG master
)
FetchContent_MakeAvailable(fmt)
find_package(fmt)

FetchContent_Declare(blaze
  GIT_REPOSITORY https://github.com/dendisuhubdy/blaze.git
  GIT_TAG master
)
FetchContent_MakeAvailable(blaze)
find_package(blaze)

add_executable(lab7
    main.cpp
)

add_executable(lab7_2
    main2.cpp
    Expression.hpp
    Matrix.hpp
)
find_package(OpenMP)
add_executable(mm1
    mm1.c
)
target_link_libraries(mm1 PUBLIC OpenMP::OpenMP_CXX)

target_link_libraries(lab7 fmt::fmt)
target_link_libraries(lab7_2 fmt::fmt blaze::blaze)

include(GNUInstallDirs)
install(TARGETS lab7
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
