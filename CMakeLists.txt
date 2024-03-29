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

find_package(TBB REQUIRED)

# FetchContent_Declare(blaze
#   GIT_REPOSITORY https://github.com/dendisuhubdy/blaze.git
#   GIT_TAG master
# )
# FetchContent_MakeAvailable(blaze)
# find_package(blaze)

add_compile_options(-mtune=native -fconstexpr-backtrace-limit=0 -fconstexpr-steps=4290000000)
# add_compile_options(-mtune=native -fexperimental-new-constant-interpreter)
add_executable(lab7
    src/main.cpp
)

add_executable(lab7_2
   src/main2.cpp
   src/Expression.hpp
   src/Matrix.hpp
)

add_executable(tracing
   src/trace.hpp
   src/trace.cpp
)

add_executable(analysis
    src/analysis.hpp
   src/test_analysis.cpp

   src/graph.hpp
   src/graph.cpp
   src/problem.hpp src/problem.cpp
)

add_executable(simple_elimination_test
    test/elimination_test.cpp
    src/Matrix.hpp
)
target_link_libraries(simple_elimination_test fmt::fmt)


find_package(OpenMP)
add_executable(mm1
    mm/mm1.c
)
target_link_libraries(mm1 PUBLIC OpenMP::OpenMP_CXX)

target_link_libraries(lab7 fmt::fmt)
target_link_libraries(tracing fmt::fmt -fsanitize=address)
# target_link_libraries(analysis fmt::fmt -fsanitize=address)
# target_link_libraries(analysis fmt::fmt TBB::tbb OpenMP::OpenMP_CXX -fsanitize=thread)
target_link_libraries(analysis fmt::fmt TBB::tbb OpenMP::OpenMP_CXX)
# set_property(TARGET analysis PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

target_link_libraries(lab7_2 fmt::fmt)

include(GNUInstallDirs)
install(TARGETS lab7
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

enable_testing()
add_test(
  NAME 3x3_elimination_test
  COMMAND $<TARGET_FILE:simple_elimination_test> ${CMAKE_CURRENT_LIST_DIR}/in.txt ${CMAKE_CURRENT_LIST_DIR}/out.txt
)

add_test(
  NAME 15x15_elimination_test
  COMMAND $<TARGET_FILE:simple_elimination_test> ${CMAKE_CURRENT_LIST_DIR}/in15.txt ${CMAKE_CURRENT_LIST_DIR}/out15.txt
)
