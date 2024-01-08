#pragma once
#include "fmt/core.h"
#include "fmt/format.h"
#include <algorithm>
#include <array>
#include <iostream>
#include <type_traits>
#include <variant>
#include <vector>

/*
 * Goals:
 * - calculate max amount of threads that can be used for the computation
 * - calculate 2D pseudo-matrix of calculations to perform
 * 
 * 2 versions:
 * - dynamic
 * - compile-time-known
 * 
 * Possible optimizations:
 * - When running on a platform supporting SIMD (like modern x86), in a foata layer, group tasks to use SIMD functionality
 * 
 * Look at https://godbolt.org/z/zdEE3fr7o
*/
