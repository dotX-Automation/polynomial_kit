cmake_minimum_required(VERSION 3.8)
project(polynomial_kit)

set(CMAKE_BUILD_TYPE "RelWithDebInfo")

# Default to C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

# Default to C++17
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic -fvisibility=hidden)
endif()

# Find dependencies
find_package(ament_cmake REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

# Dynamic systems base library configuration
add_library(polynomial_kit SHARED
  src/polynomial_kit/polynomial.cpp)
target_compile_features(polynomial_kit PUBLIC c_std_99 cxx_std_17)  # Require C99 and C++17
target_compile_definitions(polynomial_kit PRIVATE "POLYNOMIAL_KIT_BUILDING_LIBRARY")
target_include_directories(polynomial_kit PUBLIC
  ${EIGEN3_INCLUDE_DIR}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(polynomial_kit
  Eigen3::Eigen)

# Libraries installation
install(
  DIRECTORY include/
  DESTINATION include)
install(
  TARGETS polynomial_kit
  EXPORT polynomial_kitTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

# Export all dependencies and library targets for this package
ament_export_targets(polynomial_kitTargets HAS_LIBRARY_TARGET)

ament_package()
