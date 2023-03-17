cmake_minimum_required(VERSION 3.8)
project(mr_dual_arm_control)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)

find_package(robot_state_publisher REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(jsoncpp REQUIRED)
find_package(nav_msgs REQUIRED)


add_library(MR SHARED include/ModernRobotics/ModernRobotics.cpp)
add_library(MR2 SHARED include/ModernRobotics_origin/modern_robotics.cpp)
add_executable(mr_dual_arm_control src/mr_dual_arm_control.cpp include/Indy7/Indy7.cpp )
target_include_directories(mr_dual_arm_control PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/ModernRobotics>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/ModernRobotics_origin>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/Indy7>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/json>
  $<INSTALL_INTERFACE:include>)
  target_link_libraries(mr_dual_arm_control jsoncpp_lib MR MR2 )
  target_compile_features(mr_dual_arm_control PUBLIC c_std_99 cxx_std_17)  # Require C99 and C++17
ament_target_dependencies(
  mr_dual_arm_control
  "robot_state_publisher"
  "sensor_msgs"  
  "std_msgs"
  "geometry_msgs"
  "nav_msgs"
  "jsoncpp"
)
add_custom_command(TARGET mr_dual_arm_control POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/MR_info.json lib/${PROJECT_NAME})
add_custom_target(copy_MR_info ALL DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/MR_info.json)


install(TARGETS mr_dual_arm_control
  DESTINATION lib/${PROJECT_NAME})

install(
    DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/urdf/
    DESTINATION share/${PROJECT_NAME}
)
install(
    DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/meshes/visual/
    DESTINATION share/${PROJECT_NAME}/meshes/visual/
)
install(
    DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/meshes/collision/
    DESTINATION share/${PROJECT_NAME}/meshes/collision/
)
install(
    DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/launch/
    DESTINATION share/${PROJECT_NAME}
)

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

ament_package()
