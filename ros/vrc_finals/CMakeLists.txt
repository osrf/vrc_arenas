set(ROS_PACKAGE_NAME vrc_finals)

add_subdirectory(launch)
add_subdirectory(test)

set (files
  manifest.xml
)

install(FILES ${files} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/ros/vrc_finals)
