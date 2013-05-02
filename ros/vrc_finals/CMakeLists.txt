add_subdirectory(launch)

set (files
  manifest.xml
)

install(FILES ${files} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/ros/vrc_finals)
