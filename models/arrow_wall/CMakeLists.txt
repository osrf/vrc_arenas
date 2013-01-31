set (files 
  manifest.xml
  model-1_3.sdf
)

add_subdirectory(meshes)
add_subdirectory(materials)

install(FILES ${files} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/vrc_arenas/models/arrow_wall)
