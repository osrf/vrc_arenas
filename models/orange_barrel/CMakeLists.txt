set (files 
  manifest.xml
  model-1_3.sdf
)

add_subdirectory(meshes)

install(FILES ${files} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/vrc_arenas/models/orange_barrel)
