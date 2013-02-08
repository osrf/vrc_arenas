set (files 
  manifest.xml
  model.sdf
)

add_subdirectory(meshes)
add_subdirectory(materials)

install(FILES ${files} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/vrc_arenas/models/vrc_firehose)
