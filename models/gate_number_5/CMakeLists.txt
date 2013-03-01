set (files 
  model.config
  model.sdf
)

add_subdirectory(meshes)
add_subdirectory(materials)

install(FILES ${files} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/vrc_arenas/models/gate_number_5)
