set(MODEL_NAME vrc_heightmap_4)
set (files 
  model.config 
  model.sdf
)

add_subdirectory(materials)

install(FILES ${files} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/models/${MODEL_NAME}/)
