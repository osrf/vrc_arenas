set (files 
  coupling.stl
  nozzle.stl
)

install(FILES ${files} DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/models/${MODEL_NAME}/meshes)
