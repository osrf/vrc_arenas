set (files 
task1a.world
task1b.world
task1c.world
task1d.world
task1e.world
task2a.world
task2b.world
task2c.world
task2d.world
task2e.world
task3a.world
task3b.world
task3c.world
task3d.world
task3e.world
)

install(FILES ${files} DESTINATION
  ${CMAKE_INSTALL_PREFIX}/share/vrc_arenas/worlds/)
