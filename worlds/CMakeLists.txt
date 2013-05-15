set (files 
vrc_final_task1.world
vrc_final_task2.world
vrc_final_task3.world
vrc_final_task4.world
vrc_final_task5.world
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
  ${CMAKE_INSTALL_PREFIX}/vrc_arenas/worlds/)
