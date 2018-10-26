add_executable (template_example template_example.cpp)
target_link_libraries (template_example clyde-rbm ${MANDATORY_LIBRARIES})

#add_executable (provaDMD provaDMD.cpp)
#target_link_libraries (provaDMD clyde-rbm ${MANDATORY_LIBRARIES})

add_executable (stackOverflow stackOverflow.cpp)
target_link_libraries (stackOverflow clyde-rbm ${MANDATORY_LIBRARIES})

add_executable (LaminarAirfoil LaminarAirfoil.cpp)
target_link_libraries (LaminarAirfoil clyde-rbm ${MANDATORY_LIBRARIES})

add_executable (prova prova.cpp)
target_link_libraries (prova clyde-rbm ${MANDATORY_LIBRARIES})