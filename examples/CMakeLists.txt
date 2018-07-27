add_executable (template_example template_example.cpp)
target_link_libraries (template_example clyde-rbm ${MANDATORY_LIBRARIES})

add_executable (provaDMD provaDMD.cpp)
target_link_libraries (provaDMD clyde-rbm ${MANDATORY_LIBRARIES})

add_executable (stackOverflow stackOverflow.cpp)
target_link_libraries (stackOverflow clyde-rbm ${MANDATORY_LIBRARIES})