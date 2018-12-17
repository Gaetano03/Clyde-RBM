#add_executable (template_example template_example.cpp)
#target_link_libraries (template_example clyde-rbm ${MANDATORY_LIBRARIES})

#add_executable (provaDMD provaDMD.cpp)
#target_link_libraries (provaDMD clyde-rbm ${MANDATORY_LIBRARIES})

#add_executable (stackOverflow stackOverflow.cpp)
#target_link_libraries (stackOverflow clyde-rbm ${MANDATORY_LIBRARIES})

#add_executable (LaminarAirfoil LaminarAirfoil.cpp)
#target_link_libraries (LaminarAirfoil clyde-rbm ${MANDATORY_LIBRARIES})

#add_executable (prova prova.cpp)
#target_link_libraries (prova clyde-rbm ${MANDATORY_LIBRARIES})

#add_executable (Compute_best_Method Compute_best_Method.cpp)
#target_link_libraries (Compute_best_Method clyde-rbm ${MANDATORY_LIBRARIES})

#add_executable (Compute_response_surface Compute_response_surface.cpp)
#target_link_libraries (Compute_response_surface clyde-rbm ${MANDATORY_LIBRARIES})

add_executable (Compute_error_out_of_T Compute_error_out_of_T.cpp)
target_link_libraries (Compute_error_out_of_T clyde-rbm ${MANDATORY_LIBRARIES})

#add_executable (Compute_Adaptive_reconstruction Compute_Adaptive_reconstruction.cpp)
#target_link_libraries (Compute_Adaptive_reconstruction clyde-rbm ${MANDATORY_LIBRARIES})

add_executable (Adaptive_Energy_Mode Adaptive_Energy_Mode.cpp)
target_link_libraries (Adaptive_Energy_Mode clyde-rbm ${MANDATORY_LIBRARIES})

add_executable (Adaptive_reconstruction Adaptive_reconstruction.cpp)
target_link_libraries (Adaptive_reconstruction clyde-rbm ${MANDATORY_LIBRARIES})