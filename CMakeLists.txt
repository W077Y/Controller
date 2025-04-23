cmake_minimum_required (VERSION 3.19)


set(target_name "CONTROLLER")

add_library(${target_name} STATIC)
target_compile_features(${target_name} PUBLIC cxx_std_20)


message(STATUS "    -> CONTROLLER-Comp:")
message(STATUS "########################")



target_include_directories(${target_name}
 PUBLIC  "${CMAKE_CURRENT_LIST_DIR}/inc"
)

target_sources(${target_name}
 PRIVATE "${CMAKE_CURRENT_LIST_DIR}/src/controller.cpp"
)

target_link_libraries(${target_name}
 PUBLIC EXMATH
)



