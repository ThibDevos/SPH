include_directories("${PROJECT_SOURCE_DIR}/include")

add_executable(main_AoS main_AoS.cpp ../src/initialisation.cpp)
add_executable(main_SoA main_SoA.cpp ../src/initialisation.cpp)