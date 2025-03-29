cmake_minimum_required(VERSION 3.10)
project(standard_pso C)

# Find GSL package
find_package(GSL REQUIRED)

# Add executable with all source files
add_executable(standard_pso 
    main.c
    rd.c
    coil_compression_spring.c
)

# Link against GSL
target_link_libraries(standard_pso ${GSL_LIBRARIES})
target_include_directories(standard_pso PRIVATE ${GSL_INCLUDE_DIRS}) 