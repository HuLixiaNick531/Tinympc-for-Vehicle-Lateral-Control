cmake_minimum_required(VERSION 3.0.0)
project(TinyMPC VERSION 0.2.0 LANGUAGES CXX)


if(MSVC)
    # 启用 Unicode 字符集
    add_definitions(-DUNICODE -D_UNICODE)
    # 设置源文件编码为 UTF-8
    add_compile_options(/utf-8)
    # 禁用不安全函数警告（可选）
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include/Eigen)

add_subdirectory(src/tinympc)
add_subdirectory(examples)
add_subdirectory(test)