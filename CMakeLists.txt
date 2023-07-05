cmake_minimum_required(VERSION 3.10)

# プロジェクト名と使用するC++のバージョン
project(MyProject VERSION 1.0.0 LANGUAGES CXX)

# C++の標準を設定
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# IntelのC++コンパイラicpxを使用
set(CMAKE_CXX_COMPILER "icpx")

# ソースファイルの自動探索
file(GLOB_RECURSE LIB_SOURCES "src/*.cpp")

# インクルードパスを設定
include_directories(${PROJECT_SOURCE_DIR}/src)

# 出力ディレクトリを設定
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/out)

# 静的ライブラリを作成
add_library(my_library STATIC ${LIB_SOURCES})

# 実行可能ファイルを作成
add_executable(out "main_calc.cpp")

# 静的ライブラリをリンク
target_link_libraries(out my_library)