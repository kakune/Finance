cmake_minimum_required(VERSION 3.10)

# プロジェクト名と使用するC++のバージョン
project(MyProject VERSION 1.0.0 LANGUAGES CXX)

# C++の標準を設定
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# IntelのC++コンパイラicpxを使用
set(CMAKE_CXX_COMPILER "icpx")

# O3最適化を有効にする
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# 出力ディレクトリを設定
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/out)

# サブディレクトリを追加
add_subdirectory(src/core)
add_subdirectory(src/numerical)
add_subdirectory(src/analytical)
add_subdirectory(src/option)
add_subdirectory(src/monte_carlo)
add_subdirectory(src/fourier)

# インクルードパスを設定
include_directories(${PROJECT_SOURCE_DIR}/src)

# 実行可能ファイルを作成
add_executable(out "main_calc.cpp")

# ライブラリをリンク
target_link_libraries(out core numerical analytical option monte_carlo fourier)