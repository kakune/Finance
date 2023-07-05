# ソースファイルを明示的にリストアップ
set(SOURCES vanilla.cpp)

# ライブラリを作成
add_library(fourier ${SOURCES})

# ライブラリの依存関係を設定
target_link_libraries(fourier numerical analytical)

# インクルードパスを設定
target_include_directories(fourier PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)