# ソースファイルを明示的にリストアップ
set(SOURCES option_core.cpp european.cpp implied_vol.cpp)

# ライブラリを作成
add_library(option ${SOURCES})

# ライブラリの依存関係を設定
target_link_libraries(option numerical analytical)

# インクルードパスを設定
target_include_directories(option PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
