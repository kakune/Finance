# ソースファイルを明示的にリストアップ
set(SOURCES fft.cpp find_root.cpp special_functions.cpp)

# ライブラリを作成
add_library(numerical ${SOURCES})
