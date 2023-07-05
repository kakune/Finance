# ソースファイルを明示的にリストアップ
set(SOURCES black_scholes.cpp)

# ライブラリを作成
add_library(analytical ${SOURCES})