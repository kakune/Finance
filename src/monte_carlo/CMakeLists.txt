# ソースファイルを明示的にリストアップ
set(SOURCES generator.cpp process.cpp)

# ライブラリを作成
add_library(monte_carlo ${SOURCES})
