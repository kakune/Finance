# Find GSL
find_package(GSL REQUIRED)

# ソースファイルを明示的にリストアップ
set(SOURCES fft.cpp integral.cpp find_root.cpp special_functions.cpp)

# ライブラリを作成
add_library(numerical ${SOURCES})

# Include GSL headers
target_include_directories(numerical PRIVATE ${GSL_INCLUDE_DIRS})

# Link against GSL libraries
target_link_libraries(numerical ${GSL_LIBRARIES})