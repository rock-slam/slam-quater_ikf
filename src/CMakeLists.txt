rock_library(quater_ikf
    SOURCES ikf.cpp
    HEADERS ikf.h)

rock_executable(quater_ikf_bin main.cpp
    DEPS quater_ikf)

