rock_find_pkgconfig(Eigen3 eigen3 REQUIRED)

rock_library(quater_ikf
    SOURCES Ikf.cpp
    HEADERS Ikf.hpp AdaptiveAttitudeCov.hpp
    DEPS_PKGCONFIG base-lib
    )

rock_executable(quater_ikf_bin main.cpp
    DEPS quater_ikf)
