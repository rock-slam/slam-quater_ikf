configure_file(quater_ikf.pc.in quater_ikf.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/quater_ikf.pc
    DESTINATION lib/pkgconfig)
install(FILES Ikf.hpp
    AdaptiveAttitudeCov.hpp
    DESTINATION include/quater_ikf)
