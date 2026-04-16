add_executable(rcm_order ${CMAKE_CURRENT_SOURCE_DIR}/rcm_order.cc)
target_link_libraries(rcm_order sparsebase)

add_executable(rcm_order_mtx ${CMAKE_CURRENT_SOURCE_DIR}/rcm_order_mtx.cc)
target_link_libraries(rcm_order_mtx sparsebase)