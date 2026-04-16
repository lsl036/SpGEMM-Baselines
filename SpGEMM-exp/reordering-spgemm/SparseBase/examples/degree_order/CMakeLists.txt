add_executable(degree_order ${CMAKE_CURRENT_SOURCE_DIR}/degree_order.cc)
target_link_libraries(degree_order sparsebase)

add_executable(degree_order_mtx ${CMAKE_CURRENT_SOURCE_DIR}/degree_order_mtx.cc)
target_link_libraries(degree_order_mtx sparsebase)