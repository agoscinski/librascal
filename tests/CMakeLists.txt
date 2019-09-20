find_package(Boost REQUIRED COMPONENTS unit_test_framework )

add_library(test_runner STATIC test_runner.cc)
target_compile_options(test_runner PUBLIC -DBOOST_TEST_DYN_LINK)
target_link_libraries(test_runner ${Boost_LIBRARIES})

function(define_test _file_)
    get_filename_component(_name_ ${_file_} NAME_WE)

    add_executable(${_name_} ${_file_})
    target_link_libraries(${_name_} test_runner "${LIBRASCAL_NAME}")
    add_test(NAME ${_name_}
        COMMAND ${_name_} --report_level=detailed --build_info=TRUE
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_dependencies(${_name_} copy_tests_data)
endfunction()

# build c++ tests and add to test suite
file(GLOB TEST_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/test_*.cc")
list(REMOVE_ITEM TEST_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/test_runner.cc")

foreach(test_file IN LISTS TEST_SRCS)
    define_test(${test_file})
endforeach()

add_custom_target(copy_tests_data ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
    COMMENT "Copying tests data"
)

if(BUILD_BINDINGS)
    add_custom_target(copy_test_py ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/python ${CMAKE_CURRENT_BINARY_DIR}/
        COMMENT "Copying python tests"
    )

    add_test(
        NAME python_binding_test
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/python_binding_tests.py
    )

    set_tests_properties(python_binding_test
        PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}"
    )
endif()
