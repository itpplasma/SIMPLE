
if (EXISTS wout_unittests.nc)
    configure_file(
        wout_unittests.nc 
        wout.nc
        COPYONLY
    )
endif()

add_pfunit_ctest(test_example
    TEST_SOURCES test_example.f90
    LINK_LIBRARIES SIMPLE_LIB)

#add_pfunit_ctest(test_vmec
#    TEST_SOURCES test_vmec.f90
#    LINK_LIBRARIES SIMPLE_LIB)

add_pfunit_ctest(test_highlevel
    TEST_SOURCES test_highlevel.f90
    LINK_LIBRARIES SIMPLE_LIB)

add_pfunit_ctest(test_lowlevel
    TEST_SOURCES test_lowlevel.f90
    LINK_LIBRARIES SIMPLE_LIB)



