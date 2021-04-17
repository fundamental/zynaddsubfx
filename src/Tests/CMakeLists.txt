function(cp_script script_name)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${script_name}
		${CMAKE_CURRENT_BINARY_DIR}/${script_name} COPYONLY)
endfunction()

#for tests looking for files stored in the source dir
add_definitions(-DSOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

CXXTEST_ADD_TEST(ControllerTest ControllerTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/ControllerTest.h)
CXXTEST_ADD_TEST(EchoTest EchoTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/EchoTest.h)
#CXXTEST_ADD_TEST(SampleTest SampleTest.h)
CXXTEST_ADD_TEST(MicrotonalTest MicrotonalTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/MicrotonalTest.h)
CXXTEST_ADD_TEST(XMLwrapperTest XMLwrapper.cpp ${CMAKE_CURRENT_SOURCE_DIR}/XMLwrapperTest.h)
CXXTEST_ADD_TEST(ADnoteTest AdNoteTest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AdNoteTest.h)
CXXTEST_ADD_TEST(SUBnoteTest SubNoteTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/SubNoteTest.h)
CXXTEST_ADD_TEST(OscilGenTest OscilGenTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/OscilGenTest.h)
CXXTEST_ADD_TEST(RandTest RandTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/RandTest.h)
CXXTEST_ADD_TEST(TensorTest TensorTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/TensorTest.h)
CXXTEST_ADD_TEST(PADnoteTest PadNoteTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/PadNoteTest.h)
CXXTEST_ADD_TEST(PluginTest PluginTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/PluginTest.h)
CXXTEST_ADD_TEST(MsgParseTest MsgParseTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/MsgParseTest.h)
CXXTEST_ADD_TEST(MiddlewareTest MiddlewareTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/MiddlewareTest.h)
CXXTEST_ADD_TEST(MessageTest MessageTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/MessageTest.h)
CXXTEST_ADD_TEST(UnisonTest UnisonTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/UnisonTest.h)
CXXTEST_ADD_TEST(MqTest MqTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/MqTest.h)
CXXTEST_ADD_TEST(WatchTest WatchTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/WatchTest.h)
#CXXTEST_ADD_TEST(RtAllocTest RtAllocTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/RtAllocTest.h)
CXXTEST_ADD_TEST(AllocatorTest AllocatorTest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AllocatorTest.h)
CXXTEST_ADD_TEST(EffectTest EffectTest.cpp ${CMAKE_CURRENT_SOURCE_DIR}/EffectTest.h)
CXXTEST_ADD_TEST(KitTest KitTest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/KitTest.h)
CXXTEST_ADD_TEST(MemoryStressTest MemoryStressTest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MemoryStressTest.h)
CXXTEST_ADD_TEST(TriggerTest TriggerTest.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TriggerTest.h)

#Extra libraries added to make test and full compilation use the same library
#links for quirky compilers
set(test_lib zynaddsubfx_core ${GUI_LIBRARIES} ${ZLIB_LIBRARY} ${FFTW_LIBRARIES}
    ${MXML_LIBRARIES} pthread "-Wl,--no-as-needed -lpthread")

message(STATUS "Linking tests with: ${test_lib}")
target_link_libraries(ADnoteTest     ${test_lib})
target_link_libraries(SUBnoteTest    ${test_lib})
target_link_libraries(ControllerTest ${test_lib})
target_link_libraries(EchoTest       ${test_lib})
target_link_libraries(MicrotonalTest ${test_lib})
target_link_libraries(OscilGenTest   ${test_lib})
target_link_libraries(XMLwrapperTest ${test_lib})
target_link_libraries(RandTest       ${test_lib})
target_link_libraries(TensorTest     ${test_lib})
target_link_libraries(PADnoteTest    ${test_lib})
target_link_libraries(MqTest         ${test_lib})
target_link_libraries(WatchTest      ${test_lib})
target_link_libraries(TriggerTest      ${test_lib})
target_link_libraries(PluginTest     zynaddsubfx_core zynaddsubfx_nio
    zynaddsubfx_gui_bridge
    ${GUI_LIBRARIES} ${NIO_LIBRARIES} ${AUDIO_LIBRARIES})
target_link_libraries(MsgParseTest ${test_lib})
target_link_libraries(MiddlewareTest zynaddsubfx_core zynaddsubfx_nio
    zynaddsubfx_gui_bridge
    ${GUI_LIBRARIES} ${NIO_LIBRARIES} ${AUDIO_LIBRARIES})
target_link_libraries(MessageTest zynaddsubfx_core zynaddsubfx_nio
    zynaddsubfx_gui_bridge
    ${GUI_LIBRARIES} ${NIO_LIBRARIES} ${AUDIO_LIBRARIES})
target_link_libraries(UnisonTest    ${test_lib})
#target_link_libraries(RtAllocTest    ${test_lib})
target_link_libraries(AllocatorTest    ${test_lib})
target_link_libraries(KitTest    ${test_lib})
target_link_libraries(MemoryStressTest ${test_lib})
target_link_libraries(EffectTest ${test_lib})

#Testbed app
add_executable(ins-test InstrumentStats.cpp)
target_link_libraries(ins-test ${test_lib} rt)

if(LIBLO_FOUND)
    cp_script(check-ports.rb)
    add_test(PortChecker check-ports.rb)
endif()

add_executable(save-osc SaveOSC.cpp)
target_link_libraries(save-osc
                      zynaddsubfx_core zynaddsubfx_nio
                      zynaddsubfx_gui_bridge
                      ${GUI_LIBRARIES} ${NIO_LIBRARIES} ${AUDIO_LIBRARIES})
#this will be replaced with a for loop when the code will get more stable:
add_test(SaveOsc save-osc ${CMAKE_CURRENT_SOURCE_DIR}/../../instruments/examples/Arpeggio\ 1.xmz)

#message(STATUS "Plugin Test ${GUI_LIBRARIES} ${NIO_LIBRARIES} ${AUDIO_LIBRARIES}")

