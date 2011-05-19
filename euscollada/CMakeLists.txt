cmake_minimum_required(VERSION 2.4.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
#set(ROS_BUILD_TYPE RelWithDebInfo)

rosbuild_init()

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

#uncomment if you have defined messages
#rosbuild_genmsg()
#uncomment if you have defined services
#rosbuild_gensrv()

#common commands for building c++ executables and libraries
#rosbuild_add_library(${PROJECT_NAME} src/example.cpp)
#target_link_libraries(${PROJECT_NAME} another_library)
#rosbuild_add_boost_directories()
#rosbuild_link_boost(${PROJECT_NAME} thread)
#rosbuild_add_executable(example examples/example.cpp)
#target_link_libraries(example ${PROJECT_NAME})

rosbuild_add_executable(collada2eus src/collada2eus.cpp)
target_link_libraries(collada2eus qhull)
message("${CMAKE_SYSTEM_PROCESSOR}")
macro(convert_to_collada)
  set(convert_exe "convert-irtmodel-to-collada.l")
  rosbuild_invoke_rospack("euscollada" _rospack euscollada_path find)
  set(convert_exe_full_path "${_rospack_euscollada_path}/scripts/${convert_exe}")
  #message("${convert_exe_full_path}") #debug
  # search all irtmodels
  rosbuild_invoke_rospack("euslisp" _rospack euslisp_path find)
  rosbuild_invoke_rospack("${PROJECT_NAME}" _rospack current_package_path find)
  set(irtmodel_directory "${_rospack_euslisp_path}/jskeus/eus/models")
  set(collada_output_directory "${_rospack_current_package_path}/collada")
  execute_process(COMMAND "find" "${irtmodel_directory}" "-name" "*robot.l" "-o" "-name" "*object.l"
    OUTPUT_VARIABLE irtmodels)

  string(REPLACE "\n" ";" irtmodels "${irtmodels}") #space separated list
  set(irtmodels_list ${irtmodels})
  set(_autogen "")
  foreach(_model ${irtmodels_list})
    string(REGEX REPLACE ".l" "" _model_name1 ${_model})
    string(REGEX REPLACE ".*/" "" _model_name ${_model_name1})
    set(_command_arg "\"\\(progn \\(defvar *model-file* \\\"${_model}\\\"\\) \\(defvar *model-name* \\\"${_model_name}\\\"\\) \\(defvar  *output-uri* \\\"package://${PROJECT_NAME}\\\"\\)\\)\" ")
    add_custom_command(OUTPUT "${collada_output_directory}/${_model_name}.dae"
      COMMAND "roseus" "\"${_command_arg}\"" "${convert_exe_full_path}"
      DEPENDS "${convert_exe_full_path}" "${_model}"
      "${_rospack_euscollada_path}/src/collada-interface.l"
      "${_rospack_euscollada_path}/src/eusmodel-description.l")
    list(APPEND _autogen "${collada_output_directory}/${_model_name}.dae")
  endforeach(_model)
  add_custom_target(CONVERT_EUS_TO_COLLADA ALL DEPENDS ${_autogen})
endmacro(convert_to_collada)

convert_to_collada()


## test
rosbuild_download_test_data(https://openrave.svn.sourceforge.net/svnroot/openrave/data/robots/mitsubishi-pa10.zae test/mitsubishi-pa10.zae 2e995a1356befd581028b9fd0d55b10f)
rosbuild_download_test_data(https://openrave.svn.sourceforge.net/svnroot/openrave/data/robots/unimation-pumaarm.zae test/unimation-pumaarm.zae fc43ed705335e77c3e1d9634b2bf39ce)
rosbuild_download_test_data(https://openrave.svn.sourceforge.net/svnroot/openrave/data/robots/care-o-bot3.zae test/care-o-bot3.zae a745f8f4f4773cf3cedb9670619e472b)
rosbuild_download_test_data(https://openrave.svn.sourceforge.net/svnroot/openrave/data/robots/darpa-arm.zae test/darpa-arm.zae a00ba4ab0a79c484894e5e9f7cccb790)

rosbuild_add_pyunit(test/collada2eus_test.py)

