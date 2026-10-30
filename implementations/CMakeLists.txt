cmake_minimum_required(VERSION 3.16)

# Builds standalone (cmake -S implementations -B build) or as a subdirectory of a larger project
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	project(TimelineFXImplementations C CXX)
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

set(TFX_LIBRARY_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(TFX_LIBRARY_SOURCE ${TFX_LIBRARY_FOLDER}/timelinefx.cpp ${TFX_LIBRARY_FOLDER}/timelinefx.h)
set(TFX_IMPL_COMMON_SOURCE
	${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_common.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_common.h
	# Decompress only zstd for supercompressed KTX2 shapes
	${CMAKE_CURRENT_SOURCE_DIR}/third_party/zstddeclib.c
)

function(tfx_example_target target_name)
	target_include_directories(${target_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${TFX_LIBRARY_FOLDER})
	target_compile_definitions(${target_name} PRIVATE TFX_OUTPUT_NOTICE_MESSAGES TFX_IMPL_KTX2_ZSTD TFX_EXAMPLE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}")
	if (MSVC)
		target_compile_definitions(${target_name} PRIVATE _CRT_SECURE_NO_WARNINGS)
	endif()
	set_target_properties(${target_name} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	source_group("TimelineFX" FILES ${TFX_LIBRARY_SOURCE})
endfunction()

# Every shader the examples use is built from the library's one Slang source
find_program(TFX_SLANGC slangc HINTS "$ENV{VULKAN_SDK}/Bin")
set(TFX_SHADER_SOURCE ${TFX_LIBRARY_FOLDER}/Shaders/timelinefx.slang)
if (NOT TFX_SLANGC)
	message(FATAL_ERROR "slangc wasn't found. It comes with the Vulkan SDK, or with a Slang release from github.com/shader-slang/slang, put its bin folder on the PATH or set TFX_SLANGC to slangc")
endif()

if (WIN32)
	# Each entry point with its D3D11 profile, compiled to <entry>.cso. fxc can't compile Slang's #line directives, hence -line-directive-mode none
	set(TFX_DX11_SHADER_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/dx11_shaders)
	set(TFX_DX11_ENTRY_POINTS
		simulate_particles cs_5_0
		tfx_billboard_vertex vs_5_0
		tfx_billboard_baked_vertex vs_5_0
		tfx_billboard_fragment ps_5_0
		tfx_ribbon_compute cs_5_0
		tfx_ribbon_vertex vs_5_0
		tfx_ribbon_fragment ps_5_0
		tfx_sprite_playback cs_5_0
		tfx_ribbon_playback cs_5_0
	)
	set(TFX_DX11_SHADERS)
	list(LENGTH TFX_DX11_ENTRY_POINTS entry_point_list_length)
	math(EXPR last_entry_point_index "${entry_point_list_length} - 1")
	foreach(entry_point_index RANGE 0 ${last_entry_point_index} 2)
		math(EXPR profile_index "${entry_point_index} + 1")
		list(GET TFX_DX11_ENTRY_POINTS ${entry_point_index} entry_point)
		list(GET TFX_DX11_ENTRY_POINTS ${profile_index} profile)
		set(shader_output ${TFX_DX11_SHADER_FOLDER}/${entry_point}.cso)
		add_custom_command(OUTPUT ${shader_output}
			COMMAND ${CMAKE_COMMAND} -E make_directory ${TFX_DX11_SHADER_FOLDER}
			COMMAND ${TFX_SLANGC} ${TFX_SHADER_SOURCE} -entry ${entry_point} -target dxbc -profile ${profile} -D TFX_FIXED_BINDINGS -line-directive-mode none -o ${shader_output}
			DEPENDS ${TFX_SHADER_SOURCE}
			COMMENT "Compiling ${entry_point} for D3D11")
		list(APPEND TFX_DX11_SHADERS ${shader_output})
	endforeach()
	add_custom_target(tfx_dx11_shaders DEPENDS ${TFX_DX11_SHADERS})

	add_executable(example_win32_dx11
		${CMAKE_CURRENT_SOURCE_DIR}/example_win32_dx11/main.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.h
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_dx11.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_dx11.h
		${TFX_SHADER_SOURCE}
		${TFX_IMPL_COMMON_SOURCE}
		${TFX_LIBRARY_SOURCE}
	)
	# slangc builds it, so Visual Studio mustn't try to
	set_source_files_properties(${TFX_SHADER_SOURCE} PROPERTIES VS_TOOL_OVERRIDE "None")
	source_group("Shaders" FILES ${TFX_SHADER_SOURCE})
	tfx_example_target(example_win32_dx11)
	target_link_libraries(example_win32_dx11 PRIVATE d3d11 dxgi)
	add_dependencies(example_win32_dx11 tfx_dx11_shaders)
	target_compile_definitions(example_win32_dx11 PRIVATE TFX_EXAMPLE_SHADER_FOLDER="${TFX_DX11_SHADER_FOLDER}")
endif()
