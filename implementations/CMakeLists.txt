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
	target_compile_definitions(${target_name} PRIVATE TFX_IMPL_KTX2_ZSTD TFX_EXAMPLE_ROOT="${CMAKE_CURRENT_SOURCE_DIR}")
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
	# Each entry point with its D3D11 profile, compiled to <entry>.cso with TFX_ENTRY_<entry> picking its fixed bindings. fxc can't compile Slang's #line directives, hence -line-directive-mode none
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
			COMMAND ${TFX_SLANGC} ${TFX_SHADER_SOURCE} -entry ${entry_point} -target dxbc -profile ${profile} -D TFX_FIXED_BINDINGS -D TFX_ENTRY_${entry_point} -line-directive-mode none -o ${shader_output}
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

# Off unless asked for, because it downloads and builds SDL3
option(TFX_EXAMPLE_SDL3 "Build the SDL3 GPU example, fetching SDL3 with FetchContent" OFF)
if (TFX_EXAMPLE_SDL3)
	include(FetchContent)
	set(SDL_SHARED OFF CACHE BOOL "" FORCE)
	set(SDL_STATIC ON CACHE BOOL "" FORCE)
	set(SDL_TEST_LIBRARY OFF CACHE BOOL "" FORCE)
	FetchContent_Declare(SDL3
		GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
		GIT_TAG release-3.4.18
		GIT_SHALLOW TRUE)
	FetchContent_MakeAvailable(SDL3)

	# SDL_gpu's sets and spaces, SPIR-V for Vulkan and DXIL for D3D12. Each entry point keeps its name in SPIR-V so SDL can find it
	set(TFX_SDL_SHADER_DEFINES -D TFX_FIXED_BINDINGS
		-D TFX_VERTEX_UNIFORM_SPACE=1 -D TFX_VERTEX_UNIFORM_BINDING=0
		-D TFX_FRAGMENT_TEXTURE_SPACE=2 -D TFX_FRAGMENT_TEXTURE_BINDING=0 -D TFX_FRAGMENT_SAMPLER_SPACE=2 -D TFX_FRAGMENT_SAMPLER_BINDING=0
		-D TFX_COMPUTE_READ_BINDING=0 -D TFX_COMPUTE_WRITE_SPACE=1 -D TFX_COMPUTE_UNIFORM_SPACE=2 -D TFX_COMPUTE_UNIFORM_BINDING=0
		-D TFX_COMBINED_SAMPLERS)
	set(TFX_SDL_SHADER_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/sdl_shaders)
	set(TFX_SDL_ENTRY_POINTS
		simulate_particles cs_6_0
		tfx_billboard_vertex vs_6_0
		tfx_billboard_baked_vertex vs_6_0
		tfx_billboard_fragment ps_6_0
		tfx_ribbon_compute cs_6_0
		tfx_ribbon_vertex vs_6_0
		tfx_ribbon_fragment ps_6_0
		tfx_sprite_playback cs_6_0
		tfx_ribbon_playback cs_6_0
	)
	set(TFX_SDL_SHADERS)
	list(LENGTH TFX_SDL_ENTRY_POINTS entry_point_list_length)
	math(EXPR last_entry_point_index "${entry_point_list_length} - 1")
	foreach(entry_point_index RANGE 0 ${last_entry_point_index} 2)
		math(EXPR profile_index "${entry_point_index} + 1")
		list(GET TFX_SDL_ENTRY_POINTS ${entry_point_index} entry_point)
		list(GET TFX_SDL_ENTRY_POINTS ${profile_index} profile)
		set(spirv_output ${TFX_SDL_SHADER_FOLDER}/${entry_point}.spv)
		add_custom_command(OUTPUT ${spirv_output}
			COMMAND ${CMAKE_COMMAND} -E make_directory ${TFX_SDL_SHADER_FOLDER}
			COMMAND ${TFX_SLANGC} ${TFX_SHADER_SOURCE} -entry ${entry_point} -target spirv -profile spirv_1_5 -O3 -fvk-use-entrypoint-name ${TFX_SDL_SHADER_DEFINES} -D TFX_ENTRY_${entry_point} -o ${spirv_output}
			DEPENDS ${TFX_SHADER_SOURCE}
			COMMENT "Compiling ${entry_point} for SDL_gpu Vulkan")
		list(APPEND TFX_SDL_SHADERS ${spirv_output})
		if (WIN32)
			set(dxil_output ${TFX_SDL_SHADER_FOLDER}/${entry_point}.dxil)
			add_custom_command(OUTPUT ${dxil_output}
				COMMAND ${CMAKE_COMMAND} -E make_directory ${TFX_SDL_SHADER_FOLDER}
				COMMAND ${TFX_SLANGC} ${TFX_SHADER_SOURCE} -entry ${entry_point} -target dxil -profile ${profile} ${TFX_SDL_SHADER_DEFINES} -D TFX_ENTRY_${entry_point} -o ${dxil_output}
				DEPENDS ${TFX_SHADER_SOURCE}
				COMMENT "Compiling ${entry_point} for SDL_gpu D3D12")
			list(APPEND TFX_SDL_SHADERS ${dxil_output})
		endif()
	endforeach()
	add_custom_target(tfx_sdl_shaders DEPENDS ${TFX_SDL_SHADERS})

	add_executable(example_sdl3_gpu
		${CMAKE_CURRENT_SOURCE_DIR}/example_sdl3_gpu/main.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.h
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_sdlgpu.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_sdlgpu.h
		${TFX_SHADER_SOURCE}
		${TFX_IMPL_COMMON_SOURCE}
		${TFX_LIBRARY_SOURCE}
	)
	set_source_files_properties(${TFX_SHADER_SOURCE} PROPERTIES VS_TOOL_OVERRIDE "None")
	source_group("Shaders" FILES ${TFX_SHADER_SOURCE})
	tfx_example_target(example_sdl3_gpu)
	# vulkan_core.h from SDL's own copy, because the Vulkan SDK's include folder has older SDL3 headers in it
	target_include_directories(example_sdl3_gpu PRIVATE ${sdl3_SOURCE_DIR}/src/video/khronos)
	target_link_libraries(example_sdl3_gpu PRIVATE SDL3::SDL3)
	add_dependencies(example_sdl3_gpu tfx_sdl_shaders)
	target_compile_definitions(example_sdl3_gpu PRIVATE TFX_EXAMPLE_SHADER_FOLDER="${TFX_SDL_SHADER_FOLDER}")
endif()
