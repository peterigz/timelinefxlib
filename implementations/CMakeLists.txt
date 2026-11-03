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

# The OpenGL and Vulkan examples use the parent project's glfw target when it has one, otherwise they download GLFW
macro(tfx_example_require_glfw)
	if (NOT TARGET glfw)
		include(FetchContent)
		set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
		set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
		set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
		set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
		FetchContent_Declare(glfw
			GIT_REPOSITORY https://github.com/glfw/glfw.git
			GIT_TAG 3.4
			GIT_SHALLOW TRUE)
		FetchContent_MakeAvailable(glfw)
	endif()
endmacro()

# Off unless asked for, because it needs GLFW
option(TFX_EXAMPLE_GL "Build the OpenGL example with GLFW" OFF)
if (TFX_EXAMPLE_GL)
	find_program(TFX_SPIRV_CROSS spirv-cross HINTS "$ENV{VULKAN_SDK}/Bin")
	if (NOT TFX_SPIRV_CROSS)
		message(FATAL_ERROR "spirv-cross wasn't found. It comes with the Vulkan SDK, put its bin folder on the PATH or set TFX_SPIRV_CROSS to spirv-cross")
	endif()
	tfx_example_require_glfw()

	# Slang's own GLSL is Vulkan GLSL, so each entry point goes through SPIR-V and spirv-cross to GLSL 430, with the default fixed bindings
	set(TFX_GL_SHADER_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/gl_shaders)
	set(TFX_GL_ENTRY_POINTS
		simulate_particles
		tfx_billboard_vertex
		tfx_billboard_baked_vertex
		tfx_billboard_fragment
		tfx_ribbon_compute
		tfx_ribbon_vertex
		tfx_ribbon_fragment
		tfx_sprite_playback
		tfx_ribbon_playback
	)
	set(TFX_GL_SHADERS)
	foreach(entry_point ${TFX_GL_ENTRY_POINTS})
		set(spirv_output ${TFX_GL_SHADER_FOLDER}/${entry_point}.spv)
		set(glsl_output ${TFX_GL_SHADER_FOLDER}/${entry_point}.glsl)
		add_custom_command(OUTPUT ${glsl_output}
			COMMAND ${CMAKE_COMMAND} -E make_directory ${TFX_GL_SHADER_FOLDER}
			COMMAND ${TFX_SLANGC} ${TFX_SHADER_SOURCE} -entry ${entry_point} -target spirv -profile spirv_1_5 -O3 -D TFX_FIXED_BINDINGS -D TFX_COMBINED_SAMPLERS -D TFX_ENTRY_${entry_point} -o ${spirv_output}
			COMMAND ${TFX_SPIRV_CROSS} ${spirv_output} --version 430 --no-es --no-support-nonzero-baseinstance --output ${glsl_output}
			DEPENDS ${TFX_SHADER_SOURCE}
			COMMENT "Compiling ${entry_point} for OpenGL")
		list(APPEND TFX_GL_SHADERS ${glsl_output})
	endforeach()
	add_custom_target(tfx_gl_shaders DEPENDS ${TFX_GL_SHADERS})

	add_executable(example_glfw_opengl
		${CMAKE_CURRENT_SOURCE_DIR}/example_glfw_opengl/main.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.h
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_gl.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_gl.h
		${TFX_SHADER_SOURCE}
		${TFX_IMPL_COMMON_SOURCE}
		${TFX_LIBRARY_SOURCE}
	)
	set_source_files_properties(${TFX_SHADER_SOURCE} PROPERTIES VS_TOOL_OVERRIDE "None")
	source_group("Shaders" FILES ${TFX_SHADER_SOURCE})
	tfx_example_target(example_glfw_opengl)
	target_compile_definitions(example_glfw_opengl PRIVATE TFX_EXAMPLE_SHADER_FOLDER="${TFX_GL_SHADER_FOLDER}")
	target_link_libraries(example_glfw_opengl PRIVATE glfw)
	add_dependencies(example_glfw_opengl tfx_gl_shaders)
endif()

# Off unless asked for, to keep the build short. The bindless model: DXIL with no defines, each entry point reading every resource
# from one descriptor heap by the index in its root constants
option(TFX_EXAMPLE_DX12 "Build the Direct3D 12 example, with the bindless shaders" OFF)
if (WIN32 AND TFX_EXAMPLE_DX12)
	set(TFX_DX12_SHADER_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/dx12_shaders)
	set(TFX_DX12_ENTRY_POINTS
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
	set(TFX_DX12_SHADERS)
	list(LENGTH TFX_DX12_ENTRY_POINTS entry_point_list_length)
	math(EXPR last_entry_point_index "${entry_point_list_length} - 1")
	foreach(entry_point_index RANGE 0 ${last_entry_point_index} 2)
		math(EXPR profile_index "${entry_point_index} + 1")
		list(GET TFX_DX12_ENTRY_POINTS ${entry_point_index} entry_point)
		list(GET TFX_DX12_ENTRY_POINTS ${profile_index} profile)
		set(shader_output ${TFX_DX12_SHADER_FOLDER}/${entry_point}.dxil)
		add_custom_command(OUTPUT ${shader_output}
			COMMAND ${CMAKE_COMMAND} -E make_directory ${TFX_DX12_SHADER_FOLDER}
			COMMAND ${TFX_SLANGC} ${TFX_SHADER_SOURCE} -entry ${entry_point} -target dxil -profile ${profile} -o ${shader_output}
			DEPENDS ${TFX_SHADER_SOURCE}
			COMMENT "Compiling ${entry_point} for D3D12")
		list(APPEND TFX_DX12_SHADERS ${shader_output})
	endforeach()
	add_custom_target(tfx_dx12_shaders DEPENDS ${TFX_DX12_SHADERS})

	add_executable(example_win32_dx12
		${CMAKE_CURRENT_SOURCE_DIR}/example_win32_dx12/main.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.h
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_dx12.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_dx12.h
		${TFX_SHADER_SOURCE}
		${TFX_IMPL_COMMON_SOURCE}
		${TFX_LIBRARY_SOURCE}
	)
	set_source_files_properties(${TFX_SHADER_SOURCE} PROPERTIES VS_TOOL_OVERRIDE "None")
	source_group("Shaders" FILES ${TFX_SHADER_SOURCE})
	tfx_example_target(example_win32_dx12)
	target_link_libraries(example_win32_dx12 PRIVATE d3d12 dxgi)
	add_dependencies(example_win32_dx12 tfx_dx12_shaders)
	target_compile_definitions(example_win32_dx12 PRIVATE TFX_EXAMPLE_SHADER_FOLDER="${TFX_DX12_SHADER_FOLDER}")
endif()

# Off unless asked for, because it needs the Vulkan SDK and GLFW. The bindless model: SPIR-V with no defines, each entry point keeping
# its name, reading every resource from one descriptor set by the index in its push constants
option(TFX_EXAMPLE_VULKAN "Build the Vulkan example with GLFW, with the bindless shaders" OFF)
if (TFX_EXAMPLE_VULKAN)
	find_package(Vulkan REQUIRED)
	tfx_example_require_glfw()
	set(TFX_VULKAN_SHADER_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/vulkan_shaders)
	set(TFX_VULKAN_ENTRY_POINTS
		simulate_particles
		tfx_billboard_vertex
		tfx_billboard_baked_vertex
		tfx_billboard_fragment
		tfx_ribbon_compute
		tfx_ribbon_vertex
		tfx_ribbon_fragment
		tfx_sprite_playback
		tfx_ribbon_playback
	)
	set(TFX_VULKAN_SHADERS)
	foreach(entry_point ${TFX_VULKAN_ENTRY_POINTS})
		set(shader_output ${TFX_VULKAN_SHADER_FOLDER}/${entry_point}.spv)
		add_custom_command(OUTPUT ${shader_output}
			COMMAND ${CMAKE_COMMAND} -E make_directory ${TFX_VULKAN_SHADER_FOLDER}
			COMMAND ${TFX_SLANGC} ${TFX_SHADER_SOURCE} -entry ${entry_point} -target spirv -profile spirv_1_5 -O3 -fvk-use-entrypoint-name -o ${shader_output}
			DEPENDS ${TFX_SHADER_SOURCE}
			COMMENT "Compiling ${entry_point} for Vulkan")
		list(APPEND TFX_VULKAN_SHADERS ${shader_output})
	endforeach()
	add_custom_target(tfx_vulkan_shaders DEPENDS ${TFX_VULKAN_SHADERS})

	add_executable(example_glfw_vulkan
		${CMAKE_CURRENT_SOURCE_DIR}/example_glfw_vulkan/main.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.h
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_vulkan.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_vulkan.h
		${TFX_SHADER_SOURCE}
		${TFX_IMPL_COMMON_SOURCE}
		${TFX_LIBRARY_SOURCE}
	)
	set_source_files_properties(${TFX_SHADER_SOURCE} PROPERTIES VS_TOOL_OVERRIDE "None")
	source_group("Shaders" FILES ${TFX_SHADER_SOURCE})
	tfx_example_target(example_glfw_vulkan)
	target_link_libraries(example_glfw_vulkan PRIVATE glfw Vulkan::Vulkan)
	add_dependencies(example_glfw_vulkan tfx_vulkan_shaders)
	target_compile_definitions(example_glfw_vulkan PRIVATE TFX_EXAMPLE_SHADER_FOLDER="${TFX_VULKAN_SHADER_FOLDER}")
endif()

# Off unless asked for, because it downloads wgpu-native's prebuilt release and needs GLFW. The fixed bindings as WGSL, each entry
# point keeping its own name
option(TFX_EXAMPLE_WEBGPU "Build the WebGPU example with wgpu-native and GLFW" OFF)
if (TFX_EXAMPLE_WEBGPU)
	tfx_example_require_glfw()
	set(TFX_WGPU_NATIVE_DIR "" CACHE PATH "An unzipped wgpu-native release holding include/webgpu and lib, downloaded when empty")
	set(wgpu_native_folder ${TFX_WGPU_NATIVE_DIR})
	if (NOT wgpu_native_folder)
		if (WIN32 AND CMAKE_SIZEOF_VOID_P EQUAL 8)
			set(wgpu_native_archive wgpu-windows-x86_64-msvc-release.zip)
			set(wgpu_native_hash 7e67d7445c42aeb85e30f88930fd8d7d83ee769e3390aeb1ada75ebf3cf78132)
		elseif (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
			set(wgpu_native_archive wgpu-linux-x86_64-release.zip)
			set(wgpu_native_hash 95a4d90c071005a98d03eab348beaa6b07e16eb00d1dcdb9f8348f75eb97ec5a)
		else()
			message(FATAL_ERROR "There's no wgpu-native download set up for this platform. Unzip a release from github.com/gfx-rs/wgpu-native and set TFX_WGPU_NATIVE_DIR to it")
		endif()
		if (POLICY CMP0135)
			cmake_policy(SET CMP0135 NEW)
		endif()
		include(FetchContent)
		FetchContent_Declare(wgpu_native
			URL https://github.com/gfx-rs/wgpu-native/releases/download/v29.0.1.1/${wgpu_native_archive}
			URL_HASH SHA256=${wgpu_native_hash})
		FetchContent_MakeAvailable(wgpu_native)
		set(wgpu_native_folder ${wgpu_native_SOURCE_DIR})
	endif()
	if (NOT TARGET tfx_wgpu_native)
		add_library(tfx_wgpu_native SHARED IMPORTED)
		if (WIN32)
			set_target_properties(tfx_wgpu_native PROPERTIES
				IMPORTED_LOCATION ${wgpu_native_folder}/lib/wgpu_native.dll
				IMPORTED_IMPLIB ${wgpu_native_folder}/lib/wgpu_native.dll.lib)
		elseif (APPLE)
			set_target_properties(tfx_wgpu_native PROPERTIES IMPORTED_LOCATION ${wgpu_native_folder}/lib/libwgpu_native.dylib)
		else()
			set_target_properties(tfx_wgpu_native PROPERTIES IMPORTED_LOCATION ${wgpu_native_folder}/lib/libwgpu_native.so)
		endif()
		set_target_properties(tfx_wgpu_native PROPERTIES
			INTERFACE_INCLUDE_DIRECTORIES ${wgpu_native_folder}/include
			INTERFACE_COMPILE_DEFINITIONS WGPU_SHARED_LIBRARY)
	endif()

	set(TFX_WEBGPU_SHADER_FOLDER ${CMAKE_CURRENT_BINARY_DIR}/webgpu_shaders)
	set(TFX_WEBGPU_ENTRY_POINTS
		simulate_particles
		tfx_billboard_vertex
		tfx_billboard_baked_vertex
		tfx_billboard_fragment
		tfx_ribbon_compute
		tfx_ribbon_vertex
		tfx_ribbon_fragment
		tfx_sprite_playback
		tfx_ribbon_playback
	)
	set(TFX_WEBGPU_SHADERS)
	foreach(entry_point ${TFX_WEBGPU_ENTRY_POINTS})
		set(shader_output ${TFX_WEBGPU_SHADER_FOLDER}/${entry_point}.wgsl)
		add_custom_command(OUTPUT ${shader_output}
			COMMAND ${CMAKE_COMMAND} -E make_directory ${TFX_WEBGPU_SHADER_FOLDER}
			COMMAND ${TFX_SLANGC} ${TFX_SHADER_SOURCE} -entry ${entry_point} -target wgsl -D TFX_FIXED_BINDINGS -D TFX_ENTRY_${entry_point} -o ${shader_output}
			DEPENDS ${TFX_SHADER_SOURCE}
			COMMENT "Compiling ${entry_point} for WebGPU")
		list(APPEND TFX_WEBGPU_SHADERS ${shader_output})
	endforeach()
	add_custom_target(tfx_webgpu_shaders DEPENDS ${TFX_WEBGPU_SHADERS})

	add_executable(example_glfw_webgpu
		${CMAKE_CURRENT_SOURCE_DIR}/example_glfw_webgpu/main.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.h
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_webgpu.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_webgpu.h
		${TFX_SHADER_SOURCE}
		${TFX_IMPL_COMMON_SOURCE}
		${TFX_LIBRARY_SOURCE}
	)
	set_source_files_properties(${TFX_SHADER_SOURCE} PROPERTIES VS_TOOL_OVERRIDE "None")
	source_group("Shaders" FILES ${TFX_SHADER_SOURCE})
	tfx_example_target(example_glfw_webgpu)
	target_link_libraries(example_glfw_webgpu PRIVATE glfw tfx_wgpu_native)
	add_dependencies(example_glfw_webgpu tfx_webgpu_shaders)
	target_compile_definitions(example_glfw_webgpu PRIVATE TFX_EXAMPLE_SHADER_FOLDER="${TFX_WEBGPU_SHADER_FOLDER}")
	if (WIN32)
		add_custom_command(TARGET example_glfw_webgpu POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_if_different ${wgpu_native_folder}/lib/wgpu_native.dll $<TARGET_FILE_DIR:example_glfw_webgpu>)
	endif()
endif()
