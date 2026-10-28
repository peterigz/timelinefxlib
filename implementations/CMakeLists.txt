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

# The GPU particle kernel is built from the library's Slang source. Without slangc the examples run without GPU simulation
find_program(TFX_SLANGC slangc HINTS "$ENV{VULKAN_SDK}/Bin")
set(TFX_KERNEL_SOURCE ${TFX_LIBRARY_FOLDER}/Shaders/tfx_gpu_particles.slang)
if (NOT TFX_SLANGC)
	message(WARNING "slangc wasn't found (it comes with the Vulkan SDK), so the examples are built without GPU particle simulation")
endif()

if (WIN32)
	if (TFX_SLANGC)
		set(TFX_KERNEL_DX11 ${CMAKE_CURRENT_BINARY_DIR}/tfx_gpu_particles_cs_5_0.cso)
		# fxc can't compile Slang's #line directives, hence -line-directive-mode none
		add_custom_command(OUTPUT ${TFX_KERNEL_DX11}
			COMMAND ${TFX_SLANGC} ${TFX_KERNEL_SOURCE} -entry simulate_particles -target dxbc -profile cs_5_0 -D TFX_FIXED_BINDINGS -line-directive-mode none -o ${TFX_KERNEL_DX11}
			DEPENDS ${TFX_KERNEL_SOURCE}
			COMMENT "Compiling the GPU particle kernel for D3D11")
		add_custom_target(tfx_gpu_kernel_dx11 DEPENDS ${TFX_KERNEL_DX11})
	endif()

	add_executable(example_win32_dx11
		${CMAKE_CURRENT_SOURCE_DIR}/example_win32_dx11/main.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/example_scene.h
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_dx11.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/tfx_impl_dx11.h
		${CMAKE_CURRENT_SOURCE_DIR}/shaders/tfx_dx11.hlsl
		${TFX_IMPL_COMMON_SOURCE}
		${TFX_LIBRARY_SOURCE}
	)
	# Compiled at runtime with D3DCompileFromFile, so Visual Studio mustn't try to build it
	set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/shaders/tfx_dx11.hlsl PROPERTIES VS_TOOL_OVERRIDE "None")
	tfx_example_target(example_win32_dx11)
	target_link_libraries(example_win32_dx11 PRIVATE d3d11 dxgi d3dcompiler)
	if (TFX_SLANGC)
		add_dependencies(example_win32_dx11 tfx_gpu_kernel_dx11)
		target_compile_definitions(example_win32_dx11 PRIVATE TFX_EXAMPLE_KERNEL_PATH="${TFX_KERNEL_DX11}")
	endif()
endif()
