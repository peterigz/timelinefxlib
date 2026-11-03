//TimelineFX example: GLFW + Vulkan 1.3 with the bindless shaders. The window, instance, device, swap chain and main loop, see ../README.md

#include <vulkan/vulkan.h>
#ifndef GLFW_INCLUDE_NONE
#define GLFW_INCLUDE_NONE
#endif
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "timelinefx.h"
#include "tfx_impl_common.h"
#include "tfx_impl_vulkan.h"
#include "example_scene.h"

//TFX_EXAMPLE_ROOT is set by the CMake project to the implementations folder
#ifndef TFX_EXAMPLE_ROOT
#define TFX_EXAMPLE_ROOT "."
#endif
#define EXAMPLE_LIBRARY_PATH TFX_EXAMPLE_ROOT "/example.tfx"
#define EXAMPLE_SECOND_LIBRARY_PATH TFX_EXAMPLE_ROOT "/example_second.tfx"
#define EXAMPLE_SPRITE_DATA_PATH TFX_EXAMPLE_ROOT "/example.tfxsd"
//TFX_EXAMPLE_SHADER_FOLDER is where the CMake project has slangc build the shaders
#ifndef TFX_EXAMPLE_SHADER_FOLDER
#define TFX_EXAMPLE_SHADER_FOLDER "vulkan_shaders"
#endif
#define EXAMPLE_FEATURES (example_scene_feature_ribbons | example_scene_feature_gpu_simulation | example_scene_feature_baked)

#define EXAMPLE_WINDOW_WIDTH 1280
#define EXAMPLE_WINDOW_HEIGHT 800
#define EXAMPLE_MEMORY_POOL_SIZE (64ull * 1024ull * 1024ull)
#define EXAMPLE_FRAME_COUNT 2
#define EXAMPLE_MAX_SWAPCHAIN_IMAGES 8
#define EXAMPLE_DEPTH_FORMAT VK_FORMAT_D32_SFLOAT
#define EXAMPLE_VALIDATION_LAYER "VK_LAYER_KHRONOS_validation"

//The features tfx_impl_vulkan_Init lists, chained for vkGetPhysicalDeviceFeatures2 and vkCreateDevice. Never copied, it points into itself
typedef struct example_features_s {
	VkPhysicalDeviceFeatures2 features;
	VkPhysicalDeviceVulkan11Features vulkan_11;
	VkPhysicalDeviceVulkan12Features vulkan_12;
	VkPhysicalDeviceVulkan13Features vulkan_13;
} example_features_t;

typedef struct example_s {
	GLFWwindow *window;
	VkInstance instance;
	VkDebugUtilsMessengerEXT messenger;		//Null without the validation layer
	VkSurfaceKHR surface;
	VkPhysicalDevice physical_device;
	VkDevice device;
	uint32_t queue_family;
	VkQueue queue;
	VkSwapchainKHR swapchain;
	VkFormat color_format;
	VkExtent2D extent;
	uint32_t image_count;
	VkImage swapchain_images[EXAMPLE_MAX_SWAPCHAIN_IMAGES];
	VkImageView swapchain_views[EXAMPLE_MAX_SWAPCHAIN_IMAGES];
	//One per swap chain image, because presenting holds on to it until that image is acquired again
	VkSemaphore render_finished[EXAMPLE_MAX_SWAPCHAIN_IMAGES];
	VkImage depth_image;
	VkDeviceMemory depth_memory;
	VkImageView depth_view;
	VkCommandPool command_pools[EXAMPLE_FRAME_COUNT];
	VkCommandBuffer command_buffers[EXAMPLE_FRAME_COUNT];
	VkFence frame_fences[EXAMPLE_FRAME_COUNT];
	VkSemaphore image_available[EXAMPLE_FRAME_COUNT];
	bool swapchain_stale;

	tfx_impl_vulkan_t renderer;
	example_scene_t scene;
} example_t;

static example_t example;

static void ErrorCallback(int error, const char *description) {
	printf("GLFW error %i: %s\n", error, description);
}

static void KeyCallback(GLFWwindow *window, int key, int scancode, int action, int mods) {
	(void)scancode; (void)mods;
	if (action != GLFW_PRESS) {
		return;
	}
	if (key == GLFW_KEY_ESCAPE) {
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	} else if (key == GLFW_KEY_SPACE) {
		example_scene_Burst(&example.scene);
	} else if (key == GLFW_KEY_O) {
		example_scene_ToggleOrbit(&example.scene);
	} else if (key == GLFW_KEY_L) {
		example_scene_ReloadSecondLibrary(&example.scene);
	}
}

static void ScrollCallback(GLFWwindow *window, double x_offset, double y_offset) {
	(void)window; (void)x_offset;
	example_scene_Zoom(&example.scene, (float)y_offset);
}

static VKAPI_ATTR VkBool32 VKAPI_CALL DebugCallback(VkDebugUtilsMessageSeverityFlagBitsEXT severity, VkDebugUtilsMessageTypeFlagsEXT types,
	const VkDebugUtilsMessengerCallbackDataEXT *data, void *user_data) {
	(void)types; (void)user_data;
	printf("Vulkan %s: %s\n", severity >= VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT ? "error" : "warning", data->pMessage);
	return VK_FALSE;
}

static bool HasValidationLayer(void) {
	uint32_t layer_count = 0;
	vkEnumerateInstanceLayerProperties(&layer_count, NULL);
	VkLayerProperties *layers = (VkLayerProperties *)malloc(sizeof(VkLayerProperties) * (layer_count ? layer_count : 1));
	bool found = false;
	if (layers && vkEnumerateInstanceLayerProperties(&layer_count, layers) == VK_SUCCESS) {
		for (uint32_t index = 0; index != layer_count && !found; ++index) {
			found = strcmp(layers[index].layerName, EXAMPLE_VALIDATION_LAYER) == 0;
		}
	}
	free(layers);
	return found;
}

//The validation layer in debug builds when it's installed (it comes with the Vulkan SDK). GPU assisted validation checks what the shaders
//do with each descriptor, which the layer can't see from the API calls. It's slow
static bool CreateInstance(example_t *state, bool gpu_validation) {
	bool validation = false;
#ifndef NDEBUG
	validation = HasValidationLayer();
	if (!validation) {
		printf("The Vulkan validation layer isn't installed, continuing without it\n");
	}
#endif
	if (gpu_validation && !validation) {
		printf("GPU assisted validation needs the validation layer, so it's off\n");
		gpu_validation = false;
	}
	uint32_t glfw_extension_count = 0;
	const char **glfw_extensions = glfwGetRequiredInstanceExtensions(&glfw_extension_count);
	if (!glfw_extensions) {
		printf("GLFW can't find a Vulkan loader with surface support\n");
		return false;
	}
	const char *extensions[16];
	uint32_t extension_count = 0;
	for (uint32_t index = 0; index != glfw_extension_count && extension_count < 14; ++index) {
		extensions[extension_count++] = glfw_extensions[index];
	}
	if (validation) {
		extensions[extension_count++] = VK_EXT_DEBUG_UTILS_EXTENSION_NAME;
	}
	if (gpu_validation) {
		extensions[extension_count++] = VK_EXT_VALIDATION_FEATURES_EXTENSION_NAME;
	}
	const char *layers[] = { EXAMPLE_VALIDATION_LAYER };
	VkValidationFeatureEnableEXT validation_enables[] = { VK_VALIDATION_FEATURE_ENABLE_GPU_ASSISTED_EXT, VK_VALIDATION_FEATURE_ENABLE_GPU_ASSISTED_RESERVE_BINDING_SLOT_EXT };
	VkValidationFeaturesEXT validation_features = { VK_STRUCTURE_TYPE_VALIDATION_FEATURES_EXT };
	validation_features.enabledValidationFeatureCount = 2;
	validation_features.pEnabledValidationFeatures = validation_enables;

	VkApplicationInfo application = { VK_STRUCTURE_TYPE_APPLICATION_INFO };
	application.pApplicationName = "TimelineFX Vulkan Example";
	application.pEngineName = "TimelineFX";
	application.apiVersion = VK_API_VERSION_1_3;
	VkInstanceCreateInfo create_info = { VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO };
	create_info.pNext = gpu_validation ? &validation_features : NULL;
	create_info.pApplicationInfo = &application;
	create_info.enabledLayerCount = validation ? 1 : 0;
	create_info.ppEnabledLayerNames = layers;
	create_info.enabledExtensionCount = extension_count;
	create_info.ppEnabledExtensionNames = extensions;
	VkResult result = vkCreateInstance(&create_info, NULL, &state->instance);
	if (result != VK_SUCCESS) {
		printf("Could not create a Vulkan 1.3 instance (%i)\n", (int)result);
		return false;
	}
	if (validation) {
		VkDebugUtilsMessengerCreateInfoEXT messenger_info = { VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT };
		messenger_info.messageSeverity = VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT;
		messenger_info.messageType = VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT;
		messenger_info.pfnUserCallback = DebugCallback;
		PFN_vkCreateDebugUtilsMessengerEXT create_messenger = (PFN_vkCreateDebugUtilsMessengerEXT)vkGetInstanceProcAddr(state->instance, "vkCreateDebugUtilsMessengerEXT");
		if (!create_messenger || create_messenger(state->instance, &messenger_info, NULL, &state->messenger) != VK_SUCCESS) {
			printf("Could not create the debug messenger, the validation layer's messages won't be printed\n");
		}
	}
	return true;
}

static void ChainFeatures(example_features_t *chain) {
	memset(chain, 0, sizeof(example_features_t));
	chain->features.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2;
	chain->features.pNext = &chain->vulkan_11;
	chain->vulkan_11.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_1_FEATURES;
	chain->vulkan_11.pNext = &chain->vulkan_12;
	chain->vulkan_12.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES;
	chain->vulkan_12.pNext = &chain->vulkan_13;
	chain->vulkan_13.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_3_FEATURES;
}

static bool HasRequiredFeatures(const example_features_t *chain) {
	return chain->vulkan_11.shaderDrawParameters && chain->vulkan_12.runtimeDescriptorArray && chain->vulkan_12.descriptorBindingPartiallyBound
		&& chain->vulkan_12.descriptorBindingSampledImageUpdateAfterBind && chain->vulkan_12.descriptorBindingStorageBufferUpdateAfterBind
		&& chain->vulkan_12.descriptorBindingUpdateUnusedWhilePending && chain->vulkan_12.shaderSampledImageArrayNonUniformIndexing
		&& chain->vulkan_13.dynamicRendering && chain->vulkan_13.synchronization2;
}

static void RequireFeatures(example_features_t *chain) {
	ChainFeatures(chain);
	chain->vulkan_11.shaderDrawParameters = VK_TRUE;
	chain->vulkan_12.runtimeDescriptorArray = VK_TRUE;
	chain->vulkan_12.descriptorBindingPartiallyBound = VK_TRUE;
	chain->vulkan_12.descriptorBindingSampledImageUpdateAfterBind = VK_TRUE;
	chain->vulkan_12.descriptorBindingStorageBufferUpdateAfterBind = VK_TRUE;
	chain->vulkan_12.descriptorBindingUpdateUnusedWhilePending = VK_TRUE;
	chain->vulkan_12.shaderSampledImageArrayNonUniformIndexing = VK_TRUE;
	chain->vulkan_13.dynamicRendering = VK_TRUE;
	chain->vulkan_13.synchronization2 = VK_TRUE;
}

//A Vulkan 1.3 device with the features and one queue that draws, computes and presents. Discrete GPUs first
static bool PickPhysicalDevice(example_t *state) {
	uint32_t device_count = 0;
	vkEnumeratePhysicalDevices(state->instance, &device_count, NULL);
	VkPhysicalDevice devices[16];
	device_count = device_count < 16 ? device_count : 16;
	vkEnumeratePhysicalDevices(state->instance, &device_count, devices);
	int best_score = 0;
	for (uint32_t device_index = 0; device_index != device_count; ++device_index) {
		VkPhysicalDeviceProperties properties;
		vkGetPhysicalDeviceProperties(devices[device_index], &properties);
		example_features_t chain;
		ChainFeatures(&chain);
		vkGetPhysicalDeviceFeatures2(devices[device_index], &chain.features);
		if (properties.apiVersion < VK_API_VERSION_1_3 || !HasRequiredFeatures(&chain)) {
			continue;
		}
		uint32_t family_count = 0;
		vkGetPhysicalDeviceQueueFamilyProperties(devices[device_index], &family_count, NULL);
		VkQueueFamilyProperties families[16];
		family_count = family_count < 16 ? family_count : 16;
		vkGetPhysicalDeviceQueueFamilyProperties(devices[device_index], &family_count, families);
		for (uint32_t family = 0; family != family_count; ++family) {
			VkBool32 present = VK_FALSE;
			vkGetPhysicalDeviceSurfaceSupportKHR(devices[device_index], family, state->surface, &present);
			VkQueueFlags needed = VK_QUEUE_GRAPHICS_BIT | VK_QUEUE_COMPUTE_BIT;
			if (present && (families[family].queueFlags & needed) == needed) {
				int score = properties.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU ? 2 : 1;
				if (score > best_score) {
					best_score = score;
					state->physical_device = devices[device_index];
					state->queue_family = family;
				}
				break;
			}
		}
	}
	if (!state->physical_device) {
		printf("No device has Vulkan 1.3, the descriptor indexing features and a queue that draws, computes and presents\n");
		return false;
	}
	VkPhysicalDeviceProperties properties;
	vkGetPhysicalDeviceProperties(state->physical_device, &properties);
	printf("%s, Vulkan %u.%u.%u\n", properties.deviceName, VK_API_VERSION_MAJOR(properties.apiVersion), VK_API_VERSION_MINOR(properties.apiVersion), VK_API_VERSION_PATCH(properties.apiVersion));
	return true;
}

static bool CreateDevice(example_t *state) {
	float priority = 1.f;
	VkDeviceQueueCreateInfo queue_info = { VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO };
	queue_info.queueFamilyIndex = state->queue_family;
	queue_info.queueCount = 1;
	queue_info.pQueuePriorities = &priority;
	example_features_t chain;
	RequireFeatures(&chain);
	const char *extensions[] = { VK_KHR_SWAPCHAIN_EXTENSION_NAME };
	VkDeviceCreateInfo create_info = { VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO };
	create_info.pNext = &chain.features;
	create_info.queueCreateInfoCount = 1;
	create_info.pQueueCreateInfos = &queue_info;
	create_info.enabledExtensionCount = 1;
	create_info.ppEnabledExtensionNames = extensions;
	VkResult result = vkCreateDevice(state->physical_device, &create_info, NULL, &state->device);
	if (result != VK_SUCCESS) {
		printf("Could not create the Vulkan device (%i)\n", (int)result);
		return false;
	}
	vkGetDeviceQueue(state->device, state->queue_family, 0, &state->queue);
	for (uint32_t frame = 0; frame != EXAMPLE_FRAME_COUNT; ++frame) {
		VkCommandPoolCreateInfo pool_info = { VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO };
		pool_info.flags = VK_COMMAND_POOL_CREATE_TRANSIENT_BIT;
		pool_info.queueFamilyIndex = state->queue_family;
		VkCommandBufferAllocateInfo allocate_info = { VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO };
		allocate_info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
		allocate_info.commandBufferCount = 1;
		//Signalled, so the first wait on each frame slot returns straight away
		VkFenceCreateInfo fence_info = { VK_STRUCTURE_TYPE_FENCE_CREATE_INFO };
		fence_info.flags = VK_FENCE_CREATE_SIGNALED_BIT;
		VkSemaphoreCreateInfo semaphore_info = { VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO };
		if (vkCreateCommandPool(state->device, &pool_info, NULL, &state->command_pools[frame]) != VK_SUCCESS) {
			return false;
		}
		allocate_info.commandPool = state->command_pools[frame];
		if (vkAllocateCommandBuffers(state->device, &allocate_info, &state->command_buffers[frame]) != VK_SUCCESS
			|| vkCreateFence(state->device, &fence_info, NULL, &state->frame_fences[frame]) != VK_SUCCESS
			|| vkCreateSemaphore(state->device, &semaphore_info, NULL, &state->image_available[frame]) != VK_SUCCESS) {
			printf("Could not create the command buffers and frame synchronisation\n");
			return false;
		}
	}
	return true;
}

static uint32_t FindMemoryType(example_t *state, uint32_t type_bits, VkMemoryPropertyFlags properties) {
	VkPhysicalDeviceMemoryProperties memory_properties;
	vkGetPhysicalDeviceMemoryProperties(state->physical_device, &memory_properties);
	for (uint32_t type = 0; type != memory_properties.memoryTypeCount; ++type) {
		if ((type_bits & (1u << type)) && (memory_properties.memoryTypes[type].propertyFlags & properties) == properties) {
			return type;
		}
	}
	return UINT32_MAX;
}

static void ReleaseSwapchainImages(example_t *state) {
	for (uint32_t index = 0; index != state->image_count; ++index) {
		vkDestroyImageView(state->device, state->swapchain_views[index], NULL);
		vkDestroySemaphore(state->device, state->render_finished[index], NULL);
		state->swapchain_views[index] = VK_NULL_HANDLE;
		state->render_finished[index] = VK_NULL_HANDLE;
	}
	state->image_count = 0;
	vkDestroyImageView(state->device, state->depth_view, NULL);
	vkDestroyImage(state->device, state->depth_image, NULL);
	vkFreeMemory(state->device, state->depth_memory, NULL);
	state->depth_view = VK_NULL_HANDLE;
	state->depth_image = VK_NULL_HANDLE;
	state->depth_memory = VK_NULL_HANDLE;
}

static bool CreateDepthBuffer(example_t *state) {
	VkImageCreateInfo image_info = { VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO };
	image_info.imageType = VK_IMAGE_TYPE_2D;
	image_info.format = EXAMPLE_DEPTH_FORMAT;
	image_info.extent.width = state->extent.width;
	image_info.extent.height = state->extent.height;
	image_info.extent.depth = 1;
	image_info.mipLevels = 1;
	image_info.arrayLayers = 1;
	image_info.samples = VK_SAMPLE_COUNT_1_BIT;
	image_info.tiling = VK_IMAGE_TILING_OPTIMAL;
	image_info.usage = VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT;
	image_info.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
	if (vkCreateImage(state->device, &image_info, NULL, &state->depth_image) != VK_SUCCESS) {
		return false;
	}
	VkMemoryRequirements requirements;
	vkGetImageMemoryRequirements(state->device, state->depth_image, &requirements);
	VkMemoryAllocateInfo allocate_info = { VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO };
	allocate_info.allocationSize = requirements.size;
	allocate_info.memoryTypeIndex = FindMemoryType(state, requirements.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
	if (allocate_info.memoryTypeIndex == UINT32_MAX || vkAllocateMemory(state->device, &allocate_info, NULL, &state->depth_memory) != VK_SUCCESS
		|| vkBindImageMemory(state->device, state->depth_image, state->depth_memory, 0) != VK_SUCCESS) {
		return false;
	}
	VkImageViewCreateInfo view_info = { VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO };
	view_info.image = state->depth_image;
	view_info.viewType = VK_IMAGE_VIEW_TYPE_2D;
	view_info.format = EXAMPLE_DEPTH_FORMAT;
	view_info.subresourceRange.aspectMask = VK_IMAGE_ASPECT_DEPTH_BIT;
	view_info.subresourceRange.levelCount = 1;
	view_info.subresourceRange.layerCount = 1;
	return vkCreateImageView(state->device, &view_info, NULL, &state->depth_view) == VK_SUCCESS;
}

//UNORM rather than SRGB, so it blends and looks the same as the D3D examples
static VkSurfaceFormatKHR ChooseSurfaceFormat(example_t *state) {
	uint32_t format_count = 0;
	vkGetPhysicalDeviceSurfaceFormatsKHR(state->physical_device, state->surface, &format_count, NULL);
	VkSurfaceFormatKHR formats[64];
	format_count = format_count < 64 ? format_count : 64;
	vkGetPhysicalDeviceSurfaceFormatsKHR(state->physical_device, state->surface, &format_count, formats);
	for (uint32_t index = 0; index != format_count; ++index) {
		if (formats[index].format == VK_FORMAT_B8G8R8A8_UNORM || formats[index].format == VK_FORMAT_R8G8B8A8_UNORM) {
			return formats[index];
		}
	}
	return formats[0];
}

static bool CreateSwapchain(example_t *state) {
	VkSurfaceCapabilitiesKHR capabilities;
	vkGetPhysicalDeviceSurfaceCapabilitiesKHR(state->physical_device, state->surface, &capabilities);
	VkSurfaceFormatKHR surface_format = ChooseSurfaceFormat(state);
	if (state->color_format != VK_FORMAT_UNDEFINED && surface_format.format != state->color_format) {
		printf("The swap chain's format changed, and the pipelines were made for the old one\n");
		return false;
	}
	state->color_format = surface_format.format;
	if (capabilities.currentExtent.width != UINT32_MAX) {
		state->extent = capabilities.currentExtent;
	} else {
		int width = 0;
		int height = 0;
		glfwGetFramebufferSize(state->window, &width, &height);
		state->extent.width = (uint32_t)width;
		state->extent.height = (uint32_t)height;
	}
	uint32_t image_count = capabilities.minImageCount + 1;
	if (capabilities.maxImageCount && image_count > capabilities.maxImageCount) {
		image_count = capabilities.maxImageCount;
	}
	image_count = image_count < EXAMPLE_MAX_SWAPCHAIN_IMAGES ? image_count : EXAMPLE_MAX_SWAPCHAIN_IMAGES;
	VkSwapchainKHR old_swapchain = state->swapchain;
	VkSwapchainCreateInfoKHR create_info = { VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR };
	create_info.surface = state->surface;
	create_info.minImageCount = image_count;
	create_info.imageFormat = surface_format.format;
	create_info.imageColorSpace = surface_format.colorSpace;
	create_info.imageExtent = state->extent;
	create_info.imageArrayLayers = 1;
	create_info.imageUsage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
	create_info.imageSharingMode = VK_SHARING_MODE_EXCLUSIVE;
	create_info.preTransform = capabilities.currentTransform;
	create_info.compositeAlpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
	create_info.presentMode = VK_PRESENT_MODE_FIFO_KHR;
	create_info.clipped = VK_TRUE;
	create_info.oldSwapchain = old_swapchain;
	VkResult result = vkCreateSwapchainKHR(state->device, &create_info, NULL, &state->swapchain);
	vkDestroySwapchainKHR(state->device, old_swapchain, NULL);
	if (result != VK_SUCCESS) {
		printf("Could not create the swap chain (%i)\n", (int)result);
		state->swapchain = VK_NULL_HANDLE;
		return false;
	}
	vkGetSwapchainImagesKHR(state->device, state->swapchain, &state->image_count, NULL);
	state->image_count = state->image_count < EXAMPLE_MAX_SWAPCHAIN_IMAGES ? state->image_count : EXAMPLE_MAX_SWAPCHAIN_IMAGES;
	vkGetSwapchainImagesKHR(state->device, state->swapchain, &state->image_count, state->swapchain_images);
	for (uint32_t index = 0; index != state->image_count; ++index) {
		VkImageViewCreateInfo view_info = { VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO };
		view_info.image = state->swapchain_images[index];
		view_info.viewType = VK_IMAGE_VIEW_TYPE_2D;
		view_info.format = state->color_format;
		view_info.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
		view_info.subresourceRange.levelCount = 1;
		view_info.subresourceRange.layerCount = 1;
		VkSemaphoreCreateInfo semaphore_info = { VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO };
		if (vkCreateImageView(state->device, &view_info, NULL, &state->swapchain_views[index]) != VK_SUCCESS
			|| vkCreateSemaphore(state->device, &semaphore_info, NULL, &state->render_finished[index]) != VK_SUCCESS) {
			printf("Could not create the swap chain views\n");
			return false;
		}
	}
	if (!CreateDepthBuffer(state)) {
		printf("Could not create the depth buffer\n");
		return false;
	}
	state->swapchain_stale = false;
	return true;
}

static bool RecreateSwapchain(example_t *state) {
	vkDeviceWaitIdle(state->device);
	ReleaseSwapchainImages(state);
	return CreateSwapchain(state);
}

static void ReleaseDevice(example_t *state) {
	if (state->device) {
		ReleaseSwapchainImages(state);
		vkDestroySwapchainKHR(state->device, state->swapchain, NULL);
		for (uint32_t frame = 0; frame != EXAMPLE_FRAME_COUNT; ++frame) {
			vkDestroyCommandPool(state->device, state->command_pools[frame], NULL);
			vkDestroyFence(state->device, state->frame_fences[frame], NULL);
			vkDestroySemaphore(state->device, state->image_available[frame], NULL);
		}
		//The validation layer reports any object still alive here
		vkDestroyDevice(state->device, NULL);
	}
	if (state->instance) {
		vkDestroySurfaceKHR(state->instance, state->surface, NULL);
		if (state->messenger) {
			PFN_vkDestroyDebugUtilsMessengerEXT destroy_messenger = (PFN_vkDestroyDebugUtilsMessengerEXT)vkGetInstanceProcAddr(state->instance, "vkDestroyDebugUtilsMessengerEXT");
			if (destroy_messenger) {
				destroy_messenger(state->instance, state->messenger, NULL);
			}
		}
		vkDestroyInstance(state->instance, NULL);
	}
}

//Every atlas after the load, and the second library's after a reload
static bool UploadPendingAtlases(example_t *state, VkCommandBuffer command_buffer) {
	bool uploaded = true;
	tfx_impl_atlas_t *atlas = example_scene_GetPendingAtlas(&state->scene);
	while (atlas) {
		uploaded = tfx_impl_vulkan_UploadAtlas(&state->renderer, command_buffer, atlas) && uploaded;
		tfx_impl_FreeAtlasPixels(atlas);
		atlas = example_scene_GetPendingAtlas(&state->scene);
	}
	return uploaded;
}

static void Submit(example_t *state, VkCommandBuffer command_buffer, VkSemaphore wait_semaphore, VkSemaphore signal_semaphore, VkFence fence) {
	VkCommandBufferSubmitInfo command_buffer_info = { VK_STRUCTURE_TYPE_COMMAND_BUFFER_SUBMIT_INFO };
	command_buffer_info.commandBuffer = command_buffer;
	VkSemaphoreSubmitInfo wait_info = { VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO };
	wait_info.semaphore = wait_semaphore;
	wait_info.stageMask = VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT;
	VkSemaphoreSubmitInfo signal_info = { VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO };
	signal_info.semaphore = signal_semaphore;
	signal_info.stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
	VkSubmitInfo2 submit_info = { VK_STRUCTURE_TYPE_SUBMIT_INFO_2 };
	submit_info.waitSemaphoreInfoCount = wait_semaphore ? 1 : 0;
	submit_info.pWaitSemaphoreInfos = &wait_info;
	submit_info.commandBufferInfoCount = 1;
	submit_info.pCommandBufferInfos = &command_buffer_info;
	submit_info.signalSemaphoreInfoCount = signal_semaphore ? 1 : 0;
	submit_info.pSignalSemaphoreInfos = &signal_info;
	vkQueueSubmit2(state->queue, 1, &submit_info, fence);
}

//The atlas and sprite data go up in one command buffer, waited on before the first frame
static bool LoadScene(example_t *state) {
	if (!example_scene_Load(&state->scene, EXAMPLE_LIBRARY_PATH, EXAMPLE_SECOND_LIBRARY_PATH, EXAMPLE_SPRITE_DATA_PATH, EXAMPLE_FEATURES)
		|| !tfx_impl_vulkan_CreateAtlas(&state->renderer, state->scene.atlas_layer_size, state->scene.atlas_layer_count)) {
		return false;
	}
	VkCommandBuffer command_buffer = state->command_buffers[0];
	VkCommandBufferBeginInfo begin_info = { VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO };
	begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
	vkResetCommandPool(state->device, state->command_pools[0], 0);
	vkBeginCommandBuffer(command_buffer, &begin_info);
	tfx_impl_vulkan_BeginFrame(&state->renderer, 0);
	bool uploaded = UploadPendingAtlases(state, command_buffer)
		&& (!state->scene.animation_manager || tfx_impl_vulkan_UploadSpriteData(&state->renderer, command_buffer, state->scene.animation_manager, state->scene.baked_shapes));
	vkEndCommandBuffer(command_buffer);
	Submit(state, command_buffer, VK_NULL_HANDLE, VK_NULL_HANDLE, VK_NULL_HANDLE);
	vkQueueWaitIdle(state->queue);
	return uploaded;
}

static void UpdateWindowTitle(example_t *state, unsigned int frames_per_second) {
	tfx_version_t version = tfx_GetVersion();
	char title[256];
	tfx_gpu_particle_stats_t gpu_stats = tfx_GetGPUParticleStats(state->scene.stage);
	snprintf(title, sizeof(title), "TimelineFX Vulkan Example - %s %i.%i - %u fps - %u CPU particles - %u GPU particles - %llu dropped GPU spawns",
		version.name, version.major, version.minor, frames_per_second, tfx_GetParticleCount(state->scene.stage), gpu_stats.estimated_live_particles,
		(unsigned long long)gpu_stats.dropped_spawns);
	glfwSetWindowTitle(state->window, title);
}

static void TransitionAttachments(example_t *state, VkCommandBuffer command_buffer, uint32_t image_index, bool to_present) {
	VkImageMemoryBarrier2 barriers[2] = {};
	for (int index = 0; index != 2; ++index) {
		barriers[index].sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2;
		barriers[index].srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
		barriers[index].dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
		barriers[index].subresourceRange.levelCount = 1;
		barriers[index].subresourceRange.layerCount = 1;
	}
	VkImageMemoryBarrier2 *color = &barriers[0];
	color->image = state->swapchain_images[image_index];
	color->subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
	color->srcStageMask = VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT;
	if (to_present) {
		color->srcAccessMask = VK_ACCESS_2_COLOR_ATTACHMENT_WRITE_BIT;
		color->oldLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
		color->newLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
	} else {
		//The acquire semaphore is waited on at this stage, so the image is free once it starts
		color->dstStageMask = VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT;
		color->dstAccessMask = VK_ACCESS_2_COLOR_ATTACHMENT_WRITE_BIT;
		color->oldLayout = VK_IMAGE_LAYOUT_UNDEFINED;
		color->newLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
	}
	//One depth buffer for every frame, so the last frame's depth tests finish before this one clears it
	VkImageMemoryBarrier2 *depth = &barriers[1];
	depth->image = state->depth_image;
	depth->subresourceRange.aspectMask = VK_IMAGE_ASPECT_DEPTH_BIT;
	depth->srcStageMask = VK_PIPELINE_STAGE_2_EARLY_FRAGMENT_TESTS_BIT | VK_PIPELINE_STAGE_2_LATE_FRAGMENT_TESTS_BIT;
	depth->srcAccessMask = VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
	depth->dstStageMask = VK_PIPELINE_STAGE_2_EARLY_FRAGMENT_TESTS_BIT | VK_PIPELINE_STAGE_2_LATE_FRAGMENT_TESTS_BIT;
	depth->dstAccessMask = VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_READ_BIT | VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
	depth->oldLayout = VK_IMAGE_LAYOUT_UNDEFINED;
	depth->newLayout = VK_IMAGE_LAYOUT_DEPTH_ATTACHMENT_OPTIMAL;
	VkDependencyInfo dependency = { VK_STRUCTURE_TYPE_DEPENDENCY_INFO };
	dependency.imageMemoryBarrierCount = to_present ? 1 : 2;
	dependency.pImageMemoryBarriers = barriers;
	vkCmdPipelineBarrier2(command_buffer, &dependency);
}

//Usage: example_glfw_vulkan [updates per second] [--frames N] [--test-debug-output] [--gpu-validation]
int main(int argc, char **argv) {
	double updates_per_second = 0.0;
	long frame_limit = 0;
	bool test_debug_output = false;
	bool gpu_validation = false;
	for (int index = 1; index < argc; ++index) {
		if (strcmp(argv[index], "--gpu-validation") == 0) {
			gpu_validation = true;
		} else if (strcmp(argv[index], "--frames") == 0 && index + 1 < argc) {
			frame_limit = atol(argv[++index]);
		} else if (strcmp(argv[index], "--test-debug-output") == 0) {
			test_debug_output = true;
		} else if (atof(argv[index]) > 0.0) {
			updates_per_second = atof(argv[index]);
		}
	}

	glfwSetErrorCallback(ErrorCallback);
	if (!glfwInit()) {
		return 1;
	}
	glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
	example.window = glfwCreateWindow(EXAMPLE_WINDOW_WIDTH, EXAMPLE_WINDOW_HEIGHT, "TimelineFX Vulkan Example", NULL, NULL);
	bool ready = example.window && CreateInstance(&example, gpu_validation)
		&& glfwCreateWindowSurface(example.instance, example.window, NULL, &example.surface) == VK_SUCCESS
		&& PickPhysicalDevice(&example) && CreateDevice(&example) && CreateSwapchain(&example);
	if (example.window) {
		glfwSetKeyCallback(example.window, KeyCallback);
		glfwSetScrollCallback(example.window, ScrollCallback);
	}

	//Starts the worker threads and the memory pool that everything in the library allocates from
	tfx_BeginTimelineFX(tfx_GetDefaultThreadCount(), EXAMPLE_MEMORY_POOL_SIZE);
	ready = ready && tfx_impl_vulkan_Init(&example.renderer, example.physical_device, example.device, example.color_format, EXAMPLE_DEPTH_FORMAT, TFX_EXAMPLE_SHADER_FOLDER, EXAMPLE_FRAME_COUNT)
		&& LoadScene(&example);
	if (ready) {
		//A zero size buffer, to prove the validation layer's messages reach the console before trusting their silence
		if (test_debug_output) {
			if (example.messenger) {
				VkBufferCreateInfo invalid = { VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO };
				invalid.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT;
				VkBuffer buffer = VK_NULL_HANDLE;
				vkCreateBuffer(example.device, &invalid, NULL, &buffer);
				vkDestroyBuffer(example.device, buffer, NULL);
			} else {
				printf("The validation layer isn't on, so --test-debug-output has nothing to test\n");
			}
		}
		//Before the start sizes the GPU rings, so they hold a slower rate's longer updates
		if (updates_per_second > 0.0) {
			tfx_SetGPUSizingUpdateRate((float)updates_per_second);
		}
		example_scene_Start(&example.scene);
		if (updates_per_second > 0.0) {
			example.scene.timer = tfx_impl_CreateTimer(updates_per_second);
		}
		printf("Space: burst  O: pause the camera  L: free and reload the second library  Mouse wheel: zoom  Escape: quit\n");
	}

	double last_time = glfwGetTime();
	double title_timer = 0.0;
	unsigned int frame_count = 0;
	long frames_run = 0;
	uint32_t frame_number = 0;

	bool running = ready;
	while (running && !glfwWindowShouldClose(example.window)) {
		glfwPollEvents();
		int width = 0;
		int height = 0;
		glfwGetFramebufferSize(example.window, &width, &height);
		if (width == 0 || height == 0) {
			glfwWaitEvents();
			continue;
		}
		if (example.swapchain_stale || (uint32_t)width != example.extent.width || (uint32_t)height != example.extent.height) {
			if (!RecreateSwapchain(&example)) {
				break;
			}
		}

		//Wait until the GPU has finished the last frame that used this slot's command buffer and upload memory
		uint32_t frame_index = frame_number % EXAMPLE_FRAME_COUNT;
		vkWaitForFences(example.device, 1, &example.frame_fences[frame_index], VK_TRUE, UINT64_MAX);
		uint32_t image_index = 0;
		VkResult acquired = vkAcquireNextImageKHR(example.device, example.swapchain, UINT64_MAX, example.image_available[frame_index], VK_NULL_HANDLE, &image_index);
		if (acquired == VK_ERROR_OUT_OF_DATE_KHR) {
			example.swapchain_stale = true;
			continue;
		} else if (acquired != VK_SUCCESS && acquired != VK_SUBOPTIMAL_KHR) {
			printf("Could not acquire a swap chain image (%i)\n", (int)acquired);
			break;
		}
		vkResetFences(example.device, 1, &example.frame_fences[frame_index]);
		vkResetCommandPool(example.device, example.command_pools[frame_index], 0);
		VkCommandBuffer command_buffer = example.command_buffers[frame_index];
		VkCommandBufferBeginInfo begin_info = { VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO };
		begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
		vkBeginCommandBuffer(command_buffer, &begin_info);
		tfx_impl_vulkan_BeginFrame(&example.renderer, frame_index);

		double time = glfwGetTime();
		double elapsed_milliseconds = (time - last_time) * 1000.0;
		last_time = time;

		tfx_impl_frame_t frame;
		tfx_impl_uniforms_t uniforms;
		example_scene_Update(&example.scene, elapsed_milliseconds, (float)example.extent.width, (float)example.extent.height,
			(tfx_impl_clip_space_flags)(tfx_impl_clip_space_depth_zero_to_one | tfx_impl_clip_space_flip_y), &frame, &uniforms);
		UploadPendingAtlases(&example, command_buffer);
		tfx_impl_vulkan_UploadFrame(&example.renderer, command_buffer, &frame);
		tfx_impl_vulkan_Compute(&example.renderer, command_buffer, &uniforms);

		TransitionAttachments(&example, command_buffer, image_index, false);
		VkRenderingAttachmentInfo color_attachment = { VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO };
		color_attachment.imageView = example.swapchain_views[image_index];
		color_attachment.imageLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
		color_attachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
		color_attachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
		color_attachment.clearValue.color.float32[0] = 0.f;
		color_attachment.clearValue.color.float32[1] = 0.1f;
		color_attachment.clearValue.color.float32[2] = 0.2f;
		color_attachment.clearValue.color.float32[3] = 1.f;
		VkRenderingAttachmentInfo depth_attachment = { VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO };
		depth_attachment.imageView = example.depth_view;
		depth_attachment.imageLayout = VK_IMAGE_LAYOUT_DEPTH_ATTACHMENT_OPTIMAL;
		depth_attachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
		depth_attachment.storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
		depth_attachment.clearValue.depthStencil.depth = 1.f;
		VkRenderingInfo rendering_info = { VK_STRUCTURE_TYPE_RENDERING_INFO };
		rendering_info.renderArea.extent = example.extent;
		rendering_info.layerCount = 1;
		rendering_info.colorAttachmentCount = 1;
		rendering_info.pColorAttachments = &color_attachment;
		rendering_info.pDepthAttachment = &depth_attachment;
		vkCmdBeginRendering(command_buffer, &rendering_info);
		VkViewport viewport = { 0.f, 0.f, (float)example.extent.width, (float)example.extent.height, 0.f, 1.f };
		VkRect2D scissor = { { 0, 0 }, example.extent };
		vkCmdSetViewport(command_buffer, 0, 1, &viewport);
		vkCmdSetScissor(command_buffer, 0, 1, &scissor);
		tfx_impl_vulkan_Render(&example.renderer, command_buffer, &uniforms);
		vkCmdEndRendering(command_buffer);
		TransitionAttachments(&example, command_buffer, image_index, true);
		vkEndCommandBuffer(command_buffer);
		Submit(&example, command_buffer, example.image_available[frame_index], example.render_finished[image_index], example.frame_fences[frame_index]);

		VkPresentInfoKHR present_info = { VK_STRUCTURE_TYPE_PRESENT_INFO_KHR };
		present_info.waitSemaphoreCount = 1;
		present_info.pWaitSemaphores = &example.render_finished[image_index];
		present_info.swapchainCount = 1;
		present_info.pSwapchains = &example.swapchain;
		present_info.pImageIndices = &image_index;
		VkResult presented = vkQueuePresentKHR(example.queue, &present_info);
		if (presented == VK_ERROR_OUT_OF_DATE_KHR || presented == VK_SUBOPTIMAL_KHR) {
			example.swapchain_stale = true;
		}
		frame_number++;

		frame_count++;
		title_timer += elapsed_milliseconds;
		if (title_timer >= 1000.0) {
			UpdateWindowTitle(&example, frame_count);
			title_timer -= 1000.0;
			frame_count = 0;
		}
		if (frame_limit > 0 && ++frames_run >= frame_limit) {
			running = false;
		}
	}

	if (example.device) {
		vkDeviceWaitIdle(example.device);
	}
	example_scene_Shutdown(&example.scene);
	tfx_EndTimelineFX();
	tfx_impl_vulkan_Shutdown(&example.renderer);

	ReleaseDevice(&example);
	if (example.window) {
		glfwDestroyWindow(example.window);
	}
	glfwTerminate();
	return 0;
}
