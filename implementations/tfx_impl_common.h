#ifndef TFX_IMPL_COMMON_HEADER
#define TFX_IMPL_COMMON_HEADER

/*
	TimelineFX helpers shared by every graphics API backend. Nothing in here talks to a GPU, see README.md.

	Table of contents, search for the keyword to jump to the section:
	[Impl_math]			vec3 and column major mat4 types and functions
	[Impl_timer]		Fixed step timer for the stage update and the shader interpolation value
	[Impl_camera]		Look at camera with the view and projection matrices
	[Impl_shapes]		Shape loader callback that decodes every shape format to RGBA8 frames
	[Impl_atlas]		Packs every shape frame into the layers of one texture array, and the uv lookup callback
	[Impl_frame]		What a backend uploads and draws each frame: instances and the shader uniforms
*/

#include "timelinefx.h"

#ifdef __cplusplus
extern "C" {
#endif

//--- Impl_math

typedef struct tfx_impl_vec3_s {
	float x, y, z;
} tfx_impl_vec3_t;

//Column major, so columns[3] holds the translation, matching GLSL and HLSL's default constant buffer packing
typedef struct tfx_impl_mat4_s {
	float columns[4][4];
} tfx_impl_mat4_t;

tfx_impl_vec3_t tfx_impl_Vec3(float x, float y, float z);
tfx_impl_vec3_t tfx_impl_AddVec3(tfx_impl_vec3_t a, tfx_impl_vec3_t b);
tfx_impl_vec3_t tfx_impl_SubVec3(tfx_impl_vec3_t a, tfx_impl_vec3_t b);
tfx_impl_vec3_t tfx_impl_ScaleVec3(tfx_impl_vec3_t v, float scale);
tfx_impl_vec3_t tfx_impl_CrossVec3(tfx_impl_vec3_t a, tfx_impl_vec3_t b);
float tfx_impl_DotVec3(tfx_impl_vec3_t a, tfx_impl_vec3_t b);
tfx_impl_vec3_t tfx_impl_NormalizeVec3(tfx_impl_vec3_t v);
tfx_impl_mat4_t tfx_impl_IdentityMat4(void);
//result = a * b, so b is applied to a vector first
tfx_impl_mat4_t tfx_impl_MultiplyMat4(const tfx_impl_mat4_t *a, const tfx_impl_mat4_t *b);
//Right handed view matrix looking down -z, like gluLookAt
tfx_impl_mat4_t tfx_impl_LookAt(tfx_impl_vec3_t eye, tfx_impl_vec3_t target, tfx_impl_vec3_t up);

typedef enum tfx_impl_clip_space_flags {
	tfx_impl_clip_space_none = 0,
	tfx_impl_clip_space_depth_zero_to_one = 1 << 0,		//D3D, Vulkan, Metal and WebGPU. Without it depth is -1 to 1 like OpenGL
	tfx_impl_clip_space_flip_y = 1 << 1,				//Vulkan, whose clip space y points down
} tfx_impl_clip_space_flags;

//Right handed perspective projection, the field of view is vertical and in radians
tfx_impl_mat4_t tfx_impl_Perspective(float field_of_view, float aspect_ratio, float near_plane, float far_plane, tfxU32 clip_space_flags);

//--- Impl_timer

typedef struct tfx_impl_timer_s {
	double update_frequency;		//Updates per second
	double update_length;			//Milliseconds per update
	double accumulator;				//Milliseconds not yet simulated
	float lerp;						//How far between the previous and the latest update the frame is, 0 to 1
	tfxU32 max_ticks_per_frame;		//Time beyond this many updates is dropped so a long stall doesn't snowball
	tfxU32 ticks_this_frame;
} tfx_impl_timer_t;

tfx_impl_timer_t tfx_impl_CreateTimer(double update_frequency);
//Returns the number of updates to run this frame, also stored in ticks_this_frame
tfxU32 tfx_impl_AdvanceTimer(tfx_impl_timer_t *timer, double elapsed_milliseconds);

//--- Impl_camera

typedef struct tfx_impl_camera_s {
	tfx_impl_vec3_t position;
	tfx_impl_vec3_t target;
	tfx_impl_vec3_t up;
	float field_of_view;			//Vertical, in radians
	float near_plane;
	float far_plane;
	tfx_impl_mat4_t view;
	tfx_impl_mat4_t projection;
} tfx_impl_camera_t;

tfx_impl_camera_t tfx_impl_CreateCamera(tfx_impl_vec3_t position, tfx_impl_vec3_t target, float field_of_view);
//Places the camera on a sphere around its target, yaw around the y axis and pitch up from the horizon, both in radians
void tfx_impl_OrbitCamera(tfx_impl_camera_t *camera, float yaw, float pitch, float distance);
void tfx_impl_UpdateCameraMatrices(tfx_impl_camera_t *camera, float aspect_ratio, tfxU32 clip_space_flags);
//Pass to tfx_SetStageCamera so effects can face and depth sort against the camera
tfx_impl_vec3_t tfx_impl_CameraFront(const tfx_impl_camera_t *camera);

//--- Impl_shapes

typedef struct tfx_impl_atlas_region_s {
	float uv[4];					//left, top, right, bottom
	tfxU32 layer;
} tfx_impl_atlas_region_t;

//One per shape. The library keeps a pointer to it (tfx_SetImagePointer) so each one is allocated on its own and never moves
typedef struct tfx_impl_shape_s {
	tfxU32 frame_count;
	tfxU32 frame_width;
	tfxU32 frame_height;
	tfxU32 sheet_width;				//Frames are laid out left to right, top to bottom on the sheet
	tfxU32 sheet_height;
	unsigned char *pixels;			//RGBA8 sheet, freed once the atlas is built
	tfx_impl_atlas_region_t *regions;	//One per frame, filled by tfx_impl_BuildShapeAtlas
	bool unreadable;				//Couldn't be decoded, drawn as a white square
} tfx_impl_shape_t;

//Pass a pointer to this as the user data of tfx_LoadEffectLibrary
typedef struct tfx_impl_shapes_s {
	tfx_impl_shape_t **shapes;
	tfxU32 count;
	tfxU32 capacity;
	tfxU32 unreadable_count;
} tfx_impl_shapes_t;

//The tfx_shape_loader callback
void tfx_impl_ShapeLoader(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_memory_size, void *user_data);
void tfx_impl_FreeShapes(tfx_impl_shapes_t *shapes);

//--- Impl_atlas

typedef struct tfx_impl_atlas_s {
	tfxU32 layer_size;				//Width and height of every layer
	tfxU32 layer_count;
	unsigned char *pixels;			//RGBA8, layer after layer, freed with tfx_impl_FreeAtlasPixels once uploaded
} tfx_impl_atlas_t;

//Packs every frame of every shape and fills their regions. Frees the shape sheets. Call after loading, before tfx_UpdateLibraryGPUImageData
bool tfx_impl_BuildShapeAtlas(tfx_impl_shapes_t *shapes, tfx_impl_atlas_t *atlas);
void tfx_impl_FreeAtlasPixels(tfx_impl_atlas_t *atlas);
//The tfx_uv_lookup callback. The layer goes in texture_array_index, whose upper 16 bits (a bindless image index in the reference shaders) stay 0
void tfx_impl_GetUV(void *pointer, tfx_gpu_image_data_t *image_data, int offset);

//--- Impl_frame

//Mirrors the frame constants every backend's shaders read
typedef struct tfx_impl_uniforms_s {
	tfx_impl_mat4_t view;
	tfx_impl_mat4_t projection;
	float screen_size[2];
	tfxU32 milliseconds;
	float timer_lerp;
	float update_time;
	float camera_position[3];		//Faces baked ribbons to the camera, a stage's ribbons take it from tfx_SetStageCamera
} tfx_impl_uniforms_t;

//CPU memory a frame's ribbons are packed into before they're uploaded. It grows to fit the stages that exist
typedef struct tfx_impl_ribbon_staging_s {
	void *segments;
	void *ribbons;
	void *emitters;
	size_t segments_capacity;
	size_t ribbons_capacity;
	size_t emitters_capacity;
} tfx_impl_ribbon_staging_t;

typedef struct tfx_impl_ribbon_frame_s {
	tfx_stage stage;				//Null when there are no ribbons to draw
	tfx_ribbon_batch_offsets_t stage_offsets;	//Where the stage starts in the batch, for tfx_CreateRibbonBatchDispatch
	tfx_ribbon_batch_sizes_t sizes;	//The bytes to upload and the vertex and index buffer sizes the compute shader writes
	const void *segments;			//tfx_ribbon_segment_t, sizes.segment_buffer_size_in_bytes of them
	const void *ribbons;			//tfx_ribbon_t
	const void *emitters;			//tfx_gpu_ribbon_emitter_t
} tfx_impl_ribbon_frame_t;

//Baked sprite data playback, gathered every frame because animation instances advance every frame rather than on ticks
typedef struct tfx_impl_baked_frame_s {
	tfxU32 sprite_count;			//Billboards the playback compute writes this frame, nothing to do when it's 0
	tfxU32 animation_instance_count;
	bool animated_shapes;
	const void *offsets;			//A running sprite total per animation instance
	size_t offsets_bytes;
	const void *animation_instances;
	size_t animation_instances_bytes;
	tfx_animation_manager animation_manager;	//Its ribbon buckets are walked with tfx_NextAnimationRibbonDispatch
	tfxU32 ribbon_count;			//Ribbons the ribbon playback compute writes this frame
	const void *ribbon_offsets;		//A running ribbon total per animation instance
	size_t ribbon_offsets_bytes;
} tfx_impl_baked_frame_t;

typedef struct tfx_impl_frame_s {
	tfxU32 ticks;					//Updates run this frame, nothing needs uploading when it's 0
	tfx_instance_t *instances;		//Layer major, so the whole buffer goes up in one copy
	tfxU32 instance_count;
	tfxU32 layer_start[tfxLAYERS];	//The CPU instances of each layer
	tfxU32 layer_count[tfxLAYERS];
	tfx_impl_ribbon_frame_t ribbons;
	//Only fetched when asked for, because every fetch has to be dispatched. The GPU instances go after the CPU ones
	bool has_gpu_particles;
	tfx_gpu_particle_frame_t gpu_particles;
	tfx_gpu_particle_buffer_sizes_t gpu_pool_sizes;	//Fixed for the stage's life
	tfx_impl_baked_frame_t baked;
} tfx_impl_frame_t;

//Call after the stage update (or every frame with ticks 0) and waits for it to finish. Pass null ribbon staging if the backend
//doesn't draw ribbons, and only fetch the GPU particles if the backend dispatches them
void tfx_impl_GatherFrame(tfx_stage stage, tfxU32 ticks, tfx_impl_ribbon_staging_t *ribbon_staging, bool fetch_gpu_particles, tfx_impl_frame_t *frame);
void tfx_impl_FreeRibbonStaging(tfx_impl_ribbon_staging_t *ribbon_staging);
//Call every frame after tfx_UpdateAnimationManager
void tfx_impl_GatherBakedFrame(tfx_animation_manager animation_manager, tfx_impl_baked_frame_t *baked);
void tfx_impl_SetUniforms(tfx_impl_uniforms_t *uniforms, const tfx_impl_camera_t *camera, float screen_width, float screen_height, const tfx_impl_timer_t *timer);

#ifdef __cplusplus
}
#endif

#endif
