#ifndef TFX_LIBRARY_HEADER
#define TFX_LIBRARY_HEADER

#define TFX_VERSION_NAME "Alpha"
#define TFX_VERSION_MAJOR 0
#define TFX_VERSION_MINOR 33
#define TFX_VERSION_PATCH 0

#define TFX_STRINGIFY(x) #x
#if defined(__clang__)
#define TFX_DISABLE_COMPILER_WARNING(w) \
	_Pragma("clang diagnostic push") \
	_Pragma(TFX_STRINGIFY(clang diagnostic ignored w))
#define TFX_ENABLE_COMPILER_WARNING() \
	_Pragma("clang diagnostic pop")
#elif defined(__GNUC__)
#define TFX_DISABLE_COMPILER_WARNING(w)
#define TFX_ENABLE_COMPILER_WARNING()
#else
#define TFX_DISABLE_COMPILER_WARNING(w)
#define TFX_ENABLE_COMPILER_WARNING()
#endif

//Profiling is done through Tracy (https://github.com/wolfpld/tracy). Define tfxTRACY
//(along with Tracy's own TRACY_ENABLE) to turn the tfxPROFILE markers scattered through
//the library into Tracy zones. Without it they compile away to nothing, so a shipping
//build pays for no instrumentation at all.
//#define tfxTRACY
#define TFX_THREAD_SAFE
//#define TFX_EXTRA_DEBUGGING
#define SSE41		//Steam survey currently has this at 99.83% coverage 12 April 2025. I will probably make this the minimum requirement

//Override this for more layers, although currently the editor is fixed at 4
#ifndef tfxLAYERS
#define tfxLAYERS 4
#endif

//Enable this to process 8 particles at a time.
//#define tfxUSEAVX

//Enable fused multiply add in simd calculations
//#define tfxUSEFMA

/*
	Timeline FX C++ library

	This library is for implementing particle effects into your games and applications.

	This library is render agnostic, so you will have to provide your own means to render the particles. There are various API functions in the library that help you do this.

	Currently tested on Windows and MacOS, Intel and Mac based ARM processors.

	Table of contents
	Sections in this header file, you can search for the following keywords to jump to that section:

	[Header_Includes_and_Typedefs]		Base typedefs, version + platform/SIMD detection, and the tfxAPI / tfxINTERNAL markers
	[Enums]								All the definitions for enums and bit flags
	[forward_declarations]				Opaque handles, struct typedefs and the public value-type structs (tfx_random_t, tfx_stage_info_t, ...)
	[Callback_typedefs]					Callback function pointer typedefs
	[API functions]						The main functions for use by users of the library
		-[Random numbers]				Seeded random number generation
		-[Initialisation_functions]		Startup and shutdown timelinefx
		-[Global_variable_access]		Any functions that give you access to global variables relating to timelinefx
		-[Library_functions]			Functions for loading and accessing timelinefx libraries
		-[Particle_Manager_functions]	Create and update functions for effect managers where the main work is done to update particles every frame
		-[Animation_manager]			Animation manager functions for playing pre-recorded effect data
		-[Effect_templates]				Functions for working with effect templates which help modify effects in the library without actually changing the base effect in the library.
		-[Editing_graphs]				Functions to configure effect/emitter graphs
		-[General_helpers]				General math functions and other helpers.

	Everything internal - the Zest pocket allocator, containers, SIMD, work queues, hashing, profiling, file IO and all
	tfxINTERNAL / tfxAPI_EDITOR declarations - now lives in timelinefx_internal.h, included by timelinefx.cpp and the
	editor but never by a game.
*/

/*    Functions come in 3 flavours:
1) INTERNAL where they're only meant for internal use by the library and not for any use outside it. Note that these functions are declared as static.
2) API where they're meant for access within your games that you're developing. These functions are c compatible.
3) EDITOR where they can be accessed from outside the library but really they're mainly useful for editing the effects such as in in the TimelineFX Editor. These
   functions are c++ compatabile only and currently not available if you're including the library in a c project.

All functions in the library will be marked this way for clarity and naturally the API functions will all be properly documented.
*/

#ifdef __cplusplus
#define tfxAPI extern "C"
#else
#define tfxAPI 
#endif    
#if defined(__GNUC__) || defined(__clang__)
#define tfxINTERNAL static __attribute__((unused))
#else
#define tfxINTERNAL static
#endif

//----------------------------------------------------------
//Header_Includes_and_Typedefs
//----------------------------------------------------------
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64)
#define tfxX86
#include <immintrin.h>
#elif defined(__arm__) || defined(__aarch64__)
#include <arm_neon.h>
#include <mach/mach_time.h>
#define tfxARM
#endif

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#define tfxWINDOWS
#elif __APPLE__
#define tfxMAC
#elif __linux__
#define tfxLINUX
#endif

#include <stdint.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>

//type defs
typedef uint16_t tfxU16;
typedef unsigned short tfxHalf;
typedef uint32_t tfxU32;
typedef uint32_t tfxIndex;
typedef unsigned char tfxU8;
typedef unsigned int tfxEmitterID;
typedef int32_t tfxS32;
typedef uint64_t tfxU64;
typedef int64_t tfxS64;
typedef tfxU32 tfxEffectID;
typedef tfxU32 tfxAnimationID;
typedef tfxU64 tfxKey;
typedef tfxU32 tfxParticleID;
typedef short tfxShort;
typedef unsigned short tfxUShort;

<<<<<<< HEAD
//---------------------------------------
/*	Zest_Pocket_Allocator, a Two Level Segregated Fit memory allocator
	This is my own memory allocator from https://github.com/peterigz/zloc
	This is used in TimelineFX to manage memory allocation. A large pool is created and allocated from. New pools are created if it runs out of space
	(and you initialised TimelineFX to do so).
*/
//---------------------------------------
#include <assert.h>
#include <stdlib.h>
#include <stddef.h>

#define tfx__Min(a, b) (((a) < (b)) ? (a) : (b))
#define tfx__Max(a, b) (((a) > (b)) ? (a) : (b))

typedef int tfx_index;
typedef unsigned int tfx_sl_bitmap;
typedef unsigned int tfx_uint;
typedef int tfx_bool;
typedef void* tfx_pool;

#if !defined (TFX_ASSERT)
#define TFX_ASSERT assert
#endif

#define tfx__is_pow2(x) ((x) && !((x) & ((x) - 1)))
#define tfx__glue2(x, y) x ## y
#define tfx__glue(x, y) tfx__glue2(x, y)
#define tfx__static_assert(exp) \
	typedef char tfx__glue(static_assert, __LINE__) [(exp) ? 1 : -1]

#if (defined(_MSC_VER) && defined(_M_X64)) || defined(__x86_64__)
#define tfx__64BIT
typedef size_t tfx_size;
typedef size_t tfx_fl_bitmap;
typedef int32_t tfxLONG;
#define TFX_ONE 1ULL
#else
typedef size_t tfx_size;
typedef size_t tfx_fl_bitmap;
typedef int32_t tfxLONG;
#define TFX_ONE 1U
#endif

#ifndef MEMORY_ALIGNMENT_LOG2
#if defined(tfx__64BIT)
#define MEMORY_ALIGNMENT_LOG2 3		//64 bit
#else
#define MEMORY_ALIGNMENT_LOG2 2		//32 bit
#endif
#endif

#ifndef TFX_ERROR_NAME
#define TFX_ERROR_NAME "Allocator Error"
#endif

#ifndef TFX_ERROR_COLOR
#define TFX_ERROR_COLOR "\033[31m"
#endif

#ifndef TFX_NOTICE_COLOR
#define TFX_NOTICE_COLOR "\033[0m"
#endif

#ifndef TFX_NOTICE_NAME
#define TFX_NOTICE_NAME "TimelineFX Notice"
#endif

#ifdef TFX_OUTPUT_NOTICE_MESSAGES
#include <stdio.h>
#define TFX_PRINT_NOTICE(message_f, ...) printf(message_f"\033[0m", __VA_ARGS__)
#else
#define TFX_PRINT_NOTICE(message_f, ...)
#endif

#ifdef TFX_OUTPUT_ERROR_MESSAGES
#include <stdio.h>
#define TFX_PRINT_ERROR(message_f, ...) printf(message_f"\033[0m", __VA_ARGS__)
#else
#define TFX_PRINT_ERROR(message_f, ...)
#endif

#define tfx__KILOBYTE(Value) ((Value) * 1024LL)
#define tfx__MEGABYTE(Value) (tfx__KILOBYTE(Value) * 1024LL)
#define tfx__GIGABYTE(Value) (tfx__MEGABYTE(Value) * 1024LL)

#ifndef TFX_MAX_SIZE_INDEX
#if defined(tfx__64BIT)
#define TFX_MAX_SIZE_INDEX 32
#else
#define TFX_MAX_SIZE_INDEX 30
#endif
#endif

tfx__static_assert(TFX_MAX_SIZE_INDEX < 64);

#ifdef __cplusplus
extern "C" {
#endif

#define tfx__MAXIMUM_BLOCK_SIZE (TFX_ONE << TFX_MAX_SIZE_INDEX)

	enum tfx__constants {
		tfx__MEMORY_ALIGNMENT = 1 << MEMORY_ALIGNMENT_LOG2,
		tfx__SECOND_LEVEL_INDEX_LOG2 = 5,
		tfx__FIRST_LEVEL_INDEX_COUNT = TFX_MAX_SIZE_INDEX,
		tfx__SECOND_LEVEL_INDEX_COUNT = 1 << tfx__SECOND_LEVEL_INDEX_LOG2,
		tfx__BLOCK_POINTER_OFFSET = sizeof(void*) + sizeof(tfx_size),
		tfx__MINIMUM_BLOCK_SIZE = 16,
		tfx__BLOCK_SIZE_OVERHEAD = sizeof(tfx_size),
		tfx__POINTER_SIZE = sizeof(void*),
		tfx__SMALLEST_CATEGORY = (1 << (tfx__SECOND_LEVEL_INDEX_LOG2 + MEMORY_ALIGNMENT_LOG2))
	};

	typedef enum tfx__boundary_tag_flags {
		tfx__BLOCK_IS_FREE = 1 << 0,
		tfx__PREV_BLOCK_IS_FREE = 1 << 1,
	} tfx__boundary_tag_flags;

	typedef enum tfx__error_codes {
		tfx__OK,
		tfx__INVALID_FIRST_BLOCK,
		tfx__INVALID_BLOCK_FOUND,
		tfx__PHYSICAL_BLOCK_MISALIGNMENT,
		tfx__INVALID_SEGRATED_LIST,
		tfx__WRONG_BLOCK_SIZE_FOUND_IN_SEGRATED_LIST,
		tfx__SECOND_LEVEL_BITMAPS_NOT_INITIALISED
	} tfx__error_codes;

	/*
		Each block has a header that if used only has a pointer to the previous physical block
		and the size. If the block is free then the prev and next free blocks are also stored.
	*/
	typedef struct tfx_header {
		struct tfx_header *prev_physical_block;
		/*	Note that the size is either 4 or 8 bytes aligned so the boundary tag (2 flags denoting
			whether this or the previous block is free) can be stored in the first 2 least
			significant bits	*/
		tfx_size size;
		/*
		User allocation will start here when the block is used. When the block is free prev and next
		are pointers in a linked list of free blocks within the same class size of blocks
		*/
		struct tfx_header *prev_free_block;
		struct tfx_header *next_free_block;
	} tfx_header;

	/*
	A struct for making snapshots of a memory pool to get used/free memory stats
	*/
	typedef struct tfx_pool_stats_t {
		int used_blocks;
		int free_blocks;
		tfx_size free_size;
		tfx_size used_size;
	} tfx_pool_stats_t;

	typedef struct tfx_allocator {
		/*	This is basically a terminator block that free blocks can point to if they're at the end
			of a free list. */
		tfx_header null_block;
#if defined(TFX_THREAD_SAFE)
		/* Multithreading protection*/
		volatile tfxLONG access;
#endif
		tfx_size minimum_allocation_size;
		/*	Here we store all of the free block data. first_level_bitmap is either a 32bit int
		or 64bit depending on whether tfx__64BIT is set. Second_level_bitmaps are an array of 32bit
		ints. segregated_lists is a two level array pointing to free blocks or null_block if the list
		is empty. */
		tfx_fl_bitmap first_level_bitmap;
		tfx_sl_bitmap second_level_bitmaps[tfx__FIRST_LEVEL_INDEX_COUNT];
		tfx_header *segregated_lists[tfx__FIRST_LEVEL_INDEX_COUNT][tfx__SECOND_LEVEL_INDEX_COUNT];
	} tfx_allocator;

#if defined (_MSC_VER) && (_MSC_VER >= 1400) && (defined (_M_IX86) || defined (_M_X64))
	/* Microsoft Visual C++ support on x86/X64 architectures. */

#include <intrin.h>

	static inline int tfx__scan_reverse(tfx_size bitmap) {
		unsigned long index;
#if defined(tfx__64BIT)
		return _BitScanReverse64(&index, bitmap) ? index : -1;
#else
		return _BitScanReverse(&index, bitmap) ? index : -1;
#endif
	}

	static inline int tfx__scan_forward(tfx_size bitmap)
	{
		unsigned long index;
#if defined(tfx__64BIT)
		return _BitScanForward64(&index, bitmap) ? index : -1;
#else
		return _BitScanForward(&index, bitmap) ? index : -1;
#endif
	}

#ifdef _WIN32
#include <Windows.h>
	static inline tfxLONG tfx__compare_and_exchange(volatile tfxLONG* target, tfxLONG value, tfxLONG original) {
		return InterlockedCompareExchange((volatile LONG*)target, value, original);
	}

	static inline tfxLONG tfx__exchange(volatile tfxLONG* target, tfxLONG value) {
		return InterlockedExchange((volatile LONG*)target, value);
	}

	static inline uint32_t tfx__increment(uint32_t volatile *target) {
		return InterlockedIncrement(target);
	}
#endif

#define tfx__strlen strnlen_s
#define tfx__writebarrier _WriteBarrier();
#define tfx__readbarrier _ReadBarrier();
#define tfx__strcpy strcpy_s
#define tfx__fseek _fseeki64
#define tfx__ftell _ftelli64
#define TFX_ALIGN_AFFIX(v)			
#define TFX_PACKED_STRUCT

#elif defined(__GNUC__) && ((__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8)) && \
      (defined(__i386__) || defined(__x86_64__)) || defined(__clang__)
	/* GNU C/C++ or Clang support on x86/x64 architectures. */

	static inline int tfx__scan_reverse(tfx_size bitmap)
	{
#if defined(tfx__64BIT)
		return 64 - __builtin_clzll(bitmap) - 1;
#else
		return 32 - __builtin_clz((int)bitmap) - 1;
#endif
	}

	static inline int tfx__scan_forward(tfx_size bitmap)
	{
#if defined(tfx__64BIT)
		return __builtin_ffsll(bitmap) - 1;
#else
		return __builtin_ffs((int)bitmap) - 1;
#endif
	}

	static inline tfxLONG tfx__compare_and_exchange(volatile tfxLONG* target, tfxLONG value, tfxLONG original) {
		return __sync_val_compare_and_swap(target, original, value);
	}

	static inline tfxLONG tfx__exchange(volatile tfxLONG* target, tfxLONG value) {
		return __sync_lock_test_and_set(target, value);
	}

	static inline uint32_t tfx__increment(uint32_t volatile* target) {
		return __sync_add_and_fetch(target, 1);
	}

#define tfx__strlen strnlen
#define tfx__writebarrier __asm__ __volatile__ ("" : : : "memory");
#define tfx__readbarrier __asm__ __volatile__ ("" : : : "memory");
#define tfx__strcpy strcpy
#define tfx__fseek fseeko
#define tfx__ftell ftello
#define TFX_ALIGN_AFFIX(v)			__attribute__((aligned(v)))
#define TFX_PACKED_STRUCT			__attribute__((packed))

#endif

	/*
		Initialise an allocator. Pass a block of memory that you want to use to store the allocator data. This will not create a pool, only the
		necessary data structures to store the allocator.
	*/
	tfx_allocator *tfx_InitialiseAllocator(void *memory);

	/*
		Initialise an allocator and a pool at the same time. The data stucture to store the allocator will be stored at the beginning of the memory
		you pass to the function and the remaining memory will be used as the pool.
	*/
	tfx_allocator *tfx_InitialiseAllocatorWithPool(void *memory, tfx_size size, tfx_allocator **allocator);

	/*
		Add a new memory pool to the allocator. Pools don't have to all be the same size, adding a pool will create the biggest block it can within
		the pool and then add that to the segregated list of free blocks in the allocator. All the pools in the allocator will be naturally linked
		together in the segregated list because all blocks are linked together with a linked list either as physical neighbours or free blocks in
		the segregated list.
	*/
	tfx_pool *tfx_AddPool(tfx_allocator *allocator, tfx_pool *memory, tfx_size size);

	/*
		Get the structure size of an allocator. You can use this to take into account the overhead of the allocator when preparing a new allocator
		with memory pool.
	*/
	tfx_size tfx_AllocatorSize(void);

	/*
		If you initialised an allocator with a pool then you can use this function to get a pointer to the start of the pool. It won't get a pointer
		to any other pool in the allocator. You can just get that when you call tfx_AddPool.
	*/
	tfx_pool *tfx_GetPool(tfx_allocator *allocator);

	/*
		Allocate some memory within a tfx_allocator of the specified size. Minimum size is 16 bytes.
	*/
	void *tfx_Allocate(tfx_allocator *allocator, tfx_size size);

	/*
		Try to reallocate an existing memory block within the allocator. If possible the current block will be merged with the physical neigbouring
		block, otherwise a normal tfx_Allocate will take place and the data copied over to the new allocation.
	*/
	void *tfx_Reallocate(tfx_allocator *allocator, void *ptr, tfx_size size);

	/*
	Allocate some memory within a tfx_allocator of the specified size. Minimum size is 16 bytes.
*/
	void *tfx_AllocateAligned(tfx_allocator *allocator, tfx_size size, tfx_size alignment);

	/*
		Free an allocation from a tfx_allocator. When freeing a block of memory any adjacent free blocks are merged together to keep on top of
		fragmentation as much as possible. A check is also done to confirm that the block being freed is still valid and detect any memory corruption
		due to out of bounds writing of this or potentially other blocks.
	*/
	int tfx_Free(tfx_allocator *allocator, void *allocation);

	/*
		Remove a pool from an allocator. Note that all blocks in the pool must be free and therefore all merged together into one block (this happens
		automatically as all blocks are freed are merged together into bigger blocks.
	*/
	tfx_bool tfx_RemovePool(tfx_allocator *allocator, tfx_pool *pool);

	/*
	When using an allocator for managing remote memory, you need to set the bytes per block that a block storing infomation about the remote
	memory allocation will manage. For example you might set the value to 1MB so if you were to then allocate 4MB of remote memory then 4 blocks
	worth of space would be used to allocate that memory. This means that if it were to be freed and then split down to a smaller size they'd be
	enough blocks worth of space to do this.

	Note that the lower the number the more memory you need to track remote memory blocks but the more granular it will be. It will depend alot
	on the size of allocations you will need
*/
	void tfx_SetMinimumAllocationSize(tfx_allocator *allocator, tfx_size size);

	//--End of user functions

	//Private inline functions, user doesn't need to call these
	static inline void tfx__map(tfx_size size, tfx_index *fli, tfx_index *sli) {
		*fli = tfx__scan_reverse(size);
		if (*fli <= tfx__SECOND_LEVEL_INDEX_LOG2) {
			*fli = 0;
			*sli = (int)size / (tfx__SMALLEST_CATEGORY / tfx__SECOND_LEVEL_INDEX_COUNT);
			return;
		}
		size = size & ~(1 << *fli);
		*sli = (tfx_index)(size >> (*fli - tfx__SECOND_LEVEL_INDEX_LOG2)) % tfx__SECOND_LEVEL_INDEX_COUNT;
	}

	//Debug tool to make sure that if a first level bitmap has a bit set, then the corresponding second level index should contain a value
	static inline void tfx__verify_lists(tfx_allocator *allocator) {
		for (int fli = 0; fli != tfx__FIRST_LEVEL_INDEX_COUNT; ++fli) {
			if (allocator->first_level_bitmap & (1ULL << fli)) {
				//bit in first level is set but according to the second level bitmap array there are no blocks so the first level
				//bitmap bit should have been 0
				TFX_ASSERT(allocator->second_level_bitmaps[fli] > 0);
			}
		}
	}

	//Read only functions
	static inline tfx_bool tfx__has_free_block(const tfx_allocator *allocator, tfx_index fli, tfx_index sli) {
		return allocator->first_level_bitmap & (TFX_ONE << fli) && allocator->second_level_bitmaps[fli] & (1U << sli);
	}

	static inline tfx_bool tfx__is_used_block(const tfx_header *block) {
		return !(block->size & tfx__BLOCK_IS_FREE);
	}

	static inline tfx_bool tfx__is_free_block(const tfx_header *block) {
		return block->size & tfx__BLOCK_IS_FREE;
	}

	static inline tfx_bool tfx__prev_is_free_block(const tfx_header *block) {
		return block->size & tfx__PREV_BLOCK_IS_FREE;
	}

	static inline void* tfx__align_ptr(const void* ptr, tfx_size align) {
		ptrdiff_t aligned = (((ptrdiff_t)ptr) + (align - 1)) & ~(align - 1);
		TFX_ASSERT(0 == (align & (align - 1)) && "must align to a power of two");
		return (void*)aligned;
	}

	static inline tfx_bool tfx__is_aligned(tfx_size size, tfx_size alignment) {
		return (size % alignment) == 0;
	}

	static inline tfx_bool tfx__ptr_is_aligned(void *ptr, tfx_size alignment) {
		uintptr_t address = (uintptr_t)ptr;
		return (address % alignment) == 0;
	}

	static inline tfx_size tfx__align_size_down(tfx_size size, tfx_size alignment) {
		return size - (size % alignment);
	}

	static inline tfx_size tfx__align_size_up(tfx_size size, tfx_size alignment) {
		tfx_size remainder = size % alignment;
		if (remainder != 0) {
			size += alignment - remainder;
		}
		return size;
	}

	static inline tfx_size tfx__adjust_size(tfx_size size, tfx_size minimum_size, tfx_size alignment) {
		return tfx__Min(tfx__Max(tfx__align_size_up(size, alignment), minimum_size), tfx__MAXIMUM_BLOCK_SIZE);
	}

	static inline tfx_size tfx__block_size(const tfx_header *block) {
		return block->size & ~(tfx__BLOCK_IS_FREE | tfx__PREV_BLOCK_IS_FREE);
	}

	static inline tfx_header *tfx__block_from_allocation(const void *allocation) {
		return (tfx_header*)((char*)allocation - tfx__BLOCK_POINTER_OFFSET);
	}

	static inline tfx_header *tfx__null_block(tfx_allocator *allocator) {
		return &allocator->null_block;
	}

	static inline void* tfx__block_user_ptr(const tfx_header *block) {
		return (char*)block + tfx__BLOCK_POINTER_OFFSET;
	}

	static inline tfx_header* tfx__first_block_in_pool(const tfx_pool *pool) {
		return (tfx_header*)((char*)pool - tfx__POINTER_SIZE);
	}

	static inline tfx_header *tfx__next_physical_block(const tfx_header *block) {
		return (tfx_header*)((char*)tfx__block_user_ptr(block) + tfx__block_size(block));
	}

	static inline tfx_bool tfx__next_block_is_free(const tfx_header *block) {
		return tfx__is_free_block(tfx__next_physical_block(block));
	}

	static inline tfx_header *tfx__allocator_first_block(tfx_allocator *allocator) {
		return (tfx_header*)((char*)allocator + tfx_AllocatorSize() - tfx__POINTER_SIZE);
	}

	static inline tfx_bool tfx__is_last_block_in_pool(const tfx_header *block) {
		return tfx__block_size(block) == 0;
	}

	static inline tfx_index tfx__find_next_size_up(tfx_fl_bitmap map, tfx_uint start) {
		//Mask out all bits up to the start point of the scan
		map &= (~0ULL << (start + 1));
		return tfx__scan_forward(map);
	}

	//Write functions
#if defined(TFX_THREAD_SAFE)

#define tfx__lock_thread_access(alloc)										\
	do {																	\
	} while (0 != tfx__compare_and_exchange(&alloc->access, 1, 0));			\
	TFX_ASSERT(alloc->access != 0);

#define tfx__unlock_thread_access(alloc)  alloc->access = 0;

#else

#define tfx__lock_thread_access
#define tfx__unlock_thread_access 

#endif

	static inline void tfx__set_block_size(tfx_header *block, tfx_size size) {
		tfx_size boundary_tag = block->size & (tfx__BLOCK_IS_FREE | tfx__PREV_BLOCK_IS_FREE);
		block->size = size | boundary_tag;
	}

	static inline void tfx__set_prev_physical_block(tfx_header *block, tfx_header *prev_block) {
		block->prev_physical_block = prev_block;
	}

	static inline void tfx__zero_block(tfx_header *block) {
		block->prev_physical_block = 0;
		block->size = 0;
	}

	static inline void tfx__mark_block_as_used(tfx_header *block) {
		block->size &= ~tfx__BLOCK_IS_FREE;
		tfx_header *next_block = tfx__next_physical_block(block);
		next_block->size &= ~tfx__PREV_BLOCK_IS_FREE;
	}

	static inline void tfx__mark_block_as_free(tfx_header *block) {
		block->size |= tfx__BLOCK_IS_FREE;
		tfx_header *next_block = tfx__next_physical_block(block);
		next_block->size |= tfx__PREV_BLOCK_IS_FREE;
	}

	static inline void tfx__block_set_used(tfx_header *block) {
		block->size &= ~tfx__BLOCK_IS_FREE;
	}

	static inline void tfx__block_set_free(tfx_header *block) {
		block->size |= tfx__BLOCK_IS_FREE;
	}

	static inline void tfx__block_set_prev_used(tfx_header *block) {
		block->size &= ~tfx__PREV_BLOCK_IS_FREE;
	}

	static inline void tfx__block_set_prev_free(tfx_header *block) {
		block->size |= tfx__PREV_BLOCK_IS_FREE;
	}

	/*
		Push a block onto the segregated list of free blocks. Called when tfx_Free is called. Generally blocks are
		merged if possible before this is called
	*/
	static inline void tfx__push_block(tfx_allocator *allocator, tfx_header *block) {
		tfx_index fli;
		tfx_index sli;
		//Get the size class of the block
		tfx__map(tfx__block_size(block), &fli, &sli);
		tfx_header *current_block_in_free_list = allocator->segregated_lists[fli][sli];
		//If you hit this assert then it's likely that at somepoint in your code you're trying to free an allocation
		//that was already freed.
		TFX_ASSERT(block != current_block_in_free_list);
		//Insert the block into the list by updating the next and prev free blocks of
		//this and the current block in the free list. The current block in the free
		//list may well be the null_block in the allocator so this just means that this
		//block will be added as the first block in this class of free blocks.
		block->next_free_block = current_block_in_free_list;
		block->prev_free_block = &allocator->null_block;
		current_block_in_free_list->prev_free_block = block;

		allocator->segregated_lists[fli][sli] = block;
		//Flag the bitmaps to mark that this size class now contains a free block
		allocator->first_level_bitmap |= TFX_ONE << fli;
		allocator->second_level_bitmaps[fli] |= 1U << sli;
		if (allocator->first_level_bitmap & (TFX_ONE << fli)) {
			TFX_ASSERT(allocator->second_level_bitmaps[fli] > 0);
		}
		tfx__mark_block_as_free(block);
#ifdef TFX_EXTRA_DEBUGGING
		tfx__verify_lists(allocator);
#endif
	}

	/*
		Remove a block from the segregated list in the allocator and return it. If there is a next free block in the size class
		then move it down the list, otherwise unflag the bitmaps as necessary. This is only called when we're trying to allocate
		some memory with tfx_Allocate and we've determined that there's a suitable free block in segregated_lists.
	*/
	static inline tfx_header *tfx__pop_block(tfx_allocator *allocator, tfx_index fli, tfx_index sli) {
		tfx_header *block = allocator->segregated_lists[fli][sli];

		//If the block in the segregated list is actually the null_block then something went very wrong.
		//Somehow the segregated lists had the end block assigned but the first or second level bitmaps
		//did not have the masks assigned
		TFX_ASSERT(block != &allocator->null_block);
		if (block->next_free_block && block->next_free_block != &allocator->null_block) {
			//If there are more free blocks in this size class then shift the next one down and terminate the prev_free_block
			allocator->segregated_lists[fli][sli] = block->next_free_block;
			allocator->segregated_lists[fli][sli]->prev_free_block = tfx__null_block(allocator);
		}
		else {
			//There's no more free blocks in this size class so flag the second level bitmap for this class to 0.
			allocator->segregated_lists[fli][sli] = tfx__null_block(allocator);
			allocator->second_level_bitmaps[fli] &= ~(1U << sli);
			if (allocator->second_level_bitmaps[fli] == 0) {
				//And if the second level bitmap is 0 then the corresponding bit in the first lebel can be zero'd too.
				allocator->first_level_bitmap &= ~(TFX_ONE << fli);
			}
		}
		if (allocator->first_level_bitmap & (TFX_ONE << fli)) {
			TFX_ASSERT(allocator->second_level_bitmaps[fli] > 0);
		}
		tfx__mark_block_as_used(block);
#ifdef TFX_EXTRA_DEBUGGING
		tfx__verify_lists(allocator);
#endif
		return block;
	}

	/*
		Remove a block from the segregated list. This is only called when we're merging blocks together. The block is
		just removed from the list and marked as used and then merged with an adjacent block.
	*/
	static inline void tfx__remove_block_from_segregated_list(tfx_allocator *allocator, tfx_header *block) {
		tfx_index fli, sli;
		//Get the size class
		tfx__map(tfx__block_size(block), &fli, &sli);
		tfx_header *prev_block = block->prev_free_block;
		tfx_header *next_block = block->next_free_block;
		TFX_ASSERT(prev_block);
		TFX_ASSERT(next_block);
		next_block->prev_free_block = prev_block;
		prev_block->next_free_block = next_block;
		if (allocator->segregated_lists[fli][sli] == block) {
			allocator->segregated_lists[fli][sli] = next_block;
			if (next_block == tfx__null_block(allocator)) {
				allocator->second_level_bitmaps[fli] &= ~(1U << sli);
				if (allocator->second_level_bitmaps[fli] == 0) {
					allocator->first_level_bitmap &= ~(1ULL << fli);
				}
			}
		}
		if (allocator->first_level_bitmap & (TFX_ONE << fli)) {
			TFX_ASSERT(allocator->second_level_bitmaps[fli] > 0);
		}
		tfx__mark_block_as_used(block);
#ifdef TFX_EXTRA_DEBUGGING
		tfx__verify_lists(allocator);
#endif
	}

	/*
		This function is called when tfx_Allocate is called. Once a free block is found then it will be split
		if the size + header overhead + the minimum block size (16b) is greater then the size of the free block.
		If not then it simply returns the free block as it is without splitting.
		If split then the trimmed amount is added back to the segregated list of free blocks.
	*/
	static inline tfx_header *tfx__maybe_split_block(tfx_allocator *allocator, tfx_header *block, tfx_size size, tfx_size remote_size) {
		TFX_ASSERT(!tfx__is_last_block_in_pool(block));
		tfx_size size_plus_overhead = size + tfx__BLOCK_POINTER_OFFSET;
		if (size_plus_overhead + tfx__MINIMUM_BLOCK_SIZE >= tfx__block_size(block)) {
			return block;
		}
		tfx_header *trimmed = (tfx_header*)((char*)tfx__block_user_ptr(block) + size);
		trimmed->size = 0;
		tfx__set_block_size(trimmed, tfx__block_size(block) - size_plus_overhead);
		tfx_header *next_block = tfx__next_physical_block(block);
		tfx__set_prev_physical_block(next_block, trimmed);
		tfx__set_prev_physical_block(trimmed, block);
		tfx__set_block_size(block, size);
		tfx__push_block(allocator, trimmed);
		return block;
	}

	//For splitting blocks when allocating to a specific memory alignment
	static inline tfx_header *tfx__split_aligned_block(tfx_allocator *allocator, tfx_header *block, tfx_size size) {
		TFX_ASSERT(!tfx__is_last_block_in_pool(block));
		tfx_size size_minus_overhead = size - tfx__BLOCK_POINTER_OFFSET;
		tfx_header *trimmed = (tfx_header*)((char*)tfx__block_user_ptr(block) + size_minus_overhead);
		trimmed->size = 0;
		tfx__set_block_size(trimmed, tfx__block_size(block) - size);
		tfx_header *next_block = tfx__next_physical_block(block);
		tfx__set_prev_physical_block(next_block, trimmed);
		tfx__set_prev_physical_block(trimmed, block);
		tfx__set_block_size(block, size_minus_overhead);
		tfx__push_block(allocator, block);
		return trimmed;
	}

	/*
		This function is called when tfx_Free is called and the previous physical block is free. If that's the case
		then this function will merge the block being freed with the previous physical block then add that back into
		the segregated list of free blocks. Note that that happens in the tfx_Free function after attempting to merge
		both ways.
	*/
	static inline tfx_header *tfx__merge_with_prev_block(tfx_allocator *allocator, tfx_header *block) {
		TFX_ASSERT(!tfx__is_last_block_in_pool(block));
		tfx_header *prev_block = block->prev_physical_block;
		tfx__remove_block_from_segregated_list(allocator, prev_block);
		tfx__set_block_size(prev_block, tfx__block_size(prev_block) + tfx__block_size(block) + tfx__BLOCK_POINTER_OFFSET);
		tfx_header *next_block = tfx__next_physical_block(block);
		tfx__set_prev_physical_block(next_block, prev_block);
		tfx__zero_block(block);
		return prev_block;
	}

	/*
		This function might be called when tfx_Free is called to free a block. If the block being freed is not the last
		physical block then this function is called and if the next block is free then it will be merged.
	*/
	static inline void tfx__merge_with_next_block(tfx_allocator *allocator, tfx_header *block) {
		tfx_header *next_block = tfx__next_physical_block(block);
		TFX_ASSERT(next_block->prev_physical_block == block);	//could be potentional memory corruption. Check that you're not writing outside the boundary of the block size
		TFX_ASSERT(!tfx__is_last_block_in_pool(next_block));
		tfx__remove_block_from_segregated_list(allocator, next_block);
		tfx__set_block_size(block, tfx__block_size(next_block) + tfx__block_size(block) + tfx__BLOCK_POINTER_OFFSET);
		tfx_header *block_after_next = tfx__next_physical_block(next_block);
		tfx__set_prev_physical_block(block_after_next, block);
		tfx__zero_block(next_block);
	}

	static inline tfx_header *tfx__find_free_block(tfx_allocator *allocator, tfx_size size, tfx_size remote_size) {
		tfx_index fli;
		tfx_index sli;
		tfx__map(size, &fli, &sli);
		//Note that there may well be an appropriate size block in the class but that block may not be at the head of the list
		//In this situation we could opt to loop through the list of the size class to see if there is an appropriate size but instead
		//we stick to the paper and just move on to the next class up to keep a O1 speed at the cost of some extra fragmentation
		if (tfx__has_free_block(allocator, fli, sli) && tfx__block_size(allocator->segregated_lists[fli][sli]) >= size) {
			tfx_header *block = tfx__pop_block(allocator, fli, sli);
			tfx__unlock_thread_access(allocator);
			return block;
		}
		if (sli == tfx__SECOND_LEVEL_INDEX_COUNT - 1) {
			sli = -1;
		}
		else {
			sli = tfx__find_next_size_up(allocator->second_level_bitmaps[fli], sli);
		}
		if (sli == -1) {
			fli = tfx__find_next_size_up(allocator->first_level_bitmap, fli);
			if (fli > -1) {
				sli = tfx__scan_forward(allocator->second_level_bitmaps[fli]);
				tfx_header *block = tfx__pop_block(allocator, fli, sli);
				tfx_header *split_block = tfx__maybe_split_block(allocator, block, size, 0);
				tfx__unlock_thread_access(allocator);
				return split_block;
			}
		}
		else {
			tfx_header *block = tfx__pop_block(allocator, fli, sli);
			tfx_header *split_block = tfx__maybe_split_block(allocator, block, size, 0);
			tfx__unlock_thread_access(allocator);
			return split_block;
		}

		return 0;
	}
	//--End of internal functions

	//--End of header declarations

//Implementation
#if defined(TFX_ALLOCATOR_IMPLEMENTATION)

#include <limits.h>
#include <stddef.h>
#include <string.h>

//Definitions
	void* tfx_BlockUserExtensionPtr(const tfx_header *block) {
		return (char*)block + sizeof(tfx_header);
	}

	void* tfx_AllocationFromExtensionPtr(const void *block) {
		return (void*)((char*)block - tfx__MINIMUM_BLOCK_SIZE);
	}

	tfx_allocator *tfx_InitialiseAllocator(void *memory) {
		if (!memory) {
			TFX_PRINT_ERROR(TFX_ERROR_COLOR"%s: The memory pointer passed in to the initialiser was NULL, did it allocate properly?\n", TFX_ERROR_NAME);
			return 0;
		}

		tfx_allocator *allocator = (tfx_allocator*)memory;
		memset(allocator, 0, sizeof(tfx_allocator));
		allocator->null_block.next_free_block = &allocator->null_block;
		allocator->null_block.prev_free_block = &allocator->null_block;
		allocator->minimum_allocation_size = tfx__MINIMUM_BLOCK_SIZE;

		//Point all of the segregated list array pointers to the empty block
		for (tfx_uint i = 0; i < tfx__FIRST_LEVEL_INDEX_COUNT; i++) {
			for (tfx_uint j = 0; j < tfx__SECOND_LEVEL_INDEX_COUNT; j++) {
				allocator->segregated_lists[i][j] = &allocator->null_block;
			}
		}

		return allocator;
	}

	tfx_allocator *tfx_InitialiseAllocatorWithPool(void *memory, tfx_size size, tfx_allocator **allocator) {
		tfx_size array_offset = sizeof(tfx_allocator);
		if (size < array_offset + tfx__MEMORY_ALIGNMENT) {
			TFX_PRINT_ERROR(TFX_ERROR_COLOR"%s: Tried to initialise allocator with a memory allocation that is too small. Must be at least: %zi bytes\n", TFX_ERROR_NAME, array_offset + tfx__MEMORY_ALIGNMENT);
			return 0;
		}

		*allocator = tfx_InitialiseAllocator(memory);
		if (!allocator) {
			return 0;
		}
		tfx_AddPool(*allocator, tfx_GetPool(*allocator), size - tfx_AllocatorSize());
		return *allocator;
	}

	tfx_size tfx_AllocatorSize(void) {
		return sizeof(tfx_allocator);
	}

	void tfx_SetMinimumAllocationSize(tfx_allocator *allocator, tfx_size size) {
		TFX_ASSERT(allocator->minimum_allocation_size == tfx__MINIMUM_BLOCK_SIZE);		//You cannot change this once set
		TFX_ASSERT(tfx__is_pow2(size));													//Size must be a power of 2
		allocator->minimum_allocation_size = tfx__Max(tfx__MINIMUM_BLOCK_SIZE, size);
	}

	tfx_pool *tfx_GetPool(tfx_allocator *allocator) {
		return (tfx_pool*)((char*)allocator + tfx_AllocatorSize());
	}

	tfx_pool *tfx_AddPool(tfx_allocator *allocator, tfx_pool *memory, tfx_size size) {
		tfx__lock_thread_access(allocator);

		//Offset it back by the pointer size, we don't need the prev_physical block pointer as there is none
		//for the first block in the pool
		tfx_header *block = tfx__first_block_in_pool(memory);
		block->size = 0;
		//Leave room for an end block
		tfx__set_block_size(block, size - (tfx__BLOCK_POINTER_OFFSET)-tfx__BLOCK_SIZE_OVERHEAD);

		//Make sure it aligns
		tfx__set_block_size(block, tfx__align_size_down(tfx__block_size(block), tfx__MEMORY_ALIGNMENT));
		TFX_ASSERT(tfx__block_size(block) > tfx__MINIMUM_BLOCK_SIZE);
		tfx__block_set_free(block);
		tfx__block_set_prev_used(block);

		//Add a 0 sized block at the end of the pool to cap it off
		tfx_header *last_block = tfx__next_physical_block(block);
		last_block->size = 0;
		tfx__block_set_used(last_block);

		last_block->prev_physical_block = block;
		tfx__push_block(allocator, block);

		tfx__unlock_thread_access(allocator);
		return memory;
	}

	tfx_bool tfx_RemovePool(tfx_allocator *allocator, tfx_pool *pool) {
		tfx__lock_thread_access(allocator);
		tfx_header *block = tfx__first_block_in_pool(pool);

		if (tfx__is_free_block(block) && !tfx__next_block_is_free(block) && tfx__is_last_block_in_pool(tfx__next_physical_block(block))) {
			tfx__remove_block_from_segregated_list(allocator, block);
			tfx__unlock_thread_access(allocator);
			return 1;
		}
#if defined(TFX_THREAD_SAFE)
		tfx__unlock_thread_access(allocator);
		TFX_PRINT_ERROR(TFX_ERROR_COLOR"%s: In order to remove a pool there must be only 1 free block in the pool. Was possibly freed by another thread\n", TFX_ERROR_NAME);
#else
		TFX_PRINT_ERROR(TFX_ERROR_COLOR"%s: In order to remove a pool there must be only 1 free block in the pool.\n", TFX_ERROR_NAME);
#endif
		return 0;
	}

	void *tfx_Allocate(tfx_allocator *allocator, tfx_size size) {
		tfx_size remote_size = 0;
		tfx__lock_thread_access(allocator);
		size = tfx__adjust_size(size, tfx__MINIMUM_BLOCK_SIZE, tfx__MEMORY_ALIGNMENT);
		tfx_header *block = tfx__find_free_block(allocator, size, remote_size);

		if (block) {
			return tfx__block_user_ptr(block);
		}

		//Out of memory;
		TFX_PRINT_ERROR(TFX_ERROR_COLOR"%s: Not enough memory in pool to allocate %zu bytes\n", TFX_ERROR_NAME, size);
		tfx__unlock_thread_access(allocator);
		return 0;
	}

	void *tfx_Reallocate(tfx_allocator *allocator, void *ptr, tfx_size size) {
		tfx__lock_thread_access(allocator);

		if (ptr && size == 0) {
			tfx__unlock_thread_access(allocator);
			tfx_Free(allocator, ptr);
		}

		if (!ptr) {
			tfx__unlock_thread_access(allocator);
			return tfx_Allocate(allocator, size);
		}

		tfx_header *block = tfx__block_from_allocation(ptr);
		tfx_header *next_block = tfx__next_physical_block(block);
		void *allocation = 0;
		tfx_size current_size = tfx__block_size(block);
		tfx_size adjusted_size = tfx__adjust_size(size, allocator->minimum_allocation_size, tfx__MEMORY_ALIGNMENT);
		tfx_size combined_size = current_size + tfx__block_size(next_block);
		if ((!tfx__next_block_is_free(block) || adjusted_size > combined_size) && adjusted_size > current_size) {
			tfx_header *block = tfx__find_free_block(allocator, adjusted_size, 0);
			if (block) {
				allocation = tfx__block_user_ptr(block);
			}

			if (allocation) {
				tfx_size smallest_size = tfx__Min(current_size, size);
				memcpy(allocation, ptr, smallest_size);
				tfx_Free(allocator, ptr);
			}
		}
		else {
			//Reallocation is possible
			if (adjusted_size > current_size)
			{
				tfx__merge_with_next_block(allocator, block);
				tfx__mark_block_as_used(block);
			}
			tfx_header *split_block = tfx__maybe_split_block(allocator, block, adjusted_size, 0);
			allocation = tfx__block_user_ptr(split_block);
		}

		tfx__unlock_thread_access(allocator);
		return allocation;
	}

	void *tfx_AllocateAligned(tfx_allocator *allocator, tfx_size size, tfx_size alignment) {
		tfx__lock_thread_access(allocator);
		tfx_size adjusted_size = tfx__adjust_size(size, allocator->minimum_allocation_size, alignment);
		tfx_size gap_minimum = sizeof(tfx_header);
		tfx_size size_with_gap = tfx__adjust_size(adjusted_size + alignment + gap_minimum, allocator->minimum_allocation_size, alignment);
		size_t aligned_size = (adjusted_size && alignment > tfx__MEMORY_ALIGNMENT) ? size_with_gap : adjusted_size;

		tfx_header *block = tfx__find_free_block(allocator, aligned_size, 0);

		if (block) {
			void *user_ptr = tfx__block_user_ptr(block);
			void *aligned_ptr = tfx__align_ptr(user_ptr, alignment);
			tfx_size gap = (tfx_size)(((ptrdiff_t)aligned_ptr) - (ptrdiff_t)user_ptr);

			/* If gap size is too small, offset to next aligned boundary. */
			if (gap && gap < gap_minimum)
			{
				tfx_size gap_remain = gap_minimum - gap;
				tfx_size offset = tfx__Max(gap_remain, alignment);
				const void* next_aligned = (void*)((ptrdiff_t)aligned_ptr + offset);

				aligned_ptr = tfx__align_ptr(next_aligned, alignment);
				gap = (tfx_size)((ptrdiff_t)aligned_ptr - (ptrdiff_t)user_ptr);
			}

			if (gap)
			{
				TFX_ASSERT(gap >= gap_minimum && "gap size too small");
				block = tfx__split_aligned_block(allocator, block, gap);
				tfx__block_set_used(block);
			}
			TFX_ASSERT(tfx__ptr_is_aligned(tfx__block_user_ptr(block), alignment));	//pointer not aligned to requested alignment
		}
		else {
			tfx__unlock_thread_access(allocator);
			return 0;
		}

		tfx__unlock_thread_access(allocator);
		return tfx__block_user_ptr(block);
	}

	int tfx_Free(tfx_allocator *allocator, void* allocation) {
		if (!allocation) return 0;
		tfx__lock_thread_access(allocator);
		tfx_header *block = tfx__block_from_allocation(allocation);
		if (tfx__prev_is_free_block(block)) {
			TFX_ASSERT(block->prev_physical_block);		//Must be a valid previous physical block
			block = tfx__merge_with_prev_block(allocator, block);
		}
		if (tfx__next_block_is_free(block)) {
			tfx__merge_with_next_block(allocator, block);
		}
		tfx__push_block(allocator, block);
		tfx__unlock_thread_access(allocator);
		return 1;
	}

#endif

#ifdef __cplusplus
}
#endif

size_t tfxGetNextPower(size_t n);
void tfxAddHostMemoryPool(size_t size);
void* tfxAllocate(size_t size);
void* tfxReallocate(void *memory, size_t size);
void *tfxAllocateAligned(size_t size, size_t alignment);
tfx_allocator *tfxGetAllocator();

//---------------------------------------
//End of allocator code
//---------------------------------------

//----------------------------------------------------------
//Header_Includes_and_Typedefs
//----------------------------------------------------------
#if defined(_WIN32)
#include <SDKDDKVer.h>
#ifndef WIN_LEAN_AND_MEAN
#define WIN_LEAN_AND_MEAN
#endif
#include <Windows.h>
#endif

//Might possibly replace some of these in the future
#include <stdio.h>
#include <stdarg.h>					//va_list
#include <chrono>					//std::chrono::high_resolution_clock
#include <cctype>					//std::is_digit
#include <algorithm>
#include <iostream>					//temp for std::cout
#include <immintrin.h>
#include <mutex>
#include <thread>					//only using this for std::thread::hardware_ concurrency()
#include <cfloat>

#define tfxTWO63 0x8000000000000000u 
#define tfxTWO64f (tfxTWO63*2.0)
#define tfxPI 3.14159265359f
#define tfx360Radians 6.28319f
#define tfx180Radians 3.14159f
#define tfx90Radians 1.5708f
#define tfxMAXDEPTH 3

namespace tfx {
//----------------------------------------------------------
//Forward declarations

struct tfx_effect_emitter_t;
struct tfx_particle_manager_t;
struct tfx_effect_template_t;
struct tfx_compute_sprite_t;
struct tfx_compute_particle_t;
struct tfx_sprite_sheet_settings_t;
struct tfx_sprite_data_settings_t;
struct tfx_library_t;
struct tfx_str_t;
struct tfx_str16_t;
struct tfx_str32_t;
struct tfx_str64_t;
struct tfx_str128_t;
struct tfx_str256_t;
struct tfx_str512_t;

//--------------------------------------------------------------
//macros
#define TFX_VERSION "Alpha"
#define TFX_VERSION_NUMBER 3.29.2022

#define tfxMAX_FRAME 20000.f
#define tfxNullParent 0xFFFFFFFF
#define tfxINVALID 0xFFFFFFFF
#define tfxINVALID_SPRITE 0x0FFFFFFF
#define tfxEmitterPropertiesCount 26

#define tfxDel << "=" <<
#define tfxCom << "," <<
#define tfxEndLine << std::endl

#define tfxDelt "=" 
#define tfxComt ","
#define tfxEndLinet "\n"

#define tfxMin(a, b) (((a) < (b)) ? (a) : (b))
#define tfxMax(a, b) (((a) > (b)) ? (a) : (b))
#define tfxArrayCount(Array) (sizeof(Array) / sizeof((Array)[0]))

#ifndef tfxREALLOCATE
#define tfxALLOCATE(size) tfxAllocate(size)
#define tfxALLOCATE_ALIGNED(size, alignment) tfxAllocateAligned(size, alignment)
#define tfxREALLOCATE(ptr, size) tfxReallocate(ptr, size)
#endif

#ifndef tfxREALLOCATE
#define tfxALLOCATE(size) malloc(size)
#define tfxALLOCATE_ALIGNED(size, alignment) malloc(alignment)
#define tfxREALLOCATE(ptr, size) realloc(ptr, size)
#endif

#ifndef tfxFREE
#define tfxFREE(memory) tfx_Free(tfxGetAllocator(), memory)
#endif

#ifndef tfxFREE
#define tfxFREE(memory) free(memory)
#endif

#define tfxINIT_VEC_NAME 
#define tfxINIT_VEC_NAME_INIT 
#define tfxINIT_VEC_NAME_SRC_COPY 
#define tfxCONSTRUCTOR_VEC_DEF 
#define tfxCONSTRUCTOR_VEC_INIT(name) 
#define tfxCONSTRUCTOR_VEC_INIT2(name) 

typedef std::chrono::high_resolution_clock tfxClock;

/*	Functions come in 3 flavours:
1) INTERNAL where they're only meant for internal use by the library and not for any use outside it. Note that these functions are declared as static.
2) API where they're meant for access within your games that you're developing
3) EDITOR where they can be accessed from outside the library but really they're mainly useful for editing the effects such as in in the TimelineFX Editor.

All functions in the library will be marked this way for clarity and naturally the API functions will all be properly documented.
*/
//Function marker for any functions meant for external/api use
#define tfxAPI		
//Function marker for any functions meant mainly for use by the TimelineFX editor and are related to editing effects
#define tfxAPI_EDITOR		
//For internal functions
#define tfxINTERNAL static	

//Override this for more layers, although currently the editor is fixed at 4
#ifndef tfxLAYERS
#define tfxLAYERS 4
#endif 

/*
Helper macro to place inside a for loop, for example:
for(tfxEachLayer)
You can then use layer inside the loop to get the current layer
*/
#define tfxEachLayer int layer = 0; layer != tfxLAYERS; ++layer

//Internal use macro

union tfxUInt10bit
{
	struct
	{
		int x : 10;
		int y : 10;
		int z : 10;
		int w : 2;
	} data;
	tfxU32 pack;
};

//Section: OS_Specific_Functions
#ifdef _WIN32
FILE *tfx__open_file(const char *file_name, const char *mode);

tfxINTERNAL inline tfxU64 tfx_AtomicAdd64(tfxU64 volatile *value, tfxU64 amount_to_add) {
	tfxU64 result = _InterlockedExchangeAdd64((__int64 volatile*)value, amount_to_add);
	return result;
}

tfxINTERNAL inline tfxU32 tfx_AtomicAdd32(tfxU32 volatile *value, tfxU32 amount_to_add) {
	tfxU32 result = _InterlockedExchangeAdd((LONG*)value, amount_to_add);
	return result;
}
#else
FILE *tfx__open_file(const char *file_name, const char *mode);

inline tfxU64 tfx_AtomicAdd64(tfxU64 volatile* value, tfxU64 amount_to_add) {
	return __sync_fetch_and_add(value, amount_to_add);
}

inline tfxU32 tfx_AtomicAdd32(tfxU32 volatile* value, tfxU32 amount_to_add) {
	return __sync_fetch_and_add(value, amount_to_add);
}
#endif

tfxINTERNAL inline tfxU32 tfx_Millisecs() {
	auto now = tfxClock::now().time_since_epoch();
	auto m = std::chrono::duration_cast<std::chrono::milliseconds>(now).count();
	return (tfxU32)m;
}

tfxINTERNAL inline uint64_t tfx_Microsecs() {
	auto now = tfxClock::now().time_since_epoch();
	auto m = std::chrono::duration_cast<std::chrono::microseconds>(now).count();
	return m;
}

//-----------------------------------------------------------
//Section: XXHash_Implementation
//-----------------------------------------------------------

/*
	Start of xxHash code that encompasses the following license
	MIT License

	Copyright (c) 2018 Stephan Brumme

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"),
	to deal in the Software without restriction, including without limitation
	the rights to use, copy, modify, merge, publish, distribute, sublicense,
	and/or sell copies of the Software, and to permit persons to whom the Software
	is furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included
	in all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
	OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

/// XXHash (64 bit), based on Yann Collet's descriptions, see http://cyan4973.github.io/xxHash/
/** How to use:
	uint64_t myseed = 0;
	XXHash64 myhash(myseed);
	myhash.add(pointerToSomeBytes,     numberOfBytes);
	myhash.add(pointerToSomeMoreBytes, numberOfMoreBytes); // call add() as often as you like to ...
	// and compute hash:
	uint64_t result = myhash.hash();
	// or all of the above in one single line:
	uint64_t result2 = XXHash64::hash(mypointer, numBytes, myseed);
	Note: my code is NOT endian-aware !
**/
class tfxXXHash64
{
public:
	// create new XXHash (64 bit)
	/** @param seed your seed value, even zero is a valid seed **/
	explicit tfxXXHash64(uint64_t seed)
	{
		state[0] = seed + Prime1 + Prime2;
		state[1] = seed + Prime2;
		state[2] = seed;
		state[3] = seed - Prime1;
		bufferSize = 0;
		totalLength = 0;
		memset(buffer, 0, MaxBufferSize);
	}

	/// add a chunk of bytes
	/** @param  input  pointer to a continuous block of data
		@param  length number of bytes
		@return false if parameters are invalid / zero **/
	bool add(const void* input, uint64_t length)
	{
		// no data ?
		if (!input || length == 0)
			return false;

		totalLength += length;
		// byte-wise access
		const unsigned char* data = (const unsigned char*)input;

		// unprocessed old data plus new data still fit in temporary buffer ?
		if (bufferSize + length < MaxBufferSize)
		{
			// just add new data
			while (length-- > 0)
				buffer[bufferSize++] = *data++;
			return true;
		}

		// point beyond last byte
		const unsigned char* stop = data + length;
		const unsigned char* stopBlock = stop - MaxBufferSize;

		// some data left from previous update ?
		if (bufferSize > 0)
		{
			// make sure temporary buffer is full (16 bytes)
			while (bufferSize < MaxBufferSize)
				buffer[bufferSize++] = *data++;

			// process these 32 bytes (4x8)
			process(buffer, state[0], state[1], state[2], state[3]);
		}

		// copying state to local variables helps optimizer A LOT
		uint64_t s0 = state[0], s1 = state[1], s2 = state[2], s3 = state[3];
		// 32 bytes at once
		while (data <= stopBlock)
		{
			// local variables s0..s3 instead of state[0]..state[3] are much faster
			process(data, s0, s1, s2, s3);
			data += 32;
		}
		// copy back
		state[0] = s0; state[1] = s1; state[2] = s2; state[3] = s3;

		// copy remainder to temporary buffer
		bufferSize = stop - data;
		for (uint64_t i = 0; i < bufferSize; i++)
			buffer[i] = data[i];

		// done
		return true;
	}

	/// get current hash
	/** @return 64 bit XXHash **/
	uint64_t hash() const
	{
		// fold 256 bit state into one single 64 bit value
		uint64_t result;
		if (totalLength >= MaxBufferSize)
		{
			result = rotateLeft(state[0], 1) +
				rotateLeft(state[1], 7) +
				rotateLeft(state[2], 12) +
				rotateLeft(state[3], 18);
			result = (result ^ processSingle(0, state[0])) * Prime1 + Prime4;
			result = (result ^ processSingle(0, state[1])) * Prime1 + Prime4;
			result = (result ^ processSingle(0, state[2])) * Prime1 + Prime4;
			result = (result ^ processSingle(0, state[3])) * Prime1 + Prime4;
		}
		else
		{
			// internal state wasn't set in add(), therefore original seed is still stored in state2
			result = state[2] + Prime5;
		}

		result += totalLength;

		// process remaining bytes in temporary buffer
		const unsigned char* data = buffer;
		// point beyond last byte
		const unsigned char* stop = data + bufferSize;

		// at least 8 bytes left ? => eat 8 bytes per step
		for (; data + 8 <= stop; data += 8)
			result = rotateLeft(result ^ processSingle(0, *(uint64_t*)data), 27) * Prime1 + Prime4;

		// 4 bytes left ? => eat those
		if (data + 4 <= stop)
		{
			result = rotateLeft(result ^ (*(tfxU32*)data) * Prime1, 23) * Prime2 + Prime3;
			data += 4;
		}

		// take care of remaining 0..3 bytes, eat 1 byte per step
		while (data != stop)
			result = rotateLeft(result ^ (*data++) * Prime5, 11) * Prime1;

		// mix bits
		result ^= result >> 33;
		result *= Prime2;
		result ^= result >> 29;
		result *= Prime3;
		result ^= result >> 32;
		return result;
	}


	/// combine constructor, add() and hash() in one static function (C style)
	/** @param  input  pointer to a continuous block of data
		@param  length number of bytes
		@param  seed your seed value, e.g. zero is a valid seed
		@return 64 bit XXHash **/
	static uint64_t hash(const void* input, uint64_t length, uint64_t seed)
	{
		tfxXXHash64 hasher(seed);
		hasher.add(input, length);
		return hasher.hash();
	}

private:
	/// magic constants :-)
	static const uint64_t Prime1 = 11400714785074694791ULL;
	static const uint64_t Prime2 = 14029467366897019727ULL;
	static const uint64_t Prime3 = 1609587929392839161ULL;
	static const uint64_t Prime4 = 9650029242287828579ULL;
	static const uint64_t Prime5 = 2870177450012600261ULL;

	/// temporarily store up to 31 bytes between multiple add() calls
	static const uint64_t MaxBufferSize = 31 + 1;

	uint64_t      state[4];
	unsigned char buffer[MaxBufferSize];
	uint64_t      bufferSize;
	uint64_t      totalLength;

	/// rotate bits, should compile to a single CPU instruction (ROL)
	static inline uint64_t rotateLeft(uint64_t x, unsigned char bits)
	{
		return (x << bits) | (x >> (64 - bits));
	}

	/// process a single 64 bit value
	static inline uint64_t processSingle(uint64_t previous, uint64_t input)
	{
		return rotateLeft(previous + input * Prime2, 31) * Prime1;
	}

	/// process a block of 4x4 bytes, this is the main part of the XXHash32 algorithm
	static inline void process(const void* data, uint64_t& state0, uint64_t& state1, uint64_t& state2, uint64_t& state3)
	{
		const uint64_t* block = (const uint64_t*)data;
		state0 = processSingle(state0, block[0]);
		state1 = processSingle(state1, block[1]);
		state2 = processSingle(state2, block[2]);
		state3 = processSingle(state3, block[3]);
	}
};
//End of xxHash code

//----------------------------------------------------------
//Section: SIMD_defines
//----------------------------------------------------------

//Define tfxUSEAVX if you want to compile and use AVX simd operations for updating particles, otherwise SSE will be
//used by default
//Note that avx is currently only slightly faster than SSE, probably because memory bandwidth/caching becomes more of an issue at that point. But also I could be doing it wrong!
#ifdef tfxUSEAVX
#define tfxDataWidth 8	
typedef __m256 tfxWideFloat;
typedef __m256i tfxWideInt;
#define tfxWideLoad _mm256_load_ps
#define tfxWideLoadi _mm256_load_si256
#define tfxWideSet _mm256_set_ps
#define tfxWideSetSingle _mm256_set1_ps
#define tfxWideSeti _mm256_set_epi32
#define tfxWideSetSinglei _mm256_set1_epi32
#define tfxWideAdd _mm256_add_ps
#define tfxWideSub _mm256_sub_ps
#define tfxWideMul _mm256_mul_ps
#define tfxWideDiv _mm256_div_ps
#define tfxWideAddi _mm256_add_epi32
#define tfxWideSubi _mm256_sub_epi32
#define tfxWideMuli _mm256_mul_epi32
#define tfxWideSqrt _mm256_sqrt_ps
#define tfxWideMoveMask _mm256_movemask_epi8
#define tfxWideShiftRight _mm256_srli_epi32
#define tfxWideShiftLeft _mm256_slli_epi32
#define tfxWideGreaterEqual(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_GE_OS)
#define tfxWideGreater(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_GT_OS)
#define tfxWideGreateri _mm256_cmpgt_epi32
#define tfxWideLess(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_LT_OS)
#define tfxWideLessEqeual(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_LE_OS)
#define tfxWideEquals(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_EQ_OS)
#define tfxWideEqualsi _mm256_cmpeq_epi32 
#define tfxWideStore _mm256_store_ps
#define tfxWideStorei _mm256_store_si256
#define tfxWideCasti _mm256_castps_si256
#define tfxWideCast _mm256_castsi256_ps 
#define tfxWideConverti _mm256_cvttps_epi32 
#define tfxWideConvert	_mm256_cvtepi32_ps 
#define tfxWideMin _mm256_min_ps
#define tfxWideMax _mm256_max_ps
#define tfxWideMini _mm256_min_epi32
#define tfxWideMaxi _mm256_max_epi32
#define tfxWideOr _mm256_or_ps
#define tfxWideOri _mm256_or_si256
#define tfxWideXOri _mm256_xor_si256
#define tfxWideXOr _mm256_xor_ps
#define tfxWideAnd _mm256_and_ps
#define tfxWideAndi _mm256_and_si256
#define tfxWideAndNot _mm256_andnot_ps
#define tfxWideAndNoti _mm256_andnot_si256
#define tfxWideSetZero _mm256_setzero_ps
#define tfxWideSetZeroi _mm256_setzero_si256
#define tfxWideEqualsi _mm256_cmpeq_epi32 
#define tfxWideAndNot _mm256_andnot_ps
#define tfxWideLookupSet(lookup, index) tfxWideSet(lookup[index.a[7]], lookup[index.a[6]], lookup[index.a[5]], lookup[index.a[4]], lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )
#define tfxWideLookupSeti(lookup, index) tfxWideSeti(lookup[index.a[7]], lookup[index.a[6]], lookup[index.a[5]], lookup[index.a[4]], lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )
#define tfxWideLookupSetMember(lookup, member, index) tfxWideSet(lookup[index.a[7]].member, lookup[index.a[6]].member, lookup[index.a[5]].member, lookup[index.a[4]].member, lookup[index.a[3]].member, lookup[index.a[2]].member, lookup[index.a[1]].member, lookup[index.a[0]].member )
#define tfxWideLookupSetMemberi(lookup, member, index) tfxWideSeti(lookup[index.a[7]].member, lookup[index.a[6]].member, lookup[index.a[5]].member, lookup[index.a[4]].member, lookup[index.a[3]].member, lookup[index.a[2]].member, lookup[index.a[1]].member, lookup[index.a[0]].member )
#define tfxWideLookupSet2(lookup1, lookup2, index1, index2) tfxWideSet(lookup1[index1.a[7]].lookup2[index2.a[7]], lookup1[index1.a[6]].lookup2[index2.a[6]], lookup1[index1.a[5]].lookup2[index2.a[5]], lookup1[index1.a[4]].lookup2[index2.a[4]], lookup1[index1.a[3]].lookup2[index2.a[3]], lookup1[index1.a[2]].lookup2[index2.a[2]], lookup1[index1.a[1]].lookup2[index2.a[1]], lookup1[index1.a[0]].lookup2[index2.a[0]] )

const __m256 tfxWIDEF3_4 = _mm256_set1_ps(1.0f / 3.0f);
const __m256 tfxWIDEG3_4 = _mm256_set1_ps(1.0f / 6.0f);
const __m256 tfxWIDEG32_4 = _mm256_set1_ps((1.0f / 6.0f) * 2.f);
const __m256 tfxWIDEG33_4 = _mm256_set1_ps((1.0f / 6.0f) * 3.f);
const __m256i tfxWIDEONEi = _mm256_set1_epi32(1);
const __m256 tfxWIDEONE = _mm256_set1_ps(1.f);
const __m256 tfxWIDE255 = _mm256_set1_ps(255.f);
const __m256 tfxWIDEZERO = _mm256_set1_ps(0.f);
const __m256 tfxWIDETHIRTYTWO = _mm256_set1_ps(32.f);
const __m256i tfxWIDEFF = _mm256_set1_epi32(0xFF);
const __m256 tfxPWIDESIX = _mm256_set1_ps(0.6f);

typedef union {
	__m256i m;
	int a[8];
} tfxWideArrayi;

typedef union {
	__m256 m;
	float a[8];
} tfxWideArray;

#else
#define tfxDataWidth 4	
typedef __m128 tfxWideFloat;
typedef __m128i tfxWideInt;
#define tfxWideLoad _mm_load_ps
#define tfxWideLoadi _mm_load_si128
#define tfxWideSet _mm_set_ps
#define tfxWideSetSingle _mm_set_ps1
#define tfxWideSeti _mm_set_epi32
#define tfxWideSetSinglei _mm_set1_epi32
#define tfxWideAdd _mm_add_ps
#define tfxWideSub _mm_sub_ps
#define tfxWideMul _mm_mul_ps
#define tfxWideDiv _mm_div_ps
#define tfxWideAddi _mm_add_epi32
#define tfxWideSubi _mm_sub_epi32
#define tfxWideMuli _mm_mul_epu32
#define tfxWideSqrt _mm_sqrt_ps
#define tfxWideMoveMask _mm_movemask_epi8
#define tfxWideShiftRight _mm_srli_epi32
#define tfxWideShiftLeft _mm_slli_epi32
#define tfxWideGreaterEqual(v1, v2) _mm_cmpge_ps(v1, v2)
#define tfxWideGreater(v1, v2) _mm_cmpgt_ps(v1, v2)
#define tfxWideGreateri(v1, v2) _mm_cmpgt_epi32(v1, v2)
#define tfxWideLessEqual(v1, v2) _mm_cmple_ps(v1, v2)
#define tfxWideLess(v1, v2) _mm_cmplt_ps(v1, v2)
#define tfxWideLessi(v1, v2) _mm_cmplt_epi32(v1, v2)
#define tfxWideStore _mm_store_ps
#define tfxWideStorei _mm_store_si128
#define tfxWideCasti _mm_castps_si128 
#define tfxWideCast _mm_castsi128_ps
#define tfxWideConverti _mm_cvttps_epi32 
#define tfxWideConvert _mm_cvtepi32_ps 
#define tfxWideMin _mm_min_ps
#define tfxWideMax _mm_max_ps
#define tfxWideMini _mm_min_epi32
#define tfxWideMaxi _mm_max_epi32
#define tfxWideOr _mm_or_ps
#define tfxWideOri _mm_or_si128
#define tfxWideXOr _mm_xor_ps
#define tfxWideXOri _mm_xor_si128
#define tfxWideAnd _mm_and_ps
#define tfxWideAndi _mm_and_si128
#define tfxWideAndNot _mm_andnot_ps
#define tfxWideAndNoti _mm_andnot_si128
#define tfxWideSetZeroi _mm_setzero_si128
#define tfxWideSetZero _mm_setzero_ps
#define tfxWideEqualsi _mm_cmpeq_epi32 
#define tfxWideEquals _mm_cmpeq_ps
#define tfxWideShufflei _mm_shuffle_epi32

#define tfxWideLookupSet(lookup, index) tfxWideSet( lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )
#define tfxWideLookupSetMember(lookup, member, index) tfxWideSet( lookup[index.a[3]].member, lookup[index.a[2]].member, lookup[index.a[1]].member, lookup[index.a[0]].member )
#define tfxWideLookupSetMemberi(lookup, member, index) tfxWideSeti( lookup[index.a[3]].member, lookup[index.a[2]].member, lookup[index.a[1]].member, lookup[index.a[0]].member )
#define tfxWideLookupSet2(lookup1, lookup2, index1, index2) tfxWideSet( lookup1[index1.a[3]].lookup2[index2.a[3]], lookup1[index1.a[2]].lookup2[index2.a[2]], lookup1[index1.a[1]].lookup2[index2.a[1]], lookup1[index1.a[0]].lookup2[index2.a[0]] )
#define tfxWideLookupSeti(lookup, index) tfxWideSeti( lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )

const __m128 tfxWIDEF3_4 = _mm_set_ps1(1.0f / 3.0f);
const __m128 tfxWIDEG3_4 = _mm_set_ps1(1.0f / 6.0f);
const __m128 tfxWIDEG32_4 = _mm_set_ps1((1.0f / 6.0f) * 2.f);
const __m128 tfxWIDEG33_4 = _mm_set_ps1((1.0f / 6.0f) * 3.f);
const __m128i tfxWIDEONEi = _mm_set1_epi32(1);
const __m128 tfxWIDEONE = _mm_set1_ps(1.f);
const __m128 tfxWIDE255 = _mm_set1_ps(255.f);
const __m128 tfxWIDEZERO = _mm_set1_ps(0.f);
const __m128 tfxWIDETHIRTYTWO = _mm_set1_ps(32.f);
const __m128i tfxWIDEFF = _mm_set1_epi32(0xFF);
const __m128 tfxPWIDESIX = _mm_set_ps1(0.6f);

typedef union {
	__m128i m;
	int a[4];
} tfxWideArrayi;

typedef union {
	__m128 m;
	float a[4];
} tfxWideArray;

#endif

typedef __m128 tfx128;
typedef __m128i tfx128i;

typedef union {
	__m128i m;
	int a[4];
} tfx128iArray;

typedef union {
	__m128i m;
	tfxU64 a[2];
} tfx128iArray64;

typedef union {
	__m128 m;
	float a[4];
} tfx128Array;

//simd floor function thanks to Stephanie Rancourt: http://dss.stephanierct.com/DevBlog/?p=8
tfxINTERNAL inline tfx128 tfxFloor128(const tfx128& x) {
	//__m128i v0 = _mm_setzero_si128();
	//__m128i v1 = _mm_cmpeq_epi32(v0, v0);
	//__m128i ji = _mm_srli_epi32(v1, 25);
	//__m128 j = *(__m128*)&_mm_slli_epi32(ji, 23); //create vector 1.0f
	//I'm not entirely sure why original code had above lines to create a vector of 1.f. It seems to me that the below works fine 
	//Worth noting that we only need to floor small numbers for the noise algorithm so can get away with this function.
	__m128 j = _mm_set1_ps(1.f); //create vector 1.0f
	__m128i i = _mm_cvttps_epi32(x);
	__m128 fi = _mm_cvtepi32_ps(i);
	__m128 igx = _mm_cmpgt_ps(fi, x);
	j = _mm_and_ps(igx, j);
	return _mm_sub_ps(fi, j);
}

//simd mod function thanks to Stephanie Rancourt: http://dss.stephanierct.com/DevBlog/?p=8
tfxINTERNAL inline tfxWideFloat tfxWideMod(const tfxWideFloat &a, const tfxWideFloat &aDiv) {
	tfxWideFloat c = tfxWideDiv(a, aDiv);
	tfxWideInt i = tfxWideConverti(c);
	tfxWideFloat cTrunc = tfxWideConvert(i);
	tfxWideFloat base = tfxWideMul(cTrunc, aDiv);
	tfxWideFloat r = tfxWideSub(a, base);
	return r;
}

tfxINTERNAL inline tfxWideFloat tfxWideAbs(tfxWideFloat v) {
	return tfxWideAnd(tfxWideCast(tfxWideShiftRight(tfxWideSetSinglei(-1), 1)), v);
}

tfxINTERNAL inline tfxWideInt tfxWideAbsi(tfxWideInt v) {
	return tfxWideAndi(tfxWideShiftRight(tfxWideSetSinglei(-1), 1), v);
}

//----------------------------------------------------------
//Section: Enums
//----------------------------------------------------------

//tfx_graph_t presets to determine limits and scales of different graphs, mainly used for the editor
enum tfx_graph_preset {
	tfxGlobalPercentPreset,
	tfxGlobalOpacityPreset,
	tfxGlobalPercentPresetSigned,
	tfxAnglePreset,
	tfxArcPreset,
	tfxEmissionRangePreset,
	tfxDimensionsPreset,
	tfxTranslationPreset,
	tfxLifePreset,
	tfxAmountPreset,
	tfxVelocityPreset,
	tfxVelocityOvertimePreset,
	tfxWeightPreset,
	tfxWeightVariationPreset,
	tfxNoiseOffsetVariationPreset,
	tfxNoiseResolutionPreset,
	tfxWeightOvertimePreset,
	tfxSpinPreset,
	tfxSpinVariationPreset,
	tfxSpinOvertimePreset,
	tfxDirectionOvertimePreset,
	tfxDirectionVariationPreset,
	tfxFrameratePreset,
	tfxVelocityTurbulancePreset,
	tfxOpacityOvertimePreset,
	tfxColorPreset,
	tfxPercentOvertime,
	tfxIntensityOvertimePreset
};

enum tfx_graph_category : unsigned int {
	tfxGraphCategory_global,
	tfxGraphCategory_transform,
	tfxGraphCategory_property,
	tfxGraphCategory_base,
	tfxGraphCategory_variation,
	tfxGraphCategory_overtime,
	tfxGraphCategory_lifetime,
	tfxGraphCategory_spawn_rate,
	tfxGraphCategory_size,
	tfxGraphCategory_velocity,
	tfxGraphCategory_weight,
	tfxGraphCategory_spin,
	tfxGraphCategory_noise,
	tfxGraphCategory_color,
	tfxGraphCategory_max
};


#define TFX_GLOBAL_COUNT  14
#define	TFX_PROPERTY_COUNT  9
#define	TFX_BASE_COUNT  8
#define	TFX_VARIATION_COUNT  9
#define	TFX_OVERTIME_COUNT  16
#define	TFX_TRANSFORM_COUNT  6

#define TFX_GLOBAL_START 0
#define	TFX_PROPERTY_START TFX_GLOBAL_COUNT
#define	TFX_BASE_START (TFX_PROPERTY_START + TFX_PROPERTY_COUNT)
#define	TFX_VARIATION_START (TFX_BASE_START + TFX_BASE_COUNT)
#define	TFX_OVERTIME_START (TFX_VARIATION_START + TFX_VARIATION_COUNT)
#define	TFX_TRANSFORM_START (TFX_OVERTIME_START + TFX_OVERTIME_COUNT)

//All the different types of graphs, split into main type: global, property, base, variation and overtime
enum tfx_graph_type : unsigned char {
	tfxGlobal_life,
	tfxGlobal_amount,
	tfxGlobal_velocity,
	tfxGlobal_width,
	tfxGlobal_height,
	tfxGlobal_weight,
	tfxGlobal_spin,
	tfxGlobal_stretch,
	tfxGlobal_overal_scale,
	tfxGlobal_intensity,
	tfxGlobal_splatter,
	tfxGlobal_emitter_width,
	tfxGlobal_emitter_height,
	tfxGlobal_emitter_depth,

	tfxProperty_emission_pitch,
	tfxProperty_emission_yaw,
	tfxProperty_emission_range,
	tfxProperty_splatter,
	tfxProperty_emitter_width,
	tfxProperty_emitter_height,
	tfxProperty_emitter_depth,
	tfxProperty_arc_size,
	tfxProperty_arc_offset,

	tfxBase_life,
	tfxBase_amount,
	tfxBase_velocity,
	tfxBase_width,
	tfxBase_height,
	tfxBase_weight,
	tfxBase_spin,
	tfxBase_noise_offset,

	tfxVariation_life,
	tfxVariation_amount,
	tfxVariation_velocity,
	tfxVariation_width,
	tfxVariation_height,
	tfxVariation_weight,
	tfxVariation_spin,
	tfxVariation_noise_offset,
	tfxVariation_noise_resolution,

	tfxOvertime_velocity,
	tfxOvertime_width,
	tfxOvertime_height,
	tfxOvertime_weight,
	tfxOvertime_spin,
	tfxOvertime_stretch,
	tfxOvertime_red,
	tfxOvertime_green,
	tfxOvertime_blue,
	tfxOvertime_blendfactor,
	tfxOvertime_velocity_turbulance,
	tfxOvertime_direction_turbulance,
	tfxOvertime_velocity_adjuster,
	tfxOvertime_intensity,
	tfxOvertime_direction,
	tfxOvertime_noise_resolution,

	tfxTransform_roll,
	tfxTransform_pitch,
	tfxTransform_yaw,
	tfxTransform_translate_x,
	tfxTransform_translate_y,
	tfxTransform_translate_z,
	tfxGraphMaxIndex,
};

//tfx_effect_emitter_t type - effect contains emitters, and emitters spawn particles, but they both share the same struct for simplicity
enum tfx_effect_emitter_type : unsigned char {
	tfxEffectType,
	tfxEmitterType,
	tfxStage,
	tfxFolder
};

//Different ways that particles can be emitted
enum tfx_emission_type : unsigned char {
	tfxPoint,
	tfxArea,
	tfxLine,
	tfxEllipse,
	tfxCylinder,
	tfxIcosphere
};

//Determines how for area, line and ellipse emitters the direction that particles should travel
enum tfx_emission_direction : unsigned char {
	tfxInwards,
	tfxOutwards,
	tfxBothways,
	tfxSpecified
};

//For line effects where traverse line is switched on
enum tfx_line_traversal_end_behaviour : unsigned char {
	tfxLoop,
	tfxKill,
	tfxLetFree
};

//Mainly for the editor, maybe this can just be moved there instead?
enum tfx_export_color_options {
	tfxFullColor,
	tfxOneColor,
	tfxGreyScale
};

//Mainly for the editor, maybe this can just be moved there instead?
enum tfx_export_options {
	tfxSpriteSheet,
	tfxStrip,
	tfxSeparateImages
};

//tfx_graph_t data can be looked up in one of 2 ways, either by just using linear/bezier interpolation (slower), or look up the value in a pre-compiled look up table.
enum tfx_lookup_mode {
	tfxPrecise,
	tfxFast
};

//Used in file loading - for loading effects library
enum tfx_data_type {
	tfxString,
	tfxSInt,
	tfxUint,
	tfxFloat,
	tfxDouble,
	tfxBool,
	tfxUInt64,
	tfxFloat3,
	tfxFloat2
};

//Block designators for loading effects library and other files like animation sprite data
//The values of existing enums below must never change or older files won't load anymore!
enum tfx_effect_library_stream_context : tfxU32 {
	tfxStartEffect = 0x00FFFF00,
	tfxEndEffect,
	tfxStartEmitter,
	tfxEndEmitter,
	tfxStartGraphs,
	tfxEndGraphs,
	tfxStartShapes,
	tfxEndShapes,
	tfxStartAnimationSettings,
	tfxEndAnimationSettings,
	tfxStartImageData,
	tfxStartEffectData,
	tfxEndOfFile,
	tfxStartFolder,
	tfxEndFolder,
	tfxStartPreviewCameraSettings,
	tfxEndPreviewCameraSettings,
	tfxStartStage,
	tfxEndStage,
	tfxStartEffectAnimationInfo,
	tfxEndEffectAnimationInfo,
	tfxStartFrameMeta,
	tfxEndFrameMeta,
	tfxStartFrameOffsets,
	tfxEndFrameOffsets,
};

typedef tfxU32 tfxEmitterPropertyFlags;
typedef tfxU32 tfxEffectPropertyFlags;
typedef tfxU32 tfxVectorFieldFlags;
typedef tfxU32 tfxParticleFlags;
typedef tfxU32 tfxEmitterStateFlags;
typedef tfxU32 tfxEffectStateFlags;
typedef tfxU32 tfxParticleControlFlags;
typedef tfxU32 tfxAttributeNodeFlags;
typedef tfxU32 tfxAngleSettingFlags;
typedef tfxU32 tfxParticleManagerFlags;
typedef tfxU32 tfxErrorFlags;
typedef tfxU32 tfxEffectCloningFlags;
typedef tfxU32 tfxAnimationFlags;
typedef tfxU32 tfxAnimationInstanceFlags;
typedef tfxU32 tfxAnimationManagerFlags;

enum tfx_error_flag_bits {
	tfxErrorCode_success = 0,
	tfxErrorCode_incorrect_package_format = 1 << 0,
	tfxErrorCode_data_could_not_be_loaded = 1 << 1,
	tfxErrorCode_could_not_add_shape = 1 << 2,
	tfxErrorCode_error_loading_shapes = 1 << 3,
	tfxErrorCode_some_data_not_loaded = 1 << 4,
	tfxErrorCode_unable_to_open_file = 1 << 5,
	tfxErrorCode_unable_to_read_file = 1 << 6,
	tfxErrorCode_wrong_file_size = 1 << 7,
	tfxErrorCode_invalid_format = 1 << 8,
	tfxErrorCode_no_inventory = 1 << 9,
	tfxErrorCode_invalid_inventory = 1 << 10,
	tfxErrorCode_sprite_data_is_3d_but_animation_manager_is_2d = 1 << 11,
	tfxErrorCode_sprite_data_is_2d_but_animation_manager_is_3d = 1 << 12,
	tfxErrorCode_library_loaded_without_shape_loader = 1 << 13
};

enum tfx_effect_cloning_flag_bits {
	tfxEffectCloningFlags_none = 0,
	tfxEffectCloningFlags_keep_user_data = 1 << 0,
	tfxEffectCloningFlags_force_clone_global = 1 << 1,
	tfxEffectCloningFlags_clone_graphs = 1 << 2,
	tfxEffectCloningFlags_compile_graphs = 1 << 3
};

enum tfx_particle_manager_mode {
	tfxParticleManagerMode_unordered,
	tfxParticleManagerMode_ordered_by_age,
	tfxParticleManagerMode_ordered_by_depth,
	tfxParticleManagerMode_ordered_by_depth_guaranteed
};

enum tfx_billboarding_option {
	tfxBillboarding_align_to_camera = 0,			//Align to Camera only
	tfxBillboarding_free_align = 1,					//Free align
	tfxBillboarding_align_to_camera_and_vector = 2,	//Align to camera and vector
	tfxBillboarding_align_to_vector = 3				//Align to vector
};

enum tfx_particle_manager_flag_bits {
	tfxEffectManagerFlags_none = 0,
	tfxEffectManagerFlags_disable_spawning = 1,
	tfxEffectManagerFlags_force_capture = 2,			//Unused
	tfxEffectManagerFlags_use_compute_shader = 1 << 3,
	tfxEffectManagerFlags_order_by_depth = 1 << 4,
	tfxEffectManagerFlags_guarantee_order = 1 << 5,
	tfxEffectManagerFlags_update_base_values = 1 << 6,
	tfxEffectManagerFlags_dynamic_sprite_allocation = 1 << 7,
	tfxEffectManagerFlags_3d_effects = 1 << 8,
	tfxEffectManagerFlags_unordered = 1 << 9,
	tfxEffectManagerFlags_ordered_by_age = 1 << 10,
	tfxEffectManagerFlags_update_age_only = 1 << 11,
	tfxEffectManagerFlags_single_threaded = 1 << 12,
	tfxEffectManagerFlags_double_buffer_sprites = 1 << 13,
	tfxEffectManagerFlags_recording_sprites = 1 << 14,
	tfxEffectManagerFlags_using_uids = 1 << 15,
	tfxEffectManagerFlags_2d_and_3d = 1 << 16,
	tfxEffectManagerFlags_update_bounding_boxes = 1 << 17
};

enum tfx_vector_align_type {
	tfxVectorAlignType_motion,
	tfxVectorAlignType_emission,
	tfxVectorAlignType_emitter,
	tfxVectorAlignType_max,
	//Not in yet, need to think about methods of implementing
	tfxVectorAlignType_surface_normal,
};

//Particle property that defines how a particle will rotate
enum tfx_angle_setting_flag_bits {
	tfxAngleSettingFlags_none = 0,														//No flag
	tfxAngleSettingFlags_align_roll = 1 << 0,											//Align the particle with it's direction of travel in 2d
	tfxAngleSettingFlags_random_roll = 1 << 1,											//Chose a random angle at spawn time/state_flags
	tfxAngleSettingFlags_specify_roll = 1 << 2,											//Specify the angle at spawn time
	tfxAngleSettingFlags_align_with_emission = 1 << 3,									//Align the particle with the emission direction only
	tfxAngleSettingFlags_random_pitch = 1 << 4,											//3d mode allows for rotating pitch and yaw when not using billboarding (when particle always faces the camera)
	tfxAngleSettingFlags_random_yaw = 1 << 5,
	tfxAngleSettingFlags_specify_pitch = 1 << 6,
	tfxAngleSettingFlags_specify_yaw = 1 << 7
};

//All the state_flags needed by the ControlParticle function put into one enum to save space
enum tfx_particle_control_flag_bits {
	tfxParticleControlFlags_none = 0,
	tfxParticleControlFlags_random_color = 1 << 0,
	tfxParticleControlFlags_relative_position = 1 << 1,
	tfxParticleControlFlags_relative_angle = 1 << 2,
	tfxParticleControlFlags_point = 1 << 3,
	tfxParticleControlFlags_area = 1 << 4,
	tfxParticleControlFlags_line = 1 << 5,
	tfxParticleControlFlags_ellipse = 1 << 6,
	tfxParticleControlFlags_loop = 1 << 7,
	tfxParticleControlFlags_kill = 1 << 8,
	tfxParticleControlFlags_letFree = 1 << 9,
	tfxParticleControlFlags_edge_traversal = 1 << 10,
	tfxParticleControlFlags_remove = 1 << 11,
	tfxParticleControlFlags_base_uniform_size = 1 << 12,
	tfxParticleControlFlags_lifetime_uniform_size = 1 << 13,
	tfxParticleControlFlags_animate = 1 << 14,
	tfxParticleControlFlags_reverse_animation = 1 << 15,
	tfxParticleControlFlags_play_once = 1 << 16,
	tfxParticleControlFlags_align = 1 << 17,
	tfxParticleControlFlags_emission = 1 << 18,
	tfxParticleControlFlags_random_roll = 1 << 19,
	tfxParticleControlFlags_specify_roll = 1 << 20,
	tfxParticleControlFlags_random_pitch = 1 << 21,
	tfxParticleControlFlags_specify_pitch = 1 << 22,
	tfxParticleControlFlags_random_yaw = 1 << 23,
	tfxParticleControlFlags_specify_yaw = 1 << 24,
};

enum tfx_effect_property_flag_bits {
	tfxEffectPropertyFlags_none = 0,
	tfxEffectPropertyFlags_is_3d = 1 << 0,
	tfxEffectPropertyFlags_depth_draw_order = 1 << 1,
	tfxEffectPropertyFlags_guaranteed_order = 1 << 2,
	tfxEffectPropertyFlags_age_order = 1 << 3,
	tfxEffectPropertyFlags_use_keyframes = 1 << 4,
	tfxEffectPropertyFlags_include_in_sprite_data_export = 1 << 5		//In the editor you can specify which effects you want to be included in a spritedata export
};

enum tfx_emitter_property_flag_bits {
	tfxEmitterPropertyFlags_none = 0,
	tfxEmitterPropertyFlags_random_color = 1 << 0,						//Pick a random color from the color overtime gradient rather then change the color over the lifetime of the particle
	tfxEmitterPropertyFlags_relative_position = 1 << 1,					//Keep the particles position relative to the current position of the emitter
	tfxEmitterPropertyFlags_relative_angle = 1 << 2,					//Keep the angle of the particles relative to the current angle of the emitter
	tfxEmitterPropertyFlags_image_handle_auto_center = 1 << 3,			//Set the offset of the particle to the center of the image
	tfxEmitterPropertyFlags_single = 1 << 4,							//Only spawn a single particle (or number of particles specified by spawn_amount) that does not expire
	tfxEmitterPropertyFlags_specific_emission_direction = 1 << 5,		//Uses a normal vector (3d) or direction (2d) to determine emission direction
	tfxEmitterPropertyFlags_spawn_on_grid = 1 << 6,						//When using an area, line or ellipse emitter, spawn along a grid
	tfxEmitterPropertyFlags_grid_spawn_clockwise = 1 << 7,				//Spawn clockwise/left to right around the area
	tfxEmitterPropertyFlags_fill_area = 1 << 8,							//Fill the area
	tfxEmitterPropertyFlags_emitter_handle_auto_center = 1 << 9,		//Center the handle of the emitter
	tfxEmitterPropertyFlags_edge_traversal = 1 << 10,					//Line emitters only: make particles traverse the line
	tfxEmitterPropertyFlags_global_uniform_size = 1 << 11,				//Keep the global particle size uniform
	tfxEmitterPropertyFlags_base_uniform_size = 1 << 12,				//Keep the base particle size uniform
	tfxEmitterPropertyFlags_lifetime_uniform_size = 1 << 13,			//Keep the size over lifetime of the particle uniform
	tfxEmitterPropertyFlags_animate = 1 << 14,							//Animate the particle shape if it has more than one frame of animation
	tfxEmitterPropertyFlags_reverse_animation = 1 << 15,				//Make the image animation go in reverse
	tfxEmitterPropertyFlags_play_once = 1 << 16,						//Play the animation once only
	tfxEmitterPropertyFlags_random_start_frame = 1 << 17,				//Start the animation of the image from a random frame
	tfxEmitterPropertyFlags_keep_alive = 1 << 18,						//Keep the effect/emitter in the particle manager, don't remove it when it has no particles
	tfxEmitterPropertyFlags_wrap_single_sprite = 1 << 19,				//When recording sprite data, single particles can have their invalid capured index set to the current frame for better looping
	tfxEmitterPropertyFlags_is_in_folder = 1 << 20,						//This effect is located inside a folder
	tfxEmitterPropertyFlags_is_bottom_emitter = 1 << 21,				//This emitter has no child effects, so can spawn particles that could be used in a compute shader if it's enabled
	tfxEmitterPropertyFlags_use_spawn_ratio = 1 << 22,					//Option for area emitters to multiply the amount spawned by a ration of particles per pixels squared
	tfxEmitterPropertyFlags_can_grow_particle_memory = 1 << 23,			//Allows for expanding the memory used for particle emitters if the amount spawned is changed dynamically
	tfxEmitterPropertyFlags_is_3d = 1 << 24,							//Makes the effect run in 3d mode for 3d effects todo: does this need to be here, the effect dictates this?
	tfxEmitterPropertyFlags_use_dynamic = 1 << 25,						//Use a dynamic particle storage rather then a fixed one
	tfxEmitterPropertyFlags_grid_spawn_random = 1 << 26,				//Spawn on grid points but randomly rather then in sequence
	tfxEmitterPropertyFlags_area_open_ends = 1 << 27,					//Only sides of the area/cylinder are spawned on when fill area is not checked
	tfxEmitterPropertyFlags_exclude_from_hue_adjustments = 1 << 28,		//Emitter will be excluded from effect hue adjustments if this flag is checked
	tfxEmitterPropertyFlags_enabled = 1 << 29,							//The emitter is enabled or not, meaning it will or will not be added the particle manager with AddEffect
	tfxEmitterPropertyFlags_match_amount_to_grid_points = 1 << 30,		//Match the amount to spawn with a single emitter to the number of grid points in the effect
	tfxEmitterPropertyFlags_life_proportional_to_animation = 1 << 31	//When recording sprite data and animations, the life particles will be made proportional to the number of frames in the animation
};

enum tfx_particle_flag_bits : unsigned char {
	tfxParticleFlags_none = 0,
	tfxParticleFlags_fresh = 1 << 0,									//Particle has just spawned this frame	
	tfxParticleFlags_capture_after_transform = 1 << 3,					//Particle will be captured after a transfrom, used for traversing lines and looping back to the beginning to avoid lerping imbetween
	tfxParticleFlags_remove = 1 << 4,									//Particle will be removed this or next frame
	tfxParticleFlags_has_velocity = 1 << 5,								//Flagged if the particle is currently moving
	tfxParticleFlags_has_sub_effects = 1 << 6,							//Flagged if the particle has sub effects
};

enum tfx_emitter_state_flag_bits : unsigned int {
	tfxEmitterStateFlags_none = 0,
	tfxEmitterStateFlags_random_color = 1 << 0,
	tfxEmitterStateFlags_relative_position = 1 << 1,					//Keep the particles position relative to the current position of the emitter
	tfxEmitterStateFlags_relative_angle = 1 << 2,						//Keep the angle of the particles relative to the current angle of the emitter
	tfxEmitterStateFlags_stop_spawning = 1 << 3,						//Tells the emitter to stop spawning
	tfxEmitterStateFlags_remove = 1 << 4,								//Tells the effect/emitter to remove itself from the particle manager immediately
	tfxEmitterStateFlags_unused1 = 1 << 5,								//the emitter is enabled. **moved to property state_flags**
	tfxEmitterStateFlags_retain_matrix = 1 << 6,						//Internal flag about matrix usage
	tfxEmitterStateFlags_no_tween_this_update = 1 << 7,					//Internal flag generally, but you could use it if you want to teleport the effect to another location
	tfxEmitterStateFlags_is_single = 1 << 8,
	tfxEmitterStateFlags_not_line = 1 << 9,
	tfxEmitterStateFlags_is_line_traversal = 1 << 10,
	tfxEmitterStateFlags_can_spin = 1 << 11,
	tfxEmitterStateFlags_base_uniform_size = 1 << 12,
	tfxEmitterStateFlags_lifetime_uniform_size = 1 << 13,				//Keep the size over lifetime of the particle uniform
	tfxEmitterStateFlags_loop = 1 << 14,
	tfxEmitterStateFlags_kill = 1 << 15,
	tfxEmitterStateFlags_play_once = 1 << 16,							//Play the animation once only
	tfxEmitterStateFlags_single_shot_done = 1 << 17,
	tfxEmitterStateFlags_is_line_loop_or_kill = 1 << 18,
	tfxEmitterStateFlags_is_area = 1 << 19,
	tfxEmitterStateFlags_no_tween = 1 << 20,
	tfxEmitterStateFlags_align_with_velocity = 1 << 21,
	tfxEmitterStateFlags_is_sub_emitter = 1 << 22,
	tfxEmitterStateFlags_has_noise = 1 << 23
};

enum tfx_effect_state_flag_bits : unsigned int {
	tfxEffectStateFlags_none = 0,
	tfxEffectStateFlags_stop_spawning = 1 << 3,							//Tells the emitter to stop spawning
	tfxEffectStateFlags_remove = 1 << 4,								//Tells the effect/emitter to remove itself from the particle manager immediately
	tfxEffectStateFlags_retain_matrix = 1 << 6,							//Internal flag about matrix usage
	tfxEffectStateFlags_no_tween_this_update = 1 << 7,					//Internal flag generally, but you could use it if you want to teleport the effect to another location
	tfxEffectStateFlags_override_overal_scale = 1 << 8,					//Flagged when the over scale is overridden with SetEffectOveralScale
	tfxEffectStateFlags_override_orientiation = 1 << 9,					//Flagged when any of the effect angles are overridden
	tfxEffectStateFlags_override_size_multiplier = 1 << 10,				//Flagged when any of the effect size multipliers are overridden
	tfxEffectStateFlags_no_tween = 1 << 20
};

enum tfx_vector_field_flag_bits : unsigned char {
	tfxVectorFieldFlags_none = 0,
	tfxVectorFieldFlags_repeat_horizontal = 1 << 0,						//Field will repeat horizontally
	tfxVectorFieldFlags_repeat_vertical = 1 << 1						//Field will repeat vertically
};

enum tfx_attribute_node_flag_bits {
	tfxAttributeNodeFlags_none = 0,
	tfxAttributeNodeFlags_is_curve = 1 << 0,
	tfxAttributeNodeFlags_is_left_curve = 1 << 1,
	tfxAttributeNodeFlags_is_right_curve = 1 << 2,
	tfxAttributeNodeFlags_curves_initialised = 1 << 3
};

enum tfx_animation_flag_bits {
	tfxAnimationFlags_none = 0,
	tfxAnimationFlags_loop = 1 << 0,
	tfxAnimationFlags_seamless = 1 << 1,
	tfxAnimationFlags_needs_recording = 1 << 2,
	tfxAnimationFlags_export_with_transparency = 1 << 3,
	tfxAnimationFlags_auto_set_length = 1 << 4,
	tfxAnimationFlags_orthographic = 1 << 5
};

enum tfx_animation_instance_flag_bits {
	tfxAnimationInstanceFlags_none = 0,
	tfxAnimationInstanceFlags_loop = 1 << 0,
};

enum tfx_animation_manager_flag_bits {
	tfxAnimationManagerFlags_none = 0,
	tfxAnimationManagerFlags_has_animated_shapes = 1 << 0,
	tfxAnimationManagerFlags_initialised = 1 << 1,
	tfxAnimationManagerFlags_is_3d = 1 << 2,
};

//-----------------------------------------------------------
//Constants
//-----------------------------------------------------------

const float tfxMIN_FLOAT = -2147483648.f;
const float tfxMAX_FLOAT = 2147483647.f;
const tfxU32 tfxMAX_UINT = 4294967295;
const int tfxMAX_INT = INT_MAX;
const int tfxMIN_INT = INT_MIN;
const tfxS64 tfxMAX_64i = LLONG_MAX;
const tfxS64 tfxMIN_64i = LLONG_MIN;
const tfxU64 tfxMAX_64u = ULLONG_MAX;
const float tfxGAMMA = 1.f;
#if defined(__x86_64__) || defined(_M_X64)
typedef tfxU64 tfxAddress;
#else
typedef tfxU32 tfxAddress;
#endif

//These constants are the min an max levels for the emitter attribute graphs
const float tfxLIFE_MIN = 0.f;
const float tfxLIFE_MAX = 100000.f;
const float tfxLIFE_STEPS = 200.f;

const float tfxAMOUNT_MIN = 0.f;
const float tfxAMOUNT_MAX = 5000.f;
const float tfxAMOUNT_STEPS = 100.f;

const float tfxGLOBAL_PERCENT_MIN = 0.f;
const float tfxGLOBAL_PERCENT_MAX = 20.f;
const float tfxGLOBAL_PERCENT_STEPS = 100.f;

const float tfxGLOBAL_PERCENT_V_MIN = 0.f;
const float tfxGLOBAL_PERCENT_V_MAX = 10.f;
const float tfxGLOBAL_PERCENT_V_STEPS = 200.f;

const float tfxINTENSITY_MIN = 0.f;
const float tfxINTENSITY_MAX = 5.f;
const float tfxINTENSITY_STEPS = 100.f;

const float tfxANGLE_MIN = 0.f;
const float tfxANGLE_MAX = 1080.f;
const float tfxANGLE_STEPS = 54.f;

const float tfxARC_MIN = 0.f;
const float tfxARC_MAX = 6.28319f;
const float tfxARC_STEPS = .3141595f;

const float tfxEMISSION_RANGE_MIN = 0.f;
const float tfxEMISSION_RANGE_MAX = 180.f;
const float tfxEMISSION_RANGE_STEPS = 30.f;

const float tfxDIMENSIONS_MIN = 0.f;
const float tfxDIMENSIONS_MAX = 4000.f;
const float tfxDIMENSIONS_STEPS = 40.f;

const float tfxVELOCITY_MIN = 0.f;
const float tfxVELOCITY_MAX = 10000.f;
const float tfxVELOCITY_STEPS = 100.f;

const float tfxVELOCITY_OVERTIME_MIN = -20.f;
const float tfxVELOCITY_OVERTIME_MAX = 20.f;
const float tfxVELOCITY_OVERTIME_STEPS = 200.f;

const float tfxWEIGHT_MIN = -2500.f;
const float tfxWEIGHT_MAX = 2500.f;
const float tfxWEIGHT_STEPS = 200.f;

const float tfxWEIGHT_VARIATION_MIN = 0.f;
const float tfxWEIGHT_VARIATION_MAX = 2500.f;
const float tfxWEIGHT_VARIATION_STEPS = 250.f;

const float tfxSPIN_MIN = -2000.f;
const float tfxSPIN_MAX = 2000.f;
const float tfxSPIN_STEPS = 100.f;

const float tfxSPIN_VARIATION_MIN = 0.f;
const float tfxSPIN_VARIATION_MAX = 2000.f;
const float tfxSPIN_VARIATION_STEPS = 100.f;

const float tfxSPIN_OVERTIME_MIN = -20.f;
const float tfxSPIN_OVERTIME_MAX = 20.f;
const float tfxSPIN_OVERTIME_STEPS = 200.f;

const float tfxDIRECTION_OVERTIME_MIN = 0.f;
const float tfxDIRECTION_OVERTIME_MAX = 4320.f;
const float tfxDIRECTION_OVERTIME_STEPS = 216.f;

const float tfxFRAMERATE_MIN = 0.f;
const float tfxFRAMERATE_MAX = 200.f;
const float tfxFRAMERATE_STEPS = 100.f;

const float tfxMAX_DIRECTION_VARIATION = 22.5f;
const float tfxMAX_VELOCITY_VARIATION = 30.f;
const int tfxMOTION_VARIATION_INTERVAL = 30;

//Look up frequency determines the resolution of graphs that are compiled into look up arrays.
static float tfxLOOKUP_FREQUENCY = 10.f;
//Overtime frequency is for lookups that will vary in length depending on the lifetime of the particle. It should generally be a higher resolution than the base graphs
static float tfxLOOKUP_FREQUENCY_OVERTIME = 1.f;

//Look up frequency determines the resolution of graphs that are compiled into look up arrays.
static tfxWideFloat tfxLOOKUP_FREQUENCY_WIDE = tfxWideSetSingle(10.f);
//Overtime frequency is for lookups that will vary in length depending on the lifetime of the particle. It should generally be a higher resolution than the base graphs
static tfxWideFloat tfxLOOKUP_FREQUENCY_OVERTIME_WIDE = tfxWideSetSingle(1.f);

//-----------------------------------------------------------
//Section: String_Buffers
//-----------------------------------------------------------

//Very simple string builder
struct tfx_str_t {
	char *data;
	tfxU32 capacity;
	tfxU32 current_size;
	bool is_local_buffer;

	inline tfx_str_t() : current_size(0), capacity(0), data(nullptr), is_local_buffer(false) {}
	inline ~tfx_str_t() { if (data && !is_local_buffer) { tfxFREE(data); data = nullptr; } current_size = capacity = 0; }

	inline bool			empty() { return current_size == 0; }
	inline char&           operator[](tfxU32 i) { return data[i]; }
	inline const char&     operator[](tfxU32 i) const { assert(i < current_size); return data[i]; }

	inline void         free_all() { if (data) { current_size = capacity = 0; tfxFREE(data); data = nullptr; } }
	inline void         Clear() { current_size = 0; }
	inline char*           begin() { return strbuffer(); }
	inline const char*     begin() const { return strbuffer(); }
	inline char*           end() { return strbuffer() + current_size; }
	inline const char*     end() const { return strbuffer() + current_size; }
	inline char&           back() { assert(current_size > 0); return strbuffer()[current_size - 1]; }
	inline const char&     back() const { assert(current_size > 0); return strbuffer()[current_size - 1]; }
	inline void         pop() { assert(current_size > 0); current_size--; }
	inline void	        push_back(const char v) { if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); new((void*)(data + current_size)) char(v); current_size++; }

	inline tfxU32       _grow_capacity(tfxU32 sz) const { tfxU32 new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
	inline void         resize(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
	inline void         reserve(tfxU32 new_capacity) {
		if (new_capacity <= capacity) return;
		char* new_data = (char*)tfxALLOCATE((size_t)new_capacity * sizeof(char));
		assert(new_data);	//unable to allocate memory. Todo: proper handling
		if (data && !is_local_buffer) {
			memcpy(new_data, data, (size_t)current_size * sizeof(char));
			tfxFREE(data);
		}
		else if (is_local_buffer) {
			memcpy(new_data, strbuffer(), (size_t)current_size * sizeof(char));
		}
		data = new_data;
		is_local_buffer = false;
		capacity = new_capacity;
	}

	tfx_str_t(const char *text) : data(nullptr), current_size(0), capacity(0), is_local_buffer(false) { size_t length = tfx__strlen(text, 512); if (!length) { Clear(); return; }; if (capacity < length) reserve((tfxU32)length); assert(data); memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); }
	tfx_str_t(const tfx_str_t &src) : data(nullptr), current_size(0), capacity(0), is_local_buffer(false) { size_t length = src.Length(); if (!length) { Clear(); return; }; if (capacity < length) reserve((tfxU32)length); assert(data); memcpy(data, src.data, length); current_size = (tfxU32)length; NullTerminate(); }
    inline void operator=(const char *text) { if(!text) { free_all(); return;} size_t length = tfx__strlen(text, 512); if (!length) { Clear(); return; }; if (capacity < length) reserve((tfxU32)length); assert(data); memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); }
	inline void operator=(const tfx_str_t& src) { Clear(); resize(src.current_size); memcpy(data, src.strbuffer(), (size_t)current_size * sizeof(char)); }
	inline bool operator==(const char *string) { return !strcmp(string, c_str()); }
	inline bool operator==(const tfx_str_t string) { return !strcmp(c_str(), string.c_str()); }
	inline bool operator!=(const char *string) { return strcmp(string, c_str()); }
	inline bool operator!=(const tfx_str_t string) { return strcmp(c_str(), string.c_str()); }
	inline const char *strbuffer() const { return is_local_buffer ? (char*)this + sizeof(tfx_str_t) : data; }
	inline char *strbuffer() { return is_local_buffer ? (char*)this + sizeof(tfx_str_t) : data; }
	inline const char *c_str() const { return current_size ? strbuffer() : ""; }
	int Find(const char *needle);
	tfx_str_t Lower();
	inline tfxU32 Length() const { return current_size ? current_size - 1 : 0; }
	void AddLine(const char *format, ...);
	void Setf(const char *format, ...);
	void Appendf(const char* format, ...);
	void Appendv(const char* format, va_list args);
	inline void Append(char c) {
		if (current_size) {
			pop();
		}
		push_back(c);
		NullTerminate();
	}
	inline void Pop() {
		if (!Length()) return;
		if (back() == 0)
			pop();
		pop();
		NullTerminate();
	}
	inline void Trim(char c = ' ') {
		if (!Length()) return;
		if (back() == 0)
			pop();
		while (back() == c && current_size) {
			pop();
		}
		NullTerminate();
	}
	inline void TrimFront(char c = ' ') {
		if (!Length()) return;
		tfxU32 pos = 0;
		while (strbuffer()[pos] == c && pos < current_size) {
			pos++;
		}
		if (pos < current_size) {
			memcpy(strbuffer(), strbuffer() + pos, current_size - pos);
		}
		current_size -= pos;
	}
	void NullTerminate() { push_back('\0'); }
	bool SaveToFile(const char *file_name);
	const bool IsInt() const;
	const bool IsFloat() const;
} TFX_PACKED_STRUCT;

#define tfxStrType(type, size)		\
	struct type : public tfx_str_t { \
	char buffer[size];				\
	type() { memset(buffer, 0, size); data = buffer; capacity = size; current_size = 0; is_local_buffer = true; NullTerminate(); } \
	inline void operator=(const tfx_str_t& src) { \
		data = buffer; \
		is_local_buffer = true; \
		capacity = size; size_t length = src.Length(); \
		if (!length) { \
			Clear(); return; \
		}; \
		resize(src.current_size); \
		memcpy(data, src.strbuffer(), length); \
		current_size = (tfxU32)length; \
		NullTerminate(); \
	} \
	inline void operator=(const type& src) { \
		data = buffer; \
		is_local_buffer = true; \
		capacity = size; size_t length = src.Length(); \
		if (!length) { \
			Clear(); return; \
		}; \
		resize(src.current_size); \
		memcpy(data, src.strbuffer(), length); \
		current_size = (tfxU32)length; \
		NullTerminate(); \
	} \
	inline void operator=(const char *text) { data = buffer; is_local_buffer = true; capacity = size; size_t length = tfx__strlen(text, size); if (!length) { Clear(); return; } memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); } \
	type(const char *text) { memset(buffer, 0, size); data = buffer; is_local_buffer = true; capacity = size; size_t length = tfx__strlen(text, size); if (!length) { Clear(); return; } memcpy(data, text, length); current_size = (tfxU32)length; NullTerminate(); } \
	type(const tfx_str_t &src) { \
		memset(buffer, 0, size); \
		data = buffer; \
		is_local_buffer = true; \
		capacity = size; size_t length = src.Length(); \
		if (!length) { \
			Clear(); return; \
		}; \
		resize(src.current_size); \
		memcpy(data, src.strbuffer(), length); \
		current_size = (tfxU32)length; \
		NullTerminate(); \
	} \
	type(const type &src) { \
		memset(buffer, 0, size); \
		data = buffer; \
		is_local_buffer = true; \
		capacity = size; size_t length = src.Length(); \
		if (!length) { \
			Clear(); return; \
		}; \
		resize(src.current_size); \
		memcpy(data, src.strbuffer(), length); \
		current_size = (tfxU32)length; \
		NullTerminate(); \
	} \
	inline int Find(const char *needle) { type compare = needle; type lower = Lower(); compare = compare.Lower(); if (compare.Length() > Length()) return -1; tfxU32 pos = 0; int found = 0; while (compare.Length() + pos <= Length()) { if (strncmp(lower.data + pos, compare.data, compare.Length()) == 0) { return pos; } ++pos; } return -1; } \
	inline type Lower() { type convert = *this; for (auto &c : convert) { c = tolower(c); } return convert; } \
	} TFX_PACKED_STRUCT;

tfxStrType(tfx_str512_t, 512);
tfxStrType(tfx_str256_t, 256);
tfxStrType(tfx_str128_t, 128);
tfxStrType(tfx_str64_t, 64);
tfxStrType(tfx_str32_t, 32);
tfxStrType(tfx_str16_t, 16);

//-----------------------------------------------------------
//Containers_and_Memory
//-----------------------------------------------------------

//Storage
//Credit to ocornut https://github.com/ocornut/imgui/commits?author=ocornut for tfxvec
//std::vector replacement with some extra stuff and tweaks specific to TimelineFX
template<typename T>
struct tfx_vector_t {
	tfxU32 current_size;
	tfxU32 capacity;
	tfxU32 volatile locked;
	tfxU32 alignment;
	T* data;

	// Provide standard typedefs but we don't use them ourselves.
	typedef T                   value_type;
	typedef value_type*         iterator;
	typedef const value_type*   const_iterator;

	inline tfx_vector_t() { locked = false; current_size = capacity = alignment = 0; data = nullptr; tfxINIT_VEC_NAME; }
	inline tfx_vector_t(const char *name_init) { locked = false; current_size = capacity = alignment = 0; data = nullptr; tfxINIT_VEC_NAME_INIT(name_init); }
	inline tfx_vector_t(const tfx_vector_t<T> &src) { locked = false; current_size = capacity = alignment = 0; data = nullptr; tfxINIT_VEC_NAME_SRC_COPY; resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); }
	inline tfx_vector_t<T>& operator=(const tfx_vector_t<T>& src) { clear(); resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); return *this; }
	inline ~tfx_vector_t() { if (data) { tfxFREE(data); } data = nullptr; current_size = capacity = alignment = 0; }

	inline bool			empty() { return current_size == 0; }
	inline bool			full() { return current_size == capacity; }
	inline tfxU32		size() { return current_size; }
	inline const tfxU32	size() const { return current_size; }
	inline tfxU32		size_in_bytes() { return current_size * sizeof(T); }
	inline const tfxU32	size_in_bytes() const { return current_size * sizeof(T); }
	inline T&           operator[](tfxU32 i) { return data[i]; }
	inline const T&     operator[](tfxU32 i) const { assert(i < current_size); return data[i]; }
	inline T&           ts_at(tfxU32 i) { while (locked > 0); return data[i]; }

	inline void         free_all() { if (data) { current_size = capacity = alignment = 0; tfxFREE(data); data = nullptr; } }
	inline void         free() { if (data) { current_size = capacity = alignment = 0; tfxFREE(data); data = nullptr; } }
	inline void         clear() { if (data) { current_size = 0; } }
	inline T*           begin() { return data; }
	inline const T*     begin() const { return data; }
	inline T*           end() { return data + current_size; }
	inline const T*     end() const { return data + current_size; }
	inline T*           rend() { return data; }
	inline const T*     rend() const { return data; }
	inline T*           rbegin() { return data + current_size; }
	inline const T*     rbegin() const { return data + current_size; }
	inline T&           front() { assert(current_size > 0); return data[0]; }
	inline const T&     front() const { assert(current_size > 0); return data[0]; }
	inline T&           back() { assert(current_size > 0); return data[current_size - 1]; }
	inline const T&     back() const { assert(current_size > 0); return data[current_size - 1]; }
	inline T&           parent() { assert(current_size > 1); return data[current_size - 2]; }
	inline const T&     parent() const { assert(current_size > 1); return data[current_size - 2]; }
	inline tfxU32       _grow_capacity(tfxU32 sz) const { tfxU32 new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
	inline void         resize(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
	inline void         resize_bytes(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
	inline void         resize(tfxU32 new_size, const T& v) { if (new_size > capacity) reserve(_grow_capacity(new_size)); if (new_size > current_size) for (tfxU32 n = current_size; n < new_size; n++) memcpy(&data[n], &v, sizeof(v)); current_size = new_size; }
	inline void         shrink(tfxU32 new_size) { assert(new_size <= current_size); current_size = new_size; }
	inline void			set_alignment(tfxU32 align_to) { TFX_ASSERT(0 == (align_to & (align_to - 1)) && "must align to a power of two"); alignment = align_to; }
	inline void         reserve(tfxU32 new_capacity) {
		if (new_capacity <= capacity)
			return;
		T* new_data;
		if (alignment != 0) {
			new_data = (T*)tfxALLOCATE_ALIGNED((size_t)new_capacity * sizeof(T), alignment);
		}
		else {
			new_data = (T*)tfxALLOCATE((size_t)new_capacity * sizeof(T));
		}
		assert(new_data);	//Unable to allocate memory. todo: better handling
		if (data) {
			memcpy(new_data, data, (size_t)current_size * sizeof(T));
			tfxFREE(data);
		}
		data = new_data;
		capacity = new_capacity;
	}

	inline T&	        grab() {
		if (current_size == capacity) reserve(_grow_capacity(current_size + 1));
		current_size++;
		return data[current_size - 1];
	}
	inline tfxU32        locked_push_back(const T& v) {
		//suspect, just use a mutex instead?
		while (tfx__compare_and_exchange((tfxLONG volatile*)&locked, 1, 0) > 1);
		if (current_size == capacity)
			reserve(_grow_capacity(current_size + 1));
		new((void*)(data + current_size)) T(v);
		tfxU32 index = current_size++;
		tfx__exchange((tfxLONG volatile*)&locked, 0);
		return index;
	}
	inline T&	        push_back(const T& v) {
		if (current_size == capacity) {
			reserve(_grow_capacity(current_size + 1));
		}
		new((void*)(data + current_size)) T(v);
		//memcpy(&data[current_size], &v, sizeof(T));
		current_size++; return data[current_size - 1];
	}
	inline T&	        push_back_copy(const T& v) {
		if (current_size == capacity)
			reserve(_grow_capacity(current_size + 1));
		memcpy(&data[current_size], &v, sizeof(v));
		current_size++; return data[current_size - 1];
	}
	inline T&			next() {
		return push_back(T());
	}
	inline void			zero() { assert(capacity > 0); memset(data, 0, capacity * sizeof(T)); }
	inline void         pop() { assert(current_size > 0); current_size--; }
	inline T&	        pop_back() { assert(current_size > 0); current_size--; return data[current_size]; }
	inline void         push_front(const T& v) { if (current_size == 0) push_back(v); else insert(data, v); }
	inline T*           erase(const T* it) { assert(it >= data && it < data + current_size); const ptrdiff_t off = it - data; memmove(data + off, data + off + 1, ((size_t)current_size - (size_t)off - 1) * sizeof(T)); current_size--; return data + off; }
	inline T	        pop_front() { assert(current_size > 0); T front = data[0]; erase(data); return front; }
	inline T*           erase(const T* it, const T* it_last) { assert(it >= data && it < data + current_size && it_last > it && it_last <= data + current_size); const ptrdiff_t count = it_last - it; const ptrdiff_t off = it - data; memmove(data + off, data + off + count, ((size_t)current_size - (size_t)off - count) * sizeof(T)); current_size -= (tfxU32)count; return data + off; }
	inline T*           erase_unsorted(const T* it) { assert(it >= data && it < data + current_size);  const ptrdiff_t off = it - data; if (it < data + current_size - 1) memcpy(data + off, data + current_size - 1, sizeof(T)); current_size--; return data + off; }
	inline T*           insert(const T* it, const T& v) { assert(it >= data && it <= data + current_size); const ptrdiff_t off = it - data; if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); if (off < (ptrdiff_t)current_size) memmove(data + off + 1, data + off, ((size_t)current_size - (size_t)off) * sizeof(T)); new((void*)(data + off)) T(v); current_size++; return data + off; }
	inline T*           insert_after(const T* it, const T& v) { assert(it >= data && it <= data + current_size); const ptrdiff_t off = (it + 1) - data; if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); if (off < (ptrdiff_t)current_size) memmove(data + off + 1, data + off, ((size_t)current_size - (size_t)off) * sizeof(T)); new((void*)(data + off)) T(v); current_size++; return data + off; }
	inline bool         contains(const T& v) const { const T* _data = data;  const T* data_end = data + current_size; while (_data < data_end) if (*_data++ == v) return true; return false; }
	inline T*           find(const T& v) { T* _data = data;  const T* data_end = data + current_size; while (_data < data_end) if (*_data == v) break; else ++_data; return _data; }
	inline const T*     find(const T& v) const { const T* _data = data;  const T* data_end = data + current_size; while (_data < data_end) if (*_data == v) break; else ++_data; return _data; }
	inline bool         find_erase(const T& v) { const T* it = find(v); if (it < data + current_size) { erase(it); return true; } return false; }
	inline bool         find_erase_unsorted(const T& v) { const T* it = find(v); if (it < data + current_size) { erase_unsorted(it); return true; } return false; }
	inline tfxU32       index_from_ptr(const T* it) const { assert(it >= data && it < data + current_size); const ptrdiff_t off = it - data; return (tfxU32)off; }

	inline void			create_pool(tfxU32 amount) { assert(current_size == 0); T base; reserve(amount); for (tfxU32 i = 0; i != capacity; ++i) { new((void*)(data + current_size)) T(base); current_size++; } }
	inline void			create_pool_with(tfxU32 amount, const T &base) { assert(current_size == 0);  reserve(amount); for (tfxU32 i = 0; i != capacity; ++i) { new((void*)(data + current_size)) T(base); current_size++; } }

};

//Simple storage map for storing things by key/pair. The data will be in order that you add items, but the map will be in key order so just do a foreach on the data
//and use At() to retrieve data items by name use [] overload to fetch by index if you have that.
//Should not be used to constantly insert/remove things every frame, it's designed for setting up lists and fetching values in loops (by index preferably), and modifying based on user interaction or setting up new situation.
//Note that if you reference things by index and you then remove something then that index may not be valid anymore so you would need to keep checks on that.
//Not sure how efficient a hash lookup with this is, could probably be better, but isn't used much at all in any realtime particle updating.
template<typename T>
struct tfx_storage_map_t {
	struct pair {
		tfxKey key;
		tfxU32 index;
		pair(tfxKey k, tfxU32 i) : key(k), index(i) {}
	};

	tfx_vector_t<pair> map;
	tfx_vector_t<T> data;
	void(*remove_callback)(T &item) = nullptr;

	tfx_storage_map_t() : map(tfxCONSTRUCTOR_VEC_INIT("Storage Map map")), data(tfxCONSTRUCTOR_VEC_INIT("Storage Map data")) {}
	tfx_storage_map_t(const char *map_tracker, const char *data_tracker) : map(tfxCONSTRUCTOR_VEC_INIT2(map_tracker)), data(tfxCONSTRUCTOR_VEC_INIT2(data_tracker)) {}

	//Insert a new T value into the storage
	inline void Insert(const char *name, const T &value) {
		tfxKey key = tfxXXHash64::hash(name, strlen(name), 0);
		SetIndex(key, value);
	}

	//Insert a new T value into the storage
	inline void Insert(tfx_str_t name, const T &value) {
		tfxKey key = tfxXXHash64::hash(name.c_str(), name.Length(), 0);
		SetIndex(key, value);
	}

	//Insert a new T value into the storage
	void Insert(tfxKey key, const T &value) {
		SetIndex(key, value);
	}

	//Insert a new T value into the storage
	inline void InsertByInt(int name, const T &value) {
		tfxKey key = name;
		SetIndex(key, value);
	}

	inline void Clear() {
		data.clear();
		map.clear();
	}

	inline tfxKey MakeKey(const char* name) {
		return tfxXXHash64::hash(name, strlen(name), 0);
	}

	inline void FreeAll() {
		data.free_all();
		map.free_all();
	}

	inline tfxU32 Size() {
		return data.current_size;
	}

	inline tfxU32 LastIndex() {
		return data.current_size - 1;
	}

	inline bool ValidIndex(tfxU32 index) {
		return index < data.current_size;
	}

	inline bool ValidName(const char *name) {
		assert(name);	//Can't search for anything that's null
		return GetIndex(name) > -1;
	}

	inline bool ValidKey(tfxKey key) {
		return GetIndex(key) > -1;
	}

	inline bool ValidIntName(tfxU32 name) {
		return GetIntIndex(name) > -1;
	}

	inline bool ValidName(const tfx_str_t &name) {
		return GetIndex(name) > -1;
	}

	//Remove an item from the data. Slow function, 2 memmoves and then the map has to be iterated and indexes reduced by one
	//to re-align them
	inline void Remove(const char *name) {
		tfxKey key = tfxXXHash64::hash(name, strlen(name), 0);
		pair *it = LowerBound(key);
		if (remove_callback)
			remove_callback(data[it->index]);
		tfxU32 index = it->index;
		T* list_it = &data[index];
		map.erase(it);
		data.erase(list_it);
		for (auto &p : map) {
			if (p.index < index) continue;
			p.index--;
		}
	}

	//Remove an item from the data. Slow function, 2 memmoves and then the map has to be iterated and indexes reduced by one
	//to re-align them
	inline void Remove(const tfxKey &key) {
		pair *it = LowerBound(key);
		if (remove_callback)
			remove_callback(data[it->index]);
		tfxU32 index = it->index;
		T* list_it = &data[index];
		map.erase(it);
		data.erase(list_it);
		for (auto &p : map) {
			if (p.index < index) continue;
			p.index--;
		}
	}

	//Remove an item from the data. Slow function, 2 memmoves and then the map has to be iterated and indexes reduced by one
	//to re-align them
	inline void RemoveInt(int name) {
		tfxKey key = name;
		pair *it = LowerBound(key);
		if (remove_callback)
			remove_callback(data[it->index]);
		tfxU32 index = it->index;
		T* list_it = &data[index];
		map.erase(it);
		data.erase(list_it);
		for (auto &p : map) {
			if (p.index < index) continue;
			p.index--;
		}
	}

	inline T &At(const char *name) {
		int index = GetIndex(name);
		assert(index > -1);						//Key was not found
		return data[index];
	}

	inline T &At(const tfx_str_t &name) {
		int index = GetIndex(name.c_str());
		assert(index > -1);						//Key was not found
		return data[index];
	}

	inline T &AtInt(int name) {
		int index = GetIntIndex(name);
		assert(index > -1);						//Key was not found
		return data[index];
	}

	inline T &At(tfxKey key) {
		int index = GetIndex(key);
		assert(index > -1);						//Key was not found
		return data[index];
	}

	inline T &operator[](const tfxU32 index) {
		assert(index < data.current_size);		//Index was out of range
		return data[index];
	}

	void SetIndex(tfxKey key, const T &value) {
		pair* it = LowerBound(key);
		if (it == map.end() || it->key != key)
		{
			data.push_back(value);
			map.insert(it, pair(key, data.current_size - 1));
			return;
		}
		data[it->index] = value;
	}

	int GetIndex(const char *name) {
		tfxKey key = tfxXXHash64::hash(name, strlen(name), 0);
		pair* it = LowerBound(key);
		if (it == map.end() || it->key != key)
			return -1;
		return it->index;
	}

	int GetIntIndex(int name) {
		tfxKey key = name;
		pair* it = LowerBound(key);
		if (it == map.end() || it->key != key)
			return -1;
		return it->index;
	}

	int GetIndex(const tfx_str_t &name) {
		tfxKey key = tfxXXHash64::hash(name.c_str(), name.Length(), 0);
		pair* it = LowerBound(key);
		if (it == map.end() || it->key != key)
			return -1;
		return it->index;
	}

	int GetIndex(tfxKey key) {
		pair* it = LowerBound(key);
		if (it == map.end() || it->key != key)
			return -1;
		return it->index;
	}

	pair* LowerBound(tfxKey key)
	{
		tfx_storage_map_t::pair* first = map.data;
		tfx_storage_map_t::pair* last = map.data + map.current_size;
		size_t count = (size_t)(last - first);
		while (count > 0)
		{
			size_t count2 = count >> 1;
			tfx_storage_map_t::pair* mid = first + count2;
			if (mid->key < key)
			{
				first = ++mid;
				count -= count2 + 1;
			}
			else
			{
				count = count2;
			}
		}
		return first;
	}

};

#define tfxKilobyte(Value) ((Value)*1024LL)
#define tfxMegabyte(Value) (tfxKilobyte(Value)*1024LL)
#define tfxGigabyte(Value) (tfxMegabyte(Value)*1024LL)

#ifndef tfxSTACK_SIZE
#define tfxSTACK_SIZE tfxMegabyte(2)
#endif

#ifndef tfxMT_STACK_SIZE
#define tfxMT_STACK_SIZE tfxMegabyte(4)
#endif

inline tfxU32 tfxIsPowerOf2(tfxU32 v)
{
	return ((v & ~(v - 1)) == v);
}

//Used in tfx_soa_buffer_t to store pointers to arrays inside a struct of arrays
struct tfx_soa_data_t {
	void *ptr = nullptr;		//A pointer to the array in the struct
	size_t offset = 0;		//The offset to the memory location in the buffer where the array starts
	size_t unit_size = 0;	//The size of each data type in the array
};

//A buffer designed to contain structs of arrays. If the arrays need to grow then a new memory block is made and all copied over
//together. All arrays in the struct will be the same capacity but can all have different unit sizes/types.
//In order to use this you need to first prepare the buffer by calling AddStructArray for each struct member of the SoA you're setting up. 
//All members must be of the same struct.
//Then call FinishSoABufferSetup to create the memory for the struct of arrays with an initial reserve amount.
struct tfx_soa_buffer_t {
	size_t current_arena_size = 0;		//The current size of the arena that contains all the arrays
	size_t struct_size = 0;				//The size of the struct (each member unit size added up)
	tfxU32 current_size = 0;			//current size of each array
	tfxU32 start_index = 0;				//Start index if you're using the buffer as a ring buffer
	tfxU32 last_bump = 0;				//the amount the the start index was bumped by the last time Bump was called
	tfxU32 capacity = 0;				//capacity of each array
	tfxU32 block_size = tfxDataWidth;	//Keep the capacity to the nearest block size
	tfxU32 alignment = 4;				//The alignment of the memory. If you're planning on using simd for the memory, then 16 will be necessary.
	tfx_vector_t<tfx_soa_data_t> array_ptrs;		//Container for all the pointers into the arena
	void *user_data = nullptr;
	void(*resize_callback)(tfx_soa_buffer_t *ring, tfxU32 new_index_start) = nullptr;
	void *struct_of_arrays = nullptr;		//Pointer to the struct of arrays. Important that this is a stable pointer! Set with FinishSoABufferSetup
	void *data = nullptr;					//Pointer to the area in memory that contains all of the array data	
};

inline void ResetSoABuffer(tfx_soa_buffer_t *buffer) {
	buffer->current_arena_size = 0;
	buffer->struct_size = 0;
	buffer->current_size = 0;
	buffer->start_index = 0;
	buffer->last_bump = 0;
	buffer->capacity = 0;
	buffer->block_size = tfxDataWidth;
	buffer->user_data = nullptr;
	buffer->resize_callback = nullptr;
	buffer->struct_of_arrays = nullptr;
	buffer->data = nullptr;
}

inline void* GetEndOfBufferPtr(tfx_soa_buffer_t *buffer) {
	assert(buffer->data);
	return (char*)buffer->data + buffer->current_arena_size;
}

//Get the amount of free space in the buffer
inline tfxU32 FreeSpriteBufferSpace(tfx_soa_buffer_t *buffer) {
	return buffer->capacity - buffer->current_size;
}

//Get the index based on the buffer being a ring buffer
inline tfxU32 GetCircularIndex(tfx_soa_buffer_t *buffer, tfxU32 index) {
	return (buffer->start_index + index) % buffer->capacity;
}

//Get the index based on the buffer being a ring buffer
inline tfxU32 GetAbsoluteIndex(tfx_soa_buffer_t *buffer, tfxU32 circular_index) {
	return buffer->capacity - (circular_index % buffer->capacity);
}

//Add an array to a SoABuffer. parse in the size of the data type and the offset to the member variable within the struct.
//You must add all the member veriables in the struct before calling FinishSoABufferSetup
inline void AddStructArray(tfx_soa_buffer_t *buffer, size_t unit_size, size_t offset) {
	tfx_soa_data_t data;
	data.unit_size = unit_size;
	data.offset = offset;
	buffer->array_ptrs.push_back(data);
}

//In order to ensure memory alignment of all arrays in the buffer we need the following function to get the correct amount
//of memory required to contain all the data in the buffer for each array in the struct of arrays.
inline size_t GetSoACapacityRequirement(tfx_soa_buffer_t *buffer, size_t capacity) {
	size_t size_requirement = 0;
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		size_requirement += buffer->array_ptrs[i].unit_size * capacity;
		size_requirement += buffer->alignment - (size_requirement % buffer->alignment);
	}
	return size_requirement;
}

//Once you have called AddStructArray for all your member variables you must call this function in order to 
//set up the memory for all your arrays. One block of memory will be created and all your arrays will be line up
//inside the space
inline void FinishSoABufferSetup(tfx_soa_buffer_t *buffer, void *struct_of_arrays, tfxU32 reserve_amount, tfxU32 alignment = 4) {
	assert(buffer->data == nullptr && buffer->array_ptrs.current_size > 0);	//Must be an unitialised soa buffer
	assert(alignment >= 4);		//Alignment must be 4 or greater
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		buffer->struct_size += buffer->array_ptrs[i].unit_size;
	}
	reserve_amount = (reserve_amount / buffer->block_size + 1) * buffer->block_size;
	buffer->capacity = reserve_amount;
	buffer->alignment = alignment;
	buffer->current_arena_size = GetSoACapacityRequirement(buffer, reserve_amount);
	buffer->data = tfxALLOCATE_ALIGNED(buffer->current_arena_size, buffer->alignment);
	assert(buffer->data);	//Unable to allocate memory. Todo: better handling
	memset(buffer->data, 0, buffer->current_arena_size);
	buffer->struct_of_arrays = struct_of_arrays;
	size_t running_offset = 0;
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		buffer->array_ptrs[i].ptr = (char*)buffer->data + running_offset;
		memcpy((char*)buffer->struct_of_arrays + buffer->array_ptrs[i].offset, &buffer->array_ptrs[i].ptr, sizeof(void*));
		running_offset += buffer->array_ptrs[i].unit_size * buffer->capacity;
		running_offset += buffer->alignment - (running_offset % buffer->alignment);
	}
	if (buffer->resize_callback) {
		buffer->resize_callback(buffer, 0);
	}
}

//Call this function to increase the capacity of all the arrays in the buffer. Data that is already in the arrays is preserved if keep_data passed as true (default).
inline bool GrowArrays(tfx_soa_buffer_t *buffer, tfxU32 first_new_index, tfxU32 new_target_size, bool keep_data = true) {
	assert(buffer->capacity);			//buffer must already have a capacity!
	tfxU32 new_capacity = 0;
	new_capacity = new_target_size > buffer->capacity ? new_target_size + new_target_size / 2 : buffer->capacity + buffer->capacity / 2;
	new_capacity = (new_capacity / buffer->block_size + 1) * buffer->block_size;
	void *new_data = tfxALLOCATE_ALIGNED(GetSoACapacityRequirement(buffer, new_capacity), buffer->alignment);
	assert(new_data);	//Unable to allocate memory. Todo: better handling
	memset(new_data, 0, new_capacity * buffer->struct_size);
	size_t running_offset = 0;
	if (keep_data) {
		for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
			size_t capacity = buffer->capacity * buffer->array_ptrs[i].unit_size;
			size_t start_index = buffer->start_index * buffer->array_ptrs[i].unit_size;
			if ((buffer->start_index + buffer->current_size - 1) > buffer->capacity) {
				memcpy((char*)new_data + running_offset, (char*)buffer->array_ptrs[i].ptr + start_index, (size_t)(capacity - start_index));
				memcpy((char*)new_data + (capacity - start_index) + running_offset, (char*)buffer->array_ptrs[i].ptr, (size_t)(start_index));
			}
			else {
				memcpy((char*)new_data + running_offset, (char*)buffer->array_ptrs[i].ptr + start_index, (size_t)(capacity - start_index));
			}
			running_offset += buffer->array_ptrs[i].unit_size * new_capacity;
			running_offset += buffer->alignment - (running_offset % buffer->alignment);
		}
	}
	void *old_data = buffer->data;

	buffer->data = new_data;
	buffer->capacity = new_capacity;
	buffer->current_arena_size = new_capacity * buffer->struct_size;
	running_offset = 0;
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		buffer->array_ptrs[i].ptr = (char*)buffer->data + running_offset;
		memcpy((char*)buffer->struct_of_arrays + buffer->array_ptrs[i].offset, &buffer->array_ptrs[i].ptr, sizeof(void*));
		running_offset += buffer->array_ptrs[i].unit_size * buffer->capacity;
		running_offset += buffer->alignment - (running_offset % buffer->alignment);
	}
	tfxFREE(old_data);

	if (buffer->resize_callback) {
		buffer->resize_callback(buffer, first_new_index);
	}

	buffer->start_index = 0;

	return true;
}

//Increase current size of a SoA Buffer and grow if necessary.
inline void Resize(tfx_soa_buffer_t *buffer, tfxU32 new_size) {
	assert(buffer->data);			//No data allocated in buffer
	if (new_size >= buffer->capacity) {
		GrowArrays(buffer, buffer->capacity, new_size);
	}
	buffer->current_size = new_size;
}

//Increase current size of a SoA Buffer and grow if necessary. This will not shrink the capacity so if new_size is not bigger than the
//current capacity then nothing will happen
inline void SetCapacity(tfx_soa_buffer_t *buffer, tfxU32 new_size) {
	assert(buffer->data);			//No data allocated in buffer
	if (new_size >= buffer->capacity) {
		GrowArrays(buffer, buffer->capacity, new_size);
	}
}

//Increase current size of a SoA Buffer by 1 and grow if grow is true. Returns the last index.
inline tfxU32 AddRow(tfx_soa_buffer_t *buffer, bool grow = false) {
	assert(buffer->data);			//No data allocated in buffer
	tfxU32 new_size = ++buffer->current_size;
	if (grow && new_size == buffer->capacity) {
		GrowArrays(buffer, buffer->capacity, new_size);
	}
	buffer->current_size = new_size;
	assert(buffer->current_size <= buffer->capacity);	//Capacity of buffer is exceeded, set grow to true or don't exceed the capacity
	return buffer->current_size - 1;
}

//Increase current size of a SoA Buffer by a specific amount and grow if grow is true. Returns the last index.
//You can also pass in a boolean to know if the buffer had to be increased in size or not. Returns the index where the new rows start.
inline tfxU32 AddRows(tfx_soa_buffer_t *buffer, tfxU32 amount, bool grow, bool &grew) {
	assert(buffer->data);			//No data allocated in buffer
	tfxU32 first_new_index = buffer->current_size;
	tfxU32 new_size = buffer->current_size += amount;
	if (grow && new_size >= buffer->capacity) {
		grew = GrowArrays(buffer, buffer->capacity, new_size);
	}
	buffer->current_size = new_size;
	assert(buffer->current_size < buffer->capacity);	//Capacity of buffer is exceeded, set grow to true or don't exceed the capacity
	return first_new_index;
}

//Increase current size of a SoA Buffer and grow if grow is true. Returns the index where the new rows start.
inline tfxU32 AddRows(tfx_soa_buffer_t *buffer, tfxU32 amount, bool grow) {
	assert(buffer->data);			//No data allocated in buffer
	tfxU32 first_new_index = buffer->current_size;
	tfxU32 new_size = buffer->current_size + amount;
	if (grow && new_size >= buffer->capacity) {
		GrowArrays(buffer, buffer->capacity, new_size);
	}
	buffer->current_size = new_size;
	assert(buffer->current_size < buffer->capacity);	//Capacity of buffer is exceeded, set grow to true or don't exceed the capacity
	return first_new_index;
}

//Decrease the current size of a SoA Buffer by 1.
inline void PopRow(tfx_soa_buffer_t *buffer) {
	assert(buffer->data && buffer->current_size > 0);			//No data allocated in buffer
	buffer->current_size--;
}

//Bump the start index of the SoA buffer (ring buffer usage)
inline void Bump(tfx_soa_buffer_t *buffer) {
	assert(buffer->data && buffer->current_size > 0);			//No data allocated in buffer
	if (buffer->current_size == 0)
		return;
	buffer->start_index++; buffer->start_index %= buffer->capacity; buffer->current_size--;
}

//Bump the start index of the SoA buffer (ring buffer usage)
inline void Bump(tfx_soa_buffer_t *buffer, tfxU32 amount) {
	assert(buffer->data && buffer->current_size > 0);			//No data allocated in buffer
	if (buffer->current_size == 0)
		return;
	if (amount > buffer->current_size)
		amount = buffer->current_size;
	buffer->start_index += amount;
	buffer->start_index %= buffer->capacity;
	buffer->current_size -= amount;
	buffer->last_bump = amount;
}

//Free the SoA buffer
inline void FreeSoABuffer(tfx_soa_buffer_t *buffer) {
	buffer->current_arena_size = buffer->current_size = buffer->capacity = 0;
	if (buffer->data)
		tfxFREE(buffer->data);
	buffer->array_ptrs.free_all();
	ResetSoABuffer(buffer);
}

//Clear the SoA buffer
inline void ClearSoABuffer(tfx_soa_buffer_t *buffer) {
	buffer->current_size = buffer->start_index = 0;
}

//Trim an SoA buffer to the current size. This is a bit rough and ready and I just created it for trimming compressed sprite data down to size
inline void TrimSoABuffer(tfx_soa_buffer_t *buffer) {
	if (buffer->current_size == buffer->capacity) {
		return;
	}
	if (buffer->current_size == 0) {
		FreeSoABuffer(buffer);
		return;
	}
	assert(buffer->current_size < buffer->capacity);
	tfxU32 new_capacity = buffer->current_size;
	void *new_data = tfxALLOCATE_ALIGNED(GetSoACapacityRequirement(buffer, new_capacity), buffer->alignment);
	assert(new_data);	//Unable to allocate memory. Todo: better handling
	memset(new_data, 0, new_capacity * buffer->struct_size);
	size_t running_offset = 0;
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		size_t capacity = new_capacity * buffer->array_ptrs[i].unit_size;
		size_t start_index = buffer->start_index * buffer->array_ptrs[i].unit_size;
		if ((buffer->start_index + buffer->current_size - 1) > buffer->capacity) {
			memcpy((char*)new_data + running_offset, (char*)buffer->array_ptrs[i].ptr + start_index, (size_t)(capacity - start_index));
			memcpy((char*)new_data + (capacity - start_index) + running_offset, (char*)buffer->array_ptrs[i].ptr, (size_t)(start_index));
		}
		else {
			memcpy((char*)new_data + running_offset, (char*)buffer->array_ptrs[i].ptr + start_index, (size_t)(capacity - start_index));
		}
		running_offset += buffer->array_ptrs[i].unit_size * new_capacity;
		running_offset += buffer->alignment - (running_offset % buffer->alignment);

	}
	void *old_data = buffer->data;

	buffer->data = new_data;
	buffer->capacity = new_capacity;
	buffer->current_arena_size = new_capacity * buffer->struct_size;
	running_offset = 0;
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		buffer->array_ptrs[i].ptr = (char*)buffer->data + running_offset;
		memcpy((char*)buffer->struct_of_arrays + buffer->array_ptrs[i].offset, &buffer->array_ptrs[i].ptr, sizeof(void*));
		running_offset += buffer->array_ptrs[i].unit_size * buffer->capacity;
		running_offset += buffer->alignment - (running_offset % buffer->alignment);
	}
	tfxFREE(old_data);
}

#define tmpStack(type, name) tfx_vector_t<type> name

template <typename T>
struct tfx_bucket_t {
	tfx_vector_t<T> data;
	tfx_bucket_t *next;
};

template <typename T>
inline tfx_bucket_t<T> *tfxCreateBucket(tfxU32 size) {
	tfx_bucket_t<T> *bucket = (tfx_bucket_t<T>*)tfxALLOCATE(sizeof(tfx_bucket_t<T>));
	bucket->data.data = nullptr;
	bucket->data.current_size = 0;
	bucket->data.capacity = 0;
	bucket->data.locked = 0;
	bucket->data.alignment = 0;
	bucket->data.reserve(size);
	bucket->next = nullptr;
	return bucket;
}

template <typename T>
struct tfx_bucket_array_t {
	tfxU32 current_size;
	tfxU32 capacity;
	tfxU32 size_of_each_bucket;
	tfxLONG volatile locked;
	tfx_vector_t<tfx_bucket_t<T>*> bucket_list;

	inline bool			empty() { return current_size == 0; }
	inline tfxU32		size() { return current_size; }
	inline void			free_all() {
		for (tfx_bucket_t<T> *bucket : bucket_list) {
			bucket->data.free();
			tfxFREE(bucket);
		}
		bucket_list.free();
	}
	inline T&           operator[](tfxU32 i) {
		assert(i < current_size);		//Index is out of bounds
		tfxU32 bucket_index = i / size_of_each_bucket;
		tfxU32 element_index = i % size_of_each_bucket;
		return (*bucket_list[bucket_index]).data[element_index];
	}
	inline const T&     operator[](tfxU32 i) const {
		assert(i < current_size);		//Index is out of bounds
		tfxU32 bucket_index = i / size_of_each_bucket;
		tfxU32 element_index = i % size_of_each_bucket;
		return (*bucket_list[bucket_index]).data[element_index];
	}
	inline T*           begin() { return bucket_list.current_size ? (T*)bucket_list[0]->data.data : nullptr; }
	inline const T*     begin() const { return bucket_list.current_size ? (T*)bucket_list[0]->data.data : nullptr; }
	inline T*           end() { return bucket_list.current_size ? (T*)bucket_list[(current_size - 1) / size_of_each_bucket]->data.end() : nullptr; }
	inline const T*     end() const { return bucket_list.current_size ? (T*)bucket_list[(current_size - 1) / size_of_each_bucket]->data.end() : nullptr; }
	inline T&           front() { assert(current_size > 0); return bucket_list[0]->data.front(); }
	inline const T&     front() const { assert(current_size > 0); return bucket_list[0]->data.front(); }
	inline T&           back() { assert(current_size > 0); return bucket_list[(current_size - 1) / size_of_each_bucket]->data.back(); }
	inline const T&     back() const { assert(current_size > 0); return bucket_list[(current_size - 1) / size_of_each_bucket]->data.back(); }
	inline tfxU32		active_buckets() { return current_size == 0 ? 0 : current_size / size_of_each_bucket + 1; }
	inline void         clear() {
		for (tfx_bucket_t<T> *bucket : bucket_list) {
			bucket->data.clear();
		}
		current_size = 0;
	}

	inline tfx_bucket_t<T> *add_bucket(tfxU32 size_of_each_bucket) {
		//current_bucket must be the last bucket in the chain
		tfx_bucket_t<T> *bucket = tfxCreateBucket<T>(size_of_each_bucket);
		bucket_list.push_back(bucket);
		return bucket;
	}

	inline T& push_back(const T& v) {
		if (current_size == capacity) {
			add_bucket(size_of_each_bucket);
			capacity += size_of_each_bucket;
		}
		tfxU32 current_bucket = current_size / size_of_each_bucket;
		bucket_list[current_bucket]->data.push_back(v);
		current_size++;
		return bucket_list[current_bucket]->data.back();
	}

	inline tfxU32        locked_push_back(const T& v) {
		while (tfx__compare_and_exchange(&locked, 1, 0) > 1);

		push_back(v);

		tfx__exchange(&locked, 0);
		return current_size - 1;
	}

	inline T*	insert(tfxU32 insert_index, const T &v) {
		assert(insert_index < current_size);
		tfxU32 insert_bucket = insert_index / size_of_each_bucket;
		tfxU32 element_index = insert_index % size_of_each_bucket;
		if (bucket_list[insert_bucket]->data.current_size < bucket_list[insert_bucket]->data.capacity) {
			//We're inserting in the last bucket
			return bucket_list[insert_bucket]->data.insert(&bucket_list[insert_bucket]->data[element_index], v);
		}
		T end_element = bucket_list[insert_bucket]->data.pop_back();
		T end_element2;
		bool end_pushed = false;
		bool end2_pushed = true;
		bucket_list[insert_bucket]->data.insert(&bucket_list[insert_bucket]->data[element_index], v);
		tfxU32 current_insert_bucket = insert_bucket;
		tfxU32 alternator = 0;
		while (current_insert_bucket++ < active_buckets() - 1) {
			if (bucket_list[current_insert_bucket]->data.full()) {
				if (alternator == 0) {
					end_element2 = bucket_list[current_insert_bucket]->data.pop_back();
					end2_pushed = false;
					bucket_list[current_insert_bucket]->data.push_front(end_element);
					end_pushed = true;
				}
				else {
					end_element = bucket_list[current_insert_bucket]->data.pop_back();
					end_pushed = false;
					bucket_list[current_insert_bucket]->data.push_front(end_element2);
					end2_pushed = true;
				}
				alternator = !alternator;
			}
			else {
				bucket_list[current_insert_bucket]->data.push_front(alternator == 0 ? end_element : end_element2);
				end_pushed = true;
				end2_pushed = true;
				break;
			}
		}
		if (!end_pushed || !end2_pushed) {
			push_back(!end_pushed ? end_element : end_element2);
		}
		else {
			current_size++;
		}
		return &bucket_list[insert_bucket]->data[element_index];
	}

	inline T*	insert(T* position, const T &v) {
		tfxU32 index = 0;
		bool find_result = find(position, index);
		assert(find_result);	//Could not find the object to insert at, make sure it exists
		return insert(index, v);
	}

	inline void erase(tfxU32 index) {
		assert(index < current_size);
		tfxU32 bucket_index = index / size_of_each_bucket;
		tfxU32 element_index = index % size_of_each_bucket;
		bucket_list[bucket_index]->data.erase(&bucket_list[bucket_index]->data[element_index]);
		current_size--;
		if (bucket_index == bucket_list.current_size - 1) {
			//We're erasing in the last bucket
			return;
		}
		tfxU32 current_bucket_index = bucket_index;
		while (current_bucket_index < active_buckets() - 1) {
			T front = bucket_list[current_bucket_index + 1]->data.pop_front();
			bucket_list[current_bucket_index]->data.push_back(front);
			current_bucket_index++;
		}
		trim_buckets();
	}

	inline void erase(T* it) {
		tfxU32 index = 0;
		bool find_result = find(it, index);
		assert(find_result);	//pointer not found in list
		erase(index);
	}

	inline bool find(T *it, tfxU32 &index) {
		for (int i = 0; i != current_size; ++i) {
			if (it == &(*this)[i]) {
				index = i;
				return true;
			}
		}
		return false;
	}

	inline T* find(T *it) {
		for (int i = 0; i != current_size; ++i) {
			if (*it == (*this)[i]) {
				return &(*this)[i];
			}
		}
		return nullptr;
	}

	inline void trim_buckets() {
		if (active_buckets() < bucket_list.current_size) {
			for (int i = active_buckets(); i != bucket_list.current_size; ++i) {
				bucket_list[i]->data.free();
				tfxFREE(bucket_list[i]);
				capacity -= size_of_each_bucket;
			}
			bucket_list.current_size -= bucket_list.current_size - active_buckets();
		}
	}

};

template <typename T>
inline tfx_bucket_array_t<T> tfxCreateBucketArray(tfxU32 size_of_each_bucket) {
	tfx_bucket_array_t<T> bucket_array{};
	bucket_array.current_size = bucket_array.locked = bucket_array.capacity = 0;
	bucket_array.size_of_each_bucket = size_of_each_bucket;
	return bucket_array;
}

template <typename T>
inline void tfxCopyBucketArray(tfx_bucket_array_t<T> *dst, tfx_bucket_array_t<T> *src) {
	if (src == dst) {
		return;
	}
	dst->free_all();
	for (tfx_bucket_t<T> *bucket : src->bucket_list) {
		tfx_bucket_t<T> *copy = dst->add_bucket(src->size_of_each_bucket);
		copy->data = bucket->data;
	}
	dst->current_size = src->current_size;
	dst->capacity = src->capacity;
	dst->size_of_each_bucket = src->size_of_each_bucket;
}

#define tfxBucketLoop(bucket) int i = 0; i != bucket.current_size; ++i
#define tfxBucketLoopNI(bucket) int ni = 0; ni != bucket.current_size; ++ni

//A char buffer you can use to load a file into and read from
//Has no deconstructor so make sure you call FreeAll() when done
//This is meant for limited usage in timeline fx only and not recommended for use outside!
struct tfx_stream_t {
	tfxU64 size = 0;
	tfxU64 position = 0;
	char* data = nullptr;

	inline tfx_stream_t() { size = position = 0; data = nullptr; }
	inline tfx_stream_t(tfxU64 qty) { size = position = 0; data = nullptr; Resize(qty); }

	inline bool Read(char* dst, tfxU64 count) {
		if (count + position <= size) {
			memcpy(dst, data + position, count);
			position += count;
			return true;
		}
		return false;
	}
	inline tfx_str512_t ReadLine();
	inline bool Write(void *src, tfxU64 count) {
		if (count + position <= size) {
			memcpy(data + position, src, count);
			position += count;
			return true;
		}
		return false;
	}
	inline bool EoF() { return position >= size; }
	inline void Seek(tfxU64 offset) {
		if (offset < size)
			position = offset;
		else
			position = size;
	}

	inline bool			Empty() { return size == 0; }
	inline tfxU64		Size() { return size; }
	inline const tfxU64	Size() const { return size; }

	inline void			FreeAll() { if (data) { size = size = 0; tfxFREE(data); data = nullptr; } }
	inline void         Clear() { if (data) { size = 0; } }

	inline void         Resize(tfxU64 new_capacity) {
		if (new_capacity <= size)
			return;
		char* new_data = (char*)tfxALLOCATE((tfxU64)new_capacity * sizeof(char));
		assert(new_data);	//Unable to allocate memory. Todo: better handling
		if (data) {
			memcpy(new_data, data, (tfxU64)size * sizeof(char));
			tfxFREE(data);
		}
		data = new_data; size = new_capacity;
		position = 0;
	}
	inline void			NullTerminate() { *(data + size) = '\0'; }

};

//-----------------------------------------------------------
//Section: Multithreading_Work_Queues
//-----------------------------------------------------------

//Tried to keep this as simple as possible, was originally based on Casey Muratory's Hand Made Hero threading which used the Windows API for
//threading but for the sake of supporting other platforms I changed it to use std::thread which was actually a lot more simple to do then 
//I expected. I just had to swap the semaphores for condition_variable and that was pretty much it other then obviously using std::thread as well.
//There is a single thread pool created to serve multiple queues. Currently each particle manager that you create will have it's own queue.
struct tfx_work_queue_t;

#define tfxWORKQUEUECALLBACK(name) void name(tfx_work_queue_t *queue, void *data)
typedef tfxWORKQUEUECALLBACK(tfxWorkQueueCallback);

struct tfx_work_queue_entry_t {
	tfxWorkQueueCallback *call_back = nullptr;
	void *data = nullptr;
};

typedef tfxU32 tfxWorkQueueFlags;

enum tfxWorkQueueFlag_ {
	tfxWorkQueueFlag_none = 0
};

extern int tfxNumberOfThreadsInAdditionToMain;
extern bool tfxThreadUsage[];		//Used for debugging to see which threads were utilised each frame

#ifndef tfxMAX_QUEUES
#define tfxMAX_QUEUES 512
#endif

struct tfx_work_queue_t {
	tfxU32 volatile entry_completion_goal = 0;
	tfxU32 volatile entry_completion_count = 0;
	tfxLONG volatile next_read_entry = 0;
	tfxLONG volatile next_write_entry = 0;
	tfx_work_queue_entry_t entries[tfxMAX_QUEUES];
};

struct tfx_queue_processor_t {
	std::mutex mutex;
	std::condition_variable empty_condition;
	std::condition_variable full_condition;
	tfxU32 count;
	//These point to the queues stored in particle managers and anything else that needs a queue with multithreading
	tfx_work_queue_t *queues[tfxMAX_QUEUES];
};

extern tfx_queue_processor_t tfxThreadQueues;

tfxINTERNAL inline void InitialiseThreadQueues(tfx_queue_processor_t *queues) {
	queues->count = 0;
	memset(queues->queues, 0, tfxMAX_QUEUES * sizeof(void*));
}

tfxINTERNAL inline tfx_work_queue_t *tfxGetQueueWithWork(tfx_queue_processor_t *thread_processor) {
	std::unique_lock<std::mutex> lock(thread_processor->mutex);
	thread_processor->full_condition.wait(lock, [&]() { return thread_processor->count > 0; });
	tfx_work_queue_t *queue = thread_processor->queues[--thread_processor->count];
	thread_processor->empty_condition.notify_one();
	return queue;
}

tfxINTERNAL inline void tfxPushQueueWork(tfx_queue_processor_t *thread_processor, tfx_work_queue_t *queue) {
	std::unique_lock<std::mutex> lock(thread_processor->mutex);
	thread_processor->empty_condition.wait(lock, [&]() { return thread_processor->count < tfxMAX_QUEUES; });
	thread_processor->queues[thread_processor->count++] = queue;
	thread_processor->full_condition.notify_one();
}

tfxINTERNAL inline void tfxDoNextWorkQueue(tfx_queue_processor_t *queue_processor) {
	tfx_work_queue_t *queue = tfxGetQueueWithWork(queue_processor);

	if (queue) {
		tfxU32 original_read_entry = queue->next_read_entry;
		tfxU32 new_original_read_entry = (original_read_entry + 1) % tfxMAX_QUEUES;

		if (original_read_entry != queue->next_write_entry) {
			tfxU32 index = tfx__compare_and_exchange(&queue->next_read_entry, new_original_read_entry, original_read_entry);
			if (index == original_read_entry) {
				tfx_work_queue_entry_t entry = queue->entries[index];
				entry.call_back(queue, entry.data);
				tfx__increment(&queue->entry_completion_count);
			}
		}
	}
}

tfxINTERNAL inline void tfxDoNextWorkQueueEntry(tfx_work_queue_t *queue) {
	tfxU32 original_read_entry = queue->next_read_entry;
	tfxU32 new_original_read_entry = (original_read_entry + 1) % tfxMAX_QUEUES;

	if (original_read_entry != queue->next_write_entry) {
		tfxU32 index = tfx__compare_and_exchange(&queue->next_read_entry, new_original_read_entry, original_read_entry);
		if (index == original_read_entry) {
			tfx_work_queue_entry_t entry = queue->entries[index];
			entry.call_back(queue, entry.data);
			tfx__increment(&queue->entry_completion_count);
		}
	}
}

tfxINTERNAL inline void tfxAddWorkQueueEntry(tfx_work_queue_t *queue, void *data, tfxWorkQueueCallback call_back) {
	if (!tfxNumberOfThreadsInAdditionToMain) {
		call_back(queue, data);
		return;
	}

	tfxU32 new_entry_to_write = (queue->next_write_entry + 1) % tfxMAX_QUEUES;
	while (new_entry_to_write == queue->next_read_entry) {		//Not enough room in work queue
		//We can do this because we're single producer
		tfxDoNextWorkQueueEntry(queue);
	}
	queue->entries[queue->next_write_entry].data = data;
	queue->entries[queue->next_write_entry].call_back = call_back;
	tfx__increment(&queue->entry_completion_goal);

	tfx__writebarrier;

	tfxPushQueueWork(&tfxThreadQueues, queue);
	queue->next_write_entry = new_entry_to_write;

}

tfxINTERNAL inline void tfxCompleteAllWork(tfx_work_queue_t *queue) {
	tfx_work_queue_entry_t entry = {};
	while (queue->entry_completion_goal != queue->entry_completion_count) {
		tfxDoNextWorkQueueEntry(queue);
	}
	queue->entry_completion_count = 0;
	queue->entry_completion_goal = 0;
}

tfxINTERNAL inline void tfxInitialiseWorkQueue(tfx_work_queue_t *queue) {
	queue->entry_completion_count = 0;
	queue->entry_completion_goal = 0;
	queue->next_read_entry = 0;
	queue->next_write_entry = 0;
}

tfxINTERNAL inline bool tfxInitialiseThreads(tfx_queue_processor_t *thread_queues) {
	InitialiseThreadQueues(&tfxThreadQueues);

	//todo: create a function to close all the threads 

	tfxU32 threads_initialised = 0;
	for (int thread_index = 0; thread_index < tfxNumberOfThreadsInAdditionToMain; ++thread_index) {
		std::thread([thread_queues]() {
			for (;;) {
				tfxDoNextWorkQueue(thread_queues);
			}
			}).detach();

			threads_initialised++;
	}
	return true;
}

//-----------------------------------------------------------
//Section: Global_Variables
//-----------------------------------------------------------
extern const tfxU32 tfxPROFILE_COUNT;

struct tfx_data_types_dictionary_t {
	bool initialised = false;
	tfx_storage_map_t<tfx_data_type> names_and_types;
	tfx_data_types_dictionary_t() :
		names_and_types("Data Types Storage Map", "Data Types Storage Data")
	{}
	void Init();
};

//Global variables
struct tfx_storage_t {
	tfxU32 memory_pool_count;
	size_t default_memory_pool_size;
	size_t memory_pool_sizes[tfxMAX_MEMORY_POOLS];
	tfx_pool *memory_pools[tfxMAX_MEMORY_POOLS];
	tfx_data_types_dictionary_t data_types;
	//tfx_storage_map_t<tfx_particle_manager_t*> particle_managers;
	//tfx_storage_map_t<tfx_animation_manager_t*> animation_managers;
};

extern tfx_storage_t *tfxStore;
extern tfx_allocator *tfxMemoryAllocator;

//-----------------------------------------------------------
//Section: Vector_Math
//-----------------------------------------------------------

//Just the very basic vector types that we need
struct tfx_vec2_t {
	float x, y;

	tfx_vec2_t() { x = y = 0.f; }
	tfx_vec2_t(float _x, float _y) : x(_x), y(_y) {}

	inline void operator=(float v) { x = v; y = v; }

	inline tfx_vec2_t operator+(tfx_vec2_t v) { return tfx_vec2_t(x + v.x, y + v.y); }
	inline void operator+=(tfx_vec2_t v) { x += v.x; y += v.y; }
	inline tfx_vec2_t operator-(tfx_vec2_t v) { return tfx_vec2_t(x - v.x, y - v.y); }
	inline tfx_vec2_t operator-() { return tfx_vec2_t(-x, -y); }
	inline void operator-=(tfx_vec2_t v) { x -= v.x; y -= v.y; }
	inline tfx_vec2_t operator*(tfx_vec2_t v) { return tfx_vec2_t(x * v.x, y * v.y); }
	inline void operator*=(tfx_vec2_t v) { x *= v.x; y *= v.y; }
	inline tfx_vec2_t operator/(tfx_vec2_t v) { return tfx_vec2_t(x / v.x, y / v.y); }
	inline void operator/=(tfx_vec2_t v) { x /= v.x; y /= v.y; }
	inline tfx_vec2_t operator+(float v) { return tfx_vec2_t(x + v, y + v); }
	inline tfx_vec2_t operator-(float v) { return tfx_vec2_t(x - v, y - v); }
	inline tfx_vec2_t operator*(float v) { return tfx_vec2_t(x * v, y * v); }
	inline void operator*=(float v) { x *= v; y *= v; }
	inline bool operator>(tfx_vec2_t &v) { return x + y > v.x + v.y; }
	inline tfx_vec2_t operator/(float v) { return tfx_vec2_t(x / v, y / v); }

	inline tfx_vec2_t operator+(float v) const { return tfx_vec2_t(x + v, y + v); }
	inline tfx_vec2_t operator-(float v) const { return tfx_vec2_t(x - v, y - v); }
	inline tfx_vec2_t operator*(float v) const { return tfx_vec2_t(x * v, y * v); }
	inline tfx_vec2_t operator/(float v) const { return tfx_vec2_t(x / v, y / v); }

	inline tfx_vec2_t operator+(const tfx_vec2_t &v) const { return tfx_vec2_t(x + v.x, y + v.y); }
	inline tfx_vec2_t operator-(const tfx_vec2_t &v) const { return tfx_vec2_t(x - v.x, y - v.y); }
	inline tfx_vec2_t operator*(const tfx_vec2_t &v) const { return tfx_vec2_t(x * v.x, y * v.y); }
	inline tfx_vec2_t operator/(const tfx_vec2_t &v) const { return tfx_vec2_t(x / v.x, y / v.y); }

	inline float Squared() { return x * x + y * y; }
	inline bool IsNill() { return !x && !y; }
};
inline tfx_vec2_t operator*(float ls, tfx_vec2_t rs) { return tfx_vec2_t(rs.x * ls, rs.y * ls); }

typedef struct tfx_float16x4_s {
	union {
		struct {
			tfxU16 x : 16;
			tfxU16 y : 16;
			tfxU16 z : 16;
			tfxU16 w : 16;
		};
		struct {
			tfxU32 xy : 32;
			tfxU32 zw : 32;
		};
		struct { tfxU64 packed; };
	};
} tfx_float16x4_t;

typedef struct tfx_float16x2_s {
	union {
		struct {
			tfxU16 x : 16;
			tfxU16 y : 16;
		};
		struct { tfxU32 packed; };
	};
} tfx_float16x2_t;

typedef struct tfx_float8x4_s {
	union {
		struct {
			tfxU8 x : 8;
			tfxU8 y : 8;
			tfxU8 z : 8;
			tfxU8 w : 8;
		};
		struct { tfxU32 packed; };
	};
} tfx_float8x4_t;

typedef struct tfx_float32x2_s {
	float x, y;
}tfx_float32x2_t;

typedef struct tfx_float32x3_s {
	float x, y, z;
}tfx_float32x3_t;

typedef struct tfx_float32x4_s {
	float x, y, z, w;
}tfx_float32x4_t;

typedef struct tfx_instance_s {		//48 bytes
	tfx_float32x4_t position;						//The position of the billboard with stretch in w
	tfxU64 quaternion;								//Rotation of the billboard stored as a 16-bit snorm quaternion
	tfx_float16x2_t size;							//Size/Scale of the sprite
	tfx_float8x4_t alignment;						//normalised alignment vector 3 8bit floats packed into 32 bits. Free byte here.
	tfx_float16x2_t intensity_gradient_map;			//Multiplier for the color and life of particle
	tfx_float8x4_t curved_alpha_life;				//Sharpness and dissolve amount value for fading the image plus the age of the particle value packed into 3 bit unorms. Free byte here.
	tfxU32 indexes;									//[gpu property index, capture flag (1 bit << 15), image data index max 8191 images]
	tfxU32 captured_index;							//Index to the sprite in the buffer from the previous frame for interpolation
} tfx_instance_t;

typedef struct tfx_ribbon_vertex_s {
	tfx_float32x3_t position;
	tfxU32 segment_index;
	tfx_float32x2_t uv_offset_scale;
	tfxU32 ribbon_index;
	tfxU32 clipped;
} tfx_ribbon_vertex_t;

typedef struct tfx_gpu_particle_properties_s {
	tfx_float32x2_t image_handle;			
	tfxU32 color_ramp_indexes;			//[Row of color ramp bitmap, texture array]
	tfxU32 flags;						//Flags like billboard alignment type
	tfxU32 start_frame_index;
	float animation_frames;
	int padding[2];
} tfx_gpu_particle_properties_t;

typedef struct tfx_gpu_image_data_s {
	tfx_float32x4_t uv;
	tfxU64 uv_packed;
	tfx_float32x2_t image_size;
	tfxU32 texture_array_index;
	float animation_frames;
	float max_radius;
}tfx_gpu_image_data_t;

//------------------------------------------------------------
//Section: Enums
//------------------------------------------------------------

typedef enum tfx_color_format {
	tfx_color_format_rgba8,
	tfx_color_format_rgba16f,
	tfx_color_format_rgba32f,
} tfx_color_format;

typedef enum {
	tfxStageSetup_none,
	tfxStageSetup_group_sprites_by_effect,
} tfx_stage_setup;

typedef enum {
	tfxEffect_global_life_index,
	tfxEffect_global_amount_index,
	tfxEffect_global_velocity_index,
	tfxEffect_global_noise_index,
	tfxEffect_global_width_index,
	tfxEffect_global_height_index,
	tfxEffect_global_weight_index,
	tfxEffect_global_roll_spin_index,
	tfxEffect_global_pitch_spin_index,
	tfxEffect_global_yaw_spin_index,
	tfxEffect_global_stretch_index,
	tfxEffect_global_overall_scale_index,
	tfxEffect_global_intensity_index,
	tfxEffect_global_splatter_index,
	tfxEffect_global_emitter_width_index,
	tfxEffect_global_emitter_height_index,
	tfxEffect_global_emitter_depth_index,
	tfxEffectGraphs_max_index,
} tfx_global_graph_index;

typedef enum {
	tfxEmitter_property_emission_pitch_index,
	tfxEmitter_property_emission_yaw_index,
	tfxEmitter_property_emission_range_index,
	tfxEmitter_property_splatter_index,
	tfxEmitter_property_width_index,        //Also used for linear extrusion for paths as well
	tfxEmitter_property_height_index,
	tfxEmitter_property_depth_index,
	tfxEmitter_property_extrusion_index,
	tfxEmitter_property_arc_size_index,
	tfxEmitter_property_arc_offset_index,

	tfxEmitter_base_life_index,
	tfxEmitter_base_amount_index,
	tfxEmitter_base_velocity_index,
	tfxEmitter_base_width_index,
	tfxEmitter_base_height_index,
	tfxEmitter_base_weight_index,
	tfxEmitter_base_pitch_spin_index,
	tfxEmitter_base_yaw_spin_index,
	tfxEmitter_base_roll_spin_index,

	tfxEmitter_variation_life_index,
	tfxEmitter_variation_amount_index,
	tfxEmitter_variation_velocity_index,
	tfxEmitter_variation_width_index,
	tfxEmitter_variation_height_index,
	tfxEmitter_variation_weight_index,
	tfxEmitter_variation_path_trajectory_scale_index,
	tfxEmitter_variation_pitch_spin_index,
	tfxEmitter_variation_yaw_spin_index,
	tfxEmitter_variation_roll_spin_index,
	tfxEmitter_variation_noise_resolution_index,
	tfxEmitter_variation_motion_randomness_index,

	tfxEmitter_overtime_red_index,
	tfxEmitter_overtime_green_index,
	tfxEmitter_overtime_blue_index,
	tfxEmitter_overtime_blendfactor_index,
	tfxEmitter_overtime_velocity_adjuster_index,
	tfxEmitter_overtime_intensity_index,
	tfxEmitter_overtime_alpha_sharpness_index,
	tfxEmitter_overtime_curved_alpha_index,
	tfxEmitter_overtime_heat_response_index,
	tfxEmitter_overtime_gradient_mapper_index,
	tfxEmitter_overtime_velocity_index,
	tfxEmitter_overtime_width_index,
	tfxEmitter_overtime_height_index,
	tfxEmitter_overtime_weight_index,
	tfxEmitter_overtime_pitch_spin_index,
	tfxEmitter_overtime_yaw_spin_index,
	tfxEmitter_overtime_roll_spin_index,
	tfxEmitter_overtime_stretch_index,
	tfxEmitter_overtime_velocity_turbulance_index,
	tfxEmitter_overtime_direction_turbulance_index,
	tfxEmitter_overtime_direction_index,
	tfxEmitter_overtime_noise_resolution_index,
	tfxEmitter_overtime_motion_randomness_index,

	tfxEmitter_factor_life_index,
	tfxEmitter_factor_size_index,
	tfxEmitter_factor_velocity_index,
	tfxEmitter_factor_intensity_index,

	tfxEmitterGraphs_max_index,
} tfx_emitter_graph_index;

typedef enum {
	tfxRibbon_property_splatter_index,
	tfxRibbon_property_width_index,        //Also used for linear extrusion for paths as well
	tfxRibbon_property_height_index,
	tfxRibbon_property_depth_index,
	tfxRibbon_property_extrusion_index,
	tfxRibbon_property_arc_size_index,
	tfxRibbon_property_arc_offset_index,

	tfxRibbon_base_life_index,
	tfxRibbon_base_amount_index,
	tfxRibbon_base_width_index,

	tfxRibbon_variation_life_index,
	tfxRibbon_variation_amount_index,
	tfxRibbon_variation_width_index,

	tfxRibbon_overtime_red_index,
	tfxRibbon_overtime_green_index,
	tfxRibbon_overtime_blue_index,
	tfxRibbon_overtime_blendfactor_index,
	tfxRibbon_overtime_intensity_index,
	tfxRibbon_overtime_alpha_sharpness_index,
	tfxRibbon_overtime_curved_alpha_index,
	tfxRibbon_overtime_gradient_mapper_index,
	tfxRibbon_overtime_heat_response_index,
	tfxRibbon_overtime_width_index,
	tfxRibbon_overtime_overall_scale_index,
	tfxRibbon_overtime_uv_offset_y_index,
	tfxRibbon_overtime_uv_scale_y_index,
	tfxRibbon_overtime_clip_start_index,
	tfxRibbon_overtime_clip_end_index,

	tfxRibbon_overlength_intensity_index,
	tfxRibbon_overlength_alpha_sharpness_index,
	tfxRibbon_overlength_curved_alpha_index,
	tfxRibbon_overlength_gradient_map_index,
	tfxRibbon_overlength_width_index,
	tfxRibbon_overlength_fixed_angle_index,

	tfxRibbonGraphs_max_index,

	tfxRibbon_property_start_index = 0,
	tfxRibbon_base_start_index = tfxRibbon_base_life_index,
	tfxRibbon_variation_start_index = tfxRibbon_variation_life_index,
	tfxRibbon_overtime_start_index = tfxRibbon_overtime_red_index,
	tfxRibbon_property_end_index = tfxRibbon_property_arc_offset_index + 1,
	tfxRibbon_base_end_index = tfxRibbon_base_width_index + 1,
	tfxRibbon_variation_end_index = tfxRibbon_variation_width_index + 1,
	tfxRibbon_overtime_end_index = tfxRibbon_overtime_clip_end_index + 1,
	tfxRibbon_overlength_start = tfxRibbon_overlength_intensity_index,
	tfxRibbon_overlength_end = tfxRibbon_overlength_fixed_angle_index + 1,
} tfx_ribbon_graph_index;

//tfx_effect_descriptor_t type - effect contains emitters, and emitters spawn particles, but they both share the same struct for simplicity
typedef enum {
	tfxEffectType,
	tfxEmitterType,
	tfxRibbonType,
	tfxFolder,
	//Not a descriptor type an effect can have - it only ever tags a tfx_graph_list_t so that the list knows
	//how many graphs it holds and which initialiser rebuilds it. Appended rather than inserted because the
	//preceding values are saved in the file as ordinals.
	tfxForceType,
	tfxMaxDescriptorTypes
} tfx_effect_descriptor_type;

typedef enum {
	tfxColorInterpolation_linear_srgb = 0,
	tfxColorInterpolation_oklch,
	tfxColorInterpolation_hsl,
	tfxColorInterpolation_linear_rgb,
	tfxColorInterpolation_max
} tfx_color_interpolation_mode;

typedef enum {
	tfxGraphEasingType_constant                                 = 0,
	tfxGraphEasingType_smoothstep                               = 17,
	tfxGraphEasingType_out_in									= 18,
	tfxGraphEasingType_in										= 4,
	tfxGraphEasingType_out										= 5,
	tfxGraphEasingType_in_out									= 6,
	tfxGraphEasingType_linear                                   = 16,
	//Unused
	tfxGraphEasingType_ease_in_quad                             = 1,
	tfxGraphEasingType_ease_out_quad                            = 2,
	tfxGraphEasingType_ease_in_out_quad                         = 3,
	tfxGraphEasingType_ease_in_circular                         = 13,
	tfxGraphEasingType_ease_out_circular                        = 14,
	tfxGraphEasingType_ease_in_out_circular                     = 15,
} tfx_graph_easing_type;

typedef enum {
	tfxErrorCode_success                                        = 0,
	tfxErrorCode_incorrect_package_format                       = 1 << 0,
	tfxErrorCode_data_could_not_be_loaded                       = 1 << 1,
	tfxErrorCode_could_not_add_shape                            = 1 << 2,
	tfxErrorCode_error_loading_shapes                           = 1 << 3,
	tfxErrorCode_some_data_not_loaded                           = 1 << 4,
	tfxErrorCode_unable_to_open_file                            = 1 << 5,
	tfxErrorCode_unable_to_read_file                            = 1 << 6,
	tfxErrorCode_wrong_file_size                                = 1 << 7,
	tfxErrorCode_invalid_format                                 = 1 << 8,
	tfxErrorCode_no_inventory                                   = 1 << 9,
	tfxErrorCode_invalid_inventory                              = 1 << 10,
	tfxErrorCode_file_version_out_of_date                       = 1 << 11,
	tfxErrorCode_library_loaded_without_shape_loader            = 1 << 13,
	tfxErrorCode_library_object_could_not_be_created            = 1 << 14
} tfx_error_flag_bits;

typedef tfxU32 tfxErrorFlags;                   //tfx_error_flag_bits

//-----------------------------------------------------------
//Section: forward_declarations
//-----------------------------------------------------------
#define tfxMAKE_HANDLE(handle) typedef struct handle##_s* handle;

//For allocating a new object with handle. Only used internally.
#define tfxNEW(type) (type)tfxALLOCATE(sizeof(type##_t))
#define tfxNEW_ALIGNED(type, alignment) (type)tfxALLOCATE_ALIGNED(sizeof(type##_t), alignment)

typedef struct tfx_package_s tfx_package_t;

tfxMAKE_HANDLE(tfx_package)
tfxMAKE_HANDLE(tfx_library);
tfxMAKE_HANDLE(tfx_stage);
tfxMAKE_HANDLE(tfx_animation_manager);
tfxMAKE_HANDLE(tfx_effect_descriptor);
tfxMAKE_HANDLE(tfx_effect_template);
tfxMAKE_HANDLE(tfx_ribbon_buffer_requirements);
tfxMAKE_HANDLE(tfx_ribbon_dispatch);
tfxMAKE_HANDLE(tfx_gpu_shapes);

typedef struct tfx_image_data_s tfx_image_data_t;
typedef struct tfx_gpu_image_data_s tfx_gpu_image_data_t;
typedef struct tfx_bitmap_s tfx_bitmap_t;
typedef struct tfx_gpu_graph_data_s tfx_gpu_graph_data_t;
typedef struct tfx_instance_s tfx_instance_t;
typedef struct tfx_ribbon_bucket_s tfx_ribbon_bucket_t;
typedef struct tfx_sprite_data_s tfx_sprite_data_t;
typedef struct tfx_effect_index_s tfx_effect_index_t;
typedef struct tfx_effect_instance_data_s tfx_effect_instance_data_t;
typedef struct tfx_animation_instance_s tfx_animation_instance_t;
typedef struct tfx_frame_meta_s tfx_frame_meta_t;
typedef struct tfx_sprite_data_settings_s tfx_sprite_data_settings_t;
typedef struct tfx_emitter_path_s tfx_emitter_path_t;

typedef struct tfx_random_s {
	tfxU64 seeds[2];
}tfx_random_t;

typedef struct tfx_version_s {
	const char *name;
	int major;
	int minor;
	int patch;
} tfx_version_t;

typedef struct tfx_ribbon_buffer_info_s {
	tfxU32 vertices_per_segment; 
	tfxU32 triangles_per_segment; 
	tfxU32 indices_per_segment;  
	tfxU32 total_segments; 
	tfxU32 index_count; 
	tfxKey pipeline_index;
} tfx_ribbon_buffer_info_t;

//This struct is used for configuring a effect manager on creation
typedef struct tfx_stage_info_s {
	double warmup_delta_time;				//The frame length tick amount for warming up effects. Higher is more performant at the cost of accuracy.
	tfxU32 max_particles;					//The maximum number of instance_data for each layer. This setting is not relevent if dynamic_sprite_allocation is set to true or group_sprites_by_effect is true.
	tfxU32 max_effects;                     //The maximum number of effects that can be updated at the same time.
	tfxU32 max_ribbon_segments;             //All segments for ribbons are stored in a single buffer. You will need to create buffers for rendering and so whatever you decide the max segments should be your buffers
											//should be big enough to contain all ribbon segments that you might need. You can call tfx_GetSegmentBufferSizeInBytes after creating the effect manager to get the byte
											//value that you can use to create the buffers. Also note that segments are always created in multiples of 32, so whatever number you put here it will be rounded to the
											//nearest multiple of 32.
	tfxU32 max_ribbons;						//The maximum number of ribbon instances that can exist at the same time across all emitters in this effect manager.
	tfxU32 ribbon_tessellation;				//The amount of tessellation used for ribbons. Currently this is set globally. 1 is generally enough for most cases.
	tfxU32 multi_threaded_batch_size;       //The size of each batch of particles to be processed when multithreading. Must be a power of 2 and 256 or greater.
	tfxU32 sort_passes;                     //when in order by depth mode (not guaranteed order) set the number of sort passes for more accuracy. Anything above 5 and you should just be guaranteed order.
	bool double_buffer_sprites;             //Set to true to double buffer instance_data so that you can interpolate between the old and new positions for smoother animations.
	bool dynamic_sprite_allocation;         //Set to true to automatically resize the sprite buffers if they run out of space. Not applicable when grouping instance_data by effect.
	bool group_sprites_by_effect;           //Set to true to group all instance_data by effect. Effects can then be drawn in specific orders or not drawn at all on an effect by effect basis.
	bool auto_order_effects;                //When group_sprites_by_effect is true then you can set this to true to sort the effects each frame. Use tfx_SetStageCamera in 3d to set the effect depth to the distance the camera.
	void *user_data;						//User data that will get passed into the grow_staging_buffer_callback function which you can use to grow the buffer
	//If you need the staging buffer to be grown dynamically then you can use this call back to do that. It should return true if the buffer was successfully grown or false otherwise.
	bool(*grow_staging_buffer_callback)(tfxU32 new_size, tfx_stage pm, void *user_data);
} tfx_stage_info_t;

typedef struct tfx_ribbon_dispatch_s {
	tfx_ribbon_bucket_t *ribbon_data;
	tfxU32 index_offset;
	tfxU32 vertex_offset;
	tfxU32 index_count;
	tfxU32 vertex_count;
	tfxU32 ribbon_offset;
	tfxU32 segment_offset;
	tfxU32 total_segments;
	tfxU32 last_index_offset;
	tfxU32 last_vertex_offset;
	tfxU32 last_ribbon_offset;
	tfxU32 last_segment_offset;
} tfx_ribbon_dispatch_t;

typedef struct tfx_ribbon_buffer_requirements_s {
	tfxU32 segment_buffer_size_in_bytes;
	tfxU32 ribbon_buffer_size_in_bytes;
	tfxU32 emitter_buffer_size_in_bytes;
} tfx_ribbon_buffer_requirements_t;

typedef struct tfx_animation_buffer_metrics_s {
	size_t sprite_data_size;
	tfxU32 offsets_size;
	tfxU32 instances_size;
	size_t offsets_size_in_bytes;
	size_t instances_size_in_bytes;
	tfxU32 total_sprites_to_draw;
	//Ribbon Metrics
	size_t ribbon_data_size;
	tfxU32 ribbon_offsets_size;
	tfxU32 ribbon_offsets_size_in_bytes;
	tfxU32 total_ribbons_to_draw;
	size_t ribbon_segment_data_size;
}tfx_animation_buffer_metrics_t;

//This can be sent as a push constant to the gpu
typedef struct tfx_ribbon_bucket_globals_s  {
	tfx_float32x4_t camera_position;
	tfxU32 segment_count;
	tfxU32 tessellation;  
	tfxU32 index_offset;
	tfxU32 vertex_offset;
	tfxU32 ribbon_count;
	tfxU32 ribbon_offset;
	tfxU32 segment_offset;
	tfxU32 uniform_index;
	tfxU32 emitters_index;
	tfxU32 graphs_index;
	tfxU32 ribbons_index;
	tfxU32 ribbon_segments_index;
	tfxU32 vertexes_index;
	tfxU32 indexes_index;
	tfxU32 image_data_index;
	tfxU32 sampler_index;
	tfxU32 particle_texture_index;
	tfxU32 color_ramp_texture_index;
	float lerp;
	float time;
	float ndc_offset_x;
	float ndc_offset_y;
} tfx_ribbon_bucket_globals_t;

typedef struct tfx_sprite_data_push_s {
	tfxU32 animation_instances_total;
	tfxU32 billboards_total;
	tfxU32 animated_shapes;	
	tfxU32 offsets_index;
	tfxU32 animation_instances_index;
	tfxU32 billboards_index;
	tfxU32 sprite_data_index;
	tfxU32 image_data_index;
	tfxU32 emitter_properties_index;
	tfxU32 bounding_boxes_index;
} tfx_sprite_data_push_t;

typedef struct tfx_gpu_graph_data_s {
	tfx_float32x4_t node_data;
	tfx_float32x4_t oscillator;
	tfx_graph_easing_type easing_type;
	int flags;
	int padding[2];
} tfx_gpu_graph_data_t;

//------------------------------------------------------------
//Section: Callback_typedefs
//------------------------------------------------------------
typedef void(*tfx_shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data);
typedef void(*tfx_uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset);
typedef bool(*tfx_maybe_render_instance_callback)(tfx_animation_manager animation_manager, tfx_float32x3_t position, float radius, void *user_data);


tfxAPI void tfx_UpdateAnimationManagerBufferMetrics(tfx_animation_manager animation_manager);
tfxAPI float tfx_DegreesToRadians(float degrees);
tfxAPI float tfx_RadiansToDegrees(float radians);

//--------------------------------
//Random numbers
//--------------------------------
tfxAPI tfx_random_t tfx_CreateRandom(tfxU32 seed);
tfxAPI void tfx_AdvanceRandom(tfx_random_t *random);
tfxAPI void tfx_RandomReseedTime(tfx_random_t *random);
tfxAPI void tfx_RandomReseed2(tfx_random_t *random, tfxU64 seed1, tfxU64 seed2);
tfxAPI void tfx_RandomReseed(tfx_random_t *random, tfxU64 seed);
tfxAPI float tfx_GenerateRandom(tfx_random_t *random);
tfxAPI float tfx_RandomRangeZeroToMax(tfx_random_t *random, float max);
tfxAPI float tfx_RandomRangeFromTo(tfx_random_t *random, float from, float to);
tfxAPI int tfx_RandomRangeFromToInt(tfx_random_t *random, int from, int to);
tfxAPI tfxU32 tfx_RandomRangeZeroToMaxUInt(tfx_random_t *random, tfxU32 max);
tfxAPI int tfx_RandomRangeZeroToMaxInt(tfx_random_t *random, int max);
tfxAPI void tfx_AlterRandomSeedU64(tfx_random_t *random, tfxU64 amount);
tfxAPI void tfx_AlterRandomSeedU32(tfx_random_t *random, tfxU32 amount);

//[API functions]
//All the functions below represent all that you will need to call to implement TimelineFX

/*
Quickstart
----------
The minimum lifecycle to get an effect on screen. Error handling and renderer specifics are omitted; every
tfx_* function used here is documented in full below.

    // 1. Initialise the library and its memory pool (once, at startup).
    tfx_BeginTimelineFX(tfx_GetDefaultThreadCount(), 128 * 1024 * 1024);

    // 2. Load an effects library from a .tfx file. shape_loader uploads each image to your renderer;
    //    uv_lookup returns its uv rect. Then build the GPU shape/graph/color-ramp data and upload it.
    tfx_library library = tfx_LoadEffectLibrary("effects.tfx", shape_loader, uv_lookup, user_data);

    // 3. Create a stage (the runtime effect manager) to simulate effects in.
    tfx_stage stage = tfx_CreateStage(tfx_CreateStageInfo(tfxStageSetup_none));

    // 4. Create a template for the effect you want, then add an instance of it to the stage. The returned
    //    tfxEffectID is your handle for manipulating that live instance (position, scale, expiry, ...).
    tfx_effect_template explosion = tfx_CreateEffectTemplate(library, "Explosion");
    tfxEffectID effect_id = tfx_AddEffectTemplateToStage(stage, explosion);

    // 5. Per-frame update loop, ideally in a fixed timestep (see tfx_UpdateStage for the full pattern).
    tfx_UpdateStage(stage, 16.66667);   // elapsed ms since last update

    // 6. Read the results for rendering. These getters complete the stage's in-flight work for you
    //    (see the threading contract in the Particle_Manager section). Copy the instance buffer to your
    //    GPU staging buffer and draw; drive ribbons with the tfx_*RibbonDispatch iterator.
    tfxU32 count = tfx_GetInstanceCount(stage);
    tfx_instance_t *instances = tfx_GetInstanceBuffer(stage);
    // ... upload `count` instances and issue your draw call ...

    // 7. Shutdown (see "Ownership and teardown ordering" below).
    tfx_FreeEffectTemplate(explosion);
    tfx_FreeStage(stage);
    tfx_FreeLibrary(library);
    tfx_EndTimelineFX();

Ownership and teardown ordering
-------------------------------
Everything is allocated from the memory pool created by tfx_BeginTimelineFX, so tfx_EndTimelineFX must be the
very last TimelineFX call you make — it tears down the worker threads and the pool itself.

tfx_EndTimelineFX sweeps up any stages, animation managers and libraries you did not free yourself (it completes
and frees every stage, then frees every animation manager, then every library). Freeing any of them yourself
first is fine — they deregister from the sweep as they go, so there is no double free either way. It does NOT
sweep effect templates, so you must free those explicitly with tfx_FreeEffectTemplate before calling
tfx_EndTimelineFX or the built-in leak check will report leaked blocks.

If you free objects yourself rather than leaving them to the sweep, free them in dependency order: free/complete
a stage before the library whose effects it is simulating (a live stage references the library's effect data);
an effect template holds its own detached clone of the effect so it is independent of the library once created,
but freeing templates before the library keeps the ordering simple. A stage has outstanding worker-thread work
after tfx_UpdateStage, so call tfx_CompleteStageWork (or use tfx_FreeStage, which completes it for you) before
tearing anything down mid-frame.
*/

//--------------------------------
//Initialisation_functions
//--------------------------------

/*
You don't have to call this, you can just call tfx_BeginTimelineFX in order to initialise the memory, but I created this for the sake of the editor which
needs to load in an ini file before initialising timelinefx which requires the memory pool to be created before hand
* @param memory_pool_size    The size of each memory pool to contain all objects created in TimelineFX, recommended to be at least 64MB
*/
tfxAPI void tfx_InitialiseTimelineFXMemory(size_t memory_pool_size);

/*
Initialise TimelineFX. Must be called before any functionality of TimelineFX is used.
* @param max_threads        The number of worker threads to use in addition to the main thread. Pass 0 to run single threaded.
*                            The count is clamped to the number of hardware threads available on the machine.
* @param memory_pool_size    The size of each memory pool to contain all objects created in TimelineFX, recommended to be at least 64MB
*/
tfxAPI void tfx_BeginTimelineFX(int max_threads, size_t memory_pool_size);

/*
Cleanup up all threads and memory used by timelinefx
*/
tfxAPI void tfx_EndTimelineFX(void);

//--------------------------------
//Global_variable_access
//--------------------------------
/*
Set the color format used for storing color ramps. Color ramps are generated by particle emitters and dictate how the particle colors change over the lifetime
of the particle. They can be uploaded to the GPU so you can set your preference for the color format as you need. The format should be immediately after you
call tfx_BeginTimelineFX
* @param color_format        The tfx_color_format to store color ramps in: tfx_color_format_rgba8, tfx_color_format_rgba16f or tfx_color_format_rgba32f
*/
tfxAPI void tfx_SetColorRampFormat(tfx_color_format color_format);

//--------------------------------
//Library_functions
//--------------------------------

/*
Create a new, empty library handle. Normally you will use tfx_LoadEffectLibrary to create a library
already populated from a tfx file; use this only if you want to build a library up from scratch.
Free it with tfx_FreeLibrary.
* @returns tfx_library    A handle to a new, empty library.
*/
tfxAPI tfx_library tfx_CreateLibrary(void);

/*
Count the number of particle shapes stored in a tfx file without loading the whole library. Useful
if you need to reserve storage for the shapes up front before loading.
* @param filename        The name of the tfx file to inspect
* @returns int            The number of shapes in the file.
*/
tfxAPI int tfx_GetShapeCountInLibrary(const char *filename);

/*
Validate a timelinefx tfx file to make sure that it's valid.
* @param filename        The name of the file where you want to count the number of shapes
* @returns int            Returns 0 if the file successfully validated or a tfxErrorFlags if something went wrong
*/
tfxAPI int tfx_ValidateEffectPackage(const char *filename);

/**
* Loads an effect library package from the specified filename into the provided tfx_library_t object.
*
* @param filename         A pointer to a null-terminated string that contains the path and filename of the effect library package to be loaded.
* @param shape_loader     A pointer to a function that will be used to load image data into the effect library package.
*                         The function has the following signature: void shape_loader(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data).
* @param user_data        A pointer to user-defined data that will be passed to the shape_loader function. This parameter is optional and can be set to nullptr if not needed.
* @return tfx_library	  A handle to a library object
*/
tfxAPI tfx_library tfx_LoadEffectLibrary(const char *filename, tfx_shape_loader shape_loader, tfx_uv_lookup uv_lookup, void *user_data);

/**
* Loads an effect library package from memory into the provided tfx_library_t object pointer.
*
* @param data             A pointer to a memory buffer containing the library to be loaded
* @param size             The size of the memory buffer containing the library to be loaded
* @param shape_loader     A pointer to a function that will be used to load image data into the effect library package.
*                         The function has the following signature: void shape_loader(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data).
* @param user_data        A pointer to user-defined data that will be passed to the shape_loader function. This parameter is optional and can be set to nullptr if not needed.
* @return tfx_library	  A handle to a library object
*/
tfxAPI tfx_library tfx_LoadEffectLibraryFromMemory(const void *data, tfxU32 size, tfx_shape_loader shape_loader, tfx_uv_lookup uv_lookup, void *user_data);

/*
Get the error flags from a library. When you load a library from file or memory, if something goes wrong then the error status is stored in the library object and you can retrieve it with this command.
* @param lib            A handle to a tfx_library object that will hold the loaded effect library data.
* @return A tfxErrorFlags value that indicates whether the function succeeded or failed. The possible return values are:
	tfxErrorCode_success = 0
	tfxErrorCode_incorrect_package_format
	tfxErrorCode_data_could_not_be_loaded
	tfxErrorCode_could_not_add_shape
	tfxErrorCode_error_loading_shapes
	tfxErrorCode_some_data_not_loaded
	tfxErrorCode_unable_to_open_file
	tfxErrorCode_unable_to_read_file
	tfxErrorCode_wrong_file_size
	tfxErrorCode_invalid_format
	tfxErrorCode_no_inventory
	tfxErrorCode_invalid_inventory
*/
tfxAPI tfxErrorFlags tfx_GetLibraryErrorStatus(tfx_library library);

/**
* Loads a sprite data file into an animation manager
*
* @param filename        A pointer to a null-terminated string that contains the path and filename of the effect library package to be loaded.
* @param lib            A reference to a tfx_animation_manager_t object that will hold the loaded sprite data.
* @param shape_loader    A pointer to a function that will be used to load image data into the effect library package.
*                        The function has the following signature: void shape_loader(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data).
* @param user_data        A pointer to user-defined data that will be passed to the shape_loader function. This parameter is optional and can be set to nullptr if not needed.
*
* @return A tfxErrorFlags value that indicates whether the function succeeded or failed. The possible return values are:
	tfxErrorCode_success = 0
	tfxErrorCode_incorrect_package_format
	tfxErrorCode_data_could_not_be_loaded
	tfxErrorCode_could_not_add_shape
	tfxErrorCode_error_loading_shapes
	tfxErrorCode_some_data_not_loaded
	tfxErrorCode_unable_to_open_file
	tfxErrorCode_unable_to_read_file
	tfxErrorCode_wrong_file_size
	tfxErrorCode_invalid_format
	tfxErrorCode_no_inventory
	tfxErrorCode_invalid_inventory
*/
tfxAPI tfxErrorFlags tfx_LoadSpriteData(const char *filename, tfx_animation_manager animation_manager, tfx_shape_loader shape_loader, void *user_data);

/*
* Updates all the image data in the library using the uv_lookup that you set when loading a library. This allows you to add all of the uv data for
* the shapes that are loaded into the texture. You must have set the uv_lookup callback when loading the library, otherwise you can loop over the 
* shapes in the library and update the data yourself using the tfx_GetLibraryShapeArray and related functions.
* @param tfx_library                A valid pointer to a tfx_library
*/
tfxAPI void tfx_UpdateLibraryGPUImageData(tfx_library library);

/*
Get the number of shapes stored in the library
* @param tfx_library                A valid pointer to a tfx_library
* @return tfxU32					Count of shapes
*/
tfxAPI tfxU32 tfx_GetLibraryImageCount(tfx_library library);

/*
Get a particle image from a library by it's index
* @param tfx_library                A valid pointer to a tfx_library
* @return image						A tfx_image_data_t object with all the details of the image
*/
tfxAPI tfx_image_data_t *tfx_GetLibraryImage(tfx_library library, tfxU32 index);

tfxAPI void tfx_SetImagePointer(tfx_image_data_t *image, void *pointer);
tfxAPI void tfx_SetGPUImageTextureInfo(tfx_gpu_image_data_t *image, float x, float y, float z, float w, int array_index);
tfxAPI int tfx_GetImageFrameCount(tfx_image_data_t *image);
tfxAPI int tfx_GetImageWidth(tfx_image_data_t *image);
tfxAPI int tfx_GetImageHeight(tfx_image_data_t *image);
tfxAPI void* tfx_GetBitmapData(tfx_bitmap_t *bitmap);
tfxAPI size_t tfx_GetBitmapSize(tfx_bitmap_t *bitmap);
tfxAPI int tfx_GetBitmapWidth(tfx_bitmap_t *bitmap);
tfxAPI int tfx_GetBitmapHeight(tfx_bitmap_t *bitmap);

/*
Output all the effect names in a library to the console
* @param tfx_library                A valid pointer to a tfx_library
*/
tfxAPI void tfx_ListEffectNames(tfx_library library);

/*
Get an effect in the library by it's index. If you need to get an effect in a folder or an emitter then you can use tfx_GetLibraryEffectPath instead.
* @param tfx_library                A valid pointer to a tfx_library
*/
tfxAPI tfx_effect_descriptor tfx_GetEffectByIndex(tfx_library library, int index);

/*
Get an effect in the library by it's path. So for example, if you want to get a pointer to the emitter "spark" in effect "explosion" then you could do GetEffect("explosion/spark")
You will need this function to apply user data and update callbacks to effects and emitters before adding the effect to the effect manager
* @param tfx_library_t                A valid pointer to a tfx_library_t
* @param const char *path             Path to the effect or emitter
*/
tfxAPI tfx_effect_descriptor tfx_GetLibraryEffectPath(tfx_library library, const char *path);

/*
Check whether a path resolves to an effect or emitter in the library. tfx_GetLibraryEffectPath asserts
in debug builds and returns NULL in release builds when a path is missing, so use this to probe a path
safely before looking it up.
* @param library                A valid pointer to a tfx_library_t
* @param path                   Path to the effect or emitter
* @returns                      true if the path exists in the library
*/
tfxAPI bool tfx_IsValidEffectPath(tfx_library library, const char *path);

/*
Free all the memory used by a library
* param tfx_library				A pointer to the library that you want to free
*/
tfxAPI void tfx_FreeLibrary(tfx_library library);

/*
Create the image data required for shaders from a TimelineFX library. The image data will contain data such as uv coordinates. Once you have built the data you can use GetLibraryImageData to get the buffer
and upload it to the gpu.
* @param library                  A pointer to a tfx_library_t object
* @param shapes                   A pointer to a tfx_gpu_shapes_t object which will fill a buffer with all the shapes
* @param uv_lookup                A function pointer to a function that you need to set up in order to get the uv coordinates from whatever renderer you're using
*/
tfxAPI void tfx_BuildLibraryGPUShapeData(tfx_library library, tfx_gpu_shapes shapes, tfx_uv_lookup uv_lookup);

/*
Get a pointer to the particle shapes data in a library. This can be used with tfx_BuildGPUShapeData when you want to upload the data to the GPU
* @param library        A pointer to a tfx_library_t
* @param count			A pointer to an int that will be filled with the nubmer of images in the image data array that's returned
*/
tfxAPI tfx_image_data_t *tfx_GetParticleShapesLibrary(tfx_library library, int *count);

/*
Get a count of the number of color ramp bitmaps in the library. Color ramps are used to change the color of particles over time and you will need to upload them to the GPU.
* @param library        A pointer to a tfx_library_t
*/
tfxAPI tfxU32 tfx_GetColorRampBitmapCount(tfx_library library);

/*
Get a pointer to a color ramp bitmap in a library. You can use this data to upload the bitmaps to the GPU.
* @param library        A pointer to a tfx_library_t
*/
tfxAPI tfx_bitmap_t *tfx_GetColorRampBitmap(tfx_library library, tfxU32 index);

/*
Get a count of the number of color ramp bitmaps in an animation manager.
* @param animation_manager        A handle to a tfx_animation_manager
*/
tfxAPI tfxU32 tfx_GetAnimationColorRampBitmapCount(tfx_animation_manager animation_manager);

/*
Get a pointer to a color ramp bitmap in an animation manager. You can use this data to upload the bitmaps to the GPU.
* @param animation_manager        A handle to a tfx_animation_manager
* @param index                    The index of the bitmap
*/
tfxAPI tfx_bitmap_t *tfx_GetAnimationColorRampBitmap(tfx_animation_manager animation_manager, tfxU32 index);

/*
Check to see if a library has been initialised or not
* @param library        A pointer to a tfx_library_t
*/
tfxAPI bool tfx_LibraryIsInitialised(tfx_library library);

/*
Get the gpu shapes handle in library. The gpu shapes handle can be used to upload the image data for particle shapes to the gpu using functions like tfx_GetGPUShapesArray, tfx_GetGPUShapesSizeInBytes, 
tfx_GetGPUShapesCount etc.
* @param library        A handle to a tfx_library
*/
tfxAPI tfx_gpu_shapes tfx_GetLibraryGPUShapes(tfx_library library);

/*
Get a pointer to the global gpu graph lookup data. This is a single shared buffer maintained across all
libraries in tfxStore, so it takes no library argument. Shaders use this lookup data to update attributes of
particles and ribbons (currently ribbons only). Use with tfx_GetGPUGraphLookupsBufferSizeInBytes to upload it
to your GPU buffer.
* @returns tfx_gpu_graph_data_t *   A pointer to the shared gpu graph lookup data.
*/
tfxAPI tfx_gpu_graph_data_t *tfx_GetGPUGraphLookupsBuffer(void);

/*
Get the size in bytes of the global gpu graph lookup data returned by tfx_GetGPUGraphLookupsBuffer. This is the
single shared buffer in tfxStore, so it takes no library argument.
* @returns tfxU32                   The size of the shared gpu graph lookup data in bytes.
*/
tfxAPI tfxU32 tfx_GetGPUGraphLookupsBufferSizeInBytes(void);

/*
Get buffer info for ribbons based on the tessellation value. Returns a tfx_ribbon_buffer_info_t object with
vertices per segment, triangles per segment and indices per segment requirements.
* @param tfxU32        					The number of tessellations for the ribbons
* @returns tfx_ribbon_buffer_info_t		Info containing vertices, triangles and indices per segment
*/
tfxAPI tfx_ribbon_buffer_info_t tfx_GenerateRibbonBufferInfo(tfxU32 tessellation);

//--------------------------------
//Particle_Manager_functions
//--------------------------------

/*
Threading contract
------------------
A stage updates particles across the worker threads created by tfx_BeginTimelineFX. tfx_UpdateStage kicks off
that work and returns immediately while the work is still in flight on the worker threads. This means that after
tfx_UpdateStage returns, the sprite/instance and ribbon buffers are NOT yet safe to read.

Two ways to synchronise before you read that data:

  - Call tfx_CompleteStageWork(pm) explicitly. This blocks until all of the stage's outstanding update work has
    finished, after which every buffer is safe to read.
  - Or rely on the getters that synchronise implicitly. The following complete the stage's outstanding work for
    you before returning, so you can call them directly after tfx_UpdateStage without calling
    tfx_CompleteStageWork first:
        tfx_GetInstanceBuffer, tfx_GetInstanceBufferByLayer, tfx_GetInstanceCount, tfx_GetInstanceCountByLayer,
        tfx_GetNextInstanceBuffer, tfx_GetRibbonBuffers, tfx_HasRibbonsToDraw, tfx_NextRibbonDispatch,
        tfx_GetRibbonBufferRequirements, tfx_CopyRibbonDataToStagingBuffers, tfx_ClearStage and tfx_FreeStage.

Anything that reads stage state through a path NOT in that list (for example reading the instance buffer pointer
you cached last frame, or touching the raw buffers directly) must be preceded by tfx_CompleteStageWork, otherwise
you will race the worker threads. Mutators such as tfx_SetEffectPosition are also only safe once the update work
for the frame has completed so it's a good idea to set positions before the call to tfx_UpdateStage. 
*/

/*
Create a tfx_stage_info_t object which contains configuration data that you can pass to tfx_CreateStage to setup a effect manager. You can tweak the config after calling this
function if needed to fine tune the settings.
* @param setup                    A tfx_stage_setup enum which you can use to set the info based on some commonly used templates
*/
tfxAPI tfx_stage_info_t tfx_CreateStageInfo(tfx_stage_setup setup);

/*
Create a stage (effect manager) configured with a tfx_stage_info_t object. Build the info with tfx_CreateStageInfo
and tweak it as needed before passing it in. If you later want to run the stage in a different mode then call
tfx_ReconfigureStage rather than creating a new one.
* @param info                   A tfx_stage_info_t containing the configuration for the stage.
* @returns tfx_stage            A handle to the new stage, or NULL if the allocation failed (out of memory).
*/
tfxAPI tfx_stage tfx_CreateStage(tfx_stage_info_t info);

/*
Reconfigure a effect manager to make it work in a different mode. A effect manager can only run in a single mode at time like unordered, depth ordered etc so use this to change that. Also bear
in mind that you can just use more than one effect manager and utilised different modes that way as well. The modes that you need will depend on the effects that you're adding to the effect manager.
* @param pm                       A pointer to an intialised tfx_stage_t.
* @param sort_passes              The number of sort passes if you're using depth sorted effects
*/
tfxAPI void tfx_ReconfigureStage(tfx_stage pm, tfxU32 sort_passes);

/*
Block until all of the stage's in-flight update work (started by tfx_UpdateStage) has finished. Call this before
reading the stage's sprite/instance or ribbon buffers, or mutating effects, unless you are only using the getters
that synchronise implicitly (see the threading contract above). Safe to call redundantly; it is a no-op if there
is no work outstanding.
* @param pm                       A handle to an initialised tfx_stage_t.
*/
tfxAPI void tfx_CompleteStageWork(tfx_stage pm);

/*
Turn on and off whether the effect manager should sort the effects by depth order. Use tfx_SetStageCamera to set the position of the camera that the effect manager will
use to update the depth of each effect in the scene.
* @param pm                       A pointer to an intialised tfx_stage_t.
* @param yesno                    A boolean, set to true or false if you want auto ordering on or off respectively
*/
tfxAPI void tfx_ToggleStageOrderEffects(tfx_stage pm, bool yesno);

/*
Get the billboard buffer in the effect manager containing all the sprite instances that were created in the most recent frame. You can use this to copy to a staging buffer to upload to the gpu.
* @param pm                       A pointer to an intialised tfx_stage_t.
*/
tfxAPI tfx_instance_t *tfx_GetInstanceBuffer(tfx_stage  pm);

/*
Get the billboard buffer in the effect manager containing all the sprite instances for a specific layer that were created in the most recent frame. You can use this to copy to a staging buffer to upload to the gpu.
You can then use tfx_GetInstanceCountByLayer for the draw call.
* @param pm                       A pointer to an intialised tfx_stage_t.
*/
tfxAPI tfx_instance_t *tfx_GetInstanceBufferByLayer(tfx_stage pm, tfxU32 layer);

/*
Get the number of instances within the instance buffer of a effect manager
* @param pm                       A pointer to an intialised tfx_stage_t.
*/
tfxAPI int tfx_GetInstanceCount(tfx_stage pm);

/*
Get the number of instances within the instance buffer of a effect manager for a specific layer.
* @param pm                       A pointer to an intialised tfx_stage_t.
*/
tfxAPI int tfx_GetInstanceCountByLayer(tfx_stage pm, tfxU32 layer);

/*
Get the update time being used by the effect manager.
* @param pm                       A handle to an intialised tfx_stage_t.
*/
tfxAPI double tfx_GetUpdateTime(tfx_stage pm);

/*
Get the ribbon buffer for a given segment size. This will give you all the necessary info and buffer pointers for uploading the ribbon data to the GPU for processing and converting into
a vertex buffer for rendering.
* @param pm                       A pointer to an intialised tfx_stage_t.
* @param segment_count            An unsigned int specifying the ribbon length that you want the rendering info for.
* @returns						  A pointer to a tfx_ribbon_bucket_t
*/
tfxAPI tfx_ribbon_bucket_t *tfx_GetRibbonBuffers(tfx_stage pm, tfxKey bucket_id);

/*
Call this to determine whether or not any effect manager has ribbon_emitters to draw this frame.
* @returns						  True or false
*/
tfxAPI bool tfx_HasRibbonsToDraw(tfx_stage pm);

/*
Get a struct containing the info you need to compute and render ribbon_emitters of a specific length. 
* @param pm                       A pointer to an intialised tfx_stage_t.
* @param segment_count            An unsigned int specifying the ribbon length that you want the rendering info for.
* @returns						  A tfx_ribbon_buffer_info_t struct
*/
tfxAPI tfx_ribbon_buffer_info_t tfx_GetRibbonBufferInfo(tfx_stage pm, tfxKey bucket_id);

/*
--------------------------------
Ribbon rendering / GPU dispatch
--------------------------------
Ribbons are turned into geometry on the GPU by a compute shader and then drawn. Unlike billboard
sprites, which you upload wholesale from a single instance buffer, ribbons are batched into "buckets" (one per
distinct segment length / tessellation) and you drive the upload and dispatch by iterating those buckets.

Per-frame flow, after tfx_UpdateStage:

  1. Size and (re)allocate your GPU buffers. Use tfx_GetRibbonBufferRequirements() for the exact bytes needed
     this frame across every stage, or the tfx_Get*MaxSizeInBytes / tfx_GetTotal*MaxSizeInBytes helpers to size
     for the worst case up front so you never resize.
  2. Copy the CPU ribbon data into your mapped staging buffers with tfx_CopyRibbonDataToStagingBuffers, then
     transfer them to your device-local segment / ribbon / emitter buffers.
  3. Iterate the buckets to record compute dispatches: create a tfx_ribbon_dispatch_t with
     tfx_CreateRibbonDispatch(), then loop while tfx_NextRibbonDispatch(pm, &dispatch) returns true. Each call
     fills in the per-bucket offsets/counts and populates dispatch.ribbon_data->globals with the push-constant
     block; dispatch total_segments tells you how many threads to launch.
  4. Iterate the buckets again the same way to record the draw calls, using dispatch.index_count and
     dispatch.index_offset for the indexed draw.

The dispatch iterator advances an internal cursor on the stage, so if you need to walk the buckets more than
once in a frame (compute pass, then render pass) call tfx_ResetRibbonDispatchIterator(pm) between passes.
tfx_NextRibbonDispatch / tfx_GetRibbonBufferRequirements / tfx_CopyRibbonDataToStagingBuffers all implicitly
complete any in-flight stage update work first (see the threading contract above), so they are safe to call
straight after tfx_UpdateStage. See the RibbonComputeFunction / RenderRibbons functions in the editor for a
complete worked example.
*/

/*
Create a zero-initialised tfx_ribbon_dispatch_t to drive the ribbon bucket iteration. Pass its address to
tfx_NextRibbonDispatch. The struct carries the running offsets between buckets, so create a fresh one at the
start of each iteration pass (or reset the fields to zero).
* @returns tfx_ribbon_dispatch_t   A zeroed dispatch struct ready to pass to tfx_NextRibbonDispatch.
*/
tfxAPI tfx_ribbon_dispatch_t tfx_CreateRibbonDispatch(void);

/*
Advance to the next ribbon bucket that has ribbons to draw. Call in a while loop until it returns false. On each
true return, ribbon_dispatch is filled with this bucket's vertex/index offsets and counts, total_segments (the
work item count for the compute dispatch) and a pointer to the bucket in ribbon_data whose globals member is the
push-constant block for both the compute and render passes. Implicitly completes any in-flight stage update.
* @param pm                        A handle to an initialised tfx_stage_t.
* @param ribbon_dispatch           A pointer to a tfx_ribbon_dispatch_t (from tfx_CreateRibbonDispatch) updated in place.
* @returns bool                    true while there is another bucket to process, false when iteration is finished.
*/
tfxAPI bool tfx_NextRibbonDispatch(tfx_stage pm, tfx_ribbon_dispatch_t *ribbon_dispatch);

/*
Reset the stage's internal ribbon bucket iterator back to the start. Call this between two passes over the
buckets in the same frame (for example after the compute pass, before the render pass) so tfx_NextRibbonDispatch
walks them again from the beginning.
* @param pm                        A handle to an initialised tfx_stage_t.
*/
tfxAPI void tfx_ResetRibbonDispatchIterator(tfx_stage pm);

tfxAPI tfx_ribbon_bucket_globals_t *tfx_GetRibbonDispatchGlobals(tfx_ribbon_dispatch_t *ribbon_dispatch);

/*
Get the exact ribbon buffer sizes needed this frame, summed across every registered stage, for allocating a
single shared set of GPU buffers. Implicitly completes any in-flight stage update on each stage. Use this when
you would rather allocate to the exact per-frame need than to the worst case (see tfx_GetTotal*MaxSizeInBytes).
* @returns tfx_ribbon_buffer_requirements_t   Segment, ribbon and emitter buffer sizes in bytes for this frame.
*/
tfxAPI tfx_ribbon_buffer_requirements_t tfx_GetRibbonBufferRequirements(void);

/*
Copy this stage's ribbon segment, ribbon instance and emitter data into your mapped staging buffers, ready to
transfer to the GPU. Size the destinations with tfx_GetRibbonBufferRequirements or the max-size helpers.
Implicitly completes any in-flight stage update.
* @param stage                     A handle to an initialised tfx_stage_t.
* @param segments_dst              Destination for the ribbon segment data (must be non-null).
* @param ribbons_dst               Destination for the ribbon instance data (must be non-null).
* @param emitters_dst              Destination for the ribbon emitter data (must be non-null).
*/
tfxAPI void tfx_CopyRibbonDataToStagingBuffers(tfx_stage stage, void *segments_dst, void *ribbons_dst, void *emitters_dst);

/*
Worst-case size in bytes of the ribbon segment buffer for a single stage, based on its configured max_ribbon_segments.
Use to size a buffer once at startup instead of resizing per frame.
* @param pm                        A handle to an initialised tfx_stage_t.
*/
tfxAPI size_t tfx_GetSegmentBufferMaxSizeInBytes(tfx_stage pm);

/*
Worst-case size in bytes of the ribbon vertex buffer for a single stage. Pass the size of your own vertex struct,
or 0 to use the library's tfx_ribbon_vertex_t size.
* @param pm                        A handle to an initialised tfx_stage_t.
* @param vertex_size               The size of a single vertex in bytes, or 0 for the default tfx_ribbon_vertex_t.
*/
tfxAPI size_t tfx_GetSegmentVertexBufferMaxSizeInBytes(tfx_stage pm, tfxU32 vertex_size);

/*
Worst-case size in bytes of the ribbon index buffer for a single stage, based on its configured max_ribbon_segments
and tessellation.
* @param pm                        A handle to an initialised tfx_stage_t.
*/
tfxAPI size_t tfx_GetSegmentIndexBufferMaxSizeInBytes(tfx_stage pm);

/*
Worst-case size in bytes of the ribbon instance buffer for a single stage, based on its configured max_ribbons.
* @param pm                        A handle to an initialised tfx_stage_t.
*/
tfxAPI size_t tfx_GetRibbonBufferMaxSizeInBytes(tfx_stage pm);

/*
Worst-case size in bytes of the ribbon emitter buffer for a single stage.
* @param pm                        A handle to an initialised tfx_stage_t.
*/
tfxAPI size_t tfx_GetEmitterBufferMaxSizeInBytes(tfx_stage pm);

/*
Get the size in bytes of the per-emitter particle properties buffer for a library. These properties are looked up
by the particle and ribbon shaders; use with tfx_GetParticlePropertiesBuffer to upload them to the GPU.
* @param library                  A handle to a tfx_library.
*/
tfxAPI size_t tfx_GetParticlePropertiesBufferSizeInBytes(tfx_library library);

/*
Get a pointer to the per-emitter particle properties buffer for a library, for uploading to the GPU. Size it with
tfx_GetParticlePropertiesBufferSizeInBytes.
* @param library                  A handle to a tfx_library.
*/
tfxAPI void *tfx_GetParticlePropertiesBuffer(tfx_library library);

/*
Get the total worst-case buffer sizes across all registered stages, for creating a single shared set of ribbon GPU
buffers. These iterate every registered stage and sum the individual tfx_Get*MaxSizeInBytes requirements. Use these
to allocate for the worst case once; use tfx_GetRibbonBufferRequirements instead to size to the exact per-frame need.
*/
tfxAPI size_t tfx_GetTotalSegmentBufferMaxSizeInBytes(void);

tfxAPI size_t tfx_GetTotalSegmentVertexBufferMaxSizeInBytes(tfxU32 vertex_size);

tfxAPI size_t tfx_GetTotalSegmentIndexBufferMaxSizeInBytes(void);

tfxAPI size_t tfx_GetTotalRibbonBufferMaxSizeInBytes(void);

tfxAPI size_t tfx_GetTotalEmitterBufferMaxSizeInBytes(void);

/*
When a effect manager updates particles it creates work queues to handle the work. By default these each have a maximum amount of 1000 entries which should be
more than enough for most situations. However you can increase the sizes here if needed. You only need to set this manually if you hit one of the asserts when these
run out of space or you anticipate a huge amount of emitters and particles to be used (> million). On the other hand, you might be tight on memory in which case you
could reduce the numbers as well if needed (they don't take a lot of space though)
* @param pm                        A pointer to an intialised tfx_stage_t.
* @param spawn_work_max            The maximum amount of spawn work entries
* @param control_work_max          The maximum amount of control work entries
* @param age_work_max              The maximum amount of age_work work entries
*/
tfxAPI void tfx_SetStageWorkQueueSizes(tfx_stage pm, tfxU32 spawn_work_max, tfxU32 control_work_max, tfxU32 age_work_max);

/*Free the memory for a specific emitter type. When an emitter is created it creates memory to store all of the particles that it updates each frame. If you have
multiple emitters of the same type then their particle lists are resused rather then freed as they expire. When they're freed then the unused list is added to a list
of free particle banks for that emitter type so that they can then be recycled if another emitter of the same type is created. If you want to free the memory for a
specific emitter then you can call this function to do that.
NOTE: No emitters of the type passed to the function must be in use in the effect manager.
* @param pm                        A pointer to an intialised tfx_stage_t.
* @param emitter                   A pointer to a valid tfx_effect_descriptor_t of type tfxEmitterType
*/
tfxAPI void tfx_FreeParticleListsMemory(tfx_stage pm, tfx_effect_descriptor emitter);

/*
Free all the memory that is associated with an effect. Depending on the configuration of the effect manager this might be instance_data, particle lists and spawn location lists.
* @param pm                        A pointer to an intialised tfx_stage_t.
* @param emitter                   A pointer to a valid tfx_effect_descriptor_t of type tfxEffectType
*/
tfxAPI void tfx_FreeEffectListsMemory(tfx_stage pm, tfx_effect_descriptor effect);

/*
Get the current particle count for a effect manager
* @param pm                        A pointer to an tfx_stage_t
* @returns tfxU32                  The total number of particles currently being updated
*/
tfxAPI tfxU32 tfx_GetParticleCount(tfx_stage pm);

/*
Get the current ribbon count for a effect manager
* @param pm                        A pointer to an tfx_stage_t
* @returns tfxU32                  The total number of particles currently being updated
*/
tfxAPI tfxU32 tfx_GetRibbonCount(tfx_stage pm);

/*
Get the current number of effects that are currently being updated by a effect manager
* @param pm                        A pointer to an tfx_stage_t
* @returns tfxU32                  The total number of effects currently being updated
*/
tfxAPI tfxU32 tfx_GetEffectCount(tfx_stage pm);

/*
Get the current number of emitters that are currently being updated by a effect manager
* @param pm                        A pointer to an tfx_stage_t
* @returns tfxU32                  The total number of emitters currently being updated
*/
tfxAPI tfxU32 tfx_GetEmitterCount(tfx_stage pm);

/*
Get the current number of free effects in the stage that are ready to be re-used
* @param pm                        A pointer to an tfx_stage_t
* @returns tfxU32                  The total number of free effects in the stage
*/
tfxAPI tfxU32 tfx_GetFreeEffectCount(tfx_stage pm);

/*
Set the seed for the effect manager for random number generation. Setting the seed can determine how an emitters spawns particles, so if you set the seed before adding an effect to the effect manager
then the effect will look the same each time. Note that seed of 0 is invalid, it must be 1 or greater.
* @param pm                        A pointer to an initialised tfx_stage_t. 
* @param seed                      An unsigned int representing the seed (Any value other then 0)
*/
tfxAPI void tfx_SetStageSeed(tfx_stage pm, tfxU64 seed);

/*
Add an effect to a tfx_stage_t from an effect template
* @param pm                         A pointer to an initialised tfx_stage_t.
* @param effect_template			The tfx_effect_template_t object that you want to add to the effect manager. It must have already been prepared by calling tfx_CreateEffectTemplate
* @param effect_id					pointer to a tfxEffectID of the effect which will be set after it's been added to the effect manager. This index can then be used to manipulate the effect in the effect manager as it's update
									For example by calling tfx_SetEffectPosition. This will be set to tfxINVALID if the function is unable to add the effect to the effect manager if it's out of space and reached it's effect limit.
  @returns							True if the effect was succesfully added.
*/
tfxAPI tfxEffectID tfx_AddEffectTemplateToStage(tfx_stage pm, tfx_effect_template effect);

/*
Add an effect to a tfx_stage_t. Generally you should always call tfx_AddEffectTemplateToStage and use templates to organise your effects but if you want to just
test things out you can add an effect direct from a library using this command.
* @param pm							A pointer to an initialised tfx_stage_t. 
* @param effect						tfx_effect_descriptor_t object that you want to add to the effect manager.
* @param effect_id					pointer to a tfxEffectID of the effect which will be set after it's been added to the effect manager. This index can then be used to manipulate the effect in the effect manager as it's update
									For example by calling tfx_SetEffectPosition. This will be set to tfxINVALID if the function is unable to add the effect to the effect manager if it's out of space and reached it's effect limit.
  @returns							True if the effect was succesfully added.
*/
tfxAPI tfxEffectID tfx_AddRawEffectToStage(tfx_stage pm, tfx_effect_descriptor effect);

tfxAPI bool tfx_EffectIDIsValid(tfxEffectID id);
tfxAPI bool tfx_AnimationIDIsValid(tfxAnimationID id);

/*
Set the warmup time for an effect template. When the effect template is added to the stage it will be warmed up
first so that effectively it begins it's simulation from a set point in time. 
* @param tfx_effect_template		A handle to the effect template
* @param float						The number of milliseconds to warmup for
*/
tfxAPI void tfx_SetEffectTemplateWarmupTime(tfx_effect_template effect, float millisecs);

/*
Set the delta time used when warming up effects. Either match the fixed time step you're using like 1000/60 for 60fps 
or used mulitples of that value for higher performance at the cost of precision
* @param pm							A pointer to an initialised tfx_stage_t. 
* @param double						The delta time measured in milliseconds
*/
tfxAPI void tfx_SetWarmUpDeltaTime(tfx_stage pm, double delta_time);

/*
If an effect is already in the stage and you want to advance it be a given amount of time
without rendering those frame then you can use this function.
* @param pm							A pointer to an initialised tfx_stage_t. 
* @param tfxEffectID				The effect id that is in the stage
* @param float						The time to advance by measured in milliseconds
*/
tfxAPI void tfx_AdvanceEffectTime(tfx_stage pm, tfxEffectID effect_id, float time);

/*
Update a stage which manages effects. If you are interpolating particles in the vertex shader then it's important to only call this function once per frame only and idealy in a fixed step loop.
That means that if your fixed loop has to run twice to catch up (because of low frame rates) then you should still only call this function once but you can multiply the elapsed time
by the number of ticks. The ellapsed time should be the amount of time that has passed since the last frame so in a fixed step loop this will simply be the update rate in millisecs.
For example if you're updating 60 frames per second then elapsed time would be 16.666667. Psuedo code would look something like this:

	TimerAccumulate(game->timer);
	int pending_ticks = TimerPendingTicks(game->timer);	//The number of times the update loop will run this frame.

	while (tfx_TimerDoUpdate(game->timer)) {
		if (pending_ticks > 0) {
			tfx_UpdateStage(&game->pm, FrameLength * pending_ticks);
			//Set the pending ticks to 0 so we don't run the update again this frame.
			pending_ticks = 0;
		}

		TimerUnAccumulate(game->timer);
	}
	TimerSet(game->timer);	//Set the timer and calculate the interpolation value. You can pass that to a uniform or push constant for the shader

	//Only upload the sprite/billboard buffer to the gpu if the effect manager was updated.
	if (TimerUpdateWasRun(game->timer)) {
		RenderParticles(game->pm, game);
	}

* @param pm                    A pointer to an initialised tfx_stage_t.
* @param double                the amount of time that elapsed since the last frame
*/
tfxAPI void tfx_UpdateStage(tfx_stage pm, double elapsed);

/*
Get the image pointer for a sprite. Use this when rendering particles in your renderer. The pointer that is returned will be the pointer that you set in your shape loader function
used when loading an effect library. Generally you shouldn't need to use this function, simply copy the whole instance buffer in the effect manager to your staging buffer to be
copied to the gpu in one go.
* @param pm                    A pointer to an initialised tfx_stage_t. 
* @param property_indexes    The value in the instance_data->property_indexs[i] when iterating over the instance_data in your render function
  @returns                    void* pointer to the image
*/
tfxAPI void *tfx_GetSpriteImagePointer(tfx_stage pm, tfxU32 property_indexes);

/*
Get the total number of instances ready for rendering in the effect manager.
* @param pm                    A pointer to an initialised tfx_stage_t.
*/
tfxAPI tfxU32 tfx_TotalSpriteCount(tfx_stage pm);

/*
Clear all particles, instance_data and effects in a effect manager. If you don't need to use the effect manager again then call tfx_FreeStage to also
free all the memory associated with the effect manager.
* @param pm                        A pointer to an initialised tfx_stage_t.
* @param free_particle_banks    Set to true if you want to free the memory associated with the particle banks and release back to the memory pool
*/
tfxAPI void tfx_ClearStage(tfx_stage pm, bool free_particle_banks, bool free_sprite_buffers);

/*
Free all the memory used in the effect manager. You don't have to call this - tfx_EndTimelineFX frees any stages
that are still alive - but you can free them earlier if you want to reclaim the memory. Calling it is safe either
way, the stage is deregistered from the shutdown sweep so it won't be freed twice.
* @param pm                        A pointer to an initialised tfx_stage_t.
*/
tfxAPI void tfx_FreeStage(tfx_stage pm);

//[Effects functions for altering effects that are currently playing out in a effect manager]

/*
Expire an effect by telling it to stop spawning particles. This means that the effect will eventually be removed from the effect manager after all of it's remaining particles have expired.
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToStage
*/
tfxAPI void tfx_SoftExpireEffect(tfx_stage pm, tfxEffectID effect_index);

/*
Soft expire all the effects in a effect manager so that the particles complete their animation first
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
*/
tfxAPI void tfx_SoftExpireAll(tfx_stage pm);

/*
Expire an effect by telling it to stop spawning particles and remove all associated particles immediately.
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToStage
*/
tfxAPI void tfx_HardExpireEffect(tfx_stage pm, tfxEffectID effect_index);

/*
Get effect user data
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToStage
* @returns                void* pointing to the user data set in the effect. See tfx_effect_template_t::SetUserData() and tfx__set_effect_user_data()
*/
tfxAPI void *tfx_GetEffectUserData(tfx_stage pm, tfxEffectID effect_index);

/*
More for use in the editor, this function updates emitter base values for any effects that are currently running after their graph values have been changed.
*/
tfxAPI void tfx_UpdateStageBaseValues(tfx_stage pm);

/*
Set the effect manager camera. This is used to calculate particle depth if you're using depth ordered particles so it needs to be updated each frame.
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param front            An array of 3 floats representing a normalised 3d vector describing the direction that the camera is pointing
* @param position        An array of 3 floats representing the position of the camera in 3d space
*/
tfxAPI void tfx_SetStageCamera(tfx_stage pm, float front[3], float position[3]);

/*
Each effect in the effect manager can have bounding box which you can decide to keep updated or not if you wanted to do any offscreen culling of effects. Theres some
extra overhead to keep the bounding boxes updated but that can be made back if you have a number of effect particles offscreen that don't need to be drawn.
* @param pm					A pointer to a tfx_stage_t where the effect is being managed
* @param yesno				Set to true or false if you want the bounding boxes to be udpated.
*/
tfxAPI void tfx_KeepBoundingBoxesUpdated(tfx_stage pm, bool yesno);

/*
Set the effect user data for an effect already added to a effect manager
* @param pm					A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index		The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param user_data			A void* pointing to the user_data that you want to store in the effect
*/
tfxAPI void tfx_SetEffectUserData(tfx_stage pm, tfxU32 effect_index, void *data);

/*
Force a effect manager to only run in single threaded mode. In other words, only use the main thread to update particles
* @param pm                A pointer to a tfx_stage_t.
* @param switch_on        true or false to use a single thread or not
*/
tfxAPI void tfx_ForceStageSingleThreaded(tfx_stage pm, bool switch_on);

/*
Get the transform vectors for a sprite's previous position so that you can use that to interpolate between that and the current sprite position
* @param pm                A pointer to a tfx_stage_t.
* @param layer            The index of the sprite layer
* @param index            The sprite index of the sprite that you want the captured sprite for.
* @param position         This should be a pointer to a vec3 that you pass in that will get loaded with the position of the instance
*/
tfxAPI void tfx_GetCapturedInstanceTransform(tfx_stage pm, tfxU32 layer, tfxU32 index, float out_position[3]);

/*
Get the end index offset into the sprite memory for sprite data containing a pre recorded effect animation. 
a for loop to iterate over the instance_data in a pre-recorded effect
* @param sprite_data    A pointer to tfx_sprite_data_t containing all the instance_data and frame data
* @param frame            The index of the frame you want the end index for
* @param layer            The sprite layer
* @returns                tfxU32 containing the end offset
*/
tfxAPI tfxU32 tfx_SpriteDataEndIndex(tfx_sprite_data_t *sprite_data, tfxU32 frame, tfxU32 layer);

/*
Get the end index offset into the ribbon memory for ribbon data containing a pre recorded effect animation. 
a for loop to iterate over the instance_data in a pre-recorded effect
* @param sprite_data    A pointer to tfx_sprite_data_t containing all the instance_data and frame data
* @param frame            The index of the frame you want the end index for
* @returns                tfxU32 containing the end offset
*/
tfxAPI tfxU32 tfx_RibbonDataEndIndex(tfx_sprite_data_t *sprite_data, tfxU32 frame);

/*
Make a effect manager stop spawning. This will mean that all emitters in the effect manager will no longer spawn any particles so all currently running effects will expire
as the remaining particles come to the end of their life. Any single particles will also get flagged to expire
* @param pm                A pointer to a tfx_stage_t.
* @param yesno            True = disable spawning, false = enable spawning
*/
tfxAPI void tfx_DisableStageSpawning(tfx_stage pm, bool yesno);

/*
Get the buffer of effect indexes in the effect manager.
* @param pm               A pointer to a tfx_stage_t.
* @param depth            The depth of the list that you want. 0 are top level effects and anything higher are sub effects within those effects
* @param count			  A pointer to an int that you can pass in that will be filled with the count of effects in the array
* @returns                Pointer to the array of effect indexes
*/
tfxAPI tfx_effect_index_t *tfx_GetStageEffectBuffer(tfx_stage pm, int *count);

/*
Get the buffer of emitter indexes in the effect manager.
* @param pm                A pointer to a tfx_stage_t.
* @param depth            The depth of the list that you want. 0 are top level emitters and anything higher are sub emitters within those effects
* @param count			  A pointer to an int that you can pass in that will be filled with the count of emitters in the array
* @returns                Pointer to the tfxvec of effect indexes
*/
tfxAPI tfxU32 *tfx_GetStageEmitterBuffer(tfx_stage pm, int *count);

/*
Error-handling contract for the effect manager functions that take a tfxEffectID (below):
Every such function validates the id. In debug builds an invalid id trips an assert at the call site.
In release builds the call is a documented no-op instead of undefined behaviour: mutators return
without touching any state, functions that return a pointer return NULL, and functions that return a
value return 0. This means a stale effect id (for example one that expired a frame earlier than the
caller expected) can never write through a bad index and corrupt the effect manager. Effect ids are
produced by tfx_AddEffectTemplateToStage and should be treated as invalid once the effect has
expired or been removed. The same contract applies to the animation manager functions that take a
tfxAnimationID.
*/

/*
Set the position of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param x                The x value of the position
* @param y                The y value of the position
* @param z                The z value of the position
*/
tfxAPI void tfx_SetEffectPosition(tfx_stage pm, tfxEffectID effect_index, float x, float y, float z);

/*
Set the position of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param position        A float[3] array containing the x, y and z coordinates
*/
tfxAPI void tfx_SetEffectPositionVec3(tfx_stage pm, tfxEffectID effect_index, float position[3]);

/*
Move an Effect by a specified amount relative to the effect's current position
* @param pm                A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param amount            A float[3] array containing the amount to move in the x, y and z planes
*/
tfxAPI void tfx_MoveEffectVec3(tfx_stage pm, tfxEffectID effect_index, float amount[3]);

/*
Move an Effect by a specified amount relative to the effect's current position
* @param pm               A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index     The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param x                The amount to move in the x plane
* @param y                The amount to move in the y plane
* @param z                The amount to move in the z plane
*/
tfxAPI void tfx_MoveEffect(tfx_stage pm, tfxEffectID effect_index, float x, float y, float z);

/*
Get the current position of an effect
* @param pm              A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @return                float[3] array containing the effect position
*/
tfxAPI void tfx_GetEffectPositionVec3(tfx_stage pm, tfxEffectID effect_index, float out_position[3]);

/*
You can use this function to get the billboard buffer of a specific effect. 
* @param pm						A pointer to a tfx_stage_t where the effect is being managed
* @param effect_index			The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param tfxU32					Pass in a pointer to a tfxU32 which will be set to the number of instance_data in the buffer.
* @return						tfx_instance_t pointer to the buffer
*/
tfxAPI tfx_instance_t *tfx_GetEffectInstanceBuffer(tfx_stage pm, tfxEffectID effect_index, tfxU32 *sprite_count);

/*
You can use this function to get each billboard buffer for every effect that is currently active in the effect manager. Generally you would call this inside a for loop for each layer.
* @param pm						A pointer to a tfx_stage_t where the effect is being managed
* @param tfx_sprite_billboard_t	Pass in a pointer which will be set to the current sprite buffer containing all of the sprite data for this frame.
* @param tfx_effect_instance_data_t   Pass in a second pointer which will be set to the tfx_effect_instance_data_t containing all of the sprite buffer data. This can be used to gain access to all the sprite data if using double buffered instance_data (to interpolated with the previous frame).
* @param tfxU32					Pass in a pointer to a tfxU32 which will be set to the number of instance_data in the buffer.
* @return						true or false if the next billboard buffer was found. False will be returned once there are no more effect sprite buffers in the effect manager
*/
tfxAPI bool tfx_GetNextInstanceBuffer(tfx_stage pm, tfx_instance_t **sprites_soa, tfx_effect_instance_data_t **effect_sprites, tfxU32 *sprite_count);

/*After calling GetNextBillboard/SpriteBuffer in a while loop you can call this to reset the index for the next frame
* @param pm						A pointer to a tfx_stage_t
*/
tfxAPI void tfx_ResetInstanceBufferLoopIndex(tfx_stage pm);

/*
Set the roll of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current roll of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param roll            A float of the amount that you want to set the roll too
*/
tfxAPI void tfx_SetEffectRoll(tfx_stage pm, tfxEffectID effect_index, float roll);

/*
Set the pitch of a effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current pitch of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param pitch            A float of the amount that you want to set the pitch too
*/
tfxAPI void tfx_SetEffectPitch(tfx_stage pm, tfxEffectID effect_index, float pitch);

/*
Set the yaw of a effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current yaw of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param yaw            A float of the amount that you want to set the yaw too
*/
tfxAPI void tfx_SetEffectYaw(tfx_stage pm, tfxEffectID effect_index, float yaw);

/*
Set the width of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current width of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param width            A float of the amount that you want to set the width multiplier too. The width multiplier will multiply all widths of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void tfx_SetEffectWidthMultiplier(tfx_stage pm, tfxEffectID effect_index, float width);

/*
Set the height of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current height of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param height            A float of the amount that you want to set the height multiplier too. The height multiplier will multiply all heights of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void tfx_SetEffectHeightMultiplier(tfx_stage pm, tfxEffectID effect_index, float height);

/*
Set the depth of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current depth of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param depth            A float of the amount that you want to set the depth multiplier too. The depth multiplier will multiply all heights of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void tfx_SetEffectDepthMultiplier(tfx_stage pm, tfxEffectID effect_index, float depth);

/*
Set the life multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current life of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param life            A float of the amount that you want to set the life multiplier too. The life mulitplier will affect how long all particles emitted within the effect will last before expiring.
*/
tfxAPI void tfx_SetEffectLifeMultiplier(tfx_stage pm, tfxEffectID effect_index, float life);

/*
Set the particle width multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current particle width of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param width            A float of the amount that you want to set the particle width multiplier too. The particle width mulitplier will affect the width of each particle if the emitter has a non uniform particle size, otherwise
						it will uniformly size the particle
*/
tfxAPI void tfx_SetEffectParticleWidthMultiplier(tfx_stage pm, tfxEffectID effect_index, float width);

/*
Set the particle height multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current particle width of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param height            A float of the amount that you want to set the particle height multiplier too. The particle height mulitplier will affect the height of each particle if the emitter has a non uniform particle size, otherwise
						this function will have no effect.
*/
tfxAPI void tfx_SetEffectParticleHeightMultiplier(tfx_stage pm, tfxEffectID effect_index, float height);

/*
Set the velocity multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current velocity of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param velocity        A float of the amount that you want to set the particle velocity multiplier too. The particle velocity mulitplier will affect the base velocity of a particle at spawn time.
*/
tfxAPI void tfx_SetEffectVelocityMultiplier(tfx_stage pm, tfxEffectID effect_index, float velocity);

/*
Set the spin multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current spin of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param spin            A float of the amount that you want to set the particle spin multiplier too. The particle spin mulitplier will affect the base spin of a particle at spawn time.
*/
tfxAPI void tfx_SetEffectSpinMultiplier(tfx_stage pm, tfxEffectID effect_index, float spin);

/*
Set the intensity multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current intensity of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param intensity        A float of the amount that you want to set the particle intensity multiplier too. The particle intensity mulitplier will instantly affect the opacity of all particles currently emitted by the effect.
*/
tfxAPI void tfx_SetEffectIntensityMultiplier(tfx_stage pm, tfxEffectID effect_index, float intensity);

/*
Set the splatter multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current splatter of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param splatter        A float of the amount that you want to set the particle splatter multiplier too. The particle splatter mulitplier will change the amount of random offset all particles emitted in the effect will have.
*/
tfxAPI void tfx_SetEffectSplatterMultiplier(tfx_stage pm, tfxEffectID effect_index, float splatter);

/*
Set the weight multiplier of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current weight of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param weight            A float of the amount that you want to set the particle weight multiplier too. The particle weight mulitplier will change the weight applied to particles in the effect at spawn time.
*/
tfxAPI void tfx_SetEffectWeightMultiplier(tfx_stage pm, tfxEffectID effect_index, float weight);

/*
Set the overal scale of an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed. Note that this must be called after tfx_UpdateStage in order to override the current weight of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param overall_scale    A float of the amount that you want to set the overal scale to. The overal scale is an simply way to change the size of an effect
*/
tfxAPI void tfx_SetEffectOverallScale(tfx_stage pm, tfxEffectID effect_index, float overall_scale);

/*
Set the base noise offset for an effect
* @param pm                A pointer to a tfx_stage_t where the effect is being managed.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToStage
* @param noise_offset    A float of the amount that you want to set the effect noise offset to. By default when an effect is added to a effect manager a random noise offset will be set based on the Base Noise Offset Range property. Here you can override that
						value by setting it here. The most ideal time to set this would be immediately after you have added the effect to the effect manager, but you could call it any time you wanted for a constantly changing noise offset.
*/
tfxAPI void tfx_SetEffectBaseNoiseOffset(tfx_stage pm, tfxEffectID effect_index, float noise_offset);

/*
Get the name of an effect
* @param pm                A pointer to the effect
* @returns                const char * name
*/
tfxAPI const char *tfx_GetEffectName(tfx_effect_descriptor effect);


//--------------------------------
//Animation_manager
//--------------------------------

/*
Set the position of an animation instance
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @param position                A float[3] array containing the x, y and z coordinates
*/
tfxAPI void tfx_SetAnimationPosition(tfx_animation_manager animation_manager, tfxAnimationID animation_id, float position[3]);

/*
Set the scale of an animation instance
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @param scale                    A multiplier that will determine the overal size/scale of the effect
*/
tfxAPI void tfx_SetAnimationScale(tfx_animation_manager animation_manager, tfxAnimationID animation_id, float scale);

/*
Get an animation instance from an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param tfxAnimationID            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @returns pointer to instance    Pointer to a tfx_animation_instance_t
*/
tfxAPI tfx_animation_instance_t *tfx_GetAnimationInstance(tfx_animation_manager animation_manager, tfxAnimationID animation_id);

/*
Initialise an Animation Manager for use with instance data. This must be run before using an animation manager. An animation manager is used
to playback pre recorded particle effects as opposed to using a effect manager that simulates the particles in
real time. This pre-recorded data can be uploaded to the gpu for a compute shader to do all the interpolation work
to calculate the state of particles between frames for smooth animation.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param max_instances            The maximum number of animation instances that you want to be able to play at one time.
* @param initial_capacity        Optionally, you can set an initial capacity for the sprite data. The data will grow if you add
								beyond this amount but it gives you a chance to reserve a decent amount to start with to
								save too much mem copies as the data grows
*/
tfxAPI tfx_animation_manager tfx_CreateAnimationManager(tfxU32 max_instances, tfxU32 initial_sprite_data_capacity);

/*
Set the callback that you can use to determine whether or not a tfx_animation_instance_t should be added to the next frame's render queue. You can use this
to cull instances that are outside of the view frustum for example
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param callback                Pointer to the callback you want to use. It must have the following signature:
								bool(*maybe_render_instance_callback(tfx_animation_manager animation_manager, tfx_animation_instance_t *instance, tfx_frame_meta_t *meta, void *user_data))
								Values passed into the callback function are a pointer to the animation manager, a pointer to the instance being processed, a pointer to
								the frame meta of the instance, this will contain the bounding box and radius of the instance from the current frame of the instance and a pointer
								to any user data that you set that might contain the camera frustum that you want to check against.
*/
tfxAPI void tfx_SetAnimationManagerInstanceCallback(tfx_animation_manager animation_manager, tfx_maybe_render_instance_callback maybe_render_instance_callback);

/*
Get the sprite data settings for an effect in a library. Sprite data settings are the settings for an effect in the editor relating to setting up pre-baked effects
* @param library				Pointer to the tfx_library_t where the effect is stored
* @param effect					Pointer the the effect that you want the sprite settings for.
* @returns						Pointer to the tfx_sprite_data_settings
*/
tfxAPI tfx_sprite_data_settings_t *tfx_GetEffectSpriteDataSettings(tfx_library library, tfx_effect_descriptor effect);

/*
Get the sprite data settings for an effect in a library by it's path. Sprite data settings are the settings for an effect in the editor relating to setting up pre-baked effects
* @param library				Pointer to the tfx_library_t where the effect is stored
* @param path					const char* string of the path to the effect. Must be the path to a root effect.
* @returns						Pointer to the tfx_sprite_data_settings
*/
tfxAPI tfx_sprite_data_settings_t *tfx_GetEffectSpriteDataSettingsByPath(tfx_library library, const char *path);

/*
Get the index offset into the sprite memory for sprite data containing a pre recorded effect animation. Can be used along side tfx_SpriteDataEndIndex to create
a for loop to iterate over the instance_data in a pre-recorded effect
* @param sprite_data    A pointer to tfx_sprite_data_t containing all the instance_data and frame data
* @param frame            The index of the frame you want the offset for
* @param layer            The sprite layer
* @returns                tfxU32 containing the index offset
*/
tfxAPI tfxU32 tfx_SpriteDataIndexOffset(tfx_sprite_data_t *sprite_data, tfxU32 frame, tfxU32 layer);

/*
Get the index offset into the ribbon memory for ribbon data containing a pre recorded effect animation. Can be used along side tfx_RibbonDataEndIndex to create
a for loop to iterate over the instance_data in a pre-recorded effect
* @param sprite_data      A pointer to tfx_sprite_data_t containing all the instance_data and frame data
* @param frame            The index of the frame you want the offset for
* @returns                tfxU32 containing the index offset
*/
tfxAPI tfxU32 tfx_RibbonDataIndexOffset(tfx_sprite_data_t *sprite_data, tfxU32 frame);

/*
Set the user data in a tfx_animation_manager_t which can get passed through to callback functions when updated the animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param user_data                void* pointer to the data that you want to set
*/
tfxAPI void tfx_SetAnimationManagerUserData(tfx_animation_manager animation_manager, void *user_data);

/*
Add sprite data to an animation manager sprite data buffer from an effect. This will record the
animation if necessary and then convert the sprite data to tfx_sprite_data_t ready for uploading
to the GPU
* @param animation_manager       A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @param position                A float[3] array containing the x, y and z coordinates
* @param sprite_data             Optional advanced usage if you want to use sprite data you recorded in another process. Pass 0 to just use the record or lookup pre-recorded data for the effect.
*/
tfxAPI void tfx_AddSpriteData(tfx_animation_manager animation_manager, tfx_effect_descriptor effect, tfx_stage pm, float camera_position[3], tfx_sprite_data_t *sprite_data);

/*
Add an animation instance to the animation manager.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param path                    tfxKey path hash of the effect name and path: effect.path_hash
* @param start_frame            Starting frame of the animation
* @returns                        The index id of the animation instance. You can use this to reference the animation when changing position, scale etc
								Return tfxINVALID if there is no room in the animation manager
*/
tfxAPI tfxAnimationID tfx_AddAnimationInstanceByKey(tfx_animation_manager animation_manager, tfxKey path, tfxU32 start_frame);

/*
Add an animation instance to the animation manager.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param path                    const char * name of the effect. If the effect was in a folder then specify the whole path
* @param start_frame            Starting frame of the animation
* @returns                        The index id of the animation instance. You can use this to reference the animation when changing position, scale etc
								Return tfxINVALID if there is no room in the animation manager
*/
tfxAPI tfxAnimationID tfx_AddAnimationInstance(tfx_animation_manager animation_manager, const char *path, tfxU32 start_frame);

/*
Update an animation manager to advance the time and frames of all instances currently playing.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to update
* @param start_frame            Starting frame of the animation
*/
tfxAPI void tfx_UpdateAnimationManager(tfx_animation_manager animation_manager, float elapsed);

/*
Add an effect's shapes to an animation manager. You can use this function if you're manually recording particle effects and adding them to an animation
manager rather then just using the editor.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to update
* @param effect                    A pointer to the effect whose shapes you want to add
*/
tfxAPI void tfx_AddEffectShapes(tfx_animation_manager animation_manager, tfx_effect_descriptor effect);

/*
Update an animation manager so that the effects do not expire they just loop forever instead regardless of whether they're a looped effect or not.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to update
*/
tfxAPI void tfx_CycleAnimationManager(tfx_animation_manager animation_manager);

/*
Clears all animation instances currently in play in an animation manager, resulting in all currently running animations
from being drawn
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to clear
*/
tfxAPI void tfx_ClearAllAnimationInstances(tfx_animation_manager animation_manager);

/*
Clears all data from the animation manager including sprite data, metrics and instances. Essentially resetting everything back to
it's initialisation point
from being drawn
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to reset
*/
tfxAPI void tfx_ResetAnimationManager(tfx_animation_manager animation_manager);

/*
Frees all data from the animation manager including sprite data, metrics and instances and also the handle itself.
You don't have to call this - tfx_EndTimelineFX frees any animation managers that are still alive - but you can
free them earlier if you want to reclaim the memory. Calling it is safe either way, the animation manager is
deregistered from the shutdown sweep so it won't be freed twice.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to free
*/
tfxAPI void tfx_FreeAnimationManager(tfx_animation_manager animation_manager);

/*
Get the tfx_animation_buffer_metrics_t from an animation manager. This will contain the info you need to upload the sprite data,
offsets and animation instances to the GPU. Only offsets and animation instances need to be uploaded to the GPU each frame. Sprite
data can be done ahead of time.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @returns                        tfx_animation_buffer_metrics_t containing buffer sizes
*/
tfxAPI tfx_animation_buffer_metrics_t tfx_GetAnimationBufferMetrics(tfx_animation_manager animation_manager);

/*
Get the total number of instance_data that need to be drawn by an animation manager this frame. You can use this in your renderer
to draw your sprite instances
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @returns                        tfxU32 of the number of instance_data
*/
tfxAPI tfxU32 tfx_GetTotalSpritesThatNeedDrawing(tfx_animation_manager animation_manager);

/*
Get the total number of ribbons that need to be drawn by an animation manager this frame. You can use this in your renderer
to draw your sprite instances
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @returns                        tfxU32 of the number of ribbons to draw
*/
tfxAPI tfxU32 tfx_GetTotalRibbonsThatNeedDrawing(tfx_animation_manager animation_manager);

/*
Get the total number of instances being processed by an animation manager. This will not necessarily be the same number as
the instances being rendered if some are being culled in your custom callback if your using one.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to clear
* @returns int                    The number of instances being updated
*/
tfxAPI tfxU32 tfx_GetTotalInstancesBeingUpdated(tfx_animation_manager animation_manager);

/*
Create the image data required for shaders from a TimelineFX library. The image data will contain data such as uv coordinates. Once you have built the data you can use GetLibraryImageData to get the buffer
and upload it to the gpu.
* @param animation_manager		  A pointer to an tfx_animation_manager_t object
* @param shapes                   A pointer to a tfx_gpu_shapes_t object which will fill a buffer with all the shapes
* @param uv_lookup                A function pointer to a function that you need to set up in order to get the uv coordinates from whatever renderer you're using
*/
tfxAPI void tfx_BuildAnimationManagerGPUShapeData(tfx_animation_manager animation_manager, tfx_gpu_shapes shapes, tfx_uv_lookup uv_lookup);

/*
Get a pointer to the particle shapes data in the animation manager. This can be used with tfx_BuildGPUShapeData when you want to upload the data to the GPU
* @param animation_manager        A pointer the tfx_animation_manager_t
*/
tfxAPI tfx_image_data_t *tfx_GetParticleShapesAnimationManager(tfx_animation_manager animation_manager, int *count);

/*
Get the total number of instance_data in an animation manger's sprite data buffer
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32                The number of instance_data in the buffer
*/
tfxAPI tfxU32 tfx_GetTotalSpriteDataCount(tfx_animation_manager animation_manager);

/*
Get the total byte size of instance_data in an animation manger's sprite data buffer
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                 The number of instance_data in the buffer
*/
tfxAPI size_t tfx_GetSpriteDataSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the total byte size of ribbon_data in an animation manger's ribbon data buffer
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32                The number of instance_data in the buffer
*/
tfxAPI size_t tfx_GetRibbonDataSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the buffer memory address for the sprite data in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                A pointer to the sprite data memory
*/
tfxAPI void *tfx_GetSpriteDataBufferPointer(tfx_animation_manager animation_manager);

/*
Get the buffer memory address for the ribbon data in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                  A pointer to the sprite data memory
*/
tfxAPI void *tfx_GetRibbonDataBufferPointer(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the offsets buffer in an animation manager for sprite data
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                Size in bytes of the offsets buffer
*/
tfxAPI size_t tfx_GetOffsetsSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the buffer pointer for the offsets buffer in an animation manager for sprite data
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                Size in bytes of the offsets buffer
*/
tfxAPI void *tfx_GetOffsetsBufferPointer(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the offsets buffer in an animation manager for ribbons
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns size_t                Size in bytes of the offsets buffer
*/
tfxAPI size_t tfx_GetRibbonOffsetsSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the buffer pointer for the offsets buffer in an animation manager for ribbon data
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns size_t                Size in bytes of the offsets buffer
*/
tfxAPI void *tfx_GetRibbonOffsetsBufferPointer(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the render queue of animation instances buffer in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                Size in bytes of the instances buffer
*/
tfxAPI size_t tfx_GetAnimationInstancesSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the count of of animation instances in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32                 Count of the instances buffer
*/
tfxAPI tfxU32 tfx_GetAnimationInstancesCount(tfx_animation_manager animation_manager);

/*
Get the buffer pointer of of animation instances in an animation manager which you'll need to copy
the data to the GPU
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                  Pointer to the instances buffer
*/
tfxAPI void *tfx_GetAnimationInstancesBufferPointer(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the animation emitter properties list
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                Size in bytes of the properties bufffer
*/
tfxAPI size_t tfx_GetAnimationEmitterPropertySizeInBytes(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the animation ribbon properties list
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns size_t                Size in bytes of the properties bufffer
*/
tfxAPI size_t tfx_GetAnimationRibbonPropertySizeInBytes(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the animation ribbon segments list
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon segment data from
* @returns size_t                Size in bytes of the ribbon segments bufffer
*/
tfxAPI size_t tfx_GetAnimationRibbonSegmentsSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the size in bytes of the animation ribbon data list
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns size_t                Size in bytes of the ribbon segments bufffer
*/
tfxAPI size_t tfx_GetAnimationRibbonDataSizeInBytes(tfx_animation_manager animation_manager);

/*
Get the number of emitter properties being using by the animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32                Number of emitter properties
*/
tfxAPI tfxU32 tfx_GetAnimationEmitterPropertyCount(tfx_animation_manager animation_manager);

/*
Get the buffer memory address for the sprite emitter properties in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                A pointer to the sprite emitter properties data memory
*/
tfxAPI void *tfx_GetAnimationEmitterPropertiesBufferPointer(tfx_animation_manager animation_manager);

/*
Get the buffer memory address for the ribbon property data in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                A pointer to the ribbon properties data memory
*/
tfxAPI void *tfx_GetAnimationRibbonPropertiesBufferPointer(tfx_animation_manager animation_manager);

/*
Get the buffer memory address for the ribbon segments data in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                  A pointer to the ribbon segments data memory
*/
tfxAPI void *tfx_GetAnimationRibbonSegmentsBufferPointer(tfx_animation_manager animation_manager);

/*
Get the buffer memory address for the ribbon segments data in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                  A pointer to the ribbon segments data memory
*/
tfxAPI void *tfx_GetAnimationRibbonDataBufferPointer(tfx_animation_manager animation_manager);

/*
Returns true or false if the animation manager contains effects with ribbons
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns bool 	              True if the animation manager has ribbons
*/
tfxAPI bool tfx_AnimationManagerHasRibbons(tfx_animation_manager animation_manager);

/*
Returns true or false if the animation manager contains emitters that use animated shapes
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @returns bool 	              True if the animation manager does contain animated shapes
*/
tfxAPI bool tfx_HasAnimatedShapes(tfx_animation_manager animation_manager);

/*
Set the callback used by the animation manager to determin if an animation should be rendered or not
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the ribbon data from
* @param callback 	              The function callback
*/
tfxAPI void tfx_SetAnimationManagerCullCallback(tfx_animation_manager animation_manager, tfx_maybe_render_instance_callback callback);

tfxAPI size_t tfx_CalculateAnimationInstanceBufferSize(size_t instance_count);
tfxAPI size_t tfx_CalculateAnimationOffsetsBufferSize(size_t instance_count);

//--------------------------------
//Effect_templates
//--------------------------------
 
/*
Prepare a tfx_effect_template_t that you can use to customise effects in the library in various ways before adding them into a effect manager for updating and rendering. Using a template like this
means that you can tweak an effect without editing the base effect in the library.
* @param library                    A reference to a tfx_library_t that should be loaded with tfx_LoadEffectLibrary
* @param name                       The name of the effect in the library that you want to use for the template. If the effect is in a folder then use normal pathing: "My Folder/My effect"
//Returns handle					Handle to the newly created effect template or nullptr if the effect couldn't be found in the library
*/
tfxAPI tfx_effect_template tfx_CreateEffectTemplate(tfx_library library, const char *name);
 
/*
Delete an effect template and free all memory associated with it
* @param effect_template            A handle to the effect template to be deleted
//Returns handle					Handle to the newly created effect template or nullptr if the effect couldn't be found in the library
*/
tfxAPI void tfx_FreeEffectTemplate(tfx_effect_template effect_template);

/*
Reset an effect template and make it empty so you can use it to store another effect.
* @param t                        A pointer to a tfx_effect_template_t
*/
tfxAPI void tfx_ResetTemplate(tfx_effect_template t);

/*
Get the root effect from the template
* @param t                        A pointer to a tfx_effect_template_t
* @returns                        A pointer to the root effect
*/
tfxAPI tfx_effect_descriptor tfx_GetEffectFromTemplate(tfx_effect_template t);

/*
Get an emitter or sub effect from an effect template.
* @param t                        A pointer to a tfx_effect_template_t
* @param path                     A path to the emitter or sub effect that you want to retrieve. Must be a valid path. Example path might be: "Explosion/Smoke"
* @returns                        A pointer to the root effect
*/
tfxAPI tfx_effect_descriptor tfx_GetEmitterFromTemplate(tfx_effect_template t, const char *path);

/*
Get an emitter path that an emitter is using. The emitter must have the path emission type set or nullptr will be returned
* @param t                        A pointer to a tfx_effect_descriptor_t
* @param path                     A path to the emitter or sub effect that you want to retrieve. Must be a valid path. Example path might be: "Explosion/Smoke"
* @returns                        A pointer to the root effect
*/
tfxAPI tfx_emitter_path_t *tfx_GetEmitterPath(tfx_effect_descriptor e);

/*
Set the user data for any effect or emitter in the effect template. This user data will get passed through to any update callback functions
* @param t                        A pointer to a tfx_effect_template_t
* @param path                     A path to the effect or emitter in the effect template
* @param data                     A pointer to the user data
*/
tfxAPI void tfx_SetTemplateUserData(tfx_effect_template t, const char *path, void *data);

/*
Set the user data for the root effect in an effect template
* @param t                        A pointer to a tfx_effect_template_t
* @param data                     A pointer to the user data
*/
tfxAPI void tfx_SetTemplateEffectUserData(tfx_effect_template t, void *data);

/*
Set the same user data for all effects and emitters/sub effects in the effect template
* @param t                        A pointer to a tfx_effect_template_t
* @param data                     A pointer to the user data that will be set to all effects and emitters in the template
*/
tfxAPI void tfx_SetTemplateUserDataAll(tfx_effect_template t, void *data);

/*
Set an update callback for the root effect in the effect template.
* @param t                        A pointer to a tfx_effect_template_t
* @param update_callback          A pointer to the call back function
*/
tfxAPI void tfx_SetTemplateEffectUpdateCallback(tfx_effect_template t, void(*update_callback)(tfx_stage pm, tfxEffectID effect_index));

/*
Pre-record this effect into a sprite cache so that you can play the effect back without the need to actually caclulate particles in realtime.
	* @param pm					  Reference to a pm that will be used to run the particle simulation and record the sprite data
	* @param path				  const *char of a path to the emitter in the effect.Must be a valid path, for example: "My Effect/My Emitter"
	* @param camera				  Array of 3 floats with the camera position (only needed for effects that are sorted by depth)
*/
tfxAPI void tfx_RecordTemplateEffect(tfx_effect_template t, tfx_stage pm, float update_frequency, float camera_position[3]);

/*
Pre-record this effect into a sprite cache so that you can play the effect back without the need to actually caclulate particles in realtime. This version
of the function allows you to pass in specific settings
	* @param pm					  Reference to a pm that will be used to run the particle simulation and record the sprite data
	* @param path				  const *char of a path to the emitter in the effect.Must be a valid path, for example: "My Effect/My Emitter"
	* @param camera				  Array of 3 floats with the camera position (only needed for effects that are sorted by depth)
*/
tfxAPI void tfx_RecordEffect(tfx_effect_descriptor e, tfx_sprite_data_settings_t *settings, tfx_sprite_data_t *sprite_data, tfx_stage pm, float update_frequency, float camera_position[3]);

/*
Disable an emitter within an effect. Disabling an emitter will stop it being added to the effect manager when calling tfx_AddEffectTemplateToStage
* @param path					  const *char of a path to the emitter in the effect. Must be a valid path, for example: "My Effect/My Emitter"
*/
tfxAPI void tfx_DisableTemplateEmitter(tfx_effect_template t, const char *path);

/*
Enable an emitter within an effect so that it is added to the effect manager when calling tfx_AddEffectTemplateToStage. Emitters are enabled by default.
* @param path					  const *char of a path to the emitter in the effect. Must be a valid path, for example: "My Effect/My Emitter"
*/
tfxAPI void tfx_EnableTemplateEmitter(tfx_effect_template t, const char *path);

/*
Scale all nodes on a global graph graph of the effect
* @param tfx_effect_template	  Handle to a valid effect template
* @param graph_index			  tfx_global_graph_index enum of the global graph that you want to scale. 
* @param amount					  A float of the amount that you want to scale the multiplier by.
*/
tfxAPI void tfx_ScaleTemplateGlobalMultiplier(tfx_effect_template t, tfx_global_graph_index graph_index, float scale_amount);

/*
Set the single spawn amount for an emitter. Only affects emitters that have the single spawn flag set.
* @param emitter_path			 const *char of the emitter path
* @param amount					 A float of the amount that you want to set the single spawn amount to.
*/
tfxAPI void tfx_SetTemplateSingleSpawnAmount(tfx_effect_template t, const char *emitter_path, tfxU32 amount);

/*
Scale all nodes on an emitter graph
* @param tfx_effect_template	  Handle to a valid effect template
* @param emitter_path			  const *char of the emitter path
* @param global_type			  tfx_graph_type of the emitter graph that you want to scale. Must be an emitter graph or an assert will be called
* @param scale amount             A float of the amount that you want to scale the graph by.
*/
tfxAPI void tfx_ScaleTemplateEmitterGraph(tfx_effect_template t, const char *emitter_path, tfx_emitter_graph_index graph_index, float scale_amount);

//--------------------------------
//Editing_graphs
//--------------------------------
tfxAPI void tfx_ClearBaseLifetimeGraph(tfx_effect_descriptor emitter, float v );
tfxAPI void tfx_ClearVariationLifetimeGraph(tfx_effect_descriptor emitter, float v);

//--------------------------------
//General_helpers
//--------------------------------

tfxAPI void tfx_GetSpriteScale(void *instance, float out_scale[2]);

tfxAPI inline float tfx_GetDistance(float fromx, float fromy, float tox, float toy) {
	float w = tox - fromx;
	float h = toy - fromy;
	return sqrtf(w * w + h * h);
}

tfxAPI tfx_effect_descriptor tfx_CreateEffectDescriptor(tfx_effect_descriptor_type type);

/*
Create a new list for containing gpu shapes. You can use this list to upload to the GPU so that shaders can have access to particle image data like UV coordinates
* @returns tfx_gpu_shapes                A handle to the shapes list
*/
tfxAPI tfx_gpu_shapes tfx_CreateGPUShapesList(void);

/*
Delete and free all the memory for a tfx_gpu_shapes list.
* @param shapes                 A tfx_gpu_shapes handle
*/
tfxAPI void tfx_FreeGPUShapesList(tfx_gpu_shapes shapes);

/*
Clear the tfx_gpu_shapes list but keep the memory associated with it.
* @param shapes                 A tfx_gpu_shapes handle
*/
tfxAPI void tfx_ClearGPUShapesList(tfx_gpu_shapes shapes);

/*
Get a pointer to the GPU shapes array which you can use in a memcpy to a staging buffer for uploading to a GPU
* @param particle_shapes        A pointer the tfx_gpu_shapes_t
*/
tfxAPI tfx_gpu_image_data_t *tfx_GetGPUShapesArray(tfx_gpu_shapes shapes);

/*
Get the number of shapes in the GPU Shape Data buffer contained within a library. Make sure you call tfx_BuildGPUShapeData first or they'll be nothing to return
* @param shapes                 A tfx_gpu_shapes handle
* @returns tfxU32               The number of shapes in the buffer
*/
tfxAPI tfxU32 tfx_GetGPUShapesCount(tfx_gpu_shapes shapes);

/*
Get the size in bytes of shapes from a tfx_gpu_shapes handle. You can use this when uploading the shape data to the GPU along with tfx_GetGPUShapesArray
* @param shapes                 A tfx_gpu_shapes handle
* @returns tfxU32               The number of shapes in the buffer
*/
tfxAPI size_t tfx_GetGPUShapesSizeInBytes(tfx_gpu_shapes shapes);

/*
Add a new tfx_gpu_image_data_t object onto a list of gpu shapes
* @param shapes                 A tfx_gpu_shapes handle
* @param image data             The tfx_gpu_image_data_t object that you want to add to the shapes list
* @returns index				An index into the array where the new image data was added
*/
tfxAPI tfxU32 tfx_AddGPUShape(tfx_gpu_shapes shapes, tfx_gpu_image_data_t image_data);

/*
Retrieve image data point from a tfx_gpu_shapes handle containing a list of tfx_gpu_image_data_t objects
* @param shapes                 A tfx_gpu_shapes handle
* @param index					
* @returns index				An index into the array where the new image data was added
*/
tfxAPI tfx_gpu_image_data_t *tfx_GetGPUShape(tfx_gpu_shapes shapes, tfxU32 index);

tfxAPI tfx_version_t tfx_GetVersion(void);

// Helper function to get a good default thread count for thread pools
// Usually hardware threads - 1 to leave a core for the OS/main thread
tfxAPI unsigned int tfx_GetDefaultThreadCount(void);

#endif
