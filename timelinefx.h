#ifndef TFX_LIBRARY_HEADER
#define TFX_LIBRARY_HEADER

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

#define tfxENABLE_PROFILING
#define tfxPROFILER_SAMPLES 60
#define TFX_THREAD_SAFE
#define TFX_EXTRA_DEBUGGING
#define SSE41		//Steam survey current has this at 99.83% coverage 12 April 2025. I will probably make this the minimum requirement

//Enable this to process 8 particles at a time.
//#define tfxUSEAVX

//Enable fused multiply add in simd calculations
//#define tfxUSEFMA

//#define TFX_MEMORY_TRACKING

/*
	Timeline FX C++ library

	This library is for implementing particle effects into your games and applications.

	This library is render agnostic, so you will have to provide your own means to render the particles. There are various API functions in the library that help you do this.

	Currently tested on Windows and MacOS, Intel and Mac based ARM processors.

	Table of contents
	Sections in this header file, you can search for the following keywords to jump to that section:

	[Zest_Pocket_Allocator]				A single header library for allocating memory from a large pool.
	[Header_Includes_and_Typedefs]		Just your basic header stuff for setting up typedefs and some #defines
	[Macros]							Macro definitions
	[OS_Specific_Functions]				OS specific multithreading and file access
	[Pocket_Hasher]						XXHasher for the storage map.
	[SIMD_defines]						Defines for SIMD intrinsics
	[Enums]								All the definitions for enums and bit flags
		-[Graph_types]				
		-[Bit_fields]				
	[Constants]							Various constant definitions
	[Preset_spline paths] 	 			Arrays of preset path nodes
	[String_Buffers]					Basic string buffers for storing names of things in the library and reading from library files.
	[Containers_and_Memory]				Container structs and lists and defines for memory is allocated (uses Zest Pocket Allocator by default)
	[Multithreading_Work_Queues]		Implementation for work queues for threading
	[Vector_Math]						Vec2/3/4 and Matrix2/3/4 structs including wide vectors for SIMD
	[Simplex_Noise]						Some setup for implementing simplex noise.
	[Profiling]							Very basic profiling for internal use
	[File_IO]							A package manager for reading/writing files such as a tfx library effects file
	[Struct_Types]						All of the structs used for objects in TimelineFX
	[Internal_Functions]				Mainly internal functions called only by the library but also the Editor, these are marked either tfxINTERNAL or tfxAPI_EDITOR
	[Control_Position_Policies]			Inline functions for templated functions that control particle positions
	[API_Functions]						The main functions for use by users of the library
		-[Initialisation_functions]		Startup and shutdown timelinefx
		-[Global_variable_access]		Any functions that give you access to global variables relating to timelinefx
		-[Library_functions]			Functions for loading and accessing timelinefx libraries
		-[Particle_Manager_functions]	Create and update functions for effect managers where the main work is done to update particles every frame
		-[Animation_manager]			Animation manager functions for playing pre-recorded effect data
		-[Effect_templates]				Functions for working with effect templates which help modify effects in the library without actually changing the base effect in the library.
		-[General_helpers]				General math functions and other helpers.
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
#define tfxINTERNAL static    
#define tfxAPI_EDITOR 

//Override this if you'd prefer a different way to allocate the pools for sub allocation in host memory.
#ifndef tfxALLOCATE_POOL
#define tfxALLOCATE_POOL malloc
#endif

#ifndef tfxMAX_MEMORY_POOLS
#define tfxMAX_MEMORY_POOLS 32
#endif

#ifndef tfxMAX_THREADS
#define tfxMAX_THREADS 64
#endif

//---------------------------------------
/*  Zest_Pocket_Allocator, a Two Level Segregated Fit memory allocator
	This is my own memory allocator from https://github.com/peterigz/zloc
	This is used in TimelineFX to manage memory allocation. A large pool is created and allocated from. New pools are created if it runs out of space
	(and you initialised TimelineFX to do so).
*/
//---------------------------------------
#include <assert.h>
#include <stdlib.h>
#include <limits.h>
#include <stddef.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#if defined(_WIN32)
#include <process.h>
#include <SDKDDKVer.h>
#ifndef WIN_LEAN_AND_MEAN
#define WIN_LEAN_AND_MEAN
#endif
#include <Windows.h>
#else
#include <pthread.h>
#endif

#define tfx__Min(a, b) (((a) < (b)) ? (a) : (b))
#define tfx__Max(a, b) (((a) > (b)) ? (a) : (b))
#define tfx__Clamp(min, max, v) (v < min) ? min : (v > max) ? max : v

typedef int tfx_index;
typedef unsigned int tfx_sl_bitmap;
typedef unsigned int tfx_uint;
typedef unsigned char tfx_byte;
typedef void *tfx_pool;

#if !defined (TFX_ASSERT)
#define TFX_ASSERT assert
#endif

#define TFX_DEPRECATED assert(0 && "Function is deprecated");

#define TFX_ASSERT_INIT(magic) TFX_ASSERT(magic == tfxINIT_MAGIC)
#define TFX_ASSERT_UNINIT(magic) TFX_ASSERT(magic != tfxINIT_MAGIC)
#define TFX_ASSERT_HANDLE(handle) TFX_ASSERT(handle && *((tfxU32*)handle) == tfxINIT_MAGIC)
#define TFX_VALID_HANDLE(handle) (handle && *((tfxU32*)handle) == tfxINIT_MAGIC)

#define tfx__is_pow2(x) ((x) && !((x) & ((x) - 1)))
#define tfx__glue2(x, y) x ## y
#define tfx__glue(x, y) tfx__glue2(x, y)
#define tfx__static_assert(exp) \
    typedef char tfx__glue(static_assert, __LINE__) [(exp) ? 1 : -1]

#if (defined(_MSC_VER) && defined(_M_X64)) || defined(__x86_64__)
#define tfx__64BIT
typedef size_t tfx_size;
typedef size_t tfx_fl_bitmap;
typedef int tfxLONG;
#define TFX_ONE 1ULL
#else
typedef size_t tfx_size;
typedef size_t tfx_fl_bitmap;
typedef int tfxLONG;
#define TFX_ONE 1U
#endif

#ifndef MEMORY_ALIGNMENT_LOG2
#if defined(tfx__64BIT)
#define MEMORY_ALIGNMENT_LOG2 3        //64 bit
#else
#define MEMORY_ALIGNMENT_LOG2 2        //32 bit
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
#define TFX_PRINT_NOTICE(message_f, ...) printf(message_f"\033[0m", __VA_ARGS__)
#else
#define TFX_PRINT_NOTICE(message_f, ...)
#endif

#ifdef TFX_OUTPUT_ERROR_MESSAGES
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

typedef struct tfx_sync_s {
#ifdef _WIN32
	CRITICAL_SECTION mutex;
	CONDITION_VARIABLE empty_condition;
	CONDITION_VARIABLE full_condition;
#else
	pthread_mutex_t mutex;
	pthread_cond_t empty_condition;
	pthread_cond_t full_condition;
#endif
} tfx_sync_t;

#ifdef __cplusplus
extern "C" {
#endif

// Initialize synchronization primitives
tfxINTERNAL inline void tfx__sync_init(tfx_sync_t *sync) {
#ifdef _WIN32
	InitializeCriticalSection(&sync->mutex);
	InitializeConditionVariable(&sync->empty_condition);
	InitializeConditionVariable(&sync->full_condition);
#else
	pthread_mutex_init(&sync->mutex, NULL);
	pthread_cond_init(&sync->empty_condition, NULL);
	pthread_cond_init(&sync->full_condition, NULL);
#endif
}

tfxINTERNAL inline void tfx__sync_cleanup(tfx_sync_t *sync) {
#ifdef _WIN32
	DeleteCriticalSection(&sync->mutex);
#else
	pthread_mutex_destroy(&sync->mutex);
	pthread_cond_destroy(&sync->empty_condition);
	pthread_cond_destroy(&sync->full_condition);
#endif
}

tfxINTERNAL inline void tfx__sync_lock(tfx_sync_t *sync) {
#ifdef _WIN32
	EnterCriticalSection(&sync->mutex);
#else
	pthread_mutex_lock(&sync->mutex);
#endif
}

tfxINTERNAL inline void tfx__sync_unlock(tfx_sync_t *sync) {
#ifdef _WIN32
	LeaveCriticalSection(&sync->mutex);
#else
	pthread_mutex_unlock(&sync->mutex);
#endif
}

tfxINTERNAL inline void tfx__sync_wait_empty(tfx_sync_t *sync) {
#ifdef _WIN32
	SleepConditionVariableCS(&sync->empty_condition, &sync->mutex, INFINITE);
#else
	pthread_cond_wait(&sync->empty_condition, &sync->mutex);
#endif
}

tfxINTERNAL inline void tfx__sync_wait_full(tfx_sync_t *sync) {
#ifdef _WIN32
	SleepConditionVariableCS(&sync->full_condition, &sync->mutex, INFINITE);
#else
	pthread_cond_wait(&sync->full_condition, &sync->mutex);
#endif
}

tfxINTERNAL inline void tfx__sync_signal_empty(tfx_sync_t *sync) {
#ifdef _WIN32
	WakeConditionVariable(&sync->empty_condition);
#else
	pthread_cond_signal(&sync->empty_condition);
#endif
}

	tfxINTERNAL inline void tfx__sync_signal_full(tfx_sync_t *sync) {
#ifdef _WIN32
		WakeConditionVariable(&sync->full_condition);
#else
		pthread_cond_signal(&sync->full_condition);
#endif
	}

#define tfx__MAXIMUM_BLOCK_SIZE (TFX_ONE << TFX_MAX_SIZE_INDEX)

	//Note: putting this hear as a reminder for future memory tracking/features
	typedef struct tfx_header_meta_s {
		struct tfx_header *prev_physical_block;
#if defined(TFX_MEMORY_TRACKING)
		tfx_size memory_block_type;
#endif
		/*    Note that the size is either 4 or 8 bytes aligned so the boundary tag (2 flags denoting
			whether this or the previous block is free) can be stored in the first 2 least
			significant bits    */
		tfx_size size;
	} tfx_header_meta_t;

	enum tfx__constants {
		tfx__MEMORY_ALIGNMENT = 1 << MEMORY_ALIGNMENT_LOG2,
		tfx__SECOND_LEVEL_INDEX_LOG2 = 5,
		tfx__FIRST_LEVEL_INDEX_COUNT = TFX_MAX_SIZE_INDEX,
		tfx__SECOND_LEVEL_INDEX_COUNT = 1 << tfx__SECOND_LEVEL_INDEX_LOG2,
		tfx__BLOCK_POINTER_OFFSET = sizeof(void *) + sizeof(tfx_size) + sizeof(tfx_size),
		tfx__MINIMUM_BLOCK_SIZE = 16,
		tfx__BLOCK_SIZE_OVERHEAD = sizeof(tfx_size) + sizeof(tfx_size),
		tfx__POINTER_SIZE = sizeof(void *),
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
		Each block has a header that is used only has a pointer to the previous physical block
		and the size. If the block is free then the prev and next free blocks are also stored.
	*/
	typedef struct tfx_header {
		struct tfx_header *prev_physical_block;
#if defined(TFX_MEMORY_TRACKING)
		tfx_size memory_block_type;
#endif
		/*    Note that the size is either 4 or 8 bytes aligned so the boundary tag (2 flags denoting
			whether this or the previous block is free) can be stored in the first 2 least
			significant bits    */
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
		/*    This is basically a terminator block that free blocks can point to if they're at the end
			of a free list. */
		tfx_header null_block;
		//todo: just make thead safe by default and remove these conditions
#if defined(TFX_THREAD_SAFE)
		/* Multithreading protection*/
		tfx_sync_t mutex;
#endif
		tfx_size minimum_allocation_size;
		/*    Here we store all of the free block data. first_level_bitmap is either a 32bit int
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
	static inline tfxLONG tfx__compare_and_exchange(volatile tfxLONG *target, tfxLONG value, tfxLONG original) {
		return InterlockedCompareExchange((volatile LONG *)target, value, original);
	}

	static inline tfxLONG tfx__exchange(volatile tfxLONG *target, tfxLONG value) {
		return InterlockedExchange((volatile LONG *)target, value);
	}

	static inline unsigned int tfx__increment(unsigned int volatile *target) {
		return InterlockedIncrement(target);
	}
#endif

#define tfx__strlen strnlen_s
#define tfx__writebarrier _WriteBarrier();
#define tfx__readbarrier _ReadBarrier();
#define tfx__strcpy(dst, size, src) strcpy_s(dst, size, src)
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

	static inline tfxLONG tfx__compare_and_exchange(volatile tfxLONG *target, tfxLONG value, tfxLONG original) {
		return __sync_val_compare_and_swap(target, original, value);
	}

	static inline tfxLONG tfx__exchange(volatile tfxLONG *target, tfxLONG value) {
		return __sync_lock_test_and_set(target, value);
	}

	static inline unsigned int tfx__increment(unsigned int volatile *target) {
		return __sync_add_and_fetch(target, 1);
	}

#define tfx__strlen strnlen
#define tfx__writebarrier __asm__ __volatile__ ("" : : : "memory");
#define tfx__readbarrier __asm__ __volatile__ ("" : : : "memory");
#define tfx__strcpy(dst, size, src) strcpy(dst, src)
#define tfx__fseek fseeko
#define tfx__ftell ftello
#define TFX_ALIGN_AFFIX(v)            __attribute__((aligned(v)))
#define TFX_PACKED_STRUCT            __attribute__((packed))

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
	bool tfx_RemovePool(tfx_allocator *allocator, tfx_pool *pool);

	/*
	When using an allocator for managing remote memory, you need to set the bytes per block that a block storing infomation about the remote
	memory allocation will manage. For example you might set the value to 1MB so if you were to then allocate 4MB of remote memory then 4 blocks
	worth of space would be used to allocate that memory. This means that if it were to be freed and then split down to a smaller size they'd be
	enough blocks worth of space to do this.

	Note that the lower the number the more memory you need to track remote memory blocks but the more granular it will be. It will depend alot
	on the size of allocations you will need
*/
	void tfx_SetMinimumAllocationSize(tfx_allocator *allocator, tfx_size size);

	/*
	To help with debugging memory you can label a memory block with a type which gets stored in the block header.	
	*/
	void tfx_SetBlockMemoryType(void *allocation, tfx_size memory_type);

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

	//Read only functions
	static inline bool tfx__has_free_block(const tfx_allocator *allocator, tfx_index fli, tfx_index sli) {
		return allocator->first_level_bitmap & (TFX_ONE << fli) && allocator->second_level_bitmaps[fli] & (1U << sli);
	}

	static inline bool tfx__is_used_block(const tfx_header *block) {
		return !(block->size & tfx__BLOCK_IS_FREE);
	}

	static inline bool tfx__is_free_block(const tfx_header *block) {
		//Crashing here? The most likely reason is a pointer into the allocation for this block that became invalid but was still written to at some point.
		//Most likeyly cause is a tfx_vector_t or similar being resized and allocated elsewhere but you didn't account for this happening and update the pointer. Just index
		//into the array instead to fix these issues.
		//Another reason is simply that you're trying to free something that isn't actually a block of memory in the allocator, maybe you're just trying to free a struct or
		//other object?
		return block->size & tfx__BLOCK_IS_FREE;
	}

	static inline bool tfx__prev_is_free_block(const tfx_header *block) {
		return block->size & tfx__PREV_BLOCK_IS_FREE;
	}

	static inline void *tfx__align_ptr(const void *ptr, tfx_size align) {
		ptrdiff_t aligned = (((ptrdiff_t)ptr) + (align - 1)) & ~(align - 1);
		TFX_ASSERT(0 == (align & (align - 1)) && "must align to a power of two");
		return (void *)aligned;
	}

	static inline bool tfx__is_aligned(tfx_size size, tfx_size alignment) {
		return (size % alignment) == 0;
	}

	static inline bool tfx__ptr_is_aligned(void *ptr, tfx_size alignment) {
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
		return (tfx_header *)((char *)allocation - tfx__BLOCK_POINTER_OFFSET);
	}

	static inline tfx_header *tfx__null_block(tfx_allocator *allocator) {
		return &allocator->null_block;
	}

	static inline void *tfx__block_user_ptr(const tfx_header *block) {
		return (char *)block + tfx__BLOCK_POINTER_OFFSET;
	}

	static inline tfx_header *tfx__first_block_in_pool(const tfx_pool *pool) {
		return (tfx_header *)((char *)pool - tfx__POINTER_SIZE);
	}

	static inline tfx_header *tfx__next_physical_block(const tfx_header *block) {
		return (tfx_header *)((char *)tfx__block_user_ptr(block) + tfx__block_size(block));
	}

	static inline bool tfx__next_block_is_free(const tfx_header *block) {
		return tfx__is_free_block(tfx__next_physical_block(block));
	}

	static inline tfx_header *tfx__allocator_first_block(tfx_allocator *allocator) {
		return (tfx_header *)((char *)allocator + tfx_AllocatorSize() - tfx__POINTER_SIZE);
	}

	static inline bool tfx__is_last_block_in_pool(const tfx_header *block) {
		return tfx__block_size(block) == 0;
	}

	static inline tfx_index tfx__find_next_size_up(tfx_fl_bitmap map, tfx_uint start) {
		//Mask out all bits up to the start point of the scan
		map &= (~0ULL << (start + 1));
		return tfx__scan_forward(map);
	}

	//Debug tool to make sure that if a first level bitmap has a bit set, then the corresponding second level index should contain a value
	//It also checks that the blocks in the list are valid.
	//The most common cause of asserts here is where memory has been written to the wrong address. Check for buffers where they where resized
	//but the buffer pointer that was being written too was not updated after the resize for example.
	static inline void tfx__verify_lists(tfx_allocator *allocator) {
		for (int fli = 0; fli != tfx__FIRST_LEVEL_INDEX_COUNT; ++fli) {
			if (allocator->first_level_bitmap & (1ULL << fli)) {
				//bit in first level is set but according to the second level bitmap array there are no blocks so the first level
				//bitmap bit should have been 0
				TFX_ASSERT(allocator->second_level_bitmaps[fli] > 0);
				int sli = -1;
				sli = tfx__find_next_size_up(allocator->second_level_bitmaps[fli], sli);
				while (sli != -1) {
					tfx_header *block = allocator->segregated_lists[fli][sli];
					bool is_free = tfx__is_free_block(block);
					TFX_ASSERT(is_free);    //The block should be marked as free
					TFX_ASSERT(block->prev_free_block == &allocator->null_block);    //The first block in in the list should have a prev_free_block that points to the null block in the allocator
					sli = tfx__find_next_size_up(allocator->second_level_bitmaps[fli], sli);
				}
			}
		}
	}

#if defined(TFX_THREAD_SAFE)

#define tfx__lock_thread_access(alloc) tfx__sync_lock(&alloc->mutex);
#define tfx__unlock_thread_access(alloc) tfx__sync_unlock(&alloc->mutex);
#define tfx__init_allocator_mutex(alloc) tfx__sync_init(&alloc->mutex);

#else

#define tfx__lock_thread_access
#define tfx__unlock_thread_access 
#define tfx__init_allocator_mutex

#endif

	//Write functions
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
		} else {
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
		tfx_header *trimmed = (tfx_header *)((char *)tfx__block_user_ptr(block) + size);
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
		tfx_header *trimmed = (tfx_header *)((char *)tfx__block_user_ptr(block) + size_minus_overhead);
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
		TFX_ASSERT(next_block->prev_physical_block == block);    //could be potentional memory corruption. Check that you're not writing outside the boundary of the block size
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
			return block;
		}
		if (sli == tfx__SECOND_LEVEL_INDEX_COUNT - 1) {
			sli = -1;
		} else {
			sli = tfx__find_next_size_up(allocator->second_level_bitmaps[fli], sli);
		}
		if (sli == -1) {
			fli = tfx__find_next_size_up(allocator->first_level_bitmap, fli);
			if (fli > -1) {
				sli = tfx__scan_forward(allocator->second_level_bitmaps[fli]);
				tfx_header *block = tfx__pop_block(allocator, fli, sli);
				tfx_header *split_block = tfx__maybe_split_block(allocator, block, size, 0);
				return split_block;
			}
		} else {
			tfx_header *block = tfx__pop_block(allocator, fli, sli);
			tfx_header *split_block = tfx__maybe_split_block(allocator, block, size, 0);
			return split_block;
		}

		return 0;
	}
	//--End of internal functions

	//--End of header declarations

//Implementation
#if defined(TFX_ALLOCATOR_IMPLEMENTATION)

//Definitions
	void *tfx_BlockUserExtensionPtr(const tfx_header *block) {
		return (char *)block + sizeof(tfx_header);
	}

	void *tfx_AllocationFromExtensionPtr(const void *block) {
		return (void *)((char *)block - tfx__MINIMUM_BLOCK_SIZE);
	}

	tfx_allocator *tfx_InitialiseAllocator(void *memory) {
		if (!memory) {
			TFX_PRINT_ERROR(TFX_ERROR_COLOR"%s: The memory pointer passed in to the initialiser was NULL, did it allocate properly?\n", TFX_ERROR_NAME);
			return 0;
		}

		tfx_allocator *allocator = (tfx_allocator *)memory;
		memset(allocator, 0, sizeof(tfx_allocator));
		allocator->null_block.next_free_block = &allocator->null_block;
		allocator->null_block.prev_free_block = &allocator->null_block;
		allocator->minimum_allocation_size = tfx__MINIMUM_BLOCK_SIZE;
		tfx__init_allocator_mutex(allocator);

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
		TFX_ASSERT(allocator->minimum_allocation_size == tfx__MINIMUM_BLOCK_SIZE);        //You cannot change this once set
		TFX_ASSERT(tfx__is_pow2(size));                                                    //Size must be a power of 2
		allocator->minimum_allocation_size = tfx__Max(tfx__MINIMUM_BLOCK_SIZE, size);
	}

	void tfx_SetBlockMemoryType(void *allocation, tfx_size memory_type) {
		tfx_header *block = tfx__block_from_allocation(allocation);
#if defined(TFX_MEMORY_TRACKING)
		block->memory_block_type = memory_type;
#endif
	}

	tfx_pool *tfx_GetPool(tfx_allocator *allocator) {
		return (tfx_pool *)((char *)allocator + tfx_AllocatorSize());
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

	bool tfx_RemovePool(tfx_allocator *allocator, tfx_pool *pool) {
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
			tfx__unlock_thread_access(allocator);
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
		} else {
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
				const void *next_aligned = (void *)((ptrdiff_t)aligned_ptr + offset);

				aligned_ptr = tfx__align_ptr(next_aligned, alignment);
				gap = (tfx_size)((ptrdiff_t)aligned_ptr - (ptrdiff_t)user_ptr);
			}

			if (gap)
			{
				TFX_ASSERT(gap >= gap_minimum && "gap size too small");
				block = tfx__split_aligned_block(allocator, block, gap);
				tfx__block_set_used(block);
			}
			TFX_ASSERT(tfx__ptr_is_aligned(tfx__block_user_ptr(block), alignment));    //pointer not aligned to requested alignment
		} else {
			tfx__unlock_thread_access(allocator);
			return 0;
		}

		tfx__unlock_thread_access(allocator);
		return tfx__block_user_ptr(block);
	}

	int tfx_Free(tfx_allocator *allocator, void *allocation) {
		if (!allocation) return 0;
		tfx__lock_thread_access(allocator);
		tfx_header *block = tfx__block_from_allocation(allocation);
		if (tfx__prev_is_free_block(block)) {
			TFX_ASSERT(block->prev_physical_block);        //Must be a valid previous physical block
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
void *tfxAllocate(size_t size);
void *tfxReallocate(void *memory, size_t size);
void *tfxAllocateAligned(size_t size, size_t alignment);
//Do a safe copy where checks are made to ensure that the boundaries of the memory block being copied to are respected
//This assumes that dst is the start address of the block. If you're copying to a range that is offset from the beginning
//of the block then you can use tfx_SafeCopyBlock instead.
bool tfx_SafeCopy(void *dst, void *src, tfx_size size);
bool tfx_SafeCopyBlock(void *dst_block_start, void *dst, void *src, tfx_size size);
bool tfx_SafeMemset(void *allocation, void *dst, int value, tfx_size size);
tfx_allocator *tfxGetAllocator();

//---------------------------------------
//End of allocator code
//---------------------------------------

//----------------------------------------------------------
//Header_Includes_and_Typedefs
//----------------------------------------------------------
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64)
#define tfxINTEL
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

//--------------------------------------------------------------
//Macros

#define tfxTWO63 0x8000000000000000u 
#define tfxTWO64f 1.8446744073709552E19		//tfxTWO63 * 2
#define tfxPI 3.14159265359f
#define tfxHALFPI 1.570796f
#define tfxPI2 6.283185307f 
#define tfxINVTWOPI 0.1591549f
#define tfxTHREEHALFPI 4.7123889f
#define tfxQUARTERPI 0.7853982f
#define tfx720Radians 12.56638f
#define tfx360Radians 6.28319f
#define tfx180Radians 3.14159f
#define tfx90Radians 1.5708f
#define tfx270Radians 4.71239f
#define tfxNL u8"\n"
#define tfxPROPERTY_INDEX_MASK 0x00007FFF
#define tfxSPRITE_ALIGNMENT_MASK 0xFF000000
#define tfxSPRITE_IMAGE_FRAME_MASK 0x00FF0000
#define tfxEXTRACT_SPRITE_ALIGNMENT(property_index) ((property_index & tfxSPRITE_ALIGNMENT_MASK) >> 24)
#define tfxEXTRACT_SPRITE_IMAGE_FRAME(property_index) ((property_index & tfxSPRITE_IMAGE_FRAME_MASK) >> 16)
#define tfxEXTRACT_SPRITE_PROPERTY_INDEX(property_index) (property_index & tfxPROPERTY_INDEX_MASK)
#define tfxPACK_SCALE_AND_HANDLE(x, y, lib, property_index) (tfxU16)(x * 127.9960938f) | ((tfxU16)(y * 127.9960938f) << 16) | ((tfxU64)lib->emitter_properties[property_index].image_handle_packed << 32)
#define tfxPACK_SIZE_AND_HANDLE(x, y, lib, property_index) (tfxU16)(x * 7.999755859f) | ((tfxU16)(y * 7.999755859f) << 16) | ((tfxU64)lib->emitter_properties[property_index].image_handle_packed << 32)
#define tfxOSCILLATOR_SIN(t, frequency, amplitude) 0.5f + sinf((t) * (frequency) * 6.28318f) * (amplitude)
#define tfxOSCILLATOR_WIDE_SIN(t, frequency, amplitude) tfxWideAdd(tfxWIDEHALF.m, tfxWideMul(tfxWideSin(tfxWideMul(tfxWideMul(t, frequency), tfxWIDEPI2.m)), amplitude))
#define tfxCIRCLENODES 16
#define tfxMIN_SEGMENT_COUNT 32 
#define tfxMAX_SEGMENT_COUNT 1024 
#define tfxPrint(message, ...) printf(message tfxNL, ##__VA_ARGS__)

#define TFX_VERSION "Alpha"
#define TFX_VERSION_NUMBER 6.18.2024

#define tfxMAX_FRAME 20000.f
#define tfxCOLOR_RAMP_WIDTH 256
#define tfxNullParent 0xFFFFFFFF
#define tfxINVALID 0xFFFFFFFF
#define tfxINVALID_SPRITE 0x0FFFFFFF
#define tfxUNINIT_MAGIC 0xDEADBEEF
#define tfxINIT_MAGIC 0xCAFEBABE
#define tfxEmitterPropertiesCount 26

#define tfxDel << "=" <<
#define tfxCom << "," <<
#define tfxEndLine << std::endl

#define tfxDelt "=" 
#define tfxComt ","

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
#define tfxForEachLayer for (int layer = 0; layer != tfxLAYERS; ++layer)

//----------------------------------------------------------
//Forward declarations

typedef struct tfx_effect_descriptor_s tfx_effect_descriptor_t;
typedef struct tfx_effect_manager_s tfx_effect_manager_t;
typedef struct tfx_effect_template_s tfx_effect_template_t;
typedef struct tfx_compute_particle_s tfx_compute_particle_t;
typedef struct tfx_library_s tfx_library_t;
typedef struct tfx_animation_manager_s tfx_animation_manager_t;

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

union tfxUInt8bit
{
	struct
	{
		int x : 8;
		int y : 8;
		int z : 8;
		int w : 8;
	} data;
	tfxU32 pack;
};

//Section: OS_Specific_Functions
#ifdef _WIN32
FILE *tfx__open_file(const char *file_name, const char *mode);

tfxINTERNAL inline tfxU64 tfx_AtomicAdd64(tfxU64 volatile *value, tfxU64 amount_to_add) {
	tfxU64 result = _InterlockedExchangeAdd64((__int64 volatile *)value, amount_to_add);
	return result;
}

tfxINTERNAL inline tfxU32 tfx_AtomicAdd32(tfxU32 volatile *value, tfxU32 amount_to_add) {
	tfxU32 result = _InterlockedExchangeAdd((LONG *)value, amount_to_add);
	return result;
}
#else
FILE *tfx__open_file(const char *file_name, const char *mode);

inline tfxU64 tfx_AtomicAdd64(tfxU64 volatile *value, tfxU64 amount_to_add) {
	return __sync_fetch_and_add(value, amount_to_add);
}

inline tfxU32 tfx_AtomicAdd32(tfxU32 volatile *value, tfxU32 amount_to_add) {
	return __sync_fetch_and_add(value, amount_to_add);
}
#endif

#ifdef _WIN32
tfxINTERNAL tfxU32 tfx_Millisecs(void) {
	LARGE_INTEGER frequency, counter;
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&counter);
	tfxU64 ms = (tfxU64)(counter.QuadPart * 1000LL / frequency.QuadPart);
	return (tfxU32)ms;
}

tfxINTERNAL tfxU64 tfx_Microsecs(void) {
	LARGE_INTEGER frequency, counter;
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&counter);
	tfxU64 us = (tfxU64)(counter.QuadPart * 1000000LL / frequency.QuadPart);
	return (tfxU64)us;
}
#elif defined(__APPLE__)
#include <mach/mach_time.h>
tfxINTERNAL inline tfxU32 tfx_Millisecs(void) {
	static mach_timebase_info_data_t timebase_info;
	if (timebase_info.denom == 0) {
		mach_timebase_info(&timebase_info);
	}

	uint64_t time_ns = mach_absolute_time() * timebase_info.numer / timebase_info.denom;
	tfxU32 ms = (tfxU32)(time_ns / 1000000);
	return (tfxU32)ms;
}

tfxINTERNAL inline tfxU64 tfx_Microsecs(void) {
	static mach_timebase_info_data_t timebase_info;
	if (timebase_info.denom == 0) {
		mach_timebase_info(&timebase_info);
	}

	uint64_t time_ns = mach_absolute_time() * timebase_info.numer / timebase_info.denom;
	tfxU64 us = (tfxU64)(time_ns / 1000);
	return us;
}
#else
tfxINTERNAL inline tfxU32 tfx_Millisecs(void) {
	struct timespec now;
	clock_gettime(CLOCK_REALTIME, &now);
	long m = now.tv_sec * 1000 + now.tv_nsec / 1000000;
	return (tfxU32)m;
}

tfxINTERNAL inline tfxU64 tfx_Microsecs(void) {
	struct timespec now;
	clock_gettime(CLOCK_REALTIME, &now);
	tfxU64 us = now.tv_sec * 1000000ULL + now.tv_nsec / 1000;
	return (tfxU64)us;
}
#endif

// --Pocket_Hasher, converted to c from Stephen Brumme's XXHash code (https://github.com/stbrumme/xxhash) by Peter Rigby
/*
	MIT License Copyright (c) 2018 Stephan Brumme
	Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
	The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  */
#define tfx__PRIME1 11400714785074694791ULL
#define tfx__PRIME2 14029467366897019727ULL
#define tfx__PRIME3 1609587929392839161ULL
#define tfx__PRIME4 9650029242287828579ULL
#define tfx__PRIME5 2870177450012600261ULL

enum { tfx__HASH_MAX_BUFFER_SIZE = 31 + 1 };

typedef struct tfx_hasher_s {
	tfxU64      state[4];
	unsigned char buffer[tfx__HASH_MAX_BUFFER_SIZE];
	tfxU64      buffer_size;
	tfxU64      total_length;
} tfx_hasher_t;

tfxINTERNAL inline tfxU64 tfx__hash_rotate_left(tfxU64 x, unsigned char bits) {
	return (x << bits) | (x >> (64 - bits));
}
tfxINTERNAL inline tfxU64 tfx__hash_process_single(tfxU64 previous, tfxU64 input) {
	return tfx__hash_rotate_left(previous + input * tfx__PRIME2, 31) * tfx__PRIME1;
}
tfxINTERNAL inline void tfx__hasher_process(const void *data, tfxU64 *state0, tfxU64 *state1, tfxU64 *state2, tfxU64 *state3) {
	tfxU64 *block = (tfxU64 *)data;
	tfxU64 blocks[4];
	memcpy(blocks, data, sizeof(tfxU64) * 4);
	*state0 = tfx__hash_process_single(*state0, blocks[0]);
	*state1 = tfx__hash_process_single(*state1, blocks[1]);
	*state2 = tfx__hash_process_single(*state2, blocks[2]);
	*state3 = tfx__hash_process_single(*state3, blocks[3]);
}
tfxINTERNAL inline int tfx__hasher_add(tfx_hasher_t *hasher, const void *input, tfxU64 length)
{
	if (!input || length == 0) return 0;

	hasher->total_length += length;
	unsigned char *data = (unsigned char *)input;

	if (hasher->buffer_size + length < tfx__HASH_MAX_BUFFER_SIZE)
	{
		while (length-- > 0)
			hasher->buffer[hasher->buffer_size++] = *data++;
		return 1;
	}

	const unsigned char *stop = data + length;
	const unsigned char *stopBlock = stop - tfx__HASH_MAX_BUFFER_SIZE;

	if (hasher->buffer_size > 0)
	{
		while (hasher->buffer_size < tfx__HASH_MAX_BUFFER_SIZE)
			hasher->buffer[hasher->buffer_size++] = *data++;

		tfx__hasher_process(hasher->buffer, &hasher->state[0], &hasher->state[1], &hasher->state[2], &hasher->state[3]);
	}

	tfxU64 s0 = hasher->state[0], s1 = hasher->state[1], s2 = hasher->state[2], s3 = hasher->state[3];
	int test = tfx__ptr_is_aligned(&s0, 8);
	while (data <= stopBlock)
	{
		tfx__hasher_process(data, &s0, &s1, &s2, &s3);
		data += 32;
	}
	hasher->state[0] = s0; hasher->state[1] = s1; hasher->state[2] = s2; hasher->state[3] = s3;

	hasher->buffer_size = stop - data;
	for (tfxU64 i = 0; i < hasher->buffer_size; i++)
		hasher->buffer[i] = data[i];

	return 1;
}

tfxINTERNAL tfxU64 tfx__get_hash(tfx_hasher_t *hasher)
{
	tfxU64 result;
	if (hasher->total_length >= tfx__HASH_MAX_BUFFER_SIZE)
	{
		result = tfx__hash_rotate_left(hasher->state[0], 1) +
			tfx__hash_rotate_left(hasher->state[1], 7) +
			tfx__hash_rotate_left(hasher->state[2], 12) +
			tfx__hash_rotate_left(hasher->state[3], 18);
		result = (result ^ tfx__hash_process_single(0, hasher->state[0])) * tfx__PRIME1 + tfx__PRIME4;
		result = (result ^ tfx__hash_process_single(0, hasher->state[1])) * tfx__PRIME1 + tfx__PRIME4;
		result = (result ^ tfx__hash_process_single(0, hasher->state[2])) * tfx__PRIME1 + tfx__PRIME4;
		result = (result ^ tfx__hash_process_single(0, hasher->state[3])) * tfx__PRIME1 + tfx__PRIME4;
	} else
	{
		result = hasher->state[2] + tfx__PRIME5;
	}

	result += hasher->total_length;
	const unsigned char *data = hasher->buffer;
	const unsigned char *stop = data + hasher->buffer_size;
	for (; data + 8 <= stop; data += 8)
		result = tfx__hash_rotate_left(result ^ tfx__hash_process_single(0, *(tfxU64 *)data), 27) * tfx__PRIME1 + tfx__PRIME4;
	if (data + 4 <= stop)
	{
		result = tfx__hash_rotate_left(result ^ (*(tfxU32 *)data) * tfx__PRIME1, 23) * tfx__PRIME2 + tfx__PRIME3;
		data += 4;
	}
	while (data != stop)
		result = tfx__hash_rotate_left(result ^ (*data++) * tfx__PRIME5, 11) * tfx__PRIME1;
	result ^= result >> 33;
	result *= tfx__PRIME2;
	result ^= result >> 29;
	result *= tfx__PRIME3;
	result ^= result >> 32;
	return result;
}

tfxAPI_EDITOR inline void tfx__hash_initialise(tfx_hasher_t *hasher, tfxU64 seed) {
	hasher->state[0] = seed + tfx__PRIME1 + tfx__PRIME2; hasher->state[1] = seed + tfx__PRIME2; hasher->state[2] = seed; hasher->state[3] = seed - tfx__PRIME1; hasher->buffer_size = 0; hasher->total_length = 0;
}

//The only command you need for the hasher. Just used internally by the hash map.
tfxAPI tfxKey tfx_Hash(tfx_hasher_t *hasher, const void *input, tfxU64 length, tfxU64 seed);
//-- End of Pocket Hasher

//----------------------------------------------------------
//Section: SIMD_defines
//----------------------------------------------------------

#define tfx128SetConst(value) {value, value, value, value}

//Define tfxUSEAVX if you want to compile and use AVX simd operations for updating particles, otherwise SSE will be
//used by default
//Note that avx is currently only slightly faster than SSE, probably because memory bandwidth/caching becomes more of an issue at that point. But also I could be doing it wrong!
#ifdef tfxUSEAVX
#define tfxDataWidth 8    
typedef __m256 tfxWideFloat;
typedef __m256i tfxWideInt;
typedef __m256i tfxWideIntLoader;
#define tfxWideLoad _mm256_load_ps
#define tfxWideLoadi _mm256_load_si256
#define tfx128Set _mm_set_ps
#define tfx128Seti _mm_set_epi32
#define tfxWideLoadHalfs(mem_address) _mm_load_si128((tfx128i*)mem_address)
#define tfxWideStoreHalfs _mm_store_si128
#define tfx128SetSingle _mm_set_ps1
#define tfx128SetSinglei _mm_set1_epi32
#define tfx128Load64bytes _mm_loadl_epi64 
#define tfx128UnpackLo8 _mm_unpacklo_epi8
#define tfx128UnpackLo16 _mm_unpacklo_epi16
#define tfx128UnpackHi16 _mm_unpackhi_epi16
#define tfx128SetZeroi _mm_setzero_si128()
#define tfx128Packus32 _mm_packus_epi32
#define tfx128Packus16 _mm_packus_epi16
#define tfx128Convert32 _mm_cvtsi128_si32
#define tfxWideSet _mm256_set_ps
#define tfxWideSetSingle _mm256_set1_ps
#define tfxWideSeti _mm256_set_epi32
#define tfxWideSetSinglei _mm256_set1_epi32
#define tfxWideSet128i _mm256_set_m128i
#define tfxWideExtract128i _mm256_extractf128_si256
#define tfxWideBlend _mm256_blend_ps
#define tfxWideBlendv _mm256_blendv_ps
#define tfxWideAdd _mm256_add_ps
#define tfxWideSub _mm256_sub_ps
#define tfxWideMul _mm256_mul_ps
#define tfxWideDiv _mm256_div_ps
#ifdef tfxUSEFMA
#define tfxWideMulAdd(a, b, c) _mm256_fmadd_ps(a, b, c)
#else
#define tfxWideMulAdd(a, b, c) tfxWideAdd(tfxWideMul(a, b), c)
#endif
#define tfxWideAddi _mm256_add_epi32
#define tfxWideSubi _mm256_sub_epi32
#define tfxWideMuli _mm256_mullo_epi32
#define tfxWideSqrt _mm256_sqrt_ps
#define tfxWideRSqrt _mm256_rsqrt_ps
#define tfxWideMoveMask _mm256_movemask_epi8
#define tfxWideMovemasKps _mm256_movemask_ps
#define tfxWideShiftRight _mm256_srli_epi32
#define tfxWideShiftLeft _mm256_slli_epi32
#define tfxWideGreaterEqual(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_GE_OS)
#define tfxWideGreater(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_GT_OS)
#define tfxWideGreateri _mm256_cmpgt_epi32
#define tfxWideLess(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_LT_OS)
#define tfxWideLessi(v1, v2) _mm256_cmpgt_epi32(v2, v1)
#define tfxWideLessEqual(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_LE_OS)
#define tfxWideEquals(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_EQ_OS)
#define tfxWideEqualsi _mm256_cmpeq_epi32 
#define tfxWideStore _mm256_store_ps
#define tfxWideStorei _mm256_store_si256
#define tfxWideCasti _mm256_castps_si256
#define tfxWideCast _mm256_castsi256_ps 
#define tfxWideConverti _mm256_cvttps_epi32 
#define tfxWideConvert _mm256_cvtepi32_ps 
#ifdef tfxHALFFLOATS
#define tfxWideConvertHalfsToFloats _mm256_cvtph_ps 
#define tfxWideConvertFloatsToHalfs(floats) _mm256_cvtps_ph(floats, _MM_FROUND_TO_NEAREST_INT)
#endif
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
#define tfxWideSetZero _mm256_setzero_ps()
#define tfxWideSetZeroi _mm256_setzero_si256()
#define tfxWideEqualsi _mm256_cmpeq_epi32 
#define tfxWideFloor _mm256_floor_ps
#define tfxWideLookupSet(lookup, index) tfxWideSet(lookup[index.a[7]], lookup[index.a[6]], lookup[index.a[5]], lookup[index.a[4]], lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )
#define tfxWideLookupSetOffset(lookup, index, offset) tfxWideSet(lookup[index.a[7] + offset], lookup[index.a[6] + offset], lookup[index.a[5] + offset], lookup[index.a[4] + offset], lookup[index.a[3] + offset], lookup[index.a[2] + offset], lookup[index.a[1] + offset], lookup[index.a[0] + offset] )

#define tfxWideSetConst(value) {value, value, value, value, value, value, value, value}

typedef union {
	int a[8];
	__m256i m;
} tfxWideArrayi;

typedef union {
	float a[8];
	__m256 m;
} tfxWideArray;

const tfxWideArrayi tfxBASEINDEX = { 0, 1, 2, 3, 4, 5, 6, 7 };

#else

#define tfxDataWidth 4

#ifdef tfxINTEL
//Intel Intrinsics
typedef __m128 tfxWideFloat;
typedef __m128i tfxWideInt;
typedef __m128i tfxWideIntLoader;
#define tfxWideLoad _mm_load_ps
#define tfxWideLoadi _mm_load_si128
#define tfxWideSet _mm_set_ps
#define tfxWideSetSingle _mm_set_ps1
#define tfx128Load _mm_load_ps
#define tfx128Set _mm_set_ps
#define tfx128Seti _mm_set_epi32
#define tfx128SetSingle _mm_set_ps1
#define tfx128SetSinglei _mm_set1_epi32
#define tfxWideSeti _mm_set_epi32
#define tfxWideSetSinglei _mm_set1_epi32
#define tfxWideBlend _mm_blend_ps
#define tfxWideBlendv(a, b, mask) _mm_blendv_ps(a, b, mask)
#define tfxWideAdd _mm_add_ps
#define tfxWideSub _mm_sub_ps
#define tfxWideMul _mm_mul_ps
#define tfxWideDiv _mm_div_ps
#ifdef tfxUSEFMA
#define tfxWideMulAdd(a, b, c) _mm_fmadd_ps(a, b, c)
#define tfxWideMulSub(a, b, c) _mm_fmsub_ps(a, b, c)
#else
#define tfxWideMulAdd(a, b, c) tfxWideAdd(tfxWideMul(a, b), c)
#define tfxWideMulSub(a, b, c) tfxWideSub(tfxWideMul(a, b), c)
#endif
#define tfxWideLoadHalfs(mem_address) _mm_load_si128((tfx128i*)mem_address)
#define tfxWideStoreHalfs _mm_store_si128
#ifdef tfxHALFFLOATS
#define tfxWideConvertHalfsToFloats _mm_cvtph_ps 
#define tfxWideConvertFloatsToHalfs(floats) _mm_cvtps_ph(floats, _MM_FROUND_TO_NEAREST_INT)
#endif
#define tfxWideAddi _mm_add_epi32
#define tfxWideSubi _mm_sub_epi32
#define tfxWideMuli _mm_mullo_epi32
#define tfxWideSqrt _mm_sqrt_ps
#define tfxWideRSqrt _mm_rsqrt_ps
#define tfxWideMoveMask _mm_movemask_epi8
#define tfxWideMovemasKps _mm_movemask_ps
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
#define tfxWideXOr _mm_xor_ps
#define tfxWideOri _mm_or_si128
#define tfxWideXOri _mm_xor_si128
#define tfxWideAnd _mm_and_ps
#define tfxWideAndNot _mm_andnot_ps
#define tfxWideAndi _mm_and_si128
#define tfxWideAndNoti _mm_andnot_si128
#define tfxWideSetZeroi _mm_setzero_si128()
#define tfxWideSetZero _mm_setzero_ps()
#define tfxWideEqualsi _mm_cmpeq_epi32
#define tfxWideEquals _mm_cmpeq_ps
#define tfxWideShufflei _mm_shuffle_epi32
#ifdef SSE41
#define tfxWideFloor _mm_floor_ps
#else
#define tfxWideFloor tfxFloor128
#endif

#define tfxWideSetConst(value) {value, value, value, value}

typedef union {
	int a[4];
	__m128i m;
} tfxWideArrayi;

typedef union {
	float a[4];
	__m128 m;
} tfxWideArray;

#elif defined(tfxARM)
//Arm Intrinsics
typedef float32x4_t tfxWideFloat;
typedef int32x4_t tfxWideInt;
typedef int32_t tfxWideIntLoader;
#define tfxWideLoad vld1q_f32
#define tfxWideLoadi vld1q_s32
inline __attribute__((always_inline)) float32x4_t tfx__128_SET(float e3, float e2, float e1, float e0) {
	float32x4_t r;
	float data[4] __attribute__((aligned(16))) = { e0, e1, e2, e3 };
	r = vld1q_f32(data);
	return r;
}
inline __attribute__((always_inline)) int32x4_t tfx__128i_SET(int e3, int e2, int e1, int e0) {
	int32x4_t r;
	int data[4] __attribute__((aligned(16))) = { e0, e1, e2, e3 };
	r = vld1q_s32(data);
	return r;
}
#define tfx128Load vld1q_f32
#define tfx128Set tfx__128_SET
#define tfx128Seti tfx__128i_SET
#define tfx128SetSingle vdupq_n_f32
#define tfx128SetSinglei vdupq_n_s32
#define tfxWideSet vld1q_f32
#define tfxWideSetSingle vdupq_n_f32
#define tfxWideSeti vld1q_s32
#define tfxWideSetSinglei vdupq_n_s32
#define tfxWideBlendv(a, b, mask) vbslq_f32(vreinterpretq_u32_f32(mask), b, a)
#define tfxWideAdd vaddq_f32
#define tfxWideSub vsubq_f32
#define tfxWideMul vmulq_f32
#define tfxWideDiv vdivq_f32
#define tfxWideMulAdd(a, b, c) vfmaq_f32(c, a, b)
#define tfxWideMulSub(a, b, c) vsubq_f32(vmulq_f32(a, b), c)
#define tfxWideFloor vrndmq_f32
#define tfxWideAddi vaddq_s32
#define tfxWideSubi vsubq_s32
#define tfxWideMuli vmulq_s32
#define tfxWideRSqrt vrsqrteq_f32 // for reciprocal square root approximation
#define tfxWideShiftRight vshrq_n_u32
#define tfxWideShiftLeft vshlq_n_u32
#define tfxWideGreaterEqual(a, b) vreinterpretq_f32_u32(vcgeq_f32(a, b))
#define tfxWideGreater(a, b) vreinterpretq_s32_u32(vcgeq_f32(a,b))
#define tfxWideGreateri vcgtq_s32
#define tfxWideLessEqual(a, b) vreinterpretq_f32_u32(vcleq_f32(a, b))
#define tfxWideLess(a, b) vreinterpretq_f32_u32(vcltq_f32(a, b))
#define tfxWideLessi vcltq_s32
#define tfxWideStore vst1q_f32
#define tfxWideStorei vst1q_s32
#define tfxWideCasti vreinterpretq_s32_f32
#define tfxWideCast vreinterpretq_f32_s32
#define tfxWideConverti vcvtq_s32_f32
#define tfxWideConvert vcvtq_f32_s32
#define tfxWideMin vminq_f32
#define tfxWideMax vmaxq_f32
#define tfxWideMini vminq_s32
#define tfxWideMaxi vmaxq_s32
#define tfxWideOr(a, b) vreinterpretq_f32_s32(vorrq_s32(vreinterpretq_s32_f32(a), vreinterpretq_s32_f32(b)))
#define tfxWideXOr(a, b) vreinterpretq_f32_s32(veorq_s32(vreinterpretq_s32_f32(a), vreinterpretq_s32_f32(b)))
#define tfxWideOri vorrq_s32
#define tfxWideXOri veorq_s32
#define tfxWideAnd(a, b) vreinterpretq_f32_s32(vandq_s32(vreinterpretq_s32_f32(a), vreinterpretq_s32_f32(b)))
#define tfxWideAndi vandq_s32
#define tfxWideAndNoti vbicq_s32
//#define tfxWideAndNot(a, b) vreinterpretq_f32_s32(vandq_s32(vmvnq_s32(vreinterpretq_s32_f32(a)), vreinterpretq_s32_f32(b)))
#define tfxWideAndNot(a, b) vreinterpretq_f32_s32(vbicq_s32(vreinterpretq_s32_f32(b), vreinterpretq_s32_f32(a)))
#define tfxWideSetZeroi vdupq_n_s32(0)
#define tfxWideSetZero vdupq_n_f32(0.0f)
#define tfxWideEqualsi vceqq_s32
#define tfxWideEquals(a, b) vreinterpretq_f32_u32(vceqq_f32(a, b))
static inline int tfx__neon_movemask_ps(float32x4_t v) {
	uint32x4_t s = vshrq_n_u32(vreinterpretq_u32_f32(v), 31);
	return (int)(vgetq_lane_u32(s, 0) | (vgetq_lane_u32(s, 1) << 1) | (vgetq_lane_u32(s, 2) << 2) | (vgetq_lane_u32(s, 3) << 3));
}
#define tfxWideMovemasKps tfx__neon_movemask_ps

#define tfxSIMD_AND(a,b) vreinterpretq_f32_s32(vandq_s32(vreinterpretq_s32_f32(a),vreinterpretq_s32_f32(b)))
#define tfxSIMD_AND_NOT(a,b) vreinterpretq_f32_s32(vandq_s32(vmvnq_s32(vreinterpretq_s32_f32(a)),vreinterpretq_s32_f32(b)))
#define tfxSIMD_XOR(a,b) vreinterpretq_f32_s32(veorq_s32(vreinterpretq_s32_f32(a),vreinterpretq_s32_f32(b)))

#define tfxWideSetConst(value) {value, value, value, value}

typedef union {
	int a[4];
	int32x4_t m;
} tfxWideArrayi;

typedef union {
	float a[4];
	float32x4_t m;
} tfxWideArray;

#endif

#define tfxWideLookupSet(lookup, index) tfx128Set( lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )
#define tfxWideLookupSetOffset(lookup, index, offset) tfx128Set( lookup[index.a[3] + offset], lookup[index.a[2] + offset], lookup[index.a[1] + offset], lookup[index.a[0] + offset] )

const tfxWideArrayi tfxBASEINDEX = { 0, 1, 2, 3 };

#endif

const tfxWideArray tfxWIDEF3_4        = tfxWideSetConst(1.f / 3.f);
const tfxWideArray tfxWIDEG3_4        = tfxWideSetConst(1.0f / 6.0f);
const tfxWideArray tfxWIDEG32_4       = tfxWideSetConst((1.0f / 6.0f) * 2.f);
const tfxWideArray tfxWIDEG33_4       = tfxWideSetConst((1.0f / 6.0f) * 3.f);
const tfxWideArrayi tfxWIDEONEi       = tfxWideSetConst(1);
const tfxWideArray tfxWIDEMINUSONE    = tfxWideSetConst(-1.f);
const tfxWideArrayi tfxWIDEMINUSONEi  = tfxWideSetConst(-1);
const tfxWideArray tfxWIDEONE         = tfxWideSetConst(1.f);
const tfxWideArray tfxWIDETWO         = tfxWideSetConst(2.f);
const tfxWideArray tfxWIDEMINUSTWO    = tfxWideSetConst(-2.f);
const tfxWideArray tfxWIDETHREE       = tfxWideSetConst(3.f);
const tfxWideArray tfxWIDEFOUR        = tfxWideSetConst(4.f);
const tfxWideArray tfxWIDEEIGHT       = tfxWideSetConst(8.f);
const tfxWideArray tfxWIDESIXTEEN     = tfxWideSetConst(16.f);
const tfxWideArray tfxWIDEHALF        = tfxWideSetConst(0.5f);
const tfxWideArray tfxWIDE255         = tfxWideSetConst(255.f);
const tfxWideArray tfxWIDE255r         = tfxWideSetConst(1.f / 255.f);
const tfxWideArray tfxWIDEZERO        = tfxWideSetConst(0.f);
const tfxWideArray tfxWIDETHIRTYTWO   = tfxWideSetConst(32.f);
const tfxWideArray tfxPWIDESIX        = tfxWideSetConst(0.6f);
const tfxWideArray tfxMAXUINTf        = tfxWideSetConst((float)UINT32_MAX);
const tfxWideArray tfxDEGREERANGEMR   = tfxWideSetConst(0.392699f);
const tfxWideArray tfxSIGNMASK        = tfxWideSetConst(-0.f);
//const tfxWideArray tfxABSMASK		  = tfxWideSetConst(-NAN);
const tfxWideArrayi tfxABSMASKi		  = tfxWideSetConst(0x7FFFFFFF);
const tfxWideArray tfxWIDEPI          = tfxWideSetConst(3.14159265359f);
const tfxWideArray tfxWIDEHALFPI      = tfxWideSetConst(1.570796f);
const tfxWideArray tfxWIDEPI2         = tfxWideSetConst(6.283185307f);
const tfxWideArray tfxWIDEINVTWOPI    = tfxWideSetConst(0.1591549f);
const tfxWideArray tfxWIDETHREEHALFPI = tfxWideSetConst(4.7123889f);
const tfxWideArray tfxWIDEQUARTERPI   = tfxWideSetConst(0.7853982f);
const tfxWideArray tfxWIDEEPS		  = tfxWideSetConst(0.0001f);
const tfxWideArray tfxWIDEEPS2		  = tfxWideSetConst(0.0002f);
const tfxWideArray tfxWIDENOISEOFFSET = tfxWideSetConst(100.f);

#ifdef tfxINTEL
typedef __m128 tfx128;
typedef __m128i tfx128i;

typedef union {
	int a[4];
	__m128i m;
} tfx128iArray;

typedef union {
	tfxU64 a[2];
	__m128i m;
} tfx128iArray64;

typedef union {
	float a[4];
	__m128 m;
} tfx128Array;

//simd floor function thanks to Stephanie Rancourt: http://dss.stephanierct.com/DevBlog/?p=8
tfxINTERNAL inline tfx128 tfxFloor128(const tfx128 x) {
	//__m128i v0 = _mm_setzero_si128();
	//__m128i v1 = _mm_cmpeq_epi32(v0, v0);
	//__m128i ji = _mm_srli_epi32(v1, 25);
	//__m128 j = *(__m128*)&_mm_slli_epi32(ji, 23); //create vector 1.0f
	//Worth noting that we only need to floor small numbers for the noise algorithm so can get away with this function.
	__m128 j = _mm_set1_ps(1.f); //create vector 1.0f I think this is faster now on more modern CPUs
	__m128i i = _mm_cvttps_epi32(x);
	__m128 fi = _mm_cvtepi32_ps(i);
	__m128 igx = _mm_cmpgt_ps(fi, x);
	j = _mm_and_ps(igx, j);
	return _mm_sub_ps(fi, j);
}

tfxINTERNAL inline uint64_t tfx__rdtsc() {
	return __rdtsc();
}

#elif defined(tfxARM)

typedef float32x4_t tfx128;
typedef int32x4_t tfx128i;

typedef union {
	int a[4];
	tfx128i m;
} tfx128iArray;

typedef union {
	tfxU64 a[2];
	tfx128i m;
} tfx128iArray64;

typedef union {
	float a[4];
	tfx128 m;
} tfx128Array;

#define tfxFloor128 vrndmq_f32

tfxINTERNAL inline uint64_t tfx__rdtsc() {
	return mach_absolute_time();
}

#endif

inline tfxWideArrayi tfxConvertWideArray(tfxWideInt in) {
	tfxWideArrayi out;
	out.m = in;
	return out;
}

/*
Copyright (c) 2013, Robin Lobel
All rights reserved. https://github.com/divideconcept/FastTrigo

I just extracted the necessary functions that I need from the above code base and modified to use tfx SIMD_DEFINES
Also fixed a bug in atan2 function where x <= y
*/

tfxINTERNAL inline tfxWideFloat tfxWideFastSqrt(tfxWideFloat squared)
{
	//    static int csr = 0;
	//#if defined(tfxINTEL)
		//if (!csr) csr = _mm_getcsr() | 0x8040; //DAZ,FTZ (divide by zero=0)
		//_mm_setcsr(csr);
	//#endif
	return tfxWideMul(tfxWideRSqrt(squared), squared);
}

/*
possible arm function for andnot
float32x4_t andnot_ps(float32x4_t a, float32x4_t b) {
	uint32x4_t not_b = vreinterpretq_u32_f32(vmvnq_f32(b)); // Bitwise NOT of b
	return vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(a), not_b)); // Bitwise AND of a and NOT(b)
}
*/

tfxINTERNAL inline tfxWideFloat tfxWideAtan(tfxWideFloat x)
{
	//                                      tfxQUARTERPI*x
	//                                      - x*(fabs(x) - 1)
	//                                      *(0.2447f+0.0663f*fabs(x));
	const tfxWideFloat a = tfxWideSetSingle(0.2447f);
	const tfxWideFloat b = tfxWideSetSingle(0.0663f);
	return tfxWideSub(tfxWideMul(tfxWIDEQUARTERPI.m, x),
		tfxWideMul(tfxWideMul(x, tfxWideSub(tfxWideAndNot(tfxSIGNMASK.m, x), tfxWIDEONE.m)),
			(tfxWideAdd(a, tfxWideMul(b, tfxWideAndNot(tfxSIGNMASK.m, x))))));
}

/*
float atan2(float y, float x)
{
	if (fabs(x) > fabs(y)) {
		float atan = atanf(y / x);
		if (x > 0.f)
			return atan;
		else
			return y > 0.f ? atan + pi : atan - pi;
	}
	else {
		float atan = atanf(x / y);
		return y > 0.f ? halfpi - atan : (-halfpi) - atan;
	}
}
*/

tfxINTERNAL inline tfxWideFloat tfxWideAtan2(tfxWideFloat y, tfxWideFloat x)
{
	tfxWideFloat absxgreaterthanabsy = tfxWideGreater(tfxWideAndNot(tfxSIGNMASK.m, x), tfxWideAndNot(tfxSIGNMASK.m, y));
	tfxWideFloat ratio = tfxWideDiv(tfxWideAdd(tfxWideAnd(absxgreaterthanabsy, y), tfxWideAndNot(absxgreaterthanabsy, x)),
		tfxWideAdd(tfxWideAnd(absxgreaterthanabsy, x), tfxWideAndNot(absxgreaterthanabsy, y)));
	tfxWideFloat atan = tfxWideAtan(ratio);

	tfxWideFloat xgreaterthan0 = tfxWideGreater(x, tfxWideSetZero);
	tfxWideFloat ygreaterthan0 = tfxWideGreater(y, tfxWideSetZero);

	atan = tfxWideXOr(atan, tfxWideAndNot(absxgreaterthanabsy, tfxSIGNMASK.m)); //negate atan if absx<=absy & x>0

	tfxWideFloat shift = tfxWIDEPI.m;
	shift = tfxWideSub(shift, tfxWideAndNot(absxgreaterthanabsy, tfxWIDEHALFPI.m)); //substract tfxHALFPI if absx<=absy
	shift = tfxWideXOr(shift, tfxWideAndNot(ygreaterthan0, tfxSIGNMASK.m)); //negate shift if y<=0
	shift = tfxWideAndNot(tfxWideAnd(absxgreaterthanabsy, xgreaterthan0), shift); //null if abs>absy & x>0

	return tfxWideAdd(atan, shift);
}

tfxINTERNAL inline tfxWideFloat tfxWideCos52s(tfxWideFloat x)
{
	const tfxWideFloat c1 = tfxWideSetSingle(0.9999932946f);
	const tfxWideFloat c2 = tfxWideSetSingle(-0.4999124376f);
	const tfxWideFloat c3 = tfxWideSetSingle(0.0414877472f);
	const tfxWideFloat c4 = tfxWideSetSingle(-0.0012712095f);
	tfxWideFloat x2;      // The input argument squared
	x2 = tfxWideMul(x, x);
	//               (c1+           x2*          (c2+           x2*          (c3+           c4*x2)));
	return tfxWideAdd(c1, tfxWideMul(x2, tfxWideAdd(c2, tfxWideMul(x2, tfxWideAdd(c3, tfxWideMul(c4, x2))))));
}

tfxINTERNAL inline tfxWideFloat tfxWideCos(tfxWideFloat angle) {
	//clamp to the range 0..2pi

	//take absolute value
	angle = tfxWideAndNot(tfxSIGNMASK.m, angle);
	//fmod(angle,twopi)
	angle = tfxWideSub(angle, tfxWideMul(tfxWideConvert(tfxWideConverti(tfxWideMul(angle, tfxWIDEINVTWOPI.m))), tfxWIDEPI2.m)); //simplied SSE2 fmod, must always operate on absolute value
	//if SSE4.1 is always available, comment the line above and uncomment the line below
	//angle=tfxWideSub(angle,tfxWideMul(_mm_floor_ps(tfxWideMul(angle,tfxWideSetSingle(tfxINVTWOPI))),tfxWideSetSingle(tfxPI2))); //faster if SSE4.1 is always available

	tfxWideFloat cosangle = angle;
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWIDEHALFPI.m), tfxWideXOr(cosangle, tfxWideSub(tfxWIDEPI.m, angle))));
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWIDEPI.m), tfxSIGNMASK.m));
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWIDETHREEHALFPI.m), tfxWideXOr(cosangle, tfxWideSub(tfxWIDEPI2.m, angle))));

	tfxWideFloat result = tfxWideCos52s(cosangle);

	result = tfxWideXOr(result, tfxWideAnd(tfxWideAnd(tfxWideGreaterEqual(angle, tfxWIDEHALFPI.m), tfxWideLess(angle, tfxWIDETHREEHALFPI.m)), tfxSIGNMASK.m));
	return result;
}

tfxINTERNAL inline tfxWideFloat tfxWideSin(tfxWideFloat angle) {
	return tfxWideCos(tfxWideSub(tfxWIDEHALFPI.m, angle));
}

tfxINTERNAL inline void tfxWideSinCos(tfxWideFloat angle, tfxWideFloat *sin, tfxWideFloat *cos) {
	tfxWideFloat anglesign = tfxWideOr(tfxWIDEONE.m, tfxWideAnd(tfxSIGNMASK.m, angle));

	//clamp to the range 0..2pi

	//take absolute value
	angle = tfxWideAndNot(tfxSIGNMASK.m, angle);
	//fmod(angle,twopi)
	angle = tfxWideSub(angle, tfxWideMul(tfxWideConvert(tfxWideConverti(tfxWideMul(angle, tfxWIDEINVTWOPI.m))), tfxWIDEPI2.m)); //simplied SSE2 fmod, must always operate on absolute value
	//if SSE4.1 is always available, comment the line above and uncomment the line below
	//angle=tfxWideSub(angle,tfxWideMul(_mm_floor_ps(tfxWideMul(angle,tfxWideSetSingle(tfxINVTWOPI))),tfxWideSetSingle(tfxPI2))); //faster if SSE4.1 is always available

	tfxWideFloat cosangle = angle;
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWIDEHALFPI.m), tfxWideXOr(cosangle, tfxWideSub(tfxWIDEPI.m, angle))));
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWIDEPI.m), tfxSIGNMASK.m));
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWIDETHREEHALFPI.m), tfxWideXOr(cosangle, tfxWideSub(tfxWIDEPI2.m, angle))));

	tfxWideFloat result = tfxWideCos52s(cosangle);

	result = tfxWideXOr(result, tfxWideAnd(tfxWideAnd(tfxWideGreaterEqual(angle, tfxWIDEHALFPI.m), tfxWideLess(angle, tfxWIDETHREEHALFPI.m)), tfxSIGNMASK.m));
	*cos = result;

	tfxWideFloat sinmultiplier = tfxWideMul(anglesign, tfxWideOr(tfxWIDEONE.m, tfxWideAnd(tfxWideGreater(angle, tfxWIDEPI.m), tfxSIGNMASK.m)));
	*sin = tfxWideMul(sinmultiplier, tfxWideFastSqrt(tfxWideSub(tfxWideSetSingle(1.f), tfxWideMul(result, result))));
}

tfxINTERNAL inline void tfxWideSinCosAdd(tfxWideFloat angle, tfxWideFloat *sin, tfxWideFloat *cos) {
	tfxWideFloat anglesign = tfxWideOr(tfxWIDEONE.m, tfxWideAnd(tfxSIGNMASK.m, angle));

	//clamp to the range 0..2pi

	//take absolute value
	angle = tfxWideAndNot(tfxSIGNMASK.m, angle);
	//fmod(angle,twopi)
	angle = tfxWideSub(angle, tfxWideMul(tfxWideConvert(tfxWideConverti(tfxWideMul(angle, tfxWIDEINVTWOPI.m))), tfxWIDEPI2.m)); //simplied SSE2 fmod, must always operate on absolute value
	//if SSE4.1 is always available, comment the line above and uncomment the line below
	//angle=tfxWideSub(angle,tfxWideMul(_mm_floor_ps(tfxWideMul(angle,tfxWideSetSingle(tfxINVTWOPI))),tfxWideSetSingle(tfxPI2))); //faster if SSE4.1 is always available

	tfxWideFloat cosangle = angle;
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWIDEHALFPI.m), tfxWideXOr(cosangle, tfxWideSub(tfxWIDEPI.m, angle))));
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWIDEPI.m), tfxSIGNMASK.m));
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWIDETHREEHALFPI.m), tfxWideXOr(cosangle, tfxWideSub(tfxWIDEPI2.m, angle))));

	tfxWideFloat result = tfxWideCos52s(cosangle);

	result = tfxWideXOr(result, tfxWideAnd(tfxWideAnd(tfxWideGreaterEqual(angle, tfxWIDEHALFPI.m), tfxWideLess(angle, tfxWIDETHREEHALFPI.m)), tfxSIGNMASK.m));
	*cos = result;

	tfxWideFloat sinmultiplier = tfxWideMul(anglesign, tfxWideOr(tfxWIDEONE.m, tfxWideAnd(tfxWideGreater(angle, tfxWIDEPI.m), tfxSIGNMASK.m)));
	*sin = tfxWideAdd(*sin, tfxWideMul(sinmultiplier, tfxWideFastSqrt(tfxWideSub(tfxWIDEONE.m, tfxWideMul(result, result)))));
}
/*
End of Robin Lobel code
*/

//simd mod function thanks to Stephanie Rancourt: http://dss.stephanierct.com/DevBlog/?p=8
tfxINTERNAL inline tfxWideFloat tfxWideMod(const tfxWideFloat a, const tfxWideFloat aDiv) {
	tfxWideFloat c = tfxWideDiv(a, aDiv);
	tfxWideInt i = tfxWideConverti(c);
	tfxWideFloat cTrunc = tfxWideConvert(i);
	tfxWideFloat base = tfxWideMul(cTrunc, aDiv);
	tfxWideFloat r = tfxWideSub(a, base);
	return r;
}

tfxINTERNAL inline tfxWideFloat tfxWideAbs(tfxWideFloat v) {
	const tfxWideFloat tfxABSMASK		  = tfxWideCast(tfxWideSetSinglei(0x7FFFFFFF));
	return tfxWideAnd(v, tfxABSMASK);
}

tfxINTERNAL inline tfxWideFloat tfxWideCopySign(tfxWideFloat dst, tfxWideFloat src) {
	const tfxWideFloat tfxABSMASK		  = tfxWideCast(tfxWideSetSinglei(0x7FFFFFFF));
	tfxWideFloat sign_mask = tfxWideAnd(src, tfxSIGNMASK.m);
	tfxWideFloat abs_mask = tfxWideAnd(dst, tfxABSMASK);
	return tfxWideOr(abs_mask, sign_mask);
}

tfxINTERNAL inline tfxWideInt tfxWideAbsi(tfxWideInt v) {
	return tfxWideAndi(v, tfxABSMASKi.m);
}

//End of SIMD setup

//----------------------------------------------------------
//Section: Enums
//----------------------------------------------------------

//tfx_graph_t presets to determine limits and scales of different graphs, mainly used for the editor
typedef enum {
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
	tfxWeightPreset,
	tfxWeightVariationPreset,
	tfxNoiseOffsetVariationPreset,
	tfxNoiseResolutionPreset,
	tfxSpinPreset,
	tfxSpinVariationPreset,
	tfxDirectionVariationPreset,
	tfxFrameratePreset,
	tfxVelocityTurbulancePreset,
	tfxColorPreset,
	tfxSpinOvertimePreset,
	tfxVelocityOvertimePreset,
	tfxWeightOvertimePreset,
	tfxDirectionOvertimePreset,
	tfxOpacityOvertimePreset,
	tfxPercentOvertime,
	tfxIntensityOvertimePreset,
	tfxGradientMapperOvertimePreset,
	tfxPathDirectionOvertimePreset,
	tfxPathTranslationOvertimePreset,
} tfx_graph_preset;

typedef enum {
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
} tfx_graph_category;

typedef enum tfx_color_ramp_format {
	tfx_color_format_rgba,
	tfx_color_format_bgra
} tfx_color_ramp_format;

//--------------------------------------------
//--Graph_types
//--------------------------------------------

//All the different types of graphs, split into main type: global, property, base, variation and overtime
typedef enum {
	tfxGlobal_life,
	tfxGlobal_amount,
	tfxGlobal_velocity,
	tfxGlobal_noise,
	tfxGlobal_width,
	tfxGlobal_height,
	tfxGlobal_weight,
	tfxGlobal_roll_spin,
	tfxGlobal_pitch_spin,
	tfxGlobal_yaw_spin,
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
	tfxProperty_emitter_width,        //Also used for linear extrusion for paths as well
	tfxProperty_emitter_height,
	tfxProperty_emitter_depth,
	tfxProperty_extrusion,
	tfxProperty_arc_size,
	tfxProperty_arc_offset,

	tfxBase_life,
	tfxBase_amount,
	tfxBase_velocity,
	tfxBase_width,
	tfxBase_height,
	tfxBase_weight,
	tfxBase_pitch_spin,
	tfxBase_yaw_spin,
	tfxBase_roll_spin,
	tfxBase_noise_offset,

	tfxVariation_life,
	tfxVariation_amount,
	tfxVariation_velocity,
	tfxVariation_width,
	tfxVariation_height,
	tfxVariation_weight,
	tfxVariation_path_trajectory_scale,
	tfxVariation_pitch_spin,
	tfxVariation_yaw_spin,
	tfxVariation_roll_spin,
	tfxVariation_noise_offset,
	tfxVariation_noise_resolution,
	tfxVariation_motion_randomness,

	tfxOvertime_red,
	tfxOvertime_green,
	tfxOvertime_blue,
	tfxOvertime_blendfactor,
	tfxOvertime_velocity_adjuster,
	//--These compiled graph values are uploaded to the GPU
	tfxOvertime_intensity,
	tfxOvertime_alpha_sharpness,
	tfxOvertime_curved_alpha,
	tfxOvertime_gradient_mapper,
	tfxOvertime_velocity,
	tfxOvertime_width,
	tfxOvertime_height,
	tfxOvertime_overal_ribbon_scale,
	tfxOvertime_weight,
	tfxOvertime_pitch_spin,
	tfxOvertime_yaw_spin,
	tfxOvertime_roll_spin,
	tfxOvertime_stretch,
	tfxOvertime_velocity_turbulance,
	tfxOvertime_direction_turbulance,
	tfxOvertime_direction,
	tfxOvertime_noise_resolution,
	tfxOvertime_motion_randomness,
	tfxOvertime_uv_offset_y,
	tfxOvertime_uv_scale_y,
	tfxOvertime_clip_offset,
	tfxOvertime_clip_size,

	tfxOverlength_intensity,
	tfxOverlength_alpha_sharpness,
	tfxOverlength_curved_alpha,
	tfxOverlength_gradient_map,
	tfxOverlength_width,
	tfxOverlength_ribbon_fixed_angle,

	tfxFactor_life,
	tfxFactor_size,
	tfxFactor_velocity,
	tfxFactor_intensity,
	//------------------------------------------------------

	tfxTransform_roll,
	tfxTransform_pitch,
	tfxTransform_yaw,
	tfxTransform_translate_x,
	tfxTransform_translate_y,
	tfxTransform_translate_z,
	tfxEmitterGraphMaxIndex,

	tfxGraphMaxIndex
} tfx_graph_type;

#define tfxEffectGraph(graph, index_name) graph.graphs[tfxEffect_##index_name##_index]
#define tfxTransformGraph(graph, index_name) graph.graphs[tfxTransform_##index_name##_index]

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
	tfxEffect_global_overal_scale_index,
	tfxEffect_global_intensity_index,
	tfxEffect_global_splatter_index,
	tfxEffect_global_emitter_width_index,
	tfxEffect_global_emitter_height_index,
	tfxEffect_global_emitter_depth_index,
	tfxEffectGraphs_max_index,
} tfx_global_graph_index;

typedef enum {
	tfxTransform_roll_index,
	tfxTransform_pitch_index,
	tfxTransform_yaw_index,
	tfxTransform_translate_x_index,
	tfxTransform_translate_y_index,
	tfxTransform_translate_z_index,
	tfxTransformGraphs_max_index,
} tfx_transform_graph_index;

typedef enum {
	tfxPath_rotation_range_index,
	tfxPath_rotation_pitch_index,
	tfxPath_rotation_yaw_index,
} tfx_path_graph_index;

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
	tfxEmitter_base_noise_offset_index,

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
	tfxEmitter_variation_noise_offset_index,
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
	tfxRibbon_overtime_width_index,
	tfxRibbon_overtime_overal_scale_index,
	tfxRibbon_overtime_uv_offset_y_index,
	tfxRibbon_overtime_uv_scale_y_index,
	tfxRibbon_overtime_clip_offset_index,
	tfxRibbon_overtime_clip_size_index,

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
	tfxRibbon_overtime_end_index = tfxRibbon_overtime_clip_size_index + 1,
	tfxRibbon_overlength_start = tfxRibbon_overlength_intensity_index,
	tfxRibbon_overlength_end = tfxRibbon_overlength_fixed_angle_index + 1,
} tfx_ribbon_graph_index;

typedef enum {
	tfxGlobal_start = tfxGlobal_life,
	tfxGlobal_end = tfxGlobal_emitter_depth,
	tfxProperty_start = tfxProperty_emission_pitch,
	tfxProperty_end = tfxProperty_arc_offset,
	tfxBase_start = tfxBase_life,
	tfxBase_end = tfxBase_noise_offset,
	tfxVariation_start = tfxVariation_life,
	tfxVariation_end = tfxVariation_motion_randomness,
	tfxOvertime_start = tfxOvertime_red,
	tfxOvertime_end = tfxOvertime_clip_size,
	tfxOvertime_color_start = tfxOvertime_red,
	tfxOvertime_color_end = tfxOvertime_blendfactor,
	tfxOverlength_start = tfxOverlength_intensity,
	tfxOverlength_end = tfxOverlength_ribbon_fixed_angle,
	tfxFactor_start = tfxFactor_life,
	tfxFactor_end = tfxFactor_intensity,
	tfxTransform_start = tfxTransform_roll,
	tfxTransform_end = tfxTransform_translate_z,
	tfxGPU_lookup_start = tfxOvertime_intensity,
	tfxGPU_lookup_end = tfxFactor_intensity,
} tfx_graph_ranges;

typedef enum {
	tfxOscillator_sine,
	tfxOscillator_square,
	tfxOscillator_sawtooth
} tfx_oscillator_type;

//tfx_effect_descriptor_t type - effect contains emitters, and emitters spawn particles, but they both share the same struct for simplicity
typedef enum {
	tfxEffectType,
	tfxEmitterType,
	tfxRibbonType,
	tfxStage,
	tfxFolder,
	tfxMaxDescriptorTypes
} tfx_effect_descriptor_type;

//Different ways that particles can be emitted
typedef enum {
	tfxPoint,
	tfxArea,
	tfxLine,
	tfxEllipse,
	tfxCylinder,
	tfxIcosphere,
	tfxPath,
	tfxOtherEmitter,
	tfxSpawnOnRibbon,
	tfxEmissionTypeMax,
} tfx_emission_type;

typedef enum {
	tfxExtrusionArc,
	tfxExtrusionLinear
} tfx_path_extrusion_type;

//These must not change, values are used in the save file.
typedef enum {
	tfxNoNoise,
	tfxWhiteNoise,
	tfxSimplexNoise,
	tfxCurlNoise,
} tfx_noise_type;

//Determines how for area, line and ellipse emitters the direction that particles should travel when they spawn
//These must not change, values are used in the save file.
typedef enum {
	tfxInwards,
	tfxOutwards,
	tfxBothways,
	tfxSpecified,
	tfxSurface,
	tfxOrbital,
	tfxPathGradient
} tfx_emission_direction;

//For line effects where traverse line is switched on
//These must not change, values are used in the save file.
typedef enum {
	tfxLoop,
	tfxKill,
	tfxLetFree
} tfx_line_traversal_end_behaviour;

//Mainly for the editor, maybe this can just be moved there instead?
typedef enum {
	tfxFullColor,
	tfxOneColor,
	tfxGreyScale
} tfx_export_color_options;

//Mainly for the editor, maybe this can just be moved there instead?
typedef enum {
	tfxSpriteSheet,
	tfxStrip,
	tfxSeparateImages
} tfx_export_options;

//tfx_graph_t data can be looked up in one of 2 ways, either by just using linear/bezier interpolation (slower), or look up the value in a pre-compiled look up table.
typedef enum {
	tfxPrecise,
	tfxFast
} tfx_lookup_mode;

typedef enum {
	tfxCalculateFrames,
	tfxBakeSpriteData,
	tfxLinkUpSprites,
	tfxLinkUpRibbons,
	tfxCompressSpriteFrames,
	tfxCompressRibbonFrames,
	tfxBakingDone
} tfx_record_progress;

//Used in file loading - for loading effects library
typedef enum {
	tfxString,
	tfxSInt,
	tfxUInt,
	tfxFloat,
	tfxDouble,
	tfxBool,
	tfxColor,
	tfxUInt64,
	tfxFloat3,
	tfxFloat2,
	tfxAttributeGraph,
	tfxTransformGraph,
	tfxGraphProperty,
} tfx_data_type;

//Block designators for loading effects library and other files like animation sprite data
//The values of existing enums below must never change or older files won't load anymore!
typedef enum {
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
	tfxStartRibbonEmitter,
	tfxEndRibbonEmitter,
	tfxStartGraphProperties,
	tfxEndGraphProperties,
	tfxStartPathNodes,
	tfxEndPathNodes,
} tfx_effect_library_stream_context;

typedef enum {
	//Single object handles
	tfx_effect_library_mb,
	tfx_effect_manager_mb,
	//Effect Manager storage lists
	tfx_particle_array_buffers_mb,
	tfx_particle_arrays_mb,
	tfx_particle_location_buffers_mb,
	tfx_particle_location_arrays_mb,
	tfx_free_particle_lists_mb,
	tfx_free_particle_location_lists_mb,
	tfx_free_ribbon_segment_lists_mb,
	tfx_ribbon_segment_buckets_mb,
	tfx_sorting_work_entry_mb,
	tfx_spawn_work_mb,
	tfx_ribbon_work_mb,
	tfx_control_work_mb,
	tfx_ribbon_control_work_mb,
	tfx_age_work_mb,
	tfx_particle_indexes_mb,
	tfx_free_particle_indexes_mb,
	tfx_effects_in_use_mb,
	tfx_control_emitter_queue_mb,
	tfx_emitters_check_capture_mb,
	tfx_free_effects_mb,
	tfx_free_emitters_mb,
	tfx_free_gpu_emitters_mb,
	tfx_free_ribbon_emitters_mb,
	tfx_free_path_quaternions_mb,
	tfx_path_quaternions_mb,
	tfx_effects_mb,
	tfx_emitters_mb,
	tfx_gpu_emitters_mb,
	tfx_ribbon_emitters_mb,
	tfx_deffered_spawn_work_mb,
	tfx_deffered_ribbon_spawn_work_mb,
	tfx_unique_sprite_ids_mb,
	tfx_free_compute_controllers_mb,
	//Effect Library Storage lists
} tfx_memory_allocation_type;

// -- [Bit_fields]
typedef tfxU32 tfxParticleEmitterFlags;			//tfx_particle_emitter_flag_bits
typedef tfxU32 tfxRibbonEmitterFlags;			//tfx_ribbon_emitter_flag_bits
typedef tfxU32 tfxSharedEmitterFlags;			//tfx_shared_emitter_flag_bits
typedef tfxU32 tfxColorRampFlags;				//tfx_color_ramp_flag_bits
typedef tfxU32 tfxGraphFlags;			        //tfx_graph_flag_bits
typedef tfxU32 tfxEffectPropertyFlags;          //tfx_effect_property_flag_bits
typedef tfxU32 tfxParticleFlags;                 //tfx_particle_flag_bits
typedef tfxU32 tfxEmitterStateFlags;            //tfx_emitter_state_flag_bits
typedef tfxU32 tfxRibbonEmitterStateFlags;      //tfx_ribbon_emitter_state_flag_bits
typedef tfxU32 tfxRibbonFlags;		            //tfx_ribbon_flag_bits
typedef tfxU32 tfxRibbonBucketFlags;            //tfx_ribbon_bucket_flag_bits
typedef tfxU32 tfxRibbonBucketComputeShaderType;//tfx_ribbon_compute_shader_type
typedef tfxU32 tfxEffectStateFlags;             //tfx_effect_state_flag_bits
typedef tfxU32 tfxParticleControlFlags;         //tfx_particle_control_flag_bits
typedef tfxU32 tfxAttributeNodeFlags;           //tfx_attribute_node_flag_bits
typedef tfxU32 tfxAngleSettingFlags;            //tfx_angle_setting_flag_bits
typedef tfxU32 tfxEffectManagerFlags;			//tfx_effect_manager_flag_bits
typedef tfxU32 tfxErrorFlags;                   //tfx_error_flag_bits
typedef tfxU32 tfxEffectCloningFlags;           //tfx_effect_cloning_flag_bits
typedef tfxU32 tfxAnimationFlags;               //tfx_animation_flag_bits
typedef tfxU32 tfxAnimationInstanceFlags;       //tfx_animation_instance_flag_bits
typedef tfxU32 tfxAnimationManagerFlags;        //tfx_animation_manager_flag_bits
typedef tfxU32 tfxEmitterPathFlags;             //tfx_emitter_path_flag_bits
typedef tfxU32 tfxEmitterControlProfileFlags;   //tfx_emitter_control_profile_flag_bits
typedef tfxU32 tfxContextPolicyFlags;			//tfx_context_policy_flag_bits
typedef tfxU32 tfxPackageFlags;                 //tfx_package_flag_bits

typedef enum {
	tfxSharedFlag_capture_after_transform						= 1 << 8,
	tfxSharedFlag_remove										= 1 << 9
} tfx_shared_flag_bits;

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
	tfxErrorCode_library_loaded_without_shape_loader            = 1 << 13,
	tfxErrorCode_library_object_could_not_be_created            = 1 << 14
} tfx_error_flag_bits;

typedef enum {
	tfxPackageFlags_none                                        = 0,
	tfxPackageFlags_loaded_from_memory                          = 1
} tfx_package_flag_bits;

typedef enum {
	tfxEffectCloningFlags_none                                  = 0,
	tfxEffectCloningFlags_keep_user_data                        = 1 << 0,
	tfxEffectCloningFlags_clone_graphs                          = 1 << 1,
	tfxEffectCloningFlags_history                               = 1 << 2,
	tfxEffectCloningFlags_clone_camera_settings                 = 1 << 3,
	tfxEffectCloningFlags_camera_and_graphs		                = tfxEffectCloningFlags_clone_graphs | tfxEffectCloningFlags_clone_camera_settings 
} tfx_effect_cloning_flag_bits;

typedef enum {
	tfxEffectManagerSetup_none,
	tfxEffectManagerSetup_group_sprites_by_effect,
} tfx_effect_manager_setup;

typedef enum {
	tfxBillboarding_align_to_camera                             = 0,            //Align to Camera only
	tfxBillboarding_free_align                                  = 1,            //Free align
	tfxBillboarding_align_to_camera_and_vector                  = 2,            //Align to camera and vector
	tfxBillboarding_align_to_vector                             = 3,            //Align to vector
	tfxBillboarding_max                                         = 4
} tfx_billboarding_option;

typedef enum {
	tfxEmitterControlProfile_basic                              = 0,
	tfxEmitterControlProfile_simplex_noise                      = 1 << 0,
	tfxEmitterControlProfile_curl_noise                         = 1 << 1,
	tfxEmitterControlProfile_motion_randomness                  = 1 << 2,
	tfxEmitterControlProfile_orbital                            = 1 << 3,
	tfxEmitterControlProfile_path                               = 1 << 4,
	tfxEmitterControlProfile_rotated_path						= 1 << 5,
	tfxEmitterControlProfile_edge_traversal                     = 1 << 6,
	tfxEmitterControlProfile_edge_kill                          = 1 << 7,
	tfxEmitterControlProfile_edge_loop                          = 1 << 8,
	tfxEmitterControlProfile_stretch                            = 1 << 9,
	tfxEmitterControlProfile_other_ribbon_emitter_path          = 1 << 10,
	tfxEmitterControlProfile_spin3d                             = 1 << 11,
	tfxEmitterControlProfile_spin                               = 1 << 12,
	tfxEmitterControlProfile_trajectory                         = 1 << 13,
	tfxEmitterControlProfile_line		                        = 1 << 14,
	tfxEmitterControlProfile_rotated_line                       = 1 << 15,
	tfxEmitterControlProfile_has_simplex_noise_type				= tfxEmitterControlProfile_simplex_noise | tfxEmitterControlProfile_curl_noise,
	tfxEmitterControlProfile_has_any_noise						= tfxEmitterControlProfile_simplex_noise | tfxEmitterControlProfile_curl_noise | tfxEmitterControlProfile_motion_randomness,
	tfxEmitterControlProfile_has_rotated_path_or_line			= tfxEmitterControlProfile_rotated_path | tfxEmitterControlProfile_trajectory | tfxEmitterControlProfile_rotated_line,
	tfxEmitterControlProfile_any_line							= tfxEmitterControlProfile_line | tfxEmitterControlProfile_rotated_line,
} tfx_emitter_control_profile_flag_bits;

typedef enum {
	tfx_ctx_policy_flag_none                                	= 0,
	tfx_ctx_policy_flag_velocity_is_bezier_graph            	= 1,
	tfx_ctx_policy_flag_velocity_has_oscillator             	= 1 << 1,
	tfx_ctx_policy_flag_weight_is_bezier_graph              	= 1 << 2,
	tfx_ctx_policy_flag_weight_has_oscillator               	= 1 << 3,
	tfx_ctx_policy_flag_velocity_turbulance_is_bezier_graph 	= 1 << 4,
	tfx_ctx_policy_flag_velocity_turbulance_has_oscillator  	= 1 << 5,
	tfx_ctx_policy_flag_noise_resolution_is_bezier_graph    	= 1 << 6,
	tfx_ctx_policy_flag_noise_resolution_has_oscillator     	= 1 << 7,
	tfx_ctx_policy_flag_motion_randomness_is_bezier_graph   	= 1 << 8,
	tfx_ctx_policy_flag_motion_randomness_has_oscillator    	= 1 << 9,
	tfx_ctx_policy_flag_stretch_is_bezier_graph              	= 1 << 10,
	tfx_ctx_policy_flag_stretch_has_oscillator               	= 1 << 11,
	tfx_ctx_policy_flag_is_ordered				              	= 1 << 12,
	tfx_ctx_policy_flag_transform_relative						= 1 << 13,
	tfx_ctx_policy_flag_direction_is_bezier_graph            	= 1 << 14,
	tfx_ctx_policy_flag_direction_has_oscillator             	= 1 << 15,
} tfx_context_policy_flag_bits;

typedef enum {
	tfxEffectManagerFlags_none                              	= 0,
	tfxEffectManagerFlags_disable_spawning                  	= 1,
	tfxEffectManagerFlags_force_capture                     	= 2,            //Unused
	tfxEffectManagerFlags_use_compute_shader                	= 1 << 3,
	tfxEffectManagerFlags_update_base_values                	= 1 << 6,
	tfxEffectManagerFlags_dynamic_sprite_allocation         	= 1 << 7,
	tfxEffectManagerFlags_animation_loops                   	= 1 << 9,
	tfxEffectManagerFlags_update_age_only                   	= 1 << 11,
	tfxEffectManagerFlags_single_threaded                   	= 1 << 12,
	tfxEffectManagerFlags_double_buffer_sprites             	= 1 << 13,
	tfxEffectManagerFlags_recording_sprites                 	= 1 << 14,
	tfxEffectManagerFlags_using_uids                        	= 1 << 15,
	tfxEffectManagerFlags_update_bounding_boxes             	= 1 << 17,
	tfxEffectManagerFlags_auto_order_effects                	= 1 << 19,
	tfxEffectManagerFlags_direct_to_staging_buffer          	= 1 << 20,
	tfxEffectManagerFlags_has_ribbons_to_draw               	= 1 << 21,
	tfxEffectManagerFlags_record_with_compute_image_index   	= 1 << 22,
} tfx_effect_manager_flag_bits;

//These values must stay the same
typedef enum {
	tfxVectorAlignType_motion,
	tfxVectorAlignType_emission,
	tfxVectorAlignType_emitter,
	tfxVectorAlignType_velocity,
	tfxVectorAlignType_max,
} tfx_vector_align_type;

typedef enum {
	tfxPathFlags_none,
	tfxPathFlags_mode_origin                                    = 1 << 1,
	tfxPathFlags_mode_node                                      = 1 << 2,
	tfxPathFlags_reverse_direction                              = 1 << 4,
	tfxPathFlags_rotation_range_yaw_only                        = 1 << 5
} tfx_emitter_path_flag_bits;

                                                                                //Particle property that defines how a particle will rotate
typedef enum {
	tfxAngleSettingFlags_none                                   = 0,            //No flag
	tfxAngleSettingFlags_align_roll                             = 1 << 0,       //Align the particle with it's direction of travel in
	tfxAngleSettingFlags_random_roll                            = 1 << 1,       //Chose a random angle at spawn time/state_flags
	tfxAngleSettingFlags_specify_roll                           = 1 << 2,       //Specify the angle at spawn time
	tfxAngleSettingFlags_align_with_emission                    = 1 << 3,       //Align the particle with the emission direction only
	tfxAngleSettingFlags_random_pitch                           = 1 << 4,       //3d mode allows for rotating pitch and yaw when not using billboarding (when particle always faces the camera)
	tfxAngleSettingFlags_random_yaw                             = 1 << 5,
	tfxAngleSettingFlags_specify_pitch                          = 1 << 6,
	tfxAngleSettingFlags_specify_yaw                            = 1 << 7
} tfx_angle_setting_flag_bits;

                                                                                //All the state_flags needed by the ControlParticle function put into one typedef enum save typedef enum
typedef enum {
	tfxParticleControlFlags_none                                = 0,
	tfxParticleControlFlags_remove                              = tfxSharedFlag_remove,
	tfxParticleControlFlags_relative_position                   = 1 << 1,
	tfxParticleControlFlags_relative_angle                      = 1 << 2,
	tfxParticleControlFlags_point                               = 1 << 3,
	tfxParticleControlFlags_area                                = 1 << 4,
	tfxParticleControlFlags_line                                = 1 << 5,
	tfxParticleControlFlags_ellipse                             = 1 << 6,
	tfxParticleControlFlags_loop                                = 1 << 7,
	tfxParticleControlFlags_kill                                = 1 << 8,
	tfxParticleControlFlags_letFree                             = 1 << 9,
	tfxParticleControlFlags_edge_traversal                      = 1 << 10,
	tfxParticleControlFlags_random_color                        = 1 << 11,
	tfxParticleControlFlags_base_uniform_size                   = 1 << 12,
	tfxParticleControlFlags_lifetime_uniform_size               = 1 << 13,
	tfxParticleControlFlags_animate                             = 1 << 14,
	tfxParticleControlFlags_reverse_animation                   = 1 << 15,
	tfxParticleControlFlags_play_once                           = 1 << 16,
	tfxParticleControlFlags_align                               = 1 << 17,
	tfxParticleControlFlags_emission                            = 1 << 18,
	tfxParticleControlFlags_random_roll                         = 1 << 19,
	tfxParticleControlFlags_specify_roll                        = 1 << 20,
	tfxParticleControlFlags_random_pitch                        = 1 << 21,
	tfxParticleControlFlags_specify_pitch                       = 1 << 22,
	tfxParticleControlFlags_random_yaw                          = 1 << 23,
	tfxParticleControlFlags_specify_yaw                         = 1 << 24,
} tfx_particle_control_flag_bits;

typedef enum {
	tfxEffectPropertyFlags_none                                 = 0,
	tfxEffectPropertyFlags_depth_draw_order                     = 1 << 1,
	tfxEffectPropertyFlags_guaranteed_order                     = 1 << 2,
	tfxEffectPropertyFlags_age_order                            = 1 << 3,
	tfxEffectPropertyFlags_use_keyframes                        = 1 << 4,
	tfxEffectPropertyFlags_include_in_sprite_data_export        = 1 << 5,		//In the editor you can specify which effects you want to be included in a spritedata export
	tfxEffectPropertyFlags_global_uniform_size                  = 1 << 6,		//Keep the global particle size uniform
	tfxEffectPropertyFlags_is_in_folder                         = 1 << 7,		//This effect is located inside a folder. 
	tfxEffectPropertyFlags_history_effect					    = 1 << 12,		//Flagged if the effect is just a change in the editor
	tfxEffectPropertyFlags_is_ordered						    = tfxEffectPropertyFlags_depth_draw_order | tfxEffectPropertyFlags_age_order,
} tfx_effect_property_flag_bits;

typedef enum {
	tfxEmitterPropertyFlags_none							    = 0,
	tfxEmitterPropertyFlags_relative_angle					    = 1 << 1,       //Keep the angle of the particles relative to the current angle of the emitter
	tfxEmitterPropertyFlags_image_handle_auto_center		    = 1 << 2,		//Set the offset of the particle to the center of the image
	tfxEmitterPropertyFlags_edge_traversal					    = 1 << 3,       //Line and Path emitters only: make particles traverse the line/path
	tfxEmitterPropertyFlags_base_uniform_size				    = 1 << 4,       //Keep the base particle size uniform
	tfxEmitterPropertyFlags_lifetime_uniform_size			    = 1 << 5,		//Keep the size over lifetime of the particle uniform
	tfxEmitterPropertyFlags_wrap_single_sprite				    = 1 << 6,		//When recording sprite data, single particles can have their invalid capured index set to the current frame for better looping
	tfxEmitterPropertyFlags_use_spawn_ratio					    = 1 << 7,       //Option for area emitters to multiply the amount spawned by a ration of particles per pixels squared
	tfxEmitterPropertyFlags_area_open_ends					    = 1 << 8,       //Only sides of the area/cylinder are spawned on when fill area is not checked
	tfxEmitterPropertyFlags_match_amount_to_grid_points		    = 1 << 9,		//Match the amount to spawn with a single emitter to the number of grid points in the effect
	tfxEmitterPropertyFlags_run_on_gpu						    = 1 << 10,		//Makes this emitter a custom GPU emitter
	tfxEmitterPropertyFlags_alt_velocity_lifetime_sampling	    = 1 << 11,		//The point on the path dictates where on the velocity overtime graph that the particle should sample from rather then the age of the particle
	tfxEmitterPropertyFlags_alt_color_lifetime_sampling		    = 1 << 12,		//The point on the path dictates where on the color overtime graph that the particle should sample from rather then the age of the particle
	tfxEmitterPropertyFlags_alt_size_lifetime_sampling		    = 1 << 13,		//The point on the path dictates where on the size overtime graph that the particle should sample from rather then the age of the particle
	tfxEmitterPropertyFlags_use_path_as_trajectory				= 1 << 15,		//When using path emission type, all particles will be spawned at the start of the path only and travel along the path. Only available when traverse edge is active
} tfx_particle_emitter_flag_bits;

typedef enum {
	tfxSharedEmitterPropertyFlags_none                          = 0,
	tfxSharedEmitterPropertyFlags_random_color					= 1 << 0,       //Pick a random color from the color overtime gradient rather then change the color over the lifetime of the particle
	tfxSharedEmitterPropertyFlags_relative_position				= 1 << 1,       //Keep the particles position relative to the current position of the emitter
	tfxSharedEmitterPropertyFlags_single						= 1 << 2,       //Only spawn a single particle (or number of particles specified by spawn_amount) that does not expire
	tfxSharedEmitterPropertyFlags_spawn_on_grid					= 1 << 3,       //When using an area, line or ellipse emitter, spawn along a grid
	tfxSharedEmitterPropertyFlags_grid_spawn_clockwise			= 1 << 4,	    //Spawn clockwise/left to right around the area
	tfxSharedEmitterPropertyFlags_fill_area						= 1 << 5,       //Fill the area
	tfxSharedEmitterPropertyFlags_emitter_handle_auto_center	= 1 << 6,		//Center the handle of the emitter
	tfxSharedEmitterPropertyFlags_animate						= 1 << 7,       //Animate the particle shape if it has more than one frame of animation
	tfxSharedEmitterPropertyFlags_reverse_animation				= 1 << 8,       //Make the image animation go in reverse
	tfxSharedEmitterPropertyFlags_play_once						= 1 << 9,       //Play the animation once only
	tfxSharedEmitterPropertyFlags_random_start_frame			= 1 << 10,      //Start the animation of the image from a random frame
	tfxSharedEmitterPropertyFlags_grid_spawn_random				= 1 << 12,		//Spawn on grid points but randomly rather then in sequence
	tfxSharedEmitterPropertyFlags_spawn_location_source			= 1 << 13,	    //This emitter is the source for another emitter that uses it to spawn particles at the location of this emitters' particles
	tfxSharedEmitterPropertyFlags_enabled						= 1 << 14,      //The emitter is enabled or not, meaning it will or will not be added the effect manager with tfx__add_effect
	tfxSharedEmitterPropertyFlags_use_color_hint				= 1 << 15,	    //Activate a second color to tint the particles and mix between the two colors.
	tfxSharedEmitterPropertyFlags_is_in_folder					= 1 << 16,
	tfxSharedEmitterPropertyFlags_exclude_from_hue_adjustments	= 1 << 17,		//Emitter will be excluded from effect hue adjustments if this flag is checked
	tfxSharedEmitterPropertyFlags_hidden						= 1 << 18,		//Flagged when hidden from showing in the editor. This is mainly used in the undo/history system.
	tfxSharedEmitterPropertyFlags_do_not_render					= 1 << 19,		//particles will be processed but their scale will be set to 0 so that they're not rendered.
} tfx_shared_emitter_flag_bits;

typedef enum {
	tfxRibbonPropertyFlags_none                                 = 0,
	tfxRibbonPropertyFlags_use_path_from_another_emitter	    = 1 << 0,
	tfxRibbonPropertyFlags_static							    = 1 << 1,
	tfxRibbonPropertyFlags_always_face_camera				    = 1 << 2,
	tfxRibbonPropertyFlags_frenet_serret_frame				    = 1 << 3,
	tfxRibbonPropertyFlags_fixed_angle						    = 1 << 4,
} tfx_ribbon_emitter_flag_bits;

typedef enum {
	tfxColorRampFlags_none                                      = 0,
	tfxColorRampFlags_use_sinusoidal_ramp_generation            = 1 << 0,		//Use this flag to toggle between sinusoidal color ramp generation
} tfx_color_ramp_flag_bits;

typedef enum {
	tfxColorInterpolation_linear_srgb = 0,
	tfxColorInterpolation_oklch,
	tfxColorInterpolation_hsl,
	tfxColorInterpolation_linear_rgb,
	tfxColorInterpolation_max
} tfx_color_interpolation_mode;

typedef enum {
	tfxGraphFlags_none                                          = 0,
	tfxGraphFlags_use_bezier_sampling                           = 1 << 0,	
	tfxGraphFlags_enable_oscillator                             = 1 << 1,	
	tfxGraphFlags_multi_node_graph                              = 1 << 2,	
} tfx_graph_flag_bits;

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
	tfxParticleFlags_none                                       = 0,
	tfxParticleFlags_capture_after_transform                    = tfxSharedFlag_capture_after_transform,    //Particle will be captured after a transfrom, used for traversing lines and looping back to the beginning to avoid lerping imbetween
	tfxParticleFlags_remove                                     = tfxSharedFlag_remove						//Particle will be removed this or next frame
} tfx_particle_flag_bits;

typedef enum {
	tfxEmitterStateFlags_none                                   = 0,
	tfxEmitterStateFlags_stop_spawning                          = 1 << 3,       //Tells the emitter to stop spawning
	tfxEmitterStateFlags_remove                                 = 1 << 4,       //Tells the effect/emitter to remove itself from the effect manager immediately
	tfxEmitterStateFlags_unused1                                = 1 << 5,       //the emitter is enabled. **moved to property state_flags**
	tfxEmitterStateFlags_retain_matrix                          = 1 << 6,       //Internal flag about matrix usage
	tfxEmitterStateFlags_no_tween_this_update                   = 1 << 7,       //Internal flag generally, but you could use it if you want to teleport the effect to another location
	tfxEmitterStateFlags_is_single                              = 1 << 8,
	tfxEmitterStateFlags_not_line                               = 1 << 9,
	tfxEmitterStateFlags_can_spin                               = 1 << 12,
	tfxEmitterStateFlags_is_edge_traversal                      = 1 << 13,
	tfxEmitterStateFlags_loop                                   = 1 << 15,
	tfxEmitterStateFlags_kill                                   = 1 << 16,
	tfxEmitterStateFlags_single_shot_done                       = 1 << 17,
	tfxEmitterStateFlags_is_line_loop_or_kill                   = 1 << 18,
	tfxEmitterStateFlags_is_area                                = 1 << 19,
	tfxEmitterStateFlags_no_tween                               = 1 << 20,
	tfxEmitterStateFlags_align_with_velocity                    = 1 << 21,
	tfxEmitterStateFlags_can_spin_pitch_and_yaw                 = 1 << 24,      //For 3d emitters that have free alignment and not always facing the camera
	tfxEmitterStateFlags_has_path                               = 1 << 25,
	tfxEmitterStateFlags_has_rotated_path                       = 1 << 27,
	tfxEmitterStateFlags_max_active_paths_reached               = 1 << 28,
	tfxEmitterStateFlags_is_in_ordered_effect                   = 1 << 29,
	tfxEmitterStateFlags_wrap_single_sprite                     = 1 << 30,
	tfxEmitterStateFlags_src_ribbon_is_also_relative	        = 1 << 31,	    //Flagged when emission type is spawn on ribbon and both the emitter and ribbon emitter are relative.
} tfx_emitter_state_flag_bits;

typedef enum {
	tfxRibbonEmitterStateFlags_none                             = 0,
	tfxRibbonEmitterStateFlags_stop_spawning                    = 1 << 3,       //Tells the emitter to stop spawning
	tfxRibbonEmitterStateFlags_remove                           = 1 << 4,       //Tells the effect/emitter to remove itself from the effect manager immediately
	tfxRibbonEmitterStateFlags_single_shot_done                 = 1 << 17,
} tfx_ribbon_emitter_state_flag_bits;

typedef enum {
	tfxRibbonFlags_none		                                    = 0,
	tfxRibbonFlags_active	                                    = 1 << 0,
	tfxRibbonFlags_relative                                     = 1 << 1
} tfx_ribbon_flag_bits;

typedef enum {
	tfxRibbonBucketFlags_none                                   = 0,
	tfxRibbonBucketFlags_initialised                            = 1 << 1 
} tfx_ribbon_bucket_flag_bits;

typedef enum {
	tfxRibbonShader_always_face_camera,
	tfxRibbonShader_fixed_angle,
	tfxRibbonShader_frenet_serret_frame,
} tfx_ribbon_compute_shader_type;

typedef enum {
	tfxEffectStateFlags_none                                    = 0,
	tfxEffectStateFlags_stop_spawning                           = 1 << 3,       //Tells the emitter to stop spawning
	tfxEffectStateFlags_remove                                  = 1 << 4,       //Tells the effect/emitter to remove itself from the effect manager immediately
	tfxEffectStateFlags_retain_matrix                           = 1 << 6,       //Internal flag about matrix usage
	tfxEffectStateFlags_no_tween_this_update                    = 1 << 7,       //Internal flag generally, but you could use it if you want to teleport the effect to another location
	tfxEffectStateFlags_override_overal_scale                   = 1 << 8,       //Flagged when the over scale is overridden with tfx_SetEffectOveralScale
	tfxEffectStateFlags_override_orientiation                   = 1 << 9,       //Flagged when any of the effect angles are overridden
	tfxEffectStateFlags_override_size_multiplier                = 1 << 10,      //Flagged when any of the effect size multipliers are overridden
	tfxEffectStateFlags_no_tween                                = 1 << 20
} tfx_effect_state_flag_bits;

typedef enum {
	tfxAttributeNodeFlags_none                                  = 0,
	tfxAttributeNodeFlags_is_curve                              = 1 << 0,
	tfxAttributeNodeFlags_is_left_curve                         = 1 << 1,
	tfxAttributeNodeFlags_is_right_curve                        = 1 << 2,
	tfxAttributeNodeFlags_curves_initialised                    = 1 << 3,
	tfxAttributeNodeFlags_path_node_accumulate                  = 1 << 4,
} tfx_attribute_node_flag_bits;

typedef enum {
	tfxAnimationFlags_none                                      = 0,
	tfxAnimationFlags_loop                                      = 1 << 0,
	tfxAnimationFlags_seamless                                  = 1 << 1,
	tfxAnimationFlags_needs_recording                           = 1 << 2,
	tfxAnimationFlags_export_with_transparency                  = 1 << 3,
	tfxAnimationFlags_auto_set_length                           = 1 << 4,
} tfx_animation_flag_bits;

typedef enum {
	tfx_perspective_view            							= 1,
	tfx_isometric_view              							= 1 << 1,
	tfx_flat_view                   							= 1 << 2,
	tfx_3d_views                    							= tfx_isometric_view | tfx_perspective_view
} tfx_render_view_mode;

typedef enum {
	tfxAnimationInstanceFlags_none                              = 0,
	tfxAnimationInstanceFlags_loop                              = 1 << 0,
} tfx_animation_instance_flag_bits;

typedef enum {
	tfxAnimationManagerFlags_none                               = 0,
	tfxAnimationManagerFlags_has_animated_shapes                = 1 << 0,
	tfxAnimationManagerFlags_initialised                        = 1 << 1,
} tfx_animation_manager_flag_bits;

typedef enum {
	tfxPresetPath_vline  			                            = 0,
	tfxPresetPath_hline,
	tfxPresetPath_ring,
	tfxPresetPath_spiral720,
	tfxPresetPath_spiral_ball,
	tfxPresetPath_whirlpool,
	tfxPresetPath_triangle,
	tfxPresetPath_loop,
	tfxPresetPath_whispy_loop,
	tfxPresetPath_whispy_loops,
	tfxPresetPath_infinity_loop1,
	tfxPresetPath_infinity_loop2,
	tfxPresetPath_fly_path,
	tfxPresetPath_meander,
	tfxPresetPath_max_presets
} tfx_preset_path;

//Array_sizes
typedef enum {
	tfxArraySize_segment_bucket = tfxMAX_SEGMENT_COUNT / tfxMIN_SEGMENT_COUNT,
} tfx_array_sizes;

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
const tfxWideArray tfx180RadiansWide = tfxWideSetConst(3.14159f);
const float tfxGAMMA = 1.f;
#if defined(__x86_64__) || defined(_M_X64)
typedef tfxU64 tfxAddress;
#else
typedef tfxU32 tfxAddress;
#endif
#define tfxSPRITE_SIZE_SSCALE	0.25000762962736f
#define tfxSPRITE_HANDLE_SSCALE 0.00390636921292f

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

//-----------------------------------------------------------
//Section: Preset_spline_paths
//-----------------------------------------------------------
#define tfx__array_size(array) sizeof(array) / sizeof(array[0])

const float tfx__path_preset_vline[96] = {
	0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.250000f, 0.000000f, 0.000000f, 0.500000f, 0.000000f, 0.000000f, 0.750000f, 0.000000f, 0.000000f, 1.000000f, 0.000000f, 0.000000f, 1.250000f, 0.000000f, 0.000000f, 1.500000f, 0.000000f, 0.000000f, 1.750000f, 0.000000f, 
	0.000000f, 2.000000f, 0.000000f, 0.000000f, 2.250000f, 0.000000f, 0.000000f, 2.500000f, 0.000000f, 0.000000f, 2.750000f, 0.000000f, 0.000000f, 3.000000f, 0.000000f, 0.000000f, 3.250000f, 0.000000f, 0.000000f, 3.500000f, 0.000000f, 0.000000f, 3.750000f, 0.000000f, 
	0.000000f, 4.000000f, 0.000000f, 0.000000f, 4.250000f, 0.000000f, 0.000000f, 4.500000f, 0.000000f, 0.000000f, 4.750000f, 0.000000f, 0.000000f, 5.000000f, 0.000000f, 0.000000f, 5.250000f, 0.000000f, 0.000000f, 5.500000f, 0.000000f, 0.000000f, 5.750000f, 0.000000f, 
	0.000000f, 6.000000f, 0.000000f, 0.000000f, 6.250000f, 0.000000f, 0.000000f, 6.500000f, 0.000000f, 0.000000f, 6.750000f, 0.000000f, 0.000000f, 7.000000f, 0.000000f, 0.000000f, 7.250000f, 0.000000f, 0.000000f, 7.500000f, 0.000000f, 0.000000f, 7.750000f, 0.000000f
};

const float tfx__path_preset_hline[96] = {
	0.000000f, 0.000000f, 0.000000f, 0.250000f, 0.000000f, 0.000000f, 0.500000f, 0.000000f, 0.000000f, 0.750000f, 0.000000f, 0.000000f, 1.000000f, 0.000000f, 0.000000f, 1.250000f, 0.000000f, 0.000000f, 1.500000f, 0.000000f, 0.000000f, 1.750000f, 0.000000f, 0.000000f, 
	2.000000f, 0.000000f, 0.000000f, 2.250000f, 0.000000f, 0.000000f, 2.500000f, 0.000000f, 0.000000f, 2.750000f, 0.000000f, 0.000000f, 3.000000f, 0.000000f, 0.000000f, 3.250000f, 0.000000f, 0.000000f, 3.500000f, 0.000000f, 0.000000f, 3.750000f, 0.000000f, 0.000000f, 
	4.000000f, 0.000000f, 0.000000f, 4.250000f, 0.000000f, 0.000000f, 4.500000f, 0.000000f, 0.000000f, 4.750000f, 0.000000f, 0.000000f, 5.000000f, 0.000000f, 0.000000f, 5.250000f, 0.000000f, 0.000000f, 5.500000f, 0.000000f, 0.000000f, 5.750000f, 0.000000f, 0.000000f, 
	6.000000f, 0.000000f, 0.000000f, 6.250000f, 0.000000f, 0.000000f, 6.500000f, 0.000000f, 0.000000f, 6.750000f, 0.000000f, 0.000000f, 7.000000f, 0.000000f, 0.000000f, 7.250000f, 0.000000f, 0.000000f, 7.500000f, 0.000000f, 0.000000f, 7.750000f, 0.000000f, 0.000000f
};

const float tfx__path_preset_ring[102] = {
	0.000000f, 0.000000f, 0.000000f, 0.000000f, 0.496584f, 0.050498f, 0.000000f, 0.972837f, 0.199923f, 0.000000f, 1.409262f, 0.442159f, 0.000000f, 1.787992f, 0.767288f, 0.000000f, 2.093522f, 1.161999f, 0.000000f, 2.313342f, 1.610134f, 0.000000f, 2.438454f, 2.093344f, 
	0.000000f, 2.463735f, 2.591848f, 0.000000f, 2.388151f, 3.085237f, 0.000000f, 2.214795f, 3.553311f, 0.000000f, 1.950766f, 3.976907f, 0.000000f, 1.606872f, 4.338683f, 0.000000f, 1.197192f, 4.623828f, 0.000000f, 0.738499f, 4.820669f, 0.000000f, 0.249572f, 4.921146f, 
	0.000000f, -0.249572f, 4.921146f, 0.000000f, -0.738500f, 4.820669f, 0.000000f, -1.197192f, 4.623828f, 0.000000f, -1.606872f, 4.338683f, 0.000000f, -1.950766f, 3.976906f, 0.000000f, -2.214796f, 3.553310f, 0.000000f, -2.388151f, 3.085236f, 0.000000f, -2.463735f, 2.591847f, 
	0.000000f, -2.438454f, 2.093344f, 0.000000f, -2.313342f, 1.610133f, 0.000000f, -2.093522f, 1.162000f, 0.000000f, -1.787992f, 0.767288f, 0.000000f, -1.409262f, 0.442159f, 0.000000f, -0.972837f, 0.199923f, 0.000000f, -0.496583f, 0.050498f, 0.000000f, 0.000000f, 0.000000f, 
	0.000000f, 0.496584f, 0.050498f, 0.000000f, 0.972837f, 0.199923f
};

const float tfx__path_preset_spiral720[96] = {
	2.000000f, 0.000000f, 0.000000f, 1.837916f, 0.258065f, 0.788712f, 1.377934f, 0.516129f, 1.449586f, 0.694610f, 0.774194f, 1.875504f, -0.101298f, 1.032258f, 1.997433f, -0.880788f, 1.290323f, 1.795609f, -1.517516f, 1.548387f, 1.302745f, -1.908279f, 1.806452f, 0.598726f, 
	-1.989738f, 2.064516f, -0.202337f, -1.748693f, 2.322581f, -0.970605f, -1.224211f, 2.580645f, -1.581552f, -0.501304f, 2.838710f, -1.936155f, 0.302857f, 3.096774f, -1.976936f, 1.057930f, 3.354839f, -1.697288f, 1.641528f, 3.612903f, -1.142535f, 1.959060f, 3.870968f, -0.402595f, 
	1.959059f, 4.129032f, 0.402599f, 1.641525f, 4.387096f, 1.142539f, 1.057926f, 4.645161f, 1.697290f, 0.302853f, 4.903226f, 1.976937f, -0.501307f, 5.161290f, 1.936154f, -1.224214f, 5.419354f, 1.581550f, -1.748694f, 5.677419f, 0.970602f, -1.989739f, 5.935484f, 0.202334f, 
	-1.908278f, 6.193548f, -0.598729f, -1.517514f, 6.451612f, -1.302747f, -0.880785f, 6.709677f, -1.795611f, -0.101295f, 6.967742f, -1.997433f, 0.694614f, 7.225806f, -1.875503f, 1.377937f, 7.483871f, -1.449583f, 1.837917f, 7.741935f, -0.788708f, 2.000000f, 8.000000f, 0.000005f
};

const float tfx__path_preset_spiral_ball[192] = {
	0.000000f, 0.000000f, 0.000000f, 0.164738f, 0.126984f, 0.112317f, 0.145505f, 0.253968f, 0.370740f, -0.132660f, 0.380952f, 0.581222f, -0.581006f, 0.507937f, 0.539094f, -0.976005f, 0.634921f, 0.147109f, -1.062261f, 0.761905f, -0.511558f, -0.684040f, 0.888889f, -1.184793f, 
	0.116111f, 1.015873f, -1.549395f, 1.082089f, 1.142857f, -1.356896f, 1.828026f, 1.269841f, -0.563871f, 1.993077f, 1.396825f, 0.614784f, 1.404896f, 1.523810f, 1.761686f, 0.180488f, 1.650794f, 2.408464f, -1.285576f, 1.777778f, 2.226681f, -2.451258f, 1.904762f, 1.180464f, 
	-2.831485f, 2.031746f, -0.426776f, -2.198516f, 2.158730f, -2.039921f, -0.695899f, 2.285714f, -3.048916f, 1.186532f, 2.412699f, -3.023249f, 2.776244f, 2.539683f, -1.892820f, 3.464102f, 2.666667f, -0.000009f, 2.940992f, 2.793651f, 2.005122f, 1.332053f, 2.920635f, 3.393983f, 
	-0.828542f, 3.047619f, 3.630142f, -2.779509f, 3.174603f, 2.579029f, -3.807487f, 3.301588f, 0.573904f, -3.513527f, 3.428572f, -1.692004f, -1.969634f, 3.555556f, -3.411464f, 0.296578f, 3.682540f, -3.957861f, 2.486967f, 3.809524f, -3.118597f, 3.821096f, 3.936508f, -1.178679f, 
	3.821111f, 4.063492f, 1.178628f, 2.487008f, 4.190476f, 3.118563f, 0.296631f, 4.317461f, 3.957857f, -1.969588f, 4.444445f, 3.411490f, -3.513505f, 4.571429f, 1.692051f, -3.807495f, 4.698413f, -0.573853f, -2.779543f, 4.825397f, -2.578992f, -0.828590f, 4.952381f, -3.630131f, 
	1.332008f, 5.079365f, -3.394000f, 2.940965f, 5.206349f, -2.005161f, 3.464102f, 5.333333f, -0.000037f, 2.776269f, 5.460318f, 1.892783f, 1.186572f, 5.587302f, 3.023233f, -0.695859f, 5.714286f, 3.048925f, -2.198489f, 5.841270f, 2.039950f, -2.831479f, 5.968254f, 0.426814f, 
	-2.451273f, 6.095239f, -1.180431f, -1.285605f, 6.222223f, -2.226664f, 0.180456f, 6.349207f, -2.408466f, 1.404872f, 6.476191f, -1.761705f, 1.993068f, 6.603175f, -0.614812f, 1.828034f, 6.730159f, 0.563844f, 1.082106f, 6.857143f, 1.356881f, 0.116130f, 6.984128f, 1.549393f, 
	-0.684028f, 7.111112f, 1.184799f, -1.062256f, 7.238096f, 0.511567f, -0.976006f, 7.365080f, -0.147102f, -0.581009f, 7.492064f, -0.539091f, -0.132663f, 7.619048f, -0.581221f, 0.145504f, 7.746032f, -0.370741f, 0.164738f, 7.873016f, -0.112317f, -0.000000f, 8.000000f, 0.000000f
};

const float tfx__path_preset_triangle[153] = {
	2.000023f, -0.031998f, 2.149641f, 2.000017f, -0.000016f, 1.749991f, 2.000014f, -0.000011f, 1.500001f, 2.000014f, -0.000012f, 1.249997f, 2.000011f, -0.000009f, 0.999994f, 2.000007f, -0.000006f, 0.750001f, 2.000005f, -0.000003f, 0.499995f, 2.000004f, 0.000001f, 0.250000f, 
	2.000000f, 0.000000f, 0.000000f, 1.999997f, -0.000002f, -0.249999f, 1.999997f, 0.000003f, -0.499998f, 1.999990f, 0.000006f, -0.750004f, 1.999987f, 0.000009f, -0.999995f, 1.999988f, 0.000013f, -1.249999f, 1.999984f, 0.000011f, -1.500003f, 1.999980f, 0.000017f, -1.749989f, 
	1.999976f, 0.000013f, -1.999991f, 1.999982f, 0.216520f, -1.875003f, 1.999985f, 0.433024f, -1.749983f, 1.999988f, 0.649529f, -1.624991f, 1.999986f, 0.866037f, -1.499987f, 1.999988f, 1.082540f, -1.374986f, 1.999995f, 1.299042f, -1.249981f, 1.999990f, 1.515547f, -1.124981f, 
	1.999994f, 1.732054f, -0.999978f, 1.999997f, 1.948562f, -0.874979f, 1.999997f, 2.165062f, -0.749977f, 2.000000f, 2.381576f, -0.624980f, 2.000002f, 2.598070f, -0.499979f, 1.999998f, 2.814576f, -0.374981f, 2.000005f, 3.031076f, -0.249974f, 2.000004f, 3.247581f, -0.124973f, 
	2.000012f, 3.464082f, 0.000031f, 2.000006f, 3.247582f, 0.125023f, 2.000008f, 3.031068f, 0.250023f, 2.000009f, 2.814576f, 0.375022f, 2.000010f, 2.598064f, 0.500018f, 2.000013f, 2.381555f, 0.625019f, 2.000011f, 2.165051f, 0.750016f, 2.000012f, 1.948545f, 0.875012f, 
	2.000012f, 1.732035f, 1.000010f, 2.000012f, 1.515533f, 1.125007f, 2.000014f, 1.299023f, 1.250004f, 2.000015f, 1.082518f, 1.375001f, 2.000019f, 0.866011f, 1.500007f, 2.000020f, 0.649503f, 1.624998f, 2.000020f, 0.433000f, 1.749994f, 2.000018f, 0.224505f, 1.861108f, 
	2.000019f, 0.024916f, 1.934911f, 2.000017f, -0.000016f, 1.749991f, 2.000014f, -0.000011f, 1.500001f
};

const float tfx__path_preset_whirpool[144] = {
	0.400000f, 0.000000f, 0.000000f, 0.384224f, 0.170213f, 0.078118f, 0.339605f, 0.340426f, 0.144351f, 0.272871f, 0.510638f, 0.189204f, 0.193099f, 0.680851f, 0.206407f, 0.110451f, 0.851064f, 0.192974f, 0.035343f, 1.021277f, 0.148923f, -0.021956f, 1.191489f, 0.077072f, 
	-0.051617f, 1.361702f, -0.016844f, -0.044814f, 1.531915f, -0.123710f, 0.005076f, 1.702128f, -0.230557f, 0.100556f, 1.872340f, -0.320578f, 0.237647f, 2.042553f, -0.374218f, 0.403721f, 2.212766f, -0.371743f, 0.576430f, 2.382979f, -0.297264f, 0.724632f, 2.553191f, -0.143712f, 
	0.811930f, 2.723404f, 0.082297f, 0.802886f, 2.893617f, 0.357090f, 0.671199f, 3.063830f, 0.639517f, 0.408413f, 3.234042f, 0.874536f, 0.031124f, 3.404255f, 1.001168f, -0.415472f, 3.574468f, 0.964120f, -0.859555f, 3.744681f, 0.727343f, -1.211131f, 3.914893f, 0.286825f, 
	-1.376618f, 4.085106f, -0.320462f, -1.277732f, 4.255319f, -1.014406f, -0.871450f, 4.425532f, -1.678399f, -0.166918f, 4.595745f, -2.175241f, 0.764995f, 4.765957f, -2.370820f, 1.793030f, 4.936170f, -2.162077f, 2.739625f, 5.106383f, -1.504222f, 3.407241f, 5.276596f, -0.431332f, 
	3.613751f, 5.446808f, 0.935123f, 3.231270f, 5.617021f, 2.393770f, 2.221145f, 5.787234f, 3.688119f, 0.657393f, 5.957447f, 4.546606f, -1.268013f, 6.127659f, 4.732430f, -3.262012f, 6.297872f, 4.094887f, -4.969662f, 6.468085f, 2.612348f, -6.031480f, 6.638298f, 0.417148f, 
	-6.150524f, 6.808510f, -2.205142f, -5.157697f, 6.978723f, -4.845751f, -3.062493f, 7.148936f, -7.030491f, -0.077440f, 7.319149f, -8.298226f, 3.392055f, 7.489361f, -8.287976f, 6.795220f, 7.659574f, -6.819407f, 9.516871f, 7.829787f, -3.950694f, 10.981000f, 8.000000f, 0.000001f
};

const float tfx__path_preset_loop[192] = {
	0.000000f, 0.000000f, 0.000000f, 0.250000f, 0.000000f, 0.000000f, 0.500000f, 0.000000f, 0.000000f, 0.750000f, 0.000000f, 0.000000f, 1.000000f, 0.000000f, 0.000000f, 1.250000f, 0.000000f, 0.000000f, 1.500000f, 0.000000f, 0.000000f, 1.750000f, 0.000000f, 0.000000f, 
	2.000000f, 0.000000f, 0.000000f, 2.250000f, 0.000000f, 0.000000f, 2.500000f, 0.000000f, 0.000000f, 2.750000f, 0.000000f, 0.000000f, 3.000000f, 0.000000f, 0.000000f, 3.250000f, 0.000000f, 0.000000f, 3.500000f, 0.000000f, 0.000000f, 3.750000f, 0.000000f, 0.000000f, 
	4.000000f, 0.000000f, 0.000000f, 4.248292f, 0.025249f, 0.000000f, 4.486419f, 0.099962f, 0.000000f, 4.704631f, 0.221079f, 0.000000f, 4.893996f, 0.383644f, 0.000000f, 5.046761f, 0.581000f, 0.000000f, 5.156671f, 0.805067f, 0.000000f, 5.219227f, 1.046672f, 0.000000f, 
	5.231868f, 1.295924f, 0.000000f, 5.194076f, 1.542618f, 0.000000f, 5.107398f, 1.776655f, 0.000000f, 4.975383f, 1.988454f, 0.000000f, 4.803436f, 2.169342f, 0.000000f, 4.598596f, 2.311914f, 0.000000f, 4.369250f, 2.410335f, 0.000000f, 4.124786f, 2.460573f, 0.000000f, 
	3.875214f, 2.460573f, 0.000000f, 3.630750f, 2.410335f, 0.000000f, 3.401404f, 2.311914f, 0.000000f, 3.196564f, 2.169341f, 0.000000f, 3.024617f, 1.988453f, 0.000000f, 2.892602f, 1.776655f, 0.000000f, 2.805924f, 1.542618f, 0.000000f, 2.768132f, 1.295924f, 0.000000f, 
	2.780773f, 1.046672f, 0.000000f, 2.843329f, 0.805066f, 0.000000f, 2.953239f, 0.581000f, 0.000000f, 3.106004f, 0.383644f, 0.000000f, 3.295369f, 0.221079f, 0.000000f, 3.513582f, 0.099962f, 0.000000f, 3.751708f, 0.025249f, 0.000000f, 4.000000f, 0.000000f, 0.000000f, 
	4.248292f, -0.003906f, 0.000000f, 4.497859f, -0.003063f, 0.000000f, 4.747427f, -0.002220f, 0.000000f, 4.997021f, -0.001377f, 0.000000f, 5.246576f, -0.000533f, 0.000000f, 5.496140f, 0.000310f, 0.000000f, 5.745711f, 0.001153f, 0.000000f, 5.995277f, 0.001996f, 0.000000f, 
	6.244858f, 0.002839f, 0.000000f, 6.494419f, 0.003682f, 0.000000f, 6.744004f, 0.004525f, 0.000000f, 6.993561f, 0.005368f, 0.000000f, 7.243145f, 0.006211f, 0.000000f, 7.492707f, 0.007054f, 0.000000f, 7.742268f, 0.007898f, 0.000000f, 7.991854f, 0.008740f, 0.000000f
};

const float tfx__path_preset_whispy_loop[288] = {
	-17.304842f, 1.473904f, 0.000000f, -16.850300f, 1.496524f, 0.000000f, -16.396023f, 1.513536f, 0.000000f, -15.942352f, 1.520732f, 0.000000f, -15.489708f, 1.515554f, 0.000000f, -15.038718f, 1.496492f, 0.000000f, -14.590525f, 1.462688f, 0.000000f, -14.146774f, 1.414106f, 0.000000f, 
	-13.708881f, 1.351725f, 0.000000f, -13.277054f, 1.277372f, 0.000000f, -12.849972f, 1.193340f, 0.000000f, -12.425500f, 1.101992f, 0.000000f, -12.001793f, 1.005327f, 0.000000f, -11.577856f, 0.904556f, 0.000000f, -11.153130f, 0.800006f, 0.000000f, -10.726491f, 0.691534f, 0.000000f, 
	-10.295582f, 0.579239f, 0.000000f, -9.857334f, 0.464145f, 0.000000f, -9.409528f, 0.348644f, 0.000000f, -8.952070f, 0.236598f, 0.000000f, -8.486689f, 0.132937f, 0.000000f, -8.015314f, 0.042629f, 0.000000f, -7.539016f, -0.030528f, 0.000000f, -7.058889f, -0.084389f, 0.000000f, 
	-6.578042f, -0.117661f, 0.000000f, -6.102313f, -0.128383f, 0.000000f, -5.638580f, -0.113012f, 0.000000f, -5.192001f, -0.066837f, 0.000000f, -4.764553f, 0.014874f, 0.000000f, -4.355896f, 0.136099f, 0.000000f, -3.965680f, 0.299964f, 0.000000f, -3.595729f, 0.508557f, 0.000000f, 
	-3.251173f, 0.762077f, 0.000000f, -2.940368f, 1.057937f, 0.000000f, -2.673671f, 1.390961f, 0.000000f, -2.461322f, 1.754601f, 0.000000f, -2.311167f, 2.141892f, 0.000000f, -2.227362f, 2.545393f, 0.000000f, -2.210760f, 2.956773f, 0.000000f, -2.260533f, 3.366733f, 0.000000f, 
	-2.375752f, 3.764924f, 0.000000f, -2.555672f, 4.139495f, 0.000000f, -2.798370f, 4.477067f, 0.000000f, -3.098571f, 4.763999f, 0.000000f, -3.446198f, 4.988232f, 0.000000f, -3.826763f, 5.140059f, 0.000000f, -4.223272f, 5.211694f, 0.000000f, -4.618288f, 5.197400f, 0.000000f, 
	-4.995145f, 5.095395f, 0.000000f, -5.338679f, 4.910221f, 0.000000f, -5.636212f, 4.653164f, 0.000000f, -5.878634f, 4.339877f, 0.000000f, -6.060548f, 3.986773f, 0.000000f, -6.179080f, 3.608354f, 0.000000f, -6.232460f, 3.216405f, 0.000000f, -6.219676f, 2.820448f, 0.000000f, 
	-6.141235f, 2.428468f, 0.000000f, -5.999877f, 2.047423f, 0.000000f, -5.800480f, 1.683519f, 0.000000f, -5.549352f, 1.342267f, 0.000000f, -5.253495f, 1.028161f, 0.000000f, -4.920047f, 0.744088f, 0.000000f, -4.555899f, 0.490928f, 0.000000f, -4.167542f, 0.267760f, 0.000000f, 
	-3.761050f, 0.072621f, 0.000000f, -3.341917f, -0.096637f, 0.000000f, -2.914650f, -0.241756f, 0.000000f, -2.482347f, -0.363877f, 0.000000f, -2.046430f, -0.463680f, 0.000000f, -1.606515f, -0.541768f, 0.000000f, -1.160533f, -0.599191f, 0.000000f, -0.705566f, -0.637927f, 0.000000f, 
	-0.239610f, -0.660939f, 0.000000f, 0.236498f, -0.671498f, 0.000000f, 0.718370f, -0.672033f, 0.000000f, 1.199346f, -0.663326f, 0.000000f, 1.673147f, -0.644822f, 0.000000f, 2.136246f, -0.615810f, 0.000000f, 2.588543f, -0.576493f, 0.000000f, 3.032221f, -0.528202f, 0.000000f, 
	3.469941f, -0.472900f, 0.000000f, 3.903761f, -0.412569f, 0.000000f, 4.335147f, -0.348808f, 0.000000f, 4.765448f, -0.282610f, 0.000000f, 5.196054f, -0.214412f, 0.000000f, 5.628194f, -0.144536f, 0.000000f, 6.062728f, -0.073901f, 0.000000f, 6.500143f, -0.004499f, 0.000000f, 
	6.940607f, 0.060715f, 0.000000f, 7.384008f, 0.118218f, 0.000000f, 7.830073f, 0.164109f, 0.000000f, 8.278580f, 0.194239f, 0.000000f, 8.729328f, 0.204941f, 0.000000f, 9.181557f, 0.194458f, 0.000000f, 9.633205f, 0.164099f, 0.000000f, 10.080896f, 0.118246f, 0.000000f
};

const float tfx__path_preset_whispy_loops[294] = {
	-7.529158f, 1.681942f, 0.000000f, -7.142783f, 1.568462f, 0.000000f, -6.752921f, 1.464058f, 0.000000f, -6.357323f, 1.372115f, 0.000000f, -5.955426f, 1.289967f, 0.000000f, -5.548116f, 1.213203f, 0.000000f, -5.137207f, 1.139685f, 0.000000f, -4.724632f, 1.069985f, 0.000000f, 
	-4.311468f, 1.005310f, 0.000000f, -3.897368f, 0.945646f, 0.000000f, -3.480921f, 0.889453f, 0.000000f, -3.060801f, 0.834207f, 0.000000f, -2.637120f, 0.776676f, 0.000000f, -2.212366f, 0.712569f, 0.000000f, -1.791649f, 0.635771f, 0.000000f, -1.382326f, 0.537690f, 0.000000f, 
	-0.993393f, 0.407589f, 0.000000f, -0.635212f, 0.234615f, 0.000000f, -0.319977f, 0.011205f, 0.000000f, -0.062630f, -0.263470f, 0.000000f, 0.118950f, -0.581136f, 0.000000f, 0.206083f, -0.923580f, 0.000000f, 0.184580f, -1.263960f, 0.000000f, 0.051539f, -1.570868f, 0.000000f, 
	-0.180424f, -1.814339f, 0.000000f, -0.484455f, -1.970956f, 0.000000f, -0.825052f, -2.025717f, 0.000000f, -1.164404f, -1.971745f, 0.000000f, -1.467618f, -1.810746f, 0.000000f, -1.707023f, -1.554492f, 0.000000f, -1.865399f, -1.224327f, 0.000000f, -1.937004f, -0.847059f, 0.000000f, 
	-1.925545f, -0.449528f, 0.000000f, -1.839723f, -0.054994f, 0.000000f, -1.688502f, 0.317772f, 0.000000f, -1.478510f, 0.653783f, 0.000000f, -1.214314f, 0.941115f, 0.000000f, -0.900278f, 1.171363f, 0.000000f, -0.542412f, 1.341232f, 0.000000f, -0.149442f, 1.453624f, 0.000000f, 
	0.267423f, 1.516608f, 0.000000f, 0.696457f, 1.540513f, 0.000000f, 1.128368f, 1.535012f, 0.000000f, 1.558161f, 1.507956f, 0.000000f, 1.984687f, 1.465991f, 0.000000f, 2.409075f, 1.415664f, 0.000000f, 2.833758f, 1.363911f, 0.000000f, 3.261790f, 1.317908f, 0.000000f, 
	3.695247f, 1.284855f, 0.000000f, 4.133231f, 1.271930f, 0.000000f, 4.571299f, 1.286333f, 0.000000f, 5.002564f, 1.335258f, 0.000000f, 5.418790f, 1.425668f, 0.000000f, 5.810287f, 1.563913f, 0.000000f, 6.165227f, 1.755153f, 0.000000f, 6.469427f, 2.001976f, 0.000000f, 
	6.707106f, 2.301856f, 0.000000f, 6.862802f, 2.644515f, 0.000000f, 6.924239f, 3.011183f, 0.000000f, 6.885045f, 3.376563f, 0.000000f, 6.746232f, 3.712410f, 0.000000f, 6.516541f, 3.991571f, 0.000000f, 6.212396f, 4.192028f, 0.000000f, 5.857466f, 4.299921f, 0.000000f, 
	5.481092f, 4.309883f, 0.000000f, 5.115286f, 4.222881f, 0.000000f, 4.790731f, 4.044353f, 0.000000f, 4.532431f, 3.784601f, 0.000000f, 4.356110f, 3.459824f, 0.000000f, 4.267026f, 3.090841f, 0.000000f, 4.261900f, 2.699294f, 0.000000f, 4.332578f, 2.304014f, 0.000000f, 
	4.469344f, 1.919820f, 0.000000f, 4.663035f, 1.558138f, 0.000000f, 4.906178f, 1.227414f, 0.000000f, 5.193066f, 0.932661f, 0.000000f, 5.518806f, 0.675352f, 0.000000f, 5.878119f, 0.454511f, 0.000000f, 6.264873f, 0.268129f, 0.000000f, 6.672329f, 0.113816f, 0.000000f, 
	7.093729f, -0.011237f, 0.000000f, 7.523129f, -0.110162f, 0.000000f, 7.956268f, -0.185786f, 0.000000f, 8.390854f, -0.239843f, 0.000000f, 8.825981f, -0.272614f, 0.000000f, 9.261103f, -0.283194f, 0.000000f, 9.695192f, -0.270282f, 0.000000f, 10.126706f, -0.233271f, 0.000000f, 
	10.554363f, -0.173190f, 0.000000f, 10.977834f, -0.092745f, 0.000000f, 11.397480f, 0.004512f, 0.000000f, 11.813650f, 0.114647f, 0.000000f, 12.226837f, 0.232563f, 0.000000f, 12.638742f, 0.350849f, 0.000000f, 13.052684f, 0.460319f, 0.000000f, 13.472228f, 0.553206f, 0.000000f, 
	13.898808f, 0.627570f, 0.000000f, 14.330275f, 0.689663f, 0.000000f
};

const float tfx__path_preset_infinity_loop1[204] = {
	0.000000f, 0.972837f, 0.199923f, 0.000000f, 0.474194f, 0.222251f, 0.000000f, 0.972837f, 0.199923f, 0.000000f, 1.465767f, 0.278430f, 0.000000f, 1.932805f, 0.454554f, 0.000000f, 2.354820f, 0.721089f, 0.000000f, 2.714549f, 1.067117f, 0.000000f, 2.997271f, 1.478482f, 
	0.000000f, 3.191392f, 1.938331f, 0.000000f, 3.288966f, 2.427840f, 0.000000f, 3.286007f, 2.926978f, 0.000000f, 3.182640f, 3.415303f, 0.000000f, 2.983083f, 3.872824f, 0.000000f, 2.695518f, 4.280817f, 0.000000f, 2.331704f, 4.622559f, 0.000000f, 1.906557f, 4.884069f, 
	0.000000f, 1.437469f, 5.054635f, 0.000000f, 0.943642f, 5.127299f, 0.000000f, 0.445299f, 5.099059f, 0.000000f, -0.037162f, 4.971093f, 0.000000f, -0.483981f, 4.748613f, 0.000000f, -0.876878f, 4.440752f, 0.000000f, -1.199755f, 4.060121f, 0.000000f, -1.439392f, 3.622257f, 
	0.000000f, -1.585991f, 3.145127f, 0.000000f, -1.633552f, 2.648251f, 0.000000f, -1.580106f, 2.151973f, 0.000000f, -1.427860f, 1.676617f, 0.000000f, -1.183049f, 1.241638f, 0.000000f, -0.855691f, 0.864840f, 0.000000f, -0.459175f, 0.561658f, 0.000000f, -0.009749f, 0.344497f, 
	0.000000f, 0.474194f, 0.222251f, 0.000000f, 0.972837f, 0.199923f, 0.000000f, 1.466098f, 0.128941f, 0.000000f, 1.936032f, -0.036906f, 0.000000f, 2.364573f, -0.291246f, 0.000000f, 2.735264f, -0.624307f, 0.000000f, 3.033877f, -1.023285f, 0.000000f, 3.248895f, -1.472852f, 
	0.000000f, 3.372081f, -1.955724f, 0.000000f, 3.398690f, -2.453358f, 0.000000f, 3.327711f, -2.946620f, 0.000000f, 3.161860f, -3.416557f, 0.000000f, 2.907515f, -3.845119f, 0.000000f, 2.574457f, -4.215808f, 0.000000f, 2.175478f, -4.514409f, 0.000000f, 1.725917f, -4.729424f, 
	0.000000f, 1.243040f, -4.852621f, 0.000000f, 0.745410f, -4.879203f, 0.000000f, 0.252149f, -4.808245f, 0.000000f, -0.217784f, -4.642383f, 0.000000f, -0.646332f, -4.388041f, 0.000000f, -1.017027f, -4.054993f, 0.000000f, -1.315619f, -3.656011f, 0.000000f, -1.530642f, -3.206446f, 
	0.000000f, -1.653824f, -2.723568f, 0.000000f, -1.680438f, -2.225932f, 0.000000f, -1.609451f, -1.732673f, 0.000000f, -1.443605f, -1.262736f, 0.000000f, -1.189268f, -0.834186f, 0.000000f, -0.856214f, -0.463492f, 0.000000f, -0.457235f, -0.164897f, 0.000000f, -0.007671f, 0.050126f, 
	0.000000f, 0.475207f, 0.173310f, 0.000000f, 0.972837f, 0.199923f, 0.000000f, 0.474194f, 0.222251f, 0.000000f, 0.972837f, 0.199923f
};

const float tfx__path_preset_infinity_loop2[198] = {
	0.000000f, 0.144379f, -0.368948f, 0.000000f, 0.546153f, -0.013604f, 0.000000f, 0.964475f, 0.304307f, 0.000000f, 1.355363f, 0.692264f, 0.000000f, 1.711848f, 1.152925f, 0.000000f, 2.026991f, 1.680841f, 0.000000f, 2.294327f, 2.260695f, 0.000000f, 2.508458f, 2.867810f, 
	0.000000f, 2.658038f, 3.482686f, 0.000000f, 2.731355f, 4.100725f, 0.000000f, 2.725789f, 4.706430f, 0.000000f, 2.641743f, 5.272020f, 0.000000f, 2.482731f, 5.777150f, 0.000000f, 2.255292f, 6.206300f, 0.000000f, 1.968747f, 6.548832f, 0.000000f, 1.634830f, 6.798283f, 
	0.000000f, 1.267211f, 6.951271f, 0.000000f, 0.880942f, 7.006386f, 0.000000f, 0.491837f, 6.963348f, 0.000000f, 0.115825f, 6.822663f, 0.000000f, -0.231699f, 6.585835f, 0.000000f, -0.536508f, 6.256096f, 0.000000f, -0.786122f, 5.839472f, 0.000000f, -0.970319f, 5.345869f, 
	0.000000f, -1.081555f, 4.789821f, 0.000000f, -1.115262f, 4.190497f, 0.000000f, -1.070034f, 3.570651f, 0.000000f, -0.947668f, 2.954389f, 0.000000f, -0.753076f, 2.363954f, 0.000000f, -0.494056f, 1.816068f, 0.000000f, -0.180945f, 1.318688f, 0.000000f, 0.173828f, 0.869130f, 
	0.000000f, 0.556266f, 0.454306f, 0.000000f, 0.951382f, 0.053298f, 0.000000f, 1.343786f, -0.358178f, 0.000000f, 1.718285f, -0.802469f, 0.000000f, 2.060445f, -1.295167f, 0.000000f, 2.357123f, -1.841933f, 0.000000f, 2.596952f, -2.437589f, 0.000000f, 2.770757f, -3.067433f, 
	0.000000f, 2.871894f, -3.710153f, 0.000000f, 2.896503f, -4.341445f, 0.000000f, 2.843656f, -4.937447f, 0.000000f, 2.715392f, -5.477365f, 0.000000f, 2.516648f, -5.944986f, 0.000000f, 2.255062f, -6.329090f, 0.000000f, 1.940688f, -6.622950f, 0.000000f, 1.585608f, -6.823258f, 
	0.000000f, 1.203467f, -6.928858f, 0.000000f, 0.808951f, -6.939608f, 0.000000f, 0.417221f, -6.855692f, 0.000000f, 0.043330f, -6.677515f, 0.000000f, -0.298352f, -6.406212f, 0.000000f, -0.594696f, -6.044653f, 0.000000f, -0.834313f, -5.598682f, 0.000000f, -1.007993f, -5.078264f, 
	0.000000f, -1.109053f, -4.498192f, 0.000000f, -1.135694f, -3.864655f, 0.000000f, -1.088147f, -3.201942f, 0.000000f, -0.967038f, -2.543253f, 0.000000f, -0.775831f, -1.911750f, 0.000000f, -0.520616f, -1.331681f, 0.000000f, -0.209799f, -0.818895f, 0.000000f, 0.144379f, -0.368948f, 
	0.000000f, 0.546153f, -0.013604f, 0.000000f, 0.964475f, 0.304307f
};

const float tfx__path_preset_fly_path[276] = {
	-13.453751f, -0.114839f, 0.810795f, -11.865257f, 0.000000f, 0.281245f, -10.098850f, 0.000000f, 0.000000f, -9.063751f, 0.000000f, 0.000000f, -8.028650f, 0.000000f, 0.000000f, -6.993552f, 0.000000f, 0.000000f, -5.958451f, 0.000000f, 0.000000f, -5.059327f, 0.124914f, 0.497396f, 
	-4.160202f, 0.249828f, 0.994793f, -3.261076f, 0.374742f, 1.492190f, -2.361951f, 0.499656f, 1.989586f, -2.155735f, 0.367824f, 2.995332f, -1.949517f, 0.235992f, 4.001080f, -1.743299f, 0.104160f, 5.006827f, -1.537082f, -0.027672f, 6.012574f, -1.330865f, -0.159505f, 7.018320f, 
	-1.124648f, -0.291337f, 8.024067f, -0.918430f, -0.423169f, 9.029815f, 0.065629f, -0.260775f, 9.306737f, 1.049687f, -0.098382f, 9.583660f, 2.033750f, 0.064013f, 9.860583f, 3.017809f, 0.226407f, 10.137506f, 4.001871f, 0.388801f, 10.414429f, 4.863249f, 0.141254f, 9.896571f, 
	5.724630f, -0.106293f, 9.378711f, 6.586012f, -0.353840f, 8.860853f, 7.447392f, -0.601387f, 8.342998f, 8.308767f, -0.848934f, 7.825137f, 9.170145f, -1.096480f, 7.307280f, 9.054326f, -1.201880f, 6.284093f, 8.938505f, -1.307280f, 5.260907f, 8.822686f, -1.412680f, 4.237725f, 
	8.706867f, -1.518079f, 3.214543f, 8.591048f, -1.623479f, 2.191361f, 8.475229f, -1.728878f, 1.168179f, 8.359410f, -1.834278f, 0.144998f, 8.243591f, -1.939677f, -0.878184f, 7.210654f, -1.944792f, -0.811549f, 6.177718f, -1.949906f, -0.744915f, 5.144782f, -1.955021f, -0.678280f, 
	4.111845f, -1.960136f, -0.611645f, 3.427018f, -1.633621f, 0.092504f, 2.742191f, -1.307106f, 0.796652f, 2.057363f, -0.980592f, 1.500800f, 1.372537f, -0.654077f, 2.204949f, 0.687710f, -0.327562f, 2.909097f, 0.002882f, -0.001048f, 3.613246f, -0.681942f, 0.325466f, 4.317392f, 
	-1.366772f, 0.651982f, 5.021543f, -2.051598f, 0.978495f, 5.725689f, -2.736426f, 1.305011f, 6.429838f, -3.726558f, 1.242576f, 6.725073f, -4.716691f, 1.180140f, 7.020308f, -5.706825f, 1.117707f, 7.315540f, -5.880612f, 1.125261f, 6.295168f, -6.054396f, 1.132815f, 5.274794f, 
	-6.228182f, 1.140369f, 4.254415f, -6.401966f, 1.147923f, 3.234039f, -6.575748f, 1.155478f, 2.213663f, -6.749534f, 1.163033f, 1.193282f, -7.327647f, 0.776775f, 0.426466f, -7.905760f, 0.390515f, -0.340355f, -8.483875f, 0.004254f, -1.107171f, -9.061987f, -0.382006f, -1.873988f, 
	-8.926970f, -0.285420f, -2.895683f, -8.791956f, -0.188834f, -3.917377f, -8.656939f, -0.092248f, -4.939072f, -8.521921f, 0.004338f, -5.960767f, -8.386904f, 0.100924f, -6.982462f, -9.402067f, 0.069902f, -6.782714f, -10.417231f, 0.038879f, -6.582964f, -11.432394f, 0.007856f, -6.383214f, 
	-12.447557f, -0.023168f, -6.183463f, -13.462721f, -0.054190f, -5.983714f, -14.296976f, -0.082278f, -5.371631f, -15.131229f, -0.110365f, -4.759547f, -15.965483f, -0.138453f, -4.147460f, -16.799740f, -0.166540f, -3.535377f, -17.633997f, -0.194626f, -2.923287f, -18.468250f, -0.222714f, -2.311207f, 
	-18.467600f, -0.210803f, -1.276187f, -18.466948f, -0.198893f, -0.241166f, -18.466297f, -0.186982f, 0.793854f, -18.465647f, -0.175071f, 1.828873f, -18.464996f, -0.163160f, 2.863894f, -17.505991f, -0.167634f, 2.474413f, -16.546984f, -0.172109f, 2.084929f, -15.587976f, -0.176584f, 1.695447f, 
	-14.628969f, -0.181058f, 1.305965f, -13.453751f, -0.114839f, 0.810795f, -11.865257f, 0.000000f, 0.281245f, -10.098850f, 0.000000f, 0.000000f
};

const float tfx__path_preset_meander[99] = {
	-1.130700f, 0.885696f, -0.573026f, -0.671597f, 1.094656f, -0.472017f, -0.286400f, 1.344134f, -0.269407f, 0.032087f, 1.618764f, 0.006753f, 0.319323f, 1.904248f, 0.315902f, 0.636071f, 2.184167f, 0.600827f, 1.006661f, 2.446089f, 0.831720f, 1.432142f, 2.665816f, 1.000588f, 
	1.901490f, 2.827921f, 1.105987f, 2.396248f, 2.934872f, 1.156823f, 2.902027f, 2.996876f, 1.170794f, 3.411440f, 3.026704f, 1.169021f, 3.921411f, 3.038018f, 1.175556f, 4.428101f, 3.046695f, 1.221766f, 4.918196f, 3.070997f, 1.346873f, 5.363026f, 3.125611f, 1.580408f, 
	5.747606f, 3.207853f, 1.902898f, 6.088867f, 3.304912f, 2.269241f, 6.413656f, 3.406392f, 2.648317f, 6.748459f, 3.503600f, 3.019661f, 7.107994f, 3.580318f, 3.370543f, 7.414164f, 3.636668f, 3.768565f, 7.636506f, 3.696107f, 4.224109f, 7.827977f, 3.757652f, 4.692949f, 
	8.023890f, 3.829451f, 5.156253f, 8.224380f, 3.926459f, 5.613540f, 8.394526f, 4.062737f, 6.082202f, 8.514643f, 4.224124f, 6.542165f, 8.589557f, 4.406023f, 7.016800f, 8.625130f, 4.599983f, 7.492012f, 8.639706f, 4.805804f, 7.967689f, 8.698837f, 5.031759f, 8.419838f, 
	8.850214f, 5.274874f, 8.810685f
};
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
tfxMAKE_HANDLE(tfx_effect_manager);
tfxMAKE_HANDLE(tfx_animation_manager);
tfxMAKE_HANDLE(tfx_effect_descriptor);
tfxMAKE_HANDLE(tfx_effect_template);
tfxMAKE_HANDLE(tfx_ribbon_buffer_requirements);
tfxMAKE_HANDLE(tfx_ribbon_dispatch);

//-----------------------------------------------------------
//Section: String_Buffers
//-----------------------------------------------------------

int tfx_FormatString(char *buf, size_t buf_size, const char *fmt, va_list args);

#ifdef __cplusplus

//Very simple string builder for short stack based strings
template<size_t N>
struct tfx_str_t {
	char data[N];
	static constexpr tfxU32 capacity = N;
	tfxU32 current_size;

	inline tfx_str_t() : current_size(0) { memset(data, '\0', N); }
	inline tfx_str_t(const char *text) : current_size(0) {
		size_t text_len = strlen(text);
		TFX_ASSERT(text_len < capacity); //String is too big for the buffer size allowed
		memset(data, '\0', N);
		tfx__strcpy(data, N, text);
		current_size = (tfxU32)text_len + 1;
	}

	inline bool            empty() { return current_size == 0; }
	inline char &operator[](tfxU32 i) { TFX_ASSERT(i < current_size); return data[i]; }
	inline const char &operator[](tfxU32 i) const { TFX_ASSERT(i < current_size); return data[i]; }

	inline void         Clear() { current_size = 0; memset(data, '\0', N); }
	inline char *begin() { return data; }
	inline const char *begin() const { return data; }
	inline char *end() { return data + current_size; }
	inline const char *end() const { return data + current_size; }
	inline char			back() { TFX_ASSERT(current_size > 0); return data[current_size - 1]; }
	inline const char	back() const { TFX_ASSERT(current_size > 0); return data[current_size - 1]; }
	inline void         pop() { TFX_ASSERT(current_size > 0); current_size--; }
	inline void         push_back(const char v) { TFX_ASSERT(current_size < capacity); data[current_size] = char(v); current_size++; }

	inline bool operator==(const char *string) { return !strcmp(string, data); }
	inline bool operator==(const tfx_str_t string) { return !strcmp(data, string.c_str()); }
	inline bool operator!=(const char *string) { return strcmp(string, data); }
	inline bool operator!=(const tfx_str_t string) { return strcmp(data, string.c_str()); }
	inline const char *c_str() const { return data; }
	int Find(const char *needle) {
		tfx_str_t compare = needle;
		tfx_str_t lower = Lower();
		compare = compare.Lower();
		if (compare.Length() > Length()) return -1;
		tfxU32 pos = 0;
		while (compare.Length() + pos <= Length()) {
			if (strncmp(lower.data + pos, compare.data, compare.Length()) == 0) {
				return pos;
			}
			++pos;
		}
		return -1;
	}
	tfx_str_t Lower() {
		tfx_str_t convert = *this;
		for (auto &c : convert) {
			c = tolower(c);
		}
		return convert;
	}
	inline tfxU32 Length() const { return current_size ? current_size - 1 : 0; }
	void AddLine(const char *format, ...) {
		if (!current_size) {
			NullTerminate();
		}
		va_list args;
		va_start(args, format);
		Appendv(format, args);
		va_end(args);
		Append('\n');
	}
	void Set(const char *text) { size_t len = strlen(text); TFX_ASSERT(len < capacity - 1); tfx__strcpy(data, N, text); current_size = (tfxU32)len; NullTerminate(); };
	void Setf(const char *format, ...) {
		Clear();
		va_list args;
		va_start(args, format);

		va_list args_copy;
		va_copy(args_copy, args);

		int len = tfx_FormatString(NULL, 0, format, args);         // FIXME-OPT: could do a first pass write attempt, likely successful on first pass.
		if (len <= 0)
		{
			va_end(args_copy);
			return;
		}

		const int write_off = (current_size != 0) ? current_size : 1;
		const int needed_sz = write_off + len;
		TFX_ASSERT(write_off + (tfxU32)len < capacity);	//Trying to write outside of buffer space, string too small!
		tfx_FormatString(&data[write_off - 1], (size_t)len + 1, format, args_copy);
		va_end(args_copy);

		va_end(args);
		current_size = len + 1;
	}
	void Appendf(const char *format, ...) {
		va_list args;
		va_start(args, format);

		va_list args_copy;
		va_copy(args_copy, args);

		int len = tfx_FormatString(NULL, 0, format, args);         // FIXME-OPT: could do a first pass write attempt, likely successful on first pass.
		if (len <= 0)
		{
			va_end(args_copy);
			return;
		}

		const int write_off = (current_size != 0) ? current_size : 1;
		const int needed_sz = write_off + len;
		TFX_ASSERT(write_off + (tfxU32)len < capacity);	//Trying to write outside of buffer space, string too small!
		tfx_FormatString(&data[write_off - 1], (size_t)len + 1, format, args_copy);
		va_end(args_copy);

		va_end(args);
		current_size += len;
	}
	void Appendv(const char *format, va_list args) {
		va_list args_copy;
		va_copy(args_copy, args);

		int len = tfx_FormatString(NULL, 0, format, args);         // FIXME-OPT: could do a first pass write attempt, likely successful on first pass.
		if (len <= 0)
		{
			va_end(args_copy);
			return;
		}

		const int write_off = (current_size != 0) ? current_size : 1;
		const int needed_sz = write_off + len;
		TFX_ASSERT(write_off + (tfxU32)len < capacity);	//Trying to write outside of buffer space, string too small!
		tfx_FormatString(&data[write_off - 1], (size_t)len + 1, format, args_copy);

		va_end(args_copy);
		current_size += len;
	}
	inline void Append(char c) { if (current_size) { pop(); } push_back(c); NullTerminate(); }
	inline void Pop() { if (!Length()) return; if (back() == 0) pop(); pop(); NullTerminate(); }
	inline void Trim(char c = ' ') { if (!Length()) return; if (back() == 0) pop(); while (back() == c && current_size) { pop(); } NullTerminate(); }
	inline void TrimToZero() { if (current_size < capacity) { memset(data + current_size, '\0', capacity - current_size); } }
	inline void TrimFront(char c = ' ') { if (!Length()) return; tfxU32 pos = 0; while (data[pos] == c && pos < current_size) { pos++; } if (pos < current_size) { memcpy(data, data + pos, current_size - pos); } current_size -= pos; }
	inline void SanitizeLineFeeds() { if (current_size > 1) { while (current_size > 1 && (back() == '\n' || back() == '\r' || back() == '\0')) { pop(); if (current_size <= 1) { break; } } NullTerminate(); } }
	inline void NullTerminate() { push_back('\0'); }
	inline const bool IsInt() {
		if (!Length()) return false;
		for (auto c : *this) {
			if (c >= '0' && c <= '9');
			else {
				return false;
			}
		}
		return true;
	}
	inline const bool IsFloat() {
		if (!Length()) return false;
		int dot_count = 0;
		for (auto c : *this) {
			if (c >= '0' && c <= '9');
			else if (c == '.' && dot_count == 0) {
				dot_count++;
			} else {
				return false;
			}
		}
		return dot_count == 1;
	}
} TFX_PACKED_STRUCT;

#define tfx_str16_t tfx_str_t<16>
#define tfx_str32_t tfx_str_t<32>
#define tfx_str64_t tfx_str_t<64>
#define tfx_str128_t tfx_str_t<128>
#define tfx_str256_t tfx_str_t<256>
#define tfx_str512_t tfx_str_t<512>
#define tfx_str1024_t tfx_str_t<1024>

typedef struct tfx_property_pair_s {
	tfx_str256_t property_name;
	tfx_str64_t property_value;
} tfx_property_pair_t;

//-----------------------------------------------------------
//Containers_and_Memory
//-----------------------------------------------------------

#else

typedef struct tfx_vector_s {
	tfxU32 current_size;
	tfxU32 capacity;
	tfxU32 volatile locked;
	tfxU32 alignment;
	void *data;
} tfx_vector_t;

#endif

#ifdef __cplusplus

//Storage
//Credit to ocornut https://github.com/ocornut for tfxvec although it's quite a lot different now.
//std::vector replacement with some extra stuff and tweaks specific to TimelineFX
template<typename T>
struct tfx_vector_t {
	tfxU32 current_size;
	tfxU32 capacity;
	tfxU32 volatile locked;
	tfxU32 alignment;
#if defined(TFX_MEMORY_TRACKING)
	tfx_memory_allocation_type allocation_type;
#endif
	T *data;

	// Provide standard typedefs but we don't use them ourselves.
	typedef T value_type;
	typedef value_type *iterator;
	typedef const value_type *const_iterator;

	inline					tfx_vector_t(const tfx_vector_t<T> &src) { locked = false; current_size = capacity = alignment = 0; data = nullptr; resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); }
	inline					tfx_vector_t() : locked(0), current_size(0), capacity(0), alignment(0), data(nullptr) {}
	//inline					tfx_vector_t<T> &operator=(const tfx_vector_t<T> &src) { TFX_ASSERT(0); return *this; }	//Use copy instead. 
	inline					~tfx_vector_t() { TFX_ASSERT(data == nullptr); } //You must manually free containers! Call the_containter.free();

	inline void				init() { locked = false; current_size = capacity = alignment = 0; data = nullptr; }
	inline bool				empty() { return current_size == 0; }
	inline bool				full() { return current_size == capacity; }
	inline tfxU32			size() { return current_size; }
	inline const tfxU32		size() const { return current_size; }
	inline tfxU32			size_in_bytes() { return current_size * sizeof(T); }
	inline const tfxU32		size_in_bytes() const { return current_size * sizeof(T); }
	inline T &operator[](tfxU32 i) { TFX_ASSERT(i < current_size); return data[i]; }
	inline const T &operator[](tfxU32 i) const { TFX_ASSERT(i < current_size); return data[i]; }
	inline T &ts_at(tfxU32 i) { while (locked > 0)
		; 
		return data[i]; 
	}

	inline void				free() { if (data) { current_size = capacity = 0; tfxFREE(data); data = nullptr; } }
	inline void				clear() { if (data) { current_size = 0; } }
	inline T *begin() { return data; }
	inline const T *begin() const { return data; }
	inline T *end() { return data + current_size; }
	inline const T *end() const { return data + current_size; }
	inline T *rend() { return data; }
	inline const T *rend() const { return data; }
	inline T *rbegin() { return data + current_size; }
	inline const T *rbegin() const { return data + current_size; }
	inline T &front() { TFX_ASSERT(current_size > 0); return data[0]; }
	inline const T &front() const { TFX_ASSERT(current_size > 0); return data[0]; }
	inline T &back() { TFX_ASSERT(current_size > 0); return data[current_size - 1]; }
	inline const T &back() const { TFX_ASSERT(current_size > 0); return data[current_size - 1]; }
	inline T &parent() { TFX_ASSERT(current_size > 1); return data[current_size - 2]; }
	inline const T &parent() const { TFX_ASSERT(current_size > 1); return data[current_size - 2]; }
	inline void				copy(const tfx_vector_t<T> &src) { clear(); resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); }
	inline tfxU32			_grow_capacity(tfxU32 sz) const { tfxU32 new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > sz ? new_capacity : sz; }
	inline void				resize(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
	inline void				resize_bytes(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
	inline void				resize(tfxU32 new_size, const T &v) { if (new_size > capacity) reserve(_grow_capacity(new_size)); if (new_size > current_size) for (tfxU32 n = current_size; n < new_size; n++) memcpy(&data[n], &v, sizeof(v)); current_size = new_size; }
	inline void				shrink(tfxU32 new_size) { TFX_ASSERT(new_size <= current_size); current_size = new_size; }
	inline void				set_alignment(tfxU32 align_to) { TFX_ASSERT(0 == (align_to & (align_to - 1)) && "must align to a power of two"); alignment = align_to; }
	inline void				reserve(tfxU32 new_capacity) {
		if (new_capacity <= capacity)
			return;
		T *new_data;
		if (alignment > 0) {
			new_data = (T *)tfxALLOCATE_ALIGNED((size_t)new_capacity * sizeof(T), alignment);
		} else {
			new_data = (T *)tfxALLOCATE((size_t)new_capacity * sizeof(T));
		}
		TFX_ASSERT(new_data);    //Unable to allocate memory. todo: better handling
#if defined(TFX_MEMORY_TRACKING)
		tfx_SetBlockMemoryType(new_data, allocation_type);
#endif
		if (data) {
			memcpy(new_data, data, (size_t)current_size * sizeof(T));
			tfxFREE(data);
		}
		data = new_data;
		capacity = new_capacity;
	}

	inline T &grab() {
		if (current_size == capacity) reserve(_grow_capacity(current_size + 1));
		current_size++;
		return data[current_size - 1];
	}
	inline tfxU32        locked_push_back(const T &v) {
		//suspect, just use a mutex instead?
		while (tfx__compare_and_exchange((tfxLONG volatile *)&locked, 1, 0) > 1);
		if (current_size == capacity) {
			reserve(_grow_capacity(current_size + 1));
		}
		memcpy(&data[current_size], &v, sizeof(T));
		tfxU32 index = current_size++;
		tfx__exchange((tfxLONG volatile *)&locked, 0);
		return index;
	}
	inline T &push_back(const T &v) {
		if (current_size == capacity) {
			reserve(_grow_capacity(current_size + 1));
		}
		memcpy(&data[current_size], &v, sizeof(T));
		current_size++; 
		return data[current_size - 1];
	}
	inline T &push_back_copy(const T &v) {
		if (current_size == capacity) {
			reserve(_grow_capacity(current_size + 1));
		}
		memcpy(&data[current_size], &v, sizeof(v));
		current_size++; 
		return data[current_size - 1];
	}
	inline T &next() {
		return push_back(T());
	}
	inline void				zero() { TFX_ASSERT(capacity > 0); memset(data, 0, capacity * sizeof(T)); }
	inline void				pop() { TFX_ASSERT(current_size > 0); current_size--; }
	inline T &pop_back() { TFX_ASSERT(current_size > 0); current_size--; return data[current_size]; }
	inline void				push_front(const T &v) { if (current_size == 0) push_back(v); else insert(data, v); }
	inline T *erase(const T *it) { TFX_ASSERT(it >= data && it < data + current_size); const ptrdiff_t off = it - data; memmove(data + off, data + off + 1, ((size_t)current_size - (size_t)off - 1) * sizeof(T)); current_size--; return data + off; }
	inline T				pop_front() { TFX_ASSERT(current_size > 0); T front = data[0]; erase(data); return front; }
	inline T *erase(const T *it, const T *it_last) { TFX_ASSERT(it >= data && it < data + current_size && it_last > it && it_last <= data + current_size); const ptrdiff_t count = it_last - it; const ptrdiff_t off = it - data; memmove(data + off, data + off + count, ((size_t)current_size - (size_t)off - count) * sizeof(T)); current_size -= (tfxU32)count; return data + off; }
	inline T *erase_unsorted(const T *it) { TFX_ASSERT(it >= data && it < data + current_size);  const ptrdiff_t off = it - data; if (it < data + current_size - 1) memcpy(data + off, data + current_size - 1, sizeof(T)); current_size--; return data + off; }
	inline T *insert(const T *it, const T &v) { TFX_ASSERT(it >= data && it <= data + current_size); const ptrdiff_t off = it - data; if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); if (off < (ptrdiff_t)current_size) memmove(data + off + 1, data + off, ((size_t)current_size - (size_t)off) * sizeof(T)); memcpy(data + off, &v, sizeof(T)); current_size++; return data + off; }
	inline T *insert_after(const T *it, const T &v) { TFX_ASSERT(it >= data && it <= data + current_size); const ptrdiff_t off = (it + 1) - data; if (current_size == capacity) reserve(_grow_capacity(current_size + 1)); if (off < (ptrdiff_t)current_size) memmove(data + off + 1, data + off, ((size_t)current_size - (size_t)off) * sizeof(T)); memcpy(data + off, &v, sizeof(T)); current_size++; return data + off; }
	inline bool				contains(const T &v) const { const T *_data = data;  const T *data_end = data + current_size; while (_data < data_end) if (*_data++ == v) return true; return false; }
	inline T *find(const T &v) { T *_data = data;  const T *data_end = data + current_size; while (_data < data_end) if (*_data == v) break; else ++_data; return _data; }
	inline const T *find(const T &v) const { const T *_data = data;  const T *data_end = data + current_size; while (_data < data_end) if (*_data == v) break; else ++_data; return _data; }
	inline bool				find_erase(const T &v) { const T *it = find(v); if (it < data + current_size) { erase(it); return true; } return false; }
	inline bool				find_erase_unsorted(const T &v) { const T *it = find(v); if (it < data + current_size) { erase_unsorted(it); return true; } return false; }
	inline tfxU32			index_from_ptr(const T *it) const { TFX_ASSERT(it >= data && it < data + current_size); const ptrdiff_t off = it - data; return (tfxU32)off; }
};

#define tfxCastBuffer(type, buffer) static_cast<type*>(buffer->data)
#define tfxCastBufferRef(type, buffer) static_cast<type*>(buffer.data)

#else

#define tfxCastBuffer(type, buffer) (type*)buffer->data
#define tfxCastBufferRef(type, buffer) (type*)buffer.data

#endif

//Used in tfx_soa_buffer_t to store pointers to arrays inside a struct of arrays
typedef struct tfx_soa_data_s {
	void *ptr;				//A pointer to the array in the struct
	size_t offset;			//The offset to the memory location in the buffer where the array starts
	size_t unit_size;		//The size of each data type in the array
}tfx_soa_data_t;

//A buffer designed to contain structs of arrays. If the arrays need to grow then a new memory block is made and all copied over
//together. All arrays in the struct will be the same capacity but can all have different unit sizes/types.
//In order to use this you need to first prepare the buffer by calling tfx__add_struct_array for each struct member of the SoA you're setting up. 
//All members must be of the same struct.
//Then call tfx__finish_soa_buffer_setup to create the memory for the struct of arrays with an initial reserve amount.
typedef struct tfx_soa_buffer_s {
	size_t current_arena_size;				//The current size of the arena that contains all the arrays
	size_t struct_size;						//The size of the struct (each member unit size added up)
	tfxU32 current_size;					//current size of each array
	tfxU32 start_index;						//Start index if you're using the buffer as a ring buffer
	tfxU32 last_bump;						//the amount the the start index was bumped by the last time tfx__bump_soa_buffer was called
	tfxU32 capacity;						//capacity of each array
	tfxU32 block_size;						//Keep the capacity to the nearest block size
	tfxU32 buffer_amount;					//Specify an amount of extra space there should be in addition to the max capacity. This is so that in a ring buffer if you want to align to the nearest block_size amount you
											//can avoid looping over data twice once at the beginning and again at the end
	tfxU32 alignment;						//The alignment of the memory. If you're planning on using simd for the memory, then 16 will be necessary.
#ifdef __cplusplus
	tfx_vector_t<tfx_soa_data_t> array_ptrs;    //Container for all the pointers into the arena
#else
	tfx_vector_t array_ptrs;    //Container for all the pointers into the arena
#endif
	void *user_data;
	void(*resize_callback)(struct tfx_soa_buffer_s *ring, tfxU32 new_index_start);
	void *struct_of_arrays;					//Pointer to the struct of arrays. Important that this is a stable pointer! Set with tfx__finish_soa_buffer_setup
	void *data;								//Pointer to the area in memory that contains all of the array data    
}tfx_soa_buffer_t;

#ifdef __cplusplus

//This simple container struct was created for storing instance_data in the effect manager. I didn't want this templated because either 2d or 3d instance_data could be used so
//I wanted to cast as needed when writing and using the sprite data. See simple cast macros above tfxCastBuffer and tfxCastBufferRef. However, now that 2d/3d effects are now unified
//can probably just delete this and use a tfx_vector_t container instead.
struct tfx_buffer_t {
	tfxU32 current_size;
	tfxU32 capacity;
	tfxU32 struct_size;
	tfxU32 alignment;
	void *data;

	inline				tfx_buffer_t() { struct_size = current_size = capacity = alignment = 0; data = nullptr; }
	inline void         free() { if (data) { current_size = capacity = alignment = 0; tfxFREE(data); data = nullptr; } }
	inline void         clear() { if (data) { current_size = 0; } }
	inline void			reserve(tfxU32 new_capacity) {
		TFX_ASSERT(struct_size);	//Must assign a value to struct size
		if (new_capacity <= capacity) return;
		void *new_data;
		if (alignment != 0) {
			new_data = tfxALLOCATE_ALIGNED((size_t)new_capacity * struct_size, alignment);
		} else {
			new_data = tfxREALLOCATE(data, (size_t)new_capacity * struct_size);
		}
		TFX_ASSERT(new_data);    //Unable to allocate memory. todo: better handling
		if (data) {
			memcpy(new_data, data, (size_t)current_size * struct_size);
			tfxFREE(data);
		}
		data = new_data;
		capacity = new_capacity;
	}
	inline tfxU32       _grow_capacity(tfxU32 size) const { tfxU32 new_capacity = capacity ? (capacity + capacity / 2) : 8; return new_capacity > size ? new_capacity : size; }
	inline void         add_space(tfxU32 extra_size) { if (current_size + extra_size > capacity) reserve(_grow_capacity(current_size + extra_size)); current_size += extra_size; }
	inline void         resize(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
	inline void         resize_bytes(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
	inline tfxU32		size_in_bytes() { return current_size * struct_size; }
	inline const tfxU32	size_in_bytes() const { return current_size * struct_size; }
};

tfxINTERNAL inline tfx_buffer_t tfxCreateBuffer(tfxU32 struct_size, tfxU32 alignment) {
	tfx_buffer_t buffer;
	buffer.struct_size = struct_size;
	buffer.alignment = alignment;
	return buffer;
}

tfxINTERNAL inline void tfxReconfigureBuffer(tfx_buffer_t *buffer, size_t new_struct_size) {
	if ((tfxU32)new_struct_size == buffer->struct_size) return;
	tfxU32 size_in_bytes = buffer->capacity * buffer->struct_size;
	tfxU32 new_capacity = size_in_bytes / (tfxU32)new_struct_size;
	buffer->capacity = new_capacity;
	buffer->struct_size = (tfxU32)new_struct_size;
}

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

	tfx_hasher_t hasher;
	tfx_vector_t<pair> map;
	tfx_vector_t<T> data;
	tfxU32 iterator_index;
	tfxU32 last_insert_index;

	tfx_storage_map_t() {}

	inline void init() { map.init(); data.init(); iterator_index = 0; }

	inline void reserve(tfxU32 size) {
		if (size > data.capacity) {
			map.reserve(size);
			data.reserve(size);
		}
	}

	//Insert a new T value into the storage
	inline tfxKey Insert(const char *name, const T &value) {
		tfxKey key = tfx_Hash(&hasher, name, strlen(name), 0);
		SetIndex(key, value);
		return key;
	}

	//Insert a new T value into the storage
	inline T &InsertWithLength(const char *name, tfxU32 length, const T &value) {
		tfxKey key = tfx_Hash(&hasher, name, length, 0);
		return SetIndex(key, value);
	}

	//Insert a new T value into the storage
	inline T &Insert(tfxKey key, const T &value) {
		return SetIndex(key, value);
	}

	//Insert a new T value into the storage
	inline T &InsertByInt(int name, const T &value) {
		tfxKey key = name;
		return SetIndex(key, value);
	}

	inline void Clear() {
		data.clear();
		map.clear();
	}

	inline tfxKey MakeKey(const char *name) {
		return tfx_Hash(&hasher, name, strlen(name), 0);
	}

	inline void FreeAll() {
		data.free();
		map.free();
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

	inline bool ValidName(const char *name, tfxU32 length = 0) {
		TFX_ASSERT(name);    //Can't search for anything that's null
		return GetIndex(name, length) > -1;
	}

	inline bool ValidKey(tfxKey key) {
		return GetIndex(key) > -1;
	}

	inline bool ValidIntName(tfxU32 name) {
		return GetIntIndex(name) > -1;
	}

	//Remove an item from the data. Slow function, 2 memmoves and then the map has to be iterated and indexes reduced by one
	//to re-align them
	inline void Remove(const char *name) {
		tfxKey key = tfx_Hash(&hasher, name, strlen(name), 0);
		pair *it = LowerBound(key);
		tfxU32 index = it->index;
		T *list_it = &data[index];
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
		tfxU32 index = it->index;
		T *list_it = &data[index];
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
		tfxU32 index = it->index;
		T *list_it = &data[index];
		map.erase(it);
		data.erase(list_it);
		for (auto &p : map) {
			if (p.index < index) continue;
			p.index--;
		}
	}

	inline T &At(const char *name) {
		int index = GetIndex(name);
		TFX_ASSERT(index > -1);                        //Key was not found
		return data[index];
	}

	inline T *AtPtr(const char *name) {
		int index = GetIndex(name);
		return index != -1 ? &data[index] : nullptr;
	}

	inline T *AtPtr(tfxKey key) {
		int index = GetIndex(key);
		return index != -1 ? &data[index] : nullptr;
	}

	inline T &AtInt(int name) {
		int index = GetIntIndex(name);
		TFX_ASSERT(index > -1);                        //Key was not found
		return data[index];
	}

	inline T &At(tfxKey key) {
		int index = GetIndex(key);
		TFX_ASSERT(index > -1);                        //Key was not found
		return data[index];
	}

	inline T &operator[](const tfxU32 index) {
		TFX_ASSERT(index < data.current_size);        //Index was out of range
		return data[index];
	}

	T &SetIndex(tfxKey key, const T &value) {
		pair *it = LowerBound(key);
		if (it == map.end() || it->key != key)
		{
			data.push_back(value);
			map.insert(it, pair(key, data.current_size - 1));
			last_insert_index = data.current_size - 1;
			return data.back();
		}
		last_insert_index = it->index;
		data[it->index] = value;
		return data[it->index];
	}

	int GetIndex(const char *name, tfxU32 length = 0) {
		tfxKey key = tfx_Hash(&hasher, name, length ? length : strlen(name), 0);
		pair *it = LowerBound(key);
		if (it == map.end() || it->key != key)
			return -1;
		return it->index;
	}

	int GetIntIndex(int name) {
		tfxKey key = name;
		pair *it = LowerBound(key);
		if (it == map.end() || it->key != key)
			return -1;
		return it->index;
	}

	int GetIndex(tfxKey key) {
		pair *it = LowerBound(key);
		if (it == map.end() || it->key != key)
			return -1;
		return it->index;
	}

	pair *LowerBound(tfxKey key)
	{
		tfx_storage_map_t::pair *first = map.data;
		tfx_storage_map_t::pair *last = map.data + map.current_size;
		size_t count = (size_t)(last - first);
		while (count > 0)
		{
			size_t count2 = count >> 1;
			tfx_storage_map_t::pair *mid = first + count2;
			if (mid->key < key)
			{
				first = ++mid;
				count -= count2 + 1;
			} else
			{
				count = count2;
			}
		}
		return first;
	}

	T *next_item() {
		T *item = nullptr;
		if (iterator_index < data.current_size) {
			item = &data[iterator_index];
			iterator_index++;
		} else {
			iterator_index = 0;
		}
		return item;
	}

	void reset_iterator() {
		iterator_index = 0;
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

//Note this doesn't free memory, call tfx__free_soa_buffer to do that.
inline void tfx__init_soa_buffer(tfx_soa_buffer_t *buffer) {
	buffer->current_arena_size = 0;
	buffer->struct_size = 0;
	buffer->current_size = 0;
	buffer->start_index = 0;
	buffer->last_bump = 0;
	buffer->capacity = 0;
	buffer->alignment = 4;
	buffer->block_size = tfxDataWidth;
	buffer->buffer_amount = 0;
	buffer->user_data = nullptr;
	buffer->resize_callback = nullptr;
	buffer->struct_of_arrays = nullptr;
	buffer->data = nullptr;
	buffer->array_ptrs.init();
}

inline void *tfx__get_end_of_buffer_ptr(tfx_soa_buffer_t *buffer) {
	TFX_ASSERT(buffer->data);
	return (char *)buffer->data + buffer->current_arena_size;
}

//Get the amount of free space in the buffer
inline tfxU32 tfx__free_sprite_buffer_space(tfx_soa_buffer_t *buffer) {
	return buffer->capacity - buffer->current_size;
}

//Get the index based on the buffer being a ring buffer
inline tfxU32 tfx__get_circular_index(tfx_soa_buffer_t *buffer, tfxU32 index) {
	return (buffer->start_index + index) % buffer->capacity;
}

//Convert a circular index back into an index from the start of the buffer
inline tfxU32 tfx__get_absolute_index(tfx_soa_buffer_t *buffer, tfxU32 circular_index) {
	return buffer->capacity - (circular_index % buffer->capacity);
}

//Add an array to a SoABuffer. parse in the size of the data type and the offset to the member variable within the struct.
//You must add all the member veriables in the struct before calling tfx__finish_soa_buffer_setup
inline void tfx__add_struct_array(tfx_soa_buffer_t *buffer, size_t unit_size, size_t offset) {
	tfx_soa_data_t data;
	data.unit_size = unit_size;
	data.offset = offset;
	buffer->array_ptrs.push_back(data);
}

//In order to ensure memory alignment of all arrays in the buffer we need the following function to get the correct amount
//of memory required to contain all the data in the buffer for each array in the struct of arrays.
inline size_t tfx__get_soa_capacity_requirement(tfx_soa_buffer_t *buffer, size_t capacity) {
	size_t size_requirement = 0;
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		size_requirement += buffer->array_ptrs[i].unit_size * capacity;
		size_t mod = size_requirement % buffer->alignment;
		size_requirement += mod ? buffer->alignment - mod : 0;
	}
	return size_requirement;
}

//Once you have called tfx__add_struct_array for all your member variables you must call this function in order to 
//set up the memory for all your arrays. One block of memory will be created and all your arrays will be line up
//inside the space
inline void tfx__finish_soa_buffer_setup(tfx_soa_buffer_t *buffer, void *struct_of_arrays, tfxU32 reserve_amount, tfxU32 alignment = 4, tfxU32 buffer_amount = 0) {
	TFX_ASSERT(buffer->data == nullptr && buffer->array_ptrs.current_size > 0);    //Must be an unitialised soa buffer
	TFX_ASSERT(alignment >= 4);        //Alignment must be 4 or greater
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		buffer->struct_size += buffer->array_ptrs[i].unit_size;
	}
	buffer->buffer_amount = buffer_amount;
	reserve_amount = (reserve_amount / buffer->block_size + 1) * buffer->block_size;
	buffer->capacity = reserve_amount;
	TFX_ASSERT(buffer->capacity > buffer->buffer_amount);	//buffer_amount must be less than the capacity, reserver more than the buffer_amount
	buffer->alignment = alignment;
	buffer->current_arena_size = tfx__get_soa_capacity_requirement(buffer, reserve_amount);
	buffer->data = tfxALLOCATE_ALIGNED(buffer->current_arena_size, buffer->alignment);
	TFX_ASSERT(buffer->data);    //Unable to allocate memory. Todo: better handling
	memset(buffer->data, 0, buffer->current_arena_size);
	buffer->struct_of_arrays = struct_of_arrays;
	size_t running_offset = 0;
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		buffer->array_ptrs[i].ptr = (char *)buffer->data + running_offset;
		memcpy((char *)buffer->struct_of_arrays + buffer->array_ptrs[i].offset, &buffer->array_ptrs[i].ptr, sizeof(void *));
		running_offset += buffer->array_ptrs[i].unit_size * buffer->capacity;
		size_t mod = running_offset % buffer->alignment;
		running_offset += mod ? buffer->alignment - mod : 0;
	}
	if (buffer->resize_callback) {
		buffer->resize_callback(buffer, 0);
	}
}

//Call this function to increase the capacity of all the arrays in the buffer. Data that is already in the arrays is preserved if keep_data passed as true (default).
inline bool tfx__grow_soa_arrays(tfx_soa_buffer_t *buffer, tfxU32 first_new_index, tfxU32 new_target_size, bool keep_data = true) {
	TFX_ASSERT(buffer->capacity);            //buffer must already have a capacity!
	tfxU32 new_capacity = 0;
	new_capacity = new_target_size > buffer->capacity ? new_target_size + new_target_size / 2 : buffer->capacity + buffer->capacity / 2;
	new_capacity = (new_capacity / buffer->block_size + 1) * buffer->block_size;
	void *new_data = tfxALLOCATE_ALIGNED(tfx__get_soa_capacity_requirement(buffer, new_capacity), buffer->alignment);
	TFX_ASSERT(new_data);    //Unable to allocate memory. Todo: better handling
	memset(new_data, 0, new_capacity * buffer->struct_size);
	size_t running_offset = 0;
	if (keep_data) {
		for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
			size_t capacity = buffer->capacity * buffer->array_ptrs[i].unit_size;
			size_t start_index = buffer->start_index * buffer->array_ptrs[i].unit_size;
			if ((buffer->start_index + buffer->current_size - 1) > buffer->capacity) {
				memcpy((char *)new_data + running_offset, (char *)buffer->array_ptrs[i].ptr + start_index, (size_t)(capacity - start_index));
				memcpy((char *)new_data + (capacity - start_index) + running_offset, (char *)buffer->array_ptrs[i].ptr, (size_t)(start_index));
			} else {
				memcpy((char *)new_data + running_offset, (char *)buffer->array_ptrs[i].ptr + start_index, (size_t)(capacity - start_index));
			}
			running_offset += buffer->array_ptrs[i].unit_size * new_capacity;
			size_t mod = running_offset % buffer->alignment;
			running_offset += mod ? buffer->alignment - mod : 0;
		}
	}
	void *old_data = buffer->data;

	buffer->data = new_data;
	buffer->capacity = new_capacity;
	buffer->current_arena_size = new_capacity * buffer->struct_size;
	running_offset = 0;
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		buffer->array_ptrs[i].ptr = (char *)buffer->data + running_offset;
		memcpy((char *)buffer->struct_of_arrays + buffer->array_ptrs[i].offset, &buffer->array_ptrs[i].ptr, sizeof(void *));
		running_offset += buffer->array_ptrs[i].unit_size * buffer->capacity;
		size_t mod = running_offset % buffer->alignment;
		running_offset += mod ? buffer->alignment - mod : 0;
	}
	tfxFREE(old_data);

	if (buffer->resize_callback) {
		buffer->resize_callback(buffer, first_new_index);
	}

	buffer->start_index = 0;

	return true;
}

//Increase current size of a SoA Buffer and grow if necessary.
inline void tfx__resize_soa_buffer(tfx_soa_buffer_t *buffer, tfxU32 new_size) {
	TFX_ASSERT(buffer->data);            //No data allocated in buffer
	if (new_size >= buffer->capacity) {
		tfx__grow_soa_arrays(buffer, buffer->capacity, new_size);
	}
	buffer->current_size = new_size;
}

//Copy a buffer
inline void tfx__copy_soa_buffer(tfx_soa_buffer_t *dst, tfx_soa_buffer_t *src) {
	memcpy(dst, src, sizeof(tfx_soa_buffer_t));
	dst->array_ptrs.init();
	dst->array_ptrs.copy(src->array_ptrs);
}

//Increase current size of a SoA Buffer and grow if necessary. This will not shrink the capacity so if new_size is not bigger than the
//current capacity then nothing will happen
inline void tfx__set_soa_capacity(tfx_soa_buffer_t *buffer, tfxU32 new_size) {
	TFX_ASSERT(buffer->data);            //No data allocated in buffer
	if (new_size >= buffer->capacity) {
		tfx__grow_soa_arrays(buffer, buffer->capacity, new_size);
	}
}

//Increase current size of a SoA Buffer by 1 and grow if grow is true. Returns the last index.
inline tfxU32 tfx__add_soa_row(tfx_soa_buffer_t *buffer, bool grow = false) {
	TFX_ASSERT(buffer->data);            //No data allocated in buffer
	TFX_ASSERT(buffer->capacity >= buffer->buffer_amount);            //Buffer amount must not be greater than capacity of the buffer
	tfxU32 new_size = ++buffer->current_size;
	if (grow && new_size >= buffer->capacity - buffer->buffer_amount) {
		tfx__grow_soa_arrays(buffer, buffer->capacity, new_size);
	}
	buffer->current_size = new_size;
	TFX_ASSERT(buffer->current_size <= buffer->capacity);    //Capacity of buffer is exceeded, set grow to true or don't exceed the capacity
	return buffer->current_size - 1;
}

//Increase current size of a SoA Buffer by a specific amount and grow if grow is true. Returns the last index.
//You can also pass in a boolean to know if the buffer had to be increased in size or not. Returns the index where the new rows start.
inline tfxU32 tfx__add_soa_rows_grew(tfx_soa_buffer_t *buffer, tfxU32 amount, bool grow, bool &grew) {
	TFX_ASSERT(buffer->data);            //No data allocated in buffer
	TFX_ASSERT(buffer->capacity >= buffer->buffer_amount);            //Buffer amount must not be greater than capacity of the buffer
	tfxU32 first_new_index = buffer->current_size;
	tfxU32 new_size = buffer->current_size += amount;
	if (grow && new_size >= buffer->capacity - buffer->buffer_amount) {
		grew = tfx__grow_soa_arrays(buffer, buffer->capacity, new_size);
	}
	buffer->current_size = new_size;
	TFX_ASSERT(buffer->current_size < buffer->capacity);    //Capacity of buffer is exceeded, set grow to true or don't exceed the capacity
	return first_new_index;
}

//Increase current size of a SoA Buffer and grow if grow is true. Returns the index where the new rows start.
inline tfxU32 tfx__add_soa_rows(tfx_soa_buffer_t *buffer, tfxU32 amount, bool grow) {
	TFX_ASSERT(buffer->data);            //No data allocated in buffer
	TFX_ASSERT(buffer->capacity >= buffer->buffer_amount);            //Buffer amount must not be greater than capacity of the buffer
	tfxU32 first_new_index = buffer->current_size;
	tfxU32 new_size = buffer->current_size + amount;
	if (grow && new_size >= buffer->capacity - buffer->buffer_amount) {
		tfx__grow_soa_arrays(buffer, buffer->capacity, new_size);
	}
	buffer->current_size = new_size;
	TFX_ASSERT(buffer->current_size < buffer->capacity);    //Capacity of buffer is exceeded, set grow to true or don't exceed the capacity
	return first_new_index;
}

//Decrease the current size of a SoA Buffer by 1.
inline void tfx__pop_soa_row(tfx_soa_buffer_t *buffer) {
	TFX_ASSERT(buffer->data && buffer->current_size > 0);            //No data allocated in buffer
	buffer->current_size--;
}

//tfx__bump_soa_buffer the start index of the SoA buffer (ring buffer usage)
inline void tfx__bump_soa_buffer(tfx_soa_buffer_t *buffer) {
	TFX_ASSERT(buffer->data && buffer->current_size > 0);            //No data allocated in buffer
	if (buffer->current_size == 0)
		return;
	buffer->start_index++; buffer->start_index %= buffer->capacity; buffer->current_size--;
}

//tfx__bump_soa_buffer the start index of the SoA buffer (ring buffer usage)
inline void tfx__bump_soa_buffer_amount(tfx_soa_buffer_t *buffer, tfxU32 amount) {
	TFX_ASSERT(buffer->data && buffer->current_size > 0);            //No data allocated in buffer
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
inline void tfx__free_soa_buffer(tfx_soa_buffer_t *buffer) {
	buffer->current_arena_size = buffer->current_size = buffer->capacity = 0;
	if (buffer->data) {
		tfxFREE(buffer->data);
	}
	buffer->array_ptrs.free();
	tfx__init_soa_buffer(buffer);
}

//Clear the SoA buffer
inline void tfx__clear_soa_buffer(tfx_soa_buffer_t *buffer) {
	buffer->current_size = buffer->start_index = 0;
}

//Trim an SoA buffer to the current size. This is a bit rough and ready and I just created it for trimming compressed sprite data down to size
inline void tfx__trim_soa_buffer(tfx_soa_buffer_t *buffer) {
	if (buffer->current_size == buffer->capacity) {
		return;
	}
	if (buffer->current_size == 0) {
		tfx__free_soa_buffer(buffer);
		return;
	}
	TFX_ASSERT(buffer->current_size < buffer->capacity);
	tfxU32 new_capacity = buffer->current_size + buffer->buffer_amount;
	void *new_data = tfxALLOCATE_ALIGNED(tfx__get_soa_capacity_requirement(buffer, new_capacity), buffer->alignment);
	TFX_ASSERT(new_data);    //Unable to allocate memory. Todo: better handling
	memset(new_data, 0, new_capacity * buffer->struct_size);
	size_t running_offset = 0;
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		size_t capacity = new_capacity * buffer->array_ptrs[i].unit_size;
		size_t start_index = buffer->start_index * buffer->array_ptrs[i].unit_size;
		if ((buffer->start_index + buffer->current_size - 1) > buffer->capacity) {
			memcpy((char *)new_data + running_offset, (char *)buffer->array_ptrs[i].ptr + start_index, (size_t)(capacity - start_index));
			memcpy((char *)new_data + (capacity - start_index) + running_offset, (char *)buffer->array_ptrs[i].ptr, (size_t)(start_index));
		} else {
			memcpy((char *)new_data + running_offset, (char *)buffer->array_ptrs[i].ptr + start_index, (size_t)(capacity - start_index));
		}
		running_offset += buffer->array_ptrs[i].unit_size * new_capacity;
		size_t mod = running_offset % buffer->alignment;
		running_offset += mod ? buffer->alignment - mod : 0;

	}
	void *old_data = buffer->data;

	buffer->data = new_data;
	buffer->capacity = new_capacity;
	buffer->current_arena_size = new_capacity * buffer->struct_size;
	running_offset = 0;
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		buffer->array_ptrs[i].ptr = (char *)buffer->data + running_offset;
		memcpy((char *)buffer->struct_of_arrays + buffer->array_ptrs[i].offset, &buffer->array_ptrs[i].ptr, sizeof(void *));
		running_offset += buffer->array_ptrs[i].unit_size * buffer->capacity;
		size_t mod = running_offset % buffer->alignment;
		running_offset += mod ? buffer->alignment - mod : 0;
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
	tfx_vector_t<tfx_bucket_t<T> *> bucket_list;

	tfx_bucket_array_t() : size_of_each_bucket(8), current_size(0), capacity(0), locked(0) {}

	inline void			 init() { size_of_each_bucket = 8; current_size = capacity = locked = 0; bucket_list.init(); }
	inline bool          empty() { return current_size == 0; }
	inline tfxU32        size() { return current_size; }
	inline void            free() {
		for (tfx_bucket_t<T> *bucket : bucket_list) {
			bucket->data.free();
			tfxFREE(bucket);
		}
		current_size = capacity = 0;
		bucket_list.free();
	}
	inline T &operator[](tfxU32 i) {
		TFX_ASSERT(i < current_size);        //Index is out of bounds
		tfxU32 bucket_index = i / size_of_each_bucket;
		tfxU32 element_index = i % size_of_each_bucket;
		return (*bucket_list[bucket_index]).data[element_index];
	}
	inline const T &operator[](tfxU32 i) const {
		TFX_ASSERT(i < current_size);        //Index is out of bounds
		tfxU32 bucket_index = i / size_of_each_bucket;
		tfxU32 element_index = i % size_of_each_bucket;
		return (*bucket_list[bucket_index]).data[element_index];
	}
	inline T *begin() { return bucket_list.current_size ? (T *)bucket_list[0]->data.data : nullptr; }
	inline const T *begin() const { return bucket_list.current_size ? (T *)bucket_list[0]->data.data : nullptr; }
	inline T *end() { return bucket_list.current_size ? (T *)bucket_list[(current_size - 1) / size_of_each_bucket]->data.end() : nullptr; }
	inline const T *end() const { return bucket_list.current_size ? (T *)bucket_list[(current_size - 1) / size_of_each_bucket]->data.end() : nullptr; }
	inline T &front() { TFX_ASSERT(current_size > 0); return bucket_list[0]->data.front(); }
	inline const T &front() const { TFX_ASSERT(current_size > 0); return bucket_list[0]->data.front(); }
	inline T &back() { TFX_ASSERT(current_size > 0); 
		tfxU32 bucket = (current_size - 1) / size_of_each_bucket; 
		return bucket_list[bucket]->data.back(); }
	inline const T &back() const { TFX_ASSERT(current_size > 0); return bucket_list[(current_size - 1) / size_of_each_bucket]->data.back(); }
	inline tfxU32        active_buckets() { return current_size == 0 ? 0 : current_size / size_of_each_bucket + 1; }
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

	inline T &push_back(const T &v) {
		if (current_size == capacity) {
			add_bucket(size_of_each_bucket);
			capacity += size_of_each_bucket;
		}
		tfxU32 current_bucket = current_size / size_of_each_bucket;
		bucket_list[current_bucket]->data.push_back(v);
		current_size++;
		return bucket_list[current_bucket]->data.back();
	}

	inline tfxU32 locked_push_back(const T &v) {
		while (tfx__compare_and_exchange(&locked, 1, 0) > 1);

		push_back(v);

		tfx__exchange(&locked, 0);
		return current_size - 1;
	}

	inline T *insert(tfxU32 insert_index, const T &v) {
		TFX_ASSERT(insert_index < current_size);
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
				} else {
					end_element = bucket_list[current_insert_bucket]->data.pop_back();
					end_pushed = false;
					bucket_list[current_insert_bucket]->data.push_front(end_element2);
					end2_pushed = true;
				}
				alternator = alternator ^ 1;
			} else {
				bucket_list[current_insert_bucket]->data.push_front(alternator == 0 ? end_element : end_element2);
				end_pushed = true;
				end2_pushed = true;
				break;
			}
		}
		if (!end_pushed || !end2_pushed) {
			push_back(!end_pushed ? end_element : end_element2);
		} else {
			current_size++;
		}
		return &bucket_list[insert_bucket]->data[element_index];
	}

	inline T *insert(T *position, const T &v) {
		tfxU32 index = 0;
		bool find_result = find(position, index);
		TFX_ASSERT(find_result);    //Could not find the object to insert at, make sure it exists
		return insert(index, v);
	}

	inline void erase(tfxU32 index) {
		TFX_ASSERT(index < current_size);
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

	inline void erase(T *it) {
		tfxU32 index = 0;
		bool find_result = find(it, index);
		TFX_ASSERT(find_result);    //pointer not found in list
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

	inline T *find(T *it) {
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
			if (bucket_list.current_size) {
				bucket_list.back()->data.current_size = current_size == size_of_each_bucket ? current_size : (current_size % size_of_each_bucket);
			}
		} else {
			if (bucket_list.current_size) {
				bucket_list.back()->data.current_size = current_size == size_of_each_bucket ? current_size : (current_size % size_of_each_bucket);
			}
		}
	}

};

template <typename T>
inline void tfxInitBucketArray(tfx_bucket_array_t<T> *bucket_array, tfxU32 size_of_each_bucket) {
	bucket_array->current_size = bucket_array->locked = bucket_array->capacity = 0;
	bucket_array->size_of_each_bucket = size_of_each_bucket;
}

template <typename T>
inline void tfxCopyBucketArray(tfx_bucket_array_t<T> *dst, tfx_bucket_array_t<T> *src) {
	if (src == dst) {
		return;
	}
	dst->free();
	for (tfx_bucket_t<T> *bucket : src->bucket_list) {
		tfx_bucket_t<T> *copy = dst->add_bucket(src->size_of_each_bucket);
		copy->data.copy(bucket->data);
	}
	dst->current_size = src->current_size;
	dst->capacity = src->capacity;
	dst->size_of_each_bucket = src->size_of_each_bucket;
}

#define tfxBucketLoop(bucket, index) int index = 0; index != bucket.current_size; ++index

struct tfx_line_t {
	const char *start;
	const char *end;
	int length;
};

tfxAPI_EDITOR inline int tfx_FindInLine(tfx_line_t *line, const char *needle) {
	size_t needle_length = strlen(needle);
	if (needle_length > line->length) return -1;
	tfxU32 pos = 0;
	while (needle_length + pos <= line->length) {
		if (strncmp(line->start + pos, needle, needle_length) == 0) {
			return pos;
		}
		++pos;
	}
	return -1;
}

#endif

#ifdef __cplusplus
//A char buffer you can use to load a file into and read from
//Has no deconstructor so make sure you call Free() when done
//This is meant for limited usage in timeline fx only and not recommended for use outside!
typedef struct tfx_stream_s {
	tfxU32 magic;
	tfxU64 size;
	tfxU64 capacity;
	tfxU64 position;
	char *data;

	inline tfx_stream_s() { size = position = capacity = 0; data = nullptr; }
	inline tfx_stream_s(tfxU64 qty) { size = position = capacity = 0; data = nullptr; Resize(qty); }

	inline bool Read(char *dst, tfxU64 count) {
		if (count + position <= size) {
			memcpy(dst, data + position, count);
			position += count;
			return true;
		}
		return false;
	}
	tfx_line_t ReadLine();
	inline bool Write(void *src, tfxU64 count) {
		if (count + position <= size) {
			memcpy(data + position, src, count);
			position += count;
			return true;
		}
		return false;
	}
	void AddLine(const char *format, ...);
	void Appendf(const char *format, ...);
	void Appendv(const char *format, va_list args);
	void SetText(const char *text);
	inline void Append(char c) { if (size) { size--; } data[size] = c; size++; NullTerminate(); }
	inline bool EoF() { return position >= Length(); }
	inline void AddReturn() { if (size + 1 >= capacity) { tfxU64 new_capacity = capacity * 2; Reserve(new_capacity); } Append('\n'); }
	inline void Seek(tfxU64 offset) {
		if (offset < size)
			position = offset;
		else
			position = size;
	}
	inline void TrimToZero() {
		if (size < capacity) {
			bool result = tfx_SafeMemset(data, data + size, '\0', capacity - size);
			TFX_ASSERT(result);
		}
	}

	inline bool            Empty() { return size == 0; }
	inline tfxU64        Size() { return size; }
	inline const tfxU64    Size() const { return size; }
	inline tfxU64        Length() { if (!data) { return 0; } return data[size] == '\0' && size > 0 ? size - 1 : size; }
	inline const tfxU64    Length() const { if (!data) { return 0; } return data[size] == '\0' && size > 0 ? size - 1 : size; }

	inline void            Free() { if (data) { size = position = capacity = 0; tfxFREE(data); data = nullptr; } }
	inline void         Clear() { if (data) { size = 0; } }

	inline void         GrowSize(tfxU64 extra_size) {
		Resize(size + extra_size);
	}
	inline void         Reserve(tfxU64 new_capacity) {
		if (new_capacity <= capacity) {
			return;
		}
		char *new_data = (char *)tfxALLOCATE_ALIGNED((tfxU64)new_capacity * sizeof(char), 16);
		memset(new_data, 0, new_capacity);
		TFX_ASSERT(new_data);    //Unable to allocate memory. Todo: better handling
		if (data) {
			memcpy(new_data, data, (tfxU64)size * sizeof(char));
			tfxFREE(data);
		}
		data = new_data;
		capacity = new_capacity;
	}
	inline void         Resize(tfxU64 new_size) {
		if (new_size <= capacity) {
			size = new_size;
			return;
		}
		Reserve(new_size);
		size = new_size;
		position = 0;
	}
	inline void            NullTerminate() { *(data + size++) = '\0'; }

} tfx_stream_t;

tfxMAKE_HANDLE(tfx_stream)

tfxAPI_EDITOR inline tfx_stream tfx__create_stream() {
	tfx_stream_t blank_stream = { 0 };
	tfx_stream stream = tfxNEW(tfx_stream);
	*stream = blank_stream;
	stream->magic = tfxINIT_MAGIC;
	return stream;
}

tfxAPI_EDITOR inline void tfx_FreeStream(tfx_stream stream) {
	stream->Free();
	tfxFREE(stream);
}

#else

#define tfx_str_(size)	\
typedef struct tfx_str##size##_s {	\
	char data[size];	\
	tfxU32 current_size;	\
}tfx_str##size##_t

tfx_str_(32);
tfx_str_(64);
tfx_str_(128);
tfx_str_(256);
tfx_str_(512);

typedef struct tfx_storage_map_s {
	tfx_hasher_t hasher;
	tfx_vector_t map;
	tfx_vector_t data;
} tfx_storage_map_t;

typedef struct tfx_buffer_s {
	tfxU32 current_size;
	tfxU32 capacity;
	tfxU32 struct_size;
	tfxU32 alignment;
	void *data;
} tfx_buffer_t;

typedef struct tfx_bucket_array_s {
	tfxU32 current_size;
	tfxU32 capacity;
	tfxU32 size_of_each_bucket;
	tfxLONG volatile locked;
	tfx_vector_t bucket_list;
} tfx_bucket_array_t;

typedef struct tfx_str_s tfx_str_t;
typedef struct tfx_storage_map_s tfx_storage_map_t;
typedef struct tfx_soa_data_s tfx_soa_data_t;
typedef struct tfx_soa_buffer_s tfx_soa_buffer_t;
typedef struct tfx_bucket_s tfx_bucket_t;
typedef struct tfx_line_s tfx_line_t;
typedef struct tfx_stream_s {
	tfxU32 magic;
	tfxU64 size;
	tfxU64 capacity;
	tfxU64 position;
	char *data;
} tfx_stream_t;

tfxMAKE_HANDLE(tfx_stream)

#endif	//__cplusplus

inline tfxU32 tfxIsPowerOf2(tfxU32 v)
{
	return ((v & ~(v - 1)) == v);
}

//-----------------------------------------------------------
//Section: Global_Variables
//-----------------------------------------------------------
// Platform-specific synchronization wrapper
extern const tfxU32 tfxPROFILE_COUNT;

extern int tfxNumberOfThreadsInAdditionToMain;

#ifndef tfxMAX_QUEUES
#define tfxMAX_QUEUES 64
#endif

#ifndef tfxMAX_QUEUE_ENTRIES
#define tfxMAX_QUEUE_ENTRIES 512
#endif

typedef struct tfx_work_queue_s tfx_work_queue_t;

// Callback function type
typedef void (*tfx_work_queue_callback)(struct tfx_work_queue_s *queue, void *data);

typedef struct tfx_work_queue_entry_s {
	tfx_work_queue_callback call_back;
	void *data;
} tfx_work_queue_entry_t;

typedef struct tfx_work_queue_s {
	volatile tfx_uint entry_completion_goal;
	volatile tfx_uint entry_completion_count;
	volatile int next_read_entry;
	volatile int next_write_entry;
	tfx_work_queue_entry_t entries[tfxMAX_QUEUE_ENTRIES];
} tfx_work_queue_t;

typedef struct tfx_queue_processor_s {
	tfx_sync_t sync;
	tfx_uint count;
	volatile bool end_all_threads;
	tfx_work_queue_t *queues[tfxMAX_QUEUES];
} tfx_queue_processor_t;

typedef struct tfx_data_types_dictionary_s {
	int initialised;
#ifdef __cplusplus
	tfx_storage_map_t<tfx_data_type> names_and_types;
#else
	tfx_storage_map_t names_and_types;
#endif
} tfx_data_types_dictionary_t;

tfxAPI_EDITOR void tfx__initialise_dictionary(tfx_data_types_dictionary_t *dictionary);
tfxAPI_EDITOR void tfx__initialise_graph_indexes();

//Global variables
typedef struct tfx_storage_s {
	tfxU32 memory_pool_count;
	size_t default_memory_pool_size;
	size_t memory_pool_sizes[tfxMAX_MEMORY_POOLS];
	tfx_pool *memory_pools[tfxMAX_MEMORY_POOLS];
	tfx_data_types_dictionary_t data_types;
	float circle_path_x[tfxCIRCLENODES];
	float circle_path_z[tfxCIRCLENODES];
	tfx_color_ramp_format color_ramp_format;
	tfx_queue_processor_t thread_queues;
	tfx_hasher_t hasher;
#ifdef __cplusplus
	tfx_storage_map_t<tfx_effect_manager> effect_managers;
	tfx_storage_map_t<tfx_library> libraries;
	tfx_storage_map_t<tfxU32> graph_indexes;
#else
	tfx_storage_map_t effect_managers;
	tfx_storage_map_t libraries;
	tfx_storage_map_t graph_indexes;
#endif
	tfx_buffer_t gpu_graph_data;
	tfx_ribbon_dispatch last_ribbon_dispatch;
	tfx_effect_manager current_pm;
	tfx_ribbon_buffer_requirements ribbon_buffer_requirements;
#ifdef _WIN32
	HANDLE threads[tfxMAX_THREADS];
#else
	pthread_t threads[tfxMAX_THREADS];
#endif
	tfxU32 thread_count;
} tfx_storage_t;

extern tfx_storage_t *tfxStore;
extern tfx_allocator *tfxMemoryAllocator;

//-----------------------------------------------------------
//Section: Multithreading_Work_Queues
//-----------------------------------------------------------

//Tried to keep this as simple as possible, was originally based on Casey Muratory's Hand Made Hero threading which used the Windows API for
//threading but for the sake of supporting other platforms I changed it to use std::thread which was actually a lot more simple to do then 
//I expected. I just had to swap the semaphores for condition_variable and that was pretty much it other then obviously using std::thread as well.
//There is a single thread pool created to serve multiple queues. Currently each effect manager that you create will have it's own queue and then
//each emitter that the effect manager uses will be given it's own thread.

// Platform-specific atomic operations
tfxINTERNAL inline tfxU32 tfx__atomic_increment(volatile tfxU32 *value) {
#ifdef _WIN32
	return InterlockedIncrement((LONG *)value);
#else
	return __sync_fetch_and_add(value, 1) + 1;
#endif
}

tfxINTERNAL inline int tfx__atomic_compare_exchange(volatile int *dest, int exchange, int comparand) {
#ifdef _WIN32
	return InterlockedCompareExchange((LONG *)dest, exchange, comparand) == comparand;
#else
	return __sync_bool_compare_and_swap(dest, comparand, exchange);
#endif
}

tfxINTERNAL inline void tfx__memory_barrier(void) {
#ifdef _WIN32
	MemoryBarrier();
#else
	__sync_synchronize();
#endif
}

// Initialize queue processor
tfxINTERNAL inline void tfx__initialise_thread_queues(tfx_queue_processor_t *queues) {
	queues->count = 0;
	memset(queues->queues, 0, tfxMAX_QUEUES * sizeof(void *));
	tfx__sync_init(&queues->sync);
	queues->end_all_threads = 0;
}

tfxINTERNAL inline void tfx__cleanup_thread_queues(tfx_queue_processor_t *queues) {
	tfx__sync_cleanup(&queues->sync);
}

tfxINTERNAL inline tfx_work_queue_t *tfx__get_queue_with_work(tfx_queue_processor_t *thread_processor) {
	tfx__sync_lock(&thread_processor->sync);

	while (thread_processor->count == 0 && !thread_processor->end_all_threads) {
		tfx__sync_wait_full(&thread_processor->sync);
	}

	tfx_work_queue_t *queue = NULL;
	if (thread_processor->count > 0) {
		queue = thread_processor->queues[--thread_processor->count];
		tfx__sync_signal_empty(&thread_processor->sync);
	}

	tfx__sync_unlock(&thread_processor->sync);
	return queue;
}

tfxINTERNAL inline void tfx__push_queue_work(tfx_queue_processor_t *thread_processor, tfx_work_queue_t *queue) {
	tfx__sync_lock(&thread_processor->sync);

	while (thread_processor->count >= tfxMAX_QUEUES) {
		tfx__sync_wait_empty(&thread_processor->sync);
	}

	thread_processor->queues[thread_processor->count++] = queue;
	tfx__sync_signal_full(&thread_processor->sync);

	tfx__sync_unlock(&thread_processor->sync);
}

tfxINTERNAL inline bool tfx__do_next_work_queue(tfx_queue_processor_t *queue_processor) {
	tfx_work_queue_t *queue = tfx__get_queue_with_work(queue_processor);

	if (queue) {
		tfxU32 original_read_entry = queue->next_read_entry;
		tfxU32 new_original_read_entry = (original_read_entry + 1) % tfxMAX_QUEUE_ENTRIES;

		if (original_read_entry != queue->next_write_entry) {
			if (tfx__atomic_compare_exchange(&queue->next_read_entry, new_original_read_entry, original_read_entry)) {
				tfx_work_queue_entry_t entry = queue->entries[original_read_entry];
				entry.call_back(queue, entry.data);
				tfx__atomic_increment(&queue->entry_completion_count);
			}
		}
	}
	return queue_processor->end_all_threads;
}

tfxINTERNAL inline void tfx__do_next_work_queue_entry(tfx_work_queue_t *queue) {
	tfxU32 original_read_entry = queue->next_read_entry;
	tfxU32 new_original_read_entry = (original_read_entry + 1) % tfxMAX_QUEUE_ENTRIES;

	if (original_read_entry != queue->next_write_entry) {
		if (tfx__atomic_compare_exchange(&queue->next_read_entry, new_original_read_entry, original_read_entry)) {
			tfx_work_queue_entry_t entry = queue->entries[original_read_entry];
			entry.call_back(queue, entry.data);
			tfx__atomic_increment(&queue->entry_completion_count);
		}
	}
}

tfxINTERNAL inline void tfx__add_work_queue_entry(tfx_work_queue_t *queue, void *data, tfx_work_queue_callback call_back) {
	if (!tfxStore->thread_count) {
		call_back(queue, data);
		return;
	}

	tfxU32 new_entry_to_write = (queue->next_write_entry + 1) % tfxMAX_QUEUE_ENTRIES;
	while (new_entry_to_write == queue->next_read_entry) {        //Not enough room in work queue
		//We can do this because we're single producer
		tfx__do_next_work_queue_entry(queue);
	}
	queue->entries[queue->next_write_entry].data = data;
	queue->entries[queue->next_write_entry].call_back = call_back;
	tfx__atomic_increment(&queue->entry_completion_goal);

	tfx__writebarrier;

	tfx__push_queue_work(&tfxStore->thread_queues, queue);
	queue->next_write_entry = new_entry_to_write;

}

tfxINTERNAL inline void tfx__complete_all_work(tfx_work_queue_t *queue) {
	tfx_work_queue_entry_t entry = { 0 };
	while (queue->entry_completion_goal != queue->entry_completion_count) {
		tfx__do_next_work_queue_entry(queue);
	}
	queue->entry_completion_count = 0;
	queue->entry_completion_goal = 0;
}

#ifdef _WIN32
unsigned WINAPI tfx__thread_worker(void *arg);
#else
void *tfx__thread_worker(void *arg);
#endif

// Thread creation helper function
tfxINTERNAL inline int tfx__create_worker_thread(tfx_storage_t * storage, int thread_index) {
#ifdef _WIN32
	storage->threads[thread_index] = (HANDLE)_beginthreadex(
		NULL,
		0,
		tfx__thread_worker,
		&storage->thread_queues,
		0,
		NULL
	);
	return storage->threads[thread_index] != NULL;
#else
	return pthread_create(
		&storage->threads[thread_index],
		NULL,
		tfx__thread_worker,
		&storage->thread_queues
	) == 0;
#endif
}

// Thread cleanup helper function
tfxINTERNAL inline void tfx__cleanup_thread(tfx_storage_t * storage, int thread_index) {
#ifdef _WIN32
	if (storage->threads[thread_index]) {
		WaitForSingleObject(storage->threads[thread_index], INFINITE);
		CloseHandle(storage->threads[thread_index]);
		storage->threads[thread_index] = NULL;
	}
#else
	pthread_join(storage->threads[thread_index], NULL);
#endif
}

tfxAPI inline unsigned int tfx_HardwareConcurrency(void) {
#ifdef _WIN32
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	return sysinfo.dwNumberOfProcessors;
#else
#ifdef _SC_NPROCESSORS_ONLN
	long count = sysconf(_SC_NPROCESSORS_ONLN);
	return (count > 0) ? (unsigned int)count : 0;
#else
	return 0;
#endif
#endif
}

// Safe version that always returns at least 1
tfxAPI unsigned int tfx_HardwareConcurrencySafe(void);

// Helper function to get a good default thread count for thread pools
// Usually hardware threads - 1 to leave a core for the OS/main thread
tfxAPI unsigned int tfx_GetDefaultThreadCount(void);
//-----------------------------------------------------------
//Section: Vector_Math
//-----------------------------------------------------------

#ifdef __cplusplus
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

struct tfx_vec3_t {
	union {
		struct { float x, y, z; };
		struct { float pitch, yaw, roll; };
	};

	tfx_vec3_t() { x = y = z = 0.f; }
	tfx_vec3_t(float v) : x(v), y(v), z(v) {}
	tfx_vec3_t(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
	inline void operator=(const tfx_vec2_t &v) { x = v.x; y = v.y; }

	inline tfx_vec2_t xy() const { return tfx_vec2_t(x, y); }

	inline bool operator==(const tfx_vec3_t &v) const { return x == v.x && y == v.y && z == v.z; }
	inline bool operator!=(const tfx_vec3_t &v) const { return x != v.x || y != v.y || z != v.z; }

	inline tfx_vec3_t operator+(const tfx_vec3_t &v) const { return tfx_vec3_t(x + v.x, y + v.y, z + v.z); }
	inline tfx_vec3_t operator-(const tfx_vec3_t &v) const { return tfx_vec3_t(x - v.x, y - v.y, z - v.z); }
	inline tfx_vec3_t operator*(const tfx_vec3_t &v) const { return tfx_vec3_t(x * v.x, y * v.y, z * v.z); }
	inline tfx_vec3_t operator/(const tfx_vec3_t &v) const { return tfx_vec3_t(x / v.x, y / v.y, z / v.z); }

	inline tfx_vec3_t operator-() const { return tfx_vec3_t(-x, -y, -z); }

	inline void operator-=(const tfx_vec3_t &v) { x -= v.x; y -= v.y; z -= v.z; }
	inline void operator+=(const tfx_vec3_t &v) { x += v.x; y += v.y; z += v.z; }
	inline void operator*=(const tfx_vec3_t &v) { x *= v.x; y *= v.y; z *= v.z; }
	inline void operator/=(const tfx_vec3_t &v) { x /= v.x; y /= v.y; z /= v.z; }

	inline void operator-=(const tfx_vec2_t &v) { x -= v.x; y -= v.y; }
	inline void operator+=(const tfx_vec2_t &v) { x += v.x; y += v.y; }
	inline void operator*=(const tfx_vec2_t &v) { x *= v.x; y *= v.y; }
	inline void operator/=(const tfx_vec2_t &v) { x /= v.x; y /= v.y; }

	inline tfx_vec3_t operator+(float v) const { return tfx_vec3_t(x + v, y + v, z + v); }
	inline tfx_vec3_t operator-(float v) const { return tfx_vec3_t(x - v, y - v, z - v); }
	inline tfx_vec3_t operator*(float v) const { return tfx_vec3_t(x * v, y * v, z * v); }
	inline tfx_vec3_t operator/(float v) const { return tfx_vec3_t(x / v, y / v, z / v); }

	inline void operator*=(float v) { x *= v; y *= v; z *= v; }
	inline void operator/=(float v) { x /= v; y /= v; z /= v; }
	inline void operator+=(float v) { x += v; y += v; z += v; }
	inline void operator-=(float v) { x -= v; y -= v; z -= v; }

	inline float Squared() { return x * x + y * y + z * z; }
	inline bool IsNill() { return !x && !y && !z; }
};

struct tfx_vec4_t {
	union {
		struct { float x, y, z, w; };
		struct { float c0, c1, c2, c3; };
	};

	tfx_vec4_t() { x = y = z = w = 0.f; }
	tfx_vec4_t(float _x, float _y, float _z, float _w) : x(_x), y(_y), z(_z), w(_w) {}
	tfx_vec4_t(tfx_vec2_t vec1, tfx_vec2_t vec2) : x(vec1.x), y(vec1.y), z(vec2.x), w(vec2.y) {}
	tfx_vec4_t(tfx_vec3_t vec) : x(vec.x), y(vec.y), z(vec.z), w(0.f) {}
	tfx_vec4_t(tfx_vec3_t vec, float _w) : x(vec.x), y(vec.y), z(vec.z), w(_w) {}

	inline tfx_vec2_t xy() { return tfx_vec2_t(x, y); }
	inline tfx_vec2_t zw() { return tfx_vec2_t(z, w); }
	inline tfx_vec3_t xyz() { return tfx_vec3_t(x, y, z); }

	inline tfx_vec2_t xy() const { return tfx_vec2_t(x, y); }
	inline tfx_vec2_t zw() const { return tfx_vec2_t(z, w); }
	inline tfx_vec3_t xyz() const { return tfx_vec3_t(x, y, z); }

	inline void operator=(const tfx_vec2_t &v) { x = v.x; y = v.y; }
	inline void operator=(const tfx_vec3_t &v) { x = v.x; y = v.y; z = v.z; }

	inline tfx_vec4_t operator+(const tfx_vec4_t &v) { return tfx_vec4_t(x + v.x, y + v.y, z + v.z, w + v.w); }
	inline tfx_vec4_t operator-(const tfx_vec4_t &v) { return tfx_vec4_t(x - v.x, y - v.y, z - v.z, w - v.w); }
	inline tfx_vec4_t operator-() { return tfx_vec4_t(-x, -y, -z, -w); }
	inline tfx_vec4_t operator*(const tfx_vec4_t &v) { return tfx_vec4_t(x * v.x, y * v.y, z * v.z, w * v.w); }
	inline tfx_vec4_t operator/(const tfx_vec4_t &v) { return tfx_vec4_t(x / v.x, y / v.y, z / v.z, w / v.w); }

	inline void operator-=(const tfx_vec4_t &v) { x -= v.x; y -= v.y; z -= v.z; w -= v.w; }
	inline void operator+=(const tfx_vec4_t &v) { x += v.x; y += v.y; z += v.z; w += v.w; }
	inline void operator*=(const tfx_vec4_t &v) { x *= v.x; y *= v.y; z *= v.z; w *= v.w; }
	inline void operator/=(const tfx_vec4_t &v) { x /= v.x; y /= v.y; z /= v.z; w /= v.w; }

	inline void operator-=(const tfx_vec3_t &v) { x -= v.x; y -= v.y; z -= v.z; }
	inline void operator+=(const tfx_vec3_t &v) { x += v.x; y += v.y; z += v.z; }
	inline void operator*=(const tfx_vec3_t &v) { x *= v.x; y *= v.y; z *= v.z; }
	inline void operator/=(const tfx_vec3_t &v) { x /= v.x; y /= v.y; z /= v.z; }

	inline void operator-=(const tfx_vec2_t &v) { x -= v.x; y -= v.y; }
	inline void operator+=(const tfx_vec2_t &v) { x += v.x; y += v.y; }
	inline void operator*=(const tfx_vec2_t &v) { x *= v.x; y *= v.y; }
	inline void operator/=(const tfx_vec2_t &v) { x /= v.x; y /= v.y; }

	inline tfx_vec4_t operator+(float v) const { return tfx_vec4_t(x + v, y + v, z + v, w + v); }
	inline tfx_vec4_t operator-(float v) const { return tfx_vec4_t(x - v, y - v, z - v, w - v); }
	inline tfx_vec4_t operator*(float v) const { return tfx_vec4_t(x * v, y * v, z * v, w * v); }
	inline tfx_vec4_t operator/(float v) const { return tfx_vec4_t(x / v, y / v, z / v, w / v); }

	inline void operator*=(float v) { x *= v; y *= v; z *= v; w *= v; }
	inline void operator/=(float v) { x /= v; y /= v; z /= v; w /= v; }
	inline void operator+=(float v) { x += v; y += v; z += v; w += v; }
	inline void operator-=(float v) { x -= v; y -= v; z -= v; w -= v; }
};

struct tfx_quaternion_t {
	float w, x, y, z;

	tfx_quaternion_t() : w(0.f), x(0.f), y(0.f), z(0.f) {}
	tfx_quaternion_t(float w, float x, float y, float z) : w(w), x(x), y(y), z(z) {}

	tfx_vec3_t xyz() { return tfx_vec3_t(x, y, z); }
	tfx_vec4_t vec4() { return tfx_vec4_t(x, y, z, w); }
	tfx_quaternion_t operator*(const tfx_quaternion_t &q) const {
		return tfx_quaternion_t(
			w * q.w - x * q.x - y * q.y - z * q.z,
			w * q.x + x * q.w + y * q.z - z * q.y,
			w * q.y - x * q.z + y * q.w + z * q.x,
			w * q.z + x * q.y - y * q.x + z * q.w
		);
	}

};

#else

typedef struct tfx_vec2_s {
	float x, y;
} tfx_vec2_t;

typedef struct tfx_vec3_s {
	union {
		struct { float x, y, z; };
		struct { float pitch, yaw, roll; };
	};
} tfx_vec3_t;

typedef struct tfx_vec4_s {
	union {
		struct { float x, y, z, w; };
		struct { float c0, c1, c2, c3; };
	};
} tfx_vec4_t;

typedef struct tfx_quaternion_s {
	float w, x, y, z;
} tfx_quaternion_t;

#endif

tfxINTERNAL void tfx__to_quaternion2d(tfx_quaternion_t * q, float angle);
tfxINTERNAL tfx_vec2_t tfx__rotate_vector_quaternion2d(tfx_quaternion_t * q, tfx_vec2_t v);
tfxAPI_EDITOR tfx_vec3_t tfx__rotate_vector_quaternion(tfx_quaternion_t * q, tfx_vec3_t v);
tfxINTERNAL tfx_quaternion_t tfx__normalize_quaternion(tfx_quaternion_t * q);
tfxAPI_EDITOR tfx_quaternion_t tfx__euler_to_quaternion(float pitch, float yaw, float roll);
tfxAPI_EDITOR void tfx__wide_euler_to_packed_quaternion(tfxWideFloat pitch, tfxWideFloat yaw, tfxWideFloat roll, tfxWideInt *out_xy, tfxWideInt *out_zw);
tfxINTERNAL tfx_quaternion_t tfx__quaternion_from_axis_angle(float x, float y, float z, float angle);
tfxINTERNAL tfx_quaternion_t tfx__quaternion_from_direction(tfx_vec3_t * normalised_dir);

//Note, has padding for the sake of alignment on GPU compute shaders
typedef struct tfx_bounding_box_s {
	tfx_vec3_t min_corner; float padding1;
	tfx_vec3_t max_corner; float padding2;
} tfx_bounding_box_t;

const float tfxONE_DIV_255 = 1 / 255.f;
const float TFXONE_DIV_511 = 1 / 511.f;

typedef struct tfx_rgba8_s {
	union {
		struct {
			tfxU32 r : 8;
			tfxU32 g : 8;
			tfxU32 b : 8;
			tfxU32 a : 8;
		};
		struct { tfxU32 color; };
	};
} tfx_rgba8_t;
#define tfx_SWIZZLE_RGBA_TO_BGRA(rgba) ((rgba & 0xFF000000) | ((rgba & 0x00FF0000) >> 16) | (rgba & 0x0000FF00) | ((rgba & 0x000000FF) << 16))

typedef struct tfx_rgb_s {
	float r, g, b;
} tfx_rgb_t;

typedef struct tfx_hsv_s {
	float h, s, v;
} tfx_hsv_t;

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
			tfx_byte x : 8;
			tfx_byte y : 8;
			tfx_byte z : 8;
			tfx_byte w : 8;
		};
		struct { tfxU32 packed; };
	};
} tfx_float8x4_t;

const tfxWideArray one_div_127_wide = tfxWideSetConst(1 / 127.f);
const tfxWideArray one_div_511_wide = tfxWideSetConst(1 / 511.f);
const tfxWideArray one_div_32767_wide = tfxWideSetConst(1 / 32767.f);

#define tfxPACKED_X_NORMAL_3D 0x3FE7FDFF
#define tfxPACKED_Y_NORMAL_3D 0x1FFFF9FF
#define tfxPACKED_Z_NORMAL_3D 0x1FF7FFFE
#define tfxPACKED_Y_NORMAL_2D 32767
#define tfxPACKED_W_QUATERNION 0x7FFF000000000000ULL
typedef struct tfx_rgba_s {
	float r, g, b, a;
} tfx_rgba_t;

tfxAPI_EDITOR inline tfx_rgba_t tfx__create_rgba() { tfx_rgba_t color = { 1.f, 1.f, 1.f, 1.f }; return color; }
tfxAPI_EDITOR inline tfx_rgba_t tfx__create_rgba_from_rgba8(tfx_rgba8_t c) { tfx_rgba_t color = { (float)c.r * tfxONE_DIV_255, (float)c.g * tfxONE_DIV_255, (float)c.b * tfxONE_DIV_255, (float)c.a * tfxONE_DIV_255 }; return color; }

typedef struct tfx_mat4_s {
	tfx_vec4_t v[4];
} tfx_mat4_t TFX_ALIGN_AFFIX(16);

tfxINTERNAL inline void tfx__mat4_set2(tfx_mat4_t * mat, float aa, float ab, float ba, float bb) {
	mat->v[0].c0 = aa; mat->v[0].c1 = ab;
	mat->v[1].c0 = ba; mat->v[1].c1 = bb;
}

struct tfx_wide_mat4_t {
	float x[4];
	float y[4];
	float z[4];
	float w[4];
} TFX_ALIGN_AFFIX(16);;

struct tfx_mat3_t {
	tfx_vec3_t v[3];
};

//-----------------------------------------------------------
//Section: Simplex_Noise
//-----------------------------------------------------------

const float gradX[] =
{
	1,-1, 1,-1,
	1,-1, 1,-1,
	0, 0, 0, 0
};

const float gradY[] =
{
	1, 1,-1,-1,
	0, 0, 0, 0,
	1,-1, 1,-1
};

const float gradZ[] =
{
	0, 0, 0, 0,
	1, 1,-1,-1,
	1, 1,-1,-1
};

const tfxWideArray tfxF3_4 = tfxWideSetConst(1.0f / 3.0f);
const tfxWideArray tfxF2_4 = tfxWideSetConst(.366025403f);
const tfxWideArray tfxG2_4 = tfxWideSetConst(0.211324865f);
const tfxWideArray tfxG2_4x2 = tfxWideSetConst(0.42264973f);
const tfxWideArray tfxG3_4 = tfxWideSetConst(1.0f / 6.0f);
const tfxWideArray tfxG32_4 = tfxWideSetConst((1.0f / 6.0f) * 2.f);
const tfxWideArray tfxG33_4 = tfxWideSetConst((1.0f / 6.0f) * 3.f);
const tfxWideArrayi tfxONE = tfxWideSetConst(1);
const tfxWideArray tfxONEF = tfxWideSetConst(1.f);
const tfxWideArray tfxZERO = tfxWideSetConst(0.f);
const tfxWideArray tfxTHIRTYTWO = tfxWideSetConst(32.f);
const tfxWideArrayi tfxFF = tfxWideSetConst(255);
const tfxWideArray tfxPSIX = tfxWideSetConst(0.6f);

static const float tfxGRADIENTS_3D[] =
{
	0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
	1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
	1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
	0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
	1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
	1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
	0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
	1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
	1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
	0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
	1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
	1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
	0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
	1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
	1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
	1, 1, 0, 0,  0,-1, 1, 0, -1, 1, 0, 0,  0,-1,-1, 0
};

/**
* Permutation table. This is just a random jumble of all numbers 0-255.
*
* This produce a repeatable pattern of 256, but Ken Perlin stated
* that it is not a problem for graphic texture as the noise features disappear
* at a distance far enough to be able to see a repeatable pattern of 256.
*
* This needs to be exactly the same for all instances on all platforms,
* so it's easiest to just keep it as static explicit data.
* This also removes the need for any initialisation of this class.
*
* Note that making this an uint32_t[] instead of a uint8_t[] might make the
* code run faster on platforms with a high penalty for unaligned single
* byte addressing. Intel x86 is generally single-byte-friendly, but
* some other CPUs are faster with 4-aligned reads.
* However, a char[] is smaller, which avoids cache trashing, and that
* is probably the most important aspect on most architectures.
* This array is accessed a *lot* by the noise functions.
* A vector-valued noise over 3D accesses it 96 times, and a
* float-valued 4D noise 64 times. We want this to fit in the cache!
*/
const uint8_t tfx_permutation_table[] =
{
	151,160,137,91,90,15,
	131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
	190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
	88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
	77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
	102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
	135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
	5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
	223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
	129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
	251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
	49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
	138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
	151,160,137,91,90,15,
	131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
	190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
	88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
	77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
	102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
	135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
	5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
	223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
	129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
	251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
	49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
	138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
};

static const uint8_t tfx_perm_mod12[] =
{
	7, 4, 5, 7, 6, 3, 11, 1, 9, 11, 0, 5, 2, 5, 7, 9, 8, 0, 7, 6, 9, 10, 8, 3,
	1, 0, 9, 10, 11, 10, 6, 4, 7, 0, 6, 3, 0, 2, 5, 2, 10, 0, 3, 11, 9, 11, 11,
	8, 9, 9, 9, 4, 9, 5, 8, 3, 6, 8, 5, 4, 3, 0, 8, 7, 2, 9, 11, 2, 7, 0, 3, 10,
	5, 2, 2, 3, 11, 3, 1, 2, 0, 7, 1, 2, 4, 9, 8, 5, 7, 10, 5, 4, 4, 6, 11, 6,
	5, 1, 3, 5, 1, 0, 8, 1, 5, 4, 0, 7, 4, 5, 6, 1, 8, 4, 3, 10, 8, 8, 3, 2, 8,
	4, 1, 6, 5, 6, 3, 4, 4, 1, 10, 10, 4, 3, 5, 10, 2, 3, 10, 6, 3, 10, 1, 8, 3,
	2, 11, 11, 11, 4, 10, 5, 2, 9, 4, 6, 7, 3, 2, 9, 11, 8, 8, 2, 8, 10, 7, 10, 5,
	9, 5, 11, 11, 7, 4, 9, 9, 10, 3, 1, 7, 2, 0, 2, 7, 5, 8, 4, 10, 5, 4, 8, 2, 6,
	1, 0, 11, 10, 2, 1, 10, 6, 0, 0, 11, 11, 6, 1, 9, 3, 1, 7, 9, 2, 11, 11, 1, 0,
	10, 7, 1, 7, 10, 1, 4, 0, 0, 8, 7, 1, 2, 9, 7, 4, 6, 2, 6, 8, 1, 9, 6, 6, 7, 5,
	0, 0, 3, 9, 8, 3, 6, 6, 11, 1, 0, 0,
	7, 4, 5, 7, 6, 3, 11, 1, 9, 11, 0, 5, 2, 5, 7, 9, 8, 0, 7, 6, 9, 10, 8, 3,
	1, 0, 9, 10, 11, 10, 6, 4, 7, 0, 6, 3, 0, 2, 5, 2, 10, 0, 3, 11, 9, 11, 11,
	8, 9, 9, 9, 4, 9, 5, 8, 3, 6, 8, 5, 4, 3, 0, 8, 7, 2, 9, 11, 2, 7, 0, 3, 10,
	5, 2, 2, 3, 11, 3, 1, 2, 0, 7, 1, 2, 4, 9, 8, 5, 7, 10, 5, 4, 4, 6, 11, 6,
	5, 1, 3, 5, 1, 0, 8, 1, 5, 4, 0, 7, 4, 5, 6, 1, 8, 4, 3, 10, 8, 8, 3, 2, 8,
	4, 1, 6, 5, 6, 3, 4, 4, 1, 10, 10, 4, 3, 5, 10, 2, 3, 10, 6, 3, 10, 1, 8, 3,
	2, 11, 11, 11, 4, 10, 5, 2, 9, 4, 6, 7, 3, 2, 9, 11, 8, 8, 2, 8, 10, 7, 10, 5,
	9, 5, 11, 11, 7, 4, 9, 9, 10, 3, 1, 7, 2, 0, 2, 7, 5, 8, 4, 10, 5, 4, 8, 2, 6,
	1, 0, 11, 10, 2, 1, 10, 6, 0, 0, 11, 11, 6, 1, 9, 3, 1, 7, 9, 2, 11, 11, 1, 0,
	10, 7, 1, 7, 10, 1, 4, 0, 0, 8, 7, 1, 2, 9, 7, 4, 6, 2, 6, 8, 1, 9, 6, 6, 7, 5,
	0, 0, 3, 9, 8, 3, 6, 6, 11, 1, 0, 0
};

// 4 noise samples using simd
tfxWideFloat tfx__simd_noise_2d(const tfxWideFloat x4, const tfxWideFloat y4);
tfxWideFloat tfx__simd_noise_3d(const tfxWideFloat x4, const tfxWideFloat y4, const tfxWideFloat z4);

//-----------------------------------------------------------
//Section: Profiling
//-----------------------------------------------------------


typedef struct tfx_profile_stats_s {
	tfxU64 cycle_high;
	tfxU64 cycle_low;
	tfxU64 cycle_average;
	tfxU64 time_high;
	tfxU64 time_low;
	tfxU64 time_average;
	tfxU32 hit_count;
} tfx_profile_stats_t;

typedef struct tfx_profile_snapshot_s {
	tfxU32 hit_count;
	tfxU64 run_time;
	tfxU64 cycle_count;
}tfx_profile_snapshot_t;

typedef struct tfx_profile_s {
	const char *name;
	tfx_profile_snapshot_t snapshots[tfxPROFILER_SAMPLES];
} tfx_profile_t;

extern const tfxU32 tfxPROFILE_COUNT;
extern tfxU32 tfxCurrentSnapshot;
extern tfx_profile_t tfxProfileArray[];

#ifdef __cplusplus

struct tfx_profile_tag_t {
	tfx_profile_t *profile;
	tfx_profile_snapshot_t *snapshot;
	tfxU64 start_cycles;
	tfxU64 start_time;

	tfx_profile_tag_t(tfxU32 id, const char *name);

	~tfx_profile_tag_t() {
		tfx_AtomicAdd64(&snapshot->run_time, tfx_Microsecs() - start_time);
		tfx_AtomicAdd64(&snapshot->cycle_count, (tfx__rdtsc() - start_cycles));
	}

};

#else

typedef struct tfx_profile_tag_s {
	tfx_profile_t *profile;
	tfx_profile_snapshot_t *snapshot;
	tfxU64 start_cycles;
	tfxU64 start_time;
} tfx_profile_tag_t;

#endif

#ifdef tfxENABLE_PROFILING 
#define tfxPROFILE tfx_profile_tag_t tfx_tag = {(tfxU32)__COUNTER__, __FUNCTION__};
#else
#define tfxPROFILE __COUNTER__
#endif

//-----------------------------------------------------------
//Section: File_IO
//-----------------------------------------------------------

const tfxU32 tfxMAGIC_NUMBER = 559433300;                //'!XFT'
const tfxU32 tfxMAGIC_NUMBER_INVENTORY = 559304265;      //'!VNI'
const tfxU32 tfxFILE_VERSION = 2;

typedef struct tfx_package_entry_info_t {
	tfx_str512_t file_name;                     //The name of the file stored in the package
	tfxU64 offset_from_start_of_file;           //Offset from the start of the file to where the file is located
	tfxU64 file_size;                           //The size of the file
	tfx_stream_t data;                          //The file data
} tfx_package_entry_info_t;

#ifdef __cplusplus

typedef struct tfx_package_inventory_t {
	tfxU32 magic_number;                        //Magic number to confirm format of the Inventory
	tfxU32 entry_count;                         //Number of files in the inventory
	tfx_storage_map_t<tfx_package_entry_info_t> entries;        //The inventory list

	tfx_package_inventory_t() :
		magic_number(0),
		entry_count(0)
	{}
} tfx_package_inventory_t;

#else
typedef struct tfx_package_inventory_t {
	tfxU32 magic_number;                        //Magic number to confirm format of the Inventory
	tfxU32 entry_count;                         //Number of files in the inventory
	tfx_storage_map_t entries;                  //The inventory list
} tfx_package_inventory_t;
#endif

//Basic package manager used for reading/writing effects files
typedef struct tfx_package_header_s {
	tfxU32 magic_number;                        //Magic number to confirm file format
	tfxU32 file_version;                        //The version of the file
	tfxU32 flags;                               //Any state_flags for the file
	tfxU32 reserved0;                           //Reserved for future if needed
	tfxU64 offset_to_inventory;                 //Memory offset for the inventory of files
	tfxU64 user_data1;                          //Any data you might find useful
	tfxU64 user_data2;                          //Any data you might find useful
	tfxU64 reserved3;                           //More reserved space
	tfxU64 reserved4;                           //More reserved space
	tfxU64 reserved5;                           //More reserved space
} tfx_package_header_t;

typedef struct tfx_package_s {
	tfxU32 magic;
	tfx_stream_t file_path;
	tfx_package_header_t header;
	tfx_package_inventory_t inventory;
	tfxU64 file_size;                            //The total file size of the package, should match file size on disk
	tfx_stream file_data;                        //Dump of the data from the package file on disk
	tfxPackageFlags flags;
} tfx_package_t;

//------------------------------------------------------------
//Section: Struct_Types
//------------------------------------------------------------

#ifdef __cplusplus

extern tfx_vector_t<tfx_vec3_t> tfxIcospherePoints[6];

#endif

typedef struct tfx_attribute_node_s {
	float frame;
	float value;

	tfx_vec2_t left;
	tfx_vec2_t right;

	tfxAttributeNodeFlags flags;
	tfxU32 index;
	tfxU32 uid;

#ifdef __cplusplus
	tfx_attribute_node_s() : frame(0.f), value(0.f), flags(0), index(0) { }
	inline bool operator==(const struct tfx_attribute_node_s &n) { return n.frame == frame && n.value == value; }
#endif
}tfx_attribute_node_t;

typedef struct tfx_depth_index_s {
	tfxParticleID particle_id;
	float depth;
}tfx_depth_index_t;

//Used when a effect manager is grouping instances by effect. This way effects can be individually ordered and drawn/not drawn in order however you need
typedef struct tfx_effect_instance_data_s {
#ifdef __cplusplus
	tfx_vector_t<tfx_depth_index_t> depth_indexes[tfxLAYERS][2];
#else
	tfx_vector_t depth_indexes[tfxLAYERS][2];
#endif
	tfxU32 sprite_index_point[tfxLAYERS];
	tfxU32 cumulative_index_point[tfxLAYERS];
	tfxU32 depth_starting_index[tfxLAYERS];
	tfxU32 current_depth_buffer_index[tfxLAYERS];
	tfxU32 instance_start_index;
	tfxU32 instance_count;
	tfxU32 ribbon_index_point;
} tfx_effect_instance_data_t;

typedef struct tfx_face_s {
	int v[3];
} tfx_face_t;

typedef struct tfx_random_s {
	tfxU64 seeds[2];
}tfx_random_t;

typedef struct tfx_color_ramp_s {
	//These vectors are for sinusoidal color ramp generation (currently unused)
	tfx_vec3_t brightness;
	tfx_vec3_t contrast;
	tfx_vec3_t frequency;
	tfx_vec3_t offsets;
	tfxColorRampFlags flags;
	tfx_color_interpolation_mode interpolation_mode;
	tfx_rgba8_t colors[tfxCOLOR_RAMP_WIDTH];
}tfx_color_ramp_t;

typedef struct tfx_color_ramp_hash_s {
	tfxKey hash;		//A hash of the color ramp
	tfxU32 index;		//The index to lookup the ramp in the bitmap array
}tfx_color_ramp_hash_t;

typedef struct tfx_bitmap_s {
	int width;
	int height;
	int channels;
	int stride;
	tfx_size size;
	tfx_byte *data;
} tfx_bitmap_t;

typedef struct tfx_oscillator_s {
	float frequency;
	float amplitude;
	float offset_x;
	float offset_y;
	tfx_oscillator_type type;
} tfx_oscillator_t;

typedef struct tfx_oscillator_wide_s {
	tfxWideFloat frequency;
	tfxWideFloat amplitude;
	tfxWideFloat offset_x;
	tfxWideFloat offset_y;
} tfx_oscillator_wide_t;

typedef struct tfx_graph_wide_s {
	tfxWideFloat from;
	tfxWideFloat to;
	tfxWideFloat curve1;
	tfxWideFloat curve2;
} tfx_graph_wide_t;

typedef struct tfx_graph_id_s {
	tfx_graph_category category;
	tfx_graph_type type;
	tfxU32 index;
	tfxU32 graph_id;
	tfxU32 node_id;
	tfxKey path_hash;
} tfx_graph_id_t;

typedef struct tfx_graph_s {
	tfx_oscillator_wide_t wide_oscillator;
	tfx_graph_wide_t wide_graph;
	tfx_graph_preset graph_preset;
	tfx_graph_type type;
	tfx_graph_easing_type easing_type;
	tfx_effect_descriptor effector;
	tfxU32 uid_counter;
#ifdef __cplusplus
	tfx_bucket_array_t<tfx_attribute_node_t> nodes;
#else
	tfx_bucket_array_t nodes;
#endif
	tfxU32 index;
	float gamma;
	tfx_oscillator_t oscillator;
	tfxGraphFlags flags;
} tfx_graph_t TFX_ALIGN_AFFIX(16);

typedef struct tfx_graph_list_s {
#ifdef __cplusplus
	tfx_vector_t<tfx_graph_t> graphs;
#else
	tfx_vector_t graphs;
#endif
	tfx_effect_descriptor_type effect_descriptor_type;
	tfx_color_ramp_t color_ramps;
	tfx_index color_ramp_bitmap_indexes;
} tfx_graph_list_t;

tfxAPI_EDITOR void tfx__init_graph(tfx_graph_t *graph, tfxU32 node_bucket_size);

typedef struct tfx_path_nodes_soa_s {
	float *x;
	float *y;
	float *z;
	float *length;
} tfx_path_nodes_soa_t;

typedef struct tfx_path_quaternion_s {
	tfxU64 quaternion;
	float grid_coord;
	float age;
	tfxU32 cycles;
} tfx_path_quaternion_t;

typedef struct tfx_path_buffers_s {
#ifdef __cplusplus
	tfx_vector_t<tfx_vec3_t> nodes;
#else
	tfx_vector_t nodes;
#endif
	tfx_soa_buffer_t node_buffer;
	tfx_path_nodes_soa_t node_soa;
} tfx_path_buffers_t;

typedef struct tfx_path_settings_s {
	tfxKey key;
	tfx_str32_t name;
	int node_count;
	int nodes_to_commit;
	tfxEmitterPathFlags flags;
	float rotation_range;
	float rotation_pitch;      
	float rotation_yaw;
	tfxU32 maximum_active_paths;
	tfxU32 maximum_paths;
	float rotation_cycle_length;
	float rotation_stagger;
	tfx_vec3_t offset;
	tfx_vec3_t builder_parameters;
	tfx_path_extrusion_type extrusion_type;
} tfx_path_settings_t;

typedef struct tfx_emitter_path_s {
	tfx_path_settings_t settings;
	tfx_path_buffers_t buffers;
} tfx_emitter_path_t;

typedef struct tfx_base_s {
	tfx_vec2_t size;
	float velocity;
	float spin;
	float weight;
}tfx_base_t;

typedef struct tfx_camera_settings_s {
	tfx_vec3_t camera_position;
	float camera_pitch;
	float camera_yaw;
	float camera_fov;
	float camera_floor_height;
	float camera_isometric_scale;
	bool camera_isometric;
	bool camera_hide_floor;
} tfx_camera_settings_t;

//Only in the editor?
typedef struct tfx_preview_camera_settings_s {
	tfx_camera_settings_t camera_settings;
	tfx_render_view_mode view_mode;
	float effect_z_offset;
	float camera_speed;
	bool attach_effect_to_camera;
} tfx_preview_camera_settings_t;

//todo: this probably only needs to be in the editor, no use for it in the library? Maybe in the future as an alternative way to play back effects...
typedef struct tfx_sprite_sheet_settings_s {
	tfx_vec3_t position;
	tfx_vec2_t frame_size;
	float scale;
	float zoom;
	int frames;
	int current_frame;
	int frame_offset;
	int extra_frames_count;
	tfxU32 seed;
	tfxAnimationFlags animation_flags;
	tfx_render_view_mode view_mode;
	tfxU32 needs_exporting;
	float max_radius;
	tfxU32 largest_frame;
	float playback_speed;
	float effect_z_offset;
	tfx_export_color_options color_option;
	tfx_export_options export_option;
	tfx_camera_settings_t camera_settings;
} tfx_sprite_sheet_settings_t;

//This struct has the settings for recording sprite data frames so that they can be played back as an alternative to dynamic particle updating
typedef struct tfx_sprite_data_settings_s {
	int real_frames;
	int frames_after_compression;
	int current_frame;
	float current_time;
	float animation_length_in_time;
	int frame_offset;
	int extra_frames_count;
	tfxU32 seed;
	tfxAnimationFlags animation_flags;
	tfx_render_view_mode view_mode;
	tfxU32 needs_exporting;
	float max_radius;
	tfxU32 largest_frame;
	float playback_speed;
	float recording_frame_rate;
} tfx_sprite_data_settings_t;

//------------------------------------------------------------

//API structs you can access in various ways to update and render effects in realtime

//Image data for particle shapes. This is passed into your custom ShapeLoader function for loading image textures into whatever renderer you're using
typedef struct tfx_image_data_s {
	//This can be a ptr to the image texture for rendering. You must assign this in your ShapeLoader function
	void *ptr;
	//Index of the image, deprecated, image hash should be used now instead.
	tfxU32 shape_index;
	//Name of the image
	tfx_str256_t name;
	//A hash of the image data for a unique id and which can also be used to see if an image has already been loaded
	tfxU64 image_hash;
	//The size of one frame of the image
	tfx_vec2_t image_size;
	//Image index refers to any index that helps you look up the correct image to use. this could be an index in a texture array for example.
	tfxU32 image_index;
	//The number of frames in the image, can be one or more
	float animation_frames;
	//Maximum distance to the nearest transparent edge of the image
	float max_radius;
	int import_filter;
	tfxU32 compute_shape_index;
	//use this definition if you need more spefic data to point to the image texture in whatever renderer you're using
	//Just define tfxCUSTOM_IMAGE_DATA before you include timelinefx.h
#ifdef tfxCUSTOM_IMAGE_DATA
	tfxCUSTOM_IMAGE_DATA
#endif // tfxCUSTOM_IMAGE_DATA
} tfx_image_data_t;

typedef struct tfx_particle_emitter_properties_s {
	//Angle added to the rotation of the particle when spawned or random angle range if angle setting is set to tfx_random_t
	tfx_vec3_t angle_offsets;
	//When aligning the billboard along a vector, you can set the type of vector that it aligns with
	tfx_vector_align_type vector_align_type;
	//For other emitter emission types, this hash is the location of the other emitter so that it can be used to connect the two
	//The type of billboarding: 0 = use billboarding (always face camera), 1 = No billboarding, 2 = No billboarding and align with motion
	tfx_billboarding_option billboard_option;

	//The rotation of particles when they spawn, or behave overtime if tfxAlign is used
	tfxAngleSettingFlags angle_settings;
	//Should particles emit towards the center of the emitter or away, or in a specific direction
	tfx_emission_direction emission_direction;

	//The type of noise algorithm to use
	tfx_noise_type noise_algorithm;

	//How particles should behave when they reach the end of the line
	tfx_line_traversal_end_behaviour end_behaviour;
	//Bit field of various boolean state_flags
	tfxParticleControlFlags compute_flags;
	//Offset to draw particles at
	tfx_vec2_t image_handle;
	//image handle packed into 16bit floats
	tfxU32 image_handle_packed;
	//This is only used for the animation manager when sprite data is added to the animation manager. This is used to map
	//the property_index to the animation property index so the sprite data can point to a new index where some emitter properties
	//are stored on the GPU for looking up from the sprite data
	tfxU32 animation_property_index;
} tfx_particle_emitter_properties_t;

typedef struct tfx_shared_emitter_properties_s {
	//Animation frame rate
	float frame_rate;
	//The final frame index of the animation
	float end_frame;
	//The start frame index of the animation
	float start_frame;
	//Pointer to the ImageData in the EffectLibary. 
	tfx_image_data_t *image;
	//Point, area, ellipse emitter etc.
	tfx_emission_type emission_type;
	//The number of rows/columns/ellipse/line points in the grid when spawn on grid flag is used
	tfx_vec3_t grid_points;
	//Can this be removed if we're using the image hash now?
	tfxU32 image_index;
	//The shape being used for all particles spawned from the emitter
	tfxKey image_hash;
	//When single flag is set, spawn this amount of particles in one go
	tfxU32 spawn_amount;
	//When single flag is set, spawn this variable amount of particles in one go
	tfxU32 spawn_amount_variation;
	//If single shot flag is set then you can limit how many times it will loop over it's overtime graphs before expiring
	tfxU32 single_shot_limit;
	//Milliseconds to delay spawing
	float delay_spawning;
	//When relative position is set you can create a lag between the particle position and the emitter position base on the particle age
	float relative_lag;
	//When the emission type is shared emitter then this is the hash of the shared emitter.
	tfxKey paired_emitter_hash;
	//Layer of the effect manager that the particle is added to
	tfxU32 layer;
} tfx_shared_properties_t;

typedef struct tfx_ribbon_bucket_info_s {
	tfxU32 segment_count;
	tfxRibbonBucketComputeShaderType shader_type;
} tfx_ribbon_bucket_info_t;

typedef struct tfx_ribbon_emitter_properties_s {
	tfx_ribbon_bucket_info_t bucket_info;
	tfxKey ribbon_bucket_id;
	tfx_vec3_t fixed_angle_normal;
	tfxU32 animation_property_index;
} tfx_ribbon_emitter_properties_t;

//Stores the most recent parent effect (with global attributes) spawn control values to be applied to sub emitters.
typedef struct tfx_parent_spawn_controls_s {
	float life;
	float size_x;
	float size_y;
	float velocity;
	float spin;
	float pitch_spin;
	float yaw_spin;
	float intensity;
	float splatter;
	float weight;
}tfx_parent_spawn_controls_t;

typedef struct tfx_path_state_s {
	tfx_path_quaternion_t *path_quaternions;
	tfxU32 path_quaternion_index;
	tfxU32 last_path_index;
	float path_stagger_counter;
	tfxU32 path_cycle_count;
	tfxU32 active_paths;
	tfxU32 path_start_index;
	tfxU32 ribbon_index;
} tfx_path_state_t;

//This is a struct that stores an emitter state that is currently active in a effect manager.
//Todo: maybe split this up into static variables that stay the same (they're just properties copied from the emitter in the library
//      and dynamic variables that change each frame.
typedef struct tfx_particle_emitter_state_s {
	//State data
	float age;										//SoA?
	float highest_particle_age;
	float delay_spawning;
	float timeout_counter;
	float timeout;
	float amount_remainder;
	float spawn_quantity;
	float qty_step_size;
	tfx_vec3_t handle;
	tfx_bounding_box_t bounding_box;
	//Position, scale and rotation values
	tfx_vec3_t local_position;						//SoA?
	tfx_vec3_t world_position;						//SoA?
	tfx_vec3_t captured_position;					//SoA?
	tfx_vec3_t world_rotations;						//SoA?
	tfx_quaternion_t captured_rotation;				//SoA?
	tfx_quaternion_t rotation;						//SoA?

	//Static data that won't change frame by frame
	float loop_length;
	float max_life;
	float oscillator_time;
	tfxU64 image_handle_packed;

	//Indexes
	tfxIndex graph_list_index;
	tfxIndex transform_index;
	tfxIndex path_attributes;
	tfxIndex root_index;
	tfxIndex parent_index;
	tfxIndex properties_index;
	tfxIndex shared_index;
	tfxIndex info_index;
	tfxIndex sprites_index;
	tfxIndex seed_index;
	tfxIndex spawn_locations_index;    //For other_emitter emission type and storing the last known position of the particle
	tfxIndex other_emitter_index;      //For other_emitter emission type, this in the index of the other emitter in the effect manager
	tfxIndex particles_index;
	tfxIndex gpu_group_index;          //Index into pm->gpu_groups; tfxINVALID if not on the GPU path

	tfxU32 sprites_count;

	tfxParticleEmitterFlags property_flags;
	tfxSharedEmitterFlags shared_flags;
	tfx_effect_descriptor source_emitter;
	tfx_library library;

	//Control Data (May change frame by frame
	tfxKey ribbon_bucket_id;		 //For spawn on ribbon emission type and storing the last known position of the particle
	float image_frame_rate;
	float end_frame;
	tfx_vec3_t grid_coords;
	tfx_vec3_t grid_direction;
	tfx_vec3_t emitter_size;
	float emission_alternator;
	tfxEmitterStateFlags state_flags;
	tfxEmitterControlProfileFlags control_profile;
	tfx_vec2_t image_size;
	tfx_vec3_t angle_offsets;
	tfx_path_state_t path_state;
} tfx_particle_emitter_state_t TFX_ALIGN_AFFIX(16);

//This is a struct that stores an effect state that is currently active in a effect manager.
typedef struct tfx_effect_state_s {
	//State data that can change every frame
	tfx_quaternion_t rotation;
	float age;
	float highest_particle_age;
	float timeout_counter;
	float timeout;
	tfx_vec3_t handle;
	tfxParticleEmitterFlags property_flags;
	tfxEffectPropertyFlags effect_flags;
	tfxEmitterStateFlags state_flags;
	float loop_length;
	float oscillator_time;
	//Position, scale and rotation values
	tfx_vec3_t translation;
	tfx_vec3_t local_position;
	tfx_vec3_t world_position;
	tfx_vec3_t captured_position;
	tfx_vec3_t local_rotations;
	tfx_vec3_t world_rotations;
	tfx_bounding_box_t bounding_box;

	//Static data
	tfxIndex graph_list_index;
	tfxIndex transform_index;

	tfx_library library;
	tfx_effect_descriptor source_effect;
	//User Data
	void *user_data;
	void(*update_callback)(tfx_effect_manager pm, tfxEffectID effect_index);

	//Spawn controls
	tfx_parent_spawn_controls_t spawn_controls;
	tfx_vec3_t emitter_size;
	float stretch;
	float noise;
	float overal_scale;
	float noise_base_offset;
	tfxU32 sort_passes;

	//When organising instance_data per effect this is the index to the sprite buffers containing all the effects.
	tfx_effect_instance_data_t instance_data;

	//The emitters within this effect.
#ifdef __cplusplus
	tfx_vector_t<tfxU32> emitter_indexes[2];
#else
	tfx_vector_t emitter_indexes[2];
#endif
	tfxU32 emitter_start_size;

} tfx_effect_state_t TFX_ALIGN_AFFIX(16);

typedef struct tfx_ribbon_s {	//64 bytes (56 bytes data + 8 padding for std430 alignment)
	tfx_vec4_t position;
    float width;
	tfxU32 start_index;
	tfxU32 flags;
	tfxU32 captured_index;				
	tfxU64 quaternion;
	tfxU32 emitter_index;
	tfxU32 texture_indexes;
	tfxU32 intensity_gradient_map;			//Multiplier for the color of the ribbon
	tfxU32 curved_alpha;					//Sharpness and dissolve amount value for fading the image
	tfxU32 _padding[2];						//Padding to 64 bytes for std430 alignment with vec4
} tfx_ribbon_t;

typedef struct tfx_ribbon_soa_s {
	tfx_ribbon_t *ribbon_instances;
	float *age;
	float *max_age;
	float *random_age;
	float *image_frame;
	tfxU32 *path_index;
	float *grid_index;
	tfxU32 *single_loop_count;
	tfxU32 *uid;
} tfx_ribbon_soa_t;

typedef struct tfx_gpu_emitter_s {
	tfx_vec4_t quaternion;
	tfx_vec3_t position;
	tfxU32 lookup_offset;
	tfx_vec3_t captured_position;
	tfxU32 padding1;
	tfx_vec3_t scale;
	tfxU32 padding2;
	tfx_vec3_t fixed_angle_normal;
	tfxU32 padding3;
} tfx_gpu_emitter_t;

//---- GPU compute particle buffer management ----
// Uncomment to enable Phase 3 CPU-side shadow buffer that mirrors GPU ring buffer layout for validation.
// Compiles away to nothing when not defined.
// #define tfxGPU_VALIDATION

//Per-frame spawn record used by the group tracking ring to drive deterministic head bumping.
//One entry is sealed per update tick; the ring is sized to hold life_ceiling_ms worth of ticks.
typedef struct tfx_gpu_spawn_tracking_s {
	float  spawn_time_ms;    //Absolute time (ms) when this batch of particles was spawned
	tfxU32 particle_count;  //Total particles spawned into this group during that tick
} tfx_gpu_spawn_tracking_t;

//One GPU particle ring buffer shared by all emitters with the same control profile and life bucket.
//Particles from different emitters are interleaved; each particle carries its emitter param index
//so the compute shader can look up the correct graphs for that particle.
typedef struct tfx_gpu_particle_group_s {
	tfxU32  property_index;			//The property index of the emitter to identify it
	float   life_ceiling_ms;		//All particles guaranteed expired after this many ms from spawn
	tfxU32  ring_head;				//Oldest live particle position in the group ring
	tfxU32  ring_tail;				//Next write position in the group ring
	tfxU32  current_size;			//The current number of particles in the buffer.
	tfxU32  ring_capacity;			//Total particle slots; grows as emitters join during initialisation
	tfxU64  gpu_buffer_size_bytes;	//Set by renderer when the GPU buffer is allocated
	//Tracking ring: one entry per update tick; capacity = ceil(life_ceiling_ms / 8ms) + 2
#ifdef __cplusplus
	tfx_vector_t<tfx_gpu_spawn_tracking_t> tracking;
#else
	tfx_vector_t tracking;
#endif
	tfxU32  tracking_head;			//Index of the oldest valid tracking entry
	tfxU32  tracking_count;			//Number of valid entries currently in the tracking ring
	tfxU32  tracking_capacity;		//Fixed at creation; sized for minimum 8ms tick (120fps)
	tfxU32  frame_spawn_count;		//Particles spawned into this group this frame; sealed into tracking at tick
	tfxU32  active_emitter_count;	//Number of emitters currently assigned to this group
} tfx_gpu_particle_group_t;
//---- end GPU compute particle buffer management ----

//---- GPU Validation Shadow Mode (Phase 3) ----
//All particle fields that the GPU ring buffer will hold, one uint32 slot each.
typedef enum tfx_gpu_particle_field_e {
	tfx_gpu_field_position_x = 0,
	tfx_gpu_field_position_y,
	tfx_gpu_field_position_z,
	tfx_gpu_field_age,
	tfx_gpu_field_max_age,
	tfx_gpu_field_life,
	tfx_gpu_field_velocity_normal,          //packed uint (stored as float-sized slot)
	tfx_gpu_field_base_velocity,
	tfx_gpu_field_base_weight,
	tfx_gpu_field_base_size_x,
	tfx_gpu_field_base_size_y,
	tfx_gpu_field_base_roll_spin,
	tfx_gpu_field_intensity_factor,
	tfx_gpu_field_random_color,
	tfx_gpu_field_image_frame,
	tfx_gpu_field_flags_single_loop_count,  //uint, stored as float-sized slot
	tfx_gpu_field_noise_offset,
	tfx_gpu_field_noise_resolution,
	tfx_gpu_field_count
} tfx_gpu_particle_field_t;

//Index into the shadow_buffer uint32 array.
//Shadow buffer layout: [field_0: capacity uint32s][field_1: capacity uint32s]...[field_N-1: capacity uint32s]
#define tfxGPU_SHADOW_OFFSET(field, pos, capacity) ((field) * (capacity) + (pos))
//---- end GPU Validation Shadow Mode ----

typedef struct tfx_ribbon_emitter_state_s {
	//State data
	float frame;
	float age;
	float amount_remainder;
	float spawn_quantity;
	float qty_step_size;
	float timeout_counter;
	float timeout;
	float delay_spawning;
	float max_life;
	float loop_length;
	float oscillator_time;
	tfx_vec3_t handle;
	tfxRibbonEmitterFlags ribbon_property_flags;
	tfxSharedEmitterFlags shared_flags;
	//Position, scale and rotation values
	tfx_vec3_t local_position;
	tfx_vec3_t world_position;
	tfx_vec3_t captured_position;
	tfx_vec3_t local_rotations;
	tfx_vec3_t world_rotations;
	tfx_quaternion_t rotation;

	tfx_bounding_box_t bounding_box;

	tfxIndex graph_list_index;
	tfxIndex transform_index;
	tfxIndex path_attributes;
	tfxIndex seed_index;
	tfxIndex parent_index;
	tfxIndex properties_index;
	tfxIndex shared_index;
	tfxIndex info_index;
	tfxIndex spawn_locations_index;					//For other_emitter emission type and storing the last known position of the particle

	tfx_path_state_t path_state;

	tfxRibbonEmitterStateFlags state_flags;

	tfxU32 segment_count;
	tfxU32 active_ribbons;
	tfx_effect_descriptor source_ribbon;
	tfx_library library;

#ifdef __cplusplus
	tfx_vector_t<tfxU32> ribbon_indexes[2];
#else
	tfx_vector_t ribbon_indexes[2];
#endif

	//Control Data
	tfxKey ribbon_bucket_id;
	tfxU32 static_segment_start_index;				//For static paths so that we only have to build the ribbon once for all instances of it.
	float image_frame_rate;
	float end_frame;
	tfx_vec3_t emitter_size;
	tfxU32 gpu_emitter_index;
} tfx_ribbon_emitter_state_t TFX_ALIGN_AFFIX(16);

//An tfx_effect_descriptor_t can either be an effect which stores effects and global graphs for affecting all the attributes in the emitters,
//an emitter which spawns all of the particles or a ribbon for spawning ribbon segments.
//This is only for library storage, when using to update each frame aspects of this are copied to tfx_effect_state_t, tfx_particle_emitter_state_t and tfx_ribbon_emitter_state_t for realtime updates
typedef struct tfx_effect_descriptor_s {
	tfxU32 magic;
	//Name of the effect
	tfx_str64_t name;
	//The path of the effect in the library
	tfx_str512_t path;
	//Every effect and emitter in the library gets a unique id
	tfxU32 uid;
	//The max_radius of the emitter, taking into account all the particles that have spawned and active (editor only)
	float max_radius;
	//Index to sprite sheet settings stored in the effect library. 
	tfxIndex sprite_sheet_settings_index;
	//Index to sprite data settings stored in the effect library. 
	tfxIndex sprite_data_settings_index;
	//Index to preview camera settings stored in the effect library. Would like to move this at some point
	tfxIndex preview_camera_settings;
	//The current state of the effect/emitter used in the editor only at this point
	tfxEmitterStateFlags state_flags;
	//Is this an tfxEffectType or tfxEmitterType
	tfx_effect_descriptor_type type;
	//The index within the library that this exists at
	tfxIndex library_index;
	//A hash of the directory path to the effect ie Flare/spark, and also a UID for the effect/emitter
	tfxKey path_hash;
	//Pointer to the immediate parent
	tfx_effect_descriptor parent;

	//All the below fields will be used by the effect/emitter states when added to a effect manager
	//Indexes into library storage
	tfxIndex property_index;		//this will be the index to either particle emitter, ribbon emitter properties. Effects don't have any properties that aren't shared.
	//Shared properties used by all emitter/ribbon types. Doesn't apply to effects
	tfxIndex shared_index;
	//The number of millisecs before an effect or emitter will loop back round to the beginning of it's graph lookups
	float loop_length;
	//All graphs that the effect uses to lookup attribute values are stored in the library. 
	//Effects, particle emitters and ribbon emitters get their own set of graphs 
	tfxIndex graph_list_index;
	//Index to the graph list storing all of the transform graphs.
	//Transform graphs are shared by all descriptor types (except folders)
	tfxIndex transform_index;
	//If the emitter uses a path for emission then the index to the path in the library is stored here.
	tfxIndex path_attributes;
	//Graph data is compiled and uploaded to the GPU (currently for ribbons only) this is the offset into the buffer where the data starts for this emitter
	tfxU32 gpu_lookup_offset;
	//The type of function that should be called to update particle positions
	tfxEmitterControlProfileFlags control_profile;
	//Property flags for emitters
	tfxParticleEmitterFlags property_flags;
	//Property flags for ribbon_emitters
	tfxRibbonEmitterFlags ribbon_flags;
	//Shared flags for all emitter types
	tfxSharedEmitterFlags shared_flags;
	//Flags specific to effects
	tfxEffectPropertyFlags effect_flags;
	//Offset of emitters and effects
	tfx_vec3_t emitter_handle;
	//When not using insert sort to guarantee particle order, sort passes offers a more relaxed way of ordering particles over a number of frames.
	//The more passes the more quickly ordered the particles will be but at a higher cost
	tfxU32 sort_passes;
	//Base noise offset random range so that noise patterns don't repeat so much over multiple effects
	float noise_base_offset_range;
	//Custom user data, can be accessed in callback functions
	void *user_data;
	void(*update_callback)(tfx_effect_manager pm, tfxEffectID effect_index);
	//The maximum amount of life that a particle can be spawned with taking into account base + variation life values
	float max_life;
	//A link to the library that this effect/emitter belongs to
	tfx_library library;

	//List of children if this is a folder or stage, or list of particle/ribbon emitters
#ifdef __cplusplus
	tfx_vector_t<tfx_effect_descriptor> children;
#else
	tfx_vector_t children;
#endif
} tfx_effect_descriptor_t;

typedef struct tfx_unique_sprite_id_s {
	tfxU32 uid;
	tfxU32 age;
	tfxU32 property_index;
}tfx_unique_sprite_id_t;

//These all point into a tfx_soa_buffer_t, initialised with InitParticleSoA. Max Current Bandwidth: 108 bytes in total. Or if half-floats are used: 90 bytes
//Unfortunately in order make the most use of half floats the minimum requirements become AVX + F16C for half float conversion in simd which is around ~95% according to steam survey.
//Note that not all of these are used, it will depend on the emitter and which attributes it uses. So to save memory,
//when the the buffer is initialised only the fields that are needed for the emitter will be used.
typedef struct tfx_particle_soa_s {
	tfxU32 *uid;
	tfxU32 *sprite_index;
	float *age;
	float *inv_max_age;
	float *position_x;
	float *position_y;
	float *position_z;
	union {
		tfxU32 *rotation_offsets;		//Packed into 10bit ints for each axis
		float *rotation_offset;			//Just use a float if the particle always faces the camera
	};
	tfxU32 *velocity_normal;			//Packed into 10bit ints for each axis
	tfxU64 *quaternion;					//Used for paths where the path can be rotated per particle based on the emission direction
	tfxU32 *depth_index;
	float *path_position;
	float *path_offset;
	tfxU32 *flags_single_loop_count;	//Packed flags and single loop count
	union {
		float *base_velocity;
		float *path_scale_variation;
	};
	float *base_weight;
	float *base_size_x;
	float *base_size_y;
	float *noise_offset;
	float *noise_resolution;
	float *base_roll_spin;
	float *base_pitch_spin;
	float *base_yaw_spin;
	float *intensity_factor;
	float *random_color;
	float *image_frame;
} tfx_particle_soa_t;

typedef struct tfx_spawn_points_soa_s {
	float *position_x;
	float *position_y;
	float *position_z;
	float *captured_position_x;
	float *captured_position_y;
	float *captured_position_z;
	float *age;
}tfx_spawn_points_soa_t;

typedef struct tfx_sprite_transform_s {
	tfx_vec3_t position;							//The position of the sprite, x, y - world, z, w = captured for interpolating
	tfx_vec3_t rotations;							//Rotations of the sprite
	tfx_vec2_t scale;								//Scale
} tfx_sprite_transform_t;

//When exporting effects as sprite data each frame gets frame meta containing information about the frame such as bounding box and sprite count/offset into the buffer
typedef struct tfx_frame_meta_s {
	tfxU32 index_offset[tfxLAYERS];					//All sprite data is contained in a single buffer and this is the offset to the first sprite in the range
	tfxU32 sprite_count[tfxLAYERS];					//The number of instance_data in the frame for each layer
	tfxU32 cumulative_offset[tfxLAYERS];			//The cumulative number of instance_data in the frame for each layer
	tfxU32 total_sprites;							//The total number of instance_data for all layers in the frame
	tfxU32 captured_offset;							//The amount to offset the captured index by
	tfx_vec3_t min_corner;							//Bounding box min corner
	tfx_vec3_t max_corner;							//Bounding box max corner. The bounding box can be used to decide if this frame needs to be drawn
	tfx_vec3_t bb_center_point;						//The center point of the bounding box. For the fastest checking against a viewing frustum, you can combine this with radius
	float radius;									//The radius of the bounding box
} tfx_frame_meta_t;

typedef struct tfx_ribbon_frame_meta_s {
	tfxU32 index_offset;
	tfxU32 captured_offset;
	tfxU32 ribbon_count;
} tfx_ribbon_frame_meta_t;

typedef struct tfx_instance_s {		//64 bytes (52 bytes data + 12 padding for std430 alignment)
	tfx_vec4_t position;							//The position of the billboard with stretch in w
	tfxU64 quaternion;								//Rotation of the billboard stored as a 16-bit snorm quaternion
	tfx_float16x4_t size_handle;					//Size of the sprite in pixels and the handle packed into a u64 (4 16bit floats)
	tfx_float8x4_t alignment;						//normalised alignment vector 3 8bit floats packed into 32 bits. Free byte here.
	tfx_float16x2_t intensity_gradient_map;			//Multiplier for the color and life of particle
	tfx_float8x4_t curved_alpha_life;				//Sharpness and dissolve amount value for fading the image plus the age of the particle value packed into 3 bit unorms. Free byte here.
	tfxU32 indexes;									//[color ramp y index, color ramp texture array index, capture flag, image data index (1 bit << 15), billboard alignment (2 bits << 13), image data index max 8191 images]
	tfxU32 captured_index;							//Index to the sprite in the buffer from the previous frame for interpolation
	tfxU32 _padding[3];								//Padding to 64 bytes for std430 alignment with vec4
} tfx_instance_t;

//These structs are for animation sprite data that you can upload to the gpu
typedef struct tfx_sprite_instance_data_s {    //56 bytes padding to 64
	tfx_vec4_t position_stretch;                    //The position of the sprite, x, y - world, z, w = captured for interpolating
	tfxU64 quaternion;								//Rotation of the billboard stored as a 16-bit snorm quaternion
	tfx_float16x4_t size_handle;					//Size of the sprite in pixels and the handle packed into a u64 (4 16bit floats)
	tfx_float8x4_t alignment;						//normalised alignment vector 3 floats packed into 8bits
	tfx_float16x2_t intensity_gradient_map;			//Multiplier for the color and life of particle
	tfx_float8x4_t curved_alpha_life;				//Sharpness and dissolve amount value for fading the image 2 16bit floats packed
	tfxU32 indexes;									//[color ramp y index, color ramp texture array index, capture flag, image data index (1 bit << 15), billboard alignment (2 bits << 13), image data index max 8191 images]
	tfxU32 captured_index;							//Index to the sprite in the buffer from the previous frame for interpolation
	tfxU32 additional;								//Padding, but also used to pack lerp offset and property index
	tfxU32 padding[2];
} tfx_sprite_instance_data_t;

typedef struct tfx_ribbon_instance_data_s {	//64 bytes, mirrors tfx_ribbon_t layout
	tfx_vec4_t position;						//xyz position, w = normalised age
	float width;
	tfxU32 start_index;							//Into shared segment data
	tfxU32 flags;
	tfxU32 captured_index;						//Previous frame ribbon for interpolation (was _padding_pre_quat)
	tfxU64 quaternion;							//16-bit snorm packed rotation
	tfxU32 emitter_index;
	tfxU32 texture_indexes;
	tfxU32 intensity_gradient_map;
	tfxU32 curved_alpha;
	tfxU32 additional;							//lerp_offset (low 16) + property_index (high 16)
	tfxU32 padding;
} tfx_ribbon_instance_data_t;

//Animation sprite data that is used on the cpu to bake the data
typedef struct tfx_sprite_data_soa_s {
	tfx_instance_t *billboard_instance;
	tfx_unique_sprite_id_t *uid;
	float *lerp_offset;
}tfx_sprite_data_soa_t;

typedef struct tfx_ribbon_data_soa_s {
	tfx_ribbon_t *ribbon_instance;
	tfx_unique_sprite_id_t *uid;
	float *lerp_offset;
}tfx_ribbon_data_soa_t;

typedef struct tfx_wide_lerp_transform_result_s {
	float position[3];
	float rotations[3];
	float scale[2];
}tfx_wide_lerp_transform_result_t;

typedef struct tfx_wide_lerp_data_result_s {
	float stretch;
	float intensity;
	float color[4];
	float padding[2];
}tfx_wide_lerp_data_result_t;

typedef struct tfx_sprite_data_metrics_s {
	tfx_str64_t name;
	tfxKey path_hash;
	tfxU32 start_offset;    //Only applies to animation manager
	tfxU32 frames_after_compression;
	tfxU32 real_frames;
	tfxU32 frame_count;
	float animation_length_in_time;
	tfxU32 total_sprites;
	tfxU32 total_memory_for_sprites;
#ifdef __cplusplus
	tfx_vector_t<tfx_frame_meta_t> frame_meta;
	tfx_vector_t<tfx_ribbon_frame_meta_t> ribbon_frame_meta;
#else
	tfx_vector_t frame_meta;
	tfx_vector_t ribbon_frame_meta;
#endif
	tfxAnimationManagerFlags flags;
	tfxAnimationFlags animation_flags;
	//Ribbon Specific metrics
	tfxU32 total_ribbons;
	tfxU32 total_memory_for_ribbons;
	tfxU32 ribbon_start_offset;
	tfxU32 ribbon_segment_start_offset;
} tfx_sprite_data_metrics_t;

typedef struct tfx_ribbon_segment_s {
	tfx_vec3_t position;
	tfx_float16x2_t intensity_gradient_map;			//Multiplier for the color of the ribbon
	tfx_float8x4_t curved_alpha;					//Sharpness and dissolve amount value for fading the image plus the age of the particle value packed into 3 bit unorms
	tfxU32 padding[3];
} tfx_ribbon_segment_t;

typedef struct tfx_sprite_data_s {
	float frame_compression;
	tfx_sprite_data_metrics_t normal;
	tfx_sprite_data_metrics_t compressed;
	tfx_soa_buffer_t real_time_sprites_buffer;
	tfx_sprite_data_soa_t real_time_sprites;
	tfx_soa_buffer_t compressed_sprites_buffer;
	tfx_sprite_data_soa_t compressed_sprites;
	//Ribbon specific data
	tfx_soa_buffer_t real_time_ribbons_buffer;
	tfx_ribbon_data_soa_t real_time_ribbons;
	tfx_soa_buffer_t compressed_ribbons_buffer;
	tfx_ribbon_data_soa_t compressed_ribbons;
#ifdef __cplusplus
	tfx_vector_t<tfx_ribbon_segment_t> shared_segments;
#else
	tfx_vector_t shared_segments;
#endif
	tfxU32 total_ribbons;
	bool has_ribbons;
} tfx_sprite_data_t;

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

typedef struct tfx_ribbon_data_push_s {
	tfxU32 animation_instances_total;
	tfxU32 ribbons_total;
	tfxU32 offset_index;
	tfxU32 animation_instances_index;
	tfxU32 ribbon_data_index;
	tfxU32 ribbon_output_index;
	tfxU32 ribbon_properties_index;
} tfx_ribbon_data_push_t;

//This can be sent as a push constant to the gpu
typedef struct tfx_ribbon_bucket_globals_s  {
	tfx_vec4_t camera_position;
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
} tfx_ribbon_bucket_globals_t;

typedef struct tfx_ribbon_segment_soa_s {
	float *x;
	float *y;
	float *z;
	float *width;
} tfx_ribbon_segment_soa_t;

typedef struct tfx_ribbon_vertex_s {
	tfx_vec3_t position;
	tfxU32 segment_index;
	tfx_vec2_t uv_offset_scale;
	tfxU32 ribbon_index;
	tfxU32 clipped;
} tfx_ribbon_vertex_t;

typedef struct tfx_ribbon_buffer_info_s {
	tfxU32 vertices_per_segment; 
	tfxU32 triangles_per_segment; 
	tfxU32 indices_per_segment;  
	tfxU32 total_segments; 
	tfxU32 index_count; 
	tfxKey pipeline_index;
} tfx_ribbon_buffer_info_t;

typedef struct tfx_ribbon_bucket_s {
	tfx_ribbon_bucket_globals_t globals;
	tfx_ribbon_buffer_info_t buffer_info;
	tfxU32 ribbon_index_offset;
	tfxU32 active_ribbons;
	tfxU32 highest_ribbon_index;
	tfxU32 lowest_ribbon_index;
	tfxU32 highest_segment_index;
	tfxU32 lowest_segment_index;
	tfx_soa_buffer_t ribbons_buffer;
	tfx_ribbon_soa_t ribbons;
#ifdef __cplusplus
	tfx_vector_t<tfx_ribbon_segment_t> segments;
	tfx_vector_t<tfxU32> free_ribbons;
	tfx_vector_t<tfxU32> ribbon_emitter_indexes[2];
	tfx_vector_t<tfxU32> control_ribbon_queue;
	tfx_storage_map_t<tfxU32> cached_static_path_segments;
#else
	tfx_vector_t segments;
	tfx_vector_t free_ribbons;
	tfx_vector_t ribbon_emitter_indexes[2];
	tfx_vector_t control_ribbon_queue;
	tfx_storage_map_t cached_static_path_segments;
#endif
	tfxRibbonBucketFlags flags;
	tfxRibbonBucketComputeShaderType shader_type;
} tfx_ribbon_bucket_t;

typedef struct tfx_compute_fx_global_state_s {
	tfxU32 start_index;
	tfxU32 current_length;
	tfxU32 max_index;
	tfxU32 end_index;
}tfx_compute_fx_global_state_t;

typedef struct tfx_compute_controller_s {
	tfx_vec2_t position;
	float line_length;
	float angle_offset;
	tfx_vec4_t scale_rotation;              //Scale and rotation (x, y = scale, z = rotation, w = velocity_adjuster)
	float end_frame;
	tfxU32 normalised_values;				//Contains normalized values which are generally either 0 or 255, normalised in the shader to 0 and 1 (except opacity): age_rate, line_negator, spin_negator, position_negator, opacity
	tfxParticleControlFlags flags;
	tfxU32 image_data_index;				//index into the shape buffer on the gpu. CopyComputeShapeData must be called to prepare the data.
	tfx_vec2_t image_handle;
	tfx_vec2_t emitter_handle;
	float noise_offset;
	float stretch;
	float frame_rate;
	float noise_resolution;
} tfx_compute_controller_t;

typedef struct tfx_gpu_graph_data_s {
	tfx_vec4_t node_data;
	tfx_vec4_t oscillator;
	tfx_graph_easing_type easing_type;
	int flags;
	int padding[2];
} tfx_gpu_graph_data_t;

typedef struct tfx_compute_particle_s {
	tfx_vec2_t local_position;
	tfx_vec2_t base_size;

	float base_velocity;
	float base_roll_spin;
	float base_weight;

	float age;                            //The age of the particle, used by the controller to look up the current state on the graphs
	float max_age;                        //max age before the particle expires
	float emission_angle;                //Emission angle of the particle at spawn time

	float noise_offset;                    //The random velocity added each frame
	float noise_resolution;                //The random velocity added each frame
	float image_frame;
	tfxU32 control_slot_and_layer;    //index to the controller, and also stores the layer in the effect manager that the particle is on (layer << 3)
	float local_rotation;
}tfx_compute_particle_t;

typedef struct tfx_gpu_image_data_s {
	tfx_vec4_t uv;
	tfxU64 uv_packed;
	tfx_vec2_t image_size;
	tfxU32 texture_array_index;
	float animation_frames;
	float max_radius;
#ifdef tfxCUSTOM_GPU_IMAGE_DATA
	//add addition image data if needed
#endif
}tfx_gpu_image_data_t;

#ifdef __cplusplus
typedef struct tfx_gpu_shapes_s {
	tfxU32 magic;
	tfx_vector_t<tfx_gpu_image_data_t> list;
}tfx_gpu_shapes_t;
#else
typedef struct tfx_gpu_shapes_s tfx_gpu_shapes_t;
#endif
tfxMAKE_HANDLE(tfx_gpu_shapes);

typedef struct tfx_spawn_work_entry_s {
	tfx_random_t random;
	tfx_effect_manager pm;
	tfx_particle_emitter_properties_t *properties;
	tfx_shared_properties_t *shared_properties;
	tfx_parent_spawn_controls_t *parent_spawn_controls;
	tfxU32 emitter_index;
	tfxU32 parent_index;
	tfx_emission_type emission_type;
	tfxParticleEmitterFlags parent_property_flags;
	tfxEffectPropertyFlags root_effect_flags;
	tfx_particle_soa_t *particle_data;
#ifdef __cplusplus
	tfx_vector_t<tfx_depth_index_t> *depth_indexes;
#else
	tfx_vector_t *depth_indexes;
#endif
	tfxU32 depth_index_start;
	tfxU32 seed;
	float tween;
	tfxU32 max_spawn_count;
	tfxU32 amount_to_spawn;
	tfxU32 spawn_start_index;
	tfxU32 next_buffer;
	float qty_step_size;
	float highest_particle_age;
	float overal_scale;
    tfxU32 particle_uid;
}tfx_spawn_work_entry_t;

typedef struct tfx_ribbon_work_entry_s {
	tfx_random_t random;
	tfx_effect_manager pm;
	tfx_ribbon_bucket_t *ribbon_bucket;
	tfx_ribbon_emitter_properties_t *properties;
	tfx_shared_properties_t *shared_properties;
	tfxU32 parent_index;
	tfxRibbonEmitterFlags property_flags;
	tfxEffectPropertyFlags effect_flags;
	tfx_parent_spawn_controls_t *parent_spawn_controls;
	tfx_graph_list_t *graphs;
	tfxU32 new_ribbons;
	tfxU32 amount_to_spawn;
	float tween;
	float qty_step_size;
	float overal_scale;
}tfx_ribbon_work_entry_t;

typedef struct tfx_control_work_entry_s {
	float node_count;
	tfxU32 start_index;
	tfxU32 end_index;
	tfxU32 wide_end_index;
	tfxU32 start_diff;
	tfxU32 running_sprite_offset;
	tfxU32 sprites_index;
	tfxU32 cumulative_index_point;
	tfxU32 effect_instance_offset;
	tfxU32 sprite_buffer_end_index;
	tfxU32 emitter_index;
	tfx_effect_manager pm;
	tfx_graph_list_t *graphs;
	tfxU32 layer;
	tfx_particle_emitter_properties_t *properties;
	tfx_shared_properties_t *shared_properties;
	tfx_buffer_t *sprite_instances;
#ifdef __cplusplus
	tfx_vector_t<tfx_depth_index_t> *depth_indexes;
#else
	tfx_vector_t *depth_indexes;
#endif
	tfx_emitter_path_t *path;
	bool sample_path_life;
	float overal_scale;
	float global_stretch;
	float global_intensity;
	float global_noise;
}tfx_control_work_entry_t;

typedef struct tfx_control_ribbon_work_entry_s {
	tfxU32 segment_array_index;
	tfxU32 ribbon_index;
	tfx_effect_manager pm;
	tfx_random_t random;
	tfx_ribbon_bucket_t *ribbon_bucket;
} tfx_control_ribbon_work_entry_t;

typedef struct tfx_particle_age_work_entry_s {
	tfxU32 start_index;
	tfxU32 emitter_index;
	tfxU32 wide_end_index;
	tfxU32 start_diff;
	tfx_particle_emitter_properties_t *properties;
	tfx_shared_properties_t *shared_properties;
	tfx_effect_manager pm;
}tfx_particle_age_work_entry_t;

typedef struct tfx_sort_work_entry_s {
#ifdef __cplusplus
	tfx_bucket_array_t<tfx_particle_soa_t> *bank;
	tfx_vector_t<tfx_depth_index_t> *depth_indexes;
#else
	tfx_bucket_array_t *bank;
	tfx_vector_t *depth_indexes;
#endif
} tfx_sort_work_entry_t;

typedef struct tfx_compress_work_entry_s {
	tfx_sprite_data_t *sprite_data;
	tfxU32 frame;
} tfx_compress_work_entry_t;

typedef struct tfx_sprite_index_range_s {
	tfxU32 start_index;
	tfxU32 end_index;
	tfxU32 sprite_count;
} tfx_sprite_index_range_t;

typedef struct tfx_effect_data_s {
	tfxU32 *global_attributes;
	tfxU32 *transform_attributes;
	float *overal_scale;
	float *life;
	float *size_x;
	float *size_y;
	float *velocity;
	float *spin;
	float *intensity;
	float *splatter;
	float *weight;
}tfx_effect_data_t;

//An anim instance is used to let the gpu know where to draw an animation with sprite data. 48 bytes
typedef struct tfx_animation_instance_s {
	tfx_vec3_t position;                //position that the instance should be played at
	float scale;                        //Scales the overal size of the animation
	tfxU32 sprite_count;                //The number of instance_data to be drawn
	tfxU32 ribbon_count;                //The number of ribbons to be drawn
	tfxU32 frame_count;                 //The number of frames in the animation
	tfxU32 offset_into_sprite_data;     //The starting ofset in the buffer that contains all the sprite data
	tfxU32 offset_into_ribbon_data;     //The starting ofset in the buffer that contains all the ribbon data
	tfxU32 info_index;                  //Index into the effect_animation_info storage map to get at the frame meta
	float current_time;                 //Current point of time in the animation
	float animation_length_in_time;     //Total time that the animation lasts for
	float tween;                        //The point time within the frame (0..1)
	tfxAnimationInstanceFlags flags;    //Flags associated with the instance
}tfx_animation_instance_t;

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

typedef struct tfx_animation_emitter_properties_s {
	tfx_vec2_t handle;        //image handle
	tfxU32 handle_packed;
	tfxU32 flags;
	tfxU32 start_frame_index;
	tfxU32 color_ramp_index;
	float animation_frames;
	float padding;
}tfx_animation_emitter_properties_t;

typedef struct tfx_animation_ribbon_properties_s {
	tfxU32 segment_count;
	tfxU32 tessellation;
	tfxRibbonBucketComputeShaderType shader_type;
	tfxU32 color_ramp_index;
	tfxU32 segment_data_offset;
	float animation_frames;
	tfxU32 start_frame_index;
	tfxU32 flags;
	tfxU32 graph_lookup_offset;
} tfx_animation_ribbon_properties_t;

typedef struct tfx_color_ramp_bitmap_data_t {
#ifdef __cplusplus
	tfx_storage_map_t<tfxU32> color_ramp_ids;
	tfx_vector_t<tfx_bitmap_t> color_ramp_bitmaps;
#else
	tfx_storage_map_t color_ramp_ids;
	tfx_vector_t color_ramp_bitmaps;
#endif
	tfxU32 color_ramp_count;
}tfx_color_ramp_bitmap_data_t;

#ifdef __cplusplus
//Use the animation manager to control playing of pre-recorded effects
typedef struct tfx_animation_manager_s {
	tfxU32 magic;
	//All of the sprite data for all the animations that you might want to play on the GPU.
	//This could be deleted once it's uploaded to the GPU
	tfx_vector_t<tfx_sprite_instance_data_t> sprite_data;
	//List of active instances that are currently playing
	tfx_vector_t<tfx_animation_instance_t> instances;
	//List of instances in use. These index into the instances list above
	tfx_vector_t<tfxU32> instances_in_use[2];
	//List of free instance indexes
	tfx_vector_t<tfxU32> free_instances;
	//List of indexes into the instances list that will actually be sent to the GPU next frame
	//Any instances deemed not in view can be culled for example by not adding them to the queue
	tfx_vector_t<tfx_animation_instance_t> render_queue;
	//The compute shader needs to know when to switch from one animation instance to another as it
	//progresses through all the instance_data that need to be rendered. So this array of offsets tells
	//it when to do this. So 0 will always be in the first slot, then the second element will be 
	//the total number of instance_data to be drawn for the first animation instance and so on. When the
	//global index is more than or equal to the next element then we start on the next animation
	//instance and draw those instance_data.
	tfx_vector_t<tfxU32> offsets;
	//We also need to upload some emitter/ribbon properties to the GPU as well such as the sprite handle.
	//These can be looked up byt the sprite in the compute shader and the values applied to the sprite
	//before going to the vertex shader
	tfx_vector_t<tfx_animation_emitter_properties_t> emitter_properties;
	//Each animation has sprite data settings that contains properties about each animation
	tfx_vector_t<tfx_sprite_data_settings_t> sprite_data_settings;
	//Every animation that gets added to the animation manager gets info added here that describes
	//where to find the relevent sprite data in the buffer and contains other frame meta about the 
	//animation
	tfx_storage_map_t<tfx_sprite_data_metrics_t> effect_animation_info;
	//When loading in a tfxsd file the shapes are put here and can then be used to upload to the GPU
	//Other wise if you're adding sprite data from an effect library then the shapes will just be
	//referenced from there instead
	tfx_storage_map_t<tfx_image_data_t> particle_shapes;
	//Flips between 1 and 0 each frame to be used when accessing instances_in_use
	tfxU32 current_in_use_buffer;
	//This struct contains the size of the buffers that need to be uploaded to the GPU. Offsets and 
	//animation instances need to be uploaded every frame, but the sprite data only once before you
	//start drawing anything
	tfx_animation_buffer_metrics_t buffer_metrics;
	//Storage for the color ramps for all the effect sprite data in the animation manager
	tfx_color_ramp_bitmap_data_t color_ramps;
	//Bit flag field
	tfxAnimationManagerFlags flags;
	//The update frequency that the animations are recorded at. 60 is the recommended default
	float update_frequency;
	//Any pointer to user data that you want to use in callbacks such
	void *user_data;
	//Callback which you can assign in order to decide if an animation instance should be added to the render queue
	//the next frame. This callback is called inside the tfx_UpdateAnimationManager function. Set the callback
	//with SetAnimationManagerCallback
	bool((*maybe_render_instance_callback)(tfx_animation_manager animation_manager, tfx_animation_instance_t *instance, tfx_frame_meta_t *meta, void *user_data));

	//Ribbon management
	tfx_vector_t<tfxU32> ribbon_offsets;
	tfx_vector_t<tfx_ribbon_instance_data_t> ribbon_data;
	tfx_vector_t<tfx_ribbon_segment_t> ribbon_segment_data;
	tfx_vector_t<tfx_animation_ribbon_properties_t> ribbon_properties;
	tfx_vector_t<tfx_gpu_graph_data_t> ribbon_graph_data;
} tfx_animation_manager_t;
#endif

typedef struct tfx_effect_index_s {
	tfxEffectID index;
	float depth;
}tfx_effect_index_t;

//This struct is used for configuring a effect manager on creation
typedef struct tfx_effect_manager_info_s {
	tfxU32 max_particles;					//The maximum number of instance_data for each layer. This setting is not relevent if dynamic_sprite_allocation is set to true or group_sprites_by_effect is true.
	tfxU32 max_effects;                     //The maximum number of effects that can be updated at the same time.
	tfxU32 max_ribbon_segments;             //All segments for ribbons are stored in a single buffer. You will need to create buffers for rendering and so whatever you decide the max segments should be your buffers
											//should be big enough to contain all ribbon segments that you might need. You can call tfx_GetSegmentBufferSizeInBytes after creating the effect manager to get the byte
											//value that you can use to create the buffers. Also note that segments are always created in multiples of 32, so whatever number you put here it will be rounded to the
											//nearest multiple of 32.
	tfxU32 ribbon_tessellation;				//The amount of tessellation used for ribbons. Currently this is set globally. 1 is generally enough for most cases.
	//When set to false, all instance_data will be kept together in a large list.
	tfxU32 multi_threaded_batch_size;       //The size of each batch of particles to be processed when multithreading. Must be a power of 2 and 256 or greater.
	tfxU32 sort_passes;                     //when in order by depth mode (not guaranteed order) set the number of sort passes for more accuracy. Anything above 5 and you should just be guaranteed order.
	bool double_buffer_sprites;             //Set to true to double buffer instance_data so that you can interpolate between the old and new positions for smoother animations.
	bool dynamic_sprite_allocation;         //Set to true to automatically resize the sprite buffers if they run out of space. Not applicable when grouping instance_data by effect.
	bool group_sprites_by_effect;           //Set to true to group all instance_data by effect. Effects can then be drawn in specific orders or not drawn at all on an effect by effect basis.
	bool auto_order_effects;                //When group_sprites_by_effect is true then you can set this to true to sort the effects each frame. Use tfx_SetPMCamera in 3d to set the effect depth to the distance the camera.
	bool write_direct_to_staging_buffer;	//Make the effect manager write directly to the staging buffer. Use tfx_SetStagingBuffer before you call tfx_UpdateEffectManager
	void *user_data;						//User data that will get passed into the grow_staging_buffer_callback function which you can use to grow the buffer
	//If you need the staging buffer to be grown dynamically then you can use this call back to do that. It should return true if the buffer was successfully grown or false otherwise.
	bool(*grow_staging_buffer_callback)(tfxU32 new_size, tfx_effect_manager pm, void *user_data);
} tfx_effect_manager_info_t;

typedef struct tfx_ribbon_buffer_requirements_s {
	tfxU32 segment_buffer_size_in_bytes;
	tfxU32 ribbon_buffer_size_in_bytes;
	tfxU32 emitter_buffer_size_in_bytes;
} tfx_ribbon_buffer_requirements_t;

typedef struct tfx_ribbon_dispatch_s {
	tfx_ribbon_bucket_t *ribbon_data;
	tfxU32 index_offset;
	tfxU32 vertex_offset;
	tfxU32 index_count;
	tfxU32 vertex_count;
	tfxU32 ribbon_offset;
	tfxU32 segment_offset;
	tfxU32 total_segments;
} tfx_ribbon_dispatch_t;

//Use the effect manager to add multiple effects to your scene 
#ifdef __cplusplus
typedef struct tfx_effect_manager_s {
	tfxU32 magic;
	tfx_vector_t<tfx_soa_buffer_t> particle_array_buffers;
	tfx_bucket_array_t<tfx_particle_soa_t> particle_arrays;
	tfx_vector_t<tfx_soa_buffer_t> particle_location_buffers;
	tfx_bucket_array_t<tfx_spawn_points_soa_t> particle_location_arrays;
	tfx_storage_map_t<tfx_vector_t<tfxU32>> free_particle_lists;
	tfx_storage_map_t<tfx_vector_t<tfxU32>> free_particle_location_lists;
	tfx_storage_map_t<tfx_vector_t<tfxU32>> free_ribbon_segment_lists;
	tfx_storage_map_t<tfx_ribbon_bucket_t> ribbon_segment_buckets;
	//GPU compute particle buffer management
	tfx_storage_map_t<tfx_gpu_particle_group_t> gpu_groups;	//One entry per unique (profile_flags, life_bucket) combination

	//Only used when using distance from camera ordering. New particles are put in this list and then merge sorted into the particles buffer
	tfx_vector_t<tfx_sort_work_entry_t> sorting_work_entry;
	tfx_vector_t<tfx_spawn_work_entry_t> spawn_work;
	tfx_vector_t<tfx_ribbon_work_entry_t> ribbon_work;
	tfx_vector_t<tfx_control_work_entry_t> control_work;
	tfx_vector_t<tfx_control_ribbon_work_entry_t> ribbon_control_work;
	tfx_vector_t<tfx_particle_age_work_entry_t> age_work;
	tfx_vector_t<tfxParticleID> particle_indexes;
	tfx_vector_t<tfxU32> free_particle_indexes;
	tfx_vector_t<tfx_effect_index_t> effects_in_use[2];
	tfx_vector_t<tfxU32> control_emitter_queue;
	tfx_vector_t<tfxU32> emitters_check_capture;
	tfx_vector_t<tfx_effect_index_t> free_effects;
	tfx_vector_t<tfxU32> free_emitters;
	tfx_vector_t<tfxU32> free_gpu_emitters;
	tfx_vector_t<tfxU32> free_ribbon_emitters;
	tfx_vector_t<tfxU32> free_path_quaternions;
	tfx_vector_t<tfx_path_quaternion_t *> path_quaternions;
	tfx_vector_t<tfx_effect_state_t> effects;
	tfx_vector_t<tfx_particle_emitter_state_t> emitters;
	tfx_vector_t<tfx_gpu_emitter_t> gpu_emitters;
	tfx_vector_t<tfx_ribbon_emitter_state_t> ribbon_emitters;
	tfx_vector_t<tfx_spawn_work_entry_t *> deffered_spawn_work;
	tfx_vector_t<tfx_ribbon_work_entry_t *> deffered_ribbon_spawn_work;
	tfx_vector_t<tfx_unique_sprite_id_t> unique_sprite_ids[2][tfxLAYERS];
	tfx_vector_t<tfxU32> free_compute_controllers;

	tfx_work_queue_t work_queue;
	//The info config that was used to initialise the effect manager. This can be used to alter and the reconfigure the effect manager
	tfx_effect_manager_info_t info;
	//Banks of instance_data. All emitters write their sprite data to these banks. 
	tfx_buffer_t instance_buffer;
	tfx_buffer_t instance_buffer_for_recording[2][tfxLAYERS];
	tfxU32 current_sprite_buffer;
	tfxU32 highest_depth_index;

	//todo: document compute controllers once we've established this is how we'll be doing it.
	void *compute_controller_ptr;
	tfxU32 new_compute_particle_index;
	tfxU32 new_particles_count;
	void *new_compute_particle_ptr;
	//The maximum number of effects that can be updated per frame in the effect manager. If you're running effects with particles that have sub effects then this number might need 
	//to be relatively high depending on your needs. Use Init to udpate the sizes if you need to. Best to call Init at the start with the max numbers that you'll need for your application and don't adjust after.
	tfxU32 max_effects;
	//The maximum number of particles that can be updated per frame per layer. #define tfxLAYERS to set the number of allowed layers. This is currently 4 by default
	tfxU32 max_cpu_particles_per_layer[tfxLAYERS];
	//The maximum number of particles that can be updated per frame per layer in the compute shader. #define tfxLAYERS to set the number of allowed layers. This is currently 4 by default
	tfxU32 max_new_compute_particles;
	//The current effect buffer in use, can be either 0 or 1
	tfxU32 current_ebuff;
	//For looping through active effects with GetNextEffect function
	tfxU32 effect_index_position;
	tfx_ribbon_buffer_requirements_t ribbon_buffer_requirements;
	tfxU32 current_ribbon_count;

	tfxU32 effects_start_size;

	tfxU32 layer_sizes[tfxLAYERS];
	tfxU32 running_ribbon_vertex_count;

	int mt_batch_size;
	//We might not need these now.
	tfx_sync_t particle_index_mutex;
	tfx_sync_t add_effect_mutex;

	tfx_random_t random;
	tfx_random_t threaded_random;
	tfxU32 max_compute_controllers;
	tfxU32 highest_compute_controller_index;
	tfx_compute_fx_global_state_t compute_global_state;
	tfx_lookup_mode lookup_mode;
	//For when particles are ordered by distance from camera
	tfx_vec3_t camera_front;
	tfx_vec3_t camera_position;

	tfxU32 unique_particle_id;    //Used when recording sprite data
	tfxU32 unique_ribbon_id;      //Used when recording ribbon data
	//When using single particles, you can flag the emitter to set the max_age of the particle to the 
	//length in time of the animation so that it maps nicely to the animation
	float animation_length_in_time;

	//These can possibly be removed at some point, they're debugging variables
	tfxU32 particle_id;

	float gpu_current_time_ms;										//Running absolute time in ms, used for group tracking head bumps

	tfxEffectManagerFlags flags;
	//The length of time that passed since the last time Update() was called
	float frame_length;
	//You can cap the frame length to a maximum amount which I put in mainly for when stepping through 
	//with a debugger and you don't want to advance the particles too much because obviously a lot of
	//time is passing between frames because you're stepping through the code. Default is 240ms.
	float max_frame_length;
	tfxWideFloat frame_length_wide;
	float update_time;
	tfxWideFloat update_time_wide;
	float update_frequency;
} tfx_effect_manager_t;
#endif

typedef struct tfx_effect_library_stats_s {
	tfxU32 total_effects;
	tfxU32 total_particle_emitters;
	tfxU32 total_ribbon_emitters;
	tfxU32 total_node_lookup_indexes;
	tfxU32 total_shapes;
	tfxU64 required_graph_node_memory;
	tfxU64 required_graph_lookup_memory;
	tfxU32 reserved1;
	tfxU32 reserved2;
	tfxU32 reserved3;
	tfxU32 reserved4;
	tfxU32 reserved5;
	tfxU32 reserved6;
	tfxU32 reserved7;
}tfx_effect_library_stats_t;

#ifdef __cplusplus
typedef struct tfx_library_s {
	tfxU32 magic;
	tfxErrorFlags error_flags;
	tfx_storage_map_t<tfx_effect_descriptor> effect_paths;
	tfx_vector_t<tfx_effect_descriptor> effects;
	tfx_storage_map_t<tfx_image_data_t> particle_shapes;
	tfx_vector_t<tfx_particle_emitter_properties_t> emitter_properties;
	tfx_vector_t<tfx_shared_properties_t> shared_properties;
	tfx_vector_t<tfx_ribbon_emitter_properties_t> ribbon_properties;
	tfx_storage_map_t<tfx_sprite_data_t> pre_recorded_effects;

	tfx_bucket_array_t<tfx_emitter_path_t> paths;
	tfx_vector_t<tfx_graph_list_t> graphs;
	tfx_vector_t<tfx_sprite_sheet_settings_t> sprite_sheet_settings;
	tfx_vector_t<tfx_sprite_data_settings_t> sprite_data_settings;
	tfx_vector_t<tfx_preview_camera_settings_t> preview_camera_settings;

	tfx_vector_t<tfxU32> free_graph_lists;
	tfx_vector_t<tfxU32> free_animation_settings;
	tfx_vector_t<tfxU32> free_preview_camera_settings;
	tfx_vector_t<tfxU32> free_particle_emitter_properties;
	tfx_vector_t<tfxU32> free_shared_emitter_properties;
	tfx_vector_t<tfxU32> free_ribbon_emitter_properties;
	tfx_vector_t<tfxU32> free_infos;
	tfx_vector_t<tfxU32> free_keyframes;

	tfx_gpu_shapes gpu_shapes;
	tfx_color_ramp_bitmap_data_t color_ramps;
	//Get an effect from the library by index
	tfx_str256_t name;
	bool open_library;
	bool dirty;
	tfx_stream_t library_file_path;
	tfxU32 uid;
	void(*uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset);
} tfx_library_t;
#endif

#ifdef __cplusplus
typedef struct tfx_effect_template_s {
	tfxU32 magic;
	tfx_storage_map_t<tfx_effect_descriptor > paths;
	tfx_effect_descriptor effect;
	tfx_effect_descriptor original_effect;
}tfx_effect_template_t;
#endif

typedef struct tfx_data_entry_s {
	tfx_data_type type;
	tfx_str64_t key;
	tfx_str512_t str_value;
	int int_value;
	tfx_rgba8_t color_value;
	bool bool_value;
	float float_value;
	double double_value;
}tfx_data_entry_t;

//------------------------------------------------------------
//Section: Internal_Functions
//------------------------------------------------------------

#ifdef __cplusplus
tfxINTERNAL void tfx__resize_particle_soa_callback(tfx_soa_buffer_t *buffer, tfxU32 index);

//--------------------------------
//Internal functions used either by the library or editor
//--------------------------------
#define tfxMakeColorRampIndex(layer, y) ((layer & 0x00007FFF) << 16) + y
#define tfxColorRampIndex(id) (id & 0x0000FFFF)
#define tfxColorRampLayer(id) ((id & 0x7FFF0000) >> 16)
#define tfxFlagColorRampIDAsEdited(id) id |= 0x80000000
#define tfxUnFlagColorRampIDAsEdited(id) id &= ~0x80000000
#define tfxColorRampIsEdited(id) (id & 0x80000000)
#define tfxRibbonBucketID(bucket_info) (tfxKey)bucket_info.segment_count | ((tfxKey)bucket_info.shader_type << 16)
tfxINTERNAL inline tfxParticleID tfx__make_particle_id(tfxU32 bank_index, tfxU32 particle_index) { return ((bank_index & 0x00000FFF) << 20) + particle_index; }
tfxINTERNAL inline tfxU32 tfx__particle_index(tfxParticleID id) { return id & 0x000FFFFF; }
tfxINTERNAL inline tfxU32 tfx__particle_bank(tfxParticleID id) { return (id & 0xFFF00000) >> 20; }
tfxINTERNAL tfxU32 tfx__grab_particle_lists(tfx_effect_manager pm, tfxKey emitter_hash, tfxU32 reserve_amount, tfxEmitterControlProfileFlags flags);
tfxINTERNAL tfxU32 tfx__grab_gpu_emitter(tfx_effect_manager pm);
tfxINTERNAL void tfx__free_gpu_emitter(tfx_effect_manager pm, tfxU32 index);
tfxINTERNAL tfxU32 tfx__grab_ribbon(tfx_effect_manager pm, tfx_ribbon_bucket_t *bucket, tfx_ribbon_emitter_state_t *segment_count);
tfxINTERNAL void tfx__free_ribbon(tfx_effect_manager pm, tfxKey bucket_id, tfxU32 ribbon_index);
tfxINTERNAL tfxU32 tfx__grab_particle_location_lists(tfx_effect_manager pm, tfxKey emitter_hash, tfxU32 reserve_amount);
tfxINTERNAL void tfx__init_ribbon_segment_buffer(tfx_effect_manager pm, tfxKey bucket_id, tfx_ribbon_bucket_info_t *bucket_info, int tessellation);
tfxAPI_EDITOR tfx_ribbon_buffer_info_t tfx__generate_ribbon_buffer_info(tfxU32 tessellation);
tfxAPI_EDITOR void tfx__update_ribbon_bucket_id(tfx_effect_descriptor ribbon_emitter);

//--------------------------------
//Profilings
//--------------------------------
tfxAPI_EDITOR void tfx__gather_stats(tfx_profile_t *profile, tfx_profile_stats_t *stat);
tfxAPI_EDITOR void tfx__reset_snap_shot(tfx_profile_snapshot_t *snapshot);
tfxAPI_EDITOR void tfx__reset_snap_shots();
tfxAPI_EDITOR void tfx__dump_snapshots(tfx_storage_map_t<tfx_vector_t<tfx_profile_snapshot_t>> *profile_snapshots, tfxU32 amount);

//--------------------------------
//Reading/Writing files
//--------------------------------
tfxAPI_EDITOR void tfx__read_entire_file(const char *file_name, tfx_stream buffer, bool terminate = false);
tfxAPI_EDITOR tfxErrorFlags tfx__load_package_file(const char *file_name, tfx_package package);
tfxAPI_EDITOR tfxErrorFlags tfx__load_package_stream(tfx_stream stream, tfx_package package);
tfxAPI_EDITOR tfx_package tfx__create_package(const char *file_path);
tfxAPI_EDITOR bool tfx__save_package_disk(tfx_package package);
tfxAPI_EDITOR tfx_stream tfx__save_package_memory(tfx_package package);
tfxAPI_EDITOR tfx_package_entry_info_t *tfx__get_package_file(tfx_package package, const char *name);
tfxAPI_EDITOR void tfx__add_entry_to_package(tfx_package package, tfx_package_entry_info_t file);
tfxAPI_EDITOR bool tfx__file_exists_in_package(tfx_package package, const char *file_name);
tfxAPI_EDITOR void tfx__free_package(tfx_package package);
tfxAPI_EDITOR void tfx__copy_data_to_stream(tfx_stream dst, const void *src, tfxU64 size);
tfxAPI_EDITOR void tfx__copy_stream(tfx_stream dst, tfx_stream src);
tfxINTERNAL tfxU64 tfx__get_package_size(tfx_package package);
tfxINTERNAL bool tfx__validate_package(tfx_package package);
tfxINTERNAL void tfx__add_file_to_package(tfx_package package, const char *file_name, tfx_stream data);

//Some file IO functions for the editor
tfxAPI_EDITOR bool tfx__has_data_value(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key);
tfxAPI_EDITOR void tfx__add_data_value_str(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key, const char *value);
tfxAPI_EDITOR void tfx__add_data_value_int(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key, int value);
tfxAPI_EDITOR void tfx__add_color_value(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key, tfx_rgba8_t value);
tfxAPI_EDITOR void tfx__add_data_value_bool(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key, bool value);
tfxAPI_EDITOR void tfx__add_data_value_float(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key, float value);
tfxAPI_EDITOR const char* tfx__get_data_str_value(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key);
tfxAPI_EDITOR tfx_rgba8_t tfx__get_data_color_value(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key);
tfxAPI_EDITOR float tfx__get_data_float_value(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key);
tfxAPI_EDITOR bool tfx__save_data_file(tfx_storage_map_t<tfx_data_entry_t> *config, const char *path = "");
tfxAPI_EDITOR bool tfx__load_data_file(tfx_data_types_dictionary_t *data_types, tfx_storage_map_t<tfx_data_entry_t> *config, const char *path);
tfxAPI_EDITOR void tfx__stream_particle_emitter_properties(tfx_effect_descriptor emitter, tfx_shared_properties_t *shared, tfx_particle_emitter_properties_t *properties, tfxSharedEmitterFlags shared_flags, tfxParticleEmitterFlags flags, tfx_stream_t *file);
tfxAPI_EDITOR void tfx__stream_ribbon_emitter_properties(tfx_effect_descriptor emitter, tfx_shared_properties_t *shared, tfx_ribbon_emitter_properties_t *ribbon_properties, tfxSharedEmitterFlags shared_flags, tfxRibbonEmitterFlags flags, tfx_stream_t *file);
tfxAPI_EDITOR void tfx__stream_effect_properties(tfx_effect_descriptor effect, tfx_stream_t *file);
tfxAPI_EDITOR void tfx__stream_path_properties(tfx_effect_descriptor effect, tfx_stream_t *file);
tfxAPI_EDITOR void tfx__stream_path_nodes(tfx_effect_descriptor effect, tfx_stream_t *file);
tfxAPI_EDITOR void tfx__stream_graph(const char *name, tfx_effect_descriptor descriptor, tfx_graph_t *graph, tfx_stream_t *file);
tfxAPI_EDITOR void tfx__stream_graph_properties(const char *name, tfx_effect_descriptor descriptor, tfx_graph_t *graph, tfx_stream_t *file);
tfxAPI_EDITOR void tfx__split_string_stack(const char *s, int length, tfx_vector_t<tfx_str256_t> *pair, char delim = 61);
tfxAPI_EDITOR bool tfx__string_is_uint(const char *s);
tfxAPI_EDITOR bool tfx__line_is_uint(tfx_line_t *line);
tfxAPI_EDITOR tfx_str256_t tfx__get_property_as_string(tfx_effect_descriptor effect, tfx_str256_t property_name);
tfxAPI_EDITOR tfx_str256_t tfx__get_graph_property_as_string(tfx_graph_t *graph, tfx_str256_t property_name);
tfxAPI_EDITOR tfx_str64_t tfx__graph_type_to_property_string(tfx_graph_type graph_type);
tfxAPI_EDITOR tfx_stream_t tfx__get_graph_as_string(tfx_effect_descriptor effect, tfx_graph_t *graph, bool include_property_name);
tfxAPI_EDITOR void tfx__add_graph_to_stream(tfx_effect_descriptor effect, tfx_stream_t *stream, tfx_graph_t *graph, bool include_property_name);
tfxAPI_EDITOR void tfx__assign_property_from_string(tfx_effect_descriptor effect, tfx_str256_t property_name, const char *value);
tfxAPI_EDITOR void tfx__assign_property_line(tfx_effect_descriptor effect, tfx_vector_t<tfx_str256_t> *pair, tfxU32 file_version);
tfxAPI_EDITOR void tfx__assign_effector_property_u32(tfx_effect_descriptor effect, tfx_str256_t *field, tfxU32 value, tfxU32 file_version);
tfxAPI_EDITOR void tfx__assign_effector_property(tfx_effect_descriptor effect, tfx_str256_t *field, float value);
tfxAPI_EDITOR void tfx__assign_effector_property_bool(tfx_effect_descriptor effect, tfx_str256_t *field, bool value);
tfxAPI_EDITOR void tfx__assign_effector_property_int(tfx_effect_descriptor effect, tfx_str256_t *field, int value);
tfxAPI_EDITOR void tfx__assign_effector_property_str(tfx_effect_descriptor effect, tfx_str256_t *field, const char *value);
tfxAPI_EDITOR void tfx__assign_graph_node_data(tfx_effect_descriptor effect, tfx_vector_t<tfx_str256_t> *values);
tfxAPI_EDITOR void tfx__assign_graph_properties(tfx_effect_descriptor effect, tfx_vector_t<tfx_str256_t> *values);
tfxAPI_EDITOR tfx_hsv_t tfx__rgb_to_hsv(tfx_rgb_t in);
tfxAPI_EDITOR tfx_rgb_t tfx__hsv_to_rgb(tfx_hsv_t in);
tfxAPI_EDITOR tfx_rgba8_t tfx__convert_float_color(float color_array[4]);
tfxAPI_EDITOR float tfx__length_vec3(tfx_vec3_t const *v);
tfxAPI_EDITOR float tfx__dot_product_vec2(const tfx_vec2_t *a, const tfx_vec2_t *b);
tfxAPI_EDITOR void tfx__catmull_rom_spline_3d_soa(const float *p_x, const float *p_y, const float *p_z, int p0, float t, float vec[3]);
tfxAPI_EDITOR void tfx__catmull_rom_spline_3d(const tfx_vec4_t *p0, const tfx_vec4_t *p1, const tfx_vec4_t *p2, const tfx_vec4_t *p3, float t, float vec[3]);
tfxAPI_EDITOR void tfx__catmull_rom_spline_gradient_3d(const tfx_vec4_t *p0, const tfx_vec4_t *p1, const tfx_vec4_t *p2, const tfx_vec4_t *p3, float t, float vec[3]);
tfxAPI_EDITOR float tfx__vec2_length_fast(tfx_vec2_t const *v);
tfxAPI_EDITOR float tfx__vec3_length_fast(tfx_vec3_t const *v);
tfxAPI_EDITOR void tfx__wide_transform_quaternion_vec3(const tfx_quaternion_t *q, tfxWideFloat *x, tfxWideFloat *y, tfxWideFloat *z);
tfxAPI_EDITOR tfxU32 tfx__pack16bit_sscaled(float x, float y, float max_value);
tfxAPI_EDITOR tfxU32 tfx__pack16bit_unorm(float x, float y);
tfxAPI_EDITOR void tfx__transform_3d(tfx_vec3_t *out_rotations, tfx_vec3_t *out_local_rotations, float *out_scale, tfx_vec3_t *out_position, tfx_vec3_t *out_local_position, tfx_vec3_t *out_translation, tfx_quaternion_t *out_q, tfx_effect_state_t *parent);
tfxAPI_EDITOR void tfx__update_emitter_control_profile(tfx_effect_descriptor emitter);
tfxAPI_EDITOR tfx_mat3_t tfx__create_matrix3(float v = 1.f);
tfxAPI_EDITOR tfx_mat3_t tfx__rotate_matrix3(tfx_mat3_t const *m, float r);
tfxAPI_EDITOR void tfx__split_string_vec(const char *s, int length, tfx_vector_t<tfx_str256_t> *pair, char delim = 61);
tfxINTERNAL tfx_noise_type tfx__get_emitter_noise_type(tfx_effect_descriptor emitter);
tfxINTERNAL void tfx__update_library_control_profiles(tfx_library library);
tfxINTERNAL	tfx_line_t tfx__read_line(const char *s);
tfxAPI_EDITOR tfxU32 tfx__pack8bit_xyz(float const &v_x, float const &v_y, float const &v_z);
tfxINTERNAL tfxU64 tfx__pack16bit_quaternion(tfx_quaternion_t v);
tfxAPI_EDITOR tfxU64 tfx__pack16bit_quaternion_for_gpu(tfx_quaternion_t q);
tfxAPI_EDITOR tfx_quaternion_t tfx__unpack16bit_quaternion_from_gpu(tfxU64 q);
tfxINTERNAL tfxWideInt tfx__wide_pack8bitunorm_xyz(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z);
tfxINTERNAL void tfx__wide_unpack16bit(tfxWideInt xy, tfxWideInt zw, tfxWideFloat &x, tfxWideFloat &y, tfxWideFloat &z, tfxWideFloat &w);
tfxINTERNAL tfx_quaternion_t tfx__unpack16bit_quaternion(tfxU64 in);
tfxINTERNAL tfx_vec3_t tfx__get_emission_direciton_3d(tfx_effect_manager pm, tfx_library library, tfx_random_t *random, tfx_particle_emitter_state_t &emitter, float emission_pitch, float emission_yaw, tfx_vec3_t local_position, tfx_vec3_t world_position);
tfxINTERNAL tfx_quaternion_t tfx__get_path_rotation_3d(tfx_random_t *random, float range, float pitch, float yaw, bool y_axis_only);
tfxINTERNAL tfx_vec3_t tfx__cylinder_surface_normal(float x, float z, float width, float depth);
tfxINTERNAL tfx_vec3_t tfx__ellipse_surface_normal(float x, float y, float z, float width, float height, float depth);
tfxAPI_EDITOR tfx_vec3_t tfx__catmull_rom_spline_gradient_3d_soa(const float *px, const float *py, const float *pz, float t);
tfxAPI_EDITOR void tfx__wide_catmull_rom_spline_3d(tfxWideArrayi *pi, tfxWideFloat t, float *x, float *y, float *z, tfxWideFloat *vx, tfxWideFloat *vy, tfxWideFloat *vz);
tfxINTERNAL float tfx__length_vec3_nosqr(tfx_vec3_t const *v);
tfxINTERNAL void tfx__assign_effector_property_u64(tfx_effect_descriptor effect, tfx_str256_t *field, tfxU64 value, tfxU32 file_version);
tfxINTERNAL void tfx__add_data_value_double(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key, double value);
tfxINTERNAL void tfx__add_color_value_from_int(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key, tfxU32 value);
tfxINTERNAL int tfx__get_data_int_value(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key);

//--------------------------------
//Graph functions
//Mainly used by the editor to edit graphs so these are kind of API functions but you wouldn't generally use these outside of the particle editor
//--------------------------------
tfxAPI_EDITOR void tfx__init_paths_soa(tfx_soa_buffer_t *buffer, tfx_path_nodes_soa_t *soa, tfxU32 reserve_amount);
tfxAPI_EDITOR void tfx__init_emitter_properties(tfx_particle_emitter_properties_t *properties);
tfxAPI_EDITOR void tfx__init_shared_properties(tfx_shared_properties_t *shared_properties);
tfxAPI_EDITOR tfx_attribute_node_t *tfx__add_graph_node_values(tfx_graph_t *graph, float frame, float value, tfxAttributeNodeFlags flags = 0, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
tfxAPI_EDITOR tfx_attribute_node_t *tfx__append_graph_node_values(tfx_graph_t *graph, float frame, float value, tfxAttributeNodeFlags flags = 0, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
tfxAPI_EDITOR float tfx__get_graph_value_by_age(tfx_graph_t *graph, float age);
tfxAPI_EDITOR float tfx__get_linear_graph_value_by_percent_of_life(tfx_graph_t *graph, float t);
tfxAPI_EDITOR tfx_attribute_node_t *tfx__get_graph_last_node(tfx_graph_t *graph);
tfxAPI_EDITOR tfx_attribute_node_t *tfx__get_graph_first_node(tfx_graph_t *graph);
tfxAPI_EDITOR float tfx__get_graph_first_value(tfx_graph_t *graph);
tfxAPI_EDITOR tfx_attribute_node_t *tfx__insert_graph_node(tfx_graph_t *graph, float, float);
tfxAPI_EDITOR float *tfx__link_graph_first_value(tfx_graph_t *graph);
tfxAPI_EDITOR float *tfx__link_graph_last_value(tfx_graph_t *graph);
tfxAPI_EDITOR float tfx__get_graph_last_value(tfx_graph_t *graph);
tfxAPI_EDITOR float tfx__graph_value_by_index(tfx_graph_t *graph, tfxU32 index);
tfxAPI_EDITOR tfx_attribute_node_t *tfx__find_graph_node(tfx_graph_t *graph, tfx_attribute_node_t *n);
tfxAPI_EDITOR tfx_attribute_node_t *tfx__find_graph_node_by_uid(tfx_graph_t *graph, tfxU32 uid);
tfxAPI_EDITOR void tfx__validate_graph_curves(tfx_graph_t *graph);
tfxAPI_EDITOR void tfx__delete_graph_node(tfx_graph_t *graph, tfx_attribute_node_t *n);
tfxAPI_EDITOR void tfx__reset_graph(tfx_graph_t *graph, float first_node_value, tfx_graph_preset preset, bool add_node = true, float max_frames = 0);
tfxAPI_EDITOR void tfx__reset_graph_nodes(tfx_graph_t *graph, float first_node_value, tfx_graph_preset preset, bool add_node = true);
tfxAPI_EDITOR void tfx__clear_lerp_graph(tfx_graph_t *graph);
tfxAPI_EDITOR void tfx__clear_graph_to_one(tfx_graph_t *graph, float value);
tfxAPI_EDITOR void tfx__free_graph(tfx_graph_t *graph);
tfxAPI_EDITOR void tfx__copy_graph(tfx_graph_t *graph, tfx_graph_t *to, bool include_types);
tfxAPI_EDITOR void tfx__copy_graph_color(tfx_graph_list_t *from, tfx_graph_list_t *to, tfx_effect_descriptor_type from_type, tfx_effect_descriptor_type to_type);
tfxAPI_EDITOR void tfx__copy_graph_colors(tfx_graph_t *from_red, tfx_graph_t *from_blue, tfx_graph_t *from_green, tfx_graph_t *to_red, tfx_graph_t *to_green, tfx_graph_t *to_blue);
tfxAPI_EDITOR bool tfx__sort_graph(tfx_graph_t *graph);
tfxAPI_EDITOR void tfx__flip_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_blend_factor_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_lerp_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_overtime_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_overlength_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_factor_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_global_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_angle_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_translation_graph(tfx_graph_t *graph);
tfxAPI_EDITOR void tfx__multiply_all_graph_values(tfx_graph_t *graph, float scalar);
tfxAPI_EDITOR void tfx__drag_graph_values(tfx_graph_preset preset, float *frame, float *value);
tfxAPI_EDITOR void tfx__update_lerp_graph(tfx_graph_t *graph);
tfxAPI_EDITOR void tfx__update_lerp_graphs_of_effect(tfx_effect_descriptor effect, bool include_children);
tfxAPI_EDITOR void tfx__update_color_ramp(tfx_graph_list_t *graph_list, tfx_color_ramp_t *ramp, float gamma = tfxGAMMA);
tfxAPI_EDITOR bool tfx__edit_color_ramp_bitmap(tfx_library library, tfx_graph_list_t *graph_list);
tfxAPI_EDITOR void tfx__reindex_graph(tfx_graph_t *graph);
tfxAPI_EDITOR float tfx__get_graph_max_value(tfx_graph_t *graph);
tfxAPI_EDITOR tfx_vec2_t tfx__get_max_graph_values(tfx_graph_preset preset);
tfxAPI_EDITOR tfx_vec2_t tfx__get_min_graph_values(tfx_graph_preset preset);
tfxAPI_EDITOR void tfx__update_graph_wide_oscillator(tfx_graph_t *graph);
tfxINTERNAL void tfx__add_graph_node(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxINTERNAL void tfx__set_graph_node(tfx_graph_t *graph, tfxU32 index, float frame, float value, tfxAttributeNodeFlags flags = 0, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
tfxINTERNAL float tfx__get_graph_random_value(tfx_graph_t *graph, float age, tfx_random_t *seed);
tfxINTERNAL tfx_attribute_node_t *tfx__get_graph_next_node(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxINTERNAL tfx_attribute_node_t *tfx__get_graph_prev_node(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxINTERNAL tfx_attribute_node_t *tfx__add_graph_coord_node(tfx_graph_t *graph, float, float);
tfxINTERNAL float tfx__get_graph_min_value(tfx_graph_t *graph);
tfxINTERNAL float tfx__get_graph_last_frame(tfx_graph_t *graph, float udpate_frequence);
tfxINTERNAL tfx_attribute_node_t *tfx__graph_node_by_index(tfx_graph_t *graph, tfxU32 index);
tfxINTERNAL float tfx__graph_frame_by_index(tfx_graph_t *graph, tfxU32 index);
tfxINTERNAL void tfx__delete_graph_node_at_frame(tfx_graph_t *graph, float frame);
tfxINTERNAL void tfx__clear_graph(tfx_graph_t *graph);
tfxINTERNAL bool tfx__color_graph(tfx_graph_t *graph);
tfxINTERNAL bool tfx__gpu_overtime_graph(tfx_graph_t *graph);
tfxINTERNAL inline float tfx__get_vector_angle(float x, float y) { return atan2f(x, -y); }
tfxINTERNAL bool tfx__compare_nodes(tfx_attribute_node_t *left, tfx_attribute_node_t *right);
tfxINTERNAL tfxKey tfx__hash_color_ramp(tfx_color_ramp_t *ramp);
tfxINTERNAL tfx_bitmap_t tfx__create_bitmap(int width, int height, int channels);
tfxINTERNAL void tfx__plot_bitmap(tfx_bitmap_t *image, int x, int y, tfx_rgba8_t color);
tfxINTERNAL void tfx__free_bitmap(tfx_bitmap_t *bitmap);
tfxINTERNAL void tfx__plot_color_ramp(tfx_bitmap_t *bitmap, tfx_color_ramp_t *ramp, tfxU32 y);
tfxINTERNAL void tfx__create_color_ramp_bitmaps(tfx_library library);
tfxINTERNAL void tfx__maybe_insert_color_ramp_bitmap(tfx_library library, tfx_graph_list_t *list);
tfxINTERNAL tfxU32 tfx__add_color_ramp_to_bitmap(tfx_color_ramp_bitmap_data_t *ramp_data, tfx_color_ramp_t *ramp);
tfxINTERNAL void tfx__copy_emitter_color_ramp_to_animation_manager(tfx_animation_manager animation_manager, tfxU32 properties_index, tfx_color_ramp_t *ramp);
tfxINTERNAL void tfx__copy_ribbon_color_ramp_to_animation_manager(tfx_animation_manager animation_manager, tfxU32 properties_index, tfx_color_ramp_t *ramp);
tfxINTERNAL float tfx__get_max_life(tfx_effect_descriptor e);
tfxINTERNAL float tfx__sample_multi_node_graph(tfx_graph_t *graph, float frame, float osc_t);

//---- GPU compute particle buffer management functions ----
tfxINTERNAL tfxU32 tfx__compute_max_gpu_particles(tfx_effect_descriptor child);
tfxINTERNAL tfxU32 tfx__find_or_create_gpu_group(tfx_effect_manager pm, tfxU32 property_index, float max_life);
tfxINTERNAL void   tfx__assign_emitter_to_gpu_group(tfx_effect_manager pm, tfxU32 emitter_index);
tfxINTERNAL tfxU32 tfx__gpu_group_record_spawns(tfx_effect_manager pm, tfxU32 emitter_index, tfxU32 count, float current_time_ms);
tfxINTERNAL void   tfx__tick_gpu_groups(tfx_effect_manager pm, float current_time_ms);
tfxINTERNAL void   tfx__free_gpu_groups(tfx_effect_manager pm);
tfxINTERNAL void   tfx__clear_gpu_groups(tfx_effect_manager pm);
//---- end GPU compute particle buffer management functions ----

//Node Manipulation
tfxAPI_EDITOR void tfx__unset_curves(tfx_graph_t *graph, tfxU32 index);
tfxAPI_EDITOR bool tfx__set_node(tfx_graph_t *graph, tfx_attribute_node_t *node, float *frame, float *value);
tfxAPI_EDITOR void tfx__set_adjacent_node_curves(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxAPI_EDITOR void tfx__set_node_curve(tfx_graph_t *graph, tfx_attribute_node_t *node, bool is_left_curve, float *frame, float *value);
tfxAPI_EDITOR bool tfx__move_node(tfx_graph_t *graph, tfx_attribute_node_t *node, float frame, float value, bool sort = true);
tfxAPI_EDITOR void tfx__clamp_graph_nodes(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_gpu_graph_type(tfx_graph_type type);
tfxAPI_EDITOR bool tfx__is_overtime_graph_type(tfx_graph_type type);
tfxAPI_EDITOR bool tfx__is_overtime_percentage_graph_type(tfx_graph_type type);
tfxAPI_EDITOR bool tfx__is_global_graph_type(tfx_graph_type type);
tfxAPI_EDITOR bool tfx__is_transform_graph_type(tfx_graph_type type);
tfxAPI_EDITOR bool tfx__is_global_percentage_graph_type(tfx_graph_type type);
tfxAPI_EDITOR bool tfx__is_emitter_size_graph_type(tfx_graph_type type);
tfxAPI_EDITOR bool tfx__is_angle_graph_type(tfx_graph_type type);
tfxAPI_EDITOR bool tfx__is_angle_overtime_graph_type(tfx_graph_type type);
tfxAPI_EDITOR bool tfx__is_everythine_else_graph_type(tfx_graph_type type);
tfxAPI_EDITOR bool tfx__has_node_at_frame(tfx_graph_t *graph, float frame);
tfxAPI_EDITOR bool tfx__has_more_than_one_key_frame(tfx_effect_descriptor e);
tfxAPI_EDITOR bool tfx__is_node_curve(tfx_attribute_node_t *node);
tfxAPI_EDITOR bool tfx__node_curves_are_initialised(tfx_attribute_node_t *node);
tfxAPI_EDITOR bool tfx__set_node_curve_initialised(tfx_attribute_node_t *node);
tfxAPI_EDITOR bool tfx__is_color_graph_type(tfx_graph_type type);
tfxAPI_EDITOR void tfx__set_node_curve_frames(tfx_graph_t *graph);
tfxINTERNAL void tfx__clamp_node(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxINTERNAL void tfx__clamp_node_curve(tfx_graph_t *graph, tfx_vec2_t *curve, tfx_attribute_node_t *node);
tfxINTERNAL void tfx__set_left_node_curve(tfx_graph_t *graph, tfx_vec2_t *curve, tfx_attribute_node_t *node);
tfxINTERNAL void tfx__set_right_node_curve(tfx_graph_t *graph, tfx_vec2_t *curve, tfx_attribute_node_t *node);
tfxINTERNAL bool tfx__has_key_frames(tfx_effect_descriptor e);

//--------------------------------
//Grouped graph struct functions
//--------------------------------
tfxINTERNAL inline tfx_graph_t *tfx__get_descriptor_graph(tfx_effect_descriptor effect, tfxU32 graph_index) {
	return &effect->library->graphs[effect->graph_list_index].graphs[graph_index];
}
tfxAPI_EDITOR void tfx__initialise_path(tfx_emitter_path_t *path);
tfxAPI_EDITOR void tfx__space_path_nodes_evenly(tfx_emitter_path_t *path, int range_start = -1, int range_end = -1);
tfxAPI_EDITOR void tfx__build_path_nodes(tfx_emitter_path_t *path);
tfxAPI_EDITOR tfxU32 tfx__add_emitter_path_attributes(tfx_library library);
tfxAPI_EDITOR tfx_emitter_path_t *tfx__get_path(tfx_effect_descriptor descriptor);
tfxAPI_EDITOR void tfx__copy_path(tfx_emitter_path_t *src, const char *name, tfx_emitter_path_t *emitter_path);
tfxAPI_EDITOR bool tfx__has_translation_key_frames(tfx_graph_list_t *graphs);
tfxAPI_EDITOR tfx_graph_t *tfx__get_graph(tfx_library library, tfx_graph_id_t graph_id);
tfxAPI_EDITOR tfx_effect_library_stats_t tfx__create_library_stats(tfx_library lib);
tfxAPI_EDITOR bool tfx__is_library_shape_used(tfx_library library, tfxKey image_hash);
tfxAPI_EDITOR bool tfx__library_shape_exists(tfx_library library, tfxKey image_hash);
tfxAPI_EDITOR bool tfx__remove_library_shape(tfx_library library, tfxKey image_hash);
tfxAPI_EDITOR tfx_effect_descriptor tfx__insert_library_effect(tfx_library library, tfx_effect_descriptor effect, tfx_effect_descriptor position);
tfxAPI_EDITOR tfx_effect_descriptor tfx__add_library_effect(tfx_library library, tfx_effect_descriptor effect);
tfxAPI_EDITOR tfx_effect_descriptor tfx__add_new_library_folder(tfx_library library, tfx_str64_t *name);
tfxAPI_EDITOR tfx_effect_descriptor tfx__add_new_library_effect(tfx_library library, tfx_str64_t *name);
tfxAPI_EDITOR tfx_effect_descriptor tfx__add_library_stage(tfx_library library, tfx_str64_t *name);
tfxAPI_EDITOR void tfx__update_library_effect_paths(tfx_library library);
tfxAPI_EDITOR bool tfx__rename_library_effect(tfx_library library, tfx_effect_descriptor effect, const char *new_name);
tfxAPI_EDITOR bool tfx__library_name_exists(tfx_library library, tfx_effect_descriptor effect, const char *name);
tfxAPI_EDITOR void tfx__reindex_library(tfx_library library);
tfxAPI_EDITOR void tfx__update_library_particle_shape_references(tfx_library library, tfxKey default_hash);
tfxAPI_EDITOR tfx_effect_descriptor tfx__library_move_up(tfx_library library, tfx_effect_descriptor effect);
tfxAPI_EDITOR tfx_effect_descriptor tfx__library_move_down(tfx_library library, tfx_effect_descriptor effect);
tfxAPI_EDITOR void tfx__add_library_effect_graphs(tfx_library library, tfx_effect_descriptor effect);
tfxAPI_EDITOR tfxU32 tfx__add_library_transform_graphs(tfx_library library);
tfxAPI_EDITOR tfxU32 tfx__add_library_sprite_sheet_settings(tfx_library library, tfx_effect_descriptor effect);
tfxAPI_EDITOR tfxU32 tfx__add_library_sprite_data_settings(tfx_library library, tfx_effect_descriptor effect);
tfxAPI_EDITOR tfxU32 tfx__add_library_preview_camera_settings_effect(tfx_library library, tfx_effect_descriptor effect);
tfxAPI_EDITOR tfxU32 tfx__allocate_library_preview_camera_settings(tfx_library library);
tfxAPI_EDITOR tfxU32 tfx__allocate_library_particle_emitter_properties(tfx_library library);
tfxAPI_EDITOR tfxU32 tfx__allocate_library_shared_properties(tfx_library library);
tfxAPI_EDITOR tfxU32 tfx__allocate_library_ribbon_emitter_properties(tfx_library library);
tfxAPI_EDITOR void tfx__update_library_compute_nodes();
tfxAPI_EDITOR void tfx__update_library_emitter_compute_nodes(tfx_effect_descriptor_t *emitter);
tfxAPI_EDITOR void tfx__update_all_library_graphs(tfx_library library);
tfxAPI_EDITOR bool tfx__update_library_color_graphs(tfx_library library, tfxU32 index);
tfxAPI_EDITOR void tfx__init_library(tfx_library library);
tfxAPI_EDITOR bool tfx__is_valid_effect_path(tfx_library library, const char *path);
tfxAPI_EDITOR bool tfx__is_valid_effect_key(tfx_library library, tfxKey key);
tfxAPI_EDITOR tfx_effect_descriptor tfx__get_library_effect_by_key(tfx_library library, tfxKey key);
tfxAPI_EDITOR void tfx__record_sprite_data(tfx_effect_manager pm, tfx_effect_descriptor effect, tfx_sprite_data_settings_t *settings, tfx_sprite_data_t *sprite_data, float update_frequency, float camera_position[3], int *progress);
tfxAPI_EDITOR tfxU32 tfx__add_library_graphs(tfx_library library, tfx_effect_descriptor_type type);
tfxAPI_EDITOR void tfx__copy_graph_list(tfx_graph_list_t *src, tfx_graph_list_t *dst);
tfxAPI_EDITOR tfxU32 tfx__clone_library_particle_emitter_properties(tfx_library library, tfxU32 source_index, tfx_library destination_library);
tfxAPI_EDITOR tfxU32 tfx__clone_library_ribbon_emitter_properties(tfx_library library, tfxU32 source_index, tfx_library destination_library);
tfxAPI_EDITOR tfxU32 tfx__clone_library_shared_properties(tfx_library library, tfxU32 source_index, tfx_library destination_library);
tfxINTERNAL void tfx__init_graph_list(tfx_graph_list_t *graph_list);
tfxINTERNAL void tfx__copy_graph_list_range_no_lookups(tfx_graph_list_t *src, tfx_graph_list_t *dst, tfxU32 from_index, tfxU32 to_index);
tfxINTERNAL void tfx__copy_graph_list_range(tfx_graph_list_t *src, tfx_graph_list_t *dst, tfxU32 from_index, tfxU32 to_index);
tfxINTERNAL int tfx__get_effect_library_stats(const char *filename, tfx_effect_library_stats_t *stats);
tfxINTERNAL void tfx__toggle_sprites_with_uid(tfx_effect_manager pm, bool switch_on);
tfxINTERNAL void tfx__add_library_path(tfx_library library, tfx_effect_descriptor effect_emitter, const char *path, bool skip_existing);
tfxINTERNAL void tfx__free_library_graphs(tfx_graph_list_t *graph_list);
tfxINTERNAL void tfx__free_library_graph_list(tfx_library library, tfxU32 index);
tfxINTERNAL void tfx__free_library_properties(tfx_effect_descriptor descriptor);
tfxINTERNAL void tfx__free_library_emitter_properties(tfx_library library, tfxU32 index);
tfxINTERNAL void tfx__free_library_ribbon_properties(tfx_library library, tfxU32 index);
tfxINTERNAL void tfx__free_library_shared_properties(tfx_library library, tfxU32 index);
tfxINTERNAL tfxU32 tfx__clone_library_graph_list(tfx_library library, tfxU32 source_index, tfx_library destination_library);
tfxINTERNAL tfxU32 tfx__clone_library_transform_graph_list(tfx_library library, tfxU32 source_index, tfx_library destination_library);
tfxINTERNAL tfx_str256_t tfx__find_new_path_name(tfx_library library, const char *path);

//Effect/Emitter functions
tfxAPI_EDITOR void tfx__set_effect_user_data(tfx_effect_descriptor e, void *data);
tfxAPI_EDITOR void *tfx__get_effect_user_data(tfx_effect_descriptor e);
tfxAPI_EDITOR tfx_particle_emitter_properties_t *tfx__get_particle_emitter_properties(tfx_effect_descriptor e);
tfxAPI_EDITOR tfx_shared_properties_t *tfx__get_shared_emitter_properties(tfx_effect_descriptor e);
tfxAPI_EDITOR tfx_ribbon_emitter_properties_t *tfx__get_ribbon_emitter_properties(tfx_effect_descriptor e);
tfxAPI_EDITOR tfx_effect_descriptor tfx__add_emitter_to_effect(tfx_effect_descriptor effect, tfx_effect_descriptor e, tfx_effect_descriptor_type type);
tfxAPI_EDITOR tfx_effect_descriptor tfx__add_new_ribbon_to_effect(tfx_effect_descriptor effect, tfx_str64_t *name);
tfxAPI_EDITOR tfx_effect_descriptor tfx__add_effect_to_emitter(tfx_effect_descriptor effect, tfx_effect_descriptor e);
tfxAPI_EDITOR int tfx__get_effect_depth(tfx_effect_descriptor e);
tfxAPI_EDITOR tfxU32 tfx__count_all_effects(tfx_effect_descriptor effect, tfxU32 amount);
tfxAPI_EDITOR tfx_effect_descriptor tfx__get_root_effect(tfx_effect_descriptor effect);
tfxAPI_EDITOR void tfx__reindex_effect(tfx_effect_descriptor effect);
tfxAPI_EDITOR tfx_effect_descriptor tfx__move_effect_up(tfx_effect_descriptor effect_to_move);
tfxAPI_EDITOR tfx_effect_descriptor tfx__move_effect_down(tfx_effect_descriptor effect_to_move);
tfxAPI_EDITOR void tfx__free_effect(tfx_effect_descriptor effect);
tfxAPI_EDITOR void tfx__clear_effect(tfx_effect_descriptor effect);
tfxAPI_EDITOR void tfx__reset_effect_graphs(tfx_effect_descriptor effect, bool add_node = true);
tfxAPI_EDITOR void tfx__reset_transform_graphs(tfx_effect_descriptor effect, bool add_node = true);
tfxAPI_EDITOR void tfx__reset_emitter_graphs(tfx_effect_descriptor effect, bool add_node = true);
tfxAPI_EDITOR void tfx__reset_ribbon_graphs(tfx_effect_descriptor effect, bool add_node = true);
tfxAPI_EDITOR void tfx__add_emitter_color_overtime(tfx_effect_descriptor effect, float frame, tfx_rgb_t color);
tfxAPI_EDITOR void tfx__update_emitter_max_life(tfx_effect_descriptor effect);
tfxAPI_EDITOR tfx_graph_t *tfx__get_effect_graph_by_index(tfx_effect_descriptor effect, tfxU32 index);
tfxAPI_EDITOR tfx_graph_t *tfx__get_effect_transform_graph_by_index(tfx_effect_descriptor effect, tfxU32 index);
tfxAPI_EDITOR tfxU32 tfx__get_effect_graph_index_by_type(tfx_effect_descriptor effect, tfx_graph_type type);
tfxAPI_EDITOR tfx_graph_list_t *tfx__get_descriptor_graph_list(tfx_effect_descriptor emitter);
tfxAPI_EDITOR tfx_graph_list_t *tfx__get_library_graph_list(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void tfx__initialise_unitialised_graphs(tfx_effect_descriptor effect);
tfxAPI_EDITOR void tfx__set_effect_name(tfx_effect_descriptor effect, const char *n);
tfxAPI_EDITOR bool tfx__rename_child(tfx_effect_descriptor effect, const char *new_name);
tfxAPI_EDITOR bool tfx__effect_name_exists(tfx_effect_descriptor in_effect, tfx_effect_descriptor excluding_effect, const char *name);
tfxAPI_EDITOR tfx_effect_descriptor tfx__clone_effect_into_library(tfx_effect_descriptor effect_to_clone, tfx_effect_descriptor root_parent, tfx_library destination_library, tfxEffectCloningFlags flags = 0);
tfxAPI_EDITOR void tfx__overwrite_effect(tfx_effect_descriptor src, tfx_effect_descriptor *dst);
tfxAPI_EDITOR void tfx__enable_all_emitters(tfx_effect_descriptor effect);
tfxAPI_EDITOR void tfx__render_all_emitters(tfx_effect_descriptor effect);
tfxAPI_EDITOR void tfx__do_not_render_all_emitters_except(tfx_effect_descriptor effect, tfx_effect_descriptor emitter);
tfxAPI_EDITOR void tfx__disable_all_emitters(tfx_effect_descriptor effect);
tfxAPI_EDITOR void tfx__disable_all_emitters_except(tfx_effect_descriptor effect, tfx_effect_descriptor emitter);
tfxAPI_EDITOR void tfx__hide_descriptor(tfx_effect_descriptor descriptor);
tfxAPI_EDITOR void tfx__show_descriptor(tfx_effect_descriptor descriptor);
tfxAPI_EDITOR bool tfx__is_descriptor_hidden(tfx_effect_descriptor descriptor);
tfxAPI_EDITOR bool tfx__is_finite_effect(tfx_effect_descriptor effect);
tfxAPI_EDITOR bool tfx__is_finite_emitter(tfx_effect_descriptor emitter);
tfxAPI_EDITOR bool tfx__is_emitter_type(tfx_effect_descriptor emitter);
tfxAPI_EDITOR bool tfx__is_ordered_effect(tfx_effect_descriptor effect);
tfxAPI_EDITOR bool tfx__has_emission_range(tfx_effect_descriptor emitter);
tfxAPI_EDITOR tfx_preview_camera_settings_t *tfx__effect_camera_settings(tfx_effect_descriptor effect);
tfxAPI_EDITOR float tfx__get_effect_highest_loop_length(tfx_effect_descriptor effect);
tfxAPI_EDITOR void tfx__update_source_emitter_flags(tfx_effect_descriptor effect);
tfxINTERNAL void tfx__clone_effect(tfx_effect_descriptor effect_to_clone, tfx_effect_descriptor clone, tfx_library destination_library, tfxEffectCloningFlags flags = 0);
tfxINTERNAL void tfx__swap_depth_index(tfx_depth_index_t *left, tfx_depth_index_t *right);
tfxINTERNAL tfx_effect_descriptor tfx__add_effect(tfx_effect_descriptor effect);
tfxINTERNAL void tfx__enable_emitter(tfx_effect_descriptor effect);
tfxINTERNAL bool tfx__is_ordered_effect_state(tfx_effect_state_t *effect);
tfxINTERNAL void tfx__assign_stage_property_u32(tfx_effect_descriptor effect, tfx_str256_t *field, tfxU32 value);
tfxINTERNAL void tfx__assign_stage_property_float(tfx_effect_descriptor effect, tfx_str256_t *field, float value);
tfxINTERNAL void tfx__assign_stage_property_bool(tfx_effect_descriptor effect, tfx_str256_t *field, bool value);
tfxINTERNAL void tfx__assign_stage_property_int(tfx_effect_descriptor effect, tfx_str256_t *field, int value);
tfxINTERNAL void tfx__assign_stage_property_str(tfx_effect_descriptor effect, tfx_str256_t *field, tfx_str256_t *value);
tfxINTERNAL void tfx__assign_sprite_data_metrics_property_u32(tfx_sprite_data_metrics_t *metrics, tfx_str256_t *field, tfxU32 value, tfxU32 file_version);
tfxINTERNAL void tfx__assign_sprite_data_metrics_property_u64(tfx_sprite_data_metrics_t *metrics, tfx_str256_t *field, tfxU64 value, tfxU32 file_version);
tfxINTERNAL void tfx__assign_sprite_data_metrics_property_float(tfx_sprite_data_metrics_t *metrics, tfx_str256_t *field, float value, tfxU32 file_version);
tfxINTERNAL void tfx__assign_sprite_data_metrics_property_str(tfx_sprite_data_metrics_t *metrics, tfx_str256_t *field, const char *value, tfxU32 file_version);
tfxINTERNAL void tfx__assign_frame_meta_property_u32(tfx_frame_meta_t *metrics, tfx_str256_t *field, tfxU32 value, tfxU32 file_version);
tfxINTERNAL void tfx__assign_frame_meta_property_vec3(tfx_frame_meta_t *metrics, tfx_str256_t *field, tfx_vec3_t value, tfxU32 file_version);
tfxINTERNAL void tfx__assign_animation_emitter_property_u32(tfx_animation_emitter_properties_t *properties, tfx_str256_t *field, tfxU32 value, tfxU32 file_version);
tfxINTERNAL void tfx__assign_animation_emitter_property_float(tfx_animation_emitter_properties_t *properties, tfx_str256_t *field, float value, tfxU32 file_version);
tfxINTERNAL void tfx__assign_animation_emitter_property_vec2(tfx_animation_emitter_properties_t *properties, tfx_str256_t *field, tfx_vec2_t value, tfxU32 file_version);
tfxINTERNAL void tfx__assign_node_data(tfx_attribute_node_t *node, tfx_vector_t<tfx_str256_t> *values);
tfxINTERNAL tfx_vec3_t tfx__str_to_vec3(tfx_vector_t<tfx_str256_t> *str);
tfxINTERNAL tfx_vec2_t tfx__str_to_vec2(tfx_vector_t<tfx_str256_t> *str);

//Debug output functions
tfxAPI_EDITOR void tfx__print_effect(tfx_effect_descriptor effect);
tfxAPI_EDITOR tfx_str32_t tfx__descriptor_type_to_string(tfx_effect_descriptor_type type);
tfxAPI_EDITOR tfx_str32_t tfx__graph_sampling_type_to_string(tfx_graph_easing_type type);

//--------------------------------
//Math functions
//--------------------------------
tfxAPI_EDITOR float tfx__quake_sqrt(float number);
tfxINTERNAL void tfx__make_icospheres();
tfxINTERNAL int tfx__vertex_for_edge(tfx_storage_map_t<int> *point_cache, tfx_vector_t<tfx_vec3_t> *vertices, int first, int second);
tfxINTERNAL void tfx__sub_divide_icosphere(tfx_storage_map_t<int> *point_cache, tfx_vector_t<tfx_vec3_t> *vertices, tfx_vector_t<tfx_face_t> *triangles, tfx_vector_t<tfx_face_t> *result);
tfxINTERNAL int tfx__sort_icosphere_points(void const *left, void const *right);
tfxINTERNAL int tfx__sort_depth(void const *left, void const *right);
tfxINTERNAL void tfx__insertion_sort_depth(tfx_work_queue_t *queue, void *work_entry);
tfxINTERNAL float tfx__length_vec4_nosqr(tfx_vec4_t const *v);
tfxINTERNAL float tfx__length_vec4(tfx_vec4_t const *v);
tfxINTERNAL float tfx__has_length_vec3(tfx_vec3_t const *v);
tfxAPI_EDITOR tfx_vec3_t tfx__normalize_vec3(tfx_vec3_t const *v);
tfxINTERNAL tfx_vec4_t tfx__normalize_vec4(tfx_vec4_t const *v);
tfxINTERNAL tfx_vec3_t tfx__cross_product_vec3(tfx_vec3_t a, tfx_vec3_t b);
tfxINTERNAL float tfx__dot_product_vec4(const tfx_vec4_t *a, const tfx_vec4_t *b);
tfxINTERNAL float tfx__dot_product_vec3(const tfx_vec3_t *a, const tfx_vec3_t *b) { return (a->x * b->x + a->y * b->y + a->z * b->z); }
tfxINTERNAL float tfx__catmull_rom_segment(tfx_vector_t<tfx_vec4_t> *nodes, float length);
tfxINTERNAL inline tfx_vec3_t tfx__normalize_vec3_fast(tfx_vec3_t const *v) { return *v * tfx__quake_sqrt(tfx__dot_product_vec3(v, v)); }
//Quake 3 inverse square root
tfxINTERNAL tfx_mat4_t tfx__create_matrix4(float v);
tfxINTERNAL tfx_mat4_t tfx__matrix4_rotate_x(float angle);
tfxINTERNAL tfx_mat4_t tfx__matrix4_rotate_y(float angle);
tfxINTERNAL tfx_mat4_t tfx__matrix4_rotate_z(float angle);
tfxINTERNAL tfx_mat4_t tfx__transform_matrix4(const tfx_mat4_t *in, const tfx_mat4_t *m);
tfxINTERNAL tfx_vec4_t tfx__transform_matrix4_vec4(const tfx_mat4_t *mat, const tfx_vec4_t vec);
tfxINTERNAL tfxU32 tfx__pack10bit_unsigned(tfx_vec3_t const *v);
tfxINTERNAL tfxWideFloat tfx__wide_unpack10bit_y(tfxWideInt in);
tfxINTERNAL float tfx__gamma_correct(float color, float gamma = tfxGAMMA);
tfxINTERNAL inline tfx_vec2_t tfx__normalise_vec2(tfx_vec2_t v) { return v * tfx__quake_sqrt(tfx__dot_product_vec2(&v, &v)); }
tfxINTERNAL tfx_vec2_t tfx__interpolate_vec2(float tween, tfx_vec2_t from, tfx_vec2_t to);
tfxINTERNAL tfx_vec3_t tfx__interpolate_vec3(float tween, tfx_vec3_t from, tfx_vec3_t to);
tfxINTERNAL tfx_rgba8_t tfx__interpolate_rgba8(float tween, tfx_rgba8_t from, tfx_rgba8_t to);
tfxINTERNAL tfxWideFloat tfx__wide_interpolate(tfxWideFloat tween, tfxWideFloat *from, tfxWideFloat *to);
tfxINTERNAL float tfx__interpolate_float(float tween, float from, float to);

tfxINTERNAL inline tfxWideInt tfx__wide_seedgen_base(tfxWideInt base, tfxWideInt h)
{
	h = tfxWideXOri(base, h);
	tfxWideInt temp = tfxWideShiftRight(h, 15);
	h = tfxWideXOri(h, temp);
	tfxWideInt multiplier1 = tfxWideSetSinglei(0x85ebca6b);
	h = tfxWideMuli(h, multiplier1);
	temp = tfxWideShiftRight(h, 13);
	h = tfxWideXOri(h, temp);
	tfxWideInt multiplier2 = tfxWideSetSinglei(0xc2b2ae35);
	h = tfxWideMuli(h, multiplier2);
	temp = tfxWideShiftRight(h, 16);
	h = tfxWideXOri(h, temp);
	return h;
}

tfxINTERNAL inline tfxWideFloat tfx__wide_seedgen(tfxWideInt h)
{
	h = tfxWideXOri(h, tfxWideShiftRight(h, 15));
	tfxWideInt mult1 = tfxWideSetSinglei(0x85ebca6b);
	h = tfxWideMuli(h, mult1);
	h = tfxWideXOri(h, tfxWideShiftRight(h, 13));
	tfxWideInt mult2 = tfxWideSetSinglei(0xc2b2ae35);
	h = tfxWideMuli(h, mult2);
	h = tfxWideXOri(h, tfxWideShiftRight(h, 16));
	return tfxWideConvert(h);
}

//--------------------------------
//Control particle inline functions and policies
//--------------------------------
#define tfx__wide_dp_xyz(x1,y1,z1,x2,y2,z2) tfxWideMulAdd(x1, x2, tfxWideMulAdd(y1, y2, tfxWideMul(z1, z2)))
#define tfx__wide_dp_xy(x1,y1,x2,y2) tfxWideMulAdd(x1, x2, tfxWideMul(y1, y2))

tfxINTERNAL inline void tfx__wide_unpack10bit(tfxWideInt in, tfxWideFloat &x, tfxWideFloat &y, tfxWideFloat &z) {
	const tfxWideInt mask_x = tfxWideSetSinglei(0x3FF00000);
	const tfxWideInt mask_y = tfxWideSetSinglei(0x000FFC00);
	const tfxWideInt mask_z = tfxWideSetSinglei(0x000003FF);
	x = tfxWideConvert(tfxWideShiftRight(tfxWideAndi(in, mask_x), 20));
	y = tfxWideConvert(tfxWideShiftRight(tfxWideAndi(in, mask_y), 10));
	z = tfxWideConvert(tfxWideAndi(in, mask_z));
	x = tfxWideMulAdd(x, one_div_511_wide.m, tfxWIDEMINUSONE.m);
	y = tfxWideMulAdd(y, one_div_511_wide.m, tfxWIDEMINUSONE.m);
	z = tfxWideMulAdd(z, one_div_511_wide.m, tfxWIDEMINUSONE.m);
}

tfxINTERNAL inline void tfx__wide_unpack16bit(tfxWideInt xy, tfxWideInt zw, tfxWideFloat &x, tfxWideFloat &y, tfxWideFloat &z, tfxWideFloat &w) {
	const tfxWideInt mask_ffff = tfxWideSetSinglei(0xFFFF);
	const tfxWideInt sign_bit = tfxWideSetSinglei(0x8000);

	// Extract and sign-extend x (low 16 bits of xy)
	tfxWideInt x_val = tfxWideAndi(xy, mask_ffff);
	tfxWideInt x_sign = tfxWideShiftRight(tfxWideAndi(x_val, sign_bit), 15);
	x_val = tfxWideSubi(x_val, tfxWideShiftLeft(x_sign, 16));

	// Extract and sign-extend y (high 16 bits of xy)
	tfxWideInt y_val = tfxWideShiftRight(xy, 16);
	tfxWideInt y_sign = tfxWideShiftRight(tfxWideAndi(y_val, sign_bit), 15);
	y_val = tfxWideSubi(y_val, tfxWideShiftLeft(y_sign, 16));

	// Extract and sign-extend z (low 16 bits of zw)
	tfxWideInt z_val = tfxWideAndi(zw, mask_ffff);
	tfxWideInt z_sign = tfxWideShiftRight(tfxWideAndi(z_val, sign_bit), 15);
	z_val = tfxWideSubi(z_val, tfxWideShiftLeft(z_sign, 16));

	// Extract and sign-extend w (high 16 bits of zw)
	tfxWideInt w_val = tfxWideShiftRight(zw, 16);
	tfxWideInt w_sign = tfxWideShiftRight(tfxWideAndi(w_val, sign_bit), 15);
	w_val = tfxWideSubi(w_val, tfxWideShiftLeft(w_sign, 16));

	x = tfxWideMul(tfxWideConvert(x_val), one_div_32767_wide.m);
	y = tfxWideMul(tfxWideConvert(y_val), one_div_32767_wide.m);
	z = tfxWideMul(tfxWideConvert(z_val), one_div_32767_wide.m);
	w = tfxWideMul(tfxWideConvert(w_val), one_div_32767_wide.m);
}

tfxINTERNAL inline tfxWideInt tfx__wide_pack10bit_unsigned(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z) {
	const tfxWideFloat w511 = tfxWideSetSingle(511.f);
	const tfxWideInt bits10 = tfxWideSetSinglei(0x3FF);
	tfxWideInt converted_x = tfxWideConverti(tfxWideMulAdd(v_x, w511, w511));
	converted_x = tfxWideAndi(converted_x, bits10);
	converted_x = tfxWideShiftLeft(converted_x, 20);
	tfxWideInt converted_y = tfxWideConverti(tfxWideMulAdd(v_y, w511, w511));
	converted_y = tfxWideAndi(converted_y, bits10);
	converted_y = tfxWideShiftLeft(converted_y, 10);
	tfxWideInt converted_z = tfxWideConverti(tfxWideMulAdd(v_z, w511, w511));
	converted_z = tfxWideAndi(converted_z, bits10);
	return tfxWideOri(tfxWideOri(converted_x, converted_y), converted_z);
}

tfxINTERNAL inline void tfx__wide_transform_packed_quaternion_vec3(tfxWideInt *quaternion_xy, tfxWideInt *quaternion_zw, tfxWideFloat *x, tfxWideFloat *y, tfxWideFloat *z) {
	tfxWideFloat q_x, q_y, q_z, q_w;
	tfx__wide_unpack16bit(*quaternion_xy, *quaternion_zw, q_x, q_y, q_z, q_w);

	tfxWideFloat c_x = tfxWideAdd(tfxWideSub(tfxWideMul(*y, q_z), tfxWideMul(*z, q_y)), tfxWideMul(*x, q_w));
	tfxWideFloat c_y = tfxWideAdd(tfxWideSub(tfxWideMul(*z, q_x), tfxWideMul(*x, q_z)), tfxWideMul(*y, q_w));
	tfxWideFloat c_z = tfxWideAdd(tfxWideSub(tfxWideMul(*x, q_y), tfxWideMul(*y, q_x)), tfxWideMul(*z, q_w));

	*x = tfxWideAdd(tfxWideMul(tfxWideSub(tfxWideMul(c_y, q_z), tfxWideMul(c_z, q_y)), tfxWIDETWO.m), *x);
	*y = tfxWideAdd(tfxWideMul(tfxWideSub(tfxWideMul(c_z, q_x), tfxWideMul(c_x, q_z)), tfxWIDETWO.m), *y);
	*z = tfxWideAdd(tfxWideMul(tfxWideSub(tfxWideMul(c_x, q_y), tfxWideMul(c_y, q_x)), tfxWIDETWO.m), *z);
}

tfxINTERNAL inline void tfx__wide_transform_packed_quaternion_vec2(tfxWideInt *quaternion_xy, tfxWideInt *quaternion_zw, tfxWideFloat *x, tfxWideFloat *y) {
	tfxWideFloat q_x, q_y, q_z, q_w;
	tfx__wide_unpack16bit(*quaternion_xy, *quaternion_zw, q_x, q_y, q_z, q_w);

	tfxWideFloat s2 = tfxWideMul(q_z, q_z);
	tfxWideFloat c2 = tfxWideMul(q_w, q_w);
	tfxWideFloat sc = tfxWideMul(tfxWideSetSingle(2.f), tfxWideMul(q_z, q_w));

	tfxWideFloat rx = tfxWideSub(tfxWideMul(c2, *x), tfxWideAdd(tfxWideMul(sc, *y), tfxWideMul(s2, *x)));
	tfxWideFloat ry = tfxWideAdd(tfxWideMul(sc, *x), tfxWideSub(tfxWideMul(c2, *y), tfxWideMul(s2, *y)));

	*x = rx;
	*y = ry;
}

tfxINTERNAL inline tfxWideInt tfx__wide_pack16bit(tfxWideFloat v_x, tfxWideFloat v_y) {
	const tfxWideFloat w32k = tfxWideSetSingle(32767.f);
	const tfxWideInt bits16 = tfxWideSetSinglei(0xFFFF);
	tfxWideInt converted_y = tfxWideConverti(tfxWideMul(v_y, w32k));
	converted_y = tfxWideAndi(converted_y, bits16);
	converted_y = tfxWideShiftLeft(converted_y, 16);
	tfxWideInt converted_x = tfxWideConverti(tfxWideMul(v_x, w32k));
	converted_x = tfxWideAndi(converted_x, bits16);
	return tfxWideOri(converted_x, converted_y);
}

tfxINTERNAL inline tfxWideInt tfx__wide_pack16bit_2sscaled(tfxWideFloat v_x, tfxWideFloat v_y, tfxWideFloat max_value) {
	const tfxWideInt bits16 = tfxWideSetSinglei(0xFFFF);
	tfxWideInt converted_y = tfxWideConverti(tfxWideMul(v_y, max_value));
	converted_y = tfxWideAndi(converted_y, bits16);
	converted_y = tfxWideShiftLeft(converted_y, 16);
	tfxWideInt converted_x = tfxWideConverti(tfxWideMul(v_x, max_value));
	converted_x = tfxWideAndi(converted_x, bits16);
	return tfxWideOri(converted_x, converted_y);
}

tfxINTERNAL inline void tfx__wide_cross_product(tfxWideFloat ax, tfxWideFloat ay, tfxWideFloat az, tfxWideFloat *bx, tfxWideFloat *by, tfxWideFloat *bz, tfxWideFloat *rx, tfxWideFloat *ry, tfxWideFloat *rz) {
	*rx = tfxWideSub(tfxWideMul(ay, *bz), tfxWideMul(az, *by));
	*ry = tfxWideSub(tfxWideMul(az, *bx), tfxWideMul(ax, *bz));
	*rz = tfxWideSub(tfxWideMul(ax, *by), tfxWideMul(ay, *bx));
}

tfxINTERNAL inline void tfx__wide_random_vector_in_cone(tfxWideInt seed, tfxWideFloat velocity_normal_x, tfxWideFloat velocity_normal_y, tfxWideFloat velocity_normal_z, tfxWideFloat cone_angle, tfxWideFloat *random_x, tfxWideFloat *random_y, tfxWideFloat *random_z) {
	// Convert cone angle to radians

	cone_angle = tfxWideMin(cone_angle, tfx180RadiansWide.m);
	// Calculate the minimum z value for the cone
	tfxWideFloat min_z = tfxWideCos52s(cone_angle);

	tfxWideFloat max_uint = tfxWideSetSingle((float)UINT32_MAX);
	// Randomly sample z in [min_z, 1]
	tfxWideFloat z = tfxWideAdd(tfxWideDiv(tfx__wide_seedgen(seed), max_uint), tfxWideSetSingle(0.5f));
	z = tfxWideSub(tfxWIDEONE.m, tfxWideMul(tfxWideSub(tfxWIDEONE.m, min_z), z));

	// Randomly sample [0, 2pi)
	tfxWideFloat phi = tfxWideMul(tfxWideAdd(tfxWideDiv(tfx__wide_seedgen(seed), max_uint), tfxWideSetSingle(0.5f)), tfxWideSetSingle(2.f * tfxPI));

	// Calculate the corresponding x and y for the random point on the unit sphere
	tfxWideFloat sqrt_one_minus_z_squared = tfxWideSub(tfxWIDEONE.m, tfxWideMul(z, z));
	sqrt_one_minus_z_squared = tfxWideMul(tfxWideRSqrt(sqrt_one_minus_z_squared), sqrt_one_minus_z_squared);
	tfxWideFloat sin;
	tfxWideFloat cos;
	tfxWideSinCos(phi, &sin, &cos);
	tfxWideFloat x = tfxWideMul(sqrt_one_minus_z_squared, cos);
	tfxWideFloat y = tfxWideMul(sqrt_one_minus_z_squared, sin);

	// Calculate the rotation axis (cross product of (0, 0, 1) and cone_direction)
	tfx_vec3_t north_pole = { 0, 0, 1.f };
	tfxWideFloat rotation_axis_x, rotation_axis_y, rotation_axis_z;
	rotation_axis_x = tfxWideSub(tfxWideSetZero, velocity_normal_y);
	rotation_axis_y = velocity_normal_x;
	rotation_axis_z = tfxWideSetZero;
	tfxWideFloat length = tfxWideMul(rotation_axis_x, rotation_axis_x);
	length = tfxWideAdd(length, tfxWideMul(rotation_axis_y, rotation_axis_y));
	length = tfxWideMul(tfxWideRSqrt(length), length);
	rotation_axis_x = tfxWideDiv(rotation_axis_x, length);
	rotation_axis_y = tfxWideDiv(rotation_axis_y, length);

	tfxWideArray rotation_angle;
	tfxWideArray dir_z;
	dir_z.m = velocity_normal_z;
	rotation_angle.a[0] = acosf(dir_z.a[0]);
	rotation_angle.a[1] = acosf(dir_z.a[1]);
	rotation_angle.a[2] = acosf(dir_z.a[2]);
	rotation_angle.a[3] = acosf(dir_z.a[3]);
#if defined(tfxUSEAVX)
	rotation_angle.a[4] = acosf(dir_z.a[4]);
	rotation_angle.a[5] = acosf(dir_z.a[5]);
	rotation_angle.a[6] = acosf(dir_z.a[6]);
	rotation_angle.a[7] = acosf(dir_z.a[7]);
#endif

	// Rotate the random vector to align with the cone direction
	// Use Rodrigues' rotation formula
	tfxWideSinCos(rotation_angle.m, &sin, &cos);
	tfxWideFloat dot = tfxWideAdd(tfxWideMul(rotation_axis_x, x), tfxWideMul(rotation_axis_y, y));
	tfxWideFloat cx, cy, cz;
	tfx__wide_cross_product(rotation_axis_x, rotation_axis_y, rotation_axis_z, &x, &y, &z, &cx, &cy, &cz);
	//tfx_vec3_t rotated_vector = random_vector * cos + tfx__cross_product_vec3(&rotation_axis, &random_vector) * sin + rotation_axis * tfx__dot_product_vec3(&rotation_axis, &random_vector) * (1.f - cos);
	*random_x = tfxWideAdd(tfxWideAdd(tfxWideMul(x, cos), tfxWideMul(cx, sin)), tfxWideMul(tfxWideMul(rotation_axis_x, dot), tfxWideSub(tfxWIDEONE.m, cos)));
	*random_y = tfxWideAdd(tfxWideAdd(tfxWideMul(y, cos), tfxWideMul(cy, sin)), tfxWideMul(tfxWideMul(rotation_axis_y, dot), tfxWideSub(tfxWIDEONE.m, cos)));
	*random_z = tfxWideAdd(tfxWideAdd(tfxWideMul(z, cos), tfxWideMul(cz, sin)), tfxWideMul(tfxWideMul(rotation_axis_z, dot), tfxWideSub(tfxWIDEONE.m, cos)));
}

tfxINTERNAL inline tfxWideFloat tfx__wide_bezier_sampler(tfxWideFloat t, tfxWideFloat node1, tfxWideFloat curve1, tfxWideFloat curve2, tfxWideFloat node2) {
	tfxWideFloat u = tfxWideSub(tfxWIDEONE.m, t);
	tfxWideFloat w1 = tfxWideMul(tfxWideMul(u, u), u);
	tfxWideFloat w2 = tfxWideMul(tfxWideMul(tfxWideMul(tfxWIDETHREE.m, u), u), t);
	tfxWideFloat w3 = tfxWideMul(tfxWideMul(tfxWideMul(tfxWIDETHREE.m, u), t), t);
	tfxWideFloat w4 = tfxWideMul(tfxWideMul(t, t), t);
	return tfxWideAdd(tfxWideAdd(tfxWideAdd(tfxWideMul(w1, node1), tfxWideMul(w2, curve1)), tfxWideMul(w3, curve2)), tfxWideMul(w4, node2));
}

tfxINTERNAL inline tfxWideFloat tfx__wide_linear_sampler(tfxWideFloat from, tfxWideFloat to, tfxWideFloat t) {
	return tfxWideAdd(from, tfxWideMul(tfxWideSub(to, from), t));
}

tfxINTERNAL inline tfxWideInt tfx__wide_pack8bit_xyz(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z) {
	const tfxWideFloat w127 = tfxWideSetSingle(127.f);
	const tfxWideInt bits8 = tfxWideSetSinglei(0xFF);
	tfxWideInt converted_x = tfxWideConverti(tfxWideMul(v_x, w127));
	converted_x = tfxWideAndi(converted_x, bits8);
	tfxWideInt converted_y = tfxWideConverti(tfxWideMul(v_y, w127));
	converted_y = tfxWideAndi(converted_y, bits8);
	converted_y = tfxWideShiftLeft(converted_y, 8);
	tfxWideInt converted_z = tfxWideConverti(tfxWideMul(v_z, w127));
	converted_z = tfxWideAndi(converted_z, bits8);
	converted_z = tfxWideShiftLeft(converted_z, 16);
	return tfxWideOri(tfxWideOri(converted_x, converted_y), converted_z);
}

tfxINTERNAL inline tfxWideInt tfx__wide_pack16bit_unorm(tfxWideFloat v_x, tfxWideFloat v_y) {
	const tfxWideFloat w65k = tfxWideSetSingle(65535.f);
	const tfxWideInt bits16 = tfxWideSetSinglei(0xFFFF);
	tfxWideInt converted_y = tfxWideConverti(tfxWideMul(v_y, w65k));
	converted_y = tfxWideAndi(converted_y, bits16);
	converted_y = tfxWideShiftLeft(converted_y, 16);
	tfxWideInt converted_x = tfxWideConverti(tfxWideMul(v_x, w65k));
	converted_x = tfxWideAndi(converted_x, bits16);
	return tfxWideOri(converted_x, converted_y);
}

tfxINTERNAL inline float tfx__length_vec3_nosqr(tfx_vec3_t const *v) {
	return v->x * v->x + v->y * v->y + v->z * v->z;
}

typedef union 
{
	unsigned int i;
	float f;
} tfx_uint_float;

//Exponent lookup table for converting floats to halfs
static const unsigned short tfx_exponent_lookup[1 << 9] = {
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,  1024,  2048,  3072,  4096,  5120,  6144,  7168,
	 8192,  9216, 10240, 11264, 12288, 13312, 14336, 15360,
	16384, 17408, 18432, 19456, 20480, 21504, 22528, 23552,
	24576, 25600, 26624, 27648, 28672, 29696,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0, 33792, 34816, 35840, 36864, 37888, 38912, 39936,
	40960, 41984, 43008, 44032, 45056, 46080, 47104, 48128,
	49152, 50176, 51200, 52224, 53248, 54272, 55296, 56320,
	57344, 58368, 59392, 60416, 61440, 62464,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
		0,     0,     0,     0,     0,     0,     0,     0,
};

tfxAPI inline unsigned short tfx__float_to_half(float f) {
	tfx_uint_float float_int_union;

	//Use memcpy if this causes issues at some point
	float_int_union.f = f;

	if (f == 0)
	{
		//We don't care about saving the sign, just return 0;
		return 0;
	} else {
		int exponent = (float_int_union.i >> 23) & 0x000001ff;
		int e = tfx_exponent_lookup[exponent];
		if (e)
		{
			int mantissa = float_int_union.i & 0x007fffff;
			return (unsigned short)(e + ((mantissa + 0x00000fff + ((mantissa >> 13) & 1)) >> 13));
		} else if (exponent == 0xFF || exponent == 0x1FF) {
			//inf/nan
			return 0;
		} else if (exponent > 142) {
			//Overflow, return max half value 
			return 0x7BFF;
		} else {
			//all other cases
			return 0;
		}
	}
}

/*
tfxINTERNAL inline tfxWideInt tfx__load_half_ints(tfxU16 *half_ints) {
#ifdef tfxUSEAVX
	tfx128i loaded_16_bytes = _mm_loadu_si128((tfx128i*)half_ints);
	tfx128i dwords_0_3 = _mm_cvtepu16_epi32(loaded_16_bytes); 
	tfx128i packed_words_4_7_shifted = _mm_srli_si128(loaded_16_bytes, 8);
	tfx128i dwords_4_7 = _mm_cvtepu16_epi32(packed_words_4_7_shifted); 
	return tfxWideSet128i(dwords_4_7, dwords_0_3); 
#else
	tfx128i loaded_16_bytes = _mm_loadu_si64((tfx128i*)half_ints);
	return _mm_cvtepu16_epi32(loaded_16_bytes); 
#endif
}

tfxINTERNAL inline void tfx__store_half_ints(tfxU16 *dst, tfxWideInt wide_ints) {
#ifdef tfxUSEAVX
	tfx128i lo = tfxWideExtract128i(wide_ints, 0);
	tfx128i hi = tfxWideExtract128i(wide_ints, 1);
	tfx128i zero = tfx128SetZeroi;
	lo = tfx128Packus32(lo, zero);
	hi = tfx128Packus32(hi, zero);
	tfx128i packed = _mm_unpacklo_epi64(lo, hi);
	*(tfx128i*)dst = packed;
#else
*/
	/*
	wide ints contains 4 16bit ints:
	[ 0000 FFFF, 0000 FFFF, 0000 FFFF, 0000 FFFF ]
	Pack them with _mm_packus_epi32 into the lower 64 bits of the register to end up with
	[ 0000, 0000, 0000, 0000, FFFF, FFFF, FFFF, FFFF ]
	Can then use _mm_storel_epi64 to store the lower 64 bits.
	*/
/*
	tfxWideInt zero = tfxWideSetZeroi;						// SSE2
	tfxWideInt packed = _mm_packus_epi32(wide_ints, zero);	// SSE4.1

	_mm_storel_epi64((tfxWideInt*)dst, packed); // SSE2
#endif
}
*/

//Pack a vec3 into a u16 using octahedral mapping
tfxAPI_EDITOR inline tfxU16 tfx__pack_octahedral_vec3(tfx_vec3_t v) {
	float l1_norm = fabsf(v.x) + fabsf(v.y) + fabsf(v.z);
	if (l1_norm < 1e-6f) {
		return 0;
	}
	float px = v.x / l1_norm;
	float py = v.y / l1_norm;
	tfx_vec2_t p;
	if (v.z >= 0.0f) {
		p.x = px;
		p.y = py;
	} else {
		float temp_x = copysignf(1.f, px);
		float temp_y = copysignf(1.f, py);
		p.x = (1.f - fabsf(py)) * copysignf(1.f, px);
		p.y = (1.f - fabsf(px)) * copysignf(1.f, py);
	}
	p.x = fmaf(p.x, .5f, .5f);
	p.y = fmaf(p.y, .5f, .5f);
	return (tfxU16(p.x * 255.f) << 8) | tfxU16(p.y * 255.f);
}

tfxAPI_EDITOR inline tfx_vec3_t tfx__unpack_octahedral_vec3(tfxU16 uv) {
	tfx_vec3_t v;
	tfx_vec2_t packed_vec2 = { float(uv >> 8) / 255.f, float(uv & 0x00FF) / 255.f };
	packed_vec2.x = fmaf(packed_vec2.x, 2.f, -1.f);
	packed_vec2.y = fmaf(packed_vec2.y, 2.f, -1.f);
	v.x = packed_vec2.x;
	v.y = packed_vec2.y;
	v.z = 1.f - (fabsf(packed_vec2.x) + fabsf(packed_vec2.y));
	if (v.z < 0.f) {
		float temp_x = v.x;
		v.x = (1.f - fabsf(packed_vec2.y)) * copysignf(1.f, temp_x);
		v.y = (1.f - fabsf(packed_vec2.x)) * copysignf(1.f, packed_vec2.y);
	}
	return tfx__normalize_vec3_fast(&v);
}

//Pack a vec3 into a u16 using octahedral mapping
tfxAPI_EDITOR inline tfxWideInt tfx__wide_pack_octahedral_vec3(tfxWideFloat &v_x, tfxWideFloat &v_y, tfxWideFloat &v_z) {
	tfxWideFloat l1_norm = tfxWideDiv(tfxWIDEONE.m, tfxWideAdd(tfxWideAdd(tfxWideAbs(v_x), tfxWideAbs(v_y)), tfxWideAbs(v_z)));
	tfxWideFloat px = tfxWideMul(v_x, l1_norm);
	tfxWideFloat py = tfxWideMul(v_y, l1_norm);
	tfxWideFloat v_z_mask_ge_0 = tfxWideGreaterEqual(v_z, tfxWideSetZero);
	tfxWideFloat final_x = tfxWideMul(tfxWideSub(tfxWIDEONE.m, tfxWideAbs(py)), tfxWideCopySign(tfxWIDEONE.m, px));
	tfxWideFloat final_y = tfxWideMul(tfxWideSub(tfxWIDEONE.m, tfxWideAbs(px)), tfxWideCopySign(tfxWIDEONE.m, py));
	final_x = tfxWideBlendv(final_x, px, v_z_mask_ge_0);
	final_y = tfxWideBlendv(final_y, py, v_z_mask_ge_0);
	final_x = tfxWideMulAdd(final_x, tfxWIDEHALF.m, tfxWIDEHALF.m);
	final_y = tfxWideMulAdd(final_y, tfxWIDEHALF.m, tfxWIDEHALF.m);
	tfxWideInt packed_x = tfxWideShiftLeft(tfxWideConverti(tfxWideMul(final_x, tfxWIDE255.m)), 8);
	tfxWideInt packed_y = tfxWideConverti(tfxWideMul(final_y, tfxWIDE255.m));
	return tfxWideOri(packed_x, packed_y);
}

tfxAPI_EDITOR inline void tfx__wide_unpack_octahedral_vec3(tfxWideInt uv, tfxWideFloat &u_x, tfxWideFloat &u_y, tfxWideFloat &u_z) {
	const tfxWideInt y_mask = tfxWideSetSinglei(0xFF);
	tfxWideFloat packed_x = tfxWideMul(tfxWideConvert(tfxWideShiftRight(uv, 8)), tfxWIDE255r.m);
	tfxWideFloat packed_y = tfxWideMul(tfxWideConvert(tfxWideAndi(y_mask, uv)), tfxWIDE255r.m);
	packed_x = tfxWideMulAdd(packed_x, tfxWIDETWO.m, tfxWIDEMINUSONE.m);
	packed_y = tfxWideMulAdd(packed_y, tfxWIDETWO.m, tfxWIDEMINUSONE.m);
	u_x = packed_x;
	u_y = packed_y;
	u_z = tfxWideSub(tfxWIDEONE.m, tfxWideAdd(tfxWideAbs(packed_x), tfxWideAbs(packed_y)));
	tfxWideFloat u_z_mask_less_0 = tfxWideLess(u_z, tfxWideSetZero);
	tfxWideFloat temp_x = u_x;
	u_x = tfxWideBlendv(u_x, tfxWideMul(tfxWideSub(tfxWIDEONE.m, tfxWideAbs(packed_y)), tfxWideCopySign(tfxWIDEONE.m, temp_x)), u_z_mask_less_0);
	u_y = tfxWideBlendv(u_y, tfxWideMul(tfxWideSub(tfxWIDEONE.m, tfxWideAbs(packed_x)), tfxWideCopySign(tfxWIDEONE.m, packed_y)), u_z_mask_less_0);
	tfxWideFloat l = tfxWideMulAdd(u_x, u_x, tfxWideMulAdd(u_y, u_y, tfxWideMul(u_z, u_z)));
	l = tfxWideRSqrt(l);
	u_x = tfxWideMul(l, u_x);
	u_y = tfxWideMul(l, u_y);
	u_z = tfxWideMul(l, u_z);
}


typedef tfxWideFloat(*tfx_wide_easing_function)(tfxWideFloat);
typedef tfxWideFloat(*tfx_wide_bezier_function)(tfxWideFloat, tfxWideFloat, tfxWideFloat, tfxWideFloat, tfxWideFloat);

//Section: Control_Position_Policies
/*
I've gone through many iterations over how to organise the functions that update a particle's position. I went from a single function which became problematic
after introducing other features like noise, paths and such and having a single function with if statements was too large and unwieldy. Then I separated out into 
different functions which repeated a lot of code, this was harder to maintain. Then I put that repeated code into macros which was just ugly and hard to debug.
Then I switched those to inline functions which was a little better but adding yet more features made the whole process a bit of a mess. Finally I settled on template
policies. I'm not a huge fan of templates but in this case they seem to work very well and allows me to avoid repeated code and let the compiler produce all of
the different variations of the control functions that I need based on the polices I use in the templated function.

There are 2 templated functions that are run, one which sets up variables outside of the particle loop, and the other applies the policies inside the loop. So for
example the templated function would be called like this:

//Setup necessary variables outside of the loop:
tfx__setup_particles_position<
	tfx_setup_vecolity_lookup_policy, 
	tfx_setup_weight_lookup_policy
>(work_entry, ctx);

//Perform the following opperations on particles inside the loop:
tfx__update_particles_position<
	tfx_apply_life_based_on_age,
	tfx_apply_lookup_velocity,
	tfx_apply_lookup_weight,
	tfx_apply_velocity,
	tfx_apply_load_position,
	tfx_apply_position
>(work_entry, ctx);

In order to pass around variables between the different policy functions this context struct is used. Not all of these variables are used in every situation
of course.
*/
struct tfx_position_policy_context {
	tfxWideArray position_x, position_y, position_z;
	tfxWideFloat velocity_x, velocity_y, velocity_z;
	tfxWideFloat weight;
	tfxWideFloat velocity;
	tfxWideFloat age;
	tfxWideFloat life;
	tfxWideFloat velocity_adjuster;
	tfxWideFloat overal_scale_wide;
	tfxWideFloat lookup_velocity;
	tfxWideFloat lookup_direction;
	tfxWideFloat lookup_weight;
	tfxWideFloat lookup_stretch;
	tfxWideFloat lookup_velocity_turbulance;
	tfxWideFloat lookup_noise_resolution;
	tfxWideFloat node_count;
	tfxWideFloat global_noise;
	tfxWideFloat global_stretch;
	tfxWideFloat emitter_width, emitter_height, emitter_depth;
	tfxWideFloat emitter_offset_x, emitter_offset_y, emitter_offset_z;
	tfxWideFloat path_position;
	tfxWideFloat path_offset;
	tfxWideFloat path_scale_variation;
	tfxWideFloat motion_randomness_base;
	tfxWideFloat emitter_handle_x, emitter_handle_y, emitter_handle_z;
	tfxWideFloat emitter_world_x, emitter_world_y, emitter_world_z;
	tfxWideFloat emitter_scale;
	tfxWideInt capture_after_transform;
	tfx_vector_align_type vector_align_type;
	tfx_emission_type emission_type;
	tfx_control_work_entry_t *work_entry;
	tfx_instance_t *sprites;
	tfxU32 running_sprite_index;
	tfxU32 start_diff;
	tfxWideInt capture_after_transform_flag;
	tfxWideInt time_step;
	tfxWideInt time_changed_mask;
	tfxWideFloat time_step_fraction;
	tfx_particle_emitter_state_t *emitter;
	tfx_emitter_path_t *path;
	tfx_graph_t *velocity_graph;
	tfx_graph_t *direction_graph;
	tfx_graph_t *weight_graph;
	tfx_graph_t *velocity_turbulance_graph;
	tfx_graph_t *noise_resolution_graph;
	tfx_graph_t *motion_randomness_graph;
	tfx_wide_easing_function velocity_easing;
	tfx_wide_easing_function direction_easing;
	tfx_wide_easing_function weight_easing;
	tfx_wide_easing_function velocity_turbulance_easing;
	tfx_wide_easing_function noise_resolution_easing;
	tfx_wide_easing_function motion_randomness_easing;
	tfxContextPolicyFlags flags;
};

//Control setup policies
//These policy structs are run outside of the loop to set up things like graph access and emitter position etc.
struct tfx_setup_vecolity_lookup_policy {
	static void apply(tfx_control_work_entry_t *work_entry, tfx_position_policy_context &ctx);
};

struct tfx_setup_weight_lookup_policy {
	static void apply(tfx_control_work_entry_t *work_entry, tfx_position_policy_context &ctx);
};

struct tfx_setup_direction_lookup_policy {
	static void apply(tfx_control_work_entry_t *work_entry, tfx_position_policy_context &ctx);
};

struct tfx_setup_simplex_lookup_policy {
	static void apply(tfx_control_work_entry_t *work_entry, tfx_position_policy_context &ctx);
};

struct tfx_setup_path_policy {
	static void apply(tfx_control_work_entry_t *work_entry, tfx_position_policy_context &ctx);
};

struct tfx_setup_line_policy {
	static void apply(tfx_control_work_entry_t *work_entry, tfx_position_policy_context &ctx);
};

struct tfx_setup_orbital_policy {
	static void apply(tfx_control_work_entry_t *work_entry, tfx_position_policy_context &ctx);
};

struct tfx_setup_motion_randomness_policy {
	static void apply(tfx_control_work_entry_t *work_entry, tfx_position_policy_context &ctx);
};

struct tfx_setup_transform_policy {
	static void apply(tfx_control_work_entry_t *work_entry, tfx_position_policy_context &ctx);
};

//Control loop policies
//These policy structs are run inside the particle loop. They should all be inlined at compile time to build the different function
//variations based on the control profile of the emitter.
struct tfx_apply_load_life {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		ctx.age = tfxWideLoad(&bank.age[index]);
		tfxWideFloat inv_max_age = tfxWideLoad(&bank.inv_max_age[index]);
		ctx.life = tfxWideMul(ctx.age, inv_max_age);
	}
};

struct tfx_apply_life_based_on_path {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		tfxWideFloat path_position = tfxWideLoad(&bank.path_position[index]);
		ctx.age = tfxWideLoad(&bank.age[index]);
		ctx.life = tfxWideDiv(path_position, ctx.node_count);
	}
};

struct tfx_apply_lookup_velocity {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		const tfxWideFloat base_velocity = tfxWideLoad(&bank.base_velocity[index]);
		tfxWideFloat velocity_time = ctx.velocity_easing(ctx.life);
		ctx.lookup_velocity = (ctx.flags & tfx_ctx_policy_flag_velocity_is_bezier_graph) ?
			tfx__wide_bezier_sampler(velocity_time, ctx.velocity_graph->wide_graph.from, ctx.velocity_graph->wide_graph.curve1, ctx.velocity_graph->wide_graph.curve2, ctx.velocity_graph->wide_graph.to) :
			tfx__wide_linear_sampler(ctx.velocity_graph->wide_graph.from, ctx.velocity_graph->wide_graph.to, velocity_time);
		if (ctx.flags & tfx_ctx_policy_flag_velocity_has_oscillator) {
			ctx.lookup_velocity = tfxWideAdd(tfxWideMul(tfxOSCILLATOR_WIDE_SIN(velocity_time, tfxWideAdd(ctx.velocity_graph->wide_oscillator.offset_x, ctx.velocity_graph->wide_oscillator.frequency), ctx.velocity_graph->wide_oscillator.amplitude), ctx.lookup_velocity), ctx.velocity_graph->wide_oscillator.offset_y);
		}
		ctx.velocity = tfxWideMul(ctx.lookup_velocity, base_velocity);
	}
};

struct tfx_apply_lookup_weight {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		const tfxWideFloat base_weight = tfxWideLoad(&bank.base_weight[index]);
		tfxWideFloat weight_time = ctx.weight_easing(ctx.life);
		ctx.lookup_weight = (ctx.flags & tfx_ctx_policy_flag_weight_is_bezier_graph) ? 
			tfx__wide_bezier_sampler(weight_time, ctx.weight_graph->wide_graph.from, ctx.weight_graph->wide_graph.curve1, ctx.weight_graph->wide_graph.curve2, ctx.weight_graph->wide_graph.to) :
			tfx__wide_linear_sampler(ctx.weight_graph->wide_graph.from, ctx.weight_graph->wide_graph.to, weight_time);
		if (ctx.flags & tfx_ctx_policy_flag_weight_has_oscillator) {
			ctx.lookup_weight = tfxWideAdd(tfxWideMul(tfxOSCILLATOR_WIDE_SIN(weight_time, tfxWideAdd(ctx.weight_graph->wide_oscillator.offset_x, ctx.weight_graph->wide_oscillator.frequency), ctx.weight_graph->wide_oscillator.amplitude), ctx.lookup_weight), ctx.weight_graph->wide_oscillator.offset_y);
		}
		ctx.weight = tfxWideMul(ctx.lookup_weight, base_weight);
	}
};

struct tfx_apply_simplex_noise {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		tfxWideFloat velocity_turbulance_time = ctx.velocity_turbulance_easing(ctx.life);
		ctx.lookup_velocity_turbulance = (ctx.flags & tfx_ctx_policy_flag_velocity_turbulance_is_bezier_graph) ?
			tfx__wide_bezier_sampler(velocity_turbulance_time, ctx.velocity_turbulance_graph->wide_graph.from, ctx.velocity_turbulance_graph->wide_graph.curve1, ctx.velocity_turbulance_graph->wide_graph.curve2, ctx.velocity_turbulance_graph->wide_graph.to) :
			tfx__wide_linear_sampler(ctx.velocity_turbulance_graph->wide_graph.from, ctx.velocity_turbulance_graph->wide_graph.to, velocity_turbulance_time);
		if (ctx.flags & tfx_ctx_policy_flag_velocity_turbulance_has_oscillator) {
			ctx.lookup_velocity_turbulance = tfxWideAdd(tfxWideMul(tfxOSCILLATOR_WIDE_SIN(velocity_turbulance_time, tfxWideAdd(ctx.velocity_turbulance_graph->wide_oscillator.offset_x, ctx.velocity_turbulance_graph->wide_oscillator.frequency), ctx.velocity_turbulance_graph->wide_oscillator.amplitude), ctx.lookup_velocity_turbulance), ctx.velocity_turbulance_graph->wide_oscillator.offset_y);
		}

		tfxWideFloat noise_resolution_time = ctx.noise_resolution_easing(ctx.life);
		ctx.lookup_noise_resolution = (ctx.flags & tfx_ctx_policy_flag_noise_resolution_is_bezier_graph) ?
			tfx__wide_bezier_sampler(noise_resolution_time, ctx.noise_resolution_graph->wide_graph.from, ctx.noise_resolution_graph->wide_graph.curve1, ctx.noise_resolution_graph->wide_graph.curve2, ctx.noise_resolution_graph->wide_graph.to) :
			tfx__wide_linear_sampler(ctx.noise_resolution_graph->wide_graph.from, ctx.noise_resolution_graph->wide_graph.to, noise_resolution_time);
		if (ctx.flags & tfx_ctx_policy_flag_noise_resolution_has_oscillator) {
			ctx.lookup_noise_resolution = tfxWideAdd(tfxWideMul(tfxOSCILLATOR_WIDE_SIN(noise_resolution_time, tfxWideAdd(ctx.noise_resolution_graph->wide_oscillator.offset_x, ctx.noise_resolution_graph->wide_oscillator.frequency), ctx.noise_resolution_graph->wide_oscillator.amplitude), ctx.lookup_noise_resolution), ctx.noise_resolution_graph->wide_oscillator.offset_y);
		}

		const tfxWideFloat noise_resolution = tfxWideLoad(&bank.noise_resolution[index]);
		const tfxWideFloat base_noise_offset = tfxWideLoad(&bank.noise_offset[index]);

		tfxWideFloat noise_offset = tfxWideMul(base_noise_offset, ctx.overal_scale_wide);

		tfx__readbarrier;

		ctx.lookup_noise_resolution = tfxWideMul(ctx.lookup_noise_resolution, noise_resolution);
		tfxWideFloat x = tfxWideAdd(tfxWideDiv(ctx.position_x.m, ctx.lookup_noise_resolution), noise_offset);
		tfxWideFloat y = tfxWideAdd(tfxWideDiv(ctx.position_y.m, ctx.lookup_noise_resolution), noise_offset);
		tfxWideFloat z = tfxWideAdd(tfxWideDiv(ctx.position_z.m, ctx.lookup_noise_resolution), noise_offset);

		tfxWideFloat y_offset = tfxWideAdd(y, tfxWideAdd(tfxWIDENOISEOFFSET.m, ctx.life));
		tfxWideFloat z_offset = tfxWideAdd(z, tfxWideAdd(tfxWIDENOISEOFFSET.m, ctx.life));

		tfxWideFloat noise_x = tfx__simd_noise_3d(x, y, z);
		tfxWideFloat noise_y = tfx__simd_noise_3d(x, y_offset, z);
		tfxWideFloat noise_z = tfx__simd_noise_3d(x, y, z_offset);

		tfxWideFloat l = tfxWideMul(noise_x, noise_x);
		l = tfxWideAdd(l, tfxWideMul(noise_y, noise_y));
		l = tfxWideAdd(l, tfxWideMul(noise_z, noise_z));
		l = tfxWideRSqrt(l);
		noise_x = tfxWideMul(noise_x, l);
		noise_y = tfxWideMul(noise_y, l);
		noise_z = tfxWideMul(noise_z, l);

		noise_x = tfxWideMul(ctx.global_noise, tfxWideMul(ctx.lookup_velocity_turbulance, noise_x));
		noise_y = tfxWideMul(ctx.global_noise, tfxWideMul(ctx.lookup_velocity_turbulance, noise_y));
		noise_z = tfxWideMul(ctx.global_noise, tfxWideMul(ctx.lookup_velocity_turbulance, noise_z));

		ctx.velocity_x = tfxWideAdd(ctx.velocity_x, noise_x);
		ctx.velocity_y = tfxWideAdd(ctx.velocity_y, noise_y);
		ctx.velocity_z = tfxWideAdd(ctx.velocity_z, noise_z);
	}
};

struct tfx_apply_curl_noise {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		tfxWideFloat velocity_turbulance_time = ctx.velocity_turbulance_easing(ctx.life);
		ctx.lookup_velocity_turbulance = (ctx.flags & tfx_ctx_policy_flag_velocity_turbulance_is_bezier_graph) ?
			tfx__wide_bezier_sampler(velocity_turbulance_time, ctx.velocity_turbulance_graph->wide_graph.from, ctx.velocity_turbulance_graph->wide_graph.curve1, ctx.velocity_turbulance_graph->wide_graph.curve2, ctx.velocity_turbulance_graph->wide_graph.to) :
			tfx__wide_linear_sampler(ctx.velocity_turbulance_graph->wide_graph.from, ctx.velocity_turbulance_graph->wide_graph.to, velocity_turbulance_time);
		if (ctx.flags & tfx_ctx_policy_flag_velocity_turbulance_has_oscillator) {
			ctx.lookup_velocity_turbulance = tfxWideAdd(tfxWideMul(tfxOSCILLATOR_WIDE_SIN(velocity_turbulance_time, tfxWideAdd(ctx.velocity_turbulance_graph->wide_oscillator.offset_x, ctx.velocity_turbulance_graph->wide_oscillator.frequency), ctx.velocity_turbulance_graph->wide_oscillator.amplitude), ctx.lookup_velocity_turbulance), ctx.velocity_turbulance_graph->wide_oscillator.offset_y);
		}

		tfxWideFloat noise_resolution_time = ctx.noise_resolution_easing(ctx.life);
		ctx.lookup_noise_resolution = (ctx.flags & tfx_ctx_policy_flag_noise_resolution_is_bezier_graph) ?
			tfx__wide_bezier_sampler(noise_resolution_time, ctx.noise_resolution_graph->wide_graph.from, ctx.noise_resolution_graph->wide_graph.curve1, ctx.noise_resolution_graph->wide_graph.curve2, ctx.noise_resolution_graph->wide_graph.to) :
			tfx__wide_linear_sampler(ctx.noise_resolution_graph->wide_graph.from, ctx.noise_resolution_graph->wide_graph.to, noise_resolution_time);
		if (ctx.flags & tfx_ctx_policy_flag_noise_resolution_has_oscillator) {
			ctx.lookup_noise_resolution = tfxWideAdd(tfxWideMul(tfxOSCILLATOR_WIDE_SIN(noise_resolution_time, tfxWideAdd(ctx.noise_resolution_graph->wide_oscillator.offset_x, ctx.noise_resolution_graph->wide_oscillator.frequency), ctx.noise_resolution_graph->wide_oscillator.amplitude), ctx.lookup_noise_resolution), ctx.noise_resolution_graph->wide_oscillator.offset_y);
		}

		const tfxWideFloat noise_resolution = tfxWideLoad(&bank.noise_resolution[index]);
		const tfxWideFloat base_noise_offset = tfxWideLoad(&bank.noise_offset[index]);

		tfxWideFloat noise_offset = tfxWideMul(base_noise_offset, ctx.overal_scale_wide);

		tfx__readbarrier;

		ctx.lookup_noise_resolution = tfxWideMul(ctx.lookup_noise_resolution, noise_resolution);
		tfxWideFloat x = tfxWideAdd(tfxWideDiv(ctx.position_x.m, ctx.lookup_noise_resolution), noise_offset);
		tfxWideFloat y = tfxWideAdd(tfxWideDiv(ctx.position_y.m, ctx.lookup_noise_resolution), noise_offset);
		tfxWideFloat z = tfxWideAdd(tfxWideDiv(ctx.position_z.m, ctx.lookup_noise_resolution), noise_offset);

		/*
		// 2. Calculate Offset Y coordinate for the 4 particles
		tfxWideFloat y_offset = tfxWideAdd(y, tfxWIDENOISEOFFSET.m);

		// 3. Perform Noise Samples (10 calls, each processing either 4 or 8 particles for SSE or AVX)
		// --- Samples needed for Curl X component
		tfxWideFloat y_minus_eps = tfx__simd_noise_3d(x, tfxWideSub(y, tfxWIDEEPS.m), z);
		tfxWideFloat y_plus_eps = tfx__simd_noise_3d(x, tfxWideAdd(y, tfxWIDEEPS.m), z);
		tfxWideFloat z_minus_eps = tfx__simd_noise_3d(x, y, tfxWideSub(z, tfxWIDEEPS.m));
		tfxWideFloat z_plus_eps = tfx__simd_noise_3d(x, y, tfxWideAdd(z, tfxWIDEEPS.m));

		// --- Samples needed for Curl Y component
		tfxWideFloat off_z_minus_eps = tfx__simd_noise_3d(x, y_offset, tfxWideSub(z, tfxWIDEEPS.m));
		tfxWideFloat off_z_plus_eps = tfx__simd_noise_3d(x, y_offset, tfxWideAdd(z, tfxWIDEEPS.m));
		tfxWideFloat off_x_minus_eps = tfx__simd_noise_3d(tfxWideSub(x, tfxWIDEEPS.m), y_offset, z);
		tfxWideFloat off_x_plus_eps = tfx__simd_noise_3d(tfxWideAdd(x, tfxWIDEEPS.m), y_offset, z);

		// --- Samples needed for Curl Z component
		// Re-use off_x_minus_eps, off_x_plus_eps
		tfxWideFloat off_y_minus_eps = tfx__simd_noise_3d(x, tfxWideSub(y_offset, tfxWIDEEPS.m), z);
		tfxWideFloat off_y_plus_eps = tfx__simd_noise_3d(x, tfxWideAdd(y_offset, tfxWIDEEPS.m), z);

		// 4. Calculate Partial Derivatives and Curl Components using SSE math
		// Derivatives for Curl X
		tfxWideFloat dy = tfxWideDiv(tfxWideSub(y_plus_eps, y_minus_eps), tfxWIDEEPS2.m);
		tfxWideFloat dz = tfxWideDiv(tfxWideSub(z_plus_eps, z_minus_eps), tfxWIDEEPS2.m);
		// Derivatives for Curl Y (offset)
		tfxWideFloat dz_off = tfxWideDiv(tfxWideSub(off_z_plus_eps, off_z_minus_eps), tfxWIDEEPS2.m);
		tfxWideFloat dx_off = tfxWideDiv(tfxWideSub(off_x_plus_eps, off_x_minus_eps), tfxWIDEEPS2.m);
		// Derivatives for Curl Z (offset)
		// Re-use dx_off
		tfxWideFloat dy_off = tfxWideDiv(tfxWideSub(off_y_plus_eps, off_y_minus_eps), tfxWIDEEPS2.m);

		// Curl Components for 4 particles
		tfxWideFloat noise_x = tfxWideSub(dy, dz);
		tfxWideFloat noise_y = tfxWideSub(dz_off, dx_off);
		tfxWideFloat noise_z = tfxWideSub(dx_off, dy_off);

		tfxWideFloat l = tfxWideMul(noise_x, noise_x);
		l = tfxWideAdd(l, tfxWideMul(noise_y, noise_y));
		l = tfxWideAdd(l, tfxWideMul(noise_z, noise_z));
		l = tfxWideRSqrt(l);
		noise_x = tfxWideMul(noise_x, l);
		noise_y = tfxWideMul(noise_y, l);
		noise_z = tfxWideMul(noise_z, l);
		*/

		const tfxWideFloat dt = tfxWideSetSingle(1e-4f); // Epsilon for finite difference
		const tfxWideFloat tfxWIDEMINUSONE = tfxWideSetSingle(-1.0f); // Assuming you have this

		// 1. Sample noise at 4 points
		//    n0 = N(x, y, z)
		//    nx = N(x+dt, y, z)
		//    ny = N(x, y+dt, z)
		//    nz = N(x, y, z+dt)
		tfxWideFloat n0 = tfx__simd_noise_3d(x, y, z);
		tfxWideFloat nx = tfx__simd_noise_3d(tfxWideAdd(x, dt), y, z);
		tfxWideFloat ny = tfx__simd_noise_3d(x, tfxWideAdd(y, dt), z);
		tfxWideFloat nz = tfx__simd_noise_3d(x, y, tfxWideAdd(z, dt));

		// 2. Calculate partial derivatives using forward differences
		//    grad_x = (N(x+dt, y, z) - N(x, y, z)) / dt  ~= dN/dx
		//    grad_y = (N(x, y+dt, z) - N(x, y, z)) / dt  ~= dN/dy
		//    grad_z = (N(x, y, z+dt) - N(x, y, z)) / dt  ~= dN/dz
		tfxWideFloat grad_x = tfxWideDiv(tfxWideSub(nx, n0), dt);
		tfxWideFloat grad_y = tfxWideDiv(tfxWideSub(ny, n0), dt);
		tfxWideFloat grad_z = tfxWideDiv(tfxWideSub(nz, n0), dt);

		// 3. Construct the 3D noise vector (inspired by common approximations/your previous code)
		//    noise_x = dN/dy - dN/dz
		//    noise_y = dN/dz - dN/dx
		//    noise_z = dN/dx - dN/dy
		//    This specific combination is one way to generate a divergence-free field from gradients.
		tfxWideFloat noise_x = tfxWideSub(grad_y, grad_z);
		tfxWideFloat noise_y = tfxWideSub(grad_z, grad_x);
		tfxWideFloat noise_z = tfxWideSub(grad_x, grad_y);

		// 4. Normalize (Optional but recommended for consistent speed)
		tfxWideFloat l_sq = tfxWideMul(noise_x, noise_x);
		l_sq = tfxWideAdd(l_sq, tfxWideMul(noise_y, noise_y));
		l_sq = tfxWideAdd(l_sq, tfxWideMul(noise_z, noise_z));
		// Add a small epsilon before rsqrt/sqrt to avoid division by zero
		tfxWideFloat epsilon_sq = tfxWideSetSingle(1e-12f); // Example small value
		l_sq = tfxWideAdd(l_sq, epsilon_sq);
		tfxWideFloat inv_l = tfxWideRSqrt(l_sq); // Reciprocal square root if available and precise enough
		// Or: inv_l = tfxWideDiv(tfxWideSetSingle(1.0f), tfxWideSqrt(l_sq));

		noise_x = tfxWideMul(noise_x, inv_l);
		noise_y = tfxWideMul(noise_y, inv_l);
		noise_z = tfxWideMul(noise_z, inv_l);

		noise_x = tfxWideMul(ctx.global_noise, tfxWideMul(ctx.lookup_velocity_turbulance, noise_x));
		noise_y = tfxWideMul(ctx.global_noise, tfxWideMul(ctx.lookup_velocity_turbulance, noise_y));
		noise_z = tfxWideMul(ctx.global_noise, tfxWideMul(ctx.lookup_velocity_turbulance, noise_z));

		ctx.velocity_x = tfxWideAdd(ctx.velocity_x, noise_x);
		ctx.velocity_y = tfxWideAdd(ctx.velocity_y, noise_y);
		ctx.velocity_z = tfxWideAdd(ctx.velocity_z, noise_z);
	}
};

struct tfx_apply_motion_randomness {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		//----Do the random calculation

		tfxWideInt velocity_normal = tfxWideSetZeroi;
		if (!(ctx.emitter->control_profile & tfxEmitterControlProfile_orbital)) {
			velocity_normal = tfxWideLoadi((tfxWideIntLoader *)&bank.velocity_normal[index]);
			tfx__wide_unpack10bit(velocity_normal, ctx.velocity_x, ctx.velocity_y, ctx.velocity_z);
		}

		tfxWideInt uid = tfxWideLoadi((tfxWideIntLoader *)&bank.uid[index]);
		tfxWideInt seed = tfx__wide_seedgen_base(ctx.time_step, uid);
		tfxWideFloat speed = tfxWideLoad(&bank.noise_offset[index]);

		tfxWideFloat motion_randomness_time = ctx.motion_randomness_easing(ctx.life);
		tfxWideFloat lookup_motion_randomness = (ctx.flags & tfx_ctx_policy_flag_motion_randomness_is_bezier_graph) ?
			tfx__wide_bezier_sampler(motion_randomness_time, ctx.motion_randomness_graph->wide_graph.from, ctx.motion_randomness_graph->wide_graph.curve1, ctx.motion_randomness_graph->wide_graph.curve2, ctx.motion_randomness_graph->wide_graph.to) :
			tfx__wide_linear_sampler(ctx.motion_randomness_graph->wide_graph.from, ctx.motion_randomness_graph->wide_graph.to, motion_randomness_time);
		if (ctx.flags & tfx_ctx_policy_flag_motion_randomness_has_oscillator) {
			lookup_motion_randomness = tfxWideAdd(tfxWideMul(tfxOSCILLATOR_WIDE_SIN(motion_randomness_time, tfxWideAdd(ctx.motion_randomness_graph->wide_oscillator.offset_x, ctx.motion_randomness_graph->wide_oscillator.frequency), ctx.motion_randomness_graph->wide_oscillator.amplitude), lookup_motion_randomness), ctx.motion_randomness_graph->wide_oscillator.offset_y);
		}
		const tfxWideFloat influence = tfxWideMul(tfxWideMul(ctx.motion_randomness_base, ctx.global_noise), lookup_motion_randomness);
		tfxWideFloat point_one_influence = tfxWideMul(tfxWideSetSingle(0.1f), influence);
		tfxWideFloat random_speed = tfxWideMul(tfxWideDiv(tfx__wide_seedgen(seed), tfxMAXUINTf.m), tfxWideMul(tfxWideSetSingle(0.01f), influence));
		tfxWideFloat random_x, random_y, random_z;
		tfx__wide_random_vector_in_cone(seed, ctx.velocity_x, ctx.velocity_y, ctx.velocity_z, tfxWideMul(tfxDEGREERANGEMR.m, influence), &random_x, &random_y, &random_z);
		speed = tfxWideAdd(speed, random_speed);
		tfxWideFloat length = tfxWideMul(random_x, random_x);
		length = tfxWideAdd(length, tfxWideMul(random_y, random_y));
		length = tfxWideAdd(length, tfxWideMul(random_z, random_z));
		length = tfxWideMul(tfxWideRSqrt(length), length);
		tfxWideFloat length_one = tfxWideDiv(tfxWIDEONE.m, length);
		random_x = tfxWideMul(random_x, length_one);
		random_y = tfxWideMul(random_y, length_one);
		random_z = tfxWideMul(random_z, length_one);
		ctx.velocity = tfxWideAdd(ctx.velocity, tfxWideMul(speed, ctx.global_noise));
		//----

		tfxWideInt packed_normal = tfxWideSetZeroi;
		if (ctx.emitter->control_profile & tfxEmitterControlProfile_orbital) {
			tfxWideFloat vx, vy, vz;
			tfxWideInt velocity_normal = tfxWideLoadi((tfxWideIntLoader *)&bank.velocity_normal[index]);
			tfx__wide_unpack10bit(velocity_normal, vx, vy, vz);
			vx = tfxWideAdd(tfxWideMul(random_x, ctx.time_step_fraction), tfxWideMul(vx, tfxWideSub(tfxWIDEONE.m, ctx.time_step_fraction)));
			vy = tfxWideAdd(tfxWideMul(random_y, ctx.time_step_fraction), tfxWideMul(vy, tfxWideSub(tfxWIDEONE.m, ctx.time_step_fraction)));
			vz = tfxWideAdd(tfxWideMul(random_z, ctx.time_step_fraction), tfxWideMul(vz, tfxWideSub(tfxWIDEONE.m, ctx.time_step_fraction)));
			ctx.velocity_x = tfxWideAdd(ctx.velocity_x, vx);
			ctx.velocity_y = tfxWideAdd(ctx.velocity_y, vy);
			ctx.velocity_z = tfxWideAdd(ctx.velocity_z, vz);
			length = tfxWideMul(ctx.velocity_x, ctx.velocity_x);
			length = tfxWideAdd(length, tfxWideMul(ctx.velocity_y, ctx.velocity_y));
			length = tfxWideAdd(length, tfxWideMul(ctx.velocity_z, ctx.velocity_z));
			length = tfxWideMul(tfxWideRSqrt(length), length);
			ctx.velocity_x = tfxWideDiv(ctx.velocity_x, length);
			ctx.velocity_y = tfxWideDiv(ctx.velocity_y, length);
			ctx.velocity_z = tfxWideDiv(ctx.velocity_z, length);
			length = tfxWideMul(vx, vx);
			length = tfxWideAdd(length, tfxWideMul(vy, vy));
			length = tfxWideAdd(length, tfxWideMul(vz, vz));
			length = tfxWideMul(tfxWideRSqrt(length), length);
			vx = tfxWideDiv(vx, length);
			vy = tfxWideDiv(vy, length);
			vz = tfxWideDiv(vz, length);
			packed_normal = tfx__wide_pack10bit_unsigned(vx, vy, vz);
			tfxWideInt normal_to_store = tfxWideOri(tfxWideAndi(packed_normal, ctx.time_changed_mask), tfxWideAndi(velocity_normal, tfxWideXOri(ctx.time_changed_mask, tfxWIDEMINUSONEi.m)));
			tfxWideStorei((tfxWideIntLoader *)&bank.velocity_normal[index], normal_to_store);
		} else {
			//Non Orbit emission direction
			//Add the random direction to the current velocity
			ctx.velocity_x = tfxWideAdd(tfxWideMul(random_x, ctx.time_step_fraction), tfxWideMul(ctx.velocity_x, tfxWideSub(tfxWIDEONE.m, ctx.time_step_fraction)));
			ctx.velocity_y = tfxWideAdd(tfxWideMul(random_y, ctx.time_step_fraction), tfxWideMul(ctx.velocity_y, tfxWideSub(tfxWIDEONE.m, ctx.time_step_fraction)));
			ctx.velocity_z = tfxWideAdd(tfxWideMul(random_z, ctx.time_step_fraction), tfxWideMul(ctx.velocity_z, tfxWideSub(tfxWIDEONE.m, ctx.time_step_fraction)));
			length = tfxWideMul(ctx.velocity_x, ctx.velocity_x);
			length = tfxWideAdd(length, tfxWideMul(ctx.velocity_y, ctx.velocity_y));
			length = tfxWideAdd(length, tfxWideMul(ctx.velocity_z, ctx.velocity_z));
			length = tfxWideMul(tfxWideRSqrt(length), length);
			ctx.velocity_x = tfxWideDiv(ctx.velocity_x, length);
			ctx.velocity_y = tfxWideDiv(ctx.velocity_y, length);
			ctx.velocity_z = tfxWideDiv(ctx.velocity_z, length);
			packed_normal = tfx__wide_pack10bit_unsigned(ctx.velocity_x, ctx.velocity_y, ctx.velocity_z);
			tfxWideInt normal_to_store = tfxWideOri(tfxWideAndi(packed_normal, ctx.time_changed_mask), tfxWideAndi(velocity_normal, tfxWideXOri(ctx.time_changed_mask, tfxWIDEMINUSONEi.m)));
			tfxWideStorei((tfxWideIntLoader *)&bank.velocity_normal[index], normal_to_store);
			//--
		}
		tfxWideStore(&bank.noise_offset[index], speed);

		ctx.velocity_x = tfxWideMul(ctx.velocity_x, ctx.velocity);
		ctx.velocity_y = tfxWideMul(ctx.velocity_y, ctx.velocity);
		ctx.velocity_z = tfxWideMul(ctx.velocity_z, ctx.velocity);
	}
};

struct tfx_apply_velocity {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		tfxWideInt velocity_normal = tfxWideLoadi((tfxWideIntLoader *)&bank.velocity_normal[index]);
		tfx__wide_unpack10bit(velocity_normal, ctx.velocity_x, ctx.velocity_y, ctx.velocity_z);
		ctx.velocity_x = tfxWideMul(ctx.velocity_x, ctx.velocity);
		ctx.velocity_y = tfxWideMul(ctx.velocity_y, ctx.velocity);
		ctx.velocity_z = tfxWideMul(ctx.velocity_z, ctx.velocity);
	}
};

struct tfx_apply_load_position {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		ctx.position_x.m = tfxWideLoad(&bank.position_x[index]);
		ctx.position_y.m = tfxWideLoad(&bank.position_y[index]);
		ctx.position_z.m = tfxWideLoad(&bank.position_z[index]);
	}
};

struct tfx_apply_load_path {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		ctx.path_position = tfxWideLoad(&bank.path_position[index]);
		ctx.path_offset = tfxWideLoad(&bank.path_offset[index]);
		ctx.path_scale_variation = tfxWideLoad(&bank.path_scale_variation[index]);
	}
};

struct tfx_apply_update_path_position {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		ctx.path_position = tfxWideAdd(ctx.path_position, tfxWideMul(ctx.velocity, pm->update_time_wide));
	}
};

struct tfx_apply_update_path_position_trajectory {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		ctx.path_position = tfxWideMul(ctx.lookup_velocity, ctx.node_count);
	}
};

struct tfx_apply_path_end_kill {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		//Kill if the particle has reached the end of the path
		tfxWideInt remove_flag = tfxWideSetSinglei(tfxParticleFlags_remove);
		tfxWideInt remove_flags = tfxWideAndi(remove_flag, tfxWideOri(tfxWideCasti(tfxWideLess(ctx.path_position, tfxWideSetZero)), tfxWideCasti(tfxWideGreaterEqual(ctx.path_position, ctx.node_count))));
		ctx.path_position = tfxWideMax(ctx.path_position, tfxWideSetZero);
		tfxWideInt flags = tfxWideLoadi((tfxWideIntLoader *)&bank.flags_single_loop_count[index]);
		flags = tfxWideOri(flags, remove_flags);
		tfxWideStorei((tfxWideIntLoader *)&bank.flags_single_loop_count[index], flags);
	}
};

struct tfx_apply_path_end_loop {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		//Reposition if the particle is travelling along the path
		tfxWideFloat at_end = tfxWideGreaterEqual(ctx.path_position, ctx.node_count);
		ctx.path_position = tfxWideSub(ctx.path_position, tfxWideAnd(at_end, ctx.node_count));
		tfxWideInt flags = tfxWideLoadi((tfxWideIntLoader *)&bank.flags_single_loop_count[index]);
		flags = tfxWideOri(flags, tfxWideAndi(ctx.capture_after_transform_flag, tfxWideCasti(at_end)));
		at_end = tfxWideLess(ctx.path_position, tfxWideSetZero);
		ctx.path_position = tfxWideAdd(ctx.path_position, tfxWideAnd(at_end, ctx.node_count));
		flags = tfxWideOri(flags, tfxWideAndi(ctx.capture_after_transform_flag, tfxWideCasti(at_end)));
		tfxWideStorei((tfxWideIntLoader *)&bank.flags_single_loop_count[index], flags);
	}
};

struct tfx_apply_path_position {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		tfxWideArrayi node_index;
		node_index.m = tfxWideConverti(ctx.path_position);
		tfxWideArray t;
		t.m = tfxWideSub(ctx.path_position, tfxWideConvert(node_index.m));
		tfxWideArray point_x;
		tfxWideArray point_z;
		tfx__wide_catmull_rom_spline_3d(&node_index, t.m, ctx.path->buffers.node_soa.x, ctx.path->buffers.node_soa.y, ctx.path->buffers.node_soa.z, &point_x.m, &ctx.position_y.m, &point_z.m);
		if (ctx.path->settings.extrusion_type == tfxExtrusionArc) {
			tfxWideFloat radius = tfxWideAdd(tfxWideMul(point_x.m, point_x.m), tfxWideMul(point_z.m, point_z.m));
			tfxWideFloat length_mask = tfxWideGreater(radius, tfxWideSetZero);
			radius = tfxWideMul(tfxWideRSqrt(radius), radius);
			tfxWideArray angle;
			tfxWideArray rx;
			tfxWideArray rz;
			angle.m = tfxWideAtan2(point_z.m, point_x.m);
			angle.m = tfxWideAnd(length_mask, angle.m);
			angle.m = tfxWideAdd(angle.m, ctx.path_offset);
			tfxWideSinCos(angle.m, &rz.m, &rx.m);
			ctx.position_x.m = tfxWideMul(rx.m, radius);
			ctx.position_z.m = tfxWideMul(rz.m, radius);
		} else {
			ctx.position_x.m = tfxWideAdd(point_x.m, ctx.path_offset);
			ctx.position_z.m = point_z.m;
		}
		ctx.position_x.m = tfxWideAdd(ctx.position_x.m, ctx.emitter_offset_x);
		ctx.position_y.m = tfxWideAdd(ctx.position_y.m, ctx.emitter_offset_y);
		ctx.position_z.m = tfxWideAdd(ctx.position_z.m, ctx.emitter_offset_z);
		if (ctx.emitter->state_flags & tfxEmitterStateFlags_has_rotated_path) {
			tfxWideArrayi quat_xy, quat_zw;
			for (int i = 0; i < tfxDataWidth; i++) {
				tfxU64 q = bank.quaternion[index + i];
				quat_xy.a[i] = (int)(uint32_t)(q & 0xFFFFFFFF);
				quat_zw.a[i] = (int)(uint32_t)(q >> 32);
			}
			tfx__wide_transform_packed_quaternion_vec3(&quat_xy.m, &quat_zw.m, &ctx.position_x.m, &ctx.position_y.m, &ctx.position_z.m);
		}
	}
};

struct tfx_apply_path_scale_variation {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		ctx.position_x.m = tfxWideMul(ctx.position_x.m, ctx.path_scale_variation);
		ctx.position_y.m = tfxWideMul(ctx.position_y.m, ctx.path_scale_variation);
		ctx.position_z.m = tfxWideMul(ctx.position_z.m, ctx.path_scale_variation);
	}
};

struct tfx_apply_store_path_position {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		tfxWideStore(&bank.position_x[index], tfxWideMul(ctx.position_x.m, tfxWideMul(ctx.emitter_width, ctx.path_scale_variation)));
		tfxWideStore(&bank.position_y[index], tfxWideMul(ctx.position_y.m, tfxWideMul(ctx.emitter_height, ctx.path_scale_variation)));
		tfxWideStore(&bank.position_z[index], tfxWideMul(ctx.position_z.m, tfxWideMul(ctx.emitter_depth, ctx.path_scale_variation)));
		tfxWideStore(&bank.path_position[index], ctx.path_position);
	}
};

struct tfx_apply_store_path_position_edge_traversal {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		tfxWideStore(&bank.position_x[index], tfxWideMul(ctx.position_x.m, ctx.emitter_width));
		tfxWideStore(&bank.position_y[index], tfxWideMul(ctx.position_y.m, ctx.emitter_height));
		tfxWideStore(&bank.position_z[index], tfxWideMul(ctx.position_z.m, ctx.emitter_depth));
		tfxWideStore(&bank.path_position[index], ctx.path_position);
	}
};

struct tfx_apply_orbital_velocity_normal {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		ctx.velocity_x = tfxWideMul(tfxWideSub(ctx.position_z.m, ctx.emitter_offset_z), tfxWIDEMINUSONE.m);
		ctx.velocity_y = tfxWideSetZero;
		ctx.velocity_z = tfxWideSub(ctx.position_x.m, ctx.emitter_offset_x);
		tfxWideFloat l = tfxWideMul(ctx.velocity_x, ctx.velocity_x);
		l = tfxWideAdd(l, tfxWideMul(ctx.velocity_z, ctx.velocity_z));
		l = tfxWideMul(tfxWideRSqrt(l), l);
		ctx.velocity_x = tfxWideDiv(ctx.velocity_x, l);
		ctx.velocity_z = tfxWideDiv(ctx.velocity_z, l);
	}
};

struct tfx_apply_orbital_scale_velocity {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		ctx.velocity_x = tfxWideMul(ctx.velocity_x, ctx.velocity);
		ctx.velocity_z = tfxWideMul(ctx.velocity_z, ctx.velocity);
	}
};

struct tfx_apply_position {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		tfxWideFloat age = tfxWideLoad(&bank.age[index]);
		tfxWideFloat age_fraction = tfxWideMin(tfxWideDiv(age, pm->frame_length_wide), tfxWIDEONE.m);
		//tfxWideFloat age_fraction = tfxWIDEONE.m;
		ctx.velocity_y = tfxWideSub(ctx.velocity_y, ctx.weight);
		ctx.velocity_x = tfxWideMul(tfxWideMul(tfxWideMul(ctx.velocity_x, pm->update_time_wide), ctx.velocity_adjuster), age_fraction);
		ctx.velocity_y = tfxWideMul(tfxWideMul(tfxWideMul(ctx.velocity_y, pm->update_time_wide), ctx.velocity_adjuster), age_fraction);
		ctx.velocity_z = tfxWideMul(tfxWideMul(tfxWideMul(ctx.velocity_z, pm->update_time_wide), ctx.velocity_adjuster), age_fraction);
		ctx.position_x.m = tfxWideAdd(ctx.position_x.m, tfxWideMul(ctx.velocity_x, ctx.overal_scale_wide));
		ctx.position_y.m = tfxWideAdd(ctx.position_y.m, tfxWideMul(ctx.velocity_y, ctx.overal_scale_wide));
		ctx.position_z.m = tfxWideAdd(ctx.position_z.m, tfxWideMul(ctx.velocity_z, ctx.overal_scale_wide));
		tfxWideStore(&bank.position_x[index], ctx.position_x.m);
		tfxWideStore(&bank.position_y[index], ctx.position_y.m);
		tfxWideStore(&bank.position_z[index], ctx.position_z.m);
	}
};

struct tfx_apply_pack_velocity {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		tfxWideFloat length = tfxWideMul(ctx.velocity_x, ctx.velocity_x);
		length = tfxWideAdd(length, tfxWideMul(ctx.velocity_y, ctx.velocity_y));
		length = tfxWideAdd(length, tfxWideMul(ctx.velocity_z, ctx.velocity_z));
		length = tfxWideMul(tfxWideRSqrt(length), length);
		ctx.velocity_x = tfxWideDiv(ctx.velocity_x, length);
		ctx.velocity_y = tfxWideDiv(ctx.velocity_y, length);
		ctx.velocity_z = tfxWideDiv(ctx.velocity_z, length);
		tfxWideInt packed_normal = tfx__wide_pack10bit_unsigned(ctx.velocity_x, ctx.velocity_y, ctx.velocity_z);
		tfxWideStorei((tfxWideIntLoader *)&bank.velocity_normal[index], packed_normal);
	}
};

struct tfx_apply_position_line_trajectory {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		tfxWideFloat path_scale_variation = tfxWideLoad(&bank.path_scale_variation[index]);
		ctx.position_x.m = tfxWideSetZero;
		ctx.position_y.m = tfxWideMul(tfxWideMul(tfxWideMul(ctx.lookup_velocity, ctx.emitter_height), ctx.overal_scale_wide), path_scale_variation);
		ctx.position_z.m = tfxWideSetZero;
		if (ctx.emitter->control_profile & tfxEmitterControlProfile_rotated_line) {
			tfxWideArrayi quat_xy, quat_zw;
			for (int i = 0; i < tfxDataWidth; i++) {
				tfxU64 q = bank.quaternion[index + i];
				quat_xy.a[i] = (int)(uint32_t)(q & 0xFFFFFFFF);
				quat_zw.a[i] = (int)(uint32_t)(q >> 32);
			}
			tfx__wide_transform_packed_quaternion_vec3(&quat_xy.m, &quat_zw.m, &ctx.position_x.m, &ctx.position_y.m, &ctx.position_z.m);
		}
	}
};

struct tfx_apply_store_position {
	static inline void apply(tfxU32 index, tfx_effect_manager pm, tfx_particle_soa_t &bank, tfx_position_policy_context &ctx) {
		tfxWideStore(&bank.position_x[index], ctx.position_x.m);
		tfxWideStore(&bank.position_y[index], ctx.position_y.m);
		tfxWideStore(&bank.position_z[index], ctx.position_z.m);
	}
};

template <typename... Policies>
void tfx__setup_particles_position(tfx_control_work_entry_t *work_entry, tfx_position_policy_context &ctx);

template <typename... Policies>
void tfx__update_particles_position(tfx_control_work_entry_t *work_entry, tfx_position_policy_context &ctx);

//--------------------------------
//effect manager internal functions
//--------------------------------
template<typename T>
tfxINTERNAL inline void tfx__write_particle_color_sprite_data(T *sprites, tfxU32 start_diff, tfxU32 limit_index, const tfxU32 *depth_index, tfxU32 index, const tfxWideArrayi &packed_intensity_life, const tfxWideArrayi &curved_alpha_life, tfxU32 &running_sprite_index) {
	for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
		sprites[running_sprite_index].intensity_gradient_map.packed = packed_intensity_life.a[j];
		sprites[running_sprite_index].curved_alpha_life.packed = curved_alpha_life.a[j];
		running_sprite_index++;
	}
}

template<typename T>
tfxINTERNAL inline void tfx__write_particle_color_sprite_data_ordered(T *sprites, tfxU32 layer, tfxU32 start_diff, tfxU32 limit_index, const tfxU32 *depth_index, tfxU32 index, const tfxWideArrayi &packed_intensity_life, const tfxWideArrayi &curved_alpha_life, tfxU32 &running_sprite_index, tfxU32 instance_offset) {
	for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
		tfxU32 sprite_depth_index = depth_index[index + j] + instance_offset;
		sprites[sprite_depth_index].intensity_gradient_map.packed = packed_intensity_life.a[j];
		sprites[sprite_depth_index].curved_alpha_life.packed = curved_alpha_life.a[j];
		running_sprite_index++;
	}
}

template<typename T>
tfxINTERNAL inline void tfx__write_particle_image_sprite_data(T *sprites, tfx_effect_manager pm, tfxU32 layer, tfxU32 start_diff, tfxU32 limit_index, tfx_particle_soa_t &bank, tfxWideArrayi &flags, tfxWideArrayi &image_indexes, const tfxEmitterStateFlags emitter_flags, const tfx_billboarding_option billboard_option, tfxU32 index, tfxU32 &running_sprite_index) {
	for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
		int index_j = index + j;
		tfxU32 &sprites_index = bank.sprite_index[index_j];
		tfxU32 capture = flags.a[j] << 7;
		sprites[running_sprite_index].captured_index = capture == 0 ? (pm->current_sprite_buffer << 30) + running_sprite_index : (!pm->current_sprite_buffer << 30) + (sprites_index & 0x0FFFFFFF);
		sprites[running_sprite_index].captured_index |= emitter_flags & tfxEmitterStateFlags_wrap_single_sprite ? 0x80000000 : 0;
		sprites_index = layer + running_sprite_index;
		sprites[running_sprite_index].indexes = image_indexes.a[j];
		sprites[running_sprite_index].indexes |= (billboard_option << 13) | capture;
		bank.flags_single_loop_count[index_j] &= ~tfxParticleFlags_capture_after_transform;
		running_sprite_index++;
	}
}

template<typename T>
tfxINTERNAL inline void tfx__write_particle_image_sprite_data_ordered(T *sprites, tfx_effect_manager pm, tfxU32 layer, tfxU32 start_diff, tfxU32 limit_index, tfx_particle_soa_t &bank, tfxWideArrayi &flags, tfxWideArrayi &image_indexes, const tfxEmitterStateFlags emitter_flags, const tfx_billboarding_option billboard_option, tfxU32 index, tfxU32 &running_sprite_index, tfxU32 instance_offset) {
	for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
		int index_j = index + j;
		tfxU32 sprite_depth_index = bank.depth_index[index_j] + instance_offset;
		tfxU32 &sprites_index = bank.sprite_index[index_j];
		tfxU32 capture = flags.a[j] << 7;
		sprites[sprite_depth_index].captured_index = capture == 0 && (bank.flags_single_loop_count[index_j] & 0xFF) == 0 ? (pm->current_sprite_buffer << 30) + sprite_depth_index : (!pm->current_sprite_buffer << 30) + (sprites_index & 0x0FFFFFFF);
		sprites[sprite_depth_index].captured_index |= emitter_flags & tfxEmitterStateFlags_wrap_single_sprite ? 0x80000000 : 0;
		sprites_index = layer + sprite_depth_index;
		sprites[sprite_depth_index].indexes = image_indexes.a[j];
		sprites[sprite_depth_index].indexes |= (billboard_option << 13) | capture;
		bank.flags_single_loop_count[index_j] &= ~tfxParticleFlags_capture_after_transform;
		running_sprite_index++;
	}
}

template<typename T>
tfxINTERNAL inline void tfx__wrap_single_particle_instances(T *instance, tfx_sprite_data_t *sprite_data) {
	tfx_sprite_data_soa_t &sprites = sprite_data->real_time_sprites;
	for (tfxEachLayer) {
		for (int i = sprite_data->normal.frame_meta[0].index_offset[layer]; i != sprite_data->normal.frame_meta[0].index_offset[layer] + sprite_data->normal.frame_meta[0].sprite_count[layer]; ++i) {
			if (instance[i].captured_index != tfxINVALID && instance[i].captured_index & 0x80000000) {
				for (int j = sprite_data->normal.frame_meta[sprite_data->normal.frame_count - 1].index_offset[layer]; j != sprite_data->normal.frame_meta[sprite_data->normal.frame_count - 1].index_offset[layer] + sprite_data->normal.frame_meta[sprite_data->normal.frame_count - 1].sprite_count[layer]; ++j) {
					if (sprites.uid[j].uid == sprites.uid[i].uid) {
						instance[i].captured_index = j;
					}
				}
			}
		}
	}
}

template<typename T> tfxINTERNAL void tfx__clear_wrap_bit(T* instance, tfx_sprite_data_t *sprite_data);
template<typename T> tfxINTERNAL void tfx__sprite_data_offset_captured_indexes(T* instance, tfx_sprite_data_t *sprite_data, tfxU32 previous_frame, tfxU32 current_frame);
template<typename T> tfxINTERNAL void tfx__link_sprite_data_captured_indexes(T* instance, int entry_frame, tfx_sprite_data_t *sprite_data);
tfxINTERNAL void tfx__ribbon_data_set_captured_indexes(tfx_sprite_data_t *sprite_data);

template<typename T>
tfxINTERNAL inline void tfx__invalidate_new_captured_index(T* instance, tfx_vector_t<tfx_unique_sprite_id_t> &uids, tfx_effect_manager pm, tfxU32 layer) {
	for (tfxU32 i = 0; i != pm->instance_buffer_for_recording[pm->current_sprite_buffer][layer].current_size; ++i) {
		if ((uids[i].age == 0 && !(instance[i].captured_index & 0x80000000)) || ((instance[i].captured_index & 0xC0000000) >> 30 == pm->current_sprite_buffer && !(instance[i].captured_index & 0x80000000))) {
			instance[i].captured_index = tfxINVALID;
		}
	}
}

template<typename T>
tfxINTERNAL inline void tfx__invalidate_offsetted_sprite_captured_index(T* instance, tfx_vector_t<tfx_unique_sprite_id_t> &uids, tfx_effect_manager pm, tfxU32 layer) {
	for (tfxU32 i = 0; i != pm->instance_buffer_for_recording[pm->current_sprite_buffer][layer].current_size; ++i) {
		instance[i].captured_index = tfxINVALID;
	}
}

tfxINTERNAL tfx_vec3_t tfx__random_vector_in_cone(tfx_random_t *random, tfx_vec3_t cone_direction, float cone_angle);
tfxINTERNAL void tfx__transform_effect(tfx_vec3_t *world_rotations, tfx_vec3_t *local_rotations, tfx_vec3_t *world_position, tfx_vec3_t *local_position, tfx_quaternion_t *q, tfx_sprite_transform_t *parent, bool relative_position = true, bool relative_angle = false);
tfxINTERNAL void tfx__update_effect(tfx_effect_manager pm, tfxU32 index, tfxU32 parent_index = tfxINVALID);
tfxINTERNAL void tfx__update_emitter(tfx_work_queue_t *work_queue, void *data);
tfxINTERNAL void tfx__update_ribbon_bucket_emitters(tfx_work_queue_t *work_queue, void *data);
tfxINTERNAL void tfx__update_ribbon_emitter(tfxU32 ribbon_index, tfx_work_queue_t *work_queue, void *data);
tfxINTERNAL tfxU32 tfx__new_sprites_needed(tfx_effect_manager pm, tfx_random_t *random, tfxU32 index, tfx_effect_state_t *parent, tfx_shared_properties_t *shared_properties);
tfxINTERNAL tfxU32 tfx__new_ribbons_needed(tfx_effect_manager pm, tfx_random_t *random, tfxU32 index, tfx_effect_state_t *parent, tfx_shared_properties_t *shared_properties);
tfxINTERNAL void tfx__update_emitter_state(tfx_effect_manager pm, tfx_particle_emitter_state_t &emitter, tfxU32 parent_index, const tfx_parent_spawn_controls_t *parent_spawn_controls, tfx_spawn_work_entry_t *entry);
tfxINTERNAL void tfx__update_ribbon_emitter_state(tfx_effect_manager pm, tfx_ribbon_emitter_state_t &ribbon, tfxU32 parent_index, const tfx_parent_spawn_controls_t *parent_spawn_controls, tfx_ribbon_work_entry_t *entry);
tfxINTERNAL void tfx__update_effect_state(tfx_effect_manager pm, tfxU32 index);
tfxINTERNAL tfx_effect_manager tfx__next_global_effect_manager();
tfxINTERNAL tfx_library tfx__next_global_library();

tfxINTERNAL void tfx__spawn_particles(tfx_effect_manager pm, tfx_spawn_work_entry_t *work_entry);
tfxINTERNAL void tfx__spawn_particle_noise(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_motion_randomness(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_weight(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_velocity(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_image_frame(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_age(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void tfx__do_spawn_work(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void tfx__spawn_particle_point(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_other_emitter(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_other_ribbon_emitter(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_other_emitter_single(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_line(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_line_start(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_area(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_ellipsoid(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_cylinder(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_icosphere_random(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_icosphere(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_path(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_path_start(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_micro_update(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_spin(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_size(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void tfx__spawn_static_ribbons(tfxU32 ribbon_emitter_index, tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_ribbon_path(tfx_work_queue_t *queue, void *data);

typedef float(*tfx_easing_function)(float);
typedef float(*tfx_bezier_function)(float, float, float, float, float);

tfxINTERNAL tfxWideFloat tfx__wide_ease_constant(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_smoothstep(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_linear(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_in_quad(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_out_quad(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_in_out_quad(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_out_in(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_in_cubic(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_out_cubic(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_in_out_cubic(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_in_quart(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_out_quart(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_in_out_quart(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_in_quint(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_out_quint(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_in_out_quint(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_in_circular(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_out_circular(tfxWideFloat t);
tfxINTERNAL tfxWideFloat tfx__wide_ease_in_out_circular(tfxWideFloat t);

tfxINTERNAL float tfx__ease_constant(float t);
tfxINTERNAL float tfx__ease_smoothstep(float t);
tfxINTERNAL float tfx__ease_linear(float t);
tfxINTERNAL float tfx__ease_in_quad(float t);
tfxINTERNAL float tfx__ease_out_quad(float t);
tfxINTERNAL float tfx__ease_out_in(float t);
tfxINTERNAL float tfx__ease_in_out_quad(float t);
tfxINTERNAL float tfx__ease_in_cubic(float t);
tfxINTERNAL float tfx__ease_out_cubic(float t);
tfxINTERNAL float tfx__ease_in_out_cubic(float t);
tfxINTERNAL float tfx__ease_in_quart(float t);
tfxINTERNAL float tfx__ease_out_quart(float t);
tfxINTERNAL float tfx__ease_in_out_quart(float t);
tfxINTERNAL float tfx__ease_in_quint(float t);
tfxINTERNAL float tfx__ease_out_quint(float t);
tfxINTERNAL float tfx__ease_in_out_quint(float t);
tfxINTERNAL float tfx__ease_in_circular(float t);
tfxINTERNAL float tfx__ease_out_circular(float t);
tfxINTERNAL float tfx__ease_in_out_circular(float t);

tfxINTERNAL tfx_wide_easing_function tfx__get_wide_easing_function(tfx_graph_easing_type type);
tfxINTERNAL tfx_easing_function tfx__get_easing_function(tfx_graph_easing_type type);

tfxINTERNAL inline float tfx__bezier_sampler(float t, float node1, float curve1, float curve2, float node2) {
	float u = 1 - t;
	float w1 = u * u * u;
	float w2 = 3 * u * u * t;
	float w3 = 3 * u * t * t;
	float w4 = t * t * t;
	return w1 * node1 + w2 * curve1 + w3 * curve2 + w4 * node2;
}

tfxINTERNAL inline float tfx__linear_sampler(float from, float to, float t) {
	return from + (to - from) * t;
}

tfxINTERNAL inline bool tfx__graph_can_oscillate(tfx_graph_t *graph) {
	return graph->easing_type != tfxGraphEasingType_constant && graph->flags & tfxGraphFlags_enable_oscillator && graph->oscillator.amplitude != 0.f && graph->oscillator.frequency != 0.f;
}

tfxINTERNAL inline bool tfx__graph_has_bezier_curves(tfx_graph_t *graph) {
	return graph->flags & tfxGraphFlags_use_bezier_sampling && graph->easing_type != tfxGraphEasingType_constant;
}

tfxAPI_EDITOR float tfx__sample_graph(tfx_graph_t *graph, float t);

tfxINTERNAL void tfx__control_particles(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_age(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_image_frame(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_color(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_size(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_hide(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_spin_roll(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_spin_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_uid(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_capture_spawn_locations(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void tfx__control_particle_line_behaviour_kill(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_line_behaviour_loop(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void tfx__control_particle_transform(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_bounding_box(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void tfx__control_ribbons(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_ribbons_ages(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_ribbon_path_age(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_ribbon_attributes(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_ribbon_paths(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void tfx__update_ribbon_buffer_requirements(tfx_effect_manager pm);
tfxINTERNAL void tfx__reset_ribbon_buffer_requirements(tfx_effect_manager pm);
tfxAPI_EDITOR bool tfx__next_ribbon_bucket(tfx_effect_manager pm, tfx_ribbon_dispatch_t *ribbon_dispatch);

tfxINTERNAL bool tfx__control_profile_has_noise(tfxEmitterControlProfileFlags flags);
tfxINTERNAL void tfx__init_sprite_data_soa(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void tfx__init_particle_soa(tfx_soa_buffer_t *buffer, tfx_particle_soa_t *soa, tfxU32 reserve_amount, tfxEmitterControlProfileFlags control_profile);
tfxINTERNAL void tfx__init_particle_location_soa(tfx_soa_buffer_t *buffer, tfx_spawn_points_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void tfx__init_ribbons_soa(tfx_soa_buffer_t *buffer, tfx_ribbon_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void tfx__init_ribbon_data_soa(tfx_soa_buffer_t *buffer, tfx_ribbon_data_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void tfx__init_ribbon_segment_soa(tfx_soa_buffer_t *buffer, tfx_ribbon_segment_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void tfx__copy_emitter_properties(tfx_particle_emitter_properties_t *from_properties, tfx_particle_emitter_properties_t *to_properties);
tfxINTERNAL void tfx__copy_ribbon_properties(tfx_ribbon_emitter_properties_t *from_properties, tfx_ribbon_emitter_properties_t *to_properties);
tfxINTERNAL inline void tfx__free_sprite_data(tfx_sprite_data_t *sprite_data);
tfxINTERNAL inline bool tfx__is_graph_transform_rotation(tfx_graph_type type) {
	return type == tfxTransform_roll || type == tfxTransform_pitch || type == tfxTransform_yaw;
}

tfxINTERNAL inline bool tfx__is_graph_emitter_dimension(tfx_graph_type type) {
	return type == tfxProperty_emitter_width || type == tfxProperty_emitter_height || type == tfxProperty_emitter_depth;
}

tfxINTERNAL inline bool tfx__is_graph_translation(tfx_graph_type type) {
	return type == tfxTransform_translate_x || type == tfxTransform_translate_y || type == tfxTransform_translate_z;
}

tfxINTERNAL inline bool tfx__is_graph_base_spin(tfx_graph_type type) {
	return type == tfxBase_roll_spin || type == tfxBase_pitch_spin || type == tfxBase_yaw_spin;
}

tfxINTERNAL inline bool tfx__is_graph_variation_spin(tfx_graph_type type) {
	return type == tfxVariation_roll_spin || type == tfxVariation_pitch_spin || type == tfxVariation_yaw_spin;
}

tfxINTERNAL inline bool tfx__is_graph_overtime_spin(tfx_graph_type type) {
	return type == tfxOvertime_roll_spin || type == tfxOvertime_pitch_spin || type == tfxOvertime_yaw_spin;
}

tfxINTERNAL inline bool tfx__is_graph_emission(tfx_graph_type type) {
	return type == tfxProperty_emission_pitch || type == tfxProperty_emission_yaw;
}

tfxINTERNAL inline bool tfx__is_graph_particle_size(tfx_graph_type type) {
	return    type == tfxBase_width || type == tfxBase_height ||
		type == tfxVariation_width || type == tfxVariation_height ||
		type == tfxOvertime_width || type == tfxOvertime_height;
}

tfxAPI_EDITOR tfxU32 tfx__create_emitter_path_attributes(tfx_effect_descriptor emitter, bool add_node);
tfxINTERNAL void tfx__initialise_effect_graphs(tfx_graph_list_t *graph_list, tfxU32 bucket_size = 8);
tfxINTERNAL void tfx__initialise_emitter_graphs(tfx_graph_list_t *graph_list, tfxU32 bucket_size = 8);
tfxINTERNAL void tfx__initialise_ribbon_graphs(tfx_graph_list_t *graph_list, tfxU32 bucket_size = 8);
tfxINTERNAL tfxErrorFlags tfx__load_effect_library_package(tfx_package package, tfx_library lib, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset), void *user_data = nullptr);
tfxINTERNAL void tfx__build_gpu_shape_data(tfx_vector_t<tfx_image_data_t> *particle_shapes, tfx_gpu_shapes shape_data, void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset));

//--------------------------------
//Animation manager internal functions - animation manager is used to playback pre-recorded effects
//--------------------------------
tfxINTERNAL tfxAnimationID tfx__allocate_animation_instance(tfx_animation_manager animation_manager);
tfxINTERNAL void tfx__free_animation_instance(tfx_animation_manager animation_manager, tfxU32 index);
tfxINTERNAL void tfx__add_effect_emitter_properties(tfx_animation_manager animation_manager, tfx_effect_descriptor effect, bool *has_animated_shape);
tfxINTERNAL bool tfx__free_pm_effect_capacity(tfx_effect_manager pm);
tfxINTERNAL tfx_animation_manager tfx__create_animation_manager(tfxU32 max_instances);

//--------------------------------
//effect manager internal functions
//--------------------------------
tfxINTERNAL tfx_effect_index_t tfx__get_effect_slot(tfx_effect_manager pm);
tfxINTERNAL tfxU32 tfx__get_emitter_slot(tfx_effect_manager pm);
tfxINTERNAL tfxU32 tfx__get_ribbon_slot(tfx_effect_manager pm);
tfxINTERNAL tfxU32 tfx__get_particle_index_slot(tfx_effect_manager pm, tfxParticleID particle_id);
tfxINTERNAL tfxU32 tfx__allocate_path_quaternion(tfx_effect_manager pm, tfxU32 amount);
tfxINTERNAL void tfx__free_path_quaternion(tfx_effect_manager pm, tfxU32 index);
tfxINTERNAL void tfx__free_particle_index(tfx_effect_manager pm, tfxU32 *index);
tfxINTERNAL tfxU32 tfx__push_depth_index(tfx_vector_t<tfx_depth_index_t> *depth_indexes, tfx_depth_index_t depth_index);
tfxINTERNAL void tfx__reset_particle_effect_flags(tfx_effect_manager pm);
tfxINTERNAL void tfx__free_compute_slot(tfx_effect_manager pm, unsigned int slot_id);
tfxINTERNAL tfxEffectID tfx__add_effect_to_effect_manager(tfx_effect_manager pm, tfx_effect_descriptor effect, int buffer, tfxU32 root_effect_index, float add_delayed_spawning);
tfxINTERNAL void tfx__free_particle_list(tfx_effect_manager pm, tfxU32 index);
tfxINTERNAL void tfx__free_spawn_location_list(tfx_effect_manager pm, tfxU32 index);
tfxINTERNAL void tfx__free_all_particle_lists(tfx_effect_manager pm);
tfxINTERNAL void tfx__free_all_spawn_location_lists(tfx_effect_manager pm);
tfxINTERNAL void tfx__order_effect_sprites(tfx_effect_instance_data_t *sprites, tfxU32 layer, tfx_effect_manager pm);

//Compute stuff doesn't work currently. Keeping this here for now for when I get back to implementing compute shaders for TimelineFX
tfxINTERNAL void tfx__enable_compute(tfx_effect_manager pm) { pm->flags |= tfxEffectManagerFlags_use_compute_shader; }
tfxINTERNAL void tfx__disable_compute(tfx_effect_manager pm) { pm->flags &= ~tfxEffectManagerFlags_use_compute_shader; }
tfxINTERNAL int tfx__add_compute_controller(tfx_effect_manager pm);
tfxINTERNAL tfx_compute_particle_t *tfx__grab_compute_particle(tfx_effect_manager pm, unsigned int layer);
tfxINTERNAL void tfx__reset_particle_ptr(tfx_effect_manager pm, void *ptr);
tfxINTERNAL void tfx__reset_controller_ptr(tfx_effect_manager pm, void *ptr);
tfxINTERNAL void tfx__update_compute(tfx_effect_manager pm, void *sampled_particles, unsigned int sample_size = 100);
tfxINTERNAL void tfx__init_common_effect_manager(tfx_effect_manager pm, tfxU32 max_particles, unsigned int effects_limit, bool double_buffered_sprites, bool dynamic_sprite_allocation, bool group_sprites_by_effect, tfxU32 mt_batch_size);
tfxINTERNAL bool tfx__valid_effect_id(tfx_effect_manager pm, tfxEffectID id);

//--------------------------------
//Effect templates
//--------------------------------
tfxINTERNAL void tfx__add_template_path(tfx_effect_template effect_template, tfx_effect_descriptor effect_emitter, const char *path);

//--------------------------------
//Library functions, internal/Editor functions
//--------------------------------
tfxINTERNAL void tfx__prepare_library_effect_template_path(tfx_library library, const char *path, tfx_effect_template effect);
tfxINTERNAL void tfx__reset_sprite_data_lerp_offset(tfx_sprite_data_t *sprites);
tfxINTERNAL void tfx__reset_ribbon_data_lerp_offset(tfx_sprite_data_t *sprites);
tfxINTERNAL void tfx__compress_sprite_data(tfx_effect_manager pm, tfx_effect_descriptor effect, float frame_length, int *progress);
tfxINTERNAL void tfx__compress_ribbon_data(tfx_sprite_data_t *sprite_data, tfx_sprite_data_settings_t *settings, float frame_length, int *progress);
tfxINTERNAL void tfx__update_sprite_alignment_data(tfx_sprite_data_t *sprite_data, float update_time);
tfxINTERNAL void tfx__link_up_sprite_captured_indexes(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__build_all_library_paths(tfx_library library);
tfxINTERNAL tfx_str64_t tfx__get_name_from_path(const char *path);
tfxAPI_EDITOR bool tfx__is_root_effect(tfx_effect_descriptor effect);
tfxINTERNAL void tfx__reset_effect_parents(tfx_effect_descriptor effect);
tfxINTERNAL float tfx__get_effect_loop_length(tfx_effect_descriptor effect);

//------------------------------------------------------------
//Section API_Functions
//------------------------------------------------------------

#endif		//__cpluscplus

tfxAPI void tfx_UpdateAnimationManagerBufferMetrics(tfx_animation_manager animation_manager);
tfxAPI tfx_storage_t *tfx_GetGlobals();
tfxAPI tfx_pool_stats_t tfx_CreateMemorySnapshot(tfx_header *first_block);
tfxAPI float tfx_DegreesToRadians(float degrees);
tfxAPI float tfx_RadiansToDegrees(float radians);

//--------------------------------
//Random numbers
//--------------------------------
tfxAPI tfx_random_t tfx_NewRandom(tfxU32 seed);
tfxAPI void tfx_AdvanceRandom(tfx_random_t *random);
tfxAPI void tfx_RandomReSeedTime(tfx_random_t *random);
tfxAPI void tfx_RandomReSeed2(tfx_random_t *random, tfxU64 seed1, tfxU64 seed2);
tfxAPI void tfx_RandomReSeed(tfx_random_t *random, tfxU64 seed);
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

//--------------------------------
//Initialisation_functions
//--------------------------------

/*
You don't have to call this, you can just call tfx_InitialiseTimelineFX in order to initialise the memory, but I created this for the sake of the editor which
needs to load in an ini file before initialising timelinefx which requires the memory pool to be created before hand
* @param memory_pool_size    The size of each memory pool to contain all objects created in TimelineFX, recommended to be at least 64MB
*/
tfxAPI void tfx_InitialiseTimelineFXMemory(size_t memory_pool_size);

/*
Initialise TimelineFX. Must be called before any functionality of TimelineFX is used.
* @param max_threads        Pass the number of threads that you want to use in addition to the main thread.
*                            Example, if there are 12 logical cores available, 0.5 will use 6 threads. 0 means only single threaded will be used.
* @param memory_pool_size    The size of each memory pool to contain all objects created in TimelineFX, recommended to be at least 64MB
*/
tfxAPI void tfx_InitialiseTimelineFX(int max_threads, size_t memory_pool_size);

/*
Cleanup up all threads and memory used by timelinefx
*/
tfxAPI void tfx_EndTimelineFX();

//--------------------------------
//Global_variable_access
//--------------------------------
/*
Set the color format used for storing color ramps. Color ramps are generated by particle emitters and dictate how the particle colors change over the lifetime
of the particle. They can be uploaded to the GPU so you can set your preference for the color format as you need. The format should be immediately after you
call tfx_InitialiseTimelineFX
* @param color_format        Pass the tfx_color_ramp_format, either tfx_color_ramp_format_rgba or tfx_color_ramp_format_bgra
*/
tfxAPI void tfx_SetColorRampFormat(tfx_color_ramp_format color_format);

//--------------------------------
//Library_functions
//--------------------------------

/*
Create a new library handle
* @param filename        The name of the file where you want to count the number of shapes
* @returns int           The number of shapes in the library.
*/
tfxAPI tfx_library tfx_CreateLibrary();

/*
Initialise TimelineFX. Must be called before any functionality of TimelineFX is used.
* @param filename        The name of the file where you want to count the number of shapes
* @returns int            The number of shapes in the library.
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
tfxAPI tfx_library tfx_LoadEffectLibrary(const char *filename, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset), void *user_data);

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
tfxAPI tfx_library tfx_LoadEffectLibraryFromMemory(const void *data, tfxU32 size, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset), void *user_data);

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
tfxAPI tfxErrorFlags tfx_LoadSpriteData(const char *filename, tfx_animation_manager animation_manager, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void *user_data);

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
TFX_DISABLE_COMPILER_WARNING("-Wreturn-type-c-linkage")
tfxAPI tfx_image_data_t tfx_GetLibraryImage(tfx_library library, tfxU32 index);
TFX_ENABLE_COMPILER_WARNING()

/*
Output all the effect names in a library to the console
* @param tfx_library                A valid pointer to a tfx_library
*/
tfxAPI void ListEffectNames(tfx_library library);

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
tfxAPI void tfx_BuildLibraryGPUShapeData(tfx_library library, tfx_gpu_shapes shapes, void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset));

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
Get the gpu graph lookup data pointer in library. This lookup data can be used by shaders to update attributes of particles and ribbons (currently ribbons only). Use with 
tfx_GetLibraryGPUGraphLookupsBufferSizeInBytes to upload to your GPU buffer.
* @param library        A handle to a tfx_library
*/
tfxAPI tfx_gpu_graph_data_t *tfx_GetGPUGraphLookupsBuffer();

/*
Get the gpu graph lookup data size in bytes contained within the library you pass into the function.
* @param library        A handle to a tfx_library
*/
tfxAPI tfxU32 tfx_GetGPUGraphLookupsBufferSizeInBytes();

//--------------------------------
//Particle_Manager_functions
//--------------------------------
/*
Create a tfx_effect_manager_info_t object which contains configuration data that you can pass to tfx_CreateEffectManager to setup a effect manager. You can tweak the config after calling this
function if needed to fine tune the settings.
* @param setup                    A tfx_effect_manager_setup enum which you can use to set the info based on some commonly used templates
*/
tfxAPI tfx_effect_manager_info_t tfx_CreateEffectManagerInfo(tfx_effect_manager_setup setup);

/*
Initialize a effect manager with a tfx_effect_manager_info_t object which contains setup data for how to configure the effect manager. See tfx_CreateEffectManagerInfo
* @param pm						A pointer to an unitialised tfx_effect_manager_t. If you want to reconfigure a effect manager for a different usage then you can call tfx_ReconfigureEffectManager.
* @param library                A pointer to a tfx_library_t that you will be using to add all of the effects from to the effect manager.
* @param info                   A tfx_effect_manager_info_t pointer containing the configuration for the effect manager.
*/
tfxAPI tfx_effect_manager tfx_CreateEffectManager(tfx_effect_manager_info_t info);

/*
Reconfigure a effect manager to make it work in a different mode. A effect manager can only run in a single mode at time like unordered, depth ordered etc so use this to change that. Also bear
in mind that you can just use more than one effect manager and utilised different modes that way as well. The modes that you need will depend on the effects that you're adding to the effect manager.
* @param pm                       A pointer to an intialised tfx_effect_manager_t.
* @param mode                     One of the following modes:
								  tfxEffectManagerMode_unordered
								  tfxEffectManagerMode_ordered_by_age
								  tfxEffectManagerMode_ordered_by_depth
								  tfxEffectManagerMode_ordered_by_depth_guaranteed
* @param sort_passes              The number of sort passes if you're using depth sorted effects
*/
tfxAPI void tfx_ReconfigureEffectManager(tfx_effect_manager pm, tfxU32 sort_passes);

tfxAPI void tfx_CompleteEffectManagerWork(tfx_effect_manager pm);

/*
Set the staging buffer used in the effect manager. The effect manager flags must be set with tfxEffectManagerFlags_direct_to_staging_buffer when the particle
manager was created. Depending on the renderer you use you may have to call this before each time you update the effect manager so you can set the buffer to the
current frame in flight. This will probably apply in any modern renderer like vulkan, metal or dx12.
Note: It's up to you to ensure that the staging buffer has enough capacity. The effect manager will assume that the size_in_bytes that you pass to the particle
manager is correct and if tfxEffectManagerFlags_dynamic_sprite_allocation is set will attempt to grow the buffer by calling the callback you set to do this.
* @param pm                       A pointer to an intialised tfx_effect_manager_t.
* @param staging_buffer           A pointer to the staging buffer where all the instance_data/billboards are written to
* @param size_in_bytes            The size in bytes of the staging buffer
*/
tfxAPI void tfx_SetStagingBuffer(tfx_effect_manager pm, void *staging_buffer, tfxU32 size_in_bytes);

/*
Turn on and off whether the effect manager should sort the effects by depth order. Use tfx_SetPMCamera to set the position of the camera that the effect manager will
use to update the depth of each effect in the scene.
* @param pm                       A pointer to an intialised tfx_effect_manager_t.
* @param yesno                    A boolean, set to true or false if you want auto ordering on or off respectively
*/
tfxAPI void tfx_TogglePMOrderEffects(tfx_effect_manager pm, bool yesno);

/*
Get the billboard buffer in the effect manager containing all the sprite instances that were created in the most recent frame. You can use this to copy to a staging buffer to upload to the gpu.
* @param pm                       A pointer to an intialised tfx_effect_manager_t.
*/
tfxAPI tfx_instance_t *tfx_GetInstanceBuffer(tfx_effect_manager  pm);

/*
Get the number of instances within the instance buffer of a effect manager
* @param pm                       A pointer to an intialised tfx_effect_manager_t.
*/
tfxAPI int tfx_GetInstanceCount(tfx_effect_manager pm);

/*
Get the update time being used by the effect manager.
* @param pm                       A handle to an intialised tfx_effect_manager_t.
*/
tfxAPI float tfx_GetUpdateTime(tfx_effect_manager pm);

/*
Get the ribbon buffer for a given segment size. This will give you all the necessary info and buffer pointers for uploading the ribbon data to the GPU for processing and converting into
a vertex buffer for rendering.
* @param pm                       A pointer to an intialised tfx_effect_manager_t.
* @param segment_count            An unsigned int specifying the ribbon length that you want the rendering info for.
* @returns						  A pointer to a tfx_ribbon_bucket_t
*/
tfxAPI tfx_ribbon_bucket_t *tfx_GetRibbonBuffers(tfx_effect_manager pm, tfxKey bucket_id);

/*
Call this to determine whether or not any effect manager has ribbon_emitters to draw this frame.
* @returns						  True or false
*/
tfxAPI bool tfx_HasRibbonsToDraw();

/*
Get a struct containing the info you need to compute and render ribbon_emitters of a specific length. 
* @param pm                       A pointer to an intialised tfx_effect_manager_t.
* @param segment_count            An unsigned int specifying the ribbon length that you want the rendering info for.
* @returns						  A tfx_ribbon_buffer_info_t struct
*/
tfxAPI tfx_ribbon_buffer_info_t tfx_GetRibbonBufferInfo(tfx_effect_manager pm, tfxKey bucket_id);

tfxAPI tfx_ribbon_dispatch_t tfx_CreateRibbonDispatch();

tfxAPI bool tfx_NextRibbonDispatch(tfx_ribbon_dispatch_t *ribbon_dispatch);

tfxAPI void tfx_ResetRibbonDispatchIterator(tfx_effect_manager pm);

tfxAPI tfx_ribbon_buffer_requirements_t tfx_GetRibbonBufferRequirements();

tfxAPI void tfx_CopyRibbonDataToStagingBuffers(void *segments_dst, void *ribbons_dst, void *emitters_dst);

tfxAPI size_t tfx_GetSegmentBufferMaxSizeInBytes(tfx_effect_manager pm);

tfxAPI size_t tfx_GetSegmentVertexBufferMaxSizeInBytes(tfx_effect_manager pm, tfxU32 vertex_size);

tfxAPI size_t tfx_GetSegmentIndexBufferMaxSizeInBytes(tfx_effect_manager pm);

tfxAPI size_t tfx_GetRibbonBufferMaxSizeInBytes(tfx_effect_manager pm, tfxU32 max_ribbons);

tfxAPI size_t tfx_GetEmitterBufferMaxSizeInBytes(tfx_effect_manager pm);

/*
When a effect manager updates particles it creates work queues to handle the work. By default these each have a maximum amount of 1000 entries which should be
more than enough for most situations. However you can increase the sizes here if needed. You only need to set this manually if you hit one of the asserts when these
run out of space or you anticipate a huge amount of emitters and particles to be used (> million). On the other hand, you might be tight on memory in which case you
could reduce the numbers as well if needed (they don't take a lot of space though)
* @param pm                        A pointer to an intialised tfx_effect_manager_t.
* @param spawn_work_max            The maximum amount of spawn work entries
* @param control_work_max          The maximum amount of control work entries
* @param age_work_max              The maximum amount of age_work work entries
*/
tfxAPI void tfx_SetPMWorkQueueSizes(tfx_effect_manager pm, tfxU32 spawn_work_max, tfxU32 control_work_max, tfxU32 age_work_max);

/*Free the memory for a specific emitter type. When an emitter is created it creates memory to store all of the particles that it updates each frame. If you have
multiple emitters of the same type then their particle lists are resused rather then freed as they expire. When they're freed then the unused list is added to a list
of free particle banks for that emitter type so that they can then be recycled if another emitter of the same type is created. If you want to free the memory for a
specific emitter then you can call this function to do that.
NOTE: No emitters of the type passed to the function must be in use in the effect manager.
* @param pm                        A pointer to an intialised tfx_effect_manager_t.
* @param emitter                   A pointer to a valid tfx_effect_descriptor_t of type tfxEmitterType
*/
tfxAPI void tfx_FreeParticleListsMemory(tfx_effect_manager pm, tfx_effect_descriptor emitter);

/*
Free all the memory that is associated with an effect. Depending on the configuration of the effect manager this might be instance_data, particle lists and spawn location lists.
* @param pm                        A pointer to an intialised tfx_effect_manager_t.
* @param emitter                   A pointer to a valid tfx_effect_descriptor_t of type tfxEffectType
*/
tfxAPI void tfx_FreeEffectListsMemory(tfx_effect_manager pm, tfx_effect_descriptor effect);

/*
Get the current particle count for a effect manager
* @param pm                        A pointer to an tfx_effect_manager_t
* @returns tfxU32                  The total number of particles currently being updated
*/
tfxAPI tfxU32 tfx_ParticleCount(tfx_effect_manager pm);

/*
Get the current ribbon count for a effect manager
* @param pm                        A pointer to an tfx_effect_manager_t
* @returns tfxU32                  The total number of particles currently being updated
*/
tfxAPI tfxU32 tfx_RibbonCount(tfx_effect_manager pm);

/*
Get the current number of effects that are currently being updated by a effect manager
* @param pm                        A pointer to an tfx_effect_manager_t
* @returns tfxU32                  The total number of effects currently being updated
*/
tfxAPI tfxU32 tfx_EffectCount(tfx_effect_manager pm);

/*
Get the current number of emitters that are currently being updated by a effect manager
* @param pm                        A pointer to an tfx_effect_manager_t
* @returns tfxU32                  The total number of emitters currently being updated
*/
tfxAPI tfxU32 tfx_EmitterCount(tfx_effect_manager pm);

/*
Set the seed for the effect manager for random number generation. Setting the seed can determine how an emitters spawns particles, so if you set the seed before adding an effect to the effect manager
then the effect will look the same each time. Note that seed of 0 is invalid, it must be 1 or greater.
* @param pm                        A pointer to an initialised tfx_effect_manager_t. 
* @param seed                      An unsigned int representing the seed (Any value other then 0)
*/
tfxAPI void tfx_SetSeed(tfx_effect_manager pm, tfxU64 seed);

/*
Add an effect to a tfx_effect_manager_t from an effect template
* @param pm                         A pointer to an initialised tfx_effect_manager_t.
* @param effect_template			The tfx_effect_template_t object that you want to add to the effect manager. It must have already been prepared by calling tfx_PrepareEffectTemplate
* @param effect_id					pointer to a tfxEffectID of the effect which will be set after it's been added to the effect manager. This index can then be used to manipulate the effect in the effect manager as it's update
									For example by calling tfx_SetEffectPosition. This will be set to tfxINVALID if the function is unable to add the effect to the effect manager if it's out of space and reached it's effect limit.
  @returns							True if the effect was succesfully added.
*/
tfxAPI bool tfx_AddEffectTemplateToEffectManager(tfx_effect_manager pm, tfx_effect_template effect, tfxEffectID *effect_id);

/*
Add an effect to a tfx_effect_manager_t. Generally you should always call tfx_AddEffectTemplateToEffectManager and use templates to organise your effects but if you want to just
test things out you can add an effect direct from a library using this command.
* @param pm							A pointer to an initialised tfx_effect_manager_t. 
* @param effect						tfx_effect_descriptor_t object that you want to add to the effect manager.
* @param effect_id					pointer to a tfxEffectID of the effect which will be set after it's been added to the effect manager. This index can then be used to manipulate the effect in the effect manager as it's update
									For example by calling tfx_SetEffectPosition. This will be set to tfxINVALID if the function is unable to add the effect to the effect manager if it's out of space and reached it's effect limit.
  @returns							True if the effect was succesfully added.
*/
tfxAPI bool tfx_AddRawEffectToEffectManager(tfx_effect_manager pm, tfx_effect_descriptor effect, tfxEffectID *effect_id);

/*
Update a effect manager. If you are interpolating particles in the vertex shader then it's important to only call this function once per frame only and idealy in a fixed step loop.
That means that if your fixed loop has to run twice to catch up (because of low frame rates) then you should still only call this function once but you can multiply the elapsed time
by the number of ticks. The ellapsed time should be the amount of time that has passed since the last frame so in a fixed step loop this will simply be the update rate in millisecs.
For example if you're updating 60 frames per second then elapsed time would be 16.666667. Psuedo code would look something like this:

	TimerAccumulate(game->timer);
	int pending_ticks = TimerPendingTicks(game->timer);	//The number of times the update loop will run this frame.

	while (tfx_TimerDoUpdate(game->timer)) {
		if (pending_ticks > 0) {
			tfx_UpdateEffectManager(&game->pm, FrameLength * pending_ticks);
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

* @param pm                    A pointer to an initialised tfx_effect_manager_t.
*/
tfxAPI void tfx_UpdateEffectManager(tfx_effect_manager pm, float elapsed);

/*
Get the image pointer for a sprite. Use this when rendering particles in your renderer. The pointer that is returned will be the pointer that you set in your shape loader function
used when loading an effect library. Generally you shouldn't need to use this function, simply copy the whole instance buffer in the effect manager to your staging buffer to be
copied to the gpu in one go.
* @param pm                    A pointer to an initialised tfx_effect_manager_t. 
* @param property_indexes    The value in the instance_data->property_indexs[i] when iterating over the instance_data in your render function
  @returns                    void* pointer to the image
*/
tfxAPI void *tfx_GetSpriteImagePointer(tfx_effect_manager pm, tfxU32 property_indexes);

/*
Get the handle of the sprite. Use this when rendering particles in your renderer one sprite at a time.
* @param pm                    A pointer to an initialised tfx_effect_manager_t.
* @param property_indexes      The value in the instance_data->property_indexs[i] when iterating over the instance_data in your render function
  @out_handle                  Pass in a pointer to a vec2 which will be loaded with the handle values
*/
tfxAPI void tfx_GetSpriteHandle(void *instance, float out_handle[2]);

/*
Get the total number of instances ready for rendering in the effect manager.
* @param pm                    A pointer to an initialised tfx_effect_manager_t.
*/
tfxAPI tfxU32 tfx_TotalSpriteCount(tfx_effect_manager pm);

/*
Clear all particles, instance_data and effects in a effect manager. If you don't need to use the effect manager again then call tfx_FreeEffectManager to also
free all the memory associated with the effect manager.
* @param pm                        A pointer to an initialised tfx_effect_manager_t.
* @param free_particle_banks    Set to true if you want to free the memory associated with the particle banks and release back to the memory pool
*/
tfxAPI void tfx_ClearEffectManager(tfx_effect_manager pm, bool free_particle_banks, bool free_sprite_buffers);

/*
Free all the memory used in the effect manager.
* @param pm                        A pointer to an initialised tfx_effect_manager_t.
*/
tfxAPI void tfx_FreeEffectManager(tfx_effect_manager pm);

//[Effects functions for altering effects that are currently playing out in a effect manager]

/*
Expire an effect by telling it to stop spawning particles. This means that the effect will eventually be removed from the effect manager after all of it's remaining particles have expired.
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
*/
tfxAPI void tfx_SoftExpireEffect(tfx_effect_manager pm, tfxEffectID effect_index);

/*
Soft expire all the effects in a effect manager so that the particles complete their animation first
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed
*/
tfxAPI void tfx_SoftExpireAll(tfx_effect_manager pm);

/*
Expire an effect by telling it to stop spawning particles and remove all associated particles immediately.
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
*/
tfxAPI void tfx_HardExpireEffect(tfx_effect_manager pm, tfxEffectID effect_index);

/*
Get effect user data
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @returns                void* pointing to the user data set in the effect. See tfx_effect_template_t::SetUserData() and tfx__set_effect_user_data()
*/
tfxAPI void *tfx_GetEffectUserData(tfx_effect_manager pm, tfxEffectID effect_index);

/*
More for use in the editor, this function updates emitter base values for any effects that are currently running after their graph values have been changed.
*/
tfxAPI void tfx_UpdatePMBaseValues(tfx_effect_manager pm);

/*
Set the effect manager camera. This is used to calculate particle depth if you're using depth ordered particles so it needs to be updated each frame.
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed
* @param front            An array of 3 floats representing a normalised 3d vector describing the direction that the camera is pointing
* @param position        An array of 3 floats representing the position of the camera in 3d space
*/
tfxAPI void tfx_SetPMCamera(tfx_effect_manager pm, float front[3], float position[3]);

/*
Each effect in the effect manager can have bounding box which you can decide to keep updated or not if you wanted to do any offscreen culling of effects. Theres some
extra overhead to keep the bounding boxes updated but that can be made back if you have a number of effect particles offscreen that don't need to be drawn.
* @param pm					A pointer to a tfx_effect_manager_t where the effect is being managed
* @param yesno				Set to true or false if you want the bounding boxes to be udpated.
*/
tfxAPI void tfx_KeepBoundingBoxesUpdated(tfx_effect_manager pm, bool yesno);

/*
Set the effect user data for an effect already added to a effect manager
* @param pm					A pointer to a tfx_effect_manager_t where the effect is being managed
* @param effect_index		The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param user_data			A void* pointing to the user_data that you want to store in the effect
*/
tfxAPI void tfx_SetEffectUserData(tfx_effect_manager pm, tfxU32 effect_index, void *data);

/*
Force a effect manager to only run in single threaded mode. In other words, only use the main thread to update particles
* @param pm                A pointer to a tfx_effect_manager_t.
* @param switch_on        true or false to use a single thread or not
*/
tfxAPI void tfx_ForcePMSingleThreaded(tfx_effect_manager pm, bool switch_on);

/*
Get the transform vectors for a sprite's previous position so that you can use that to interpolate between that and the current sprite position
* @param pm                A pointer to a tfx_effect_manager_t.
* @param layer            The index of the sprite layer
* @param index            The sprite index of the sprite that you want the captured sprite for.
* @param position         This should be a pointer to a vec3 that you pass in that will get loaded with the position of the instance
*/
tfxAPI void tfx_GetCapturedInstanceTransform(tfx_effect_manager pm, tfxU32 layer, tfxU32 index, float out_position[3]);

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
* @param pm                A pointer to a tfx_effect_manager_t.
* @param yesno            True = disable spawning, false = enable spawning
*/
tfxAPI void tfx_DisablePMSpawning(tfx_effect_manager pm, bool yesno);

/*
Get the buffer of effect indexes in the effect manager.
* @param pm               A pointer to a tfx_effect_manager_t.
* @param depth            The depth of the list that you want. 0 are top level effects and anything higher are sub effects within those effects
* @param count			  A pointer to an int that you can pass in that will be filled with the count of effects in the array
* @returns                Pointer to the array of effect indexes
*/
tfxAPI tfx_effect_index_t *tfx_GetPMEffectBuffer(tfx_effect_manager pm, int *count);

/*
Get the buffer of emitter indexes in the effect manager.
* @param pm                A pointer to a tfx_effect_manager_t.
* @param depth            The depth of the list that you want. 0 are top level emitters and anything higher are sub emitters within those effects
* @param count			  A pointer to an int that you can pass in that will be filled with the count of emitters in the array
* @returns                Pointer to the tfxvec of effect indexes
*/
tfxAPI tfxU32 *tfx_GetPMEmitterBuffer(tfx_effect_manager pm, int *count);

/*
Set the position of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param x                The x value of the position
* @param y                The y value of the position
* @param z                The z value of the position
*/
tfxAPI void tfx_SetEffectPosition(tfx_effect_manager pm, tfxEffectID effect_index, float x, float y, float z);

/*
Set the position of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param position        A tfx_vec3_t vector object containing the x, y and z coordinates
*/
tfxAPI void tfx_SetEffectPositionVec3(tfx_effect_manager pm, tfxEffectID effect_index, tfx_vec3_t position);

/*
Move an Effect by a specified amount relative to the effect's current position
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param amount            A tfx_vec3_t vector object containing the amount to move in the x, y and z planes
*/
tfxAPI void tfx_MoveEffectVec3(tfx_effect_manager pm, tfxEffectID effect_index, tfx_vec3_t amount);

/*
Move an Effect by a specified amount relative to the effect's current position
* @param pm               A pointer to a tfx_effect_manager_t where the effect is being managed
* @param effect_index     The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param x                The amount to move in the x plane
* @param y                The amount to move in the y plane
* @param z                The amount to move in the z plane
*/
tfxAPI void tfx_MoveEffect(tfx_effect_manager pm, tfxEffectID effect_index, float x, float y, float z);

/*
Get the current position of an effect
* @param pm              A pointer to a tfx_effect_manager_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @return                tfx_vec3_t containing the effect position
*/
tfxAPI void tfx_GetEffectPositionVec3(tfx_effect_manager pm, tfxEffectID effect_index, float out_position[3]);

/*
You can use this function to get the billboard buffer of a specific effect. 
* @param pm						A pointer to a tfx_effect_manager_t where the effect is being managed
* @param effect_index			The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param tfxU32					Pass in a pointer to a tfxU32 which will be set to the number of instance_data in the buffer.
* @return						tfx_instance_t pointer to the buffer
*/
tfxAPI tfx_instance_t *tfx_GetEffectInstanceBuffer(tfx_effect_manager pm, tfxEffectID effect_index, tfxU32 *sprite_count);

/*
You can use this function to get each billboard buffer for every effect that is currently active in the effect manager. Generally you would call this inside a for loop for each layer.
* @param pm						A pointer to a tfx_effect_manager_t where the effect is being managed
* @param tfx_sprite_billboard_t	Pass in a pointer which will be set to the current sprite buffer containing all of the sprite data for this frame.
* @param tfx_effect_instance_data_t   Pass in a second pointer which will be set to the tfx_effect_instance_data_t containing all of the sprite buffer data. This can be used to gain access to all the sprite data if using double buffered instance_data (to interpolated with the previous frame).
* @param tfxU32					Pass in a pointer to a tfxU32 which will be set to the number of instance_data in the buffer.
* @return						true or false if the next billboard buffer was found. False will be returned once there are no more effect sprite buffers in the effect manager
*/
tfxAPI bool tfx_GetNextInstanceBuffer(tfx_effect_manager pm, tfx_instance_t **sprites_soa, tfx_effect_instance_data_t **effect_sprites, tfxU32 *sprite_count);

/*After calling GetNextBillboard/SpriteBuffer in a while loop you can call this to reset the index for the next frame
* @param pm						A pointer to a tfx_effect_manager_t
*/
tfxAPI void tfx_ResetInstanceBufferLoopIndex(tfx_effect_manager pm);

/*
Set the roll of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current roll of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param roll            A float of the amount that you want to set the roll too
*/
tfxAPI void tfx_SetEffectRoll(tfx_effect_manager pm, tfxEffectID effect_index, float roll);

/*
Set the pitch of a effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current pitch of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param pitch            A float of the amount that you want to set the pitch too
*/
tfxAPI void tfx_SetEffectPitch(tfx_effect_manager pm, tfxEffectID effect_index, float pitch);

/*
Set the yaw of a effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current yaw of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param yaw            A float of the amount that you want to set the yaw too
*/
tfxAPI void tfx_SetEffectYaw(tfx_effect_manager pm, tfxEffectID effect_index, float yaw);

/*
Set the width of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current width of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param width            A float of the amount that you want to set the width multiplier too. The width multiplier will multiply all widths of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void tfx_SetEffectWidthMultiplier(tfx_effect_manager pm, tfxEffectID effect_index, float width);

/*
Set the height of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current height of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param height            A float of the amount that you want to set the height multiplier too. The height multiplier will multiply all heights of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void tfx_SetEffectHeightMultiplier(tfx_effect_manager pm, tfxEffectID effect_index, float height);

/*
Set the depth of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current depth of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param depth            A float of the amount that you want to set the depth multiplier too. The depth multiplier will multiply all heights of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void tfx_SetEffectDepthMultiplier(tfx_effect_manager pm, tfxEffectID effect_index, float depth);

/*
Set the life multiplier of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current life of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param life            A float of the amount that you want to set the life multiplier too. The life mulitplier will affect how long all particles emitted within the effect will last before expiring.
*/
tfxAPI void tfx_SetEffectLifeMultiplier(tfx_effect_manager pm, tfxEffectID effect_index, float life);

/*
Set the particle width multiplier of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current particle width of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param width            A float of the amount that you want to set the particle width multiplier too. The particle width mulitplier will affect the width of each particle if the emitter has a non uniform particle size, otherwise
						it will uniformly size the particle
*/
tfxAPI void tfx_SetEffectParticleWidthMultiplier(tfx_effect_manager pm, tfxEffectID effect_index, float width);

/*
Set the particle height multiplier of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current particle width of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param height            A float of the amount that you want to set the particle height multiplier too. The particle height mulitplier will affect the height of each particle if the emitter has a non uniform particle size, otherwise
						this function will have no effect.
*/
tfxAPI void tfx_SetEffectParticleHeightMultiplier(tfx_effect_manager pm, tfxEffectID effect_index, float height);

/*
Set the velocity multiplier of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current velocity of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param velocity        A float of the amount that you want to set the particle velocity multiplier too. The particle velocity mulitplier will affect the base velocity of a particle at spawn time.
*/
tfxAPI void tfx_SetEffectVelocityMultiplier(tfx_effect_manager pm, tfxEffectID effect_index, float velocity);

/*
Set the spin multiplier of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current spin of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param spin            A float of the amount that you want to set the particle spin multiplier too. The particle spin mulitplier will affect the base spin of a particle at spawn time.
*/
tfxAPI void tfx_SetEffectSpinMultiplier(tfx_effect_manager pm, tfxEffectID effect_index, float spin);

/*
Set the intensity multiplier of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current intensity of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param intensity        A float of the amount that you want to set the particle intensity multiplier too. The particle intensity mulitplier will instantly affect the opacity of all particles currently emitted by the effect.
*/
tfxAPI void tfx_SetEffectIntensityMultiplier(tfx_effect_manager pm, tfxEffectID effect_index, float intensity);

/*
Set the splatter multiplier of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current splatter of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param splatter        A float of the amount that you want to set the particle splatter multiplier too. The particle splatter mulitplier will change the amount of random offset all particles emitted in the effect will have.
*/
tfxAPI void tfx_SetEffectSplatterMultiplier(tfx_effect_manager pm, tfxEffectID effect_index, float splatter);

/*
Set the weight multiplier of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current weight of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param weight            A float of the amount that you want to set the particle weight multiplier too. The particle weight mulitplier will change the weight applied to particles in the effect at spawn time.
*/
tfxAPI void tfx_SetEffectWeightMultiplier(tfx_effect_manager pm, tfxEffectID effect_index, float weight);

/*
Set the overal scale of an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateEffectManager in order to override the current weight of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param overal_scale    A float of the amount that you want to set the overal scale to. The overal scale is an simply way to change the size of an effect
*/
tfxAPI void tfx_SetEffectOveralScale(tfx_effect_manager pm, tfxEffectID effect_index, float overal_scale);

/*
Set the base noise offset for an effect
* @param pm                A pointer to a tfx_effect_manager_t where the effect is being managed.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToEffectManager
* @param noise_offset    A float of the amount that you want to set the effect noise offset to. By default when an effect is added to a effect manager a random noise offset will be set based on the Base Noise Offset Range property. Here you can override that
						value by setting it here. The most ideal time to set this would be immediately after you have added the effect to the effect manager, but you could call it any time you wanted for a constantly changing noise offset.
*/
tfxAPI void tfx_SetEffectBaseNoiseOffset(tfx_effect_manager pm, tfxEffectID effect_index, float noise_offset);

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
* @param position                A tfx_vec3_t vector object containing the x, y and z coordinates
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
tfxAPI void tfx_SetAnimationManagerInstanceCallback(tfx_animation_manager animation_manager, bool((*maybe_render_instance_callback)(tfx_animation_manager animation_manager, tfx_animation_instance_t *instance, tfx_frame_meta_t *meta, void *user_data)));

/*
Get the sprite data settings for an effect in a library. Sprite data settings are the settings for an effect in the editor relating to setting up pre-baked effects
* @param library				Pointer to the tfx_library_t where the effect is stored
* @param effect					Pointer the the effect that you want the sprite settings for.
* @returns						Pointer to the tfx_sprite_data_settings
*/
tfx_sprite_data_settings_t *tfx_GetEffectSpriteDataSettings(tfx_library library, tfx_effect_descriptor effect);

/*
Get the sprite data settings for an effect in a library by it's path. Sprite data settings are the settings for an effect in the editor relating to setting up pre-baked effects
* @param library				Pointer to the tfx_library_t where the effect is stored
* @param path					const char* string of the path to the effect. Must be the path to a root effect.
* @returns						Pointer to the tfx_sprite_data_settings
*/
tfx_sprite_data_settings_t *tfx_GetEffectSpriteDataSettingsByPath(tfx_library library, const char *path);

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
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @param position                A tfx_vec3_t vector object containing the x, y and z coordinates
*/
tfxAPI void tfx_AddSpriteData(tfx_animation_manager animation_manager, tfx_effect_descriptor effect, tfx_effect_manager pm, tfx_vec3_t camera_position);

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
Frees all data from the animation manager including sprite data, metrics and instances and also the handle itself
from being drawn
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to reset
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
tfxAPI void tfx_BuildAnimationManagerGPUShapeData(tfx_animation_manager animation_manager, tfx_gpu_shapes shapes, void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset));

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
tfxAPI void tfx_SetTemplateEffectUpdateCallback(tfx_effect_template t, void(*update_callback)(tfx_effect_manager pm, tfxEffectID effect_index));

/*
Pre-record this effect into a sprite cache so that you can play the effect back without the need to actually caclulate particles in realtime.
	* @param pm					  Reference to a pm that will be used to run the particle simulation and record the sprite data
	* @param path				  const *char of a path to the emitter in the effect.Must be a valid path, for example: "My Effect/My Emitter"
	* @param camera				  Array of 3 floats with the camera position (only needed for effects that are sorted by depth)
*/
tfxAPI void tfx_RecordTemplateEffect(tfx_effect_template t, tfx_effect_manager pm, float update_frequency, float camera_position[3]);

/*
Pre-record this effect into a sprite cache so that you can play the effect back without the need to actually caclulate particles in realtime. This version
of the function allows you to pass in specific settings
	* @param pm					  Reference to a pm that will be used to run the particle simulation and record the sprite data
	* @param path				  const *char of a path to the emitter in the effect.Must be a valid path, for example: "My Effect/My Emitter"
	* @param camera				  Array of 3 floats with the camera position (only needed for effects that are sorted by depth)
*/
tfxAPI void tfx_RecordEffect(tfx_effect_descriptor e, tfx_sprite_data_settings_t *settings, tfx_sprite_data_t *sprite_data, tfx_effect_manager pm, float update_frequency, float camera_position[3]);

/*
Disable an emitter within an effect. Disabling an emitter will stop it being added to the effect manager when calling tfx_AddEffectTemplateToEffectManager
* @param path					  const *char of a path to the emitter in the effect. Must be a valid path, for example: "My Effect/My Emitter"
*/
tfxAPI void tfx_DisableTemplateEmitter(tfx_effect_template t, const char *path);

/*
Enable an emitter within an effect so that it is added to the effect manager when calling tfx_AddEffectTemplateToEffectManager. Emitters are enabled by default.
* @param path					  const *char of a path to the emitter in the effect. Must be a valid path, for example: "My Effect/My Emitter"
*/
tfxAPI void tfx_EnableTemplateEmitter(tfx_effect_template t, const char *path);

/*
Scale all nodes on a global graph graph of the effect
* @param global_type			  tfx_graph_type of the global graph that you want to scale. Must be a global graph or an assert will be called
* @param amount					  A float of the amount that you want to scale the multiplier by.
*/
tfxAPI void tfx_ScaleTemplateGlobalMultiplier(tfx_effect_template t, tfx_global_graph_index graph_index, float amount);

/*
Scale all nodes on an emitter graph
* @param emitter_path			  const *char of the emitter path
* @param global_type			  tfx_graph_type of the emitter graph that you want to scale. Must be an emitter graph or an assert will be called
* @param amount                   A float of the amount that you want to scale the graph by.
*/
tfxAPI void tfx_ScaleTemplateEmitterGraph(tfx_effect_template t, const char *emitter_path, tfx_emitter_graph_index graph_index, float amount);

/*
Set the single spawn amount for an emitter. Only affects emitters that have the single spawn flag set.
* @param emitter_path			 const *char of the emitter path
* @param amount					 A float of the amount that you want to set the single spawn amount to.
*/
tfxAPI void tfx_SetTemplateSingleSpawnAmount(tfx_effect_template t, const char *emitter_path, tfxU32 amount);

//--------------------------------
//General_helpers
//--------------------------------

/*
Interpolate between 2 tfxVec3s. You can make use of this in your render function when rendering instance_data and interpolating between captured and current positions
* @param tween                The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param world                The current tvxVec3 position
* @param captured            The captured tvxVec3 position
* @returns tfx_vec3_t            The interpolated tfx_vec3_t
*/
tfxAPI void tfx_Lerp3d(float lerp, const tfx_vec3_t *world, const tfx_vec3_t *captured, float out_lerp[3]);

/*
Interpolate between 2 tfxVec2s. You can make use of this in your render function when rendering instance_data and interpolating between captured and current positions
* @param tween        The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param world        The current tvxVec2 position
* @param captured    The captured tvxVec2 position
* @returns tfx_vec2_t    The interpolated tfx_vec2_t
*/
tfxAPI void tfx_Lerp2d(float lerp, const tfx_vec2_t *world, const tfx_vec2_t *captured, float out_lerp[2]);

/*
Interpolate between 2 float. You can make use of this in your render function when rendering instance_data and interpolating between captured and current float values like intensity
* @param tween        The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param world        The current tvxVec2 position
* @param captured    The captured tvxVec2 position
* @returns tfx_vec2_t    The interpolated tfx_vec2_t
*/
tfxAPI inline float tfx_LerpFloat(float lerp, const float current, const float captured) {
	return current * lerp + captured * (1.f - lerp);
}

tfxAPI void tfx_GetSpriteScale(void *instance, float out_scale[2]);

#ifdef tfxINTEL

/*
Interpolate all sprite transform data in a single function. This will interpolate position, scale and rotation.
* @param tween                The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param current            The current transform struct of the sprite
* @param captured            The captured transform struct of the sprite
* @returns tfx_wide_lerp_transform_result_t            The interpolated transform data in a tfx_wide_lerp_transform_result_t
*/
tfxAPI inline tfx_wide_lerp_transform_result_t tfx_InterpolateSpriteTransform(const __m128 *tween, const tfx_sprite_transform_t *current, const tfx_sprite_transform_t *captured) {
	tfx128 to1 = _mm_load_ps(&current->position.x);
	tfx128 from1 = _mm_load_ps(&captured->position.x);
	tfx128 to2 = _mm_load_ps(&current->rotations.y);
	tfx128 from2 = _mm_load_ps(&captured->rotations.y);
	tfx128 one_minus_tween = _mm_sub_ps(_mm_set1_ps(1.f), *tween);
	tfx128 to_lerp1 = _mm_mul_ps(to1, *tween);
	tfx128 from_lerp1 = _mm_mul_ps(from1, one_minus_tween);
	tfx128 result = _mm_add_ps(from_lerp1, to_lerp1);
	tfx_wide_lerp_transform_result_t out;
	_mm_store_ps(out.position, result);
	to_lerp1 = _mm_mul_ps(to2, *tween);
	from_lerp1 = _mm_mul_ps(from2, one_minus_tween);
	result = _mm_add_ps(from_lerp1, to_lerp1);
	_mm_store_ps(&out.rotations[1], result);
	return out;
}

#elif defined(tfxARM)

/*
Interpolate all sprite transform data in a single function. This will interpolate position, scale and rotation.
* @param tween                The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param current            The current transform struct of the sprite
* @param captured            The captured transform struct of the sprite
* @returns tfx_wide_lerp_transform_result_t            The interpolated transform data in a tfx_wide_lerp_transform_result_t
*/
tfxAPI inline tfx_wide_lerp_transform_result_t tfx_InterpolateSpriteTransform(const tfx128 *tween, const tfx_sprite_transform_t *current, const tfx_sprite_transform_t *captured) {
	tfx128 to1 = vld1q_f32(&current->position.x);
	tfx128 from1 = vld1q_f32(&captured->position.x);
	tfx128 to2 = vld1q_f32(&current->rotations.y);
	tfx128 from2 = vld1q_f32(&captured->rotations.y);
	tfx128 one_minus_tween = vsubq_f32(vdupq_n_f32(1.f), *tween);
	tfx128 to_lerp1 = vmulq_f32(to1, *tween);
	tfx128 from_lerp1 = vmulq_f32(from1, one_minus_tween);
	tfx128 result = vaddq_f32(from_lerp1, to_lerp1);
	tfx_wide_lerp_transform_result_t out;
	vst1q_f32(out.position, result);
	to_lerp1 = vmulq_f32(to2, *tween);
	from_lerp1 = vmulq_f32(from2, one_minus_tween);
	result = vaddq_f32(from_lerp1, to_lerp1);
	vst1q_f32(&out.rotations[1], result);
	return out;
}

#endif

tfxAPI inline float tfx_GetDistance(float fromx, float fromy, float tox, float toy) {
	float w = tox - fromx;
	float h = toy - fromy;
	return sqrtf(w * w + h * h);
}

tfxAPI tfx_effect_descriptor tfx_NewEffectDescriptor(tfx_effect_descriptor_type type);

/*
Create a new list for containing gpu shapes. You can use this list to upload to the GPU so that shaders can have access to particle image data like UV coordinates
* @returns tfx_gpu_shapes                A handle to the shapes list
*/
tfxAPI tfx_gpu_shapes tfx_CreateGPUShapesList();

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

#endif
