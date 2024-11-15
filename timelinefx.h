#ifndef TFX_LIBRARY_HEADER
#define TFX_LIBRARY_HEADER

#define tfxENABLE_PROFILING
#define tfxPROFILER_SAMPLES 60
#define TFX_THREAD_SAFE
#define TFX_EXTRA_DEBUGGING

/*
	Timeline FX C++ library

	This library is for implementing particle effects into your games and applications.

	This library is render agnostic, so you will have to provide your own means to render the particles. There are various API functions in the library that help you do this.

	Currently tested on Windows and MacOS, Intel and Mac based ARM processors.

	Table of contents
	Sections in this header file, you can search for the following keywords to jump to that section:

	[Zest_Pocket_Allocator]				A single header library for allocating memory from a large pool.
	[Header_Includes_and_Typedefs]		Just your basic header stuff for setting up typedefs and some #defines
	[OS_Specific_Functions]				OS specific multithreading and file access
	[XXHash_Implementation]				XXHasher for the storage map.
	[SIMD_defines]						Defines for SIMD intrinsics
	[Enums]								All the definitions for enums and bit flags
	[Constants]							Various constant definitions
	[String_Buffers]					Basic string buffers for storing names of things in the library and reading from library files.
	[Containers_and_Memory]				Container structs and lists and defines for memory is allocated (uses Zest Pocket Allocator by default)
	[Multithreading_Work_Queues]		Implementation for work queues for threading
	[Vector_Math]						Vec2/3/4 and Matrix2/3/4 structs including wide vectors for SIMD
	[Simplex_Noise]						Some setup for implementing simplex noise.
	[Profiling]							Very basic profiling for internal use
	[File_IO]							A package manager for reading/writing files such as a tfx library effects file
	[Struct_Types]						All of the structs used for objects in TimelineFX
	[Internal_Functions]				Mainly internal functions called only by the library but also the Editor, these are marked either tfxINTERNAL or tfxAPI_EDITOR
	[API_Functions]						The main functions for use by users of the library
		-[Initialisation_functions]		Startup and shutdown timelinefx
		-[Global_variable_access]		Any functions that give you access to global variables relating to timelinefx
		-[Library_functions]			Functions for loading and accessing timelinefx libraries
		-[Particle_Manager_functions]	Create and update functions for particle managers where the main work is done to update particles every frame
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
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>

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

#ifdef __cplusplus
extern "C" {
#endif

#define tfx__MAXIMUM_BLOCK_SIZE (TFX_ONE << TFX_MAX_SIZE_INDEX)

	enum tfx__constants {
		tfx__MEMORY_ALIGNMENT = 1 << MEMORY_ALIGNMENT_LOG2,
		tfx__SECOND_LEVEL_INDEX_LOG2 = 5,
		tfx__FIRST_LEVEL_INDEX_COUNT = TFX_MAX_SIZE_INDEX,
		tfx__SECOND_LEVEL_INDEX_COUNT = 1 << tfx__SECOND_LEVEL_INDEX_LOG2,
		tfx__BLOCK_POINTER_OFFSET = sizeof(void *) + sizeof(tfx_size),
		tfx__MINIMUM_BLOCK_SIZE = 16,
		tfx__BLOCK_SIZE_OVERHEAD = sizeof(tfx_size),
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
		volatile tfxLONG access;
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
#define tfx__strcpy(left, right, size) strcpy(left, right)
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

	//Write functions
#if defined(TFX_THREAD_SAFE)

#define tfx__lock_thread_access(alloc)                                        \
    do {                                                                    \
    } while (0 != tfx__compare_and_exchange(&alloc->access, 1, 0));            \
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
			tfx__unlock_thread_access(allocator);
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
				tfx__unlock_thread_access(allocator);
				return split_block;
			}
		} else {
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
typedef uint32_t tfxU32;
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

#if defined(_WIN32)
#include <SDKDDKVer.h>
#ifndef WIN_LEAN_AND_MEAN
#define WIN_LEAN_AND_MEAN
#endif
#include <Windows.h>
#endif

#define tfxTWO63 0x8000000000000000u 
#define tfxTWO64f (tfxTWO63*2.0)
#define tfxPI 3.14159265359f
#define tfxHALFPI 1.570796f
#define tfxPI2 6.283185307f 
#define tfxINVTWOPI 0.1591549f
#define tfxTHREEHALFPI 4.7123889f
#define tfxQUARTERPI 0.7853982f
#define tfx360Radians 6.28319f
#define tfx180Radians 3.14159f
#define tfx90Radians 1.5708f
#define tfx270Radians 4.71239f
#define tfxMAXDEPTH 3
#define tfxNL u8"\n"
#define tfxPROPERTY_INDEX_MASK 0x00007FFF
#define tfxSPRITE_ALIGNMENT_MASK 0xFF000000
#define tfxSPRITE_IMAGE_FRAME_MASK 0x00FF0000
#define tfxEXTRACT_SPRITE_ALIGNMENT(property_index) ((property_index & tfxSPRITE_ALIGNMENT_MASK) >> 24)
#define tfxEXTRACT_SPRITE_IMAGE_FRAME(property_index) ((property_index & tfxSPRITE_IMAGE_FRAME_MASK) >> 16)
#define tfxEXTRACT_SPRITE_PROPERTY_INDEX(property_index) (property_index & tfxPROPERTY_INDEX_MASK)
#define tfxPACK_SCALE_AND_HANDLE(x, y, lib, property_index) (tfxU16)(x * 127.9960938f) | ((tfxU16)(y * 127.9960938f) << 16) | ((tfxU64)lib->emitter_properties[property_index].image_handle_packed << 32)
#define tfxPACK_SIZE_AND_HANDLE(x, y, lib, property_index) (tfxU16)(x * 7.999755859f) | ((tfxU16)(y * 7.999755859f) << 16) | ((tfxU64)lib->emitter_properties[property_index].image_handle_packed << 32)
#define tfxCIRCLENODES 16
#define tfxPrint(message, ...) printf(message tfxNL, ##__VA_ARGS__)

//----------------------------------------------------------
//Forward declarations

typedef struct tfx_effect_emitter_s tfx_effect_emitter_t;
typedef struct tfx_particle_manager_s tfx_particle_manager_t;
typedef struct tfx_effect_template_s tfx_effect_template_t;
typedef struct tfx_compute_sprite_s tfx_compute_sprite_t;
typedef struct tfx_compute_particle_s tfx_compute_particle_t;
typedef struct tfx_library_s tfx_library_t;
typedef struct tfx_animation_manager_s tfx_animation_manager_t;

//--------------------------------------------------------------
//macros
#define TFX_VERSION "Alpha"
#define TFX_VERSION_NUMBER 6.18.2024

#define tfxMAX_FRAME 20000.f
#define tfxCOLOR_RAMP_WIDTH 256
#define tfxNullParent 0xFFFFFFFF
#define tfxINVALID 0xFFFFFFFF
#define tfxINVALID_SPRITE 0x0FFFFFFF
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

tfxINTERNAL inline tfxU64 tfx__get_hash(tfx_hasher_t *hasher)
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
tfxAPI_EDITOR inline tfxKey tfx_Hash(tfx_hasher_t *hasher, const void *input, tfxU64 length, tfxU64 seed) {
	tfx__hash_initialise(hasher, seed); tfx__hasher_add(hasher, input, length); return (tfxKey)tfx__get_hash(hasher);
}
//-- End of Pocket Hasher

//----------------------------------------------------------
//Section: SIMD_defines
//----------------------------------------------------------

//Currently there's no advantage to using avx so I have some work to do optimising there, probably to do with cache and general memory bandwidth
//Hiding this define here for now until I can test and improve AVX more
//#define tfxUSEAVX

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
#define tfx128SetSingle _mm_set_ps1
#define tfx128SetSinglei _mm_set1_epi32
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
#define tfxWideMuli _mm256_mullo_epi32
#define tfxWideSqrt _mm256_sqrt_ps
#define tfxWideRSqrt _mm256_rsqrt_ps
#define tfxWideMoveMask _mm256_movemask_epi8
#define tfxWideShiftRight _mm256_srli_epi32
#define tfxWideShiftLeft _mm256_slli_epi32
#define tfxWideGreaterEqual(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_GE_OS)
#define tfxWideGreater(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_GT_OS)
#define tfxWideGreateri _mm256_cmpgt_epi32
#define tfxWideLess(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_LT_OS)
#define tfxWideLessi(v1, v2) _mm256_cmpgt_epi32(v2, v1)
#define tfxWideLessEqeual(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_LE_OS)
#define tfxWideEquals(v1, v2) _mm256_cmp_ps(v1, v2, _CMP_EQ_OS)
#define tfxWideEqualsi _mm256_cmpeq_epi32 
#define tfxWideStore _mm256_store_ps
#define tfxWideStorei _mm256_store_si256
#define tfxWideCasti _mm256_castps_si256
#define tfxWideCast _mm256_castsi256_ps 
#define tfxWideConverti _mm256_cvttps_epi32 
#define tfxWideConvert    _mm256_cvtepi32_ps 
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
#define tfxWideLookupSet(lookup, index) tfxWideSet(lookup[index.a[7]], lookup[index.a[6]], lookup[index.a[5]], lookup[index.a[4]], lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )
#define tfxWideLookupSeti(lookup, index) tfxWideSeti(lookup[index.a[7]], lookup[index.a[6]], lookup[index.a[5]], lookup[index.a[4]], lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )
#define tfxWideLookupSetColor(lookup, index) tfxWideSeti(lookup[index.a[7]].color, lookup[index.a[6]].color, lookup[index.a[5]].color, lookup[index.a[4]].color, lookup[index.a[3]].color, lookup[index.a[2]].color, lookup[index.a[1]].color, lookup[index.a[0]].color )
#define tfxWideLookupSetMember(lookup, member, index) tfxWideSet(lookup[index.a[7]].member, lookup[index.a[6]].member, lookup[index.a[5]].member, lookup[index.a[4]].member, lookup[index.a[3]].member, lookup[index.a[2]].member, lookup[index.a[1]].member, lookup[index.a[0]].member )
#define tfxWideLookupSetMemberi(lookup, member, index) tfxWideSeti(lookup[index.a[7]].member, lookup[index.a[6]].member, lookup[index.a[5]].member, lookup[index.a[4]].member, lookup[index.a[3]].member, lookup[index.a[2]].member, lookup[index.a[1]].member, lookup[index.a[0]].member )
#define tfxWideLookupSet2(lookup1, lookup2, index1, index2) tfxWideSet(lookup1[index1.a[7]].lookup2[index2.a[7]], lookup1[index1.a[6]].lookup2[index2.a[6]], lookup1[index1.a[5]].lookup2[index2.a[5]], lookup1[index1.a[4]].lookup2[index2.a[4]], lookup1[index1.a[3]].lookup2[index2.a[3]], lookup1[index1.a[2]].lookup2[index2.a[2]], lookup1[index1.a[1]].lookup2[index2.a[1]], lookup1[index1.a[0]].lookup2[index2.a[0]] )
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

const __m256 tfxWIDEF3_4 = tfxWideSetConst(1.0f / 3.0f);
const __m256 tfxWIDEG3_4 = tfxWideSetConst(1.0f / 6.0f);
const __m256 tfxWIDEG32_4 = tfxWideSetConst((1.0f / 6.0f) * 2.f);
const __m256 tfxWIDEG33_4 = tfxWideSetConst((1.0f / 6.0f) * 3.f);
const __m256i tfxWIDEONEi = tfxWideSetConst(1);
const __m256 tfxWIDEMINUSONE = tfxWideSetConst(-1.f);
const __m256i tfxWIDEMINUSONEi = tfxWideSetConst(-1);
const __m256 tfxWIDEONE = tfxWideSetConst(1.f);
const __m256 tfxWIDE255 = tfxWideSetConst(255.f);
const __m256 tfxWIDEZERO = tfxWideSetConst(0.f);
const __m256 tfxWIDETHIRTYTWO = tfxWideSetConst(32.f);
const __m256 tfxPWIDESIX = tfxWideSetConst(0.6f);
const __m256 tfxMAXUINTf = tfxWideSetConst((float)UINT32_MAX);
const __m256 tfxDEGREERANGEMR = tfxWideSetConst(0.392699f);
tfxINTERNAL const __m256 SIGNMASK = tfxWideSetConst(-0.f);

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
#define tfxWideAdd _mm_add_ps
#define tfxWideSub _mm_sub_ps
#define tfxWideMul _mm_mul_ps
#define tfxWideDiv _mm_div_ps
#define tfxWideAddi _mm_add_epi32
#define tfxWideSubi _mm_sub_epi32
#define tfxWideMuli _mm_mullo_epi32
#define tfxWideSqrt _mm_sqrt_ps
#define tfxWideRSqrt _mm_rsqrt_ps
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

#define tfxWideSetConst(value) {value, value, value, value}

typedef union {
	int a[4];
	__m128i m;
} tfxWideArrayi;

typedef union {
	float a[4];
	__m128 m;
} tfxWideArray;

const tfxWideArray tfxWIDEF3_4 = tfxWideSetConst(1.f / 3.f);
const tfxWideArray tfxWIDEG3_4 = tfxWideSetConst(1.0f / 6.0f);
const tfxWideArray tfxWIDEG32_4 = tfxWideSetConst((1.0f / 6.0f) * 2.f);
const tfxWideArray tfxWIDEG33_4 = tfxWideSetConst((1.0f / 6.0f) * 3.f);
const tfxWideArrayi tfxWIDEONEi = tfxWideSetConst(1);
const tfxWideArray tfxWIDEMINUSONE = tfxWideSetConst(-1.f);
const tfxWideArrayi tfxWIDEMINUSONEi = tfxWideSetConst(-1);
const tfxWideArray tfxWIDEONE = tfxWideSetConst(1.f);
const tfxWideArray tfxWIDE255 = tfxWideSetConst(255.f);
const tfxWideArray tfxWIDEZERO = tfxWideSetConst(0.f);
const tfxWideArray tfxWIDETHIRTYTWO = tfxWideSetConst(32.f);
const tfxWideArray tfxPWIDESIX = tfxWideSetConst(0.6f);
const tfxWideArray tfxMAXUINTf = tfxWideSetConst((float)UINT32_MAX);
const tfxWideArray tfxDEGREERANGEMR = tfxWideSetConst(0.392699f);
const tfxWideArray SIGNMASK = tfxWideSetConst(-0.f);

#elif defined(tfxARM)
//Arm Intrinsics
typedef float32x4_t tfxWideFloat;
typedef int32x4_t tfxWideInt;
typedef int32_t tfxWideIntLoader;
#define tfxWideLoad vld1q_f32
#define tfxWideLoadi vld1q_s32
inline __attribute__((always_inline)) float32x4_t tfx__128_SET(float e3, float e2, float e1, float e0) {
	float32x4_t r;
	alignas(16) float data[4] = { e0, e1, e2, e3 };
	r = vld1q_f32(data);
	return r;
}
inline __attribute__((always_inline)) int32x4_t tfx__128i_SET(int e3, int e2, int e1, int e0) {
	int32x4_t r;
	alignas(16) int data[4] = { e0, e1, e2, e3 };
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
#define tfxWideAdd vaddq_f32
#define tfxWideSub vsubq_f32
#define tfxWideMul vmulq_f32
#define tfxWideDiv vdivq_f32
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

const float32x4_t tfxWIDEF3_4 = tfxWideSetConst(1.0f / 3.0f);
const float32x4_t tfxWIDEG3_4 = tfxWideSetConst(1.0f / 6.0f);
const float32x4_t tfxWIDEG32_4 = tfxWideSetConst((1.0f / 6.0f) * 2.f);
const float32x4_t tfxWIDEG33_4 = tfxWideSetConst((1.0f / 6.0f) * 3.f);
const int32x4_t tfxWIDEONEi = tfxWideSetConst(1);
const float32x4_t tfxWIDEMINUSONE = tfxWideSetConst(-1.f);
const int32x4_t tfxWIDEMINUSONEi = tfxWideSetConst(-1);
const float32x4_t tfxWIDEONE = tfxWideSetConst(1.f);
const float32x4_t tfxWIDE255 = tfxWideSetConst(255.f);
const float32x4_t tfxWIDEZERO = tfxWideSetConst(0.f);
const float32x4_t tfxWIDETHIRTYTWO = tfxWideSetConst(32.f);
const float32x4_t tfxPWIDESIX = tfxWideSetConst(0.6f);
const float32x4_t tfxMAXUINTf = tfxWideSetConst((float)UINT32_MAX);
const float32x4_t tfxDEGREERANGEMR = tfxWideSetConst(0.392699f);

tfxINTERNAL const float32x4_t SIGNMASK = tfxWideSetConst(-0.f);

#endif

#define tfxWideLookupSet(lookup, index) tfx128Set( lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )
#define tfxWideLookupSetMember(lookup, member, index) tfx128Set( lookup[index.a[3]].member, lookup[index.a[2]].member, lookup[index.a[1]].member, lookup[index.a[0]].member )
#define tfxWideLookupSetMemberi(lookup, member, index) tfx128Seti( lookup[index.a[3]].member, lookup[index.a[2]].member, lookup[index.a[1]].member, lookup[index.a[0]].member )
#define tfxWideLookupSet2(lookup1, lookup2, index1, index2) tfx128Set( lookup1[index1.a[3]].lookup2[index2.a[3]], lookup1[index1.a[2]].lookup2[index2.a[2]], lookup1[index1.a[1]].lookup2[index2.a[1]], lookup1[index1.a[0]].lookup2[index2.a[0]] )
#define tfxWideLookupSeti(lookup, index) tfx128Seti( lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )
#define tfxWideLookupSetColor(lookup, index) tfx128Seti( lookup[index.a[3]].color, lookup[index.a[2]].color, lookup[index.a[1]].color, lookup[index.a[0]].color )
#define tfxWideLookupSetOffset(lookup, index, offset) tfx128Set( lookup[index.a[3] + offset], lookup[index.a[2] + offset], lookup[index.a[1] + offset], lookup[index.a[0] + offset] )

#endif

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

#define tfx128SetConst(value) {value, value, value, value}

//simd floor function thanks to Stephanie Rancourt: http://dss.stephanierct.com/DevBlog/?p=8
tfxINTERNAL inline tfx128 tfxFloor128(const tfx128 x) {
	//__m128i v0 = _mm_setzero_si128();
	//__m128i v1 = _mm_cmpeq_epi32(v0, v0);
	//__m128i ji = _mm_srli_epi32(v1, 25);
	//__m128 j = *(__m128*)&_mm_slli_epi32(ji, 23); //create vector 1.0f
	//Worth noting that we only need to floor small numbers for the noise algorithm so can get away with this function.
	__m128 j = _mm_set1_ps(1.f); //create vector 1.0f
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
	return tfxWideSub(tfxWideMul(tfxWideSetSingle(tfxQUARTERPI), x),
		tfxWideMul(tfxWideMul(x, tfxWideSub(tfxWideAndNot(SIGNMASK.m, x), tfxWideSetSingle(1.f))),
			(tfxWideAdd(tfxWideSetSingle(0.2447f), tfxWideMul(tfxWideSetSingle(0.0663f), tfxWideAndNot(SIGNMASK.m, x))))));
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
	tfxWideFloat absxgreaterthanabsy = tfxWideGreater(tfxWideAndNot(SIGNMASK.m, x), tfxWideAndNot(SIGNMASK.m, y));
	tfxWideFloat ratio = tfxWideDiv(tfxWideAdd(tfxWideAnd(absxgreaterthanabsy, y), tfxWideAndNot(absxgreaterthanabsy, x)),
		tfxWideAdd(tfxWideAnd(absxgreaterthanabsy, x), tfxWideAndNot(absxgreaterthanabsy, y)));
	tfxWideFloat atan = tfxWideAtan(ratio);

	tfxWideFloat xgreaterthan0 = tfxWideGreater(x, tfxWideSetSingle(0.f));
	tfxWideFloat ygreaterthan0 = tfxWideGreater(y, tfxWideSetSingle(0.f));

	atan = tfxWideXOr(atan, tfxWideAndNot(absxgreaterthanabsy, SIGNMASK.m)); //negate atan if absx<=absy & x>0

	tfxWideFloat shift = tfxWideSetSingle(tfxPI);
	shift = tfxWideSub(shift, tfxWideAndNot(absxgreaterthanabsy, tfxWideSetSingle(tfxHALFPI))); //substract tfxHALFPI if absx<=absy
	shift = tfxWideXOr(shift, tfxWideAndNot(ygreaterthan0, SIGNMASK.m)); //negate shift if y<=0
	shift = tfxWideAndNot(tfxWideAnd(absxgreaterthanabsy, xgreaterthan0), shift); //null if abs>absy & x>0

	return tfxWideAdd(atan, shift);
}

inline tfxWideFloat tfxWideCos52s(tfxWideFloat x)
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

inline void tfxWideSinCos(tfxWideFloat angle, tfxWideFloat *sin, tfxWideFloat *cos) {
	tfxWideFloat anglesign = tfxWideOr(tfxWideSetSingle(1.f), tfxWideAnd(SIGNMASK.m, angle));

	//clamp to the range 0..2pi

	//take absolute value
	angle = tfxWideAndNot(SIGNMASK.m, angle);
	//fmod(angle,tfxPI2)
	angle = tfxWideSub(angle, tfxWideMul(tfxWideConvert(tfxWideConverti(tfxWideMul(angle, tfxWideSetSingle(tfxINVTWOPI)))), tfxWideSetSingle(tfxPI2))); //simplied SSE2 fmod, must always operate on absolute value
	//if SSE4.1 is always available, comment the line above and uncomment the line below
	//angle=tfxWideSub(angle,tfxWideMul(_mm_floor_ps(tfxWideMul(angle,tfxWideSetSingle(tfxINVTWOPI))),tfxWideSetSingle(tfxPI2))); //faster if SSE4.1 is always available

	tfxWideFloat cosangle = angle;
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWideSetSingle(tfxHALFPI)), tfxWideXOr(cosangle, tfxWideSub(tfxWideSetSingle(tfxPI), angle))));
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWideSetSingle(tfxPI)), SIGNMASK.m));
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWideSetSingle(tfxTHREEHALFPI)), tfxWideXOr(cosangle, tfxWideSub(tfxWideSetSingle(tfxPI2), angle))));

	tfxWideFloat result = tfxWideCos52s(cosangle);

	result = tfxWideXOr(result, tfxWideAnd(tfxWideAnd(tfxWideGreaterEqual(angle, tfxWideSetSingle(tfxHALFPI)), tfxWideLess(angle, tfxWideSetSingle(tfxTHREEHALFPI))), SIGNMASK.m));
	*cos = result;

	tfxWideFloat sinmultiplier = tfxWideMul(anglesign, tfxWideOr(tfxWideSetSingle(1.f), tfxWideAnd(tfxWideGreater(angle, tfxWideSetSingle(tfxPI)), SIGNMASK.m)));
	*sin = tfxWideMul(sinmultiplier, tfxWideFastSqrt(tfxWideSub(tfxWideSetSingle(1.f), tfxWideMul(result, result))));
}

inline void tfxWideSinCosAdd(tfxWideFloat angle, tfxWideFloat *sin, tfxWideFloat *cos) {
	tfxWideFloat anglesign = tfxWideOr(tfxWideSetSingle(1.f), tfxWideAnd(SIGNMASK.m, angle));

	//clamp to the range 0..2pi

	//take absolute value
	angle = tfxWideAndNot(SIGNMASK.m, angle);
	//fmod(angle,tfxPI2)
	angle = tfxWideSub(angle, tfxWideMul(tfxWideConvert(tfxWideConverti(tfxWideMul(angle, tfxWideSetSingle(tfxINVTWOPI)))), tfxWideSetSingle(tfxPI2))); //simplied SSE2 fmod, must always operate on absolute value
	//if SSE4.1 is always available, comment the line above and uncomment the line below
	//angle=tfxWideSub(angle,tfxWideMul(_mm_floor_ps(tfxWideMul(angle,tfxWideSetSingle(tfxINVTWOPI))),tfxWideSetSingle(tfxPI2))); //faster if SSE4.1 is always available

	tfxWideFloat cosangle = angle;
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWideSetSingle(tfxHALFPI)), tfxWideXOr(cosangle, tfxWideSub(tfxWideSetSingle(tfxPI), angle))));
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWideSetSingle(tfxPI)), SIGNMASK.m));
	cosangle = tfxWideXOr(cosangle, tfxWideAnd(tfxWideGreaterEqual(angle, tfxWideSetSingle(tfxTHREEHALFPI)), tfxWideXOr(cosangle, tfxWideSub(tfxWideSetSingle(tfxPI2), angle))));

	tfxWideFloat result = tfxWideCos52s(cosangle);

	result = tfxWideXOr(result, tfxWideAnd(tfxWideAnd(tfxWideGreaterEqual(angle, tfxWideSetSingle(tfxHALFPI)), tfxWideLess(angle, tfxWideSetSingle(tfxTHREEHALFPI))), SIGNMASK.m));
	*cos = tfxWideAdd(*cos, result);

	tfxWideFloat sinmultiplier = tfxWideMul(anglesign, tfxWideOr(tfxWideSetSingle(1.f), tfxWideAnd(tfxWideGreater(angle, tfxWideSetSingle(tfxPI)), SIGNMASK.m)));
	*sin = tfxWideAdd(*sin, tfxWideMul(sinmultiplier, tfxWideFastSqrt(tfxWideSub(tfxWideSetSingle(1.f), tfxWideMul(result, result)))));
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
	return tfxWideAnd(tfxWideCast(tfxWideShiftRight(tfxWideSetSinglei(-1), 1)), v);
}

tfxINTERNAL inline tfxWideInt tfxWideAbsi(tfxWideInt v) {
	return tfxWideAndi(tfxWideShiftRight(tfxWideSetSinglei(-1), 1), v);
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
	tfxIntensityOvertimePreset,
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

#define TFX_GLOBAL_COUNT     17
#define TFX_PROPERTY_COUNT   10
#define TFX_BASE_COUNT       10
#define TFX_VARIATION_COUNT  12
#define TFX_OVERTIME_COUNT   25
#define TFX_FACTOR_COUNT     4
#define TFX_TRANSFORM_COUNT  6

enum {
	TFX_GLOBAL_START = 0,
	TFX_PROPERTY_START = TFX_GLOBAL_COUNT,
	TFX_BASE_START = (TFX_PROPERTY_START + TFX_PROPERTY_COUNT),
	TFX_VARIATION_START = (TFX_BASE_START + TFX_BASE_COUNT),
	TFX_OVERTIME_START = (TFX_VARIATION_START + TFX_VARIATION_COUNT),
	TFX_FACTOR_START = (TFX_OVERTIME_START + TFX_OVERTIME_COUNT),
	TFX_TRANSFORM_START = (TFX_FACTOR_START + TFX_FACTOR_COUNT)
};

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
	tfxVariation_pitch_spin,
	tfxVariation_yaw_spin,
	tfxVariation_roll_spin,
	tfxVariation_noise_offset,
	tfxVariation_noise_resolution,
	tfxVariation_motion_randomness,

	tfxOvertime_velocity,
	tfxOvertime_width,
	tfxOvertime_height,
	tfxOvertime_weight,
	tfxOvertime_pitch_spin,
	tfxOvertime_yaw_spin,
	tfxOvertime_roll_spin,
	tfxOvertime_stretch,
	tfxOvertime_red,
	tfxOvertime_green,
	tfxOvertime_blue,
	tfxOvertime_blendfactor,
	tfxOvertime_red_hint,
	tfxOvertime_green_hint,
	tfxOvertime_blue_hint,
	tfxOvertime_blendfactor_hint,
	tfxOvertime_velocity_turbulance,
	tfxOvertime_direction_turbulance,
	tfxOvertime_velocity_adjuster,
	tfxOvertime_intensity,
	tfxOvertime_alpha_sharpness,
	tfxOvertime_curved_alpha,
	tfxOvertime_direction,
	tfxOvertime_noise_resolution,
	tfxOvertime_motion_randomness,

	tfxFactor_life,
	tfxFactor_size,
	tfxFactor_velocity,
	tfxFactor_intensity,

	tfxTransform_roll,
	tfxTransform_pitch,
	tfxTransform_yaw,
	tfxTransform_translate_x,
	tfxTransform_translate_y,
	tfxTransform_translate_z,
	tfxEmitterGraphMaxIndex,

	tfxPath_angle_x,
	tfxPath_angle_y,
	tfxPath_angle_z,
	tfxPath_offset_x,
	tfxPath_offset_y,
	tfxPath_offset_z,
	tfxPath_distance,
	tfxPath_rotation_range,
	tfxPath_rotation_pitch,
	tfxPath_rotation_yaw,
	tfxGraphMaxIndex
} tfx_graph_type;

//tfx_effect_emitter_t type - effect contains emitters, and emitters spawn particles, but they both share the same struct for simplicity
typedef enum {
	tfxEffectType,
	tfxEmitterType,
	tfxStage,
	tfxFolder
} tfx_effect_emitter_type;

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
	tfxEmissionTypeMax,
} tfx_emission_type;

typedef enum {
	tfxExtrusionArc,
	tfxExtrusionLinear
} tfx_path_extrusion_type;

//Determines how for area, line and ellipse emitters the direction that particles should travel when they spawn
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
	tfxCompressFrames,
	tfxBakingDone
} tfx_record_progress;

//Used in file loading - for loading effects library
typedef enum {
	tfxString,
	tfxSInt,
	tfxUint,
	tfxFloat,
	tfxDouble,
	tfxBool,
	tfxColor,
	tfxUInt64,
	tfxFloat3,
	tfxFloat2
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
} tfx_effect_library_stream_context;

typedef tfxU32 tfxEmitterPropertyFlags;         //tfx_emitter_property_flag_bits
typedef tfxU32 tfxColorRampFlags;		        //tfx_color_ramp_flag_bits
typedef tfxU32 tfxEffectPropertyFlags;          //tfx_effect_property_flag_bits
typedef tfxU32 tfxVectorFieldFlags;             //tfx_vector_field_flag_bits
typedef tfxU32 tfxParticleFlags;                //tfx_particle_flag_bits
typedef tfxU32 tfxEmitterStateFlags;            //tfx_emitter_state_flag_bits
typedef tfxU32 tfxEffectStateFlags;             //tfx_effect_state_flag_bits
typedef tfxU32 tfxParticleControlFlags;         //tfx_particle_control_flag_bits
typedef tfxU32 tfxAttributeNodeFlags;           //tfx_attribute_node_flag_bits
typedef tfxU32 tfxAngleSettingFlags;            //tfx_angle_setting_flag_bits
typedef tfxU32 tfxParticleManagerFlags;         //tfx_particle_manager_flag_bits
typedef tfxU32 tfxErrorFlags;                   //tfx_error_flag_bits
typedef tfxU32 tfxEffectCloningFlags;           //tfx_effect_cloning_flag_bits
typedef tfxU32 tfxAnimationFlags;               //tfx_animation_flag_bits
typedef tfxU32 tfxAnimationInstanceFlags;       //tfx_animation_instance_flag_bits
typedef tfxU32 tfxAnimationManagerFlags;        //tfx_animation_manager_flag_bits
typedef tfxU32 tfxEmitterPathFlags;             //tfx_emitter_path_flag_bits
typedef tfxU32 tfxEmitterControlProfileFlags;   //tfx_emitter_control_profile_flag_bits
typedef tfxU32 tfxPackageFlags;                 //tfx_package_flag_bits

typedef enum {
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
} tfx_error_flag_bits;

typedef enum {
	tfxPackageFlags_none = 0,
	tfxPackageFlags_loaded_from_memory = 1
} tfx_package_flag_bits;

typedef enum {
	tfxEffectCloningFlags_none = 0,
	tfxEffectCloningFlags_keep_user_data = 1 << 0,
	tfxEffectCloningFlags_force_clone_global = 1 << 1,
	tfxEffectCloningFlags_clone_graphs = 1 << 2,
	tfxEffectCloningFlags_compile_graphs = 1 << 3
} tfx_effect_cloning_flag_bits;

typedef enum {
	tfxParticleManagerMode_unordered,
	tfxParticleManagerMode_ordered_by_age,
	tfxParticleManagerMode_ordered_by_depth,
	tfxParticleManagerMode_ordered_by_depth_guaranteed
} tfx_particle_manager_mode;

typedef enum {
	tfxParticleManagerSetup_none,
	tfxParticleManagerSetup_2d_unordered,
	tfxParticleManagerSetup_2d_ordered_by_age,
	tfxParticleManagerSetup_2d_group_sprites_by_effect,
	tfxParticleManagerSetup_3d_unordered,
	tfxParticleManagerSetup_3d_ordered_by_age,
	tfxParticleManagerSetup_3d_ordered_by_depth,
	tfxParticleManagerSetup_3d_ordered_by_depth_guaranteed,
	tfxParticleManagerSetup_3d_group_sprites_by_effect,
} tfx_particle_manager_setup;

typedef enum {
	tfxBillboarding_align_to_camera = 0,            //Align to Camera only
	tfxBillboarding_free_align = 1,                    //Free align
	tfxBillboarding_align_to_camera_and_vector = 2,    //Align to camera and vector
	tfxBillboarding_align_to_vector = 3,            //Align to vector
	tfxBillboarding_max = 4
} tfx_billboarding_option;

typedef enum {
	tfxEmitterControlProfile_basic = 0,
	tfxEmitterControlProfile_noise = 1 << 0,
	tfxEmitterControlProfile_orbital = 1 << 1,
	tfxEmitterControlProfile_motion_randomness = 1 << 2,
	tfxEmitterControlProfile_path = 1 << 3,
	tfxEmitterControlProfile_path_rotated_path = 1 << 4,
	tfxEmitterControlProfile_edge_traversal = 1 << 5,
	tfxEmitterControlProfile_edge_kill = 1 << 6,
	tfxEmitterControlProfile_edge_loop = 1 << 7,
	tfxEmitterControlProfile_stretch = 1 << 8
} tfx_emitter_control_profile_flag_bits;

typedef enum {
	tfxParticleManagerFlags_none = 0,
	tfxParticleManagerFlags_disable_spawning = 1,
	tfxParticleManagerFlags_force_capture = 2,            //Unused
	tfxParticleManagerFlags_use_compute_shader = 1 << 3,
	tfxParticleManagerFlags_order_by_depth = 1 << 4,
	tfxParticleManagerFlags_guarantee_order = 1 << 5,
	tfxParticleManagerFlags_update_base_values = 1 << 6,
	tfxParticleManagerFlags_dynamic_sprite_allocation = 1 << 7,
	tfxParticleManagerFlags_3d_effects = 1 << 8,
	tfxParticleManagerFlags_unordered = 1 << 9,
	tfxParticleManagerFlags_ordered_by_age = 1 << 10,
	tfxParticleManagerFlags_update_age_only = 1 << 11,
	tfxParticleManagerFlags_single_threaded = 1 << 12,
	tfxParticleManagerFlags_double_buffer_sprites = 1 << 13,
	tfxParticleManagerFlags_recording_sprites = 1 << 14,
	tfxParticleManagerFlags_using_uids = 1 << 15,
	tfxParticleManagerFlags_2d_and_3d = 1 << 16,
	tfxParticleManagerFlags_update_bounding_boxes = 1 << 17,
	tfxParticleManagerFlags_use_effect_sprite_buffers = 1 << 18,
	tfxParticleManagerFlags_auto_order_effects = 1 << 19,
	tfxParticleManagerFlags_direct_to_staging_buffer = 1 << 20,
} tfx_particle_manager_flag_bits;

typedef enum {
	tfxVectorAlignType_motion,
	tfxVectorAlignType_emission,
	tfxVectorAlignType_emitter,
	tfxVectorAlignType_max,
} tfx_vector_align_type;

typedef enum {
	tfxPathGenerator_spiral,
	tfxPathGenerator_loop,
	tfxPathGenerator_arc,
	tfxPathGenerator_s_curve,
	tfxPathGenerator_bend,
	tfxPathGenerator_free_mode_origin,
	tfxPathGenerator_free_mode_distance,
	tfxPathGenerator_max,
} tfx_path_generator_type;

typedef enum {
	tfxPathFlags_none,
	tfxPathFlags_2d = 1 << 0,
	tfxPathFlags_mode_origin = 1 << 1,
	tfxPathFlags_mode_node = 1 << 2,
	tfxPathFlags_space_nodes_evenly = 1 << 3,
	tfxPathFlags_reverse_direction = 1 << 4,
	tfxPathFlags_rotation_range_yaw_only = 1 << 5
} tfx_emitter_path_flag_bits;

//Particle property that defines how a particle will rotate
typedef enum {
	tfxAngleSettingFlags_none = 0,                                                        //No flag
	tfxAngleSettingFlags_align_roll = 1 << 0,                                            //Align the particle with it's direction of travel in 2d
	tfxAngleSettingFlags_random_roll = 1 << 1,                                            //Chose a random angle at spawn time/state_flags
	tfxAngleSettingFlags_specify_roll = 1 << 2,                                            //Specify the angle at spawn time
	tfxAngleSettingFlags_align_with_emission = 1 << 3,                                    //Align the particle with the emission direction only
	tfxAngleSettingFlags_random_pitch = 1 << 4,                                            //3d mode allows for rotating pitch and yaw when not using billboarding (when particle always faces the camera)
	tfxAngleSettingFlags_random_yaw = 1 << 5,
	tfxAngleSettingFlags_specify_pitch = 1 << 6,
	tfxAngleSettingFlags_specify_yaw = 1 << 7
} tfx_angle_setting_flag_bits;

//All the state_flags needed by the ControlParticle function put into one typedef enum save typedef enum
typedef enum {
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
} tfx_particle_control_flag_bits;

typedef enum {
	tfxEffectPropertyFlags_none = 0,
	tfxEffectPropertyFlags_depth_draw_order = 1 << 1,
	tfxEffectPropertyFlags_guaranteed_order = 1 << 2,
	tfxEffectPropertyFlags_age_order = 1 << 3,
	tfxEffectPropertyFlags_use_keyframes = 1 << 4,
	tfxEffectPropertyFlags_include_in_sprite_data_export = 1 << 5,        //In the editor you can specify which effects you want to be included in a spritedata export
	tfxEffectPropertyFlags_global_uniform_size = 1 << 6,                //Keep the global particle size uniform
} tfx_effect_property_flag_bits;

typedef enum {
	tfxEmitterPropertyFlags_none = 0,
	tfxEmitterPropertyFlags_random_color = 1 << 0,                      //Pick a random color from the color overtime gradient rather then change the color over the lifetime of the particle
	tfxEmitterPropertyFlags_relative_position = 1 << 1,                 //Keep the particles position relative to the current position of the emitter
	tfxEmitterPropertyFlags_relative_angle = 1 << 2,                    //Keep the angle of the particles relative to the current angle of the emitter
	tfxEmitterPropertyFlags_image_handle_auto_center = 1 << 3,          //Set the offset of the particle to the center of the image
	tfxEmitterPropertyFlags_single = 1 << 4,                            //Only spawn a single particle (or number of particles specified by spawn_amount) that does not expire
	tfxEmitterPropertyFlags_spawn_on_grid = 1 << 5,                     //When using an area, line or ellipse emitter, spawn along a grid
	tfxEmitterPropertyFlags_grid_spawn_clockwise = 1 << 6,              //Spawn clockwise/left to right around the area
	tfxEmitterPropertyFlags_fill_area = 1 << 7,                         //Fill the area
	tfxEmitterPropertyFlags_emitter_handle_auto_center = 1 << 8,        //Center the handle of the emitter
	tfxEmitterPropertyFlags_edge_traversal = 1 << 9,                    //Line and Path emitters only: make particles traverse the line/path
	tfxEmitterPropertyFlags_base_uniform_size = 1 << 10,                //Keep the base particle size uniform
	tfxEmitterPropertyFlags_lifetime_uniform_size = 1 << 11,            //Keep the size over lifetime of the particle uniform
	tfxEmitterPropertyFlags_animate = 1 << 12,                          //Animate the particle shape if it has more than one frame of animation
	tfxEmitterPropertyFlags_reverse_animation = 1 << 13,                //Make the image animation go in reverse
	tfxEmitterPropertyFlags_play_once = 1 << 14,                        //Play the animation once only
	tfxEmitterPropertyFlags_random_start_frame = 1 << 15,               //Start the animation of the image from a random frame
	tfxEmitterPropertyFlags_wrap_single_sprite = 1 << 16,               //When recording sprite data, single particles can have their invalid capured index set to the current frame for better looping
	tfxEmitterPropertyFlags_is_in_folder = 1 << 17,                     //This effect is located inside a folder. Todo: move this to effect properties
	tfxEmitterPropertyFlags_use_spawn_ratio = 1 << 18,                  //Option for area emitters to multiply the amount spawned by a ration of particles per pixels squared
	tfxEmitterPropertyFlags_effect_is_3d = 1 << 19,                     //Makes the effect run in 3d mode for 3d effects todo: does this need to be here, the effect dictates this?
	tfxEmitterPropertyFlags_grid_spawn_random = 1 << 20,                //Spawn on grid points but randomly rather then in sequence
	tfxEmitterPropertyFlags_area_open_ends = 1 << 21,                   //Only sides of the area/cylinder are spawned on when fill area is not checked
	tfxEmitterPropertyFlags_exclude_from_hue_adjustments = 1 << 22,     //Emitter will be excluded from effect hue adjustments if this flag is checked
	tfxEmitterPropertyFlags_enabled = 1 << 23,                          //The emitter is enabled or not, meaning it will or will not be added the particle manager with tfx__add_effect
	tfxEmitterPropertyFlags_match_amount_to_grid_points = 1 << 24,      //Match the amount to spawn with a single emitter to the number of grid points in the effect
	tfxEmitterPropertyFlags_use_path_for_direction = 1 << 25,           //Make the particles use a path to dictate their direction of travel
	tfxEmitterPropertyFlags_alt_velocity_lifetime_sampling = 1 << 26,	//The point on the path dictates where on the velocity overtime graph that the particle should sample from rather then the age of the particle
	tfxEmitterPropertyFlags_alt_color_lifetime_sampling = 1 << 27,      //The point on the path dictates where on the color overtime graph that the particle should sample from rather then the age of the particle
	tfxEmitterPropertyFlags_alt_size_lifetime_sampling = 1 << 28,       //The point on the path dictates where on the size overtime graph that the particle should sample from rather then the age of the particle
	tfxEmitterPropertyFlags_use_simple_motion_randomness = 1 << 29,     //Use a simplified way to generate random particle movement which is much less computationally intensive than simplex noise
	tfxEmitterPropertyFlags_spawn_location_source = 1 << 30,            //This emitter is the source for another emitter that uses it to spawn particles at the location of this emitters' particles
	tfxEmitterPropertyFlags_use_color_hint = 1 << 31	                //Activate a second color to tint the particles and mix between the two colors.
} tfx_emitter_property_flag_bits;

typedef enum {
	tfxColorRampFlags_none = 0,
	tfxColorRampFlags_use_sinusoidal_ramp_generation = 1 << 0,			//Use this flag to toggle between sinusoidal color ramp generation
} tfx_color_ramp_flag_bits;

typedef enum {
	tfxParticleFlags_none = 0,
	tfxParticleFlags_fresh = 1 << 0,                                    //Particle has just spawned this frame    
	tfxParticleFlags_remove = 1 << 4,                                   //Particle will be removed this or next frame
	tfxParticleFlags_has_velocity = 1 << 5,                             //Flagged if the particle is currently moving
	tfxParticleFlags_has_sub_effects = 1 << 6,                          //Flagged if the particle has sub effects
	tfxParticleFlags_capture_after_transform = 1 << 15,                 //Particle will be captured after a transfrom, used for traversing lines and looping back to the beginning to avoid lerping imbetween
} tfx_particle_flag_bits;

typedef enum {
	tfxEmitterStateFlags_none = 0,
	tfxEmitterStateFlags_random_color = 1 << 0,
	tfxEmitterStateFlags_relative_position = 1 << 1,                    //Keep the particles position relative to the current position of the emitter
	tfxEmitterStateFlags_relative_angle = 1 << 2,                       //Keep the angle of the particles relative to the current angle of the emitter
	tfxEmitterStateFlags_stop_spawning = 1 << 3,                        //Tells the emitter to stop spawning
	tfxEmitterStateFlags_remove = 1 << 4,                               //Tells the effect/emitter to remove itself from the particle manager immediately
	tfxEmitterStateFlags_unused1 = 1 << 5,                              //the emitter is enabled. **moved to property state_flags**
	tfxEmitterStateFlags_retain_matrix = 1 << 6,                        //Internal flag about matrix usage
	tfxEmitterStateFlags_no_tween_this_update = 1 << 7,                 //Internal flag generally, but you could use it if you want to teleport the effect to another location
	tfxEmitterStateFlags_is_single = 1 << 8,
	tfxEmitterStateFlags_not_line = 1 << 9,
	tfxEmitterStateFlags_base_uniform_size = 1 << 10,
	tfxEmitterStateFlags_lifetime_uniform_size = 1 << 11,               //Keep the size over lifetime of the particle uniform
	tfxEmitterStateFlags_can_spin = 1 << 12,
	tfxEmitterStateFlags_is_edge_traversal = 1 << 13,
	tfxEmitterStateFlags_play_once = 1 << 14,                           //Play the animation once only
	tfxEmitterStateFlags_loop = 1 << 15,
	tfxEmitterStateFlags_kill = 1 << 16,
	tfxEmitterStateFlags_single_shot_done = 1 << 17,
	tfxEmitterStateFlags_is_line_loop_or_kill = 1 << 18,
	tfxEmitterStateFlags_is_area = 1 << 19,
	tfxEmitterStateFlags_no_tween = 1 << 20,
	tfxEmitterStateFlags_align_with_velocity = 1 << 21,
	tfxEmitterStateFlags_is_sub_emitter = 1 << 22,
	tfxEmitterStateFlags_unused2 = 1 << 23,
	tfxEmitterStateFlags_can_spin_pitch_and_yaw = 1 << 24,              //For 3d emitters that have free alignment and not always facing the camera
	tfxEmitterStateFlags_has_path = 1 << 25,
	tfxEmitterStateFlags_is_bottom_emitter = 1 << 26,                   //This emitter has no child effects, so can spawn particles that could be used in a compute shader if it's enabled
	tfxEmitterStateFlags_has_rotated_path = 1 << 27,
	tfxEmitterStateFlags_max_active_paths_reached = 1 << 28,
	tfxEmitterStateFlags_is_in_ordered_effect = 1 << 29,
	tfxEmitterStateFlags_wrap_single_sprite = 1 << 30
} tfx_emitter_state_flag_bits;

typedef enum {
	tfxEffectStateFlags_none = 0,
	tfxEffectStateFlags_stop_spawning = 1 << 3,                            //Tells the emitter to stop spawning
	tfxEffectStateFlags_remove = 1 << 4,                                //Tells the effect/emitter to remove itself from the particle manager immediately
	tfxEffectStateFlags_retain_matrix = 1 << 6,                            //Internal flag about matrix usage
	tfxEffectStateFlags_no_tween_this_update = 1 << 7,                    //Internal flag generally, but you could use it if you want to teleport the effect to another location
	tfxEffectStateFlags_override_overal_scale = 1 << 8,                    //Flagged when the over scale is overridden with tfx_SetEffectOveralScale
	tfxEffectStateFlags_override_orientiation = 1 << 9,                    //Flagged when any of the effect angles are overridden
	tfxEffectStateFlags_override_size_multiplier = 1 << 10,                //Flagged when any of the effect size multipliers are overridden
	tfxEffectStateFlags_no_tween = 1 << 20
} tfx_effect_state_flag_bits;

typedef enum {
	tfxVectorFieldFlags_none = 0,
	tfxVectorFieldFlags_repeat_horizontal = 1 << 0,                        //Field will repeat horizontally
	tfxVectorFieldFlags_repeat_vertical = 1 << 1                        //Field will repeat vertically
} tfx_vector_field_flag_bits;

typedef enum {
	tfxAttributeNodeFlags_none = 0,
	tfxAttributeNodeFlags_is_curve = 1 << 0,
	tfxAttributeNodeFlags_is_left_curve = 1 << 1,
	tfxAttributeNodeFlags_is_right_curve = 1 << 2,
	tfxAttributeNodeFlags_curves_initialised = 1 << 3,
	tfxAttributeNodeFlags_path_node_accumulate = 1 << 4,
} tfx_attribute_node_flag_bits;

typedef enum {
	tfxAnimationFlags_none = 0,
	tfxAnimationFlags_loop = 1 << 0,
	tfxAnimationFlags_seamless = 1 << 1,
	tfxAnimationFlags_needs_recording = 1 << 2,
	tfxAnimationFlags_export_with_transparency = 1 << 3,
	tfxAnimationFlags_auto_set_length = 1 << 4,
	tfxAnimationFlags_orthographic = 1 << 5
} tfx_animation_flag_bits;

typedef enum {
	tfxAnimationInstanceFlags_none = 0,
	tfxAnimationInstanceFlags_loop = 1 << 0,
} tfx_animation_instance_flag_bits;

typedef enum {
	tfxAnimationManagerFlags_none = 0,
	tfxAnimationManagerFlags_has_animated_shapes = 1 << 0,
	tfxAnimationManagerFlags_initialised = 1 << 1,
	tfxAnimationManagerFlags_is_3d = 1 << 2,
} tfx_animation_manager_flag_bits;

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

//Look up frequency determines the resolution of graphs that are compiled into look up arrays.
static float tfxLOOKUP_FREQUENCY = 10.f;
//Overtime frequency is for lookups that will vary in length depending on the lifetime of the particle. It should generally be a higher resolution than the base graphs
static float tfxLOOKUP_FREQUENCY_OVERTIME = 1.f;

//Look up frequency determines the resolution of graphs that are compiled into look up arrays.
static tfxWideArray tfxLOOKUP_FREQUENCY_WIDE = tfxWideSetConst(10.f);
//Overtime frequency is for lookups that will vary in length depending on the lifetime of the particle. It should generally be a higher resolution than the base graphs
//Experiment with lower resolution and use interpolation instead? Could be a lot better on the cache.
static tfxWideArray tfxLOOKUP_FREQUENCY_OVERTIME_WIDE = tfxWideSetConst(1.f);

//-----------------------------------------------------------
//Section: forward_declarations
//-----------------------------------------------------------
#define tfxMAKE_HANDLE(handle) typedef struct handle##_s* handle;

//For allocating a new object with handle. Only used internally.
#define tfxNEW(type) (type)tfxALLOCATE(sizeof(type##_t))
#define tfxNEW_ALIGNED(type, alignment) (type)tfxALLOCATE_ALIGNED(sizeof(type##_t), alignment)

typedef struct tfx_package_s tfx_package_t;

tfxMAKE_HANDLE(tfx_package)

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
		current_size = text_len + 1;
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
	inline void SanitizeLineFeeds() { if (current_size > 1) { while (current_size > 1 && back() == '\n' || back() == '\r' || back() == '\0') { pop(); if (current_size <= 1) { break; } } NullTerminate(); } }
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
//Credit to ocornut https://github.com/ocornut/imgui/commits?author=ocornut for tfxvec although it's quite a lot different now.
//std::vector replacement with some extra stuff and tweaks specific to TimelineFX
template<typename T>
struct tfx_vector_t {
	tfxU32 current_size;
	tfxU32 capacity;
	tfxU32 volatile locked;
	tfxU32 alignment;
	T *data;

	// Provide standard typedefs but we don't use them ourselves.
	typedef T value_type;
	typedef value_type *iterator;
	typedef const value_type *const_iterator;

	inline					tfx_vector_t(const tfx_vector_t<T> &src) { locked = false; current_size = capacity = alignment = 0; data = nullptr; resize(src.current_size); memcpy(data, src.data, (size_t)current_size * sizeof(T)); }
	inline					tfx_vector_t() : locked(0), current_size(0), capacity(0), alignment(0), data(nullptr) {}
	inline					tfx_vector_t<T> &operator=(const tfx_vector_t<T> &src) { TFX_ASSERT(0); return *this; }	//Use copy instead. 
	inline					~tfx_vector_t() { TFX_ASSERT(data == nullptr); } //You must manually free containers!

	inline void				init() { locked = false; current_size = capacity = alignment = 0; data = nullptr; }
	inline bool				empty() { return current_size == 0; }
	inline bool				full() { return current_size == capacity; }
	inline tfxU32			size() { return current_size; }
	inline const tfxU32		size() const { return current_size; }
	inline tfxU32			size_in_bytes() { return current_size * sizeof(T); }
	inline const tfxU32		size_in_bytes() const { return current_size * sizeof(T); }
	inline T &operator[](tfxU32 i) { TFX_ASSERT(i < current_size); return data[i]; }
	inline const T &operator[](tfxU32 i) const { TFX_ASSERT(i < current_size); return data[i]; }
	inline T &ts_at(tfxU32 i) { while (locked > 0); return data[i]; }

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
		if (alignment != 0) {
			new_data = (T *)tfxALLOCATE_ALIGNED((size_t)new_capacity * sizeof(T), alignment);
		} else {
			new_data = (T *)tfxALLOCATE((size_t)new_capacity * sizeof(T));
		}
		TFX_ASSERT(new_data);    //Unable to allocate memory. todo: better handling
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
		current_size++; return data[current_size - 1];
	}
	inline T &push_back_copy(const T &v) {
		if (current_size == capacity) {
			reserve(_grow_capacity(current_size + 1));
		}
		memcpy(&data[current_size], &v, sizeof(v));
		current_size++; return data[current_size - 1];
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
//This simple container struct was created for storing instance_data in the particle manager. I didn't want this templated because either 2d or 3d instance_data could be used so
//I wanted to cast as needed when writing and using the sprite data. See simple cast macros above tfxCastBuffer and tfxCastBufferRef
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
	inline void         resize(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
	inline void         resize_bytes(tfxU32 new_size) { if (new_size > capacity) reserve(_grow_capacity(new_size)); current_size = new_size; }
	inline tfxU32		size_in_bytes() { return current_size * struct_size; }
	inline const tfxU32	size_in_bytes() const { return current_size * struct_size; }
};
inline tfx_buffer_t tfxCreateBuffer(tfxU32 struct_size, tfxU32 alignment) {
	tfx_buffer_t buffer;
	buffer.struct_size = struct_size;
	buffer.alignment = alignment;
	return buffer;
}
inline void tfxReconfigureBuffer(tfx_buffer_t *buffer, size_t new_struct_size) {
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
	void(*remove_callback)(T &item) = nullptr;

	tfx_storage_map_t() {}

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
	inline void InsertWithLength(const char *name, tfxU32 length, const T &value) {
		tfxKey key = tfx_Hash(&hasher, name, length, 0);
		SetIndex(key, value);
	}

	//Insert a new T value into the storage
	inline void Insert(tfxKey key, const T &value) {
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
		if (remove_callback)
			remove_callback(data[it->index]);
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
		if (remove_callback)
			remove_callback(data[it->index]);
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
		if (remove_callback)
			remove_callback(data[it->index]);
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

	void SetIndex(tfxKey key, const T &value) {
		pair *it = LowerBound(key);
		if (it == map.end() || it->key != key)
		{
			data.push_back(value);
			map.insert(it, pair(key, data.current_size - 1));
			return;
		}
		data[it->index] = value;
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
inline void tfx__reset_soa_buffer(tfx_soa_buffer_t *buffer) {
	buffer->current_arena_size = 0;
	buffer->struct_size = 0;
	buffer->current_size = 0;
	buffer->start_index = 0;
	buffer->last_bump = 0;
	buffer->capacity = 0;
	buffer->alignment = 4;
	buffer->block_size = tfxDataWidth;
	buffer->user_data = nullptr;
	buffer->resize_callback = nullptr;
	buffer->struct_of_arrays = nullptr;
	buffer->data = nullptr;
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
inline void tfx__finish_soa_buffer_setup(tfx_soa_buffer_t *buffer, void *struct_of_arrays, tfxU32 reserve_amount, tfxU32 alignment = 4) {
	TFX_ASSERT(buffer->data == nullptr && buffer->array_ptrs.current_size > 0);    //Must be an unitialised soa buffer
	TFX_ASSERT(alignment >= 4);        //Alignment must be 4 or greater
	for (int i = 0; i != buffer->array_ptrs.current_size; ++i) {
		buffer->struct_size += buffer->array_ptrs[i].unit_size;
	}
	reserve_amount = (reserve_amount / buffer->block_size + 1) * buffer->block_size;
	buffer->capacity = reserve_amount;
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
	tfxU32 new_size = ++buffer->current_size;
	if (grow && new_size == buffer->capacity) {
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
	tfxU32 first_new_index = buffer->current_size;
	tfxU32 new_size = buffer->current_size += amount;
	if (grow && new_size >= buffer->capacity) {
		grew = tfx__grow_soa_arrays(buffer, buffer->capacity, new_size);
	}
	buffer->current_size = new_size;
	TFX_ASSERT(buffer->current_size < buffer->capacity);    //Capacity of buffer is exceeded, set grow to true or don't exceed the capacity
	return first_new_index;
}

//Increase current size of a SoA Buffer and grow if grow is true. Returns the index where the new rows start.
inline tfxU32 tfx__add_soa_rows(tfx_soa_buffer_t *buffer, tfxU32 amount, bool grow) {
	TFX_ASSERT(buffer->data);            //No data allocated in buffer
	tfxU32 first_new_index = buffer->current_size;
	tfxU32 new_size = buffer->current_size + amount;
	if (grow && new_size >= buffer->capacity) {
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
	if (buffer->data)
		tfxFREE(buffer->data);
	buffer->array_ptrs.free();
	tfx__reset_soa_buffer(buffer);
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
	tfxU32 new_capacity = buffer->current_size;
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
	inline T &back() { TFX_ASSERT(current_size > 0); return bucket_list[(current_size - 1) / size_of_each_bucket]->data.back(); }
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
	inline bool EoF() { return position >= size; }
	inline void AddReturn() { if (size + 1 >= capacity) { tfxU64 new_capacity = capacity * 2; Reserve(new_capacity); } data[size] = '\n'; size++; }
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
	struct pair {
		tfxKey key;
		tfxU32 index;
	};

	tfx_hasher_t hasher;
	tfx_vector_t map;
	tfx_vector_t data;
	void(*remove_callback)(void *item);
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
extern const tfxU32 tfxPROFILE_COUNT;

extern int tfxNumberOfThreadsInAdditionToMain;

#ifndef tfxMAX_QUEUES
#define tfxMAX_QUEUES 64
#endif

#ifndef tfxMAX_QUEUE_ENTRIES
#define tfxMAX_QUEUE_ENTRIES 512
#endif

#ifdef _WIN32
#include <process.h>
#else
#include <pthread.h>
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

// Platform-specific synchronization wrapper
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

typedef struct tfx_queue_processor_s {
	tfx_sync_t sync;
	tfx_uint count;
	volatile bool end_all_threads;
	tfx_work_queue_t *queues[tfxMAX_QUEUES];
} tfx_queue_processor_t;

typedef struct tfx_data_types_dictionary_s {
	int initialised;
	void *names_and_types;
} tfx_data_types_dictionary_t;

tfxAPI_EDITOR void tfx__initialise_dictionary(tfx_data_types_dictionary_t *dictionary);

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
//There is a single thread pool created to serve multiple queues. Currently each particle manager that you create will have it's own queue and then
//each emitter that the particle manager uses will be given it's own thread.

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
tfxINTERNAL inline unsigned WINAPI tfx__thread_worker(void *arg) {
#else
inline void *tfx__thread_worker(void *arg) {
#endif
	tfx_queue_processor_t *queue_processor = (tfx_queue_processor_t *)arg;
	while (!tfx__do_next_work_queue(queue_processor)) {
		// Continue processing
	}
	return 0;
}

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
tfxAPI inline unsigned int tfx_HardwareConcurrencySafe(void) {
	unsigned int count = tfx_HardwareConcurrency();
	return count > 0 ? count : 1;
}

// Helper function to get a good default thread count for thread pools
// Usually hardware threads - 1 to leave a core for the OS/main thread
tfxAPI inline unsigned int tfx_GetDefaultThreadCount(void) {
	unsigned int count = tfx_HardwareConcurrency();
	return count > 1 ? count - 1 : 1;
}
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
tfxINTERNAL tfx_vec2_t tfx__rotate_vector_quaternion2d(const tfx_quaternion_t * q, const tfx_vec2_t v);
tfxAPI_EDITOR tfx_vec3_t tfx__rotate_vector_quaternion(const tfx_quaternion_t * q, tfx_vec3_t v);
tfxINTERNAL tfx_quaternion_t tfx__normalize_quaternion(tfx_quaternion_t * q);
tfxINTERNAL tfx_quaternion_t tfx__euler_to_quaternion(float pitch, float yaw, float roll);
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

#define tfxPACKED_X_NORMAL_3D 0x3FE7FDFF
#define tfxPACKED_Y_NORMAL_3D 0x1FFFF9FF
#define tfxPACKED_Z_NORMAL_3D 0x1FF7FFFE
#define tfxPACKED_Y_NORMAL_2D 32767

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

#ifdef tfxINTEL
const tfx128Array tfxF3_4 = tfx128SetConst(1.0f / 3.0f);
const tfx128Array tfxF2_4 = tfx128SetConst(.366025403f);
const tfx128Array tfxG2_4 = tfx128SetConst(0.211324865f);
const tfx128Array tfxG2_4x2 = tfx128SetConst(0.42264973f);
const tfx128Array tfxG3_4 = tfx128SetConst(1.0f / 6.0f);
const tfx128Array tfxG32_4 = tfx128SetConst((1.0f / 6.0f) * 2.f);
const tfx128Array tfxG33_4 = tfx128SetConst((1.0f / 6.0f) * 3.f);
const tfx128iArray tfxONE = tfx128SetConst(1);
const tfx128Array tfxONEF = tfx128SetConst(1.f);
const tfx128Array tfxZERO = tfx128SetConst(0.f);
const tfx128Array tfxTHIRTYTWO = tfx128SetConst(32.f);
const tfx128iArray tfxFF = tfx128SetConst(0xFF);
const tfx128Array tfxPSIX = tfx128SetConst(0.6f);
#elif defined(tfxARM)
const tfx128 tfxF3_4 = vdupq_n_f32(1.0f / 3.0f);
const tfx128 tfxF2_4 = vdupq_n_f32(.366025403f);
const tfx128 tfxG2_4 = vdupq_n_f32(0.211324865f);
const tfx128 tfxG2_4x2 = vdupq_n_f32(0.42264973f);
const tfx128 tfxG3_4 = vdupq_n_f32(1.0f / 6.0f);
const tfx128 tfxG32_4 = vdupq_n_f32((1.0f / 6.0f) * 2.f);
const tfx128 tfxG33_4 = vdupq_n_f32((1.0f / 6.0f) * 3.f);
const tfx128i tfxONE = vdupq_n_s32(1);
const tfx128 tfxONEF = vdupq_n_f32(1.f);
const tfx128 tfxZERO = vdupq_n_f32(0.f);
const tfx128 tfxTHIRTYTWO = vdupq_n_f32(32.f);
const tfx128i tfxFF = vdupq_n_s32(0xFF);
const tfx128 tfxPSIX = vdupq_n_f32(0.6f);
#endif

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
tfx128Array tfxNoise4_2d(const tfx128Array x4, const tfx128Array y4);
tfx128Array tfxNoise4_3d(const tfx128Array x4, const tfx128Array y4, const tfx128Array z4);

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

#ifdef __cplusplus
	tfx_attribute_node_s() : frame(0.f), value(0.f), flags(0), index(0) { }
	inline bool operator==(const struct tfx_attribute_node_s &n) { return n.frame == frame && n.value == value; }
#endif
}tfx_attribute_node_t;

typedef struct tfx_depth_index_s {
	tfxParticleID particle_id;
	float depth;
}tfx_depth_index_t;

typedef struct tfx_graph_lookup_t {
#ifdef __cplusplus
	tfx_vector_t<float> values;
#else
	tfx_vector_t values;
#endif
	tfxU32 last_frame;
	float life;
} tfx_graph_lookup_t;

//Used when a particle manager is grouping instances by effect. This way effects can be individually ordered and drawn/not drawn in order however you need
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
} tfx_effect_instance_data_t;

typedef struct tfx_face_s {
	int v[3];
} tfx_face_t;

typedef struct tfx_random_s {
	tfxU64 seeds[2];
}tfx_random_t;

typedef struct tfx_color_ramp_s {
	//These vectors are for sinusoidal color ramp generation
	tfx_vec3_t brightness;
	tfx_vec3_t contrast;
	tfx_vec3_t frequency;
	tfx_vec3_t offsets;
	tfxColorRampFlags flags;
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
}tfx_bitmap_t;

typedef struct tfx_graph_id_s {
	tfx_graph_category category;
	tfx_graph_type type;
	tfxU32 graph_id;
	tfxU32 node_id;
	tfxKey path_hash;
}tfx_graph_id_t;

typedef struct tfx_graph_lookup_index_s {
	tfxU32 start_index;
	tfxU32 length;
	float max_life;
	float padding1;
}tfx_graph_lookup_index_t;

//This typedef struct is used to store indexing data in order to index into large lists containing either the node data of graphs
//or the lookup data of compiled graphs. This is so that we can upload that data into a buffer on the GPU to get the particles
//updating in a compute shader.
typedef struct tfx_effect_lookup_data_s {
	tfx_graph_lookup_index_t overtime_velocity;
	tfx_graph_lookup_index_t overtime_width;
	tfx_graph_lookup_index_t overtime_height;
	tfx_graph_lookup_index_t overtime_weight;
	tfx_graph_lookup_index_t overtime_roll_spin;
	tfx_graph_lookup_index_t overtime_pitch_spin;
	tfx_graph_lookup_index_t overtime_yaw_spin;
	tfx_graph_lookup_index_t overtime_stretch;
	tfx_graph_lookup_index_t overtime_red;
	tfx_graph_lookup_index_t overtime_green;
	tfx_graph_lookup_index_t overtime_blue;
	tfx_graph_lookup_index_t overtime_blendfactor;
	tfx_graph_lookup_index_t overtime_red_hint;
	tfx_graph_lookup_index_t overtime_green_hint;
	tfx_graph_lookup_index_t overtime_blue_hint;
	tfx_graph_lookup_index_t overtime_blendfactor_hint;
	tfx_graph_lookup_index_t overtime_velocity_turbulance;
	tfx_graph_lookup_index_t overtime_direction_turbulance;
	tfx_graph_lookup_index_t overtime_velocity_adjuster;
	tfx_graph_lookup_index_t overtime_intensity;
	tfx_graph_lookup_index_t overtime_hint_intensity;
	tfx_graph_lookup_index_t overtime_color_mix_overtime;
	tfx_graph_lookup_index_t overtime_direction;
	tfx_graph_lookup_index_t overtime_noise_resolution;
	tfx_graph_lookup_index_t overtime_motion_randomness;
}tfx_effect_lookup_data_t;

typedef struct tfx_graph_s {
	//The ratio to transalte graph frame/value to grid x/y coords on a graph editor
	tfx_vec2_t min;
	tfx_vec2_t max;
	tfx_graph_preset graph_preset;
	tfx_graph_type type;
	tfx_effect_emitter_t *effector;
#ifdef __cplusplus
	tfx_bucket_array_t<tfx_attribute_node_t> nodes;
#else
	tfx_bucket_array_t nodes;
#endif
	tfx_graph_lookup_t lookup;
	tfxU32 index;
	float gamma;
} tfx_graph_t;

tfxAPI_EDITOR void tfx__init_graph(tfx_graph_t * graph, tfxU32 node_bucket_size);

//The following structs group graphs together under the attribute categories Global, Transform, Properties, Base, Variation and Overtime
typedef struct tfx_global_attributes_s {
	tfx_graph_t life;
	tfx_graph_t amount;
	tfx_graph_t velocity;
	tfx_graph_t noise;
	tfx_graph_t width;
	tfx_graph_t height;
	tfx_graph_t weight;
	tfx_graph_t spin;
	tfx_graph_t pitch_spin;
	tfx_graph_t yaw_spin;
	tfx_graph_t stretch;
	tfx_graph_t overal_scale;
	tfx_graph_t intensity;
	tfx_graph_t splatter;
	tfx_graph_t emitter_width;
	tfx_graph_t emitter_height;
	tfx_graph_t emitter_depth;
}tfx_global_attributes_t;

typedef struct tfx_transform_attributes_s {
	tfx_graph_t roll;
	tfx_graph_t pitch;
	tfx_graph_t yaw;
	tfx_graph_t translation_x;
	tfx_graph_t translation_y;
	tfx_graph_t translation_z;
}tfx_transform_attributes_t;

typedef struct tfx_property_attributes_s {
	tfx_graph_t emission_pitch;
	tfx_graph_t emission_yaw;
	tfx_graph_t emission_range;
	tfx_graph_t splatter;
	tfx_graph_t emitter_width;
	tfx_graph_t emitter_height;
	tfx_graph_t emitter_depth;
	tfx_graph_t extrusion;
	tfx_graph_t arc_size;
	tfx_graph_t arc_offset;
}tfx_property_attributes_t;

typedef struct tfx_base_attributes_s {
	tfx_graph_t life;
	tfx_graph_t amount;
	tfx_graph_t velocity;
	tfx_graph_t width;
	tfx_graph_t height;
	tfx_graph_t weight;
	tfx_graph_t pitch_spin;
	tfx_graph_t yaw_spin;
	tfx_graph_t spin;
	tfx_graph_t noise_offset;
}tfx_base_attributes_t;

typedef struct tfx_variation_attributes_s {
	tfx_graph_t life;
	tfx_graph_t amount;
	tfx_graph_t velocity;
	tfx_graph_t width;
	tfx_graph_t height;
	tfx_graph_t weight;
	tfx_graph_t pitch_spin;
	tfx_graph_t yaw_spin;
	tfx_graph_t spin;
	tfx_graph_t noise_offset;
	tfx_graph_t noise_resolution;
	tfx_graph_t motion_randomness;
} tfx_variation_attributes_t;

typedef struct tfx_overtime_attributes_s {
	tfx_graph_t velocity;
	tfx_graph_t width;
	tfx_graph_t height;
	tfx_graph_t weight;
	tfx_graph_t pitch_spin;
	tfx_graph_t yaw_spin;
	tfx_graph_t spin;
	tfx_graph_t stretch;
	tfx_graph_t red;
	tfx_graph_t green;
	tfx_graph_t blue;
	tfx_graph_t blendfactor;
	tfx_graph_t red_hint;
	tfx_graph_t green_hint;
	tfx_graph_t blue_hint;
	tfx_graph_t blendfactor_hint;
	tfx_graph_t velocity_turbulance;
	tfx_graph_t direction_turbulance;
	tfx_graph_t velocity_adjuster;
	tfx_graph_t intensity;
	tfx_graph_t alpha_sharpness;
	tfx_graph_t curved_alpha;
	tfx_graph_t direction;
	tfx_graph_t noise_resolution;
	tfx_graph_t motion_randomness;
	tfx_color_ramp_t color_ramps[2];
	tfx_index color_ramp_bitmap_indexes[2];
} tfx_overtime_attributes_t;

typedef struct tfx_factor_attributes_s {
	tfx_graph_t life;
	tfx_graph_t size;
	tfx_graph_t velocity;
	tfx_graph_t intensity;
} tfx_factor_attributes_t;

typedef struct tfx_path_nodes_soa_s {
	float *x;
	float *y;
	float *z;
	float *length;
}tfx_path_nodes_soa_t;

typedef struct tfx_path_quaternion_s {
	tfxU32 quaternion;
	float grid_coord;
	float age;
	tfxU32 cycles;
}tfx_path_quaternion_t;

typedef struct tfx_emitter_path_s {
	tfxKey key;
	tfx_str32_t name;
	int node_count;
	tfxEmitterPathFlags flags;
	tfx_path_generator_type generator_type;
	tfx_graph_t angle_x;
	tfx_graph_t angle_y;
	tfx_graph_t angle_z;
	tfx_graph_t offset_x;
	tfx_graph_t offset_y;
	tfx_graph_t offset_z;
	tfx_graph_t distance;
	float rotation_range;
	union {
		float rotation_pitch;    //3d paths
		float rotation_offset;    //For 2d paths
	};
	float rotation_yaw;
	tfxU32 maximum_active_paths;
	tfxU32 maximum_paths;
	float rotation_cycle_length;
	float rotation_stagger;
	tfx_vec3_t offset;
	tfx_vec3_t builder_parameters;
#ifdef __cplusplus
	tfx_vector_t<tfx_vec4_t> nodes;
#else
	tfx_vector_t nodes;
#endif
	tfx_soa_buffer_t node_buffer;
	tfx_path_nodes_soa_t node_soa;
	tfx_path_extrusion_type extrusion_type;
} tfx_emitter_path_t;

typedef struct tfx_emitter_attributes_s {
	tfx_property_attributes_t properties;
	tfx_base_attributes_t base;
	tfx_variation_attributes_t variation;
	tfx_overtime_attributes_t overtime;
	tfx_factor_attributes_t factor;
}tfx_emitter_attributes_t;

static float(*lookup_overtime_callback)(tfx_graph_t * graph, float age, float lifetime);
static float(*lookup_callback)(tfx_graph_t * graph, float age);
static float(*lookup_random_callback)(tfx_graph_t * graph, float age, tfx_random_t * random);

typedef struct tfx_shape_data_s {
	tfx_str64_t name;
	tfxU32 frame_count;
	tfxU32 width;
	tfxU32 height;
	tfxU32 shape_index;
	tfxKey image_hash;
	int import_filter;
}tfx_shape_data_t;

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
}tfx_camera_settings_t;

typedef struct tfx_preview_camera_settings_s {
	tfx_camera_settings_t camera_settings;
	float effect_z_offset;
	float camera_speed;
	bool attach_effect_to_camera;
}tfx_preview_camera_settings_t;

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
	tfxU32 needs_exporting;
	float max_radius;
	tfxU32 largest_frame;
	float playback_speed;
	float effect_z_offset;
	tfx_export_color_options color_option;
	tfx_export_options export_option;
	tfx_camera_settings_t camera_settings;
	tfx_camera_settings_t camera_settings_orthographic;
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
	tfxU32 needs_exporting;
	float max_radius;
	tfxU32 largest_frame;
	float playback_speed;
	float recording_frame_rate;
}tfx_sprite_data_settings_t;

//------------------------------------------------------------

//API structs you can access in various ways to update and render effects in realtime

//Image data for particle shapes. This is passed into your custom ShapeLoader function for loading image textures into whatever renderer you're using
typedef struct tfx_image_data_s {
	//This can be a ptr to the image texture for rendering. You must assign this in your ShapeLoader function
	void *ptr;
	//Index of the image, deprecated, image hash should be used now instead.
	tfxU32 shape_index;
	//Name of the image
	tfx_str64_t name;
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
}tfx_image_data_t;

typedef struct tfx_emitter_properties_s {
	//Angle added to the rotation of the particle when spawned or random angle range if angle setting is set to tfx_random_t
	tfx_vec3_t angle_offsets;
	//When aligning the billboard along a vector, you can set the type of vector that it aligns with
	tfx_vector_align_type vector_align_type;
	//Point, area, ellipse emitter etc.
	tfx_emission_type emission_type;
	//For other emitter emission types, this hash is the location of the other emitter so that it can be used to connect the two
	//emitters together when added to a particle manager.
	tfxKey paired_emitter_hash;
	//If single shot flag is set then you can limit how many times it will loop over it's overtime graphs before expiring
	tfxU32 single_shot_limit;
	//Animation frame rate
	float frame_rate;
	//The final frame index of the animation
	float end_frame;
	//Pointer to the ImageData in the EffectLibary. 
	tfx_image_data_t *image;
	//For 3d effects, the type of billboarding: 0 = use billboarding (always face camera), 1 = No billboarding, 2 = No billboarding and align with motion
	tfx_billboarding_option billboard_option;

	//The number of rows/columns/ellipse/line points in the grid when spawn on grid flag is used
	tfx_vec3_t grid_points;
	//The rotation of particles when they spawn, or behave overtime if tfxAlign is used
	tfxAngleSettingFlags angle_settings;
	//Layer of the particle manager that the particle is added to
	tfxU32 layer;
	//Milliseconds to delay spawing
	float delay_spawning;
	//Should particles emit towards the center of the emitter or away, or in a specific direction
	tfx_emission_direction emission_direction;

	//How particles should behave when they reach the end of the line
	tfx_line_traversal_end_behaviour end_behaviour;
	//Bit field of various boolean state_flags
	tfxParticleControlFlags compute_flags;
	//Offset to draw particles at
	tfx_vec2_t image_handle;
	//image handle packed into 16bit floats
	tfxU32 image_handle_packed;
	//Offset of emitters and effects
	tfx_vec3_t emitter_handle;
	//When single flag is set, spawn this amount of particles in one go
	tfxU32 spawn_amount;
	//When single flag is set, spawn this variable amount of particles in one go
	tfxU32 spawn_amount_variation;
	//The shape being used for all particles spawned from the emitter (deprecated, hash now used instead)
	tfxU32 image_index;
	//The shape being used for all particles spawned from the emitter
	tfxKey image_hash;
	//The number of millisecs before an effect or emitter will loop back round to the beginning of it's graph lookups
	float loop_length;
	//The start frame index of the animation
	float start_frame;
	//Base noise offset random range so that noise patterns don't repeat so much over multiple effects
	float noise_base_offset_range;
	//This is only used for the animation manager when sprite data is added to the animation manager. This is used to map
	//the property_index to the animation property index so the sprite data can point to a new index where some emitter properties
	//are stored on the GPU for looking up from the sprite data
	tfxU32 animation_property_index;
}tfx_emitter_properties_t;

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

//This is a struct that stores an emitter state that is currently active in a particle manager.
//Todo: maybe split this up into static variables that stay the same (they're just properties copied from the emitter in the library
//        and dynamic variables that change each frame.
typedef struct tfx_emitter_state_s {
	//State data
	float frame;
	float age;
	float highest_particle_age;
	float delay_spawning;
	float timeout_counter;
	float timeout;
	float amount_remainder;
	float spawn_quantity;
	float qty_step_size;
	tfx_vec3_t handle;
	tfxEmitterPropertyFlags property_flags;
	//Position, scale and rotation values
	tfx_vec3_t local_position;
	tfx_vec3_t world_position;
	tfx_vec3_t captured_position;
	tfx_vec3_t world_rotations;

	float loop_length;
	tfx_quaternion_t rotation;
	tfxU64 image_handle_packed;
	tfx_bounding_box_t bounding_box;

	tfxU32 emitter_attributes;
	tfxU32 transform_attributes;
	tfxU32 overtime_attributes;
	tfxU32 path_attributes;
	tfx_path_quaternion_t *path_quaternions;
	tfxU32 path_quaternion_index;
	tfxU32 last_path_index;
	float path_stagger_counter;
	tfxU32 path_cycle_count;
	tfxU32 active_paths;
	tfxU32 path_start_index;

	tfxU32 root_index;
	tfxU32 parent_index;
	tfxU32 properties_index;
	tfxU32 info_index;
	tfxU32 hierarchy_depth;
	tfxU32 sprites_count;
	tfxU32 sprites_index;
	tfxU32 seed_index;
	tfxKey path_hash;

	//Control Data
	tfxU32 particles_index;
	tfxU32 spawn_locations_index;    //For other_emitter emission type and storing the last known position of the particle
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
}tfx_emitter_state_t TFX_ALIGN_AFFIX(16);

//This is a struct that stores an effect state that is currently active in a particle manager.
typedef struct tfx_effect_state_s {
	tfx_quaternion_t rotation;
	//State data
	float frame;
	float age;
	float highest_particle_age;
	float timeout_counter;
	float timeout;
	tfx_vec3_t handle;
	tfxEmitterPropertyFlags property_flags;
	tfxEffectPropertyFlags effect_flags;
	float loop_length;
	//Position, scale and rotation values
	tfx_vec3_t translation;
	tfx_vec3_t local_position;
	tfx_vec3_t world_position;
	tfx_vec3_t captured_position;
	tfx_vec3_t local_rotations;
	tfx_vec3_t world_rotations;
	tfx_bounding_box_t bounding_box;

	tfxU32 global_attributes;
	tfxU32 transform_attributes;

	tfxU32 properties_index;
	tfxU32 info_index;
	tfxU32 parent_particle_index;
	tfx_library_t *library;
	tfxKey path_hash;

	//Spawn controls
	tfx_parent_spawn_controls_t spawn_controls;
	tfx_vec3_t emitter_size;
	float stretch;
	float noise;
	float overal_scale;
	float noise_base_offset;
	tfxEmitterStateFlags state_flags;
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

	//User Data
	void *user_data;
	void(*update_callback)(tfx_particle_manager_t *pm, tfxEffectID effect_index);
}tfx_effect_state_t TFX_ALIGN_AFFIX(16);

//An tfx_effect_emitter_t can either be an effect which stores effects and global graphs for affecting all the attributes in the emitters
//Or it can be an emitter which spawns all of the particles. 
//This is only for library storage, when using to update each frame this is copied to tfx_effect_state_t and tfx_emitter_state_t for realtime updates
typedef struct tfx_effect_emitter_s {
	//Required for frame by frame updating
	//The current state of the effect/emitter used in the editor only at this point
	tfxEmitterStateFlags state_flags;
	//Property flags for emitters
	tfxEmitterPropertyFlags property_flags;
	//Flags specific to effects
	tfxEffectPropertyFlags effect_flags;
	//A link to the library that this effect/emitter belongs to
	tfx_library_t *library;
	//Is this an tfxEffectType or tfxEmitterType
	tfx_effect_emitter_type type;
	//The index within the library that this exists at
	tfxU32 library_index;
	//A hash of the directory path to the effect ie Flare/spark, and also a UID for the effect/emitter
	tfxKey path_hash;
	//All graphs that the effect uses to lookup attribute values are stored in the library. These variables here are indexes to the array where they're stored
	tfxU32 global;
	tfxU32 emitter_attributes;
	tfxU32 transform_attributes;
	tfxU32 path_attributes;
	//The type of function that should be called to update particle positions
	tfxEmitterControlProfileFlags control_profile;
	//Pointer to the immediate parent
	struct tfx_effect_emitter_s *parent;
	//When not using insert sort to guarantee particle order, sort passes offers a more relaxed way of ordering particles over a number of frames.
	//The more passes the more quickly ordered the particles will be but at a higher cost
	tfxU32 sort_passes;
	//Custom user data, can be accessed in callback functions
	void *user_data;
	void(*update_callback)(tfx_particle_manager_t *pm, tfxEffectID effect_index);

	tfxU32 buffer_index;

	//Indexes into library storage
	tfxU32 info_index;
	tfxU32 property_index;
} tfx_effect_emitter_t;

typedef struct tfx_effect_emitter_info_s {
	//Name of the effect
	tfx_str64_t name;
	//The path of the effect in the library
	tfx_str512_t path;
	//Every effect and emitter in the library gets a unique id
	tfxU32 uid;
	//The max_radius of the emitter, taking into account all the particles that have spawned and active (editor only)
	float max_radius;
	//Experiment: index into the lookup index data in the effect library
	tfxU32 lookup_node_index;
	tfxU32 lookup_value_index;
	//Index to sprite sheet settings stored in the effect library. 
	tfxU32 sprite_sheet_settings_index;
	//Index to sprite data settings stored in the effect library. 
	tfxU32 sprite_data_settings_index;
	//Index to preview camera settings stored in the effect library. Would like to move this at some point
	tfxU32 preview_camera_settings;
	//The maximum amount of life that a particle can be spawned with taking into account base + variation life values
	float max_life;
	//List of sub_effects ( effects contain emitters, emitters contain sub effects )
#ifdef __cplusplus
	tfx_vector_t<tfx_effect_emitter_t> sub_effectors;
#else
	tfx_vector_t sub_effectors;
#endif
}tfx_effect_emitter_info_t;

typedef struct tfx_compute_sprite_s {    //64 bytes
	tfx_vec4_t bounds;                //the min/max x,y coordinates of the image being drawn
	tfx_vec4_t uv;                    //The UV coords of the image in the texture
	tfx_vec4_t scale_rotation;        //Scale and rotation (x, y = scale, z = rotation, w = multiply blend factor)
	tfx_vec2_t position;            //The position of the sprite
	tfx_rgba8_t color;                //The color tint of the sprite
	tfxU32 parameters;    //4 extra parameters packed into a tfxU32: blend_mode, image layer index, shader function index, blend type
}tfx_compute_sprite_t;

typedef struct tfx_unique_sprite_id_s {
	tfxU32 uid;
	tfxU32 age;
	tfxU32 property_index;
}tfx_unique_sprite_id_t;

//These all point into a tfx_soa_buffer_t, initialised with InitParticleSoA. Current Bandwidth: 108 bytes
//Note that not all of these are used, it will depend on the emitter and which attributes it uses. So to save memory,
//when the the buffer is initialised only the fields that are needed for the emitter will be used.
typedef struct tfx_particle_soa_s {
	tfxU32 *uid;
	tfxU32 *sprite_index;
	tfxU32 *particle_index;
	tfxParticleFlags *flags;
	float *age;
	float *max_age;
	float *position_x;
	float *position_y;
	float *position_z;
	float *captured_position_x;
	float *captured_position_y;
	float *captured_position_z;
	float *local_rotations_x;    //In 2d this is the direction
	float *local_rotations_y;
	float *local_rotations_z;    //In 2d this is the roll
	tfxU32 *velocity_normal;
	tfxU32 *quaternion;            //Used for paths where the path can be rotated per particle based on the emission direction
	tfxU32 *depth_index;
	float *path_position;
	float *path_offset;
	float *base_weight;
	float *base_velocity;
	float *base_spin;
	float *base_pitch_spin;
	float *base_yaw_spin;
	float *base_size_x;
	float *base_size_y;
	float *noise_offset;
	float *noise_resolution;
	float *intensity_factor;
	tfx_rgba8_t *color;
	float *image_frame;
	tfxU32 *single_loop_count;
}tfx_particle_soa_t;

typedef struct tfx_spawn_points_soa_s {
	float *position_x;
	float *position_y;
	float *position_z;
	float *captured_position_x;
	float *captured_position_y;
	float *captured_position_z;
	float *age;
}tfx_spawn_points_soa_t;

typedef struct tfx_sprite_transform2d_s {
	tfx_vec2_t position;							//The position of the sprite, x, y - world, z, w = captured for interpolating
	tfx_vec2_t scale;								//Scale
	float rotation;
}tfx_sprite_transform2d_t;

typedef struct tfx_sprite_transform3d_s {
	tfx_vec3_t position;							//The position of the sprite, x, y - world, z, w = captured for interpolating
	tfx_vec3_t rotations;							//Rotations of the sprite
	tfx_vec2_t scale;								//Scale
}tfx_sprite_transform3d_t;

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
}tfx_frame_meta_t;

//This is the exact typedef struct to upload to the GPU for 2d instance_data, so timelinefx will prepare buffers so that they're ready to just
//upload to the GPU in one go. Of course you don't *have* to do this you could loop over the buffer and draw the instance_data
//in a different way if you don't have this option for some reason, but the former way is by far the most efficient.
typedef struct tfx_2d_instance_s {			//44 bytes + padding to 48
	tfx_vec4_t position;							//The position of the sprite, rotation in w, stretch in z
	tfx_float16x4_t size_handle;					//Size of the sprite in pixels and the handle packed into a u64 (4 16bit floats)
	tfx_float16x2_t alignment;						//normalised alignment vector 2 floats packed into 16bits or 3 8bit floats for 3d
	tfx_float16x2_t intensity_life;					//Multiplier for the color and life of particle
	tfx_float16x2_t curved_alpha;					//Sharpness and dissolve amount value for fading the image 2 16bit floats packed
	tfxU32 indexes;									//[color ramp y index, color ramp texture array index, capture flag, image data index (1 bit << 15), billboard alignment (2 bits << 13), image data index max 8191 images]
	tfxU32 captured_index;							//Index to the sprite in the buffer from the previous frame for interpolation
	float lerp_offset;
}tfx_2d_instance_t;

typedef struct tfx_3d_instance_s {		//56 bytes + padding to 64
	tfx_vec4_t position;							//The position of the billboard
	tfx_vec3_t rotations;				            //Rotations of the billboard with stretch in w
	tfx_float8x4_t alignment;						//normalised alignment vector 2 floats packed into 16bits or 3 8bit floats for 3d
	tfx_float16x4_t size_handle;					//Size of the sprite in pixels and the handle packed into a u64 (4 16bit floats)
	tfx_float16x2_t intensity_life;					//Multiplier for the color and life of particle
	tfx_float16x2_t curved_alpha;					//Sharpness and dissolve amount value for fading the image 2 16bit floats packed
	tfxU32 indexes;									//[color ramp y index, color ramp texture array index, capture flag, image data index (1 bit << 15), billboard alignment (2 bits << 13), image data index max 8191 images]
	tfxU32 captured_index;							//Index to the sprite in the buffer from the previous frame for interpolation
	float lerp_offset;
	tfxU32 padding;
}tfx_3d_instance_t;

//This typedef struct of arrays is used for both 2d and 3d instance_data, but obviously the transform_3d data is either 2d or 3d depending on which effects you're using in the particle manager.
//InitSprite3dSoA is called to initialise 3d instance_data and InitSprite2dArray for 2d instance_data. This is all managed internally by the particle manager. It's convenient to have both 2d and
//3d in one typedef struct like this as it makes it a lot easier to use the same control functions where we can. Also note that stretch and alignment for 3d instance_data are packed into
//stretch_alignment_x and alignment_yz as 16bit floats. 2d uses the float for stretch and packs xy alignment into alignment_yz
typedef struct tfx_sprite_soa_s {                       //3d takes 56 bytes of bandwidth, 2d takes 40 bytes of bandwidth
	tfxU32 *property_indexes;                   //The image frame of animation index packed with alignment option flag and property_index
	tfxU32 *captured_index;                     //The index of the sprite in the previous frame so that it can be looked up and interpolated with
	tfx_unique_sprite_id_t *uid;                //Unique particle id of the sprite, only used when recording sprite data
	tfx_sprite_transform3d_t *transform_3d;     //Transform data for 3d instance_data
	tfx_sprite_transform2d_t *transform_2d;     //Transform data for 2d instance_data
	tfxU32 *intensity_life;                     //The multiplier for the sprite color and the lifetime of the particle (0..1)
	float *stretch;                             //Multiplier for how much the particle is stretched in the shader
	tfxU32 *alignment;                          //The alignment of the particle. 2 16bit floats for 2d and 3 8bit floats for 3d
	tfxU32 *curved_alpha;						//Alpha sharpness for the texture dissolve/alpha calculations
	tfxU32 *indexes;							//The indexes to lookup the color in the color ramp textures and the image texture data
}tfx_sprite_soa_t;

typedef enum {
	tfxSpriteBufferMode_2d,
	tfxSpriteBufferMode_3d,
	tfxSpriteBufferMode_both,
} tfxSpriteBufferMode;

//These structs are for animation sprite data that you can upload to the gpu
typedef struct tfx_sprite_data3d_s {    //60 bytes aligning to 64
	tfx_vec4_t position_stretch;                    //The position of the sprite, x, y - world, z, w = captured for interpolating
	tfx_vec3_t rotations;				            //Rotations of the sprite
	tfx_float8x4_t alignment;						//normalised alignment vector 3 floats packed into 8bits
	tfx_float16x4_t size_handle;					//Size of the sprite in pixels and the handle packed into a u64 (4 16bit floats)
	tfx_float16x2_t intensity_life;					//Multiplier for the color and life of particle
	tfx_float16x2_t curved_alpha;					//Sharpness and dissolve amount value for fading the image 2 16bit floats packed
	tfxU32 indexes;									//[color ramp y index, color ramp texture array index, capture flag, image data index (1 bit << 15), billboard alignment (2 bits << 13), image data index max 8191 images]
	tfxU32 captured_index;							//Index to the sprite in the buffer from the previous frame for interpolation
	tfxU32 additional;								//Padding, but also used to pack lerp offset and property index
	tfxU32 padding;
}tfx_sprite_data3d_t;

typedef struct tfx_sprite_data2d_s {    //48 bytes
	tfx_vec4_t position_stretch_rotation;           //The position of the sprite, rotation in w, stretch in z
	tfx_float16x4_t size_handle;					//Size of the sprite in pixels and the handle packed into a u64 (4 16bit floats)
	tfx_float16x2_t alignment;						//normalised alignment vector 2 floats packed into 16bits or 3 8bit floats for 3d
	tfx_float16x2_t intensity_life;					//Multiplier for the color and life of particle
	tfx_float16x2_t curved_alpha;					//Sharpness and dissolve amount value for fading the image 2 16bit floats packed
	tfxU32 indexes;									//[color ramp y index, color ramp texture array index, capture flag, image data index (1 bit << 15), billboard alignment (2 bits << 13), image data index max 8191 images]
	tfxU32 captured_index;							//Index to the sprite in the buffer from the previous frame for interpolation
	tfxU32 additional;								//Padding, but also used to pack lerp offset and property index
}tfx_sprite_data2d_t;

//Animation sprite data that is used on the cpu to bake the data
typedef struct tfx_sprite_data_soa_s {
	tfx_2d_instance_t *sprite_instance;
	tfx_3d_instance_t *billboard_instance;
	tfx_unique_sprite_id_t *uid;
	float *lerp_offset;
}tfx_sprite_data_soa_t;

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
#else
	tfx_vector_t frame_meta;
#endif
	tfxAnimationManagerFlags flags;
	tfxAnimationFlags animation_flags;
}tfx_sprite_data_metrics_t;

typedef struct tfx_sprite_data_s {
	float frame_compression;
	tfx_sprite_data_metrics_t normal;
	tfx_sprite_data_metrics_t compressed;
	tfx_soa_buffer_t real_time_sprites_buffer;
	tfx_sprite_data_soa_t real_time_sprites;
	tfx_soa_buffer_t compressed_sprites_buffer;
	tfx_sprite_data_soa_t compressed_sprites;
}tfx_sprite_data_t;

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
	tfx_vec4_t scale_rotation;                //Scale and rotation (x, y = scale, z = rotation, w = velocity_adjuster)
	float end_frame;
	tfxU32 normalised_values;        //Contains normalized values which are generally either 0 or 255, normalised in the shader to 0 and 1 (except opacity): age_rate, line_negator, spin_negator, position_negator, opacity
	tfxParticleControlFlags flags;
	tfxU32 image_data_index;        //index into the shape buffer on the gpu. CopyComputeShapeData must be called to prepare the data.
	tfx_vec2_t image_handle;
	tfx_vec2_t emitter_handle;
	float noise_offset;
	float stretch;
	float frame_rate;
	float noise_resolution;
}tfx_compute_controller_t;

typedef struct tfx_compute_particle_s {
	tfx_vec2_t local_position;
	tfx_vec2_t base_size;

	float base_velocity;
	float base_spin;
	float base_weight;

	float age;                            //The age of the particle, used by the controller to look up the current state on the graphs
	float max_age;                        //max age before the particle expires
	float emission_angle;                //Emission angle of the particle at spawn time

	float noise_offset;                    //The random velocity added each frame
	float noise_resolution;                //The random velocity added each frame
	float image_frame;
	tfxU32 control_slot_and_layer;    //index to the controller, and also stores the layer in the particle manager that the particle is on (layer << 3)
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

typedef struct tfx_gpu_shapes_s {
#ifdef __cplusplus
	tfx_vector_t<tfx_gpu_image_data_t> list;
#else
	tfx_vector_t list;
#endif
}tfx_gpu_shapes_t;

typedef struct tfx_spawn_work_entry_s {
	tfx_random_t random;
	tfx_particle_manager_t *pm;
	tfx_emitter_properties_t *properties;
	tfx_parent_spawn_controls_t *parent_spawn_controls;
	tfxU32 emitter_index;
	tfxU32 parent_index;
	tfx_emission_type emission_type;
	tfxEmitterPropertyFlags property_flags;
	tfxEmitterPropertyFlags parent_property_flags;
	tfxEffectPropertyFlags root_effect_flags;
	tfx_particle_soa_t *particle_data;
#ifdef __cplusplus
	tfx_vector_t<tfx_effect_emitter_t> *sub_effects;
	tfx_vector_t<tfx_depth_index_t> *depth_indexes;
#else
	tfx_vector_t *sub_effects;
	tfx_vector_t *depth_indexes;
#endif
	tfxU32 depth_index_start;
	tfxU32 seed;
	float tween;
	tfxU32 max_spawn_count;
	tfxU32 amount_to_spawn;
	tfxU32 spawn_start_index;
	tfxU32 next_buffer;
	int depth;
	float qty_step_size;
	float highest_particle_age;
	float overal_scale;
}tfx_spawn_work_entry_t;

typedef struct tfx_control_work_entry_s {
	float node_count;
	tfxU32 start_index;
	tfxU32 end_index;
	tfxU32 wide_end_index;
	tfxU32 start_diff;
	tfxU32 sprites_index;
	tfxU32 cumulative_index_point;
	tfxU32 effect_instance_offset;
	tfxU32 sprite_buffer_end_index;
	tfxU32 emitter_index;
	tfx_particle_manager_t *pm;
	tfx_overtime_attributes_t *graphs;
	tfxU32 layer;
	tfx_emitter_properties_t *properties;
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

typedef struct tfx_particle_age_work_entry_s {
	tfxU32 start_index;
	tfxU32 emitter_index;
	tfxU32 wide_end_index;
	tfxU32 start_diff;
	tfx_emitter_properties_t *properties;
	tfx_particle_manager_t *pm;
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
	bool is_3d;
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
	tfxU32 frame_count;                    //The number of frames in the animation
	tfxU32 offset_into_sprite_data;        //The starting ofset in the buffer that contains all the sprite data
	tfxU32 info_index;                    //Index into the effect_animation_info storage map to get at the frame meta
	float current_time;                    //Current point of time in the animation
	float animation_length_in_time;        //Total time that the animation lasts for
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

//Use the animation manager to control playing of pre-recorded effects
typedef struct tfx_animation_manager_s {
#ifdef __cplusplus
	//All of the sprite data for all the animations that you might want to play on the GPU.
	//This could be deleted once it's uploaded to the GPU
	//An animation manager can only be used for either 2d or 3d not both
	tfx_vector_t<tfx_sprite_data3d_t> sprite_data_3d;
	tfx_vector_t<tfx_sprite_data2d_t> sprite_data_2d;
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
	//We also need to upload some emitter properties to the GPU as well such as the sprite handle.
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
#else
	tfx_vector_t sprite_data_3d;
	tfx_vector_t sprite_data_2d;
	tfx_vector_t instances;
	tfx_vector_t instances_in_use[2];
	tfx_vector_t free_instances;
	tfx_vector_t render_queue;
	tfx_vector_t offsets;
	tfx_vector_t emitter_properties;
	tfx_vector_t sprite_data_settings;
	tfx_storage_map_t effect_animation_info;
	tfx_storage_map_t particle_shapes;
#endif
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
	bool((*maybe_render_instance_callback)(tfx_animation_manager_t *animation_manager, tfx_animation_instance_t *instance, tfx_frame_meta_t *meta, void *user_data));
} tfx_animation_manager_t;

typedef struct tfx_effect_index_s {
	tfxEffectID index;
	float depth;
}tfx_effect_index_t;

typedef struct tfx_particle_manager_info_s {
	tfxU32 max_particles;					//The maximum number of instance_data for each layer. This setting is not relevent if dynamic_sprite_allocation is set to true or group_sprites_by_effect is true.
	tfxU32 max_effects;                     //The maximum number of effects that can be updated at the same time.
	tfx_particle_manager_mode order_mode;   //When not grouping instance_data by effect, you can set the mode of the particle manager to order instance_data or not.
	//When set to false, all instance_data will be kept together in a large list.
	tfxU32 multi_threaded_batch_size;       //The size of each batch of particles to be processed when multithreading. Must be a power of 2 and 256 or greater.
	tfxU32 sort_passes;                     //when in order by depth mode (not guaranteed order) set the number of sort passes for more accuracy. Anything above 5 and you should just be guaranteed order.
	bool double_buffer_sprites;             //Set to true to double buffer instance_data so that you can interpolate between the old and new positions for smoother animations.
	bool dynamic_sprite_allocation;         //Set to true to automatically resize the sprite buffers if they run out of space. Not applicable when grouping instance_data by effect.
	bool group_sprites_by_effect;           //Set to true to group all instance_data by effect. Effects can then be drawn in specific orders or not drawn at all on an effect by effect basis.
	bool auto_order_effects;                //When group_sprites_by_effect is true then you can set this to true to sort the effects each frame. Use tfx_SetPMCamera in 3d to set the effect depth to the distance the camera, in 2d the depth is set to the effect y position.
	bool is_3d;                             //All effects are 3d
	bool write_direct_to_staging_buffer;	//Make the particle manager write directly to the staging buffer. Use tfx_SetStagingBuffer before you call tfx_UpdateParticleManager
	void *user_data;						//User data that will get passed into the grow_staging_buffer_callback function which you can use to grow the buffer
	//If you need the staging buffer to be grown dynamically then you can use this call back to do that. It should return true if the buffer was successfully grown or false otherwise.
	bool(*grow_staging_buffer_callback)(tfxU32 new_size, tfx_particle_manager_t *pm, void *user_data);
}tfx_particle_manager_info_t;

//Use the particle manager to add multiple effects to your scene 
typedef struct tfx_particle_manager_s {
#ifdef __cplusplus
	tfx_vector_t<tfx_soa_buffer_t> particle_array_buffers;
	tfx_bucket_array_t<tfx_particle_soa_t> particle_arrays;
	tfx_vector_t<tfx_soa_buffer_t> particle_location_buffers;
	tfx_bucket_array_t<tfx_spawn_points_soa_t> particle_location_arrays;
	tfx_storage_map_t<tfx_vector_t<tfxU32>> free_particle_lists;
	tfx_storage_map_t<tfx_vector_t<tfxU32>> free_particle_location_lists;
	//Only used when using distance from camera ordering. New particles are put in this list and then merge sorted into the particles buffer
	tfx_vector_t<tfx_sort_work_entry_t> sorting_work_entry;
	tfx_vector_t<tfx_spawn_work_entry_t> spawn_work;
	tfx_vector_t<tfx_control_work_entry_t> control_work;
	tfx_vector_t<tfx_particle_age_work_entry_t> age_work;
	tfx_vector_t<tfxParticleID> particle_indexes;
	tfx_vector_t<tfxU32> free_particle_indexes;
	tfx_vector_t<tfx_effect_index_t> effects_in_use[tfxMAXDEPTH][2];
	tfx_vector_t<tfxU32> control_emitter_queue;
	tfx_vector_t<tfxU32> emitters_check_capture;
	tfx_vector_t<tfx_effect_index_t> free_effects;
	tfx_vector_t<tfxU32> free_emitters;
	tfx_vector_t<tfxU32> free_path_quaternions;
	tfx_vector_t<tfx_path_quaternion_t *> path_quaternions;
	tfx_vector_t<tfx_effect_state_t> effects;
	tfx_vector_t<tfx_emitter_state_t> emitters;
	tfx_vector_t<tfx_spawn_work_entry_t *> deffered_spawn_work;
	tfx_vector_t <tfx_unique_sprite_id_t> unique_sprite_ids[2][tfxLAYERS];
	tfx_vector_t<unsigned int> free_compute_controllers;
#else
	tfx_vector_t particle_array_buffers;
	tfx_bucket_array_t particle_arrays;
	tfx_vector_t particle_location_buffers;
	tfx_bucket_array_t particle_location_arrays;
	tfx_storage_map_t free_particle_lists;
	tfx_storage_map_t free_particle_location_lists;
	//Only used when using distance from camera ordering. New particles are put in this list and then merge sorted into the particles buffer
	tfx_vector_t sorting_work_entry;
	tfx_vector_t spawn_work;
	tfx_vector_t control_work;
	tfx_vector_t age_work;
	tfx_vector_t particle_indexes;
	tfx_vector_t free_particle_indexes;
	tfx_vector_t effects_in_use[tfxMAXDEPTH][2];
	tfx_vector_t control_emitter_queue;
	tfx_vector_t emitters_check_capture;
	tfx_vector_t free_effects;
	tfx_vector_t free_emitters;
	tfx_vector_t free_path_quaternions;
	tfx_vector_t path_quaternions;
	tfx_vector_t effects;
	tfx_vector_t emitters;
	tfx_vector_t deffered_spawn_work;
	tfx_vector_t unique_sprite_ids[2][tfxLAYERS];
	tfx_vector_t free_compute_controllers;
#endif

	tfx_library_t *library;
	tfx_work_queue_t work_queue;
	//The info config that was used to initialise the particle manager. This can be used to alter and the reconfigure the particle manager
	tfx_particle_manager_info_t info;
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
	//The maximum number of effects that can be updated per frame in the particle manager. If you're running effects with particles that have sub effects then this number might need 
	//to be relatively high depending on your needs. Use Init to udpate the sizes if you need to. Best to call Init at the start with the max numbers that you'll need for your application and don't adjust after.
	tfxU32 max_effects;
	//The maximum number of particles that can be updated per frame per layer. #define tfxLAYERS to set the number of allowed layers. This is currently 4 by default
	tfxU32 max_cpu_particles_per_layer[tfxLAYERS];
	//The maximum number of particles that can be updated per frame per layer in the compute shader. #define tfxLAYERS to set the number of allowed layers. This is currently 4 by default
	tfxU32 max_new_compute_particles;
	//The current effect buffer in use, can be either 0 or 1
	tfxU32 current_ebuff;
	tfxU32 next_ebuff;
	//For looping through active effects with GetNextEffect function
	tfxU32 effect_index_position;

	tfxU32 effects_start_size[tfxMAXDEPTH];

	tfxU32 sprite_index_point[tfxLAYERS];
	tfxU32 cumulative_index_point[tfxLAYERS];
	tfxU32 layer_sizes[tfxLAYERS];

	int mt_batch_size;
	tfx_sync_t particle_index_mutex;
	tfx_sync_t add_effect_mutex;

	tfx_random_t random;
	tfx_random_t threaded_random;
	tfxU32 max_compute_controllers;
	tfxU32 highest_compute_controller_index;
	tfx_compute_fx_global_state_t compute_global_state;
	tfxU32 sort_passes;
	tfx_lookup_mode lookup_mode;
	//For when particles are ordered by distance from camera (3d effects)
	tfx_vec3_t camera_front;
	tfx_vec3_t camera_position;

	tfxU32 unique_particle_id;    //Used when recording sprite data
	//When using single particles, you can flag the emitter to set the max_age of the particle to the 
	//length in time of the animation so that it maps nicely to the animation
	float animation_length_in_time;

	//These can possibly be removed at some point, they're debugging variables
	tfxU32 particle_id;
	tfxParticleManagerFlags flags;
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
} tfx_particle_manager_t;

typedef struct tfx_effect_library_stats_s {
	tfxU32 total_effects;
	tfxU32 total_sub_effects;
	tfxU32 total_emitters;
	tfxU32 total_attribute_nodes;
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

typedef struct tfx_library_s {
#ifdef __cplusplus
	tfx_storage_map_t<tfx_effect_emitter_t *> effect_paths;
	tfx_vector_t<tfx_effect_emitter_t> effects;
	tfx_storage_map_t<tfx_image_data_t> particle_shapes;
	tfx_vector_t<tfx_effect_emitter_info_t> effect_infos;
	tfx_vector_t<tfx_emitter_properties_t> emitter_properties;
	tfx_storage_map_t<tfx_sprite_data_t> pre_recorded_effects;

	tfx_bucket_array_t<tfx_emitter_path_t> paths;
	tfx_vector_t<tfx_global_attributes_t> global_graphs;
	tfx_vector_t<tfx_emitter_attributes_t> emitter_attributes;
	tfx_vector_t<tfx_transform_attributes_t> transform_attributes;
	tfx_vector_t<tfx_sprite_sheet_settings_t> sprite_sheet_settings;
	tfx_vector_t<tfx_sprite_data_settings_t> sprite_data_settings;
	tfx_vector_t<tfx_preview_camera_settings_t> preview_camera_settings;
	tfx_vector_t<tfx_attribute_node_t> all_nodes;
	tfx_vector_t<tfx_effect_lookup_data_t> node_lookup_indexes;
	tfx_vector_t<float> compiled_lookup_values;
	tfx_vector_t<tfx_graph_lookup_index_t> compiled_lookup_indexes;
	//This could probably be stored globally
	tfx_vector_t<tfx_vec4_t> graph_min_max;

	tfx_vector_t<tfxU32> free_global_graphs;
	tfx_vector_t<tfxU32> free_keyframe_graphs;
	tfx_vector_t<tfxU32> free_emitter_attributes;
	tfx_vector_t<tfxU32> free_animation_settings;
	tfx_vector_t<tfxU32> free_preview_camera_settings;
	tfx_vector_t<tfxU32> free_properties;
	tfx_vector_t<tfxU32> free_infos;
	tfx_vector_t<tfxU32> free_keyframes;
#else
	tfx_storage_map_t effect_paths;
	tfx_vector_t effects;
	tfx_storage_map_t particle_shapes;
	tfx_vector_t effect_infos;
	tfx_vector_t emitter_properties;
	tfx_storage_map_t pre_recorded_effects;

	tfx_bucket_array_t paths;
	tfx_vector_t global_graphs;
	tfx_vector_t emitter_attributes;
	tfx_vector_t transform_attributes;
	tfx_vector_t sprite_sheet_settings;
	tfx_vector_t sprite_data_settings;
	tfx_vector_t preview_camera_settings;
	tfx_vector_t all_nodes;
	tfx_vector_t node_lookup_indexes;
	tfx_vector_t compiled_lookup_values;
	tfx_vector_t compiled_lookup_indexes;
	tfx_vector_t graph_min_max;

	tfx_vector_t free_global_graphs;
	tfx_vector_t free_keyframe_graphs;
	tfx_vector_t free_emitter_attributes;
	tfx_vector_t free_animation_settings;
	tfx_vector_t free_preview_camera_settings;
	tfx_vector_t free_properties;
	tfx_vector_t free_infos;
	tfx_vector_t free_keyframes;
#endif

	tfx_gpu_shapes_t gpu_shapes;
	tfx_color_ramp_bitmap_data_t color_ramps;
	//Get an effect from the library by index
	tfx_str64_t name;
	bool open_library;
	bool dirty;
	tfx_stream_t library_file_path;
	tfxU32 uid;
	void(*uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset);
} tfx_library_t;

typedef struct tfx_effect_template_s {
#ifdef __cplusplus
	tfx_storage_map_t<tfx_effect_emitter_t *> paths;
#else
	tfx_storage_map_t paths;
#endif
	tfx_effect_emitter_t effect;
	tfxKey original_effect_hash;
}tfx_effect_template_t;

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
tfxINTERNAL inline tfxParticleID tfx__make_particle_id(tfxU32 bank_index, tfxU32 particle_index) { return ((bank_index & 0x00000FFF) << 20) + particle_index; }
tfxINTERNAL inline tfxU32 tfx__particle_index(tfxParticleID id) { return id & 0x000FFFFF; }
tfxINTERNAL inline tfxU32 tfx__particle_bank(tfxParticleID id) { return (id & 0xFFF00000) >> 20; }
tfxINTERNAL tfxU32 tfx__grab_particle_lists(tfx_particle_manager_t *pm, tfxKey emitter_hash, bool is_3d, tfxU32 reserve_amount, tfxEmitterControlProfileFlags flags);
tfxINTERNAL tfxU32 tfx__grab_particle_location_lists(tfx_particle_manager_t *pm, tfxKey emitter_hash, bool is_3d, tfxU32 reserve_amount);

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
tfxAPI_EDITOR void tfx__stream_emitter_properties(tfx_emitter_properties_t *property, tfxEmitterPropertyFlags flags, tfx_stream_t *file);
tfxAPI_EDITOR void tfx__stream_effect_properties(tfx_effect_emitter_t *effect, tfx_stream_t *file);
tfxAPI_EDITOR void tfx__stream_path_properties(tfx_effect_emitter_t *effect, tfx_stream_t *file);
tfxAPI_EDITOR void tfx__stream_graph(const char *name, tfx_graph_t *graph, tfx_stream_t *file);
tfxAPI_EDITOR void tfx__stream_graph_line(const char *name, tfx_graph_t *graph, tfx_str512_t *line);
tfxAPI_EDITOR void tfx__split_string_stack(const char *s, int length, tfx_vector_t<tfx_str256_t> *pair, char delim = 61);
tfxAPI_EDITOR bool tfx__string_is_uint(const char *s);
tfxAPI_EDITOR bool tfx__line_is_uint(tfx_line_t *line);
tfxAPI_EDITOR void tfx__assign_property_line(tfx_effect_emitter_t *effect, tfx_vector_t<tfx_str256_t> *pair, tfxU32 file_version);
tfxAPI_EDITOR void tfx__assign_effector_property_u32(tfx_effect_emitter_t *effect, tfx_str256_t *field, tfxU32 value, tfxU32 file_version);
tfxAPI_EDITOR void tfx__assign_effector_property(tfx_effect_emitter_t *effect, tfx_str256_t *field, float value);
tfxAPI_EDITOR void tfx__assign_effector_property_bool(tfx_effect_emitter_t *effect, tfx_str256_t *field, bool value);
tfxAPI_EDITOR void tfx__assign_effector_property_int(tfx_effect_emitter_t *effect, tfx_str256_t *field, int value);
tfxAPI_EDITOR void tfx__assign_effector_property_str(tfx_effect_emitter_t *effect, tfx_str256_t *field, const char *value);
tfxAPI_EDITOR void tfx__assign_graph_data(tfx_effect_emitter_t *effect, tfx_vector_t<tfx_str256_t> *values);
tfxAPI_EDITOR tfx_hsv_t tfx__rgb_to_hsv(tfx_rgb_t in);
tfxAPI_EDITOR tfx_rgb_t tfx__hsv_to_rgb(tfx_hsv_t in);
tfxAPI_EDITOR tfx_rgba8_t tfx__convert_float_color(float color_array[4]);
tfxAPI_EDITOR float tfx__length_vec3(tfx_vec3_t const *v);
tfxAPI_EDITOR float tfx__dot_product_vec2(const tfx_vec2_t *a, const tfx_vec2_t *b);
tfxAPI_EDITOR void tfx__catmull_rom_spline_2d_soa(const float *p_x, const float *p_y, int p0, float t, float vec[2]);
tfxAPI_EDITOR void tfx__catmull_rom_spline_3d_soa(const float *p_x, const float *p_y, const float *p_z, int p0, float t, float vec[3]);
tfxAPI_EDITOR void tfx__catmull_rom_spline_3d(const tfx_vec4_t *p0, const tfx_vec4_t *p1, const tfx_vec4_t *p2, const tfx_vec4_t *p3, float t, float vec[3]);
tfxAPI_EDITOR void tfx__catmull_rom_spline_gradient_3d(const tfx_vec4_t *p0, const tfx_vec4_t *p1, const tfx_vec4_t *p2, const tfx_vec4_t *p3, float t, float vec[3]);
tfxAPI_EDITOR float tfx__vec2_length_fast(tfx_vec2_t const *v);
tfxAPI_EDITOR float tfx__vec3_length_fast(tfx_vec3_t const *v);
tfxAPI_EDITOR void tfx__wide_transform_quaternion_vec3(const tfx_quaternion_t *q, tfxWideFloat *x, tfxWideFloat *y, tfxWideFloat *z);
tfxAPI_EDITOR void tfx__wide_transform_quaternion_vec2(const tfx_quaternion_t *q, tfxWideFloat *x, tfxWideFloat *y);
tfxAPI_EDITOR tfxU32 tfx__pack16bit_sscaled(float x, float y, float max_value);
tfxAPI_EDITOR float tfx__distance_2d(float fromx, float fromy, float tox, float toy);
tfxAPI_EDITOR void tfx__transform_3d(tfx_vec3_t *out_rotations, tfx_vec3_t *out_local_rotations, float *out_scale, tfx_vec3_t *out_position, tfx_vec3_t *out_local_position, tfx_vec3_t *out_translation, tfx_quaternion_t *out_q, const tfx_effect_state_t *parent);
tfxAPI_EDITOR void tfx__update_emitter_control_profile(tfx_effect_emitter_t *emitter);
tfxAPI_EDITOR void tfx__complete_particle_manager_work(tfx_particle_manager_t *pm);
tfxAPI_EDITOR tfx_mat3_t tfx__create_matrix3(float v = 1.f);
tfxAPI_EDITOR tfx_mat3_t tfx__rotate_matrix3(tfx_mat3_t const *m, float r);
tfxAPI_EDITOR void tfx__split_string_vec(const char *s, int length, tfx_vector_t<tfx_str256_t> *pair, char delim = 61);
tfxINTERNAL	tfx_line_t tfx__read_line(const char *s);
tfxINTERNAL void tfx__wide_transform_packed_quaternion_vec2(tfxWideInt *quaternion, tfxWideFloat *x, tfxWideFloat *y);
tfxINTERNAL void tfx__wide_transform_packed_quaternion_vec3(tfxWideInt *quaternion, tfxWideFloat *x, tfxWideFloat *y, tfxWideFloat *z);
tfxINTERNAL tfxU32 tfx__pack8bit_quaternion(tfx_quaternion_t v);
tfxINTERNAL tfxWideInt tfx__wide_pack16bit(tfxWideFloat v_x, tfxWideFloat v_y);
tfxINTERNAL tfxWideInt tfx__wide_pack16bit_2sscaled(tfxWideFloat v_x, tfxWideFloat  v_y, float max_value);
tfxINTERNAL tfxWideInt tfx__wide_pack8bit_xyz(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z);
tfxINTERNAL void tfx__wide_unpack8bit(tfxWideInt in, tfxWideFloat &x, tfxWideFloat &y, tfxWideFloat &z, tfxWideFloat &w);
tfxINTERNAL tfx_quaternion_t tfx__unpack8bit_quaternion(tfxU32 in);
tfxINTERNAL tfx_vec3_t tfx__get_emission_direciton_3d(tfx_particle_manager_t *pm, tfx_library_t *library, tfx_random_t *random, tfx_emitter_state_t &emitter, float emission_pitch, float emission_yaw, tfx_vec3_t local_position, tfx_vec3_t world_position);
tfxINTERNAL tfx_quaternion_t tfx__get_path_rotation_3d(tfx_random_t *random, float range, float pitch, float yaw, bool y_axis_only);
tfxINTERNAL tfx_quaternion_t tfx__get_path_rotation_2d(tfx_random_t *random, float range, float angle);
tfxINTERNAL tfx_vec3_t tfx__normalize_vec3_fast(tfx_vec3_t const *v);
tfxINTERNAL tfx_vec3_t tfx__cylinder_surface_normal(float x, float z, float width, float depth);
tfxINTERNAL tfx_vec3_t tfx__ellipse_surface_normal(float x, float y, float z, float width, float height, float depth);
tfxINTERNAL tfx_vec2_t tfx__catmull_rom_spline_gradient_2d_soa(const float *px, const float *py, float t);
tfxINTERNAL tfx_vec3_t tfx__catmull_rom_spline_gradient_3d_soa(const float *px, const float *py, const float *pz, float t);
tfxINTERNAL void tfx__wide_catmull_rom_spline_2d(tfxWideArrayi *pi, tfxWideFloat t, float *x, float *y, tfxWideFloat *vx, tfxWideFloat *vy);
tfxINTERNAL void tfx__wide_catmull_tom_spline_3d(tfxWideArrayi *pi, tfxWideFloat t, float *x, float *y, float *z, tfxWideFloat *vx, tfxWideFloat *vy, tfxWideFloat *vz);
tfxINTERNAL void tfx__catmull_rom_spline_2d(const tfx_vec4_t *p0, const tfx_vec4_t *p1, const tfx_vec4_t *p2, const tfx_vec4_t *p3, float t, float vec[2]);
tfxINTERNAL float tfx__length_vec3_nosqr(tfx_vec3_t const *v);
tfxINTERNAL void tfx__assign_effector_property_u64(tfx_effect_emitter_t *effect, tfx_str256_t *field, tfxU64 value, tfxU32 file_version);
tfxINTERNAL void tfx__add_data_value_double(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key, double value);
tfxINTERNAL void tfx__add_color_value_from_int(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key, tfxU32 value);
tfxINTERNAL int tfx__get_data_int_value(tfx_storage_map_t<tfx_data_entry_t> *config, const char *key);

//--------------------------------
//Graph functions
//Mainly used by the editor to edit graphs so these are kind of API functions but you wouldn't generally use these outside of the particle editor
//--------------------------------
tfxAPI_EDITOR void tfx__init_paths_soa_3d(tfx_soa_buffer_t *buffer, tfx_path_nodes_soa_t *soa, tfxU32 reserve_amount);
tfxAPI_EDITOR void tfx__init_emitter_properties(tfx_emitter_properties_t *properties);
tfxAPI_EDITOR tfx_attribute_node_t *tfx__add_graph_node_values(tfx_graph_t *graph, float frame, float value, tfxAttributeNodeFlags flags = 0, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
tfxAPI_EDITOR float tfx__get_graph_value_by_age(tfx_graph_t *graph, float age);
tfxAPI_EDITOR float tfx__get_graph_value_by_percent_of_life(tfx_graph_t *graph, float age, float life);
tfxAPI_EDITOR tfx_attribute_node_t *tfx__get_graph_last_node(tfx_graph_t *graph);
tfxAPI_EDITOR float tfx__get_graph_first_value(tfx_graph_t *graph);
tfxAPI_EDITOR tfx_attribute_node_t *tfx__insert_graph_node(tfx_graph_t *graph, float, float);
tfxAPI_EDITOR float *tfx__link_graph_first_value(tfx_graph_t *graph);
tfxAPI_EDITOR float *tfx__link_graph_last_value(tfx_graph_t *graph);
tfxAPI_EDITOR float tfx__get_graph_last_value(tfx_graph_t *graph);
tfxAPI_EDITOR float tfx__graph_value_by_index(tfx_graph_t *graph, tfxU32 index);
tfxAPI_EDITOR tfx_attribute_node_t *tfx__find_graph_node(tfx_graph_t *graph, tfx_attribute_node_t *n);
tfxAPI_EDITOR void tfx__validate_graph_curves(tfx_graph_t *graph);
tfxAPI_EDITOR void tfx__delete_graph_node(tfx_graph_t *graph, tfx_attribute_node_t *n);
tfxAPI_EDITOR void tfx__reset_graph(tfx_graph_t *graph, float first_node_value, tfx_graph_preset preset, bool add_node = true, float max_frames = 0);
tfxAPI_EDITOR void tfx__reset_graph_nodes(tfx_graph_t *graph, float first_node_value, tfx_graph_preset preset, bool add_node = true);
tfxAPI_EDITOR void tfx__clear_graph_to_one(tfx_graph_t *graph, float value);
tfxAPI_EDITOR void tfx__free_graph(tfx_graph_t *graph);
tfxAPI_EDITOR void tfx__copy_graph(tfx_graph_t *graph, tfx_graph_t *to, bool compile = true);
tfxAPI_EDITOR void tfx__copy_graph_color(tfx_overtime_attributes_t *from, tfx_overtime_attributes_t *to);
tfxAPI_EDITOR void tfx__copy_graph_color_hint(tfx_overtime_attributes_t *from, tfx_overtime_attributes_t *to);
tfxAPI_EDITOR void tfx__copy_graph_colors(tfx_graph_t *from_red, tfx_graph_t *from_blue, tfx_graph_t *from_green, tfx_graph_t *to_red, tfx_graph_t *to_green, tfx_graph_t *to_blue);
tfxAPI_EDITOR bool tfx__sort_graph(tfx_graph_t *graph);
tfxAPI_EDITOR void tfx__glip_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_blend_factor_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_overtime_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_factor_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_global_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_angle_graph(tfx_graph_t *graph);
tfxAPI_EDITOR bool tfx__is_translation_graph(tfx_graph_t *graph);
tfxAPI_EDITOR void tfx__multiply_all_graph_values(tfx_graph_t *graph, float scalar);
tfxAPI_EDITOR void tfx__copy_graph_no_lookups(tfx_graph_t *src_graph, tfx_graph_t *dst_graph);
tfxAPI_EDITOR void tfx__drag_graph_values(tfx_graph_preset preset, float *frame, float *value);
tfxAPI_EDITOR void tfx__compile_graph(tfx_graph_t *graph);
tfxAPI_EDITOR void tfx__compile_graph_overtime(tfx_graph_t *graph);
tfxAPI_EDITOR void tfx__compile_color_ramp(tfx_overtime_attributes_t *attributes, tfx_color_ramp_t *ramp, float gamma = tfxGAMMA);
tfxAPI_EDITOR void tfx__compile_color_ramp_hint(tfx_overtime_attributes_t *attributes, tfx_color_ramp_t *ramp, float gamma = tfxGAMMA);
tfxAPI_EDITOR void tfx__edit_color_ramp_bitmap(tfx_library_t *library, tfx_overtime_attributes_t *a, tfxU32 ramp_id);
tfxINTERNAL void tfx__init_paths_soa_2d(tfx_soa_buffer_t *buffer, tfx_path_nodes_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void tfx__add_graph_node(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxINTERNAL void tfx__set_graph_node(tfx_graph_t *graph, tfxU32 index, float frame, float value, tfxAttributeNodeFlags flags = 0, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
tfxINTERNAL float tfx__get_graph_random_value(tfx_graph_t *graph, float age, tfx_random_t *seed);
tfxINTERNAL tfx_attribute_node_t *tfx__get_graph_next_node(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxINTERNAL tfx_attribute_node_t *tfx__get_graph_prev_node(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxINTERNAL tfx_attribute_node_t *tfx__add_graph_coord_node(tfx_graph_t *graph, float, float);
tfxINTERNAL float tfx__get_graph_max_value(tfx_graph_t *graph);
tfxINTERNAL float tfx__get_graph_min_value(tfx_graph_t *graph);
tfxINTERNAL float tfx__get_graph_last_frame(tfx_graph_t *graph, float udpate_frequence);
tfxINTERNAL tfx_attribute_node_t *tfx__graph_node_by_index(tfx_graph_t *graph, tfxU32 index);
tfxINTERNAL float tfx__graph_frame_by_index(tfx_graph_t *graph, tfxU32 index);
tfxINTERNAL void tfx__delete_graph_node_at_frame(tfx_graph_t *graph, float frame);
tfxINTERNAL void tfx__clear_graph(tfx_graph_t *graph);
tfxINTERNAL void tfx__reindex_graph(tfx_graph_t *graph);
tfxINTERNAL bool tfx__color_graph(tfx_graph_t *graph);
tfxINTERNAL tfx_vec4_t tfx__get_min_max_graph_values(tfx_graph_preset preset);
tfxINTERNAL tfx_vec2_t tfx__get_quad_bezier_clamp(tfx_vec2_t p0, tfx_vec2_t p1, tfx_vec2_t p2, float t, float ymin, float ymax);
tfxINTERNAL tfx_vec2_t tfx__get_cubic_bezier_clamp(tfx_vec2_t p0, tfx_vec2_t p1, tfx_vec2_t p2, tfx_vec2_t p3, float t, float ymin, float ymax);
tfxINTERNAL tfx_vec3_t tfx__get_cubic_bezier_3d(tfx_vec4_t *p0, tfx_vec4_t *p1, tfx_vec4_t *p2, tfx_vec4_t *p3, float t);
tfxINTERNAL float tfx__get_bezier_value(const tfx_attribute_node_t *lastec, const tfx_attribute_node_t *a, float t, float ymin, float ymax);
tfxINTERNAL inline float tfx__get_vector_angle(float x, float y) { return atan2f(x, -y); }
tfxINTERNAL bool tfx__compare_nodes(tfx_attribute_node_t *left, tfx_attribute_node_t *right);
tfxINTERNAL void tfx__compile_graph_ramp_overtime(tfx_graph_t *graph);
tfxINTERNAL void tfx__compile_color_overtime(tfx_graph_t *graph, float gamma = tfxGAMMA);
tfxINTERNAL tfxKey tfx__hash_color_ramp(tfx_color_ramp_t *ramp);
tfxINTERNAL tfx_bitmap_t tfx__create_bitmap(int width, int height, int channels);
tfxINTERNAL void tfx__plot_bitmap(tfx_bitmap_t *image, int x, int y, tfx_rgba8_t color);
tfxINTERNAL void tfx__free_bitmap(tfx_bitmap_t *bitmap);
tfxINTERNAL void tfx__plot_color_ramp(tfx_bitmap_t *bitmap, tfx_color_ramp_t *ramp, tfxU32 y);
tfxINTERNAL void tfx__create_color_ramp_bitmaps(tfx_library_t *library);
tfxINTERNAL void tfx__maybe_insert_color_ramp_bitmap(tfx_library_t *library, tfx_overtime_attributes_t *a, tfxU32 ramp_id);
tfxINTERNAL tfxU32 tfx__add_color_ramp_to_bitmap(tfx_color_ramp_bitmap_data_t *ramp_data, tfx_color_ramp_t *ramp);
tfxINTERNAL void tfx__copy_color_ramp_to_animation_manager(tfx_animation_manager_t *animation_manager, tfxU32 properties_index, tfx_color_ramp_t *ramp);
tfxINTERNAL float tfx__get_max_life(tfx_effect_emitter_t *e);
tfxINTERNAL float tfx__lookup_fast_overtime(tfx_graph_t *graph, float age, float lifetime);
tfxINTERNAL float tfx__lookup_fast(tfx_graph_t *graph, float frame);
tfxINTERNAL float tfx__lookup_precise_overtime(tfx_graph_t *graph, float age, float lifetime);
tfxINTERNAL float tfx__lookup_precise(tfx_graph_t *graph, float frame);
tfxINTERNAL float tfx__get_random_fast(tfx_graph_t *graph, float frame, tfx_random_t *random);
tfxINTERNAL float tfx__get_random_precise(tfx_graph_t *graph, float frame, tfx_random_t *random);

//Node Manipulation
tfxAPI_EDITOR bool tfx__set_node(tfx_graph_t *graph, tfx_attribute_node_t *node, float *frame, float *value);
tfxAPI_EDITOR void tfx__set_node_curve(tfx_graph_t *graph, tfx_attribute_node_t *node, bool is_left_curve, float *frame, float *value);
tfxAPI_EDITOR bool tfx__move_node(tfx_graph_t *graph, tfx_attribute_node_t *node, float frame, float value, bool sort = true);
tfxAPI_EDITOR void tfx__clamp_graph_nodes(tfx_graph_t *graph);
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
tfxAPI_EDITOR bool tfx__has_more_than_one_key_frame(tfx_effect_emitter_t *e);
tfxAPI_EDITOR bool tfx__is_node_curve(tfx_attribute_node_t *node);
tfxAPI_EDITOR bool tfx__node_curves_are_initialised(tfx_attribute_node_t *node);
tfxAPI_EDITOR bool tfx__set_node_curve_initialised(tfx_attribute_node_t *node);
tfxINTERNAL void tfx__clamp_node(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxINTERNAL void tfx__clamp_node_curve(tfx_graph_t *graph, tfx_vec2_t *curve, tfx_attribute_node_t *node);
tfxINTERNAL bool tfx__is_color_graph_type(tfx_graph_type type);
tfxINTERNAL bool tfx__is_emitter_graph_type(tfx_graph_type type);
tfxINTERNAL bool tfx__has_key_frames(tfx_effect_emitter_t *e);
tfxINTERNAL void tfx__push_translation_points(tfx_effect_emitter_t *e, tfx_vector_t<tfx_vec3_t> *points, float frame);

//--------------------------------
//Grouped graph struct functions
//--------------------------------
tfxAPI_EDITOR void tfx__initialise_path(tfx_emitter_path_t *path);
tfxAPI_EDITOR void tfx__initialise_path_graphs(tfx_emitter_path_t *path, tfxU32 bucket_size = 8);
tfxAPI_EDITOR void tfx__reset_path_graphs(tfx_emitter_path_t *path, tfx_path_generator_type generator);
tfxAPI_EDITOR void tfx__build_path_nodes_3d(tfx_emitter_path_t *path);
tfxAPI_EDITOR void tfx__build_path_nodes_2d(tfx_emitter_path_t *path);
tfxAPI_EDITOR tfxU32 tfx__add_emitter_path_attributes(tfx_library_t *library);
tfxAPI_EDITOR void tfx__copy_path(tfx_emitter_path_t *src, const char *name, tfx_emitter_path_t *emitter_path);
tfxAPI_EDITOR bool tfx__has_translation_key_frames(tfx_transform_attributes_t *graphs);
tfxAPI_EDITOR void tfx__add_translation_nodes(tfx_transform_attributes_t *keyframes, float frame);
tfxAPI_EDITOR tfx_graph_t *tfx__get_graph(tfx_library_t *library, tfx_graph_id_t graph_id);
tfxAPI_EDITOR tfx_effect_library_stats_t tfx__create_library_stats(tfx_library_t *lib);
tfxAPI_EDITOR bool tfx__is_library_shape_used(tfx_library_t *library, tfxKey image_hash);
tfxAPI_EDITOR bool tfx__library_shape_exists(tfx_library_t *library, tfxKey image_hash);
tfxAPI_EDITOR bool tfx__remove_library_shape(tfx_library_t *library, tfxKey image_hash);
tfxAPI_EDITOR tfx_effect_emitter_t *tfx__insert_library_effect(tfx_library_t *library, tfx_effect_emitter_t *effect, tfx_effect_emitter_t *position);
tfxAPI_EDITOR tfx_effect_emitter_t *tfx__add_library_effect(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfx_effect_emitter_t *tfx__add_new_library_effect(tfx_library_t *library, tfx_str64_t *name);
tfxAPI_EDITOR tfx_effect_emitter_t *tfx__add_library_stage(tfx_library_t *library, tfx_str64_t *name);
tfxAPI_EDITOR void tfx__update_library_effect_paths(tfx_library_t *library);
tfxAPI_EDITOR bool tfx__rename_library_effect(tfx_library_t *library, tfx_effect_emitter_t *effect, const char *new_name);
tfxAPI_EDITOR bool tfx__library_name_exists(tfx_library_t *library, tfx_effect_emitter_t *effect, const char *name);
tfxAPI_EDITOR void tfx__reindex_library(tfx_library_t *library);
tfxAPI_EDITOR void tfx__update_library_particle_shape_references(tfx_library_t *library, tfxKey default_hash);
tfxAPI_EDITOR tfx_effect_emitter_t *tfx__library_move_up(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfx_effect_emitter_t *tfx__library_move_down(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void tfx__add_library_emitter_graphs(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void tfx__add_library_effect_graphs(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void tfx__add_library_transform_graphs(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfxU32 tfx__add_library_sprite_sheet_settings(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfxU32 tfx__add_library_sprite_data_settings(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void tfx__add_library_sprite_sheet_settings_sub(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void tfx__add_library_sprite_data_settings_sub(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfxU32 tfx__add_library_preview_camera_settings_effect(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void tfx__add_library_preview_camera_settings_sub_effects(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfxU32 tfx__allocate_library_preview_camera_settings(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 tfx__allocate_library_effect_emitter_info(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 tfx__allocate_library_emitter_properties(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 tfx__allocate_library_key_frames(tfx_library_t *library);
tfxAPI_EDITOR void tfx__update_library_compute_nodes(tfx_library_t *library);
tfxAPI_EDITOR void tfx__compile_all_library_graphs(tfx_library_t *library);
tfxAPI_EDITOR void tfx__compile_library_property_graphs(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void tfx__compile_library_base_graphs(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void tfx__compile_library_variation_graphs(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void tfx__compile_library_overtime_graph(tfx_library_t *library, tfxU32 index, bool including_color_ramps = true);
tfxAPI_EDITOR void tfx__compile_library_color_graphs(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void tfx__compile_library_graphs_of_effect(tfx_library_t *library, tfx_effect_emitter_t *effect, tfxU32 depth = 0);
tfxAPI_EDITOR void tfx__set_library_min_max_data(tfx_library_t *library);
tfxAPI_EDITOR void tfx__init_library(tfx_library_t *library);
tfxAPI_EDITOR bool tfx__is_valid_effect_path(tfx_library_t *library, const char *path);
tfxAPI_EDITOR bool tfx__is_valid_effect_key(tfx_library_t *library, tfxKey key);
tfxAPI_EDITOR tfx_effect_emitter_t *tfx__get_library_effect_by_key(tfx_library_t *library, tfxKey key);
tfxAPI_EDITOR void tfx__record_sprite_data(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect, float update_frequency, float camera_position[3], int *progress);
tfxINTERNAL void tfx__build_path_nodes_complex(tfx_emitter_path_t *path);
tfxINTERNAL void tfx__free_overtime_attributes(tfx_overtime_attributes_t *attributes);
tfxINTERNAL void tfx__copy_overtime_attributes_no_lookups(tfx_overtime_attributes_t *src, tfx_overtime_attributes_t *dst);
tfxINTERNAL void tfx__copy_overtime_attributes(tfx_overtime_attributes_t *src, tfx_overtime_attributes_t *dst);
tfxINTERNAL void tfx__free_factor_attributes(tfx_factor_attributes_t *attributes);
tfxINTERNAL void tfx__copy_factor_attributes_no_lookups(tfx_factor_attributes_t *src, tfx_factor_attributes_t *dst);
tfxINTERNAL void tfx__copy_factor_attributes(tfx_factor_attributes_t *src, tfx_factor_attributes_t *dst);
tfxINTERNAL void tfx__free_variation_attributes(tfx_variation_attributes_t *attributes);
tfxINTERNAL void tfx__copy_variation_attributes_no_lookups(tfx_variation_attributes_t *src, tfx_variation_attributes_t *dst);
tfxINTERNAL void tfx__copy_variation_attributes(tfx_variation_attributes_t *src, tfx_variation_attributes_t *dst);
tfxINTERNAL void tfx__free_base_attributes(tfx_base_attributes_t *attributes);
tfxINTERNAL void tfx__copy_base_attributes_no_lookups(tfx_base_attributes_t *src, tfx_base_attributes_t *dst);
tfxINTERNAL void tfx__copy_base_attributes(tfx_base_attributes_t *src, tfx_base_attributes_t *dst);
tfxINTERNAL void tfx__free_property_attributes(tfx_property_attributes_t *attributes);
tfxINTERNAL void tfx__copy_property_attributes_no_lookups(tfx_property_attributes_t *src, tfx_property_attributes_t *dst);
tfxINTERNAL void tfx__copy_property_attributes(tfx_property_attributes_t *src, tfx_property_attributes_t *dst);
tfxINTERNAL void tfx__free_transform_attributes(tfx_transform_attributes_t *attributes);
tfxINTERNAL void tfx__copy_transfrom_attributes_no_lookups(tfx_transform_attributes_t *src, tfx_transform_attributes_t *dst);
tfxINTERNAL void tfx__copy_transform_attributes(tfx_transform_attributes_t *src, tfx_transform_attributes_t *dst);
tfxINTERNAL void tfx__copy_global_attributes_no_lookups(tfx_global_attributes_t *src, tfx_global_attributes_t *dst);
tfxINTERNAL void tfx__copy_global_attributes(tfx_global_attributes_t *src, tfx_global_attributes_t *dst);
tfxINTERNAL int tfx__get_effect_library_stats(const char *filename, tfx_effect_library_stats_t *stats);
tfxINTERNAL void tfx__toggle_sprites_with_uid(tfx_particle_manager_t *pm, bool switch_on);
tfxINTERNAL tfxU32 tfx__get_library_lookup_indexes_size_in_bytes(tfx_library_t *library);
tfxINTERNAL tfxU32 tfx__get_library_lookup_values_size_in_bytes(tfx_library_t *library);
tfxINTERNAL void tfx__add_library_path(tfx_library_t *library, tfx_effect_emitter_t *effect_emitter, const char *path, bool skip_existing);
tfxINTERNAL tfxU32 tfx__add_library_global(tfx_library_t *library);
tfxINTERNAL tfxU32 tfx__add_library_emitter_attributes(tfx_library_t *library);
tfxINTERNAL void tfx__free_library_global(tfx_library_t *library, tfxU32 index);
tfxINTERNAL void tfx__free_library_key_frames(tfx_library_t *library, tfxU32 index);
tfxINTERNAL void tfx__free_library_emitter_attributes(tfx_library_t *library, tfxU32 index);
tfxINTERNAL void tfx__free_library_properties(tfx_library_t *library, tfxU32 index);
tfxINTERNAL void tfx_free_library_info(tfx_library_t *library, tfxU32 index);
tfxINTERNAL tfxU32 tfx__clone_library_global(tfx_library_t *library, tfxU32 source_index, tfx_library_t *destination_library);
tfxINTERNAL tfxU32 tfx__clone_library_key_frames(tfx_library_t *library, tfxU32 source_index, tfx_library_t *destination_library);
tfxINTERNAL tfxU32 tfx__clone_library_emitter_attributes(tfx_library_t *library, tfxU32 source_index, tfx_library_t *destination_library);
tfxINTERNAL tfxU32 tfx__clone_library_info(tfx_library_t *library, tfxU32 source_index, tfx_library_t *destination_library);
tfxINTERNAL tfxU32 tfx__clone_library_properties(tfx_library_t *library, tfx_emitter_properties_t *source, tfx_library_t *destination_library);
tfxINTERNAL void tfx__compile_library_global_graphs(tfx_library_t *library, tfxU32 index);
tfxINTERNAL void tfx__compile_library_key_frame_graphs(tfx_library_t *library, tfxU32 index);
tfxINTERNAL void tfx__compile_library_emitter_graphs(tfx_library_t *library, tfxU32 index);
tfxINTERNAL void tfx__compile_library_factor_graphs(tfx_library_t *library, tfxU32 index);
tfxINTERNAL tfx_str256_t tfx__find_new_path_name(tfx_library_t *library, const char *path);

//Effect/Emitter functions
tfxAPI_EDITOR void tfx__set_effect_user_data(tfx_effect_emitter_t *e, void *data);
tfxAPI_EDITOR void *tfx__get_effect_user_data(tfx_effect_emitter_t *e);
tfxAPI_EDITOR tfx_emitter_properties_t *tfx__get_effect_properties(tfx_effect_emitter_t *e);
tfxAPI_EDITOR tfx_effect_emitter_t *tfx__add_emitter_to_effect(tfx_effect_emitter_t *effect, tfx_effect_emitter_t *e);
tfxAPI_EDITOR tfx_effect_emitter_t *tfx__add_effect_to_emitter(tfx_effect_emitter_t *effect, tfx_effect_emitter_t *e);
tfxAPI_EDITOR int tfx__get_effect_depth(tfx_effect_emitter_t *e);
tfxAPI_EDITOR tfxU32 tfx__count_all_effects(tfx_effect_emitter_t *effect, tfxU32 amount);
tfxAPI_EDITOR tfx_effect_emitter_t *tfx__get_root_effect(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void tfx__reindex_effect(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfx_effect_emitter_t *tfx__move_effect_up(tfx_effect_emitter_t *effect_to_move);
tfxAPI_EDITOR tfx_effect_emitter_t *tfx__move_effect_down(tfx_effect_emitter_t *effect_to_move);
tfxAPI_EDITOR void tfx__delete_emitter_from_effect(tfx_effect_emitter_t *emitter_to_delete);
tfxAPI_EDITOR void tfx__clean_up_effect(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void tfx__reset_effect_graphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxAPI_EDITOR void tfx__reset_transform_graphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxAPI_EDITOR void tfx__reset_emitter_graphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxAPI_EDITOR void tfx__add_emitter_color_overtime(tfx_effect_emitter_t *effect, float frame, tfx_rgb_t color);
tfxAPI_EDITOR void tfx__add_emitter_color_hint_overtime(tfx_effect_emitter_t *effect, float frame, tfx_rgb_t color);
tfxAPI_EDITOR void tfx__update_effect_max_life(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfx_graph_t *tfx__get_effect_graph_by_type(tfx_effect_emitter_t *effect, tfx_graph_type type);
tfxAPI_EDITOR tfxU32 tfx__get_effect_graph_index_by_type(tfx_effect_emitter_t *effect, tfx_graph_type type);
tfxAPI_EDITOR tfx_emitter_attributes_t *tfx__get_emitter_attributes(tfx_effect_emitter_t *emitter);
tfxAPI_EDITOR void tfx__initialise_unitialised_graphs(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void tfx__set_effect_name(tfx_effect_emitter_t *effect, const char *n);
tfxAPI_EDITOR bool tfx__rename_sub_effector(tfx_effect_emitter_t *effect, const char *new_name);
tfxAPI_EDITOR bool tfx__effect_name_exists(tfx_effect_emitter_t *in_effect, tfx_effect_emitter_t *excluding_effect, const char *name);
tfxAPI_EDITOR void tfx__clone_effect(tfx_effect_emitter_t *effect_to_clone, tfx_effect_emitter_t *clone, tfx_effect_emitter_t *root_parent, tfx_library_t *destination_library, tfxEffectCloningFlags flags = 0);
tfxAPI_EDITOR void tfx__enable_all_emitters(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void tfx__disable_all_emitters(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void tfx__disable_all_emitters_except(tfx_effect_emitter_t *effect, tfx_effect_emitter_t *emitter);
tfxAPI_EDITOR bool tfx__is_finite_effect(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR bool tfx__is_finite_emitter(tfx_effect_emitter_t *emitter);
tfxAPI_EDITOR void tfx__flag_effect_as_3d(tfx_effect_emitter_t *effect, bool flag);
tfxAPI_EDITOR void tfx__flag_effects_as_3d(tfx_library_t *library);
tfxAPI_EDITOR bool tfx__is_3d_effect(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR bool tfx__is_ordered_effect(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfx_particle_manager_mode tfx__get_required_particle_manager_mode(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfx_preview_camera_settings_t *tfx__effect_camera_settings(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR float tfx__get_effect_highest_loop_length(tfx_effect_emitter_t *effect);
tfxINTERNAL void tfx__swap_effects(tfx_effect_emitter_t *left, tfx_effect_emitter_t *right);
tfxINTERNAL void tfx__swap_depth_index(tfx_depth_index_t *left, tfx_depth_index_t *right);
tfxINTERNAL tfx_effect_emitter_t *tfx__add_effect(tfx_effect_emitter_t *effect);
tfxINTERNAL void tfx__reset_emitter_base_graphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxINTERNAL void tfx__emitter_property_graphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxINTERNAL void tfx__reset_emitter_variation_graphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxINTERNAL void tfx__reset_emitter_overtime_graphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxINTERNAL void tfx__reset_emitter_factor_graphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxINTERNAL void tfx__enable_emitter(tfx_effect_emitter_t *effect);
tfxINTERNAL bool tfx__is_ordered_effect_state(tfx_effect_state_t *effect);
tfxINTERNAL void tfx__assign_stage_property_u32(tfx_effect_emitter_t *effect, tfx_str256_t *field, tfxU32 value);
tfxINTERNAL void tfx__assign_stage_property_float(tfx_effect_emitter_t *effect, tfx_str256_t *field, float value);
tfxINTERNAL void tfx__assign_stage_property_bool(tfx_effect_emitter_t *effect, tfx_str256_t *field, bool value);
tfxINTERNAL void tfx__assign_stage_property_int(tfx_effect_emitter_t *effect, tfx_str256_t *field, int value);
tfxINTERNAL void tfx__assign_stage_property_str(tfx_effect_emitter_t *effect, tfx_str256_t *field, tfx_str256_t *value);
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
tfxINTERNAL tfx128 tfx__dot128_xyz(const tfx128 *x1, const tfx128 *y1, const tfx128 *z1, const tfx128 *x2, const tfx128 *y2, const tfx128 *z);
tfxINTERNAL tfx128 tfx__dot128_xy(const tfx128 *x1, const tfx128 *y1, const tfx128 *x2, const tfx128 *y2);
tfxINTERNAL float tfx__length_vec4_nosqr(tfx_vec4_t const *v);
tfxINTERNAL float tfx__length_vec4(tfx_vec4_t const *v);
tfxINTERNAL float tfx__has_length_vec3(tfx_vec3_t const *v);
tfxINTERNAL tfx_vec3_t tfx__normalize_vec3(tfx_vec3_t const *v);
tfxINTERNAL tfx_vec4_t tfx__normalize_vec4(tfx_vec4_t const *v);
tfxINTERNAL tfx_vec3_t tfx__cross_product_vec3(tfx_vec3_t *a, tfx_vec3_t *b);
tfxINTERNAL void tfx__wide_cross_product(tfxWideFloat ax, tfxWideFloat ay, tfxWideFloat az, tfxWideFloat *bx, tfxWideFloat *by, tfxWideFloat *bz, tfxWideFloat *rx, tfxWideFloat *ry, tfxWideFloat *rz);
tfxINTERNAL float tfx__dot_product_vec4(const tfx_vec4_t *a, const tfx_vec4_t *b);
tfxINTERNAL float tfx__dot_product_vec3(const tfx_vec3_t *a, const tfx_vec3_t *b);
tfxINTERNAL float tfx__catmull_rom_segment(tfx_vector_t<tfx_vec4_t> *nodes, float length);
//Quake 3 inverse square root
tfxINTERNAL tfx_mat4_t tfx__create_matrix4(float v);
tfxINTERNAL tfx_mat4_t tfx__matrix4_rotate_x(float angle);
tfxINTERNAL tfx_mat4_t tfx__matrix4_rotate_y(float angle);
tfxINTERNAL tfx_mat4_t tfx__matrix4_rotate_z(float angle);
tfxINTERNAL tfx_mat4_t tfx__transform_matrix4(const tfx_mat4_t *in, const tfx_mat4_t *m);
tfxINTERNAL tfx_vec4_t tfx__transform_matrix4_vec4(const tfx_mat4_t *mat, const tfx_vec4_t vec);
tfxINTERNAL tfxU32 tfx__pack10bit_unsigned(tfx_vec3_t const *v);
tfxINTERNAL tfxWideInt tfx__wide_pack10bit_unsigned(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z);
tfxINTERNAL void tfx__wide_unpack10bit(tfxWideInt in, tfxWideFloat &x, tfxWideFloat &y, tfxWideFloat &z);
tfxINTERNAL tfxWideFloat tfx__wide_unpack10bit_y(tfxWideInt in);
tfxINTERNAL void tfx__transform_2d(tfx_vec3_t *out_rotations, tfx_vec3_t *out_local_rotations, float *out_scale, tfx_vec3_t *out_position, tfx_vec3_t *out_local_position, tfx_vec3_t *out_translation, tfx_quaternion_t *out_q, tfx_effect_state_t *parent);
tfxINTERNAL float tfx__gamma_correct(float color, float gamma = tfxGAMMA);
tfxINTERNAL tfx_vec2_t tfx__normalise_vec2(tfx_vec2_t const *v);
tfxINTERNAL tfx_vec2_t tfx__interpolate_vec2(float tween, tfx_vec2_t from, tfx_vec2_t to);
tfxINTERNAL tfx_vec3_t tfx__interpolate_vec3(float tween, tfx_vec3_t from, tfx_vec3_t to);
tfxINTERNAL tfx_rgba8_t tfx__interpolate_rgba8(float tween, tfx_rgba8_t from, tfx_rgba8_t to);
tfxINTERNAL tfxWideFloat tfx__wide_interpolate(tfxWideFloat tween, tfxWideFloat *from, tfxWideFloat *to);
tfxINTERNAL float tfx__interpolate_float(float tween, float from, float to);

//-------------------------------------------------
//--New transform_3d particle functions for SoA data--
//--------------------------2d---------------------
tfxINTERNAL void tfx__transform_particle_position(const float local_position_x, const float local_position_y, const float roll, tfx_vec2_t *world_position, float *world_rotations);

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
//Particle manager internal functions
//--------------------------------
template<typename T>
tfxINTERNAL inline void tfx__write_particle_color_sprite_data(T *sprites, tfxU32 start_diff, tfxU32 limit_index, const tfxU32 *depth_index, tfxU32 index, const tfxWideArrayi &packed_intensity_life, const tfxWideArrayi &curved_alpha, tfxU32 &running_sprite_index) {
	for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
		sprites[running_sprite_index].intensity_life.packed = packed_intensity_life.a[j];
		sprites[running_sprite_index].curved_alpha.packed = curved_alpha.a[j];
		running_sprite_index++;
	}
}

template<typename T>
tfxINTERNAL inline void tfx__write_particle_color_sprite_data_ordered(T *sprites, tfx_particle_manager_t &pm, tfxU32 layer, tfxU32 start_diff, tfxU32 limit_index, const tfxU32 *depth_index, tfxU32 index, const tfxWideArrayi &packed_intensity_life, const tfxWideArrayi &curved_alpha, tfxU32 &running_sprite_index, tfxU32 instance_offset) {
	for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
		tfxU32 sprite_depth_index = depth_index[index + j] + instance_offset;
		sprites[sprite_depth_index].intensity_life.packed = packed_intensity_life.a[j];
		sprites[sprite_depth_index].curved_alpha.packed = curved_alpha.a[j];
		running_sprite_index++;
	}
}

template<typename T>
tfxINTERNAL inline void tfx__write_particle_image_sprite_data(T *sprites, tfx_particle_manager_t &pm, tfxU32 layer, tfxU32 start_diff, tfxU32 limit_index, tfx_particle_soa_t &bank, tfxWideArrayi &flags, tfxWideArrayi &image_indexes, const tfxEmitterStateFlags emitter_flags, const tfx_billboarding_option billboard_option, tfxU32 index, tfxU32 &running_sprite_index) {
	for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
		int index_j = index + j;
		tfxU32 &sprites_index = bank.sprite_index[index_j];
		tfxU32 capture = flags.a[j];
		sprites[running_sprite_index].captured_index = capture == 0 ? (pm.current_sprite_buffer << 30) + running_sprite_index : (!pm.current_sprite_buffer << 30) + (sprites_index & 0x0FFFFFFF);
		sprites[running_sprite_index].captured_index |= emitter_flags & tfxEmitterStateFlags_wrap_single_sprite ? 0x80000000 : 0;
		sprites_index = layer + running_sprite_index;
		sprites[running_sprite_index].indexes = image_indexes.a[j];
		sprites[running_sprite_index].indexes |= (billboard_option << 13) | capture;
		bank.flags[index_j] &= ~tfxParticleFlags_capture_after_transform;
		running_sprite_index++;
	}
}

template<typename T>
tfxINTERNAL inline void tfx__write_particle_image_sprite_data_ordered(T *sprites, tfx_particle_manager_t &pm, tfxU32 layer, tfxU32 start_diff, tfxU32 limit_index, tfx_particle_soa_t &bank, tfxWideArrayi &flags, tfxWideArrayi &image_indexes, const tfxEmitterStateFlags emitter_flags, const tfx_billboarding_option billboard_option, tfxU32 index, tfxU32 &running_sprite_index, tfxU32 instance_offset) {
	for (tfxU32 j = start_diff; j < tfxMin(limit_index + start_diff, tfxDataWidth); ++j) {
		int index_j = index + j;
		tfxU32 sprite_depth_index = bank.depth_index[index_j] + instance_offset;
		tfxU32 &sprites_index = bank.sprite_index[index_j];
		tfxU32 capture = flags.a[j];
		sprites[sprite_depth_index].captured_index = capture == 0 && bank.single_loop_count[index_j] == 0 ? (pm.current_sprite_buffer << 30) + sprite_depth_index : (!pm.current_sprite_buffer << 30) + (sprites_index & 0x0FFFFFFF);
		sprites[sprite_depth_index].captured_index |= emitter_flags & tfxEmitterStateFlags_wrap_single_sprite ? 0x80000000 : 0;
		sprites_index = layer + sprite_depth_index;
		sprites[sprite_depth_index].indexes = image_indexes.a[j];
		sprites[sprite_depth_index].indexes |= (billboard_option << 13) | capture;
		bank.flags[index_j] &= ~tfxParticleFlags_capture_after_transform;
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

template<typename T>
tfxINTERNAL inline void tfx__invalidate_new_captured_index(T* instance, tfx_vector_t<tfx_unique_sprite_id_t> &uids, tfx_particle_manager_t *pm, tfxU32 layer) {
	for (tfxU32 i = 0; i != pm->instance_buffer_for_recording[pm->current_sprite_buffer][layer].current_size; ++i) {
		if ((uids[i].age == 0 && !(instance[i].captured_index & 0x80000000)) || (instance[i].captured_index & 0xC0000000) >> 30 == pm->current_sprite_buffer && !(instance[i].captured_index & 0x80000000)) {
			instance[i].captured_index = tfxINVALID;
		}
	}
}

template<typename T>
tfxINTERNAL inline void tfx__invalidate_offsetted_sprite_captured_index(T* instance, tfx_vector_t<tfx_unique_sprite_id_t> &uids, tfx_particle_manager_t *pm, tfxU32 layer) {
	for (tfxU32 i = 0; i != pm->instance_buffer_for_recording[pm->current_sprite_buffer][layer].current_size; ++i) {
		instance[i].captured_index = tfxINVALID;
	}
}

tfxINTERNAL float tfx__get_emission_direciton_2d(tfx_particle_manager_t *pm, tfx_library_t *library, tfx_random_t *random, tfx_emitter_state_t &emitter, tfx_vec2_t local_position, tfx_vec2_t world_position);
tfxINTERNAL tfx_vec3_t tfx__random_vector_in_cone(tfx_random_t *random, tfx_vec3_t cone_direction, float cone_angle);
tfxINTERNAL void tfx__wide_random_vector_in_cone(tfxWideInt seed, tfxWideFloat dx, tfxWideFloat dy, tfxWideFloat dz, tfxWideFloat cone_angle, tfxWideFloat *result_x, tfxWideFloat *result_y, tfxWideFloat *result_z);
tfxINTERNAL void tfx__transform_effector_2d(tfx_vec3_t *world_rotations, tfx_vec3_t *local_rotations, tfx_vec3_t *world_position, tfx_vec3_t *local_position, tfx_quaternion_t *q, tfx_sprite_transform2d_t *parent, bool relative_position = true, bool relative_angle = false);
tfxINTERNAL void tfx__transform_effector_3d(tfx_vec3_t *world_rotations, tfx_vec3_t *local_rotations, tfx_vec3_t *world_position, tfx_vec3_t *local_position, tfx_quaternion_t *q, tfx_sprite_transform3d_t *parent, bool relative_position = true, bool relative_angle = false);
tfxINTERNAL void tfx__update_effect(tfx_particle_manager_t *pm, tfxU32 index, tfxU32 parent_index = tfxINVALID);
tfxINTERNAL void tfx__update_emitter(tfx_work_queue_t *work_queue, void *data);
tfxINTERNAL tfxU32 tfx__new_sprites_needed(tfx_particle_manager_t *pm, tfx_random_t *random, tfxU32 index, tfx_effect_state_t *parent, tfx_emitter_properties_t *properties);
tfxINTERNAL void tfx__update_emitter_state(tfx_particle_manager_t *pm, tfx_emitter_state_t &emitter, tfxU32 parent_index, const tfx_parent_spawn_controls_t *parent_spawn_controls, tfx_spawn_work_entry_t *entry);
tfxINTERNAL void tfx__update_effect_state(tfx_particle_manager_t *pm, tfxU32 index);

tfxINTERNAL tfxU32 tfx__spawn_particles(tfx_particle_manager_t *pm, tfx_spawn_work_entry_t *work_entry);
tfxINTERNAL void tfx__spawn_particle_point_2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_line_2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_area_2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_ellipse_2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_path_2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_micro_update_2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_noise(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_motion_randomness(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_weight(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_velocity(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_roll(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_image_frame(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_age(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_size_2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_spin_2d(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void tfx__do_spawn_work_2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__do_spawn_work_3d(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void tfx__spawn_particle_point_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_other_emitter_2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_other_emitter_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_line_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_area_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_ellipsoid(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_cylinder(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_icosphere_random(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_icosphere(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_path_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_micro_update_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_spin_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__spawn_particle_size_3d(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void tfx__control_particles(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_age(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_image_frame(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_color(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_size(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_spin_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_spin(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_uid(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_position_2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_transform_2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_position_path_2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_position_basic_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_position_orbital_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_noise_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_orbital_noise_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_motion_randomness_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_motion_randomness_orbital_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_line_behaviour_kill(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_line_behaviour_loop(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_position_path_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_transform_3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__control_particle_bounding_box(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__init_sprite_data_soa_compression_3d(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void tfx__init_sprite_data_soa_3d(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void tfx__init_sprite_data_soa_compression_2d(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void tfx__init_sprite_data_soa_2d(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void tfx__init_particle_soa_2d(tfx_soa_buffer_t *buffer, tfx_particle_soa_t *soa, tfxU32 reserve_amount, tfxEmitterControlProfileFlags control_profile);
tfxINTERNAL void tfx__init_particle_soa_3d(tfx_soa_buffer_t *buffer, tfx_particle_soa_t *soa, tfxU32 reserve_amount, tfxEmitterControlProfileFlags control_profile);
tfxINTERNAL void tfx__init_particle_location_soa_3d(tfx_soa_buffer_t *buffer, tfx_spawn_points_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void tfx__init_particle_location_soa_2d(tfx_soa_buffer_t *buffer, tfx_spawn_points_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void tfx__copy_emitter_properties(tfx_emitter_properties_t *from_properties, tfx_emitter_properties_t *to_properties);
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

tfxAPI_EDITOR void tfx__free_path_graphs(tfx_emitter_path_t *path);
tfxINTERNAL void tfx__copy_path_graphs(tfx_emitter_path_t *src, tfx_emitter_path_t *dst);
tfxAPI_EDITOR tfxU32 tfx__create_emitter_path_attributes(tfx_effect_emitter_t *emitter, bool add_node);
tfxINTERNAL void tfx__initialise_global_attributes(tfx_global_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void tfx__initialise_overtime_attributes(tfx_overtime_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void tfx__initialise_factor_attributes(tfx_factor_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void tfx__initialise_variation_attributes(tfx_variation_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void tfx__initialise_base_attributes(tfx_base_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void tfx__initialise_property_attributes(tfx_property_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void tfx__initialise_transform_attributes(tfx_transform_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void tfx__initialise_emitter_attributes(tfx_emitter_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void tfx__free_emitter_attributes(tfx_emitter_attributes_t *attributes);
tfxINTERNAL void tfx__free_global_attributes(tfx_global_attributes_t *attributes);
tfxINTERNAL tfxErrorFlags tfx__load_effect_library_package(tfx_package package, tfx_library_t *lib, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset), void *user_data = nullptr);
tfxINTERNAL void tfx__build_gpu_shape_data(tfx_vector_t<tfx_image_data_t> *particle_shapes, tfx_gpu_shapes_t *shape_data, void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset));

//--------------------------------
//Animation manager internal functions - animation manager is used to playback pre-recorded effects
//--------------------------------
tfxINTERNAL tfxAnimationID tfx__allocate_animation_instance(tfx_animation_manager_t *animation_manager);
tfxINTERNAL void tfx__free_animation_instance(tfx_animation_manager_t *animation_manager, tfxU32 index);
tfxINTERNAL void tfx__add_effect_emitter_properties(tfx_animation_manager_t *animation_manager, tfx_effect_emitter_t *effect, bool *has_animated_shape);
tfxINTERNAL bool tfx__free_pm_effect_capacity(tfx_particle_manager_t *pm);
tfxINTERNAL void tfx__initialise_animation_manager(tfx_animation_manager_t *animation_manager, tfxU32 max_instances);

//--------------------------------
//Particle manager internal functions
//--------------------------------
tfxINTERNAL tfx_effect_index_t tfx__get_effect_slot(tfx_particle_manager_t *pm);
tfxINTERNAL tfxU32 tfx__get_emitter_slot(tfx_particle_manager_t *pm);
tfxINTERNAL tfxU32 tfx__get_particle_index_slot(tfx_particle_manager_t *pm, tfxParticleID particle_id);
tfxINTERNAL tfxU32 tfx__allocate_path_quaternion(tfx_particle_manager_t *pm, tfxU32 amount);
tfxINTERNAL void tfx__free_path_quaternion(tfx_particle_manager_t *pm, tfxU32 index);
tfxINTERNAL void tfx__free_particle_index(tfx_particle_manager_t *pm, tfxU32 *index);
tfxINTERNAL tfxU32 tfx__push_depth_index(tfx_vector_t<tfx_depth_index_t> *depth_indexes, tfx_depth_index_t depth_index);
tfxINTERNAL void tfx__reset_particle_manager_flags(tfx_particle_manager_t *pm);
tfxINTERNAL tfxU32 tfx__get_particle_sprite_index(tfx_particle_manager_t *pm, tfxParticleID id);
tfxINTERNAL void tfx__free_compute_slot(tfx_particle_manager_t *pm, unsigned int slot_id);
tfxINTERNAL tfxEffectID tfx__add_effect_to_particle_manager(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect, int buffer, int hierarchy_depth, bool is_sub_emitter, tfxU32 root_effect_index, float add_delayed_spawning);
tfxINTERNAL void tfx__free_particle_list(tfx_particle_manager_t *pm, tfxU32 index);
tfxINTERNAL void tfx__free_spawn_location_list(tfx_particle_manager_t *pm, tfxU32 index);
tfxINTERNAL void tfx__free_all_particle_lists(tfx_particle_manager_t *pm);
tfxINTERNAL void tfx__free_all_spawn_location_lists(tfx_particle_manager_t *pm);
tfxINTERNAL void tfx__order_effect_sprites(tfx_effect_instance_data_t *sprites, tfxU32 layer, tfx_particle_manager_t *pm);

//Compute stuff doesn't work currently. Keeping this here for now for when I get back to implementing compute shaders for TimelineFX
tfxINTERNAL void tfx__enable_compute(tfx_particle_manager_t *pm) { pm->flags |= tfxParticleManagerFlags_use_compute_shader; }
tfxINTERNAL void tfx__disable_compute(tfx_particle_manager_t *pm) { pm->flags &= ~tfxParticleManagerFlags_use_compute_shader; }
tfxINTERNAL int tfx__add_compute_controller(tfx_particle_manager_t *pm);
tfxINTERNAL tfx_compute_particle_t *tfx__grab_compute_particle(tfx_particle_manager_t *pm, unsigned int layer);
tfxINTERNAL void tfx__reset_particle_ptr(tfx_particle_manager_t *pm, void *ptr);
tfxINTERNAL void tfx__reset_controller_ptr(tfx_particle_manager_t *pm, void *ptr);
tfxINTERNAL void tfx__update_compute(tfx_particle_manager_t *pm, void *sampled_particles, unsigned int sample_size = 100);
tfxINTERNAL void tfx__init_common_particle_manager(tfx_particle_manager_t *pm, tfx_library_t *library, tfxU32 max_particles, unsigned int effects_limit, tfx_particle_manager_mode mode, bool double_buffered_sprites, bool dynamic_sprite_allocation, bool group_sprites_by_effect, tfxU32 mt_batch_size);
tfxINTERNAL bool tfx__valid_effect_id(tfx_particle_manager_t *pm, tfxEffectID id);
tfxINTERNAL tfxU32 tfx__count_library_global_lookup_values(tfx_library_t *library, tfxU32 index);
tfxINTERNAL tfxU32 tfx__count_library_emitter_lookup_values(tfx_library_t *library, tfxU32 index);

//--------------------------------
//Effect templates
//--------------------------------
tfxINTERNAL void tfx__add_template_path(tfx_effect_template_t *effect_template, tfx_effect_emitter_t *effect_emitter, const char *path);

//--------------------------------
//Library functions, internal/Editor functions
//--------------------------------
tfxINTERNAL void tfx__prepare_library_effect_template_path(tfx_library_t *library, const char *path, tfx_effect_template_t *effect);
tfxINTERNAL void tfx__reset_sprite_data_lerp_offset(tfx_sprite_data_t *sprites);
tfxINTERNAL void tfx__compress_sprite_data(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect, bool is_3d, float frame_length, int *progress);
tfxINTERNAL void tfx__link_up_sprite_captured_indexes(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void tfx__build_all_library_paths(tfx_library_t *library);
tfxINTERNAL tfx_str64_t tfx__get_name_from_path(const char *path);
tfxINTERNAL bool tfx__is_root_effect(tfx_effect_emitter_t *effect);
tfxINTERNAL void tfx__reset_effect_parents(tfx_effect_emitter_t *effect);
tfxINTERNAL void tfx__free_effect_graphs(tfx_effect_emitter_t *effect);
tfxINTERNAL tfxU32 tfx__count_all_effect_lookup_values(tfx_effect_emitter_t *effect);
tfxINTERNAL float tfx__get_effect_loop_length(tfx_effect_emitter_t *effect);

//------------------------------------------------------------
//Section API_Functions
//------------------------------------------------------------

#endif		//__cpluscplus

tfxAPI void tfx_UpdateAnimationManagerBufferMetrics(tfx_animation_manager_t *animation_manager);
tfxAPI tfx_storage_t *tfx_GetGlobals();
tfxAPI tfx_pool_stats_t tfx_CreateMemorySnapshot(tfx_header *first_block);
tfxAPI float tfx_DegreesToRadians(float degrees);
tfxAPI float tfx_RadiansToDegrees(float radians);
tfxAPI tfx_effect_emitter_info_t *tfx_GetEffectInfo(tfx_effect_emitter_t *e);

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
* @param filename        A pointer to a null-terminated string that contains the path and filename of the effect library package to be loaded.
* @param lib            A reference to a tfx_library_t object that will hold the loaded effect library data.
* @param shape_loader    A pointer to a function that will be used to load image data into the effect library package.
*                        The function has the following signature: void shape_loader(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data).
* @param user_data        A pointer to user-defined data that will be passed to the shape_loader function. This parameter is optional and can be set to nullptr if not needed.
* @param read_only        A boolean value that determines whether the effect library data will be loaded in read-only mode. (Maybe removed in the future).
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
tfxAPI tfxErrorFlags tfx_LoadEffectLibrary(const char *filename, tfx_library_t *lib, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset), void *user_data);

/**
* Loads an effect library package from memory into the provided tfx_library_t object pointer.
*
* @param data            A pointer to a memory buffer containing the library to be loaded
* @param size            The size of the memory buffer containing the library to be loaded
* @param lib            A reference to a tfx_library_t object that will hold the loaded effect library data.
* @param shape_loader    A pointer to a function that will be used to load image data into the effect library package.
*                        The function has the following signature: void shape_loader(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data).
* @param user_data        A pointer to user-defined data that will be passed to the shape_loader function. This parameter is optional and can be set to nullptr if not needed.
* @param read_only        A boolean value that determines whether the effect library data will be loaded in read-only mode. (Maybe removed in the future).
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
tfxAPI tfxErrorFlags tfx_LoadEffectLibraryFromMemory(const void *data, tfxU32 size, tfx_library_t *lib, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset), void *user_data);

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
tfxAPI tfxErrorFlags tfx_LoadSpriteData(const char *filename, tfx_animation_manager_t *animation_manager, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void *user_data);

/*
* Updates all the image data in the library using the uv_lookup that you set when loading a library. This allows you to add all of the uv data for
* the shapes that are loaded into the texture.
* @param tfx_library_t                A valid pointer to a tfx_library_t
*/
tfxAPI void tfx_UpdateLibraryGPUImageData(tfx_library_t *library);

/*
Output all the effect names in a library to the console
* @param tfx_library_t                A valid pointer to a tfx_library_t
*/
tfxAPI void ListEffectNames(tfx_library_t *library);

/*
Get an effect in the library by it's index. If you need to get an effect in a folder or an emitter then you can use tfx_GetLibraryEffectPath instead.
* @param tfx_library_t                A valid pointer to a tfx_library_t
*/
tfxAPI tfx_effect_emitter_t *tfx_GetEffectByIndex(tfx_library_t *library, int index);

/*
Get an effect in the library by it's path. So for example, if you want to get a pointer to the emitter "spark" in effect "explosion" then you could do GetEffect("explosion/spark")
You will need this function to apply user data and update callbacks to effects and emitters before adding the effect to the particle manager
* @param tfx_library_t                A valid pointer to a tfx_library_t
* @param const char *path             Path to the effect or emitter
*/
tfxAPI tfx_effect_emitter_t *tfx_GetLibraryEffectPath(tfx_library_t *library, const char *path);

/*
Free all the memory used by a library
* param tfx_library_t				A pointer to the library that you want to free
*/
tfxAPI void tfx_FreeLibrary(tfx_library_t *library);

/*
Create the image data required for shaders from a TimelineFX library. The image data will contain data such as uv coordinates. Once you have built the data you can use GetLibraryImageData to get the buffer
and upload it to the gpu.
* @param library                  A pointer to a tfx_library_t object
* @param shapes                   A pointer to a tfx_gpu_shapes_t object which will fill a buffer with all the shapes
* @param uv_lookup                A function pointer to a function that you need to set up in order to get the uv coordinates from whatever renderer you're using
*/
tfxAPI void tfx_BuildLibraryGPUShapeData(tfx_library_t *library, tfx_gpu_shapes_t *shapes, void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset));

/*
Get a pointer to the particle shapes data in a library. This can be used with tfx_BuildGPUShapeData when you want to upload the data to the GPU
* @param animation_manager        A pointer the tfx_animation_manager_t
*/
tfxAPI inline tfx_image_data_t *tfx_GetParticleShapesLibrary(tfx_library_t *library, int *count) {
	*count = library->particle_shapes.data.current_size;
	return library->particle_shapes.data.data;
}

//--------------------------------
//Particle_Manager_functions
//--------------------------------
/*
Create a tfx_particle_manager_info_t object which contains configuration data that you can pass to tfx_InitializeParticleManager to setup a particle manager. You can tweak the config after calling this
function if needed to fine tune the settings.
* @param setup                    A tfx_particle_manager_setup enum which you can use to set the info based on some commonly used templates
*/
tfxAPI tfx_particle_manager_info_t tfx_CreateParticleManagerInfo(tfx_particle_manager_setup setup);

/*
Initialize a particle manager with a tfx_particle_manager_info_t object which contains setup data for how to configure the particle manager. See tfx_CreateParticleManagerInfo
* @param pm						A pointer to an unitialised tfx_particle_manager_t. If you want to reconfigure a particle manager for a different usage then you can call tfx_ReconfigureParticleManager.
* @param library                A pointer to a tfx_library_t that you will be using to add all of the effects from to the particle manager.
* @param info                   A tfx_particle_manager_info_t pointer containing the configuration for the particle manager.
*/
tfxAPI void tfx_InitializeParticleManager(tfx_particle_manager_t *pm, tfx_library_t *library, tfx_particle_manager_info_t info);

/*
Reconfigure a particle manager to make it work in a different mode. A particle manager can only run in a single mode at time like unordered, depth ordered etc so use this to change that. Also bear
in mind that you can just use more than one particle manager and utilised different modes that way as well. The modes that you need will depend on the effects that you're adding to the particle manager.
* @param pm                       A pointer to an intialised tfx_particle_manager_t.
* @param mode                     One of the following modes:
								  tfxParticleManagerMode_unordered
								  tfxParticleManagerMode_ordered_by_age
								  tfxParticleManagerMode_ordered_by_depth
								  tfxParticleManagerMode_ordered_by_depth_guaranteed
* @param sort_passes              The number of sort passes if you're using depth sorted effects
* @param is_3d                    True if the particle manager should be configured for 3d effects.
*/
tfxAPI void tfx_ReconfigureParticleManager(tfx_particle_manager_t *pm, tfx_particle_manager_mode mode, tfxU32 sort_passes, bool is_3d);

/*
Set the staging buffer used in the particle manager. The particle manager flags must be set with tfxParticleManagerFlags_direct_to_staging_buffer when the particle
manager was created. Depending on the renderer you use you may have to call this before each time you update the particle manager so you can set the buffer to the
current frame in flight. This will probably apply in any modern renderer like vulkan, metal or dx12.
Note: It's up to you to ensure that the staging buffer has enough capacity. The particle manager will assume that the size_in_bytes that you pass to the particle
manager is correct and if tfxParticleManagerFlags_dynamic_sprite_allocation is set will attempt to grow the buffer by calling the callback you set to do this.
* @param pm                       A pointer to an intialised tfx_particle_manager_t.
* @param staging_buffer           A pointer to the staging buffer where all the instance_data/billboards are written to
* @param size_in_bytes            The size in bytes of the staging buffer
*/
tfxAPI void tfx_SetStagingBuffer(tfx_particle_manager_t *pm, void *staging_buffer, tfxU32 size_in_bytes);

/*
Turn on and off whether the particle manager should sort the effects by depth order. Use tfx_SetPMCamera to set the position of the camera that the particle manager will
use to update the depth of each effect in the scene (3d mode). In 2d mode the depth will be auto set to the y position of the effect.
* @param pm                       A pointer to an intialised tfx_particle_manager_t.
* @param yesno                    A boolean, set to true or false if you want auto ordering on or off respectively
*/
tfxAPI inline void tfx_TogglePMOrderEffects(tfx_particle_manager_t *pm, bool yesno) {
	if (yesno) {
		pm->flags |= tfxParticleManagerFlags_auto_order_effects;
	}
	else {
		pm->flags &= ~tfxParticleManagerFlags_auto_order_effects;
	}
}

/*
Get the sprite buffer in the particle manager containing all the 2d instance_data that were created the last frame. You can use this to copy to a staging buffer to upload to the gpu.
This will be a pointer to the start of the buffer for uploading all the instance_data. If you want to do this for each effect then you can call tfx_GetEffect2dInstanceBuffer.
* @param pm                       A pointer to an intialised tfx_particle_manager_t.
*/
tfxAPI inline tfx_2d_instance_t *tfx_Get2dInstanceBuffer(tfx_particle_manager_t *pm) {
	return tfxCastBufferRef(tfx_2d_instance_t, pm->instance_buffer);
}

/*
Get the billboard buffer in the particle manager containing all the 3d billboards that were created the last frame. You can use this to copy to a staging buffer to upload to the gpu.
* @param pm                       A pointer to an intialised tfx_particle_manager_t.
*/
tfxAPI inline tfx_3d_instance_t *tfx_Get3dInstanceBuffer(tfx_particle_manager_t *pm) {
	return tfxCastBufferRef(tfx_3d_instance_t, pm->instance_buffer);
}

/*
When a particle manager updates particles it creates work queues to handle the work. By default these each have a maximum amount of 1000 entries which should be
more than enough for most situations. However you can increase the sizes here if needed. You only need to set this manually if you hit one of the asserts when these
run out of space or you anticipate a huge amount of emitters and particles to be used (> million). On the other hand, you might be tight on memory in which case you
could reduce the numbers as well if needed (they don't take a lot of space though)
* @param pm                        A pointer to an intialised tfx_particle_manager_t.
* @param spawn_work_max            The maximum amount of spawn work entries
* @param control_work_max        The maximum amount of control work entries
* @param age_work_max            The maximum amount of age_work work entries
*/
tfxAPI void tfx_SetPMWorkQueueSizes(tfx_particle_manager_t *pm, tfxU32 spawn_work_max, tfxU32 control_work_max, tfxU32 age_work_max);

/*Free the memory for a specific emitter type. When an emitter is created it creates memory to store all of the particles that it updates each frame. If you have
multiple emitters of the same type then their particle lists are resused rather then freed as they expire. When they're freed then the unused list is added to a list
of free particle banks for that emitter type so that they can then be recycled if another emitter of the same type is created. If you want to free the memory for a
specific emitter then you can call this function to do that.
NOTE: No emitters of the type passed to the function must be in use in the particle manager.
* @param pm                        A pointer to an intialised tfx_particle_manager_t.
* @param emitter                A pointer to a valid tfx_effect_emitter_t of type tfxEmitterType
*/
tfxAPI void tfx_FreeParticleListsMemory(tfx_particle_manager_t *pm, tfx_effect_emitter_t *emitter);

/*
Free all the memory that is associated with an effect. Depending on the configuration of the particle manager this might be instance_data, particle lists and spawn location lists.
* @param pm                        A pointer to an intialised tfx_particle_manager_t.
* @param emitter                A pointer to a valid tfx_effect_emitter_t of type tfxEffectType
*/
tfxAPI void tfx_FreeEffectListsMemory(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect);

/*
Get the current particle count for a particle manager
* @param pm                        A pointer to an tfx_particle_manager_t
* @returns tfxU32                The total number of particles currently being updated
*/
tfxAPI tfxU32 tfx_ParticleCount(tfx_particle_manager_t *pm);

/*
Get the current number of effects that are currently being updated by a particle manager
* @param pm                        A pointer to an tfx_particle_manager_t
* @returns tfxU32                The total number of effects currently being updated
*/
tfxAPI tfxU32 tfx_EffectCount(tfx_particle_manager_t *pm);

/*
Get the current number of emitters that are currently being updated by a particle manager
* @param pm                        A pointer to an tfx_particle_manager_t
* @returns tfxU32                The total number of emitters currently being updated
*/
tfxAPI tfxU32 tfx_EmitterCount(tfx_particle_manager_t *pm);

/*
Set the seed for the particle manager for random number generation. Setting the seed can determine how an emitters spawns particles, so if you set the seed before adding an effect to the particle manager
then the effect will look the same each time. Note that seed of 0 is invalid, it must be 1 or greater.
* @param pm                            A pointer to an initialised tfx_particle_manager_t. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
* @param seed                        An unsigned int representing the seed (Any value other then 0)
*/
tfxAPI inline void tfx_SetSeed(tfx_particle_manager_t *pm, tfxU64 seed) {
	tfx_RandomReSeed(&pm->random, seed == 0 ? tfxMAX_UINT : seed);
	tfx_RandomReSeed(&pm->threaded_random, seed == 0 ? tfxMAX_UINT : seed);
}
 
/*
Prepare a tfx_effect_template_t that you can use to customise effects in the library in various ways before adding them into a particle manager for updating and rendering. Using a template like this
means that you can tweak an effect without editing the base effect in the library.
* @param library                    A reference to a tfx_library_t that should be loaded with tfx_LoadEffectLibrary
* @param name                        The name of the effect in the library that you want to use for the template. If the effect is in a folder then use normal pathing: "My Folder/My effect"
* @param effect_template            The empty tfx_effect_template_t object that you want the effect loading into
//Returns true on success.
*/
tfxAPI bool tfx_PrepareEffectTemplate(tfx_library_t *library, const char *name, tfx_effect_template_t *effect_template);

/*
Add an effect to a tfx_particle_manager_t from an effect template
* @param pm                    A pointer to an initialised tfx_particle_manager_t. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
* @param effect_template    The tfx_effect_template_t object that you want to add to the particle manager. It must have already been prepared by calling tfx_PrepareEffectTemplate
* @param effect_id            pointer to a tfxEffectID of the effect which will be set after it's been added to the particle manager. This index can then be used to manipulate the effect in the particle manager as it's update
							For example by calling tfx_SetEffectPosition2d. This will be set to tfxINVALID if the function is unable to add the effect to the particle manager if it's out of space and reached it's effect limit.
  @returns                    True if the effect was succesfully added.
*/
tfxAPI bool tfx_AddEffectTemplateToParticleManager(tfx_particle_manager_t *pm, tfx_effect_template_t *effect, tfxEffectID *effect_id);

/*
Add an effect to a tfx_particle_manager_t. Generally you should always call tfx_AddEffectTemplateToParticleManager and use templates to organise your effects but if you want to just
test things out you can add an effect direct from a library using this command.
* @param pm                    A pointer to an initialised tfx_particle_manager_t. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
* @param effect                tfx_effect_emitter_t object that you want to add to the particle manager.
* @param effect_id            pointer to a tfxEffectID of the effect which will be set after it's been added to the particle manager. This index can then be used to manipulate the effect in the particle manager as it's update
							For example by calling tfx_SetEffectPosition2d. This will be set to tfxINVALID if the function is unable to add the effect to the particle manager if it's out of space and reached it's effect limit.
  @returns                    True if the effect was succesfully added.
*/
tfxAPI bool tfx_AddRawEffectToParticleManager(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect, tfxEffectID *effect_id);

/*
Update a particle manager. If you are interpolating particles in the vertex shader then it's important to only call this function once per frame only and idealy in a fixed step loop.
That means that if your fixed loop has to run twice to catch up (because of low frame rates) then you should still only call this function once but you can multiply the elapsed time
by the number of ticks. The ellapsed time should be the amount of time that has passed since the last frame so in a fixed step loop this will simply be the update rate in millisecs.
For example if you're updating 60 frames per second then elapsed time would be 16.666667. Psuedo code would look something like this:

	TimerAccumulate(game->timer);
	int pending_ticks = TimerPendingTicks(game->timer);	//The number of times the update loop will run this frame.

	while (tfx_TimerDoUpdate(game->timer)) {
		if (pending_ticks > 0) {
			tfx_UpdateParticleManager(&game->pm, FrameLength * pending_ticks);
			//Set the pending ticks to 0 so we don't run the update again this frame.
			pending_ticks = 0;
		}

		TimerUnAccumulate(game->timer);
	}
	TimerSet(game->timer);	//Set the timer and calculate the interpolation value. You can pass that to a uniform or push constant for the shader

	//Only upload the sprite/billboard buffer to the gpu if the particle manager was updated.
	if (TimerUpdateWasRun(game->timer)) {
		RenderParticles(game->pm, game);
	}

* @param pm                    A pointer to an initialised tfx_particle_manager_t. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
*/
tfxAPI void tfx_UpdateParticleManager(tfx_particle_manager_t *pm, float elapsed);

/*
Get the image pointer for a sprite. Use this when rendering particles in your renderer. The pointer that is returned will be the pointer that you set in your shape loader function
used when loading an effect library. Generally you shouldn't need to use this function, simply copy the whole instance buffer in the particle manager to your staging buffer to be
copied to the gpu in one go.
* @param pm                    A pointer to an initialised tfx_particle_manager_t. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
* @param property_indexes    The value in the instance_data->property_indexs[i] when iterating over the instance_data in your render function
  @returns                    void* pointer to the image
*/
tfxAPI void *tfx_GetSpriteImagePointer(tfx_particle_manager_t *pm, tfxU32 property_indexes);

/*
Get the handle of the sprite. Use this when rendering particles in your renderer one sprite at a time.
* @param pm                    A pointer to an initialised tfx_particle_manager_t. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
* @param property_indexes      The value in the instance_data->property_indexs[i] when iterating over the instance_data in your render function
  @out_handle                  Pass in a pointer to a vec2 which will be loaded with the handle values
*/
tfxAPI void tfx_GetSpriteHandle(void *instance, float out_handle[2]);

/*
Get the total number of instances ready for rendering in the particle manager.
* @param pm                    A pointer to an initialised tfx_particle_manager_t.
*/
tfxAPI inline tfxU32 tfx_TotalSpriteCount(tfx_particle_manager_t *pm) {
	return pm->instance_buffer.current_size;
}

/*
Clear all particles, instance_data and effects in a particle manager. If you don't need to use the particle manager again then call tfx_FreeParticleManager to also
free all the memory associated with the particle manager.
* @param pm                        A pointer to an initialised tfx_particle_manager_t.
* @param free_particle_banks    Set to true if you want to free the memory associated with the particle banks and release back to the memory pool
*/
tfxAPI void tfx_ClearParticleManager(tfx_particle_manager_t *pm, bool free_particle_banks, bool free_sprite_buffers);

/*
Free all the memory used in the particle manager.
* @param pm                        A pointer to an initialised tfx_particle_manager_t.
*/
tfxAPI void tfx_FreeParticleManager(tfx_particle_manager_t *pm);

//[Effects functions for altering effects that are currently playing out in a particle manager]

/*
Expire an effect by telling it to stop spawning particles. This means that the effect will eventually be removed from the particle manager after all of it's remaining particles have expired.
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
*/
tfxAPI void tfx_SoftExpireEffect(tfx_particle_manager_t *pm, tfxEffectID effect_index);

/*
Soft expire all the effects in a particle manager so that the particles complete their animation first
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
*/
tfxAPI void tfx_SoftExpireAll(tfx_particle_manager_t *pm);

/*
Expire an effect by telling it to stop spawning particles and remove all associated particles immediately.
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
*/
tfxAPI void tfx_HardExpireEffect(tfx_particle_manager_t *pm, tfxEffectID effect_index);

/*
Get effect user data
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @returns                void* pointing to the user data set in the effect. See tfx_effect_template_t::SetUserData() and tfx__set_effect_user_data()
*/
tfxAPI void *tfx_GetEffectUserData(tfx_particle_manager_t *pm, tfxEffectID effect_index);

/*
More for use in the editor, this function updates emitter base values for any effects that are currently running after their graph values have been changed.
*/
tfxAPI void tfx_UpdatePMBaseValues(tfx_particle_manager_t *pm);

/*
Set the tfx_library_t that the particle manager will use to render instance_data and lookup all of the various properties required to update emitters and particles.
This is also set when you initialise a particle manager
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param lib            A pointer to a tfx_library_t
*/
tfxAPI void tfx_SetPMLibrary(tfx_particle_manager_t *pm, tfx_library_t *library);

/*
Set the particle manager camera. This is used to calculate particle depth if you're using depth ordered particles so it needs to be updated each frame.
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param front            An array of 3 floats representing a normalised 3d vector describing the direction that the camera is pointing
* @param position        An array of 3 floats representing the position of the camera in 3d space
*/
tfxAPI void tfx_SetPMCamera(tfx_particle_manager_t *pm, float front[3], float position[3]);

/*
Each effect in the particle manager can have bounding box which you can decide to keep updated or not if you wanted to do any offscreen culling of effects. Theres some
extra overhead to keep the bounding boxes updated but that can be made back if you have a number of effect particles offscreen that don't need to be drawn.
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param yesno            Set to true or false if you want the bounding boxes to be udpated.
*/
tfxAPI void tfx_KeepBoundingBoxesUpdated(tfx_particle_manager_t *pm, bool yesno);

/*
Set the effect user data for an effect already added to a particle manager
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index    The index of the effect that you want to expire. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param user_data        A void* pointing to the user_data that you want to store in the effect
*/
tfxAPI void tfx_SetEffectUserData(tfx_particle_manager_t *pm, tfxEffectID effect_index, void *user_data);

/*
Force a particle manager to only run in single threaded mode. In other words, only use the main thread to update particles
* @param pm                A pointer to a tfx_particle_manager_t.
* @param switch_on        true or false to use a single thread or not
*/
tfxAPI inline void tfx_ForcePMSingleThreaded(tfx_particle_manager_t *pm, bool switch_on) {
	if (switch_on) pm->flags |= tfxParticleManagerFlags_single_threaded; else pm->flags &= ~tfxParticleManagerFlags_single_threaded;
}

/*
Get the transform vectors for a 3d sprite's previous position so that you can use that to interpolate between that and the current sprite position
* @param pm                A pointer to a tfx_particle_manager_t.
* @param layer            The index of the sprite layer
* @param index            The sprite index of the sprite that you want the captured sprite for.
* @param position         This should be a pointer to a vec3 that you pass in that will get loaded with the position of the instance
*/
tfxAPI void tfx_GetCapturedInstance3dTransform(tfx_particle_manager_t *pm, tfxU32 layer, tfxU32 index, float out_position[3]);

/*
Get the transform vectors for a 2d sprite's previous position so that you can use that to interpolate between that and the current sprite position
* @param pm                A pointer to a tfx_particle_manager_t.
* @param layer            The index of the sprite layer
* @param index            The sprite index of the sprite that you want the captured sprite for.
*/
tfxAPI void tfx_GetCapturedInstance2dTransform(tfx_particle_manager_t *pm, tfxU32 layer, tfxU32 index, float out_position[3]);

/*
Get the index offset into the sprite memory for sprite data containing a pre recorded effect animation. Can be used along side tfx_SpriteDataEndIndex to create
a for loop to iterate over the instance_data in a pre-recorded effect
* @param sprite_data    A pointer to tfx_sprite_data_t containing all the instance_data and frame data
* @param frame            The index of the frame you want the end index for
* @param layer            The sprite layer
* @returns                tfxU32 containing the end offset
*/
tfxAPI tfxU32 tfx_SpriteDataEndIndex(tfx_sprite_data_t *sprite_data, tfxU32 frame, tfxU32 layer);

/*
Make a particle manager stop spawning. This will mean that all emitters in the particle manager will no longer spawn any particles so all currently running effects will expire
as the remaining particles come to the end of their life. Any single particles will also get flagged to expire
* @param pm                A pointer to a tfx_particle_manager_t.
* @param yesno            True = disable spawning, false = enable spawning
*/
tfxAPI inline void tfx_DisablePMSpawning(tfx_particle_manager_t *pm, bool yesno) {
	if (yesno) {
		pm->flags |= tfxParticleManagerFlags_disable_spawning;
	}
	else {
		pm->flags &= ~tfxParticleManagerFlags_disable_spawning;
	}
}

/*
Get the buffer of effect indexes in the particle manager.
* @param pm               A pointer to a tfx_particle_manager_t.
* @param depth            The depth of the list that you want. 0 are top level effects and anything higher are sub effects within those effects
* @param count			  A pointer to an int that you can pass in that will be filled with the count of effects in the array
* @returns                Pointer to the array of effect indexes
*/
tfxAPI tfx_effect_index_t *tfx_GetPMEffectBuffer(tfx_particle_manager_t *pm, tfxU32 depth, int *count);

/*
Get the buffer of emitter indexes in the particle manager.
* @param pm                A pointer to a tfx_particle_manager_t.
* @param depth            The depth of the list that you want. 0 are top level emitters and anything higher are sub emitters within those effects
* @param count			  A pointer to an int that you can pass in that will be filled with the count of emitters in the array
* @returns                Pointer to the tfxvec of effect indexes
*/
tfxAPI tfxU32 *tfx_GetPMEmitterBuffer(tfx_particle_manager_t *pm, tfxU32 depth, int *count);

/*
Set the position of a 2d effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param x                The x value of the position
* @param y                The y value of the position
*/
tfxAPI void tfx_SetEffectPosition2d(tfx_particle_manager_t *pm, tfxEffectID effect_index, float x, float y);

/*
Set the position of a 3d effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param x                The x value of the position
* @param y                The y value of the position
* @param z                The z value of the position
*/
tfxAPI void tfx_SetEffectPosition3d(tfx_particle_manager_t *pm, tfxEffectID effect_index, float x, float y, float z);

/*
Set the position of a 2d effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param position        A tfx_vec2_t vector object containing the x and y coordinates
*/
tfxAPI void tfx_SetEffectPositionVec2(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfx_vec2_t position);

/*
Set the position of a 3d effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param position        A tfx_vec3_t vector object containing the x, y and z coordinates
*/
tfxAPI void tfx_SetEffectPositionVec3(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfx_vec3_t position);

/*
Move an Effect by a specified amount relative to the effect's current position
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param amount            A tfx_vec3_t vector object containing the amount to move in the x, y and z planes
*/
tfxAPI void tfx_MoveEffectVec3(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfx_vec3_t amount);

/*
Move an Effect by a specified amount relative to the effect's current position
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param x                The amount to move in the x plane
* @param y                The amount to move in the y plane
* @param z                The amount to move in the z plane
*/
tfxAPI void tfx_MoveEffect3d(tfx_particle_manager_t *pm, tfxEffectID effect_index, float x, float y, float z);

/*
Get the current position of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @return                tfx_vec3_t containing the effect position
*/
tfxAPI void tfx_GetEffectPositionVec3(tfx_particle_manager_t *pm, tfxEffectID effect_index, float out_position[3]);

/*
You can use this function to get the sprite buffer of a specific effect. 
* @param pm						A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index			The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param tfxU32					Pass in a pointer to a tfxU32 which will be set to the number of instance_data in the buffer.
* @return						tfx_2d_instance_t pointer to the buffer
*/
tfxAPI tfx_2d_instance_t *tfx_GetEffect2dInstanceBuffer(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfxU32 *sprite_count);

/*
You can use this function to get the billboard buffer of a specific effect. 
* @param pm						A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index			The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param tfxU32					Pass in a pointer to a tfxU32 which will be set to the number of instance_data in the buffer.
* @return						tfx_3d_instance_t pointer to the buffer
*/
tfxAPI tfx_3d_instance_t *tfx_GetEffect3dInstanceBuffer(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfxU32 *sprite_count);

/*
You can use this function to get each sprite buffer for every effect that is currently active in the particle manager. Generally you would call this inside a for loop for each layer.
* @param pm						A pointer to a tfx_particle_manager_t where the effect is being managed
* @param tfxU32					The index of the sprite layer that you want
* @param tfx_2d_instance_t	Pass in a pointer which will be set to the current sprite buffer containing all of the sprite data for this frame.
* @param tfx_effect_instance_data_t   Pass in a second pointer which will be set to the tfx_effect_instance_data_t containing all of the sprite buffer data. This can be used to gain access to all the sprite data if using double buffered instance_data (to interpolated with the previous frame).
*                               You can use this with functions like GetCapturedEffectSprite3dTransform.
* @param tfxU32					Pass in a pointer to a tfxU32 which will be set to the number of instance_data in the buffer.
* @return						true or false if the next sprite buffer was found. False will be returned once there are no more effect sprite buffers in the particle manager
*/
tfxAPI bool tfx_GetNext2dInstanceBuffer(tfx_particle_manager_t *pm, tfx_2d_instance_t **sprites_soa, tfx_effect_instance_data_t **effect_sprites, tfxU32 *sprite_count);

/*
You can use this function to get each billboard buffer for every effect that is currently active in the particle manager. Generally you would call this inside a for loop for each layer.
* @param pm						A pointer to a tfx_particle_manager_t where the effect is being managed
* @param tfx_sprite_billboard_t	Pass in a pointer which will be set to the current sprite buffer containing all of the sprite data for this frame.
* @param tfx_effect_instance_data_t   Pass in a second pointer which will be set to the tfx_effect_instance_data_t containing all of the sprite buffer data. This can be used to gain access to all the sprite data if using double buffered instance_data (to interpolated with the previous frame).
*                               You can use this with functions like GetCapturedEffectSprite3dTransform.
* @param tfxU32					Pass in a pointer to a tfxU32 which will be set to the number of instance_data in the buffer.
* @return						true or false if the next billboard buffer was found. False will be returned once there are no more effect sprite buffers in the particle manager
*/
tfxAPI bool tfx_GetNext3dInstanceBuffer(tfx_particle_manager_t *pm, tfx_3d_instance_t **sprites_soa, tfx_effect_instance_data_t **effect_sprites, tfxU32 *sprite_count);

/*After calling GetNextBillboard/SpriteBuffer in a while loop you can call this to reset the index for the next frame
* @param pm						A pointer to a tfx_particle_manager_t
*/
tfxAPI void tfx_ResetInstanceBufferLoopIndex(tfx_particle_manager_t *pm);

/*
Set the rotation of a 2d effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current rotation of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param rotation        A float of the amount that you want to set the rotation too
*/
tfxAPI void tfx_SetEffect2dRotation(tfx_particle_manager_t *pm, tfxEffectID effect_index, float rotation);

/*
Set the roll of a 3d effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current roll of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param roll            A float of the amount that you want to set the roll too
*/
tfxAPI void tfx_SetEffectRoll(tfx_particle_manager_t *pm, tfxEffectID effect_index, float roll);

/*
Set the pitch of a 3d effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current pitch of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param pitch            A float of the amount that you want to set the pitch too
*/
tfxAPI void tfx_SetEffectPitch(tfx_particle_manager_t *pm, tfxEffectID effect_index, float pitch);

/*
Set the yaw of a 3d effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current yaw of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param yaw            A float of the amount that you want to set the yaw too
*/
tfxAPI void tfx_SetEffectYaw(tfx_particle_manager_t *pm, tfxEffectID effect_index, float yaw);

/*
Set the width of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current width of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param width            A float of the amount that you want to set the width multiplier too. The width multiplier will multiply all widths of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void tfx_SetEffectWidthMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float width);

/*
Set the height of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current height of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param height            A float of the amount that you want to set the height multiplier too. The height multiplier will multiply all heights of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void tfx_SetEffectHeightMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float height);

/*
Set the depth of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current depth of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param depth            A float of the amount that you want to set the depth multiplier too. The depth multiplier will multiply all heights of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void tfx_SetEffectDepthMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float depth);

/*
Set the life multiplier of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current life of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param life            A float of the amount that you want to set the life multiplier too. The life mulitplier will affect how long all particles emitted within the effect will last before expiring.
*/
tfxAPI void tfx_SetEffectLifeMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float life);

/*
Set the particle width multiplier of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current particle width of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param width            A float of the amount that you want to set the particle width multiplier too. The particle width mulitplier will affect the width of each particle if the emitter has a non uniform particle size, otherwise
						it will uniformly size the particle
*/
tfxAPI void tfx_SetEffectParticleWidthMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float width);

/*
Set the particle height multiplier of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current particle width of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param height            A float of the amount that you want to set the particle height multiplier too. The particle height mulitplier will affect the height of each particle if the emitter has a non uniform particle size, otherwise
						this function will have no effect.
*/
tfxAPI void tfx_SetEffectParticleHeightMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float height);

/*
Set the velocity multiplier of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current velocity of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param velocity        A float of the amount that you want to set the particle velocity multiplier too. The particle velocity mulitplier will affect the base velocity of a particle at spawn time.
*/
tfxAPI void tfx_SetEffectVelocityMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float velocity);

/*
Set the spin multiplier of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current spin of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param spin            A float of the amount that you want to set the particle spin multiplier too. The particle spin mulitplier will affect the base spin of a particle at spawn time.
*/
tfxAPI void tfx_SetEffectSpinMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float spin);

/*
Set the intensity multiplier of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current intensity of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param intensity        A float of the amount that you want to set the particle intensity multiplier too. The particle intensity mulitplier will instantly affect the opacity of all particles currently emitted by the effect.
*/
tfxAPI void tfx_SetEffectIntensityMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float intensity);

/*
Set the splatter multiplier of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current splatter of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param splatter        A float of the amount that you want to set the particle splatter multiplier too. The particle splatter mulitplier will change the amount of random offset all particles emitted in the effect will have.
*/
tfxAPI void tfx_SetEffectSplatterMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float splatter);

/*
Set the weight multiplier of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current weight of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param weight            A float of the amount that you want to set the particle weight multiplier too. The particle weight mulitplier will change the weight applied to particles in the effect at spawn time.
*/
tfxAPI void tfx_SetEffectWeightMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float weight);

/*
Set the overal scale of an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after tfx_UpdateParticleManager in order to override the current weight of the effect that was
*                        set in the TimelineFX editor.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param overal_scale    A float of the amount that you want to set the overal scale to. The overal scale is an simply way to change the size of an effect
*/
tfxAPI void tfx_SetEffectOveralScale(tfx_particle_manager_t *pm, tfxEffectID effect_index, float overal_scale);

/*
Set the base noise offset for an effect
* @param pm                A pointer to a tfx_particle_manager_t where the effect is being managed.
* @param effect_index    The index of the effect. This is the index returned when calling tfx_AddEffectTemplateToParticleManager
* @param noise_offset    A float of the amount that you want to set the effect noise offset to. By default when an effect is added to a particle manager a random noise offset will be set based on the Base Noise Offset Range property. Here you can override that
						value by setting it here. The most ideal time to set this would be immediately after you have added the effect to the particle manager, but you could call it any time you wanted for a constantly changing noise offset.
*/
tfxAPI void tfx_SetEffectBaseNoiseOffset(tfx_particle_manager_t *pm, tfxEffectID effect_index, float noise_offset);

/*
Get the name of an effect
* @param pm                A pointer to the effect
* @returns                const char * name
*/
tfxAPI inline const char *tfx_GetEffectName(tfx_effect_emitter_t *effect);


//--------------------------------
//Animation_manager
//--------------------------------

/*
Set the position of a 3d animation
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @param position                A tfx_vec3_t vector object containing the x, y and z coordinates
*/
tfxAPI void tfx_SetAnimationPosition3d(tfx_animation_manager_t *animation_manager, tfxAnimationID animation_id, float position[3]);

/*
Set the position of a 2d animation
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @param x                        A float of the x position
* @param y                        A float of the y position
*/
tfxAPI void tfx_SetAnimationPosition2d(tfx_animation_manager_t *animation_manager, tfxAnimationID animation_id, float x, float y);

/*
Set the scale of a 3d animation
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @param scale                    A multiplier that will determine the overal size/scale of the effect
*/
tfxAPI void tfx_SetAnimationScale(tfx_animation_manager_t *animation_manager, tfxAnimationID animation_id, float scale);

/*
Get an animation instance from an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param tfxAnimationID            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @returns pointer to instance    Pointer to a tfx_animation_instance_t
*/
tfxAPI tfx_animation_instance_t *tfx_GetAnimationInstance(tfx_animation_manager_t *animation_manager, tfxAnimationID animation_id);

/*
Initialise an Animation Manager for use with 3d instance_data. This must be run before using an animation manager. An animation manager is used
to playback pre recorded particle effects as opposed to using a particle manager that simulates the particles in
real time. This pre-recorded data can be uploaded to the gpu for a compute shader to do all the interpolation work
to calculate the state of particles between frames for smooth animation.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param max_instances            The maximum number of animation instances that you want to be able to play at one time.
* @param initial_capacity        Optionally, you can set an initial capacity for the sprite data. The data will grow if you add
								beyond this amount but it gives you a chance to reserve a decent amount to start with to
								save too much mem copies as the data grows
*/
tfxAPI void tfx_InitialiseAnimationManagerFor3d(tfx_animation_manager_t *animation_manager, tfxU32 max_instances, tfxU32 initial_sprite_data_capacity);

/*
Initialise an Animation Manager for use with 2d instance_data. This must be run before using an animation manager. An animation manager is used
to playback pre recorded particle effects as opposed to using a particle manager that simulates the particles in
real time. This pre-recorded data can be uploaded to the gpu for a compute shader to do all the interpolation work
to calculate the state of particles between frames for smooth animation.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param max_instances            The maximum number of animation instances that you want to be able to play at one time.
* @param initial_capacity        Optionally, you can set an initial capacity for the sprite data. The data will grow if you add
								beyond this amount but it gives you a chance to reserve a decent amount to start with to
								save too much mem copies as the data grows
*/
tfxAPI void tfx_InitialiseAnimationManagerFor2d(tfx_animation_manager_t *animation_manager, tfxU32 max_instances, tfxU32 initial_sprite_data_capacity);

/*
Set the callback that you can use to determine whether or not a tfx_animation_instance_t should be added to the next frame's render queue. You can use this
to cull instances that are outside of the view frustum for example
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param callback                Pointer to the callback you want to use. It must have the following signature:
								bool(*maybe_render_instance_callback(tfx_animation_manager_t *animation_manager, tfx_animation_instance_t *instance, tfx_frame_meta_t *meta, void *user_data))
								Values passed into the callback function are a pointer to the animation manager, a pointer to the instance being processed, a pointer to
								the frame meta of the instance, this will contain the bounding box and radius of the instance from the current frame of the instance and a pointer
								to any user data that you set that might contain the camera frustum that you want to check against.
*/
tfxAPI void tfx_SetAnimationManagerInstanceCallback(tfx_animation_manager_t *animation_manager, bool((*maybe_render_instance_callback)(tfx_animation_manager_t *animation_manager, tfx_animation_instance_t *instance, tfx_frame_meta_t *meta, void *user_data)));

/*
Get the sprite data settings for an effect in a library. Sprite data settings are the settings for an effect in the editor relating to setting up pre-baked effects
* @param library				Pointer to the tfx_library_t where the effect is stored
* @param effect					Pointer the the effect that you want the sprite settings for.
* @returns						Pointer to the tfx_sprite_data_settings
*/
tfx_sprite_data_settings_t *tfx_GetEffectSpriteDataSettings(tfx_library_t *library, tfx_effect_emitter_t *effect);

/*
Get the sprite data settings for an effect in a library by it's path. Sprite data settings are the settings for an effect in the editor relating to setting up pre-baked effects
* @param library				Pointer to the tfx_library_t where the effect is stored
* @param path					const char* string of the path to the effect. Must be the path to a root effect.
* @returns						Pointer to the tfx_sprite_data_settings
*/
tfx_sprite_data_settings_t *tfx_GetEffectSpriteDataSettingsByPath(tfx_library_t *library, const char *path);

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
Set the user data in a tfx_animation_manager_t which can get passed through to callback functions when updated the animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param user_data                void* pointer to the data that you want to set
*/
tfxAPI void tfx_SetAnimationManagerUserData(tfx_animation_manager_t *animation_manager, void *user_data);

/*
Add sprite data to an animation manager sprite data buffer from an effect. This will record the
animation if necessary and then convert the sprite data to tfx_sprite_data3d_t ready for uploading
to the GPU
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index            The index of the effect. This is the index returned when calling tfx_AddAnimationInstance
* @param position                A tfx_vec3_t vector object containing the x, y and z coordinates
*/
tfxAPI void tfx_AddSpriteData(tfx_animation_manager_t *animation_manager, tfx_effect_emitter_t *effect, tfx_particle_manager_t *pm, tfx_vec3_t camera_position);

/*
Add an animation instance to the animation manager.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param path                    tfxKey path hash of the effect name and path: effect.path_hash
* @param start_frame            Starting frame of the animation
* @returns                        The index id of the animation instance. You can use this to reference the animation when changing position, scale etc
								Return tfxINVALID if there is no room in the animation manager
*/
tfxAPI tfxAnimationID tfx_AddAnimationInstanceByKey(tfx_animation_manager_t *animation_manager, tfxKey path, tfxU32 start_frame);

/*
Add an animation instance to the animation manager.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param path                    const char * name of the effect. If the effect was in a folder then specify the whole path
* @param start_frame            Starting frame of the animation
* @returns                        The index id of the animation instance. You can use this to reference the animation when changing position, scale etc
								Return tfxINVALID if there is no room in the animation manager
*/
tfxAPI tfxAnimationID tfx_AddAnimationInstance(tfx_animation_manager_t *animation_manager, const char *path, tfxU32 start_frame);

/*
Update an animation manager to advance the time and frames of all instances currently playing.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to update
* @param start_frame            Starting frame of the animation
*/
tfxAPI void tfx_UpdateAnimationManager(tfx_animation_manager_t *animation_manager, float elapsed);

/*
Add an effect's shapes to an animation manager. You can use this function if you're manually recording particle effects and adding them to an animation
manager rather then just using the editor.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to update
* @param effect                    A pointer to the effect whose shapes you want to add
*/
tfxAPI void tfx_AddEffectShapes(tfx_animation_manager_t *animation_manager, tfx_effect_emitter_t *effect);

/*
Update an animation manager so that the effects do not expire they just loop forever instead regardless of whether they're a looped effect or not.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to update
*/
tfxAPI void tfx_CycleAnimationManager(tfx_animation_manager_t *animation_manager);

/*
Clears all animation instances currently in play in an animation manager, resulting in all currently running animations
from being drawn
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to clear
*/
tfxAPI void tfx_ClearAllAnimationInstances(tfx_animation_manager_t *animation_manager);

/*
Clears all data from the animation manager including sprite data, metrics and instances. Essentially resetting everything back to
it's initialisation point
from being drawn
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to reset
*/
tfxAPI void tfx_ResetAnimationManager(tfx_animation_manager_t *animation_manager);

/*
Frees all data from the animation manager including sprite data, metrics and instances.
from being drawn
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to reset
*/
tfxAPI void tfx_FreeAnimationManager(tfx_animation_manager_t *animation_manager);

/*
Get the tfx_animation_buffer_metrics_t from an animation manager. This will contain the info you need to upload the sprite data,
offsets and animation instances to the GPU. Only offsets and animation instances need to be uploaded to the GPU each frame. Sprite
data can be done ahead of time.
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @returns                        tfx_animation_buffer_metrics_t containing buffer sizes
*/
tfxAPI inline tfx_animation_buffer_metrics_t tfx_GetAnimationBufferMetrics(tfx_animation_manager_t *animation_manager) {
	return animation_manager->buffer_metrics;
}

/*
Get the total number of instance_data that need to be drawn by an animation manager this frame. You can use this in your renderer
to draw your sprite instances
* @param animation_manager        A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @returns                        tfxU32 of the number of instance_data
*/
tfxAPI inline tfxU32 tfx_GetTotalSpritesThatNeedDrawing(tfx_animation_manager_t *animation_manager) {
	return animation_manager->buffer_metrics.total_sprites_to_draw;
}

/*
Get the total number of instances being processed by an animation manager. This will not necessarily be the same number as
the instances being rendered if some are being culled in your custom callback if your using one.
* @param animation_manager        A pointer to a tfx_animation_manager_t that you want to clear
* @returns int                    The number of instances being updated
*/
tfxAPI tfxU32 tfx_GetTotalInstancesBeingUpdated(tfx_animation_manager_t *animation_manager);

/*
Create the image data required for shaders from a TimelineFX library. The image data will contain data such as uv coordinates. Once you have built the data you can use GetLibraryImageData to get the buffer
and upload it to the gpu.
* @param animation_manager		  A pointer to an tfx_animation_manager_t object
* @param shapes                   A pointer to a tfx_gpu_shapes_t object which will fill a buffer with all the shapes
* @param uv_lookup                A function pointer to a function that you need to set up in order to get the uv coordinates from whatever renderer you're using
*/
tfxAPI void tfx_BuildAnimationManagerGPUShapeData(tfx_animation_manager_t *animation_manager, tfx_gpu_shapes_t *shapes, void(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset));

/*
Get a pointer to the GPU shapes which you can use in a memcpy
* @param particle_shapes        A pointer the tfx_gpu_shapes_t
*/
tfxAPI inline void *tfx_GetGPUShapesPointer(tfx_gpu_shapes_t *particle_shapes) {
	return particle_shapes->list.data;
}

/*
Get a pointer to the particle shapes data in the animation manager. This can be used with tfx_BuildGPUShapeData when you want to upload the data to the GPU
* @param animation_manager        A pointer the tfx_animation_manager_t
*/
tfxAPI inline tfx_image_data_t *tfx_GetParticleShapesAnimationManager(tfx_animation_manager_t *animation_manager, int *count) {
	*count = animation_manager->particle_shapes.data.current_size;
	return animation_manager->particle_shapes.data.data;
}

/*
Get the number of shapes in the GPU Shape Data buffer. Make sure you call tfx_BuildGPUShapeData first or they'll be nothing to return
* @param library                A pointer to a tfx_animation_manager_t where the image data will be created.
* @returns tfxU32                The number of shapes in the buffer
*/
tfxAPI tfxU32 tfx_GetGPUShapeCount(tfx_gpu_shapes_t *particle_shapes);

/*
Get the size in bytes of the GPU image data in a tfx_library_t
* @param library                A pointer to a tfx_library_t where the image data exists.
* @returns size_t                The size in bytes of the image data
*/
tfxAPI size_t tfx_GetGPUShapesSizeInBytes(tfx_gpu_shapes_t *particle_shapes);

/*
Get the total number of instance_data in an animation manger's sprite data buffer
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32                The number of instance_data in the buffer
*/
tfxAPI inline tfxU32 tfx_GetTotalSpriteDataCount(tfx_animation_manager_t *animation_manager) {
	if (animation_manager->flags & tfxAnimationManagerFlags_is_3d) {
		return animation_manager->sprite_data_3d.current_size;
	}
	return animation_manager->sprite_data_2d.current_size;
}

/*
Get the total number of instance_data in an animation manger's sprite data buffer
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32                The number of instance_data in the buffer
*/
tfxAPI inline size_t tfx_GetSpriteDataSizeInBytes(tfx_animation_manager_t *animation_manager);

/*
Get the buffer memory address for the sprite data in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                A pointer to the sprite data memory
*/
tfxAPI inline void *tfx_GetSpriteDataBufferPointer(tfx_animation_manager_t *animation_manager) {
	if (animation_manager->flags & tfxAnimationManagerFlags_is_3d) {
		return animation_manager->sprite_data_3d.data;
	}
	return animation_manager->sprite_data_2d.data;
}

/*
Get the size in bytes of the offsets buffer in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                Size in bytes of the offsets buffer
*/
tfxAPI inline size_t tfx_GetOffsetsSizeInBytes(tfx_animation_manager_t *animation_manager) {
	return animation_manager->offsets.current_size * sizeof(tfxU32);
}

/*
Get the size in bytes of the render queue of animation instances buffer in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                Size in bytes of the instances buffer
*/
tfxAPI inline size_t tfx_GetAnimationInstancesSizeInBytes(tfx_animation_manager_t *animation_manager) {
	return animation_manager->render_queue.current_size * sizeof(tfx_animation_instance_t);
}

/*
Get the size in bytes of the animation emitter properties list
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t                Size in bytes of the properties bufffer
*/
tfxAPI inline size_t tfx_GetAnimationEmitterPropertySizeInBytes(tfx_animation_manager_t *animation_manager) {
	return animation_manager->emitter_properties.current_size * sizeof(tfx_animation_emitter_properties_t);
}

/*
Get the number of emitter properties being using by the animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32                Number of emitter properties
*/
tfxAPI inline tfxU32 tfx_GetAnimationEmitterPropertyCount(tfx_animation_manager_t *animation_manager) {
	return animation_manager->emitter_properties.current_size;
}

/*
Get the buffer memory address for the sprite data in an animation manager
* @param animation_manager        A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*                A pointer to the sprite data memory
*/
tfxAPI void *tfx_GetAnimationEmitterPropertiesBufferPointer(tfx_animation_manager_t *animation_manager);

//--------------------------------
//Effect_templates
//--------------------------------

/*
Reset an effect template and make it empty so you can use it to store another effect.
* @param t                        A pointer to a tfx_effect_template_t
*/
tfxAPI void tfx_ResetTemplate(tfx_effect_template_t *t);

/*
Get the root effect from the template
* @param t                        A pointer to a tfx_effect_template_t
* @returns                        A pointer to the root effect
*/
tfxAPI tfx_effect_emitter_t *tfx_GetEffectFromTemplate(tfx_effect_template_t *t);

/*
Get an emitter or sub effect from an effect template.
* @param t                        A pointer to a tfx_effect_template_t
* @param path                    A path to the emitter or sub effect that you want to retrieve. Must be a valid path. Example path might be: "Explosion/Smoke"
* @returns                        A pointer to the root effect
*/
tfxAPI tfx_effect_emitter_t *tfx_GetEmitterFromTemplate(tfx_effect_template_t *t, const char *path);

/*
Get an emitter path that an emitter is using. The emitter must have the path emission type set or nullptr will be returned
* @param t                        A pointer to a tfx_effect_emitter_t
* @param path                    A path to the emitter or sub effect that you want to retrieve. Must be a valid path. Example path might be: "Explosion/Smoke"
* @returns                        A pointer to the root effect
*/
tfxAPI tfx_emitter_path_t *tfx_GetEmitterPath(tfx_effect_emitter_t *e);

/*
Set the user data for any effect or emitter in the effect template. This user data will get passed through to any update callback functions
* @param t                        A pointer to a tfx_effect_template_t
* @param path                    A path to the effect or emitter in the effect template
* @param data                    A pointer to the user data
*/
tfxAPI void tfx_SetTemplateUserData(tfx_effect_template_t *t, const char *path, void *data);

/*
Set the user data for the root effect in an effect template
* @param t                        A pointer to a tfx_effect_template_t
* @param data                    A pointer to the user data
*/
tfxAPI void tfx_SetTemplateEffectUserData(tfx_effect_template_t *t, void *data);

/*
Set the same user data for all effects and emitters/sub effects in the effect template
* @param t                        A pointer to a tfx_effect_template_t
* @param data                    A pointer to the user data that will be set to all effects and emitters in the template
*/
tfxAPI void tfx_SetTemplateUserDataAll(tfx_effect_template_t *t, void *data);

/*
Set an update callback for the root effect in the effect template.
* @param t                        A pointer to a tfx_effect_template_t
* @param update_callback        A pointer to the call back function
*/
tfxAPI void tfx_SetTemplateEffectUpdateCallback(tfx_effect_template_t *t, void(*update_callback)(tfx_particle_manager_t *pm, tfxEffectID effect_index));

/*
Pre-record this effect into a sprite cache so that you can play the effect back without the need to actually caclulate particles in realtime.
	* @param pm            Reference to a pm that will be used to run the particle simulation and record the sprite data
	* @param path        const *char of a path to the emitter in the effect.Must be a valid path, for example: "My Effect/My Emitter"
	* @param camera        Array of 3 floats with the camera position (only needed for 3d effects that are sorted by depth
*/
tfxAPI void tfx_RecordTemplateEffect(tfx_effect_template_t *t, tfx_particle_manager_t *pm, float update_frequency, float camera_position[3]);

/*
Disable an emitter within an effect. Disabling an emitter will stop it being added to the particle manager when calling tfx_AddEffectTemplateToParticleManager
* @param path        const *char of a path to the emitter in the effect. Must be a valid path, for example: "My Effect/My Emitter"
*/
tfxAPI void tfx_DisableTemplateEmitter(tfx_effect_template_t *t, const char *path);

/*
Enable an emitter within an effect so that it is added to the particle manager when calling tfx_AddEffectTemplateToParticleManager. Emitters are enabled by default.
* @param path        const *char of a path to the emitter in the effect. Must be a valid path, for example: "My Effect/My Emitter"
*/
tfxAPI void tfx_EnableTemplateEmitter(tfx_effect_template_t *t, const char *path);

/*
Scale all nodes on a global graph graph of the effect
* @param global_type        tfx_graph_type of the global graph that you want to scale. Must be a global graph or an assert will be called
* @param amount                A float of the amount that you want to scale the multiplier by.
*/
tfxAPI void tfx_ScaleTemplateGlobalMultiplier(tfx_effect_template_t *t, tfx_graph_type global_type, float amount);

/*
Scale all nodes on an emitter graph
* @param emitter_path        const *char of the emitter path
* @param global_type        tfx_graph_type of the emitter graph that you want to scale. Must be an emitter graph or an assert will be called
* @param amount                A float of the amount that you want to scale the graph by.
*/
tfxAPI void tfx_ScaleTemplateEmitterGraph(tfx_effect_template_t *t, const char *emitter_path, tfx_graph_type graph_type, float amount);

/*
Set the single spawn amount for an emitter. Only affects emitters that have the single spawn flag set.
* @param emitter_path        const *char of the emitter path
* @param amount                A float of the amount that you want to set the single spawn amount to.
*/
tfxAPI void tfx_SetTemplateSingleSpawnAmount(tfx_effect_template_t *t, const char *emitter_path, tfxU32 amount);

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
tfxAPI inline void tfx_Lerp3d(float lerp, const tfx_vec3_t *world, const tfx_vec3_t *captured, float out_lerp[3]);

/*
Interpolate between 2 tfxVec2s. You can make use of this in your render function when rendering instance_data and interpolating between captured and current positions
* @param tween        The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param world        The current tvxVec2 position
* @param captured    The captured tvxVec2 position
* @returns tfx_vec2_t    The interpolated tfx_vec2_t
*/
tfxAPI inline void tfx_Lerp2d(float lerp, const tfx_vec2_t *world, const tfx_vec2_t *captured, float out_lerp[2]);

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

/*
Check if a particle sprite is newly spawned. This means that there will be no captured index to interpolate with so if you want you can opt to not draw the sprite or
draw the sprite but with 0 alpha. A float is returned, either 0.f or 1.f so you can use that to multiply the alpha value or scale of the sprite to not draw it.
* @param instance_data    A pointer to a tfx_sprite_soa_t
* @param index        The index of the sprite that you're checking
* @returns float    0.f if it IS the first frame of the sprite otherwise 1.f.
*/
tfxAPI float tfx_IsFirstFrame(tfx_sprite_soa_t *sprites, tfxU32 sprite_index);

tfxAPI void tfx_GetSpriteScale(void *instance, float out_scale[2]);

#ifdef tfxINTEL
/*
Interpolate between 2 colors in tfx_rgba8_t format. You can make use of this in your render function when rendering instance_data and interpolating between captured and current colors
* @param tween						The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param current					The current tfx_rgba8_t color
* @param captured					The captured tfx_rgba8_t color
* @returns tfx_rgba8_t				The interpolated tfx_rgba8_t
*/
tfxAPI inline tfx_rgba8_t tfx_TweenColor(float tween, const tfx_rgba8_t current, const tfx_rgba8_t captured) {
	tfx128 color1 = _mm_set_ps((float)current.a, (float)current.b, (float)current.g, (float)current.r);
	tfx128 color2 = _mm_set_ps((float)captured.a, (float)captured.b, (float)captured.g, (float)captured.r);
	tfx128 wide_tween = _mm_set1_ps(tween);
	tfx128 wide_tween_m1 = _mm_sub_ps(_mm_set1_ps(1.f), wide_tween);
	color1 = _mm_div_ps(color1, _mm_set1_ps(255.f));
	color2 = _mm_div_ps(color2, _mm_set1_ps(255.f));
	color1 = _mm_mul_ps(color1, wide_tween);
	color2 = _mm_mul_ps(color2, wide_tween_m1);
	color1 = _mm_add_ps(color1, color2);
	color1 = _mm_mul_ps(color1, _mm_set1_ps(255.f));
	tfx128iArray packed;
	packed.m = _mm_cvtps_epi32(color1);
	tfx_rgba8_t color = { (unsigned char)packed.a[0], (unsigned char)packed.a[1], (unsigned char)packed.a[2], (unsigned char)packed.a[3] };
	return color;
}

/*
Interpolate all sprite transform data in a single function. This will interpolate position, scale and rotation.
* @param tween                The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param current            The current transform struct of the sprite
* @param captured            The captured transform struct of the sprite
* @returns tfx_wide_lerp_transform_result_t            The interpolated transform data in a tfx_wide_lerp_transform_result_t
*/
tfxAPI inline tfx_wide_lerp_transform_result_t tfx_InterpolateSpriteTransform(const __m128 *tween, const tfx_sprite_transform3d_t *current, const tfx_sprite_transform3d_t *captured) {
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
Interpolate between 2 colors in tfx_rgba8_t format. You can make use of this in your render function when rendering instance_data and interpolating between captured and current colors
* @param tween                The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param current            The current tfx_rgba8_t color
* @param captured            The captured tfx_rgba8_t color
* @returns tfx_rgba8_t            The interpolated tfx_rgba8_t
*/
tfxAPI inline tfx_rgba8_t tfx_TweenColor(float tween, const tfx_rgba8_t current, const tfx_rgba8_t captured) {
	tfx128 color1 = { (float)current.a, (float)current.b, (float)current.g, (float)current.r };
	tfx128 color2 = { (float)captured.a, (float)captured.b, (float)captured.g, (float)captured.r };
	tfx128 wide_tween = vdupq_n_f32(tween);
	tfx128 wide_tween_m1 = vsubq_f32(vdupq_n_f32(1.f), wide_tween);
	color1 = vmulq_f32(vdivq_f32(color1, vdupq_n_f32(255.f)), wide_tween);
	color2 = vmulq_f32(vdivq_f32(color2, vdupq_n_f32(255.f)), wide_tween_m1);
	color1 = vaddq_f32(color1, color2);
	color1 = vmulq_f32(color1, vdupq_n_f32(255.f));
	tfx128i packed = vcvtq_s32_f32(color1);
	return tfx_rgba8_t(vgetq_lane_s32(packed, 0), vgetq_lane_s32(packed, 1), vgetq_lane_s32(packed, 2), vgetq_lane_s32(packed, 3));
}


/*
Interpolate all sprite transform data in a single function. This will interpolate position, scale and rotation.
* @param tween                The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param current            The current transform struct of the sprite
* @param captured            The captured transform struct of the sprite
* @returns tfx_wide_lerp_transform_result_t            The interpolated transform data in a tfx_wide_lerp_transform_result_t
*/
tfxAPI inline tfx_wide_lerp_transform_result_t tfx_InterpolateSpriteTransform(const tfx128 *tween, const tfx_sprite_transform3d_t *current, const tfx_sprite_transform3d_t *captured) {
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

tfxAPI tfx_effect_emitter_t tfx_NewEffect();

#endif
