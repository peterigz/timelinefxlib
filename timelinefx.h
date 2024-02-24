#ifndef TFX_LIBRARY_HEADER
#define TFX_LIBRARY_HEADER

#define tfxENABLE_PROFILING
#define tfxPROFILER_SAMPLES 60
#define TFX_THREAD_SAFE
//#define TFX_EXTRA_DEBUGGING
//#define tfxUSEAVX

/*
	Timeline FX C++ library

	This library is for implementing particle effects into your games and applications.

	This library is render agnostic, so you will have to provide your own means to render the particles. You will use ParticleManager::GetParticleBuffer() to get all of the active particles in the particle manager
	and then use the values in Particle struct to draw a correctly scaled and rotated particle.

	Sections in this header file, you can search for the following keywords to jump to that section:

	[Zest_Pocket_Allocator]			A single header library for allocating memory from a large pool.
	[Header_Includes_and_Typedefs]	Just your basic header stuff for setting up typedefs and some #defines
	[OS_Specific_Functions]			OS specific multithreading and file access
	[XXHash_Implementation]			XXHasher for the storage map.
	[SIMD_defines]					Defines for SIMD intrinsics
	[Enums]							All the definitions for enums and bit flags
	[Constants]						Various constant definitions
	[String_Buffers]				Basic string buffers for storing names of things in the library and reading from library files.
	[Containers_and_Memory]			Container structs and lists and defines for memory is allocated (uses Zest Pocket Allocator by default)
	[Multithreading_Work_Queues]	Implementation for work queues for threading
	[Vector_Math]					Vec2/3/4 and Matrix2/3/4 structs including wide vectors for SIMD
	[Simplex_Noise]					Some setup for implementing simplex noise.
	[Profiling]						Very basic profiling for internal use
	[File_IO]						A package manager for reading/writing files such as a tfx library effects file
	[Struct_Types]					All of the structs used for objects in TimelineFX
	[Internal_Functions]			Mainly internal functions called only by the library but also the Editor, these are marked either tfxINTERNAL or tfxAPI_EDITOR
	[API_Functions]					The main functions for use by users of the library
*/

//Override this if you'd prefer a different way to allocate the pools for sub allocation in host memory.
#ifndef tfxALLOCATE_POOL
#define tfxALLOCATE_POOL malloc
#endif

#ifndef tfxMAX_MEMORY_POOLS
#define tfxMAX_MEMORY_POOLS 32
#endif

#if defined(__x86_64__) || defined(__i386__)
#define tfxINTEL
#include <immintrin.h>
#elif defined(__arm__) || defined(__aarch64__)
#include <arm_neon.h>
#define tfxARM
#endif

#include <stdint.h>
#include <math.h>

//type defs
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
		//Note to self: Crashing here? The most likely reason is a pointer into the allocation for this block that became invalid but was still written to at some point.
		//Most likeyly cause is a tfx_vector_t or similar being resized and allocated elsewhere but you didn't account for this happening and update the pointer. Just index
		//into the array instead to fix these issues.
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
    How to use:
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

#ifdef tfxINTEL
//Intel Intrinsics
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
#define tfxWideOri _mm_or_si128
#define tfxWideXOri _mm_xor_si128
#define tfxWideAnd _mm_and_ps
#define tfxWideAndi _mm_and_si128
#define tfxWideAndNoti _mm_andnot_si128
#define tfxWideSetZeroi _mm_setzero_si128
#define tfxWideSetZero _mm_setzero_ps
#define tfxWideEqualsi _mm_cmpeq_epi32 
#define tfxWideEquals _mm_cmpeq_ps
#define tfxWideShufflei _mm_shuffle_epi32

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

#elifdef tfxARM
//Arm Intrinsics
typedef float32x4_t tfxWideFloat;
typedef int32x4_t tfxWideInt;
#define tfxWideLoad vld1q_f32
#define tfxWideLoadi vld1q_s32
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
#define tfxWideSqrt vrsqrteq_f32 // for reciprocal square root approximation
#define tfxWideShiftRight vshrq_n_s32
#define tfxWideShiftLeft vshlq_n_s32
#define tfxWideGreaterEqual vreinterpretq_f32_u32(vcgeq_f32)
#define tfxWideGreater vreinterpretq_f32_u32(vcgtq_f32)
#define tfxWideGreateri vcgtq_s32
#define tfxWideLessEqual vreinterpretq_f32_u32(vcleq_f32)
#define tfxWideLess vreinterpretq_f32_u32(vcltq_f32)
#define tfxWideLessi vcltq_s32
#define tfxWideStore vst1q_f32
#define tfxWideStorei vst1q_s32
#define tfxWideCasti vreinterpretq_s32_f32
#define tfxWideCast vreinterpretq_f32_s32
#define tfxWideConverti vcvtnq_s32_f32
#define tfxWideConvert vcvtq_f32_s32
#define tfxWideMin vminq_f32
#define tfxWideMax vmaxq_f32
#define tfxWideMini vminq_s32
#define tfxWideMaxi vmaxq_s32
#define tfxWideOri vorrq_s32
#define tfxWideXOri veorq_s32
#define tfxWideAnd(a, b) vreinterpretq_f32_s32(vandq_s32(vreinterpretq_s32_f32(a), vreinterpretq_s32_f32(b)))
#define tfxWideAndi vandq_s32
#define tfxWideAndNoti vbicq_s32
#define tfxWideSetZeroi vdupq_n_s32(0)
#define tfxWideSetZero vdupq_n_f32(0.0f)
#define tfxWideEqualsi vceqq_s32
#define tfxWideEquals vreinterpretq_f32_u32(vceqq_f32)

const float32x4_t tfxWIDEF3_4 = vdupq_n_f32(1.0f / 3.0f);
const float32x4_t tfxWIDEG3_4 = vdupq_n_f32(1.0f / 6.0f);
const float32x4_t tfxWIDEG32_4 = vdupq_n_f32((1.0f / 6.0f) * 2.f);
const float32x4_t tfxWIDEG33_4 = vdupq_n_f32((1.0f / 6.0f) * 3.f);
const int32x4_t tfxWIDEONEi = vdupq_n_s32(1);
const float32x4_t tfxWIDEONE = vdupq_n_f32(1.f);
const float32x4_t tfxWIDE255 = vdupq_n_f32(255.f);
const float32x4_t tfxWIDEZERO = vdupq_n_f32(0.f);
const float32x4_t tfxWIDETHIRTYTWO = vdupq_n_f32(32.f);
const int32x4_t tfxWIDEFF = vdupq_n_s32(0xFF);
const float32x4_t tfxPWIDESIX = vdupq_n_f32(0.6f);

typedef union {
    int32x4_t m;
    int a[4];
} tfxWideArrayi;

typedef union {
    float32x4_t m;
    float a[4];
} tfxWideArray;


#endif

#define tfxWideLookupSet(lookup, index) tfxWideSet( lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )
#define tfxWideLookupSetMember(lookup, member, index) tfxWideSet( lookup[index.a[3]].member, lookup[index.a[2]].member, lookup[index.a[1]].member, lookup[index.a[0]].member )
#define tfxWideLookupSetMemberi(lookup, member, index) tfxWideSeti( lookup[index.a[3]].member, lookup[index.a[2]].member, lookup[index.a[1]].member, lookup[index.a[0]].member )
#define tfxWideLookupSet2(lookup1, lookup2, index1, index2) tfxWideSet( lookup1[index1.a[3]].lookup2[index2.a[3]], lookup1[index1.a[2]].lookup2[index2.a[2]], lookup1[index1.a[1]].lookup2[index2.a[1]], lookup1[index1.a[0]].lookup2[index2.a[0]] )
#define tfxWideLookupSeti(lookup, index) tfxWideSeti( lookup[index.a[3]], lookup[index.a[2]], lookup[index.a[1]], lookup[index.a[0]] )

#endif

#ifdef tfxINTEL
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

#elifdef tfxARM

typedef float32x4_t tfx128;
typedef int32x4_t tfx128i;

typedef union {
    tfx128i m;
    int a[4];
} tfx128iArray;

typedef union {
    tfx128i m;
    tfxU64 a[2];
} tfx128iArray64;

typedef union {
    tfx128 m;
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
    tfx128 j = vdupq_n_f32(1.f); //create vector 1.0f
    tfx128i i = vcvtnq_s32_f32(x);
    tfx128 fi = vreinterpretq_f32_s32(i);
    tfx128 igx = vreinterpretq_f32_s32(vcgtq_f32(fi, x));
    j = vreinterpretq_f32_s32(vandq_s32(vreinterpretq_s32_f32(igx), vreinterpretq_s32_f32(j)));
    return vsubq_f32(fi, j);
}

#endif

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
	tfxBillboarding_align_to_vector = 3,			//Align to vector
	tfxBillboarding_max = 4				
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
	tfx_str512_t ReadLine();
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

//Wide simd versions of tfx_vec2_t/3
struct tfx_wide_vec3_t {
	union {
		struct { tfxWideFloat x, y, z; };
		struct { tfxWideFloat pitch, yaw, roll; };
	};

	tfx_wide_vec3_t() { x = tfxWideSetSingle(0.f); y = tfxWideSetSingle(0.f); z = tfxWideSetSingle(0.f); }
	tfx_wide_vec3_t(tfxWideFloat _x, tfxWideFloat _y, tfxWideFloat _z) { x = _x; y = _y; z = _z; }

	inline tfx_wide_vec3_t operator+(const tfx_wide_vec3_t &v) const { return tfx_wide_vec3_t(tfxWideAdd(x, v.x), tfxWideAdd(y, v.y), tfxWideAdd(z, v.z)); }
	inline tfx_wide_vec3_t operator-(const tfx_wide_vec3_t &v) const { return tfx_wide_vec3_t(tfxWideSub(x, v.x), tfxWideSub(y, v.y), tfxWideSub(z, v.z)); }
	inline tfx_wide_vec3_t operator*(const tfx_wide_vec3_t &v) const { return tfx_wide_vec3_t(tfxWideMul(x, v.x), tfxWideMul(y, v.y), tfxWideMul(z, v.z)); }
	inline tfx_wide_vec3_t operator/(const tfx_wide_vec3_t &v) const { return tfx_wide_vec3_t(tfxWideDiv(x, v.x), tfxWideDiv(y, v.y), tfxWideDiv(z, v.z)); }

	inline tfx_wide_vec3_t operator-() const { return tfx_wide_vec3_t(tfxWideSub(tfxWideSetSingle(0.f), x), tfxWideSub(tfxWideSetSingle(0.f), y), tfxWideSub(tfxWideSetSingle(0.f), z)); }

	inline void operator-=(const tfx_wide_vec3_t &v) { x = tfxWideSub(x, v.x); y = tfxWideSub(y, v.y); z = tfxWideSub(z, v.z); }
	inline void operator+=(const tfx_wide_vec3_t &v) { x = tfxWideAdd(x, v.x); y = tfxWideAdd(y, v.y); z = tfxWideAdd(z, v.z); }
	inline void operator*=(const tfx_wide_vec3_t &v) { x = tfxWideMul(x, v.x); y = tfxWideMul(y, v.y); z = tfxWideMul(z, v.z); }
	inline void operator/=(const tfx_wide_vec3_t &v) { x = tfxWideDiv(x, v.x); y = tfxWideDiv(y, v.y); z = tfxWideDiv(z, v.z); }

	inline tfx_wide_vec3_t operator+(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfx_wide_vec3_t(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v), tfxWideAdd(z, wide_v)); }
	inline tfx_wide_vec3_t operator-(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfx_wide_vec3_t(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v), tfxWideAdd(z, wide_v)); }
	inline tfx_wide_vec3_t operator*(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfx_wide_vec3_t(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v), tfxWideAdd(z, wide_v)); }
	inline tfx_wide_vec3_t operator/(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfx_wide_vec3_t(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v), tfxWideAdd(z, wide_v)); }

	inline tfx_wide_vec3_t operator+(tfxWideFloat v) const { return tfx_wide_vec3_t(tfxWideAdd(x, v), tfxWideAdd(y, v), tfxWideAdd(z, v)); }
	inline tfx_wide_vec3_t operator-(tfxWideFloat v) const { return tfx_wide_vec3_t(tfxWideAdd(x, v), tfxWideAdd(y, v), tfxWideAdd(z, v)); }
	inline tfx_wide_vec3_t operator*(tfxWideFloat v) const { return tfx_wide_vec3_t(tfxWideAdd(x, v), tfxWideAdd(y, v), tfxWideAdd(z, v)); }
	inline tfx_wide_vec3_t operator/(tfxWideFloat v) const { return tfx_wide_vec3_t(tfxWideAdd(x, v), tfxWideAdd(y, v), tfxWideAdd(z, v)); }

	inline void operator*=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideMul(x, wide_v); y = tfxWideMul(y, wide_v); z = tfxWideMul(z, wide_v); }
	inline void operator/=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideDiv(x, wide_v); y = tfxWideDiv(y, wide_v); z = tfxWideDiv(z, wide_v); }
	inline void operator+=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideAdd(x, wide_v); y = tfxWideAdd(y, wide_v); z = tfxWideAdd(z, wide_v); }
	inline void operator-=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideSub(x, wide_v); y = tfxWideSub(y, wide_v); z = tfxWideSub(z, wide_v); }

	inline tfxWideFloat Squared() { return tfxWideAdd(tfxWideMul(x, x), tfxWideAdd(tfxWideMul(y, y), tfxWideMul(z, z))); }
};

struct tfx_wide_vec2_t {
	tfxWideFloat x, y;

	tfx_wide_vec2_t() { x = tfxWideSetSingle(0.f); y = tfxWideSetSingle(0.f); }
	tfx_wide_vec2_t(tfxWideFloat _x, tfxWideFloat _y) { x = _x; y = _y; }

	inline tfx_wide_vec2_t operator+(const tfx_wide_vec2_t &v) const { return tfx_wide_vec2_t(tfxWideAdd(x, v.x), tfxWideAdd(y, v.y)); }
	inline tfx_wide_vec2_t operator-(const tfx_wide_vec2_t &v) const { return tfx_wide_vec2_t(tfxWideSub(x, v.x), tfxWideSub(y, v.y)); }
	inline tfx_wide_vec2_t operator*(const tfx_wide_vec2_t &v) const { return tfx_wide_vec2_t(tfxWideMul(x, v.x), tfxWideMul(y, v.y)); }
	inline tfx_wide_vec2_t operator/(const tfx_wide_vec2_t &v) const { return tfx_wide_vec2_t(tfxWideDiv(x, v.x), tfxWideDiv(y, v.y)); }

	inline tfx_wide_vec2_t operator-() const { return tfx_wide_vec2_t(tfxWideSub(tfxWideSetSingle(0.f), x), tfxWideSub(tfxWideSetSingle(0.f), y)); }

	inline void operator-=(const tfx_wide_vec2_t &v) { x = tfxWideSub(x, v.x); y = tfxWideSub(y, v.y); }
	inline void operator+=(const tfx_wide_vec2_t &v) { x = tfxWideAdd(x, v.x); y = tfxWideAdd(y, v.y); }
	inline void operator*=(const tfx_wide_vec2_t &v) { x = tfxWideMul(x, v.x); y = tfxWideMul(y, v.y); }
	inline void operator/=(const tfx_wide_vec2_t &v) { x = tfxWideDiv(x, v.x); y = tfxWideDiv(y, v.y); }

	inline tfx_wide_vec2_t operator+(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfx_wide_vec2_t(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v)); }
	inline tfx_wide_vec2_t operator-(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfx_wide_vec2_t(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v)); }
	inline tfx_wide_vec2_t operator*(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfx_wide_vec2_t(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v)); }
	inline tfx_wide_vec2_t operator/(float v) const { tfxWideFloat wide_v = tfxWideSetSingle(v); return tfx_wide_vec2_t(tfxWideAdd(x, wide_v), tfxWideAdd(y, wide_v)); }

	inline tfx_wide_vec2_t operator+(tfxWideFloat v) const { return tfx_wide_vec2_t(tfxWideAdd(x, v), tfxWideAdd(y, v)); }
	inline tfx_wide_vec2_t operator-(tfxWideFloat v) const { return tfx_wide_vec2_t(tfxWideAdd(x, v), tfxWideAdd(y, v)); }
	inline tfx_wide_vec2_t operator*(tfxWideFloat v) const { return tfx_wide_vec2_t(tfxWideAdd(x, v), tfxWideAdd(y, v)); }
	inline tfx_wide_vec2_t operator/(tfxWideFloat v) const { return tfx_wide_vec2_t(tfxWideAdd(x, v), tfxWideAdd(y, v)); }

	inline void operator*=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideMul(x, wide_v); y = tfxWideMul(y, wide_v); }
	inline void operator/=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideDiv(x, wide_v); y = tfxWideDiv(y, wide_v); }
	inline void operator+=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideAdd(x, wide_v); y = tfxWideAdd(y, wide_v); }
	inline void operator-=(float v) { tfxWideFloat wide_v = tfxWideSetSingle(v); x = tfxWideSub(x, wide_v); y = tfxWideSub(y, wide_v); }

	inline tfxWideFloat Squared() { return tfxWideAdd(tfxWideMul(x, x), tfxWideMul(y, y)); }
};

//Note, has padding for the sake of alignment on GPU compute shaders
struct tfx_bounding_box_t {
	tfx_vec3_t min_corner; float padding1;
	tfx_vec3_t max_corner; float padding2;
};

inline tfx_wide_vec3_t InterpolateWideVec3(tfxWideFloat &tween, tfx_wide_vec3_t &from, tfx_wide_vec3_t &to) {
	return to * tween + from * (tfxWideSub(tfxWIDEONE, tween));
}

static inline void ScaleVec4xyz(tfx_vec4_t &v, float scalar) {
	v.x *= scalar;
	v.y *= scalar;
	v.z *= scalar;
}

const float tfxONE_DIV_255 = 1 / 255.f;
const float TFXONE_DIV_511 = 1 / 511.f;

struct tfx_rgba8_t {
	union {
		struct {	tfxU32 r : 8; 
					tfxU32 g : 8; 
					tfxU32 b : 8; 
					tfxU32 a : 8; 
		};
		struct {	tfxU32 color; };
	};

	tfx_rgba8_t() { r = g = b = a = 0; }
	tfx_rgba8_t(unsigned char _r, unsigned char _g, unsigned char _b, unsigned char _a) : r(_r), g(_g), b(_b), a(_a) { }
	tfx_rgba8_t(float _r, float _g, float _b, float _a) : r((char)_r), g((char)_g), b((char)_b), a((char)_a) { }
	tfx_rgba8_t(tfxU32 _r, tfxU32 _g, tfxU32 _b, tfxU32 _a) : r((char)_r), g((char)_g), b((char)_b), a((char)_a) { }
	tfx_rgba8_t(int _r, int _g, int _b, int _a) : r((char)_r), g((char)_g), b((char)_b), a((char)_a) { }
	tfx_rgba8_t(tfx_rgba8_t _c, char _a) : r(_c.r), g(_c.g), b(_c.b), a((char)_a) { }
};

struct tfx_rgb_t {
	float r, g, b;
	tfx_rgb_t() { r = g = b = 0.f; }
	tfx_rgb_t(float _r, float _g, float _b) : r(_r), g(_g), b(_b) { }
};

struct tfx_hsv_t {
	float h, s, v;
	tfx_hsv_t() { h = s = v = 0.f; }
	tfx_hsv_t(float _h, float _s, float _v) : h(_h), s(_s), v(_v) { }
};

const tfxWideFloat one_div_511_wide = tfxWideSetSingle(1 / 511.f);
const tfxWideFloat one_div_32k_wide = tfxWideSetSingle(1 / 32767.f);
#define tfxPACKED_Y_NORMAL_3D 0x1FFFF9FF
#define tfxPACKED_Y_NORMAL_2D 32767

struct tfx_rgba_t {
	float r, g, b, a;
	tfx_rgba_t() { r = g = b = a = 1.f; }
	tfx_rgba_t(float _r, float _g, float _b, float _a) : r(_r), g(_g), b(_b), a(_a) { }
	tfx_rgba_t(tfx_rgba8_t c) : r((float)c.r * tfxONE_DIV_255), g((float)c.g * tfxONE_DIV_255), b((float)c.b * tfxONE_DIV_255), a((float)c.a * tfxONE_DIV_255) { }
};

struct tfx_mat4_t {

	tfx_vec4_t v[4];

	inline void Set2(float aa, float ab, float ba, float bb) {
		v[0].c0 = aa; v[0].c1 = ab;
		v[1].c0 = ba; v[1].c1 = bb;
	}

} TFX_ALIGN_AFFIX(16);

struct tfx_wide_mat4_t {
	float x[4];
	float y[4];
	float z[4];
	float w[4];
} TFX_ALIGN_AFFIX(16);;

struct tfx_mat3_t {

	tfx_vec3_t v[3];

	inline tfx_vec3_t operator*(const tfx_vec3_t &vec) const {
		return tfx_vec3_t(
			v[0].x * vec.x + v[1].x * vec.y + v[2].x * vec.z,
			v[0].y * vec.x + v[1].y * vec.y + v[2].y * vec.z,
			v[0].z * vec.x + v[1].z * vec.y + v[2].z * vec.z);
	}

};

//Very simple 2D Matix
struct tfx_mat2_t {

	float aa, ab, ba, bb;

	tfx_mat2_t() :aa(1.f), ab(0.f), ba(0.f), bb(1.f) {}

	void Set(float _aa = 1.f, float _ab = 1.f, float _ba = 1.f, float _bb = 1.f) {
		aa = _aa;
		ab = _ab;
		ba = _ba;
		bb = _bb;
	}

	void Transpose() {
		float abt = ab;
		ab = ba;
		ba = abt;
	}

	void Scale(float s) {
		aa *= s;
		ab *= s;
		ba *= s;
		bb *= s;
	}

	tfx_mat2_t Transform(const tfx_mat2_t &m) {
		tfx_mat2_t r;
		r.aa = aa * m.aa + ab * m.ba; r.ab = aa * m.ab + ab * m.bb;
		r.ba = ba * m.aa + bb * m.ba; r.bb = ba * m.ab + bb * m.bb;
		return r;
	}

	tfx_mat2_t Transform(const tfx_mat4_t &m) {
		tfx_mat2_t r;
		r.aa = aa * m.v[0].x + ab * m.v[1].x; r.ab = aa * m.v[0].y + ab * m.v[1].y;
		r.ba = ba * m.v[0].x + bb * m.v[1].x; r.bb = ba * m.v[0].y + bb * m.v[1].y;
		return r;
	}

	tfx_vec2_t TransformVector(const tfx_vec2_t v) {
		tfx_vec2_t tv = tfx_vec2_t(0.f, 0.f);
		tv.x = v.x * aa + v.y * ba;
		tv.y = v.x * ab + v.y * bb;
		return tv;
	}

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
const tfx128 tfxF3_4 = _mm_set_ps1(1.0f / 3.0f);
const tfx128 tfxF2_4 = _mm_set_ps1(.366025403f);
const tfx128 tfxG2_4 = _mm_set_ps1(0.211324865f);
const tfx128 tfxG2_4x2 = _mm_set_ps1(0.42264973f);
const tfx128 tfxG3_4 = _mm_set_ps1(1.0f / 6.0f);
const tfx128 tfxG32_4 = _mm_set_ps1((1.0f / 6.0f) * 2.f);
const tfx128 tfxG33_4 = _mm_set_ps1((1.0f / 6.0f) * 3.f);
const tfx128i tfxONE = _mm_set1_epi32(1);
const tfx128 tfxONEF = _mm_set1_ps(1.f);
const tfx128 tfxZERO = _mm_set1_ps(0.f);
const tfx128 tfxTHIRTYTWO = _mm_set1_ps(32.f);
const tfx128i tfxFF = _mm_set1_epi32(0xFF);
const tfx128 tfxPSIX = _mm_set_ps1(0.6f);
#elifdef tfxARM
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
const uint8_t perm[] =
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

static const uint8_t permMOD12[] =
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
tfx128Array tfxNoise4_2d(const tfx128 &x4, const tfx128 &y4);
tfx128Array tfxNoise4_3d(const tfx128 &x4, const tfx128 &y4, const tfx128 &z4);

//-----------------------------------------------------------
//Section: Profiling
//-----------------------------------------------------------

struct tfx_profile_stats_t {
	tfxU64 cycle_high;
	tfxU64 cycle_low;
	tfxU64 cycle_average;
	tfxU64 time_high;
	tfxU64 time_low;
	tfxU64 time_average;
	tfxU32 hit_count;
};

struct tfx_profile_snapshot_t {
	tfxU32 hit_count = 0;
	tfxU64 run_time = 0;
	tfxU64 cycle_count = 0;
};

struct tfx_profile_t {
	const char *name;
	tfx_profile_snapshot_t snapshots[tfxPROFILER_SAMPLES];
};

extern const tfxU32 tfxPROFILE_COUNT;
extern tfxU32 tfxCurrentSnapshot;
extern tfx_profile_t tfxProfileArray[];

struct tfx_profile_tag_t {
	tfx_profile_t *profile;
	tfx_profile_snapshot_t *snapshot;
	tfxU64 start_cycles;
	tfxU64 start_time;

	tfx_profile_tag_t(tfxU32 id, const char *name);

	~tfx_profile_tag_t() {
		tfx_AtomicAdd64(&snapshot->run_time, tfx_Microsecs() - start_time);
		tfx_AtomicAdd64(&snapshot->cycle_count, (__rdtsc() - start_cycles));
	}

};

#ifdef tfxENABLE_PROFILING 
#define tfxPROFILE tfx_profile_tag_t tfx_tag((tfxU32)__COUNTER__, __FUNCTION__)
#else
#define tfxPROFILE __COUNTER__
#endif

//-----------------------------------------------------------
//Section: File_IO
//-----------------------------------------------------------

const tfxU32 tfxMAGIC_NUMBER = 559433300;				//'!XFT'
const tfxU32 tfxMAGIC_NUMBER_INVENTORY = 559304265;		//'!VNI'
const tfxU32 tfxFILE_VERSION = 2;

//Basic package manager used for reading/writing effects files
struct tfx_package_header_t {
	tfxU32 magic_number;						//Magic number to confirm file format
	tfxU32 file_version;						//The version of the file
	tfxU32 flags;								//Any state_flags for the file
	tfxU32 reserved0;							//Reserved for future if needed
	tfxU64 offset_to_inventory;					//Memory offset for the inventory of files
	tfxU64 user_data1;							//Any data you might find useful
	tfxU64 user_data2;							//Any data you might find useful
	tfxU64 reserved3;							//More reserved space
	tfxU64 reserved4;							//More reserved space
	tfxU64 reserved5;							//More reserved space
};

struct tfx_package_entry_info_t {
	tfx_str_t file_name;						//The name of the file stored in the package
	tfxU64 offset_from_start_of_file = 0;		//Offset from the start of the file to where the file is located
	tfxU64 file_size = 0;						//The size of the file
	tfx_stream_t data;							//The file data

	void FreeData();
};

struct tfx_package_inventory_t {
	tfxU32 magic_number;						//Magic number to confirm format of the Inventory
	tfxU32 entry_count;							//Number of files in the inventory
	tfx_storage_map_t<tfx_package_entry_info_t> entries;		//The inventory list

	tfx_package_inventory_t() :
		entries("Inventory Map", "Inventory Data"),
		magic_number(0),
		entry_count(0)
	{}
};

struct tfx_package_t {
	tfx_str_t file_path;
	tfx_package_header_t header;
	tfx_package_inventory_t inventory;
	tfxU64 file_size = 0;						//The total file size of the package, should match file size on disk
	tfx_stream_t file_data;						//Dump of the data from the package file on disk

	~tfx_package_t();

};

//------------------------------------------------------------
//Section: Struct_Types
//------------------------------------------------------------

struct tfx_face_t {
	int v[3];
};
extern tfx_vector_t<tfx_vec3_t> tfxIcospherePoints[6];

struct tfx_attribute_node_t {
	float frame;
	float value;

	tfx_vec2_t left;
	tfx_vec2_t right;

	tfxAttributeNodeFlags flags;
	tfxU32 index;

	tfx_attribute_node_t() : frame(0.f), value(0.f), flags(0), index(0) { }
	inline bool operator==(const tfx_attribute_node_t& n) { return n.frame == frame && n.value == value; }
};

struct tfx_random_t {
	tfxU64 seeds[2];
};

struct tfx_graph_lookup_t {
	tfx_vector_t<float> values;
	tfxU32 last_frame;
	float life;

	tfx_graph_lookup_t() : last_frame(0), life(0) {}
};

struct tfx_graph_id_t {
	tfx_graph_category category;
	tfx_graph_type type = tfxGraphMaxIndex;
	tfxU32 graph_id = 0;
	tfxU32 node_id = 0;
};

struct tfx_graph_lookup_index_t {
	tfxU32 start_index;
	tfxU32 length;
	float max_life;
	float padding1;
};

//This struct is used to store indexing data in order to index into large lists containing either the node data of graphs
//or the lookup data of compiled graphs. This is so that we can upload that data into a buffer on the GPU to get the particles
//updating in a compute shader.
struct tfx_effect_lookup_data_t {
	tfx_graph_lookup_index_t overtime_velocity;
	tfx_graph_lookup_index_t overtime_width;
	tfx_graph_lookup_index_t overtime_height;
	tfx_graph_lookup_index_t overtime_weight;
	tfx_graph_lookup_index_t overtime_spin;
	tfx_graph_lookup_index_t overtime_stretch;
	tfx_graph_lookup_index_t overtime_red;
	tfx_graph_lookup_index_t overtime_green;
	tfx_graph_lookup_index_t overtime_blue;
	tfx_graph_lookup_index_t overtime_opacity;
	tfx_graph_lookup_index_t overtime_velocity_turbulance;
	tfx_graph_lookup_index_t overtime_direction_turbulance;
	tfx_graph_lookup_index_t overtime_velocity_adjuster;
	tfx_graph_lookup_index_t overtime_intensity;
	tfx_graph_lookup_index_t overtime_direction;
	tfx_graph_lookup_index_t overtime_noise_resolution;
};

struct tfx_graph_t {
	//The ratio to transalte graph frame/value to grid x/y coords on a graph editor

	tfx_vec2_t min;
	tfx_vec2_t max;
	tfx_graph_preset graph_preset;
	tfx_graph_type type;
	tfx_effect_emitter_t *effector;
	tfx_bucket_array_t<tfx_attribute_node_t> nodes;
	tfx_graph_lookup_t lookup;
	tfxU32 index;
	float gamma;

	tfx_graph_t();
	tfx_graph_t(tfxU32 bucket_size);
	~tfx_graph_t();

};

//The following structs group graphs together under the attribute categories Global, Transform, Properties, Base, Variation and Overtime
struct tfx_global_attributes_t {
	tfx_graph_t life;
	tfx_graph_t amount;
	tfx_graph_t velocity;
	tfx_graph_t width;
	tfx_graph_t height;
	tfx_graph_t weight;
	tfx_graph_t spin;
	tfx_graph_t stretch;
	tfx_graph_t overal_scale;
	tfx_graph_t intensity;
	tfx_graph_t splatter;
	tfx_graph_t emitter_width;
	tfx_graph_t emitter_height;
	tfx_graph_t emitter_depth;
};

struct tfx_transform_attributes_t {
	tfx_graph_t roll;
	tfx_graph_t pitch;
	tfx_graph_t yaw;
	tfx_graph_t translation_x;
	tfx_graph_t translation_y;
	tfx_graph_t translation_z;
};

struct tfx_property_attributes_t {
	tfx_graph_t emission_pitch;
	tfx_graph_t emission_yaw;
	tfx_graph_t emission_range;
	tfx_graph_t splatter;
	tfx_graph_t emitter_width;
	tfx_graph_t emitter_height;
	tfx_graph_t emitter_depth;
	tfx_graph_t arc_size;
	tfx_graph_t arc_offset;
};

struct tfx_base_attributes_t {
	tfx_graph_t life;
	tfx_graph_t amount;
	tfx_graph_t velocity;
	tfx_graph_t width;
	tfx_graph_t height;
	tfx_graph_t weight;
	tfx_graph_t spin;
	tfx_graph_t noise_offset;
};

struct tfx_variation_attributes_t {
	tfx_graph_t life;
	tfx_graph_t amount;
	tfx_graph_t velocity;
	tfx_graph_t width;
	tfx_graph_t height;
	tfx_graph_t weight;
	tfx_graph_t spin;
	tfx_graph_t noise_offset;
	tfx_graph_t noise_resolution;
};

struct tfx_overtime_attributes_t {
	tfx_graph_t velocity;
	tfx_graph_t width;
	tfx_graph_t height;
	tfx_graph_t weight;
	tfx_graph_t spin;
	tfx_graph_t stretch;
	tfx_graph_t red;
	tfx_graph_t green;
	tfx_graph_t blue;
	tfx_graph_t blendfactor;
	tfx_graph_t velocity_turbulance;
	tfx_graph_t direction_turbulance;
	tfx_graph_t velocity_adjuster;
	tfx_graph_t intensity;
	tfx_graph_t direction;
	tfx_graph_t noise_resolution;
};

struct tfx_emitter_attributes_t {
	tfx_property_attributes_t properties;
	tfx_base_attributes_t base;
	tfx_variation_attributes_t variation;
	tfx_overtime_attributes_t overtime;
};

static float(*lookup_overtime_callback)(tfx_graph_t *graph, float age, float lifetime);
static float(*lookup_callback)(tfx_graph_t *graph, float age);
static float(*lookup_random_callback)(tfx_graph_t *graph, float age, tfx_random_t *random);

struct tfx_shape_data_t {
	char name[64];
	tfxU32 frame_count = 0;
	tfxU32 width = 0;
	tfxU32 height = 0;
	tfxU32 shape_index = 0;
	tfxKey image_hash = 0;
	int import_filter = 0;
};

struct tfx_base_t {
	tfx_vec2_t size;
	float velocity;
	float spin;
	float weight;
};

struct tfx_camera_settings_t {
	tfx_vec3_t camera_position;
	float camera_pitch;
	float camera_yaw;
	float camera_fov;
	float camera_floor_height;
	float camera_isometric_scale;
	bool camera_isometric;
	bool camera_hide_floor;
};

struct tfx_preview_camera_settings_t {
	tfx_camera_settings_t camera_settings;
	float effect_z_offset;
	float camera_speed;
	bool attach_effect_to_camera;
};

//this probably only needs to be in the editor, no use for it in the library? Maybe in the future as an alternative way to play back effects...
struct tfx_sprite_sheet_settings_t {
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
	float effect_z_offset = 5.f;
	tfx_export_color_options color_option;
	tfx_export_options export_option;
	tfx_camera_settings_t camera_settings;
	tfx_camera_settings_t camera_settings_orthographic;
};

//This struct has the settings for recording sprite data frames so that they can be played back as an alternative to dynamic particle updating
struct tfx_sprite_data_settings_t {
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
};

//------------------------------------------------------------

//API structs you can access in various ways to update and render effects in realtime

//Image data for particle shapes. This is passed into your custom ShapeLoader function for loading image textures into whatever renderer you're using
struct tfx_image_data_t {
	//This can be a ptr to the image texture for rendering. You must assign this in your ShapeLoader function
	void *ptr;
	//Index of the image, deprecated, image hash should be used now instead.
	tfxU32 shape_index;
	//Name of the image
	tfx_str_t name;
	//A hash of the image data for a unique and which can also be used to see if an image has already been loaded
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

		tfx_image_data_t() :
		image_index(0),
		shape_index(0),
		ptr(nullptr),
		animation_frames(1.f),
		image_hash(0),
		max_radius(0),
		import_filter(0),
		compute_shape_index(0)
	{ }
};

struct tfx_emitter_properties_t {
	//Angle added to the rotation of the particle when spawned or random angle range if angle setting is set to tfx_random_t
	tfx_vec3_t angle_offsets;
	//When aligning the billboard along a vector, you can set the type of vector that it aligns with
	tfx_vector_align_type vector_align_type;
	//Point, area, ellipse emitter etc.
	tfx_emission_type emission_type;
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
	//Offset of emitters
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

	tfx_emitter_properties_t() { memset(this, 0, sizeof(tfx_emitter_properties_t)); }
};

//Stores the most recent parent effect (with global attributes) spawn control values to be applied to sub emitters.
struct tfx_parent_spawn_controls_t {
	float life;
	float size_x;
	float size_y;
	float velocity;
	float spin;
	float intensity;
	float splatter;
	float weight;
};

struct tfx_effect_emitter_info_t {
	//Name of the effect
	tfx_str64_t name;
	//Every effect and emitter in the library gets a unique id
	tfxU32 uid;
	//The max_radius of the emitter, taking into account all the particles that have spawned and active (editor only)
	float max_radius;
	//List of sub_effects ( effects contain emitters, emitters contain sub effects )
	tfx_vector_t<tfx_effect_emitter_t> sub_effectors;
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
	//The estimated maximum time that the sub emitter might last for, taking into account the parent particle lifetime
	float max_sub_emitter_life;
	//The maximum amount of particles that this effect can spawn (root effects and emitters only)
	tfxU32 max_particles[tfxLAYERS];
	tfxU32 max_sub_emitters;

	tfx_effect_emitter_info_t() :
		lookup_node_index(0),
		lookup_value_index(0),
		sprite_data_settings_index(0),
		uid(0),
		max_particles{2500, 2500, 2500, 2500},
		max_radius(0),
		sprite_sheet_settings_index(0),
		preview_camera_settings(0),
		max_sub_emitters(0),
		max_life(0),
		max_sub_emitter_life(0.f),
		sub_effectors(tfxCONSTRUCTOR_VEC_INIT("sub_effectors"))
	{
		for (int i = 0; i != tfxLAYERS; ++i) {
			max_particles[i] = 0;
		}
	}
};

//This is a struct that stores an emitter state that is currently active in a particle manager.
struct tfx_emitter_state_t {
	//State data
	float frame;
	float age;
	float highest_particle_age;
	float delay_spawning;
	float timeout_counter;
	float timeout;
	tfx_vec3_t handle;
	tfxEmitterPropertyFlags property_flags;
	float loop_length;
	//Position, scale and rotation values
	tfx_vec3_t local_position;
	tfx_vec3_t world_position;
	tfx_vec3_t captured_position;
	tfx_vec3_t world_rotations;
	//Todo: save space and use a quaternion here... maybe
	tfx_mat4_t matrix;
	tfx_vec2_t image_handle;
	tfx_bounding_box_t bounding_box;

	float amount_remainder;
	float spawn_quantity;
	float qty_step_size;

	tfxU32 emitter_attributes;
	tfxU32 transform_attributes;
	tfxU32 overtime_attributes;

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
	float image_frame_rate;
	float end_frame;
	tfx_vec3_t grid_coords;
	tfx_vec3_t grid_direction;
	tfx_vec3_t emitter_size;
	float emission_alternator;
	tfxEmitterStateFlags state_flags;
	tfx_vec2_t image_size;
	tfx_vec3_t angle_offsets;
} TFX_ALIGN_AFFIX(16);

//This is a struct that stores an effect state that is currently active in a particle manager.
struct tfx_effect_state_t {
	tfx_mat4_t matrix;
	//State data
	float frame;
	float age;
	float highest_particle_age;
	float timeout_counter;
	float timeout;
	tfx_vec3_t handle;
	tfxEmitterPropertyFlags property_flags;
	float loop_length;
	//Position, scale and rotation values
	tfx_vec3_t translation;
	tfx_vec3_t local_position;
	tfx_vec3_t world_position;
	tfx_vec3_t captured_position;
	tfx_vec3_t local_rotations;
	tfx_vec3_t world_rotations;
	//Todo: save space and use a quaternion here?
	tfx_bounding_box_t bounding_box;

	tfxU32 global_attributes;
	tfxU32 transform_attributes;

	tfxU32 properties_index;
	tfxU32 info_index;
	tfxU32 parent_particle_index;
	tfx_library_t *library;

	//Spawn controls
	tfx_parent_spawn_controls_t spawn_controls;
	tfx_vec3_t emitter_size;
	float stretch;
	float overal_scale;
	float noise_base_offset;
	tfxEmitterStateFlags state_flags;

	//User Data
	void *user_data;
	void(*update_callback)(tfx_particle_manager_t *pm, tfxEffectID effect_index);
} TFX_ALIGN_AFFIX(16);

//An tfx_effect_emitter_t can either be an effect which stores emitters and global graphs for affecting all the attributes in the emitters
//Or it can be an emitter which spawns all of the particles. Effectors are stored in the particle manager effects list buffer.
//This is only for library storage, when using to update each frame this is copied to tfx_effect_state_t and tfx_emitter_state_t for realtime updates
//suited for realtime use.
struct tfx_effect_emitter_t {
	//Required for frame by frame updating
	//The current state of the effect/emitter used in the editor only at this point
	tfxEmitterStateFlags state_flags;
	tfxEmitterPropertyFlags property_flags;
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
	//Pointer to the immediate parent
	tfx_effect_emitter_t *parent;
	//State state_flags for emitters and effects
	tfxEffectPropertyFlags effect_flags;
	//When not using insert sort to guarantee particle order, sort passes offers a more lax way of ordering particles over a number of frames.
	//The more passes the more quickly ordered the particles will be but at a higher cost
	tfxU32 sort_passes;
	//Custom user data, can be accessed in callback functions
	void *user_data;
	void(*update_callback)(tfx_particle_manager_t *pm, tfxEffectID effect_index);

	tfxU32 buffer_index;

	//Indexes into library storage
	tfxU32 info_index;
	tfxU32 property_index;
	tfxU32 pm_index;

	tfx_effect_emitter_t() :
		buffer_index(0),
		path_hash(0),
		pm_index(0),
		parent(nullptr),
		user_data(nullptr),
		update_callback(nullptr),
		effect_flags(tfxEffectPropertyFlags_none),
		sort_passes(1),
		info_index(tfxINVALID),
		property_index(tfxINVALID),
		global(tfxINVALID),
		emitter_attributes(tfxINVALID),
		transform_attributes(tfxINVALID),
		property_flags(tfxEmitterPropertyFlags_image_handle_auto_center |
			tfxEmitterPropertyFlags_grid_spawn_clockwise |
			tfxEmitterPropertyFlags_emitter_handle_auto_center |
			tfxEmitterPropertyFlags_global_uniform_size |
			tfxEmitterPropertyFlags_base_uniform_size |
			tfxEmitterPropertyFlags_lifetime_uniform_size),
		state_flags(0)
	{ }
	~tfx_effect_emitter_t();

};

struct tfx_compute_sprite_t {	//64 bytes
	tfx_vec4_t bounds;				//the min/max x,y coordinates of the image being drawn
	tfx_vec4_t uv;					//The UV coords of the image in the texture
	tfx_vec4_t scale_rotation;		//Scale and rotation (x, y = scale, z = rotation, w = multiply blend factor)
	tfx_vec2_t position;			//The position of the sprite
	tfx_rgba8_t color;				//The color tint of the sprite
	tfxU32 parameters;	//4 extra parameters packed into a tfxU32: blend_mode, image layer index, shader function index, blend type
};

struct tfx_depth_index_t {
	tfxParticleID particle_id;
	float depth;
};

struct tfx_unique_sprite_id_t {
	tfxU32 uid;
	tfxU32 age;
};

//These all point into a tfx_soa_buffer_t, initialised with InitParticleSoA. Current Bandwidth: 108 bytes
struct tfx_particle_soa_t {
	tfxU32 *uid;		//Only used for recording sprite data
	tfxU32 *parent_index;
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
	float *local_rotations_x;
	float *local_rotations_y;
	float *local_rotations_z;
	tfxU32 *velocity_normal;
	tfxU32 *depth_index;
	float *base_weight;
	float *base_velocity;
	float *base_spin;
	float *base_size_x;
	float *base_size_y;
	float *noise_offset;
	float *noise_resolution;
	tfx_rgba8_t *color;
	float *image_frame;
	tfxU32 *single_loop_count;
};

struct tfx_sprite_transform2d_t {
	tfx_vec2_t position;					//The position of the sprite, x, y - world, z, w = captured for interpolating
	tfx_vec2_t scale;						//Scale
	float rotation;
};

struct tfx_sprite_transform3d_t {
	tfx_vec3_t position;					//The position of the sprite, x, y - world, z, w = captured for interpolating
	tfx_vec3_t rotations;					//Rotations of the sprite
	tfx_vec2_t scale;						//Scale
};

//When exporting effects as sprite data each frame gets frame meta containing information about the frame such as bounding box and sprite count/offset into the buffer
struct tfx_frame_meta_t {
	tfxU32 index_offset[tfxLAYERS];		//All sprite data is contained in a single buffer and this is the offset to the first sprite in the range
	tfxU32 sprite_count[tfxLAYERS];		//The number of sprites in the frame for each layer
	tfxU32 total_sprites;				//The total number of sprites for all layers in the frame
	tfx_vec3_t min_corner;				//Bounding box min corner
	tfx_vec3_t max_corner;				//Bounding box max corner. The bounding box can be used to decide if this frame needs to be drawn
	tfx_vec3_t bb_center_point;			//The center point of the bounding box. For the fastest checking against a viewing frustum, you can combine this with radius
	float radius;						//The radius of the bounding box
};

//This struct of arrays is used for both 2d and 3d sprites, but obviously the transform_3d data is either 2d or 3d depending on which effects you're using in the particle manager.
//InitSprite3dSoA is called to initialise 3d sprites and InitSprite2dArray for 2d sprites. This is all managed internally by the particle manager. It's convenient to have both 2d and
//3d in one struct like this as it makes it a lot easier to use the same control functions where we can. Also note that stretch and alignment for 3d sprites are packed into
//stretch_alignment_x and alignment_yz as 16bit floats. 2d uses the float float for stretch and packs xy alignment into alignment_yz
struct tfx_sprite_soa_t {						//3d takes 56 bytes of bandwidth, 2d takes 40 bytes of bandwidth
	tfxU32 *property_indexes;					//The image frame of animation index packed with alignment option flag and property_index
	tfxU32 *captured_index;						//The index of the sprite in the previous frame so that it can be looked up and interpolated with
	tfx_unique_sprite_id_t *uid;				//Unique particle id of the sprite, only used when recording sprite data
	tfx_sprite_transform3d_t *transform_3d;		//Transform data for 3d sprites
	tfx_sprite_transform2d_t *transform_2d;		//Transform data for 2d sprites
	tfx_rgba8_t *color;							//The color tint of the sprite and blend factor in alpha channel
	float *intensity;							//The multiplier for the sprite color
	float *stretch;								//Multiplier for how much the particle is stretched in the shader
	tfxU32 *alignment;							//The alignment of the particle. 2 16bit floats for 2d and 3 8bit floats for 3d
};

enum tfxSpriteBufferMode {
	tfxSpriteBufferMode_2d,
	tfxSpriteBufferMode_3d,
	tfxSpriteBufferMode_both,
};

//These structs are for animation sprite data that you can upload to the gpu
struct alignas(16) tfx_sprite_data3d_t {	//60 bytes aligning to 64
	tfx_vec3_t position;
	float lerp_offset;
	tfx_vec3_t rotations;
	float stretch;
	tfx_vec2_t scale;
	tfxU32 property_indexes;
	tfxU32 captured_index;
	tfxU32 alignment;
	tfx_rgba8_t color;
	float intensity;
	//Free space for extra 4 bytes if needed
};

struct alignas(16) tfx_sprite_data2d_t {	//48 bytes
	tfx_vec2_t position;
	tfx_vec2_t scale;
	float rotation;
	tfxU32 property_indexes;
	tfxU32 captured_index;
	tfxU32 alignment;
	tfx_rgba8_t color;
	float lerp_offset;
	float stretch;
	float intensity;
};

//Animation sprite data that is used on the cpu to bake the data
struct tfx_sprite_data_soa_t {	//64 bytes or 60 after uid is removed as it's only needed for compressing the sprite data down to size.
	tfxU32 *property_indexes;	//The image frame of animation index packed with alignment option flag and property_index
	tfxU32 *captured_index;
	tfx_unique_sprite_id_t *uid;
	float *lerp_offset;
	tfx_sprite_transform3d_t *transform_3d;
	tfx_sprite_transform2d_t *transform_2d;
	tfx_rgba8_t *color;			//The color tint of the sprite and blend factor in a
	float *intensity;
	float *stretch;
	tfxU32 *alignment;			//normalised alignment vector 3 floats packed into 10bits each with 2 bits left over
};

struct tfx_wide_lerp_transform_result_t {
	float position[3];
	float rotations[3];
	float scale[2];
};

struct tfx_wide_lerp_data_result_t {
	float stretch;
	float intensity;
	float color[4];
	float padding[2];
};

struct tfx_sprite_data_metrics_t {
	tfx_str64_t name;
	tfxKey path_hash;
	tfxU32 start_offset;	//Only applies to animation manager
	tfxU32 frames_after_compression;
	tfxU32 real_frames;
	tfxU32 frame_count;
	float animation_length_in_time;
	tfxU32 total_sprites;
	tfxU32 total_memory_for_sprites;
	tfx_vector_t<tfx_frame_meta_t> frame_meta;
	tfxAnimationManagerFlags flags;
	tfxAnimationFlags animation_flags;
};

struct tfx_sprite_data_t {
	float frame_compression;
	tfx_sprite_data_metrics_t normal;
	tfx_sprite_data_metrics_t compressed;
	tfx_soa_buffer_t real_time_sprites_buffer;
	tfx_sprite_data_soa_t real_time_sprites;
	tfx_soa_buffer_t compressed_sprites_buffer;
	tfx_sprite_data_soa_t compressed_sprites;
};

struct tfx_compute_fx_global_state_t {
	tfxU32 start_index = 0;
	tfxU32 current_length = 0;
	tfxU32 max_index = 0;
	tfxU32 end_index = 0;
};

struct tfx_compute_controller_t {
	tfx_vec2_t position;
	float line_length;
	float angle_offset;
	tfx_vec4_t scale_rotation;				//Scale and rotation (x, y = scale, z = rotation, w = velocity_adjuster)
	float end_frame;
	tfxU32 normalised_values;		//Contains normalized values which are generally either 0 or 255, normalised in the shader to 0 and 1 (except opacity): age_rate, line_negator, spin_negator, position_negator, opacity
	tfxParticleControlFlags flags;
	tfxU32 image_data_index;		//index into the shape buffer on the gpu. CopyComputeShapeData must be called to prepare the data.
	tfx_vec2_t image_handle;
	tfx_vec2_t emitter_handle;
	float noise_offset;
	float stretch;
	float frame_rate;
	float noise_resolution;
};

struct tfx_compute_particle_t {
	tfx_vec2_t local_position;
	tfx_vec2_t base_size;

	float base_velocity = 1;
	float base_spin = 1;
	float base_weight = 1;

	float age = 1;							//The age of the particle, used by the controller to look up the current state on the graphs
	float max_age = 1;						//max age before the particle expires
	float emission_angle = 1;				//Emission angle of the particle at spawn time

	float noise_offset = 1;					//The random velocity added each frame
	float noise_resolution = 1;				//The random velocity added each frame
	float image_frame = 0;
	tfxU32 control_slot_and_layer;	//index to the controller, and also stores the layer in the particle manager that the particle is on (layer << 3)
	float local_rotation;
};

struct alignas(16) tfx_gpu_image_data_t {
	tfx_vec4_t uv;
	tfxU32 uv_xy;
	tfxU32 uv_zw;
	tfx_vec2_t image_size;
	tfxU32 texture_array_index;
	float animation_frames;
#ifdef tfxCUSTOM_GPU_IMAGE_DATA
	//add addition image data if needed
#endif
};

struct tfx_gpu_shapes_t {
	tfx_vector_t<tfx_gpu_image_data_t> list;
};

//Struct to contain a static state of a particle in a frame of animation. Used in the editor for recording frames of animation so probably not needed here really!
struct tfx_particle_frame_t {
	tfxU32 property_indexes;	//The image frame of animation index packed with alignment option flag and property_index
	tfxU32 captured_index;
	tfx_sprite_transform3d_t transform;
	tfx_rgba8_t color;				//The color tint of the sprite and blend factor in a
	float intensity;
	float depth;
	tfxU32 alignment;		
	float stretch;
};

struct tfx_spawn_work_entry_t {
	tfx_particle_manager_t *pm;
	tfx_emitter_properties_t *properties;
	tfx_parent_spawn_controls_t *parent_spawn_controls;
	tfxU32 emitter_index;
	tfxU32 parent_index;
	tfx_emission_type emission_type;
	tfxEmitterPropertyFlags property_flags;
	tfxEmitterPropertyFlags parent_property_flags;
	tfx_particle_soa_t *particle_data;
	tfx_vector_t<tfx_effect_emitter_t> *sub_effects;
	tfxU32 seed;
	float tween;
	tfxU32 max_spawn_count;
	tfxU32 amount_to_spawn = 0;
	tfxU32 spawn_start_index;
	tfxU32 next_buffer;
	int depth;
	float qty_step_size;
	float highest_particle_age;
	float overal_scale;
};

struct tfx_control_work_entry_t {
	tfxU32 start_index;
	tfxU32 end_index;
	tfxU32 wide_end_index;
	tfxU32 start_diff;
	tfxU32 sprites_index;
	tfxU32 sprite_buffer_end_index;
	tfxU32 emitter_index;
	tfx_particle_manager_t *pm;
	tfx_overtime_attributes_t *graphs;
	tfxU32 layer;
	tfx_emitter_properties_t *properties;
	tfx_sprite_soa_t *sprites;
	float overal_scale;
	float stretch;
	float intensity;
};

struct tfx_particle_age_work_entry_t {
	tfxU32 start_index;
	tfxU32 emitter_index;
	tfxU32 wide_end_index;
	tfxU32 start_diff;
	tfx_emitter_properties_t *properties;
	tfx_particle_manager_t *pm;
};

struct tfx_sort_work_entry_t {
	tfx_bucket_array_t<tfx_particle_soa_t> *bank;
	tfx_vector_t<tfx_depth_index_t> *depth_indexes;
};

struct tfx_compress_work_entry_t {
	tfx_sprite_data_t *sprite_data;
	tfxU32 frame;
};

struct tfx_effect_data_t {
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
};

//An anim instance is used to let the gpu know where to draw an animation with sprite data. 48 bytes
struct tfx_animation_instance_t {
	tfx_vec3_t position;				//position that the instance should be played at
	float scale;						//Scales the overal size of the animation
	tfxU32 sprite_count;				//The number of sprites to be drawn
	tfxU32 frame_count;					//The number of frames in the animation
	tfxU32 offset_into_sprite_data;		//The starting ofset in the buffer that contains all the sprite data
	tfxU32 info_index;					//Index into the effect_animation_info storage map to get at the frame meta
	float current_time;					//Current point of time in the animation
	float animation_length_in_time;		//Total time that the animation lasts for
	float tween;						//The point time within the frame (0..1)
	tfxAnimationInstanceFlags flags;	//Flags associated with the instance
};

struct tfx_animation_buffer_metrics_t {
	size_t sprite_data_size;
	tfxU32 offsets_size;
	tfxU32 instances_size;
	size_t offsets_size_in_bytes;
	size_t instances_size_in_bytes;
	tfxU32 total_sprites_to_draw;

	tfx_animation_buffer_metrics_t() : sprite_data_size(0), offsets_size(0), instances_size(0), total_sprites_to_draw(0), instances_size_in_bytes(0), offsets_size_in_bytes(0) {}
};

struct alignas(16) tfx_animation_emitter_properties_t {
	tfx_vec2_t handle;
	tfxU32 flags;
	tfxU32 start_frame_index;
	float animation_frames;
	void *image_ptr;		//Note: not needed on the GPU, only used if you interpolate and render on the cpu for whatever reason
};

//Use the animation manager to control playing of pre-recorded effects
struct tfx_animation_manager_t {
	//All of the sprite data for all the animations that you might want to play on the GPU.
	//This could be deleted once it's uploaded to the GPU
	//An animation manager can only be used for either 2d or 3d not both
	tfx_vector_t<tfx_sprite_data3d_t> sprite_data_3d;
	tfx_vector_t<tfx_sprite_data2d_t> sprite_data_2d;
	//List of active instances that are currently playing
	tfx_vector_t<tfx_animation_instance_t> instances;
	//List of instances in use. These index into the instances list above
	tfx_vector_t<tfxU32> instances_in_use[2];
	//Flips between 1 and 0 each frame to be used when accessing instances_in_use
	tfxU32 current_in_use_buffer;
	//List of free instance indexes
	tfx_vector_t<tfxU32> free_instances;
	//List of indexes into the instances list that will actually be sent to the GPU next frame
	//Any instances deemed not in view can be culled for example by not adding them to the queue
	tfx_vector_t<tfx_animation_instance_t> render_queue;
	//The compute shader needs to know when to switch from one animation instance to another as it
	//progresses through all the sprites that need to be rendered. So this array of offsets tells
	//it when to do this. So 0 will always be in the first slot, then the second element will be 
	//the total number of sprites to be drawn for the first animation instance and so on. When the
	//global index is more than or equal to the next element then we start on the next animation
	//instance and draw those sprites.
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
	//This struct contains the size of the buffers that need to be uploaded to the GPU. Offsets and 
	//animation instances need to be uploaded every frame, but the sprite data only once before you
	//start drawing anything
	tfx_animation_buffer_metrics_t buffer_metrics;
	//Bit flag field
	tfxAnimationManagerFlags flags;
	//The update frequency that the animations are recorded at. 60 is the recommended default
	float update_frequency;
	//Any pointer to user data that you want to use in callbacks such
	void *user_data;
    //Callback which you can assign in order to decide if an animation instance should be added to the render queue
    //the next frame. This callback is called inside the UpdateAnimationManager function. Set the callback
    //with SetAnimationManagerCallback
	bool((*maybe_render_instance_callback)(tfx_animation_manager_t *animation_manager, tfx_animation_instance_t *instance, tfx_frame_meta_t *meta, void *user_data));
};

//Use the particle manager to add multiple effects to your scene 
struct tfx_particle_manager_t {
	tfx_vector_t<tfx_soa_buffer_t> particle_array_buffers;
	tfx_bucket_array_t<tfx_particle_soa_t> particle_arrays;

	//In unordered mode emitters that expire have their particle banks added here to be reused
	tfx_storage_map_t<tfx_vector_t<tfxU32>> free_particle_lists;
	//Only used when using distance from camera ordering. New particles are put in this list and then merge sorted into the particles buffer
	tfx_sort_work_entry_t sorting_work_entry[tfxLAYERS];

	tfx_vector_t<tfx_spawn_work_entry_t> spawn_work;
	tfx_vector_t<tfx_control_work_entry_t> control_work;
	tfx_vector_t<tfx_particle_age_work_entry_t> age_work;
	tfx_vector_t<tfxParticleID> particle_indexes;
	tfx_vector_t<tfxU32> free_particle_indexes;
	tfx_vector_t<tfx_depth_index_t> depth_indexes[tfxLAYERS][2];
	tfx_vector_t<tfxU32> effects_in_use[tfxMAXDEPTH][2];
	tfx_vector_t<tfxU32> emitters_in_use[tfxMAXDEPTH][2];
	tfx_vector_t<tfxU32> emitters_check_capture;
	tfx_vector_t<tfxU32> free_effects;
	tfx_vector_t<tfxU32> free_emitters;
	tfx_vector_t<tfx_effect_state_t> effects;
	tfx_vector_t<tfx_emitter_state_t> emitters;
	tfx_library_t *library;

	tfx_work_queue_t work_queue;

	//Banks of sprites for drawing in unordered mode
	tfx_soa_buffer_t sprite_buffer[2][tfxLAYERS];
	tfx_sprite_soa_t sprites[2][tfxLAYERS];
	tfxU32 active_particles_count[tfxLAYERS];
	tfxU32 current_sprite_buffer;
	tfxU32 current_depth_index_buffer;

	//todo: document compute controllers once we've established this is how we'll be doing it.
	void *compute_controller_ptr;
	tfx_vector_t<unsigned int> free_compute_controllers;
	unsigned int new_compute_particle_index;
	unsigned int new_particles_count;
	void *new_compute_particle_ptr;
	//The maximum number of effects that can be updated per frame in the particle manager. If you're running effects with particles that have sub effects then this number might need 
	//to be relatively high depending on your needs. Use Init to udpate the sizes if you need to. Best to call Init at the start with the max numbers that you'll need for your application and don't adjust after.
	unsigned int max_effects;
	//The maximum number of particles that can be updated per frame per layer. #define tfxLAYERS to set the number of allowed layers. This is currently 4 by default
	unsigned int max_cpu_particles_per_layer[tfxLAYERS];
	//The maximum number of particles that can be updated per frame per layer in the compute shader. #define tfxLAYERS to set the number of allowed layers. This is currently 4 by default
	unsigned int max_new_compute_particles;
	//The current effect buffer in use, can be either 0 or 1
	unsigned int current_ebuff;
	unsigned int next_ebuff;

	tfxU32 effects_start_size[tfxMAXDEPTH];
	tfxU32 emitter_start_size[tfxMAXDEPTH];

	tfxU32 sprite_index_point[tfxLAYERS];

	int mt_batch_size;
	std::mutex particle_index_mutex;

	tfx_random_t random;
	unsigned int max_compute_controllers;
	unsigned int highest_compute_controller_index;
	tfx_compute_fx_global_state_t compute_global_state;
	tfxU32 sort_passes;
	tfx_lookup_mode lookup_mode;
	//For when particles are ordered by distance from camera (3d effects)
	tfx_vec3_t camera_front;
	tfx_vec3_t camera_position;

	tfxU32 unique_particle_id = 0;	//Used when recording sprite data
	//When using single particles, you can flag the emitter to set the max_age of the particle to the 
	//length in time of the animation so that it maps nicely to the animation
	float animation_length_in_time;

	//These can possibly be removed at some point, they're debugging variables
	unsigned int particle_id;
	tfxParticleManagerFlags flags;
	//The length of time that passed since the last time Update() was called
	float frame_length;
	tfxWideFloat frame_length_wide;
	float update_time;
	tfxWideFloat update_time_wide;
	float update_frequency;

	tfx_particle_manager_t() :
		flags(0),
		lookup_mode(tfxFast),
		max_effects(10000),
		current_ebuff(0),
		highest_compute_controller_index(0),
		new_compute_particle_ptr(nullptr),
		compute_controller_ptr(nullptr),
		sorting_work_entry{ 0 },
		max_compute_controllers(10000),
		max_new_compute_particles(10000),
		new_compute_particle_index(0),
		new_particles_count(0),
		mt_batch_size(512),
		current_sprite_buffer(0),
		current_depth_index_buffer(0),
		free_compute_controllers(tfxCONSTRUCTOR_VEC_INIT(pm "free_compute_controllers")),
		library(nullptr),
		sort_passes(0)
	{
	}
	~tfx_particle_manager_t();
};

struct tfx_effect_library_stats_t {
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
};

struct tfx_library_t {
	tfx_storage_map_t<tfx_effect_emitter_t*> effect_paths;
	tfx_vector_t<tfx_effect_emitter_t> effects;
	tfx_storage_map_t<tfx_image_data_t> particle_shapes;
	tfx_vector_t<tfx_effect_emitter_info_t> effect_infos;
	tfx_vector_t<tfx_emitter_properties_t> emitter_properties;
	tfx_storage_map_t<tfx_sprite_data_t> pre_recorded_effects;

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

	//Get an effect from the library by index
	tfx_effect_emitter_t& operator[] (tfxU32 index);
	tfx_str64_t name;
	bool open_library = false;
	bool dirty = false;
	tfx_str_t library_file_path;
	tfxU32 uid;

	tfx_library_t() :
		uid(0),
		effect_paths("EffectLib effect paths map", "EffectLib effect paths data"),
		particle_shapes("EffectLib shapes map", "EffectLib shapes data"),
		effects(tfxCONSTRUCTOR_VEC_INIT("effects")),
		effect_infos(tfxCONSTRUCTOR_VEC_INIT("effect_infos")),
		global_graphs(tfxCONSTRUCTOR_VEC_INIT("global_graphs")),
		emitter_attributes(tfxCONSTRUCTOR_VEC_INIT("emitter_attributes")),
		sprite_sheet_settings(tfxCONSTRUCTOR_VEC_INIT("animation_settings")),
		preview_camera_settings(tfxCONSTRUCTOR_VEC_INIT("preview_camera_settings")),
		all_nodes(tfxCONSTRUCTOR_VEC_INIT("all_nodes")),
		node_lookup_indexes(tfxCONSTRUCTOR_VEC_INIT("nodes_lookup_indexes")),
		compiled_lookup_values(tfxCONSTRUCTOR_VEC_INIT("compiled_lookup_values")),
		compiled_lookup_indexes(tfxCONSTRUCTOR_VEC_INIT("compiled_lookup_indexes")),
		graph_min_max(tfxCONSTRUCTOR_VEC_INIT("graph_min_max")),
		free_global_graphs(tfxCONSTRUCTOR_VEC_INIT("free_global_graphs")),
		free_keyframe_graphs(tfxCONSTRUCTOR_VEC_INIT("free_keyframe_graphs")),
		free_emitter_attributes(tfxCONSTRUCTOR_VEC_INIT("free_emitter_attributes")),
		free_animation_settings(tfxCONSTRUCTOR_VEC_INIT("free_animation_settings")),
		free_preview_camera_settings(tfxCONSTRUCTOR_VEC_INIT("free_preview_camera_settings")),
		free_properties(tfxCONSTRUCTOR_VEC_INIT("free_properties")),
		free_infos(tfxCONSTRUCTOR_VEC_INIT("free_infos"))
	{}

	//Free everything in the library

};

struct tfx_effect_template_t {
	tfx_storage_map_t<tfx_effect_emitter_t*> paths;
	tfx_effect_emitter_t effect;
	tfxKey original_effect_hash;

	tfx_effect_template_t() :
		original_effect_hash(0),
		paths("Effect template paths map", "Effect template paths data")
	{}
};

struct tfx_data_entry_t {
	tfx_data_type type = tfxSInt;
	tfx_str32_t key;
	tfx_str_t str_value;
	int int_value = 0;
	bool bool_value = 0;
	float float_value = 0;
	double double_value = 0;
};

//------------------------------------------------------------
//Section: Internal_Functions
//------------------------------------------------------------

tfxAPI tfx_storage_t *GetGlobals();
tfxAPI tfx_pool_stats_t CreateMemorySnapshot(tfx_header *first_block);

tfxINTERNAL void ResizeParticleSoACallback(tfx_soa_buffer_t *buffer, tfxU32 index);

//--------------------------------
//Internal functions used either by the library or editor
//--------------------------------
tfxINTERNAL inline tfxParticleID MakeParticleID(tfxU32 bank_index, tfxU32 particle_index);
tfxINTERNAL inline tfxU32 ParticleIndex(tfxParticleID id);
tfxINTERNAL inline tfxU32 ParticleBank(tfxParticleID id);
//Dump sprites for Debugging
tfxAPI inline void DumpSprites(tfx_particle_manager_t *pm, tfxU32 layer);
tfxINTERNAL tfxU32 GrabParticleLists(tfx_particle_manager_t *pm, tfxKey emitter_hash, tfxU32 reserve_amount = 100);

//--------------------------------
//Profilings
//--------------------------------
tfxAPI_EDITOR void GatherStats(tfx_profile_t *profile, tfx_profile_stats_t *stat);
tfxAPI_EDITOR void ResetSnapshot(tfx_profile_snapshot_t *snapshot);
tfxAPI_EDITOR void ResetSnapshots();
tfxAPI_EDITOR void DumpSnapshots(tfx_storage_map_t<tfx_vector_t<tfx_profile_snapshot_t>> *profile_snapshots, tfxU32 amount);

//--------------------------------
//Reading/Writing files
//--------------------------------
tfxAPI_EDITOR tfx_stream_t ReadEntireFile(const char *file_name, bool terminate = false);
tfxAPI_EDITOR tfxErrorFlags LoadPackage(const char *file_name, tfx_package_t *package);
tfxAPI_EDITOR tfxErrorFlags LoadPackage(tfx_stream_t *stream, tfx_package_t *package);
tfxAPI_EDITOR tfx_package_t CreatePackage(const char *file_path);
tfxAPI_EDITOR bool SavePackageDisk(tfx_package_t *package);
tfxAPI_EDITOR tfx_stream_t SavePackageMemory(tfx_package_t *package);
tfxAPI_EDITOR tfxU64 GetPackageSize(tfx_package_t *package);
tfxAPI_EDITOR bool ValidatePackage(tfx_package_t *package);
tfxAPI_EDITOR tfx_package_entry_info_t *GetPackageFile(tfx_package_t *package, const char *name);
tfxAPI_EDITOR void AddEntryToPackage(tfx_package_t *package, tfx_package_entry_info_t file);
tfxAPI_EDITOR void AddFileToPackage(tfx_package_t *package, const char *file_name, tfx_stream_t *data);
tfxAPI_EDITOR bool FileExists(tfx_package_t *package, const char *file_name);
tfxAPI_EDITOR void FreePackage(tfx_package_t *package);

//Some file IO functions for the editor
tfxAPI_EDITOR bool HasDataValue(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key);
tfxAPI_EDITOR void AddDataValue(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key, const char *value);
tfxAPI_EDITOR void AddDataValue(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key, int value);
tfxAPI_EDITOR void AddDataValue(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key, bool value);
tfxAPI_EDITOR void AddDataValue(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key, double value);
tfxAPI_EDITOR void AddDataValue(tfx_storage_map_t<tfx_data_entry_t> *config, tfx_str32_t key, float value);
tfxAPI_EDITOR tfx_str_t GetDataStrValue(tfx_storage_map_t<tfx_data_entry_t> *config, const char* key);
tfxAPI_EDITOR int GetDataIntValue(tfx_storage_map_t<tfx_data_entry_t> *config, const char* key);
tfxAPI_EDITOR float GetDataFloatValue(tfx_storage_map_t<tfx_data_entry_t> *config, const char* key);
tfxAPI_EDITOR bool SaveDataFile(tfx_storage_map_t<tfx_data_entry_t> *config, const char* path = "");
tfxAPI_EDITOR bool LoadDataFile(tfx_data_types_dictionary_t *data_types, tfx_storage_map_t<tfx_data_entry_t> *config, const char* path);
tfxAPI_EDITOR void StreamProperties(tfx_emitter_properties_t *property, tfxEmitterPropertyFlags flags, tfx_str_t *file);
tfxAPI_EDITOR void StreamProperties(tfx_effect_emitter_t *effect, tfx_str_t *file);
tfxAPI_EDITOR void StreamGraph(const char * name, tfx_graph_t *graph, tfx_str_t *file);
tfxAPI_EDITOR void SplitStringStack(const tfx_str_t s, tfx_vector_t<tfx_str256_t> *pair, char delim = 61);
tfxAPI_EDITOR bool StringIsUInt(const tfx_str_t s);
tfxAPI_EDITOR void AssignEffectorProperty(tfx_effect_emitter_t *effect, tfx_str_t *field, tfxU64 value, tfxU32 file_version);
tfxAPI_EDITOR void AssignEffectorProperty(tfx_effect_emitter_t *effect, tfx_str_t *field, tfxU32 value, tfxU32 file_version);
tfxAPI_EDITOR void AssignEffectorProperty(tfx_effect_emitter_t *effect, tfx_str_t *field, float value);
tfxAPI_EDITOR void AssignEffectorProperty(tfx_effect_emitter_t *effect, tfx_str_t *field, bool value);
tfxAPI_EDITOR void AssignEffectorProperty(tfx_effect_emitter_t *effect, tfx_str_t *field, int value);
tfxAPI_EDITOR void AssignEffectorProperty(tfx_effect_emitter_t *effect, tfx_str_t *field, tfx_str_t &value);
tfxAPI_EDITOR void AssignGraphData(tfx_effect_emitter_t *effect, tfx_vector_t<tfx_str256_t> *values);
tfxINTERNAL void SplitStringVec(const tfx_str_t s, tfx_vector_t<tfx_str256_t> *pair, char delim = 61);
tfxINTERNAL int GetDataType(const tfx_str_t &s);
tfxINTERNAL void AssignStageProperty(tfx_effect_emitter_t *effect, tfx_str_t *field, tfxU32 value);
tfxINTERNAL void AssignStageProperty(tfx_effect_emitter_t *effect, tfx_str_t *field, float value);
tfxINTERNAL void AssignStageProperty(tfx_effect_emitter_t *effect, tfx_str_t *field, bool value);
tfxINTERNAL void AssignStageProperty(tfx_effect_emitter_t *effect, tfx_str_t *field, int value);
tfxINTERNAL void AssignStageProperty(tfx_effect_emitter_t *effect, tfx_str_t *field, tfx_str_t *value);
tfxINTERNAL void AssignSpriteDataMetricsProperty(tfx_sprite_data_metrics_t *metrics, tfx_str_t *field, tfxU32 value, tfxU32 file_version);
tfxINTERNAL void AssignSpriteDataMetricsProperty(tfx_sprite_data_metrics_t *metrics, tfx_str_t *field, tfxU64 value, tfxU32 file_version);
tfxINTERNAL void AssignSpriteDataMetricsProperty(tfx_sprite_data_metrics_t *metrics, tfx_str_t *field, float value, tfxU32 file_version);
tfxINTERNAL void AssignSpriteDataMetricsProperty(tfx_sprite_data_metrics_t *metrics, tfx_str_t *field, tfx_str_t value, tfxU32 file_version);
tfxINTERNAL void AssignFrameMetaProperty(tfx_frame_meta_t *metrics, tfx_str_t *field, tfxU32 value, tfxU32 file_version);
tfxINTERNAL void AssignFrameMetaProperty(tfx_frame_meta_t *metrics, tfx_str_t *field, tfx_vec3_t value, tfxU32 file_version);
tfxINTERNAL void AssignAnimationEmitterProperty(tfx_animation_emitter_properties_t *properties, tfx_str_t *field, tfxU32 value, tfxU32 file_version);
tfxINTERNAL void AssignAnimationEmitterProperty(tfx_animation_emitter_properties_t *properties, tfx_str_t *field, float value, tfxU32 file_version);
tfxINTERNAL void AssignAnimationEmitterProperty(tfx_animation_emitter_properties_t *properties, tfx_str_t *field, tfx_vec2_t value, tfxU32 file_version);
tfxINTERNAL void AssignNodeData(tfx_attribute_node_t *node, tfx_vector_t<tfx_str256_t> *values);
tfxINTERNAL tfx_vec3_t StrToVec3(tfx_vector_t<tfx_str256_t> *str);
tfxINTERNAL tfx_vec2_t StrToVec2(tfx_vector_t<tfx_str256_t> *str);

//--------------------------------
//Inline Math functions
//--------------------------------
tfxINTERNAL void MakeIcospheres();
tfxINTERNAL int VertexForEdge(tfx_storage_map_t<int> *point_cache, tfx_vector_t<tfx_vec3_t> *vertices, int first, int second);
tfxINTERNAL tfx_vector_t<tfx_face_t> SubDivideIcosphere(tfx_storage_map_t<int> *point_cache, tfx_vector_t<tfx_vec3_t> *vertices, tfx_vector_t<tfx_face_t> *triangles);

tfxINTERNAL int SortIcospherePoints(void const *left, void const *right);
tfxINTERNAL int SortDepth(void const *left, void const *right);
tfxINTERNAL void InsertionSortDepth(tfx_work_queue_t *queue, void *work_entry);
tfxAPI_EDITOR void InsertionSortParticleFrame(tfx_vector_t<tfx_particle_frame_t> *particles);
tfxINTERNAL tfx128 Dot128XYZ(const tfx128 *x1, const tfx128 *y1, const tfx128 *z1, const tfx128 *x2, const tfx128 *y2, const tfx128 *z);
tfxINTERNAL tfx128 Dot128XY(const tfx128 *x1, const tfx128 *y1, const tfx128 *x2, const tfx128 *y2);
tfxAPI_EDITOR tfx_hsv_t RGBtoHSV(tfx_rgb_t in);
tfxAPI_EDITOR tfx_rgb_t HSVtoRGB(tfx_hsv_t in);
tfxAPI_EDITOR float DegreesToRadians(float degrees);
tfxAPI_EDITOR float RadiansToDegrees(float radians);

tfxAPI_EDITOR float LengthVec3NoSqR(tfx_vec3_t const *v);
tfxINTERNAL float LengthVec4NoSqR(tfx_vec4_t const *v);
tfxAPI_EDITOR float LengthVec(tfx_vec3_t const *v);
tfxINTERNAL float LengthVec(tfx_vec4_t const *v);
tfxINTERNAL float HasLength(tfx_vec3_t const *v);
tfxINTERNAL tfx_vec3_t NormalizeVec3(tfx_vec3_t const *v);
tfxINTERNAL tfx_vec4_t NormalizeVec4(tfx_vec4_t const *v);
tfxINTERNAL tfx_vec3_t Cross(tfx_vec3_t *a, tfx_vec3_t *b);
tfxINTERNAL float DotProductVec4(const tfx_vec4_t *a, const tfx_vec4_t *b);
tfxINTERNAL float DotProductVec3(const tfx_vec3_t *a, const tfx_vec3_t *b);
tfxAPI_EDITOR float DotProductVec2(const tfx_vec2_t *a, const tfx_vec2_t *b);
//Quake 3 inverse square root
tfxINTERNAL float QuakeSqrt(float number);
tfxINTERNAL tfxU32 GetLayerFromID(tfxU32 index);
tfxINTERNAL tfxU32 GetIndexFromID(tfxU32 index);
//Todo: can delete this now?
tfxINTERNAL tfxU32 SetNibbleID(tfxU32 nibble, tfxU32 index);
tfxAPI_EDITOR float Vec2LengthFast(tfx_vec2_t const *v);
tfxAPI_EDITOR float Vec3FastLength(tfx_vec3_t const *v);
tfxAPI_EDITOR tfx_vec3_t NormalizeVec3Fast(tfx_vec3_t const *v);
tfxINTERNAL tfx_vec2_t NormalizeVec2(tfx_vec2_t const *v);
tfxAPI_EDITOR tfx_mat3_t CreateMatrix3(float v = 1.f);
tfxINTERNAL tfx_mat3_t TranslateMatrix3Vec3(tfx_mat3_t const *m, tfx_vec3_t const *v);
tfxAPI_EDITOR tfx_mat3_t RotateMatrix3(tfx_mat3_t const *m, float r);
tfxINTERNAL tfx_mat3_t ScaleMatrix3Vec2(tfx_mat3_t const *m, tfx_vec2_t const &v);
tfxAPI_EDITOR tfx_vec2_t TransformVec2Matrix4(const tfx_mat4_t *mat, const tfx_vec2_t v);
tfxINTERNAL tfx_mat4_t CreateMatrix4(float v);
tfxINTERNAL tfx_mat4_t Matrix4FromVecs(tfx_vec4_t a, tfx_vec4_t b, tfx_vec4_t c, tfx_vec4_t d);
tfxINTERNAL tfx_mat4_t Matrix4RotateX(float angle);
tfxINTERNAL tfx_mat4_t Matrix4RotateY(float angle);
tfxINTERNAL tfx_mat4_t Matrix4RotateZ(float angle);
tfxINTERNAL tfx_mat4_t TransposeMatrix4(tfx_mat4_t *mat);
tfxINTERNAL tfx_mat4_t TransformMatrix42d(const tfx_mat4_t *in, const tfx_mat4_t *m);
tfxINTERNAL tfx_mat4_t TransformMatrix4ByMatrix2(const tfx_mat4_t *in, const tfx_mat2_t *m);
tfxINTERNAL tfx_mat4_t TransformMatrix4(const tfx_mat4_t *in, const tfx_mat4_t *m);
tfxINTERNAL void TransformMatrix4Vec3(const tfx_mat4_t *mat, tfxWideFloat *x, tfxWideFloat *y, tfxWideFloat *z);
tfxINTERNAL void TransformMatrix4Vec2(const tfx_mat4_t *mat, tfxWideFloat *x, tfxWideFloat *y);
tfxINTERNAL void MaskedTransformMatrix2(const tfxWideFloat *r0c, const tfxWideFloat *r1c, tfxWideFloat *x, tfxWideFloat *y, tfxWideFloat *mask, tfxWideFloat *xor_mask);
tfxINTERNAL void MaskedTransformMatrix42d(const tfx_mat4_t *mat, tfxWideFloat *x, tfxWideFloat *y, tfxWideFloat *mask, tfxWideFloat *xor_mask);
tfxINTERNAL void MaskedTransformMatrix4Vec3(const tfxWideFloat *r0c, const tfxWideFloat *r1c, const tfxWideFloat *r2c, tfxWideFloat *x, tfxWideFloat *y, tfxWideFloat *z, tfxWideFloat *mask, tfxWideFloat *xor_mask);
tfxINTERNAL tfx_vec4_t TransformVec4Matrix4(const tfx_mat4_t *mat, const tfx_vec4_t vec);
tfxAPI_EDITOR tfx_vec4_t WideTransformVec4Matrix4(const tfx128 *row1, const tfx128 *row2, const tfx128 *row3, const tfx128 *row4, const tfx_vec4_t vec);
tfxINTERNAL tfx_vec3_t TransformVec3Matrix4(const tfx_mat4_t *mat, const tfx_vec4_t *vec);
tfxINTERNAL tfx_mat4_t Matrix4RotateAxis(tfx_mat4_t const *m, float r, tfx_vec3_t const *v);
tfxAPI_EDITOR int tfxClampi(int lower, int upper, int value);
tfxAPI_EDITOR float tfxClampf(float lower, float upper, float value);
tfxAPI_EDITOR tfxU32 Pack10bit(tfx_vec3_t const *v, tfxU32 extra);
tfxINTERNAL tfxU32 Pack10bitUnsigned(tfx_vec3_t const *v);
tfxAPI_EDITOR tfxU32 Pack16bit(float x, float y);
tfxAPI_EDITOR tfxU32 Pack8bit(tfx_vec3_t v);
tfxINTERNAL tfxU32 Pack16bitUnsigned(float x, float y);
tfxAPI_EDITOR tfx_vec2_t UnPack16bit(tfxU32 in);
tfxINTERNAL tfx_vec2_t UnPack16bitUnsigned(tfxU32 in);
tfxINTERNAL tfxWideInt PackWide16bitStretch(tfxWideFloat &v_x, tfxWideFloat &v_y);
tfxINTERNAL tfxWideInt PackWide16bit(tfxWideFloat &v_x, tfxWideFloat &v_y);
tfxINTERNAL void UnPackWide16bit(tfxWideInt in, tfxWideFloat &x, tfxWideFloat &y);
tfxAPI tfxWideInt PackWide8bitXYZ(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z);
tfxINTERNAL tfxWideInt PackWide10bit(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z);
tfxINTERNAL tfxWideInt PackWide10bit(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z, tfxU32 extra);
tfxINTERNAL tfxWideInt PackWide10bitUnsigned(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z, tfxU32 extra);
tfxINTERNAL void UnPackWide10bit(tfxWideInt in, tfxWideFloat &x, tfxWideFloat &y, tfxWideFloat &z);
tfxINTERNAL tfxWideFloat UnPackWide10bitY(tfxWideInt in);
tfxINTERNAL tfxWideInt PackWideColor(tfxWideFloat const &v_r, tfxWideFloat const &v_g, tfxWideFloat const &v_b, tfxWideFloat v_a);
tfxINTERNAL tfxWideInt PackWide10bit(tfxWideFloat const &v_x, tfxWideFloat const &v_y, tfxWideFloat const &v_z, tfxWideInt extra);
tfxAPI_EDITOR tfx_vec4_t UnPack10bit(tfxU32 in);
tfxAPI_EDITOR tfx_vec3_t UnPack8bit(tfxU32 in);
tfxINTERNAL tfx_vec3_t UnPack10bitVec3(tfxU32 in);
tfxINTERNAL tfxU32 Get2bitFromPacked10bit(tfxU32 in);
tfxINTERNAL size_t ClampStringSize(size_t compare, size_t string_size);
tfxAPI_EDITOR float Distance2d(float fromx, float fromy, float tox, float toy);
tfxINTERNAL tfxUInt10bit UintToPacked10bit(tfxU32 in);
tfxINTERNAL tfx_vec2_t InterpolateVec2(float tween, tfx_vec2_t from, tfx_vec2_t to);
tfxINTERNAL tfx_vec3_t InterpolateVec3(float tween, tfx_vec3_t from, tfx_vec3_t to);
tfxINTERNAL tfx_rgba8_t InterpolateRGBA(float tween, tfx_rgba8_t from, tfx_rgba8_t to);
tfxINTERNAL float GammaCorrect(float color, float gamma = tfxGAMMA);
tfxINTERNAL tfxU32 InterpolateAlignment(float tween, tfxU32 from, tfxU32 to);
tfxINTERNAL tfx_vec4_t InterpolateVec4(float tween, tfx_vec4_t *from, tfx_vec4_t *to);
tfxINTERNAL tfxWideFloat WideInterpolate(tfxWideFloat tween, tfxWideFloat *from, tfxWideFloat *to);
tfxINTERNAL float Interpolatef(float tween, float from, float to);
tfxINTERNAL void Transform2d(tfx_vec3_t *out_rotations, tfx_vec3_t *out_local_rotations, float *out_scale, tfx_vec3_t *out_position, tfx_vec3_t *out_local_position, tfx_vec3_t *out_translation, tfx_mat4_t *out_matrix, tfx_effect_state_t *parent);
tfxAPI_EDITOR void Transform3d(tfx_vec3_t *out_rotations, tfx_vec3_t *out_local_rotations, float *out_scale, tfx_vec3_t *out_position, tfx_vec3_t *out_local_position, tfx_vec3_t *out_translation, tfx_mat4_t *out_matrix, const tfx_effect_state_t *parent);
//-------------------------------------------------
//--New transform_3d particle functions for SoA data--
//--------------------------2d---------------------
tfxINTERNAL void TransformParticlePosition(const float local_position_x, const float local_position_y, const float roll, tfx_vec2_t *world_position, float *world_rotations, const tfx_vec3_t *parent_rotations, const tfx_mat4_t *matrix, const tfx_vec3_t *handle, const float *scale, const tfx_vec3_t *from_position);

//--------------------------------
//Random numbers
//--------------------------------
tfx_random_t NewRandom(tfxU32 seed);

void AdvanceRandom(tfx_random_t *random);
void RandomReSeed(tfx_random_t *random);
void RandomReSeed(tfx_random_t *random, tfxU64 seed1, tfxU64 seed2);
void RandomReSeed(tfx_random_t *random, tfxU64 seed);
float GenerateRandom(tfx_random_t *random);
float RandomRange(tfx_random_t *random, float max);
float RandomRange(tfx_random_t *random, float from, float to);
int RandomRange(tfx_random_t *random, int from, int to);
tfxU32 RandomRange(tfx_random_t *random, tfxU32 max);
void AlterRandomSeed(tfx_random_t *random, tfxU64 amount);
void AlterRandomSeed(tfx_random_t *random, tfxU32 amount);

//--------------------------------
//Particle manager internal functions
//--------------------------------
tfxINTERNAL float GetEmissionDirection2d(tfx_particle_manager_t *pm, tfx_library_t *library, tfx_random_t *random, tfx_emitter_state_t &emitter, tfx_vec2_t local_position, tfx_vec2_t world_position);
tfxINTERNAL tfx_vec3_t GetEmissionDirection3d(tfx_particle_manager_t *pm, tfx_library_t *library, tfx_random_t *random, tfx_emitter_state_t &emitter, float emission_pitch, float emission_yaw, tfx_vec3_t local_position, tfx_vec3_t world_position);
tfxINTERNAL void TransformEffector2d(tfx_vec3_t *world_rotations, tfx_vec3_t *local_rotations, tfx_vec3_t *world_position, tfx_vec3_t *local_position, tfx_mat4_t *matrix, tfx_sprite_transform2d_t *parent, bool relative_position = true, bool relative_angle = false);
tfxINTERNAL void TransformEffector3d(tfx_vec3_t *world_rotations, tfx_vec3_t *local_rotations, tfx_vec3_t *world_position, tfx_vec3_t *local_position, tfx_mat4_t *matrix, tfx_sprite_transform3d_t *parent, bool relative_position = true, bool relative_angle = false);
tfxINTERNAL void UpdatePMEffect(tfx_particle_manager_t *pm, tfxU32 index, tfxU32 parent_index = tfxINVALID);
tfxINTERNAL void UpdatePMEmitter(tfx_work_queue_t *work_queue, void *data);
tfxINTERNAL tfxU32 NewSpritesNeeded(tfx_particle_manager_t *pm, tfxU32 index, tfx_effect_state_t *parent, tfx_emitter_properties_t *properties);
tfxINTERNAL void UpdateEmitterState(tfx_particle_manager_t *pm, tfx_emitter_state_t &emitter, tfxU32 parent_index, const tfx_parent_spawn_controls_t *parent_spawn_controls, tfx_spawn_work_entry_t *entry);
tfxINTERNAL void UpdateEffectState(tfx_particle_manager_t *pm, tfxU32 index);

tfxAPI_EDITOR void CompletePMWork(tfx_particle_manager_t *pm);

tfxINTERNAL tfxU32 SpawnParticles2d(tfx_particle_manager_t *pm, tfx_spawn_work_entry_t *spawn_work_entry, tfxU32 max_spawn_count);
tfxINTERNAL void SpawnParticlePoint2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleLine2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleArea2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleEllipse2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleMicroUpdate2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleNoise(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void SpawnParticleWeight(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleVelocity(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleRoll(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleImageFrame(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleAge(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleSize2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleSpin2d(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void DoSpawnWork3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void DoSpawnWork2d(tfx_work_queue_t *queue, void *data);

tfxINTERNAL tfxU32 SpawnParticles3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticlePoint3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleLine3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleArea3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleEllipse3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleCylinder3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleIcosphereRandom3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleIcosphere3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleMicroUpdate3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleSpin3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void SpawnParticleSize3d(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void ControlParticles(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void ControlParticleAge(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void ControlParticleImageFrame(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void ControlParticleColor(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void ControlParticleSize(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void ControlParticleUID(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void ControlParticleCaptureFlag(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void ControlParticlePosition2d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void ControlParticleTransform2d(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void ControlParticlePosition3d(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void ControlParticleTransform3d(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void ControlParticleBoundingBox(tfx_work_queue_t *queue, void *data);

tfxINTERNAL void InitSpriteData3dSoACompression(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void InitSpriteData3dSoA(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void InitSpriteData2dSoACompression(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void InitSpriteData2dSoA(tfx_soa_buffer_t *buffer, tfx_sprite_data_soa_t *soa, tfxU32 reserve_amount);
tfxINTERNAL void InitSpriteBufferSoA(tfx_soa_buffer_t *buffer, tfx_sprite_soa_t *soa, tfxU32 reserve_amount, tfxSpriteBufferMode mode, bool use_uid = false);
tfxINTERNAL void InitParticleSoA(tfx_soa_buffer_t *buffer, tfx_particle_soa_t *soa, tfxU32 reserve_amount);

tfxAPI_EDITOR void InitEmitterProperites(tfx_emitter_properties_t *properties);
tfxINTERNAL void CopyEmitterProperites(tfx_emitter_properties_t *from_properties, tfx_emitter_properties_t *to_properties);

tfxINTERNAL inline void FreeSpriteData(tfx_sprite_data_t *sprite_data);

//--------------------------------
//Graph functions
//Mainly used by the editor to edit graphs so these are kind of API functions but you wouldn't generally use these outside of the particle editor
//--------------------------------
tfxAPI_EDITOR tfx_attribute_node_t* AddGraphNode(tfx_graph_t *graph, float frame, float value, tfxAttributeNodeFlags flags = 0, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
tfxAPI_EDITOR void AddGraphNode(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxAPI_EDITOR void SetGraphNode(tfx_graph_t *graph, tfxU32 index, float frame, float value, tfxAttributeNodeFlags flags = 0, float x1 = 0, float y1 = 0, float x2 = 0, float y2 = 0);
tfxAPI_EDITOR float GetGraphValue(tfx_graph_t *graph, float age);
tfxAPI_EDITOR float GetGraphRandomValue(tfx_graph_t *graph, float age, tfx_random_t *seed);
tfxAPI_EDITOR float GetGraphValue(tfx_graph_t *graph, float age, float life);
tfxAPI_EDITOR tfx_attribute_node_t *GetGraphNextNode(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxAPI_EDITOR tfx_attribute_node_t *GetGraphPrevNode(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxAPI_EDITOR tfx_attribute_node_t *GetGraphLastNode(tfx_graph_t *graph);
tfxAPI_EDITOR float GetGraphFirstValue(tfx_graph_t *graph);
tfxAPI_EDITOR tfx_attribute_node_t* AddGraphCoordNode(tfx_graph_t *graph, float, float);
tfxAPI_EDITOR tfx_attribute_node_t* InsertGraphCoordNode(tfx_graph_t *graph, float, float);
tfxAPI_EDITOR tfx_attribute_node_t* InsertGraphNode(tfx_graph_t *graph, float, float);
tfxAPI_EDITOR float *LinkGraphFirstValue(tfx_graph_t *graph);
tfxAPI_EDITOR float GetGraphLastValue(tfx_graph_t *graph);
tfxAPI_EDITOR float GetGraphMaxValue(tfx_graph_t *graph);
tfxAPI_EDITOR float GetGraphMinValue(tfx_graph_t *graph);
tfxAPI_EDITOR float GetGraphLastFrame(tfx_graph_t *graph, float udpate_frequence);
tfxAPI_EDITOR tfx_attribute_node_t* GraphNodeByIndex(tfx_graph_t *graph, tfxU32 index);
tfxAPI_EDITOR float GraphValueByIndex(tfx_graph_t *graph, tfxU32 index);
tfxAPI_EDITOR float GraphFrameByIndex(tfx_graph_t *graph, tfxU32 index);
tfxAPI_EDITOR tfx_attribute_node_t* FindGraphNode(tfx_graph_t *graph, tfx_attribute_node_t *n);
tfxAPI_EDITOR void ValidateGraphCurves(tfx_graph_t *graph);
tfxAPI_EDITOR void DeleteGraphNode(tfx_graph_t *graph, tfx_attribute_node_t *n);
tfxAPI_EDITOR void DeleteGraphNodeAtFrame(tfx_graph_t *graph, float frame);
tfxAPI_EDITOR void ResetGraph(tfx_graph_t *graph, float first_node_value, tfx_graph_preset preset, bool add_node = true);
tfxAPI_EDITOR void ClearGraphToOne(tfx_graph_t *graph, float value);
tfxAPI_EDITOR void ClearGraph(tfx_graph_t *graph);
tfxAPI_EDITOR void FreeGraph(tfx_graph_t *graph);
tfxAPI_EDITOR void CopyGraph(tfx_graph_t *graph, tfx_graph_t *to, bool compile = true);
tfxAPI_EDITOR bool SortGraph(tfx_graph_t *graph);
tfxAPI_EDITOR void ReIndexGraph(tfx_graph_t *graph);
tfxAPI_EDITOR tfx_vec2_t GetGraphInitialZoom(tfx_graph_t *graph);
tfxAPI_EDITOR tfx_vec2_t GetGraphInitialZoom3d(tfx_graph_t *graph);
tfxAPI_EDITOR bool IsColorGraph(tfx_graph_t *graph);
tfxAPI_EDITOR bool IsOvertimeGraph(tfx_graph_t *graph);
tfxAPI_EDITOR bool IsGlobalGraph(tfx_graph_t *graph);
tfxAPI_EDITOR bool IsAngleGraph(tfx_graph_t *graph);
tfxAPI_EDITOR bool IsTranslationGraph(tfx_graph_t *graph);
tfxAPI_EDITOR void MultiplyAllGraphValues(tfx_graph_t *graph, float scalar);
tfxAPI_EDITOR void CopyGraphNoLookups(tfx_graph_t *src_graph, tfx_graph_t *dst_graph);
tfxAPI_EDITOR void DragGraphValues(tfx_graph_preset preset, float *frame, float *value);
tfxAPI_EDITOR tfx_vec4_t GetMinMaxGraphValues(tfx_graph_preset preset);
tfxAPI_EDITOR tfx_vec2_t GetQuadBezier(tfx_vec2_t p0, tfx_vec2_t p1, tfx_vec2_t p2, float t, float ymin, float ymax, bool clamp = true);
tfxAPI_EDITOR tfx_vec2_t GetCubicBezier(tfx_vec2_t p0, tfx_vec2_t p1, tfx_vec2_t p2, tfx_vec2_t p3, float t, float ymin, float ymax, bool clamp = true);
tfxAPI_EDITOR float GetBezierValue(const tfx_attribute_node_t *lastec, const tfx_attribute_node_t *a, float t, float ymin, float ymax);
tfxAPI_EDITOR float GetDistance(float fromx, float fromy, float tox, float toy);
tfxAPI_EDITOR float inline GetVectorAngle(float x, float y) { return atan2(x, -y); }
tfxAPI_EDITOR bool CompareNodes(tfx_attribute_node_t *left, tfx_attribute_node_t *right);
tfxAPI_EDITOR void CompileGraph(tfx_graph_t *graph);
tfxAPI_EDITOR void CompileGraphOvertime(tfx_graph_t *graph);
tfxAPI_EDITOR void CompileColorOvertime(tfx_graph_t *graph, float gamma = tfxGAMMA);
tfxAPI_EDITOR float GetMaxLife(tfx_effect_emitter_t *e);
tfxAPI_EDITOR float LookupFastOvertime(tfx_graph_t *graph, float age, float lifetime);
tfxAPI_EDITOR float LookupFast(tfx_graph_t *graph, float frame);
tfxAPI_EDITOR float LookupPreciseOvertime(tfx_graph_t *graph, float age, float lifetime);
tfxAPI_EDITOR float LookupPrecise(tfx_graph_t *graph, float frame);
tfxAPI_EDITOR float GetRandomFast(tfx_graph_t *graph, float frame, tfx_random_t *random);
tfxAPI_EDITOR float GetRandomPrecise(tfx_graph_t *graph, float frame, tfx_random_t *random);

//Node Manipulation
tfxAPI_EDITOR bool SetNode(tfx_graph_t *graph, tfx_attribute_node_t *node, float, float, tfxAttributeNodeFlags flags, float = 0, float = 0, float = 0, float = 0);
tfxAPI_EDITOR bool SetNode(tfx_graph_t *graph, tfx_attribute_node_t *node, float *frame, float *value);
tfxAPI_EDITOR void SetCurve(tfx_graph_t *graph, tfx_attribute_node_t *node, bool is_left_curve, float *frame, float *value);
tfxAPI_EDITOR bool MoveNode(tfx_graph_t *graph, tfx_attribute_node_t *node, float frame, float value, bool sort = true);
tfxAPI_EDITOR bool SetNodeFrame(tfx_graph_t *graph, tfx_attribute_node_t *node, float *frame);
tfxAPI_EDITOR bool SetNodeValue(tfx_graph_t *graph, tfx_attribute_node_t *node, float *value);
tfxAPI_EDITOR void ClampNode(tfx_graph_t *graph, tfx_attribute_node_t *node);
tfxAPI_EDITOR void ClampCurve(tfx_graph_t *graph, tfx_vec2_t *curve, tfx_attribute_node_t *node);
tfxAPI_EDITOR void ClampGraph(tfx_graph_t *graph);
tfxAPI_EDITOR bool IsOvertimeGraph(tfx_graph_type type);
tfxAPI_EDITOR bool IsColorGraph(tfx_graph_type type);
tfxAPI_EDITOR bool IsOvertimePercentageGraph(tfx_graph_type type);
tfxAPI_EDITOR bool IsGlobalGraph(tfx_graph_type type);
tfxAPI_EDITOR bool IsEmitterGraph(tfx_graph_type type);
tfxAPI_EDITOR bool IsTransformGraph(tfx_graph_type type);
tfxAPI_EDITOR bool IsGlobalPercentageGraph(tfx_graph_type type);
tfxAPI_EDITOR bool IsAngleGraph(tfx_graph_type type);
tfxAPI_EDITOR bool IsAngleOvertimeGraph(tfx_graph_type type);
tfxAPI_EDITOR bool IsEverythingElseGraph(tfx_graph_type type);
tfxAPI_EDITOR bool HasNodeAtFrame(tfx_graph_t *graph, float frame);
tfxAPI_EDITOR bool HasKeyframes(tfx_effect_emitter_t *e);
tfxAPI_EDITOR bool HasMoreThanOneKeyframe(tfx_effect_emitter_t *e);
tfxAPI_EDITOR void PushTranslationPoints(tfx_effect_emitter_t *e, tfx_vector_t<tfx_vec3_t> *points, float frame);

tfxAPI_EDITOR bool IsNodeCurve(tfx_attribute_node_t *node);
tfxAPI_EDITOR bool NodeCurvesAreInitialised(tfx_attribute_node_t *node);
tfxAPI_EDITOR bool SetNodeCurveInitialised(tfx_attribute_node_t *node);

tfxINTERNAL inline bool IsGraphTransformRotation(tfx_graph_type type) {
	return type == tfxTransform_roll || type == tfxTransform_pitch || type == tfxTransform_yaw;
}

tfxINTERNAL inline bool IsGraphEmitterDimension(tfx_graph_type type) {
	return type == tfxProperty_emitter_width || type == tfxProperty_emitter_height || type == tfxProperty_emitter_depth;
}

tfxINTERNAL inline bool IsGraphTranslation(tfx_graph_type type) {
	return type == tfxTransform_translate_x || type == tfxTransform_translate_y || type == tfxTransform_translate_z;
}

tfxINTERNAL inline bool IsGraphEmission(tfx_graph_type type) {
	return type == tfxProperty_emission_pitch || type == tfxProperty_emission_yaw;
}

tfxINTERNAL inline bool IsGraphParticleSize(tfx_graph_type type) {
	return	type == tfxBase_width || type == tfxBase_height ||
		type == tfxVariation_width || type == tfxVariation_height ||
		type == tfxOvertime_width || type == tfxOvertime_height;
}

//--------------------------------
//Grouped graph struct functions
//--------------------------------
tfxINTERNAL void InitialiseGlobalAttributes(tfx_global_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void InitialiseOvertimeAttributes(tfx_overtime_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void InitialiseVariationAttributes(tfx_variation_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void InitialiseBaseAttributes(tfx_base_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void InitialisePropertyAttributes(tfx_property_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void InitialiseTransformAttributes(tfx_transform_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void InitialiseEmitterAttributes(tfx_emitter_attributes_t *attributes, tfxU32 bucket_size = 8);
tfxINTERNAL void FreeEmitterAttributes(tfx_emitter_attributes_t *attributes);
tfxINTERNAL void FreeGlobalAttributes(tfx_global_attributes_t *attributes);
tfxAPI_EDITOR void FreeOvertimeAttributes(tfx_overtime_attributes_t *attributes);
tfxAPI_EDITOR void CopyOvertimeAttributesNoLookups(tfx_overtime_attributes_t *src, tfx_overtime_attributes_t *dst);
tfxAPI_EDITOR void CopyOvertimeAttributes(tfx_overtime_attributes_t *src, tfx_overtime_attributes_t *dst);
tfxAPI_EDITOR void FreeVariationAttributes(tfx_variation_attributes_t *attributes);
tfxAPI_EDITOR void CopyVariationAttributesNoLookups(tfx_variation_attributes_t *src, tfx_variation_attributes_t *dst);
tfxAPI_EDITOR void CopyVariationAttributes(tfx_variation_attributes_t *src, tfx_variation_attributes_t *dst);
tfxAPI_EDITOR void FreeBaseAttributes(tfx_base_attributes_t *attributes);
tfxAPI_EDITOR void CopyBaseAttributesNoLookups(tfx_base_attributes_t *src, tfx_base_attributes_t *dst);
tfxAPI_EDITOR void CopyBaseAttributes(tfx_base_attributes_t *src, tfx_base_attributes_t *dst);
tfxAPI_EDITOR void FreePropertyAttributes(tfx_property_attributes_t *attributes);
tfxAPI_EDITOR void CopyPropertyAttributesNoLookups(tfx_property_attributes_t *src, tfx_property_attributes_t *dst);
tfxAPI_EDITOR void CopyPropertyAttributes(tfx_property_attributes_t *src, tfx_property_attributes_t *dst);
tfxAPI_EDITOR void FreeTransformAttributes(tfx_transform_attributes_t *attributes);
tfxAPI_EDITOR void CopyTransformAttributesNoLookups(tfx_transform_attributes_t *src, tfx_transform_attributes_t *dst);
tfxAPI_EDITOR void CopyTransformAttributes(tfx_transform_attributes_t *src, tfx_transform_attributes_t *dst);
tfxAPI_EDITOR bool HasTranslationKeyframes(tfx_transform_attributes_t *graphs);
tfxAPI_EDITOR void AddTranslationNodes(tfx_transform_attributes_t *keyframes, float frame);
tfxAPI_EDITOR void CopyGlobalAttributesNoLookups(tfx_global_attributes_t *src, tfx_global_attributes_t *dst);
tfxAPI_EDITOR void CopyGlobalAttributes(tfx_global_attributes_t *src, tfx_global_attributes_t *dst);

//Get a graph by tfx_graph_id_t
tfxAPI_EDITOR tfx_graph_t *GetGraph(tfx_library_t *library, tfx_graph_id_t graph_id);

tfxAPI_EDITOR int GetEffectLibraryStats(const char *filename, tfx_effect_library_stats_t *stats);
tfxAPI_EDITOR tfx_effect_library_stats_t CreateLibraryStats(tfx_library_t *lib);
tfxINTERNAL tfxErrorFlags LoadEffectLibraryPackage(tfx_package_t *package, tfx_library_t *lib, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void *user_data = nullptr, bool read_only = true);

//--------------------------------
//Animation manager internal functions - animation manager is used to playback pre-recorded effects
//--------------------------------
tfxINTERNAL tfxAnimationID AddAnimationInstance(tfx_animation_manager_t *animation_manager);
tfxINTERNAL void FreeAnimationInstance(tfx_animation_manager_t *animation_manager, tfxU32 index);
tfxINTERNAL void AddEffectEmitterProperties(tfx_animation_manager_t *animation_manager, tfx_effect_emitter_t *effect, bool *has_animated_shape);
tfxAPI void UpdateAnimationManagerBufferMetrics(tfx_animation_manager_t *animation_manager);
tfxINTERNAL bool FreePMEffectCapacity(tfx_particle_manager_t *pm);
tfxINTERNAL void InitialiseAnimationManager(tfx_animation_manager_t *animation_manager, tfxU32 max_instances);

//--------------------------------
//Particle manager internal functions
//--------------------------------
tfxINTERNAL tfxU32 GetPMEffectSlot(tfx_particle_manager_t *pm);
tfxINTERNAL tfxU32 GetPMEmitterSlot(tfx_particle_manager_t *pm);
tfxINTERNAL tfxU32 GetPMParticleIndexSlot(tfx_particle_manager_t *pm, tfxParticleID particle_id);
tfxINTERNAL void FreePMParticleIndex(tfx_particle_manager_t *pm, tfxU32 *index);
tfxINTERNAL tfxU32 PushPMDepthIndex(tfx_particle_manager_t *pm, tfxU32 layer, tfx_depth_index_t depth_index);
tfxINTERNAL void ResetPMFlags(tfx_particle_manager_t *pm);
tfxINTERNAL tfxU32 GetParticleSpriteIndex(tfx_particle_manager_t *pm, tfxParticleID id);
tfxINTERNAL unsigned int GetControllerMemoryUsage(tfx_particle_manager_t *pm);
tfxINTERNAL unsigned int GetParticleMemoryUsage(tfx_particle_manager_t *pm);
tfxINTERNAL void FreeComputeSlot(tfx_particle_manager_t *pm, unsigned int slot_id);
tfxINTERNAL tfxEffectID AddEffectToParticleManager(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect, int buffer, int hierarchy_depth, bool is_sub_emitter, tfxU32 root_effect_index, float add_delayed_spawning);
tfxINTERNAL void ToggleSpritesWithUID(tfx_particle_manager_t *pm, bool switch_on);
tfxINTERNAL void FreeParticleList(tfx_particle_manager_t *pm, tfxU32 index);
tfxINTERNAL void FreeParticleBanks(tfx_particle_manager_t *pm);

//Compute stuff doesn't work currently
tfxINTERNAL void EnableCompute(tfx_particle_manager_t *pm) { pm->flags |= tfxEffectManagerFlags_use_compute_shader; }
tfxINTERNAL void DisableCompute(tfx_particle_manager_t *pm) { pm->flags &= ~tfxEffectManagerFlags_use_compute_shader; }
tfxINTERNAL int AddComputeController(tfx_particle_manager_t *pm);
tfxINTERNAL tfx_compute_particle_t *GrabComputeParticle(tfx_particle_manager_t *pm, unsigned int layer);
tfxINTERNAL void ResetParticlePtr(tfx_particle_manager_t *pm, void *ptr);
tfxINTERNAL void ResetControllerPtr(tfx_particle_manager_t *pm, void *ptr);
tfxINTERNAL void UpdateCompute(tfx_particle_manager_t *pm, void *sampled_particles, unsigned int sample_size = 100);
tfxINTERNAL void InitCommonParticleManager(tfx_particle_manager_t *pm, tfx_library_t *library, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit, tfx_particle_manager_mode mode, bool double_buffered_sprites, bool dynamic_sprite_allocation, tfxU32 mt_batch_size);
tfxINTERNAL bool ValidEffectID(tfx_particle_manager_t *pm, tfxEffectID id);

//--------------------------------
//Effect templates
//--------------------------------
tfxINTERNAL void AddTemplatePath(tfx_effect_template_t *effect_template, tfx_effect_emitter_t *effect_emitter, tfx_str256_t path);

//--------------------------------
//Library functions, internal/Editor functions
//--------------------------------
tfxAPI_EDITOR tfx_effect_emitter_info_t *GetEffectInfo(tfx_effect_emitter_t *e);					//Required by editor
tfxINTERNAL void PrepareLibraryEffectTemplate(tfx_library_t *library, tfx_str256_t path, tfx_effect_template_t *effect);
tfxINTERNAL void PrepareLibraryEffectTemplate(tfx_library_t *library, tfx_effect_emitter_t *effect, tfx_effect_template_t *effect_template);
//Copy the shape data to a memory location, like a staging buffer ready to be uploaded to the GPU for use in a compute shader
tfxINTERNAL void CopyLibraryLookupIndexesData(tfx_library_t *library, void* dst);
tfxINTERNAL void CopyLibraryLookupValuesData(tfx_library_t *library, void* dst);
tfxINTERNAL tfxU32 CountLibraryKeyframeLookUpValues(tfx_library_t *library, tfxU32 index);
tfxINTERNAL tfxU32 CountLibraryGlobalLookUpValues(tfx_library_t *library, tfxU32 index);
tfxINTERNAL tfxU32 CountLibraryEmitterLookUpValues(tfx_library_t *library, tfxU32 index);
tfxINTERNAL float LookupLibraryPreciseOvertimeNodeList(tfx_library_t *library, tfx_graph_type graph_type, int index, float age, float life);
tfxINTERNAL float LookupLibraryPreciseNodeList(tfx_library_t *library, tfx_graph_type graph_type, int index, float age);
tfxINTERNAL float LookupLibraryFastOvertimeValueList(tfx_library_t *library, tfx_graph_type graph_type, int index, float age, float life);
tfxINTERNAL float LookupLibraryFastValueList(tfx_library_t *library, tfx_graph_type graph_type, int index, float age);
tfxINTERNAL void InvalidateNewSpriteCapturedIndex(tfx_particle_manager_t *pm);
tfxINTERNAL void ResetSpriteDataLerpOffset(tfx_sprite_data_t *sprites);
tfxINTERNAL void CompressSpriteData(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect, bool is_3d, float frame_lengt);
tfxINTERNAL void LinkUpSpriteCapturedIndexes(tfx_work_queue_t *queue, void *data);
tfxINTERNAL void WrapSingleParticleSprites(tfx_sprite_data_t *sprite_data);
tfxINTERNAL void ClearWrapBit(tfx_sprite_data_t *sprite_data);
tfxINTERNAL void MaybeGrowLibraryInfos(tfx_library_t *library);

tfxAPI_EDITOR void MaybeGrowLibraryProperties(tfx_library_t *library, tfxU32 size_offset);	//Required by editor
tfxAPI_EDITOR tfxU32 GetLibraryComputeShapeDataSizeInBytes(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 GetLibraryComputeShapeCount(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 GetLibraryLookupIndexCount(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 GetLibraryLookupValueCount(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 GetLibraryLookupIndexesSizeInBytes(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 GetLibraryLookupValuesSizeInBytes(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 CountOfGraphsInUse(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 CountOfFreeGraphs(tfx_library_t *library);
tfxAPI_EDITOR bool IsLibraryShapeUsed(tfx_library_t *library, tfxKey image_hash);
tfxAPI_EDITOR bool LibraryShapeExists(tfx_library_t *library, tfxKey image_hash);
tfxAPI_EDITOR bool RemoveLibraryShape(tfx_library_t *library, tfxKey image_hash);
tfxAPI_EDITOR tfx_effect_emitter_t *InsertLibraryEffect(tfx_library_t *library, tfx_effect_emitter_t *effect, tfx_effect_emitter_t *position);
tfxAPI_EDITOR tfx_effect_emitter_t *AddLibraryEffect(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfx_effect_emitter_t *AddLibraryFolder(tfx_library_t *library, tfx_str64_t *name);
tfxAPI_EDITOR tfx_effect_emitter_t *AddLibraryFolder(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfx_effect_emitter_t *AddLibraryStage(tfx_library_t *library, tfx_str64_t *name);
tfxAPI_EDITOR void UpdateLibraryEffectPaths(tfx_library_t *library);
tfxAPI_EDITOR void AddLibraryPath(tfx_library_t *library, tfx_effect_emitter_t *effect_emitter, tfx_str256_t *path);
tfxAPI_EDITOR void DeleteLibraryEffect(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR bool RenameLibraryEffect(tfx_library_t *library, tfx_effect_emitter_t *effect, const char *new_name);
tfxAPI_EDITOR bool LibraryNameExists(tfx_library_t *library, tfx_effect_emitter_t *effect, const char *name);
tfxAPI_EDITOR void ReIndexLibrary(tfx_library_t *library);
tfxAPI_EDITOR void UpdateLibraryParticleShapeReferences(tfx_library_t *library, tfxKey default_hash);
tfxAPI_EDITOR tfx_effect_emitter_t* LibraryMoveUp(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfx_effect_emitter_t* LibraryMoveDown(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfxU32 AddLibraryGlobal(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 AddLibraryEmitterAttributes(tfx_library_t *library);
tfxAPI_EDITOR void FreeLibraryGlobal(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void FreeLibraryKeyframes(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void FreeLibraryEmitterAttributes(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void FreeLibraryProperties(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void FreeLibraryInfo(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR tfxU32 CloneLibraryGlobal(tfx_library_t *library, tfxU32 source_index, tfx_library_t *destination_library);
tfxAPI_EDITOR tfxU32 CloneLibraryKeyframes(tfx_library_t *library, tfxU32 source_index, tfx_library_t *destination_library);
tfxAPI_EDITOR tfxU32 CloneLibraryEmitterAttributes(tfx_library_t *library, tfxU32 source_index, tfx_library_t *destination_library);
tfxAPI_EDITOR tfxU32 CloneLibraryInfo(tfx_library_t *library, tfxU32 source_index, tfx_library_t *destination_library);
tfxAPI_EDITOR tfxU32 CloneLibraryProperties(tfx_library_t *library, tfx_emitter_properties_t *source, tfx_library_t *destination_library);
tfxAPI_EDITOR void AddLibraryEmitterGraphs(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void AddLibraryEffectGraphs(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void AddLibraryTransformGraphs(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfxU32 AddLibrarySpriteSheetSettings(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfxU32 AddLibrarySpriteDataSettings(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void AddLibrarySpriteSheetSettingsSub(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void AddLibrarySpriteDataSettingsSub(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfxU32 AddLibraryPreviewCameraSettings(tfx_library_t *library, tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfxU32 AddLibraryPreviewCameraSettings(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 AddLibraryEffectEmitterInfo(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 AddLibraryEmitterProperties(tfx_library_t *library);
tfxAPI_EDITOR tfxU32 AddLibraryKeyframes(tfx_library_t *library);
tfxAPI_EDITOR void UpdateLibraryComputeNodes(tfx_library_t *library);
tfxAPI_EDITOR void CompileAllLibraryGraphs(tfx_library_t *library);
tfxAPI_EDITOR void CompileLibraryGlobalGraph(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void CompileLibraryKeyframeGraph(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void CompileLibraryEmitterGraphs(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void CompileLibraryPropertyGraph(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void CompileLibraryBaseGraph(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void CompileLibraryVariationGraph(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void CompileLibraryOvertimeGraph(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void CompileLibraryColorGraphs(tfx_library_t *library, tfxU32 index);
tfxAPI_EDITOR void CompileLibraryGraphsOfEffect(tfx_library_t *library, tfx_effect_emitter_t *effect, tfxU32 depth = 0);
tfxAPI_EDITOR void SetLibraryMinMaxData(tfx_library_t *library);
tfxAPI_EDITOR void ClearLibrary(tfx_library_t *library);
tfxAPI_EDITOR void InitLibrary(tfx_library_t *library);
//Get an effect in the library by it's path. So for example, if you want to get a pointer to the emitter "spark" in effect "explosion" then you could do GetEffect("explosion/spark")
//You will need this function to apply user data and update callbacks to effects and emitters before adding the effect to the particle manager
//These are mainly for use by the editor, use effect templates instead, see PrepareEffectTemplate.
tfxAPI_EDITOR tfx_effect_emitter_t *GetLibraryEffect(tfx_library_t *library, const char *path);
//Get an effect by it's path hash key
tfxAPI_EDITOR tfx_effect_emitter_t *GetLibraryEffect(tfx_library_t *library, tfxKey key);
tfxAPI_EDITOR void RecordSpriteData(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect, float update_frequency, float camera_position[3]);

//Effect/Emitter functions
void SetEffectUserData(tfx_effect_emitter_t *e, void *data);
void *GetEffectUserData(tfx_effect_emitter_t *e);

tfxINTERNAL bool IsRootEffect(tfx_effect_emitter_t *effect);
tfxINTERNAL void ResetEffectParents(tfx_effect_emitter_t *effect);
tfxINTERNAL void CompileEffectGraphs(tfx_effect_emitter_t *effect);
tfxINTERNAL void FreeEffectGraphs(tfx_effect_emitter_t *effect);
tfxINTERNAL tfxU32 CountAllEffectLookupValues(tfx_effect_emitter_t *effect);
tfxINTERNAL float GetEffectLoopLength(tfx_effect_emitter_t *effect);

tfxAPI_EDITOR tfx_emitter_properties_t *GetEffectProperties(tfx_effect_emitter_t *e);
tfxAPI_EDITOR tfx_effect_emitter_t* AddEmitterToEffect(tfx_effect_emitter_t *effect, tfx_effect_emitter_t *e);
tfxAPI_EDITOR tfx_effect_emitter_t* AddEffectToEmitter(tfx_effect_emitter_t *effect, tfx_effect_emitter_t *e);
tfxAPI_EDITOR tfx_effect_emitter_t* AddEffect(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR int GetEffectDepth(tfx_effect_emitter_t *e);
tfxAPI_EDITOR tfxU32 CountAllEffects(tfx_effect_emitter_t *effect, tfxU32 amount);
tfxAPI_EDITOR tfx_effect_emitter_t* tfx_GetRootEffect(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void ReIndexEffect(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void CountEffectChildren(tfx_effect_emitter_t *effect, int *emitters, int *effects);
tfxAPI_EDITOR tfx_effect_emitter_t* MoveEffectUp(tfx_effect_emitter_t *effect_to_move);
tfxAPI_EDITOR tfx_effect_emitter_t* MoveEffectDown(tfx_effect_emitter_t *effect_to_move);
tfxAPI_EDITOR void DeleteEmitterFromEffect(tfx_effect_emitter_t *emitter_to_delete);
tfxAPI_EDITOR void CleanUpEffect(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void ResetEffectGraphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxAPI_EDITOR void ResetTransformGraphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxAPI_EDITOR void ResetEmitterBaseGraphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxAPI_EDITOR void ResetEmitterPropertyGraphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxAPI_EDITOR void ResetEmitterVariationGraphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxAPI_EDITOR void ResetEmitterOvertimeGraphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);
tfxAPI_EDITOR void ResetEmitterGraphs(tfx_effect_emitter_t *effect, bool add_node = true, bool compile = true);

tfxAPI_EDITOR void AddEmitterColorOvertime(tfx_effect_emitter_t *effect, float frame, tfx_rgb_t color);
tfxAPI_EDITOR void UpdateEffectMaxLife(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfx_graph_t* GetEffectGraphByType(tfx_effect_emitter_t *effect, tfx_graph_type type);
tfxAPI_EDITOR tfxU32 GetEffectGraphIndexByType(tfx_effect_emitter_t *effect, tfx_graph_type type);
tfxAPI_EDITOR void InitialiseUninitialisedGraphs(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void SetEffectName(tfx_effect_emitter_t *effect, const char *n);
tfxAPI_EDITOR bool RenameSubEffector(tfx_effect_emitter_t *effect, const char *new_name);
tfxAPI_EDITOR bool EffectNameExists(tfx_effect_emitter_t *in_effect, tfx_effect_emitter_t *excluding_effect, const char *name);
tfxAPI_EDITOR void CloneEffect(tfx_effect_emitter_t *effect_to_clone, tfx_effect_emitter_t *clone, tfx_effect_emitter_t *root_parent, tfx_library_t *destination_library, tfxEffectCloningFlags flags = 0);
tfxAPI_EDITOR void EnableAllEmitters(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void EnableEmitter(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void DisableAllEmitters(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void DisableAllEmittersExcept(tfx_effect_emitter_t *effect, tfx_effect_emitter_t *emitter);
tfxAPI_EDITOR bool IsFiniteEffect(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR void FlagEffectAs3D(tfx_effect_emitter_t *effect, bool flag);
tfxAPI_EDITOR bool Is3DEffect(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfx_particle_manager_mode GetRequiredParticleManagerMode(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR tfx_preview_camera_settings_t *GetEffectCameraSettings(tfx_effect_emitter_t *effect);
tfxAPI_EDITOR float GetEffectHighestLoopLength(tfx_effect_emitter_t *effect);

//------------------------------------------------------------
//Section API_Functions
//------------------------------------------------------------

//[API functions]
//All the functions below represent all that you will need to call to implement TimelineFX

/*
You don't have to call this, you can just call InitialiseTimelineFX in order to initialise the memory, but I created this for the sake of the editor which
needs to load in an ini file before initialising timelinefx which requires the memory pool to be created before hand
* @param memory_pool_size	The size of each memory pool to contain all objects created in TimelineFX
*/
tfxAPI void InitialiseTimelineFXMemory(size_t memory_pool_size = tfxMegabyte(128));

/*
Initialise TimelineFX. Must be called before any functionality of TimelineFX is used.
* @param max_threads	Pass the number of threads that you want to use in addition to the main thread.
*						Example, if there are 12 logical cores available, 0.5 will use 6 threads. 0 means only single threaded will be used.
*/
tfxAPI void InitialiseTimelineFX(int max_threads = 0, size_t memory_pool_size = tfxMegabyte(128));

/*
Initialise TimelineFX. Must be called before any functionality of TimelineFX is used.
* @param filename		The name of the file where you want to count the number of shapes
* @returns int			The number of shapes in the library.
*/
tfxAPI int GetShapeCountInLibrary(const char *filename);

/*
Validate a timelinefx tfx file to make sure that it's valid.
* @param filename		The name of the file where you want to count the number of shapes
* @returns int			Returns 0 if the file successfully validated or a tfxErrorFlags if something went wrong
*/
tfxAPI int ValidateEffectPackage(const char *filename);

/**
* Loads an effect library package from the specified filename into the provided tfx_library_t object.
*
* @param filename		A pointer to a null-terminated string that contains the path and filename of the effect library package to be loaded.
* @param lib			A reference to a tfx_library_t object that will hold the loaded effect library data.
* @param shape_loader	A pointer to a function that will be used to load image data into the effect library package.
*						The function has the following signature: void shape_loader(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data).
* @param user_data		A pointer to user-defined data that will be passed to the shape_loader function. This parameter is optional and can be set to nullptr if not needed.
* @param read_only		A boolean value that determines whether the effect library data will be loaded in read-only mode. (Maybe removed in the future).
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
tfxAPI tfxErrorFlags LoadEffectLibraryPackage(const char *filename, tfx_library_t *lib, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void *user_data = nullptr, bool read_only = true);

/**
* Loads a sprite data file into an animation manager
*
* @param filename		A pointer to a null-terminated string that contains the path and filename of the effect library package to be loaded.
* @param lib			A reference to a tfx_animation_manager_t object that will hold the loaded sprite data.
* @param shape_loader	A pointer to a function that will be used to load image data into the effect library package.
*						The function has the following signature: void shape_loader(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data).
* @param user_data		A pointer to user-defined data that will be passed to the shape_loader function. This parameter is optional and can be set to nullptr if not needed.
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
tfxAPI tfxErrorFlags LoadSpriteData(const char *filename, tfx_animation_manager_t *animation_manager, void(*shape_loader)(const char *filename, tfx_image_data_t *image_data, void *raw_image_data, int image_size, void *user_data), void *user_data = nullptr);

/*
Output all the effect names in a library to the console
* @param tfx_library_t				A valid pointer to a tfx_library_t
*/
inline tfxAPI void ListEffectNames(tfx_library_t *library) {
	tfxU32 index = 0;
	for (auto &effect : library->effects) {
		printf("%i) %s\n", index++, GetEffectInfo(&effect)->name.c_str());
	}
}

//[Particle Manager functions]

/*
Initialise a tfx_particle_manager_t for 3d usage
* @param pm						A pointer to an unitialised tfx_particle_manager_t. If you want to reconfigure a particle manager for a different usage then you can call ReconfigureParticleManager.
* @param library				A pointer to a tfx_library_t that you will be using to add all of the effects from to the particle manager.
* @param layer_max_values		An array of unsigned ints representing the maximum amount of particles you want available for each layer. This will allocate the appropriate amount of memory ahead of time.
* @param effects_limit			The maximum amount of effects and emitters that can be updated in a single frame. This will allocate the appropriate amount of memory ahead of time. Default: 1000.
* @param mode					The operation mode of the particle manager regarding how particles are ordered. Default value: tfxParticleManagerMode_unordered. Possible modes are:
	tfxParticleManagerMode_unordered					Particles will be updated by emitter. No ordering is maintained, each emitter will spawn and update their particles in turn and sprites will be ordered
														according to that sequence.
	tfxParticleManagerMode_ordered_by_age				Particles will be kept in age order, older particles will be drawn first and newer ones last
	tfxParticleManagerMode_ordered_by_depth				Particles will be drawn in depth order or distance from the camera. You can specify the number of sort passes when setting up the effects in TimelineFX editor
	tfxParticleManagerMode_ordered_by_depth_guaranteed	Particles will be sorted each update and kept in depth order
* @param double_buffer_sprites	True or False, whether the last frame of sprites is kept so that you can use to do interpolations for smoother animation
* @param dynamic_allocation		If set to true then when the layer_max_values is hit for a layer the sprite and particle memory allocation will be grown dynamically. This can be useful when you're unsure of how
								many particles you will need to display while developing you're game/app. Default is false.
* @param mt_batch_size			When using multithreading you can alter the size of each batch of particles that each thread will update. The default is 512

*/
tfxAPI void InitParticleManagerFor3d(tfx_particle_manager_t *pm, tfx_library_t *library, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit = 1000, tfx_particle_manager_mode mode = tfxParticleManagerMode_unordered, bool double_buffer_sprites = true, bool dynamic_allocation = false, tfxU32 mt_batch_size = 512);

/*
Initialise a tfx_particle_manager_t for 2d usage
* @param pm						A pointer to an unitialised tfx_particle_manager_t. If you want to reconfigure a particle manager for a different usage then you can call ReconfigureParticleManager.
* @param library				A pointer to a tfx_library_t that you will be using to add all of the effects from to the particle manager.
* @param layer_max_values		An array of unsigned ints representing the maximum amount of particles you want available for each layer. This will allocate the appropriate amount of memory ahead of time.
* @param effects_limit			The maximum amount of effects and emitters that can be updated in a single frame. This will allocate the appropriate amount of memory ahead of time. Default: 1000.
* @param mode					The operation mode of the particle manager regarding how particles are ordered. Default value: tfxParticleManagerMode_unordered. Possible modes are:
	tfxParticleManagerMode_unordered					Particles will be updated by emitter. No ordering is maintained, each emitter will spawn and update their particles in turn and sprites will be ordered
														according to that sequence.
	tfxParticleManagerMode_ordered_by_age				Particles will be kept in age order, older particles will be drawn first and newer ones last
* @param double_buffer_sprites	True or False, whether the last frame of sprites is kept so that you can use to do interpolations for smoother animation
* @param dynamic_allocation		If set to true then when the layer_max_values is hit for a layer the sprite and particle memory allocation will be grown dynamically. This can be useful when you're unsure of how
								many particles you will need to display while developing you're game/app. Default is false.
* @param mt_batch_size			When using multithreading you can alter the size of each batch of particles that each thread will update. The default is 512.

*/
tfxAPI void InitParticleManagerFor2d(tfx_particle_manager_t *pm, tfx_library_t *library, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit = 1000, tfx_particle_manager_mode mode = tfxParticleManagerMode_unordered, bool double_buffer_sprites = true, bool dynamic_allocation = false, tfxU32 mt_batch_size = 512);

/*
Initialise a tfx_particle_manager_t for both 2d and 3d. This just allocates buffers for both 2d and 3d anticipating that you'll be using ReconfigureParticleManager to switch between 2d/3d modes. If you want to update
both 2d and 3d particles at the same time then just use 2 separate particle managers instead as a particle manager can only update one type of particle 2d or 3d.
* @param pm						A pointer to an unitialised tfx_particle_manager_t. If you want to reconfigure a particle manager for a different usage then you can call ReconfigureParticleManager.
* @param library				A pointer to a tfx_library_t that you will be using to add all of the effects from to the particle manager.
* @param layer_max_values		An array of unsigned ints representing the maximum amount of particles you want available for each layer. This will allocate the appropriate amount of memory ahead of time.
* @param effects_limit			The maximum amount of effects and emitters that can be updated in a single frame. This will allocate the appropriate amount of memory ahead of time. Default: 1000.
* @param mode					The operation mode of the particle manager regarding how particles are ordered. Default value: tfxParticleManagerMode_unordered. Possible modes are:
	tfxParticleManagerMode_unordered					Particles will be updated by emitter. No ordering is maintained, each emitter will spawn and update their particles in turn and sprites will be ordered
														according to that sequence.
	tfxParticleManagerMode_ordered_by_age				Particles will be kept in age order, older particles will be drawn first and newer ones last
* @param double_buffer_sprites	True or False, whether the last frame of sprites is kept so that you can use to do interpolations for smoother animation
* @param dynamic_allocation		If set to true then when the layer_max_values is hit for a layer the sprite and particle memory allocation will be grown dynamically. This can be useful when you're unsure of how
								many particles you will need to display while developing you're game/app. Default is false.
* @param mt_batch_size			When using multithreading you can alter the size of each batch of particles that each thread will update. The default is 512.

*/
tfxAPI void InitParticleManagerForBoth(tfx_particle_manager_t *pm, tfx_library_t *library, tfxU32 layer_max_values[tfxLAYERS], unsigned int effects_limit = 1000, tfx_particle_manager_mode mode = tfxParticleManagerMode_unordered, bool double_buffer_sprites = true, bool dynamic_sprite_allocation = false, tfxU32 multi_threaded_batch_size = 512);

/*
Reconfigure a particle manager to make it work in a different mode. A particle manager can only run in a single mode at time like unordered, depth ordered etc so use this to change that. Also bear
in mind that you can just use more than one particle manager and utilised different modes that way as well. The modes that you need will depend on the effects that you're adding to the particle manager.
* @param pm						A pointer to an intialised tfx_particle_manager_t.
* @param mode					One of the following modes:
								tfxParticleManagerMode_unordered
								tfxParticleManagerMode_ordered_by_age
								tfxParticleManagerMode_ordered_by_depth
								tfxParticleManagerMode_ordered_by_depth_guaranteed
* @param sort_passes			The number of sort passes if you're using depth sorted effects
* @param is_3d					True if the particle manager should be configured for 3d effects.
*/
void ReconfigureParticleManager(tfx_particle_manager_t *pm, tfx_particle_manager_mode mode, tfxU32 sort_passes, bool is_3d);

/*
When a particle manager updates particles it creates work queues to handle the work. By default these each have a maximum amount of 1000 entries which should be
more than enough for most situations. However you can increase the sizes here if needed. You only need to set this manually if you hit one of the asserts when these
run out of space or you anticipate a huge amount of emitters and particles to be used (> million). On the other hand, you might be tight on memory in which case you
could reduce the numbers as well if needed (they don't take a lot of space though)
* @param pm						A pointer to an intialised tfx_particle_manager_t.
* @param spawn_work_max			The maximum amount of spawn work entries
* @param control_work_max		The maximum amount of control work entries
* @param age_work_max			The maximum amount of age_work work entries
*/
void SetPMWorkQueueSizes(tfx_particle_manager_t *pm, tfxU32 spawn_work_max = 1000, tfxU32 control_work_max = 1000, tfxU32 age_work_max = 1000);

/*
Get the current particle count for a particle manager
* @param pm						A pointer to an tfx_particle_manager_t
* @returns tfxU32				The total number of particles currently being updated
*/
tfxU32 ParticleCount(tfx_particle_manager_t *pm);

/*
Get the current number of effects that are currently being updated by a particle manager
* @param pm						A pointer to an tfx_particle_manager_t
* @returns tfxU32				The total number of effects currently being updated
*/
tfxU32 EffectCount(tfx_particle_manager_t *pm);

/*
Get the current number of emitters that are currently being updated by a particle manager
* @param pm						A pointer to an tfx_particle_manager_t
* @returns tfxU32				The total number of emitters currently being updated
*/
tfxU32 EmitterCount(tfx_particle_manager_t *pm);

/*
Set the seed for the particle manager for random number generation. Setting the seed can determine how an emitters spawns particles, so if you set the seed before adding an effect to the particle manager
then the effect will look the same each time. Note that seed of 0 is invalid, it must be 1 or greater.
* @param pm							A pointer to an initialised tfx_particle_manager_t. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
* @param seed						An unsigned int representing the seed (Any value other then 0)
*/
tfxAPI inline void SetSeed(tfx_particle_manager_t *pm, tfxU64 seed) {
	RandomReSeed(&pm->random, seed == 0 ? tfxMAX_UINT : seed);
}

/*
Prepare a tfx_effect_template_t that you can use to customise effects in the library in various ways before adding them into a particle manager for updating and rendering. Using a template like this
means that you can tweak an effect without editing the base effect in the library.
* @param library					A reference to a tfx_library_t that should be loaded with LoadEffectLibraryPackage
* @param name						The name of the effect in the library that you want to use for the template. If the effect is in a folder then use normal pathing: "My Folder/My effect"
* @param effect_template			The empty tfx_effect_template_t object that you want the effect loading into
//Returns true on success.
*/
tfxAPI bool PrepareEffectTemplate(tfx_library_t *library, const char *name, tfx_effect_template_t *effect_template);

/*
Add an effect to a tfx_particle_manager_t from an effect template
* @param pm					A pointer to an initialised tfx_particle_manager_t. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
* @param effect_template	The tfx_effect_template_t object that you want to add to the particle manager. It must have already been prepared by calling PrepareEffectTemplate
* @param effect_id			pointer to a tfxEffectID of the effect which will be set after it's been added to the particle manager. This index can then be used to manipulate the effect in the particle manager as it's update
							For example by calling SetEffectPosition. This will be set to tfxINVALID if the function is unable to add the effect to the particle manager if it's out of space and reached it's effect limit.
  @returns					True if the effect was succesfully added.
*/
tfxAPI bool AddEffectToParticleManager(tfx_particle_manager_t *pm, tfx_effect_template_t *effect, tfxEffectID *effect_id = nullptr);

/*
Add an effect to a tfx_particle_manager_t.
* @param pm					A pointer to an initialised tfx_particle_manager_t. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
* @param effect				tfx_effect_emitter_t object that you want to add to the particle manager.
* @param effect_id			pointer to a tfxEffectID of the effect which will be set after it's been added to the particle manager. This index can then be used to manipulate the effect in the particle manager as it's update
							For example by calling SetEffectPosition. This will be set to tfxINVALID if the function is unable to add the effect to the particle manager if it's out of space and reached it's effect limit.
  @returns					True if the effect was succesfully added.
*/
tfxAPI bool AddEffectToParticleManager(tfx_particle_manager_t *pm, tfx_effect_emitter_t *effect, tfxEffectID *effect_id = nullptr);

/*
Update a particle manager. Call this function each frame in your update loop. It should be called the same number of times per second as set with SetUpdateFrequency.
* @param pm					A pointer to an initialised tfx_particle_manager_t. The particle manager must have already been initialised by calling InitFor3d or InitFor2d
*/
tfxAPI void UpdateParticleManager(tfx_particle_manager_t *pm, float elapsed);

/*
Get the total number of sprites within the layer of the particle manager
* @param pm					A pointer to an initialised tfx_particle_manager_t.
* @param layer				The layer of the sprites to the count of
*/
tfxAPI inline tfxU32 SpritesInLayerCount(tfx_particle_manager_t *pm, tfxU32 layer) {
	return pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
}

/*
Get the total number of sprites within the layer of the particle manager
* @param pm					A pointer to an initialised tfx_particle_manager_t.
* @param layer				The layer of the sprites to the count of
*/
tfxAPI inline tfx_sprite_soa_t *SpritesInLayer(tfx_particle_manager_t *pm, tfxU32 layer) {
	return &pm->sprites[pm->current_sprite_buffer][layer];
}

/*
Get the total number of 3d sprites ready for rendering in the particle manager
* @param pm					A pointer to an initialised tfx_particle_manager_t.
*/
tfxAPI inline tfxU32 TotalSpriteCount(tfx_particle_manager_t *pm) {
	tfxU32 count = 0;
	for (tfxEachLayer) {
		count += pm->sprite_buffer[pm->current_sprite_buffer][layer].current_size;
	}
	return count;
}

/*
Clear all particles, sprites and effects in a particle manager. If you don't need to use the particle manager again then call FreeParticleManager to also
free all the memory associated with the particle manager.
* @param pm						A pointer to an initialised tfx_particle_manager_t.
* @param free_particle_banks	Set to true if you want to free the memory associated with the particle banks and release back to the memory pool
*/
tfxAPI void ClearParticleManager(tfx_particle_manager_t *pm, bool free_particle_banks);

/*
Free all the memory used in the particle manager.
* @param pm						A pointer to an initialised tfx_particle_manager_t.
*/
tfxAPI void FreeParticleManager(tfx_particle_manager_t *pm);

//[Effects functions for altering effects that are currently playing out in a particle manager]

/*
Expire an effect by telling it to stop spawning particles. This means that the effect will eventually be removed from the particle manager after all of it's remaining particles have expired.
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index	The index of the effect that you want to expire. This is the index returned when calling AddEffectToParticleManager
*/
tfxAPI inline void SoftExpireEffect(tfx_particle_manager_t *pm, tfxEffectID effect_index) {
	pm->effects[effect_index].state_flags |= tfxEmitterStateFlags_stop_spawning;
}

/*
Soft expire all the effects in a particle manager so that the particles complete their animation first
*/
tfxAPI void SoftExpireAll(tfx_particle_manager_t *pm);

/*
Expire an effect by telling it to stop spawning particles and remove all associated particles immediately.
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index	The index of the effect that you want to expire. This is the index returned when calling AddEffectToParticleManager
*/
tfxAPI inline void HardExpireEffect(tfx_particle_manager_t *pm, tfxEffectID effect_index) {
	pm->effects[effect_index].state_flags |= tfxEmitterStateFlags_stop_spawning;
	pm->effects[effect_index].state_flags |= tfxEmitterStateFlags_remove;
}

/*
Get effect user data
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index	The index of the effect that you want to expire. This is the index returned when calling AddEffectToParticleManager
* @returns				void* pointing to the user data set in the effect. See tfx_effect_template_t::SetUserData() and SetEffectUserData()
*/
tfxAPI inline void* GetEffectUserData(tfx_particle_manager_t *pm, tfxEffectID effect_index) {
	return pm->effects[effect_index].user_data;
}

/*
More for use in the editor, this function updates emitter base values for any effects that are currently running after their graph values have been changed.
*/
tfxAPI void UpdatePMBaseValues(tfx_particle_manager_t *pm);

/*
Set the tfx_library_t that the particle manager will use to render sprites and lookup all of the various properties required to update emitters and particles.
This is also set when you initialise a particle manager
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param lib			A pointer to a tfx_library_t
*/
tfxAPI void SetPMLibrary(tfx_particle_manager_t *pm, tfx_library_t *library);

/*
Set the particle manager camera. This is used to calculate particle depth if you're using depth ordered particles so it needs to be updated each frame.
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param front			An array of 3 floats representing a normalised 3d vector describing the direction that the camera is pointing
* @param position		An array of 3 floats representing the position of the camera in 3d space
*/
tfxAPI void SetPMCamera(tfx_particle_manager_t *pm, float front[3], float position[3]);

/*
Set the lookup mode for the particle manager. Lookup modes are either calculating the emitter attribute graphs in real time, or looking up arrays of the graphs
in pre-compiled lookup tables. Generally tfxFast will be a little bit faster due to less math involved. Having said that more graph nodes could be cached in memory
so maybe the extra math wouldn't make much difference but either way it needs more testing and profiling!
todo: callbacks should be moved into the particle manager, currently they're global callbacks so having this function is a bit pointless!
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param mode			The look up mode you want to set. tfxFast is the default mode.
*/
tfxAPI void SetPMLookUpMode(tfx_particle_manager_t *pm, tfx_lookup_mode mode);

/*
Each effect in the particle manager can have bounding box which you can decide to keep updated or not if you wanted to do any offscreen culling of effects. Theres some
extra overhead to keep the bounding boxes updated but that can be made back if you have a number of effect particles offscreen that don't need to be drawn.
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param yesno			Set to true or false if you want the bounding boxes to be udpated.
*/
tfxAPI void KeepBoundingBoxesUpdated(tfx_particle_manager_t *pm, bool yesno);

/*
Set the effect user data for an effect already added to a particle manager
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index	The index of the effect that you want to expire. This is the index returned when calling AddEffectToParticleManager
* @param user_data		A void* pointing to the user_data that you want to store in the effect
*/
tfxAPI inline void SetEffectUserData(tfx_particle_manager_t *pm, tfxEffectID effect_index, void* user_data) {
	pm->effects[effect_index].user_data = user_data;
}

/*
Force a particle manager to only run in single threaded mode. In other words, only use the main thread to update particles
* @param pm				A pointer to a tfx_particle_manager_t.
* @param switch_on		true or false to use a single thread or not
*/
tfxAPI inline void ForcePMSingleThreaded(tfx_particle_manager_t *pm, bool switch_on) {
	if (switch_on) pm->flags |= tfxEffectManagerFlags_single_threaded; else pm->flags &= ~tfxEffectManagerFlags_single_threaded;
}

/*
Get the transform vectors for a 3d sprite's previous position so that you can use that to interpolate between that and the current sprite position
* @param pm				A pointer to a tfx_particle_manager_t.
* @param layer			The index of the sprite layer
* @param index			The sprite index of the sprite that you want the captured sprite for.
*/
tfxAPI inline tfx_sprite_transform3d_t *GetCapturedSprite3dTransform(tfx_particle_manager_t *pm, tfxU32 layer, tfxU32 index) {
	return &pm->sprites[(index & 0x40000000) >> 30][layer].transform_3d[index & 0x0FFFFFFF];
}

/*
Get the transform vectors for a 2d sprite's previous position so that you can use that to interpolate between that and the current sprite position
* @param pm				A pointer to a tfx_particle_manager_t.
* @param layer			The index of the sprite layer
* @param index			The sprite index of the sprite that you want the captured sprite for.
*/
tfxAPI inline tfx_sprite_transform2d_t *GetCapturedSprite2dTransform(tfx_particle_manager_t *pm, tfxU32 layer, tfxU32 index) {
	return &pm->sprites[(index & 0x40000000) >> 30][layer].transform_2d[index & 0x0FFFFFFF];
}

/*
Get the intensity for a sprite's previous frame so that you can use that to interpolate between that and the current sprite intensity
* @param pm				A pointer to a tfx_particle_manager_t.
* @param layer			The index of the sprite layer
* @param index			The sprite index of the sprite that you want the captured sprite for.
*/
tfxAPI inline float *GetCapturedSprite3dIntensity(tfx_particle_manager_t *pm, tfxU32 layer, tfxU32 index) {
	return &pm->sprites[(index & 0x40000000) >> 30][layer].intensity[index & 0x0FFFFFFF];
}

/*
Get the index offset into the sprite memory for sprite data containing a pre recorded effect animation. Can be used along side SpriteDataEndIndex to create
a for loop to iterate over the sprites in a pre-recorded effect
* @param sprite_data	A pointer to tfx_sprite_data_t containing all the sprites and frame data
* @param frame			The index of the frame you want the offset for
* @param layer			The sprite layer
* @returns				tfxU32 containing the index offset
*/
tfxAPI inline tfxU32 SpriteDataIndexOffset(tfx_sprite_data_t *sprite_data, tfxU32 frame, tfxU32 layer) {
	assert(frame < sprite_data->normal.frame_meta.size());			//frame is outside index range
	assert(layer < tfxLAYERS);								//layer is outside index range
	return sprite_data->normal.frame_meta[frame].index_offset[layer];
}

/*
Make a particle manager stop spawning. This will mean that all emitters in the particle manager will no longer spawn any particles so all currently running effects will expire
as the remaining particles come to the end of their life. Any single particles will also get flagged to expire
* @param pm				A pointer to a tfx_particle_manager_t.
* @param yesno			True = disable spawning, false = enable spawning
*/
tfxAPI inline void DisablePMSpawning(tfx_particle_manager_t *pm, bool yesno) {
	if (yesno) {
		pm->flags |= tfxEffectManagerFlags_disable_spawning;
	}
	else {
		pm->flags &= ~tfxEffectManagerFlags_disable_spawning;
	}
}

/*
Get the buffer of effect indexes in the particle manager.
* @param pm				A pointer to a tfx_particle_manager_t.
* @param depth			The depth of the list that you want. 0 are top level effects and anything higher are sub effects within those effects
* @returns				Pointer to the tfxvec of effect indexes
*/
tfxAPI tfx_vector_t<tfxU32> *GetPMEffectBuffer(tfx_particle_manager_t *pm, tfxU32 depth);

/*
Get the buffer of emitter indexes in the particle manager.
* @param pm				A pointer to a tfx_particle_manager_t.
* @param depth			The depth of the list that you want. 0 are top level emitters and anything higher are sub emitters within those effects
* @returns				Pointer to the tfxvec of effect indexes
*/
tfxAPI tfx_vector_t<tfxU32> *GetPMEmitterBuffer(tfx_particle_manager_t *pm, tfxU32 depth);

/*
Get the end index offset into the sprite memory for sprite data containing a pre recorded effect animation that has been compresssed into fewer frames. Can be used along side SpriteDataIndexOffset to create
a for loop to iterate over the sprites in a pre-recorded effect
* @param sprite_data	A pointer to tfx_sprite_data_t containing all the sprites and frame data
* @param frame			The index of the frame you want the offset for
* @param layer			The sprite layer
* @returns				tfxU32 containing the index offset
*/
tfxAPI inline tfxU32 CompressedSpriteDataIndexOffset(tfx_sprite_data_t *sprite_data, tfxU32 frame, tfxU32 layer) {
	assert(frame < sprite_data->compressed.frame_meta.size());			//frame is outside index range
	assert(layer < tfxLAYERS);								//layer is outside index range
	return sprite_data->compressed.frame_meta[frame].index_offset[layer];
}

/*
Get the index offset into the sprite memory for sprite data containing a pre recorded effect animation. Can be used along side SpriteDataEndIndex to create
a for loop to iterate over the sprites in a pre-recorded effect
* @param sprite_data	A pointer to tfx_sprite_data_t containing all the sprites and frame data
* @param frame			The index of the frame you want the end index for
* @param layer			The sprite layer
* @returns				tfxU32 containing the end offset
*/
tfxAPI inline tfxU32 SpriteDataEndIndex(tfx_sprite_data_t *sprite_data, tfxU32 frame, tfxU32 layer) {
	assert(frame < sprite_data->normal.frame_meta.size());			//frame is outside index range
	assert(layer < tfxLAYERS);								//layer is outside index range
	return sprite_data->normal.frame_meta[frame].index_offset[layer] + sprite_data->normal.frame_meta[frame].sprite_count[layer];
}

/*
Get the end index offset into the sprite memory for sprite data containing a pre recorded effect animation that has been compressed into fewer frames. Can be used along side CompressedSpriteDataIndexOffset to create
a for loop to iterate over the sprites in a pre-recorded effect
* @param sprite_data	A pointer to tfx_sprite_data_t containing all the sprites and frame data
* @param frame			The index of the frame you want the end index for
* @param layer			The sprite layer
* @returns				tfxU32 containing the end offset
*/
tfxAPI inline tfxU32 CompressedSpriteDataEndIndex(tfx_sprite_data_t *sprite_data, tfxU32 frame, tfxU32 layer) {
	assert(frame < sprite_data->compressed.frame_meta.size());			//frame is outside index range
	assert(layer < tfxLAYERS);								//layer is outside index range
	return sprite_data->compressed.frame_meta[frame].index_offset[layer] + sprite_data->compressed.frame_meta[frame].sprite_count[layer];
}

/*
Get the 3d transform struct of a sprite data by its index in the sprite data struct of arrays
* @param sprite_data	A pointer to tfx_sprite_data_t containing all the sprites and frame data
* @param index			The index of the sprite you want to retrieve
* @returns				tfx_sprite_transform3d_t reference
*/
tfxAPI inline tfx_sprite_transform3d_t *GetSpriteData3dTransform(tfx_sprite_data_soa_t *sprites, tfxU32 index) {
	return &sprites->transform_3d[index];
}

/*
Get the 2d transform struct of a sprite data by its index in the sprite data struct of arrays
* @param sprite_data	A pointer to tfx_sprite_data_t containing all the sprites and frame data
* @param index			The index of the sprite you want to retrieve
* @returns				tfx_sprite_transform2d_t reference
*/
tfxAPI inline tfx_sprite_transform2d_t *GetSpriteData2dTransform(tfx_sprite_data_soa_t *sprites, tfxU32 index) {
	return &sprites->transform_2d[index];
}

/*
Get the intensity of a sprite data by its index in the sprite data struct of arrays
* @param sprite_data	A pointer to tfx_sprite_data_t containing all the sprites and frame data
* @param index			The index of the sprite you want to retrieve
* @returns				float of the intensity value
*/
tfxAPI inline float GetSpriteDataIntensity(tfx_sprite_data_soa_t *sprites, tfxU32 index) {
	return sprites->intensity[index];
}

/*
Get the alignment of a sprite data by its index in the sprite data struct of arrays
* @param sprite_data	A pointer to tfx_sprite_data_t containing all the sprites and frame data
* @param index			The index of the sprite you want to retrieve
* @returns				tfxU32 
*/
tfxAPI inline tfxU32 GetSpriteDataAlignment(tfx_sprite_data_soa_t *sprites, tfxU32 index) {
	return sprites->alignment[index];
}

/*
Get the image frame of a sprite data by its index in the sprite data struct of arrays
* @param sprite_data	A pointer to tfx_sprite_data_t containing all the sprites and frame data
* @param index			The index of the sprite you want to retrieve
* @returns				tfxU32 of the frame value
*/
tfxAPI inline tfxU32 GetSpriteDataFrame(tfx_sprite_data_soa_t *sprites, tfxU32 index) {
	return (sprites->property_indexes[index] & 0x00FF0000) >> 16;
}

/*
Get the color of a sprite data by its index in the sprite data struct of arrays
* @param sprite_data	A pointer to tfx_sprite_data_t containing all the sprites and frame data
* @param index			The index of the sprite you want to retrieve
* @returns				tfx_rgba8_t of the frame value
*/
tfxAPI inline tfx_rgba8_t GetSpriteDataColor(tfx_sprite_data_soa_t *sprites, tfxU32 index) {
	return sprites->color[index];
}

/*
Set the position of a 2d effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param x				The x value of the position
* @param y				The y value of the position
*/
tfxAPI void SetEffectPosition(tfx_particle_manager_t *pm, tfxEffectID effect_index, float x, float y);

/*
Set the position of a 3d effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param x				The x value of the position
* @param y				The y value of the position
* @param z				The z value of the position
*/
tfxAPI void SetEffectPosition(tfx_particle_manager_t *pm, tfxEffectID effect_index, float x, float y, float z);

/*
Set the position of a 2d effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param position		A tfx_vec2_t vector object containing the x and y coordinates
*/
tfxAPI void SetEffectPosition(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfx_vec2_t position);

/*
Set the position of a 3d effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param position		A tfx_vec3_t vector object containing the x, y and z coordinates
*/
tfxAPI void SetEffectPosition(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfx_vec3_t position);

/*
Move an Effect by a specified amount relative to the effect's current position
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param amount			A tfx_vec3_t vector object containing the amount to move in the x, y and z planes
*/
tfxAPI void MoveEffect(tfx_particle_manager_t *pm, tfxEffectID effect_index, tfx_vec3_t amount);

/*
Move an Effect by a specified amount relative to the effect's current position
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param x				The amount to move in the x plane
* @param y				The amount to move in the y plane
* @param z				The amount to move in the z plane
*/
tfxAPI void MoveEffect(tfx_particle_manager_t *pm, tfxEffectID effect_index, float x, float y, float z);

/*
Get the current position of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @return				tfx_vec3_t containing the effect position
*/
tfxAPI tfx_vec3_t GetEffectPosition(tfx_particle_manager_t *pm, tfxEffectID effect_index);

/*
Set the rotation of a 2d effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current rotation of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param rotation		A float of the amount that you want to set the rotation too
*/
tfxAPI void SetEffectRotation(tfx_particle_manager_t *pm, tfxEffectID effect_index, float rotation);

/*
Set the roll of a 3d effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current roll of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param roll			A float of the amount that you want to set the roll too
*/
tfxAPI void SetEffectRoll(tfx_particle_manager_t *pm, tfxEffectID effect_index, float roll);

/*
Set the pitch of a 3d effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current pitch of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param pitch			A float of the amount that you want to set the pitch too
*/
tfxAPI void SetEffectPitch(tfx_particle_manager_t *pm, tfxEffectID effect_index, float pitch);

/*
Set the yaw of a 3d effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current yaw of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param yaw			A float of the amount that you want to set the yaw too
*/
tfxAPI void SetEffectYaw(tfx_particle_manager_t *pm, tfxEffectID effect_index, float yaw);

/*
Set the width of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current width of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param width			A float of the amount that you want to set the width multiplier too. The width multiplier will multiply all widths of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void SetEffectWidthMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float width);

/*
Set the height of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current height of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param height			A float of the amount that you want to set the height multiplier too. The height multiplier will multiply all heights of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void SetEffectHeightMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float height);

/*
Set the depth of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current depth of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param depth			A float of the amount that you want to set the depth multiplier too. The depth multiplier will multiply all heights of emitters within the effect so it can be an easy way to alter the size
						of area, line, ellipse etc., emitters.
*/
tfxAPI void SetEffectDepthMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float depth);

/*
Set the life multiplier of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current life of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param life			A float of the amount that you want to set the life multiplier too. The life mulitplier will affect how long all particles emitted within the effect will last before expiring.
*/
tfxAPI void SetEffectLifeMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float life);

/*
Set the particle width multiplier of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current particle width of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param width			A float of the amount that you want to set the particle width multiplier too. The particle width mulitplier will affect the width of each particle if the emitter has a non uniform particle size, otherwise
						it will uniformly size the particle
*/
tfxAPI void SetEffectParticleWidthMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float width);

/*
Set the particle height multiplier of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current particle width of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param height			A float of the amount that you want to set the particle height multiplier too. The particle height mulitplier will affect the height of each particle if the emitter has a non uniform particle size, otherwise
						this function will have no effect.
*/
tfxAPI void SetEffectParticleHeightMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float height);

/*
Set the velocity multiplier of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current velocity of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param velocity		A float of the amount that you want to set the particle velocity multiplier too. The particle velocity mulitplier will affect the base velocity of a particle at spawn time.
*/
tfxAPI void SetEffectVelocityMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float velocity);

/*
Set the spin multiplier of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current spin of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param spin			A float of the amount that you want to set the particle spin multiplier too. The particle spin mulitplier will affect the base spin of a particle at spawn time.
*/
tfxAPI void SetEffectSpinMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float spin);

/*
Set the intensity multiplier of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current intensity of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param intensity		A float of the amount that you want to set the particle intensity multiplier too. The particle intensity mulitplier will instantly affect the opacity of all particles currently emitted by the effect.
*/
tfxAPI void SetEffectIntensityMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float intensity);

/*
Set the splatter multiplier of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current splatter of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param splatter		A float of the amount that you want to set the particle splatter multiplier too. The particle splatter mulitplier will change the amount of random offset all particles emitted in the effect will have.
*/
tfxAPI void SetEffectSplatterMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float splatter);

/*
Set the weight multiplier of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current weight of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param weight			A float of the amount that you want to set the particle weight multiplier too. The particle weight mulitplier will change the weight applied to particles in the effect at spawn time.
*/
tfxAPI void SetEffectWeightMultiplier(tfx_particle_manager_t *pm, tfxEffectID effect_index, float weight);

/*
Set the overal scale of an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed. Note that this must be called after UpdateParticleManager in order to override the current weight of the effect that was
*						set in the TimelineFX editor.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param overal_scale	A float of the amount that you want to set the overal scale to. The overal scale is an simply way to change the size of an effect
*/
tfxAPI void SetEffectOveralScale(tfx_particle_manager_t *pm, tfxEffectID effect_index, float overal_scale);

/*
Set the base noise offset for an effect
* @param pm				A pointer to a tfx_particle_manager_t where the effect is being managed.
* @param effect_index	The index of the effect. This is the index returned when calling AddEffectToParticleManager
* @param noise_offset	A float of the amount that you want to set the effect noise offset to. By default when an effect is added to a particle manager a random noise offset will be set based on the Base Noise Offset Range property. Here you can override that
						value by setting it here. The most ideal time to set this would be immediately after you have added the effect to the particle manager, but you could call it any time you wanted for a constantly changing noise offset.
*/
tfxAPI void SetEffectBaseNoiseOffset(tfx_particle_manager_t *pm, tfxEffectID effect_index, float noise_offset);

/*
Get the name of an effect
* @param pm				A pointer to the effect
* @returns				const char * name
*/
inline tfxAPI const char *GetEffectName(tfx_effect_emitter_t *effect) {
	return GetEffectInfo(effect)->name.c_str();
}

//-------Functions related to tfx_animation_manager_t--------

/*
Set the position of a 3d animation
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index			The index of the effect. This is the index returned when calling AddAnimationInstance
* @param position				A tfx_vec3_t vector object containing the x, y and z coordinates
*/
tfxAPI void SetAnimationPosition(tfx_animation_manager_t *animation_manager, tfxAnimationID animation_id, float position[3]);

/*
Set the position of a 2d animation
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index			The index of the effect. This is the index returned when calling AddAnimationInstance
* @param x						A float of the x position
* @param y						A float of the y position
*/
tfxAPI void SetAnimationPosition(tfx_animation_manager_t *animation_manager, tfxAnimationID animation_id, float x, float y);

/*
Set the scale of a 3d animation
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index			The index of the effect. This is the index returned when calling AddAnimationInstance
* @param scale					A multiplier that will determine the overal size/scale of the effect
*/
tfxAPI void SetAnimationScale(tfx_animation_manager_t *animation_manager, tfxAnimationID animation_id, float scale);

/*
Get an animation instance from an animation manager
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param tfxAnimationID			The index of the effect. This is the index returned when calling AddAnimationInstance
* @returns pointer to instance	Pointer to a tfx_animation_instance_t
*/
tfxAPI tfx_animation_instance_t *GetAnimationInstance(tfx_animation_manager_t *animation_manager, tfxAnimationID animation_id);

/*
Initialise an Animation Manager for use with 3d sprites. This must be run before using an animation manager. An animation manager is used
to playback pre recorded particle effects as opposed to using a particle manager that simulates the particles in
real time. This pre-recorded data can be uploaded to the gpu for a compute shader to do all the interpolation work
to calculate the state of particles between frames for smooth animation.
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param max_instances			The maximum number of animation instances that you want to be able to play at one time.
* @param initial_capacity		Optionally, you can set an initial capacity for the sprite data. The data will grow if you add
								beyond this amount but it gives you a chance to reserve a decent amount to start with to
								save too much mem copies as the data grows
*/
tfxAPI void InitialiseAnimationManagerFor3d(tfx_animation_manager_t *animation_manager, tfxU32 max_instances, tfxU32 initial_sprite_data_capacity = 100000);

/*
Initialise an Animation Manager for use with 2d sprites. This must be run before using an animation manager. An animation manager is used
to playback pre recorded particle effects as opposed to using a particle manager that simulates the particles in
real time. This pre-recorded data can be uploaded to the gpu for a compute shader to do all the interpolation work
to calculate the state of particles between frames for smooth animation.
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param max_instances			The maximum number of animation instances that you want to be able to play at one time.
* @param initial_capacity		Optionally, you can set an initial capacity for the sprite data. The data will grow if you add
								beyond this amount but it gives you a chance to reserve a decent amount to start with to
								save too much mem copies as the data grows
*/
tfxAPI void InitialiseAnimationManagerFor2d(tfx_animation_manager_t *animation_manager, tfxU32 max_instances, tfxU32 initial_sprite_data_capacity = 100000);

/*
Set the callback that you can use to determine whether or not a tfx_animation_instance_t should be added to the next frame's render queue. You can use this
to cull instances that are outside of the view frustum for example
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param callback				Pointer to the callback you want to use. It must have the following signature:
								bool(*maybe_render_instance_callback(tfx_animation_manager_t *animation_manager, tfx_animation_instance_t *instance, tfx_frame_meta_t *meta, void *user_data))
								Values passed into the callback function are a pointer to the animation manager, a pointer to the instance being processed, a pointer to
								the frame meta of the instance, this will contain the bounding box and radius of the instance from the current frame of the instance and a pointer
								to any user data that you set that might contain the camera frustum that you want to check against.
*/
tfxAPI void SetAnimationManagerInstanceCallback(tfx_animation_manager_t *animation_manager, bool((*maybe_render_instance_callback)(tfx_animation_manager_t *animation_manager, tfx_animation_instance_t *instance, tfx_frame_meta_t *meta, void *user_data)));

/*
Set the user data in a tfx_animation_manager_t which can get passed through to callback functions when updated the animation manager
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param user_data				void* pointer to the data that you want to set
*/
tfxAPI void SetAnimationManagerUserData(tfx_animation_manager_t *animation_manager, void *user_data);

/*
Add sprite data to an animation manager sprite data buffer from an effect. This will record the
animation if necessary and then convert the sprite data to tfx_sprite_data3d_t ready for uploading
to the GPU
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param effect_index			The index of the effect. This is the index returned when calling AddAnimationInstance
* @param position				A tfx_vec3_t vector object containing the x, y and z coordinates
*/
tfxAPI void AddSpriteData(tfx_animation_manager_t *animation_manager, tfx_effect_emitter_t *effect, tfx_particle_manager_t *pm = nullptr, tfx_vec3_t camera_position = { 0.f, 0.f, 0.f });

/*
Add an animation instance to the animation manager.
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param path					tfxKey path hash of the effect name and path: effect.path_hash
* @param start_frame			Starting frame of the animation
* @returns						The index id of the animation instance. You can use this to reference the animation when changing position, scale etc
								Return tfxINVALID if there is no room in the animation manager
*/
tfxAPI tfxAnimationID AddAnimationInstance(tfx_animation_manager_t *animation_manager, tfxKey path, tfxU32 start_frame = 0);

/*
Add an animation instance to the animation manager.
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @param path					const char * name of the effect. If the effect was in a folder then specify the whole path
* @param start_frame			Starting frame of the animation
* @returns						The index id of the animation instance. You can use this to reference the animation when changing position, scale etc
								Return tfxINVALID if there is no room in the animation manager
*/
tfxAPI tfxAnimationID AddAnimationInstance(tfx_animation_manager_t *animation_manager, const char *path, tfxU32 start_frame = 0);

/*
Update an animation manager to advance the time and frames of all instances currently playing.
* @param animation_manager		A pointer to a tfx_animation_manager_t that you want to update
* @param start_frame			Starting frame of the animation
*/
tfxAPI void UpdateAnimationManager(tfx_animation_manager_t *animation_manager, float elapsed);

/*
Add an effect's shapes to an animation manager. You can use this function if you're manually recording particle effects and adding them to an animation
manager rather then just using the editor.
* @param animation_manager		A pointer to a tfx_animation_manager_t that you want to update
* @param effect					A pointer to the effect whose shapes you want to add
*/
tfxAPI void AddEffectShapes(tfx_animation_manager_t *animation_manager, tfx_effect_emitter_t *effect);

/*
Update an animation manager so that the effects do not expire they just loop forever instead regardless of whether they're a looped effect or not.
* @param animation_manager		A pointer to a tfx_animation_manager_t that you want to update
*/
tfxAPI void CycleAnimationManager(tfx_animation_manager_t *animation_manager);

/*
Clears all animation instances currently in play in an animation manager, resulting in all currently running animations
from being drawn
* @param animation_manager		A pointer to a tfx_animation_manager_t that you want to clear
*/
tfxAPI void ClearAllAnimationInstances(tfx_animation_manager_t *animation_manager);

/*
Clears all data from the animation manager including sprite data, metrics and instances. Essentially resetting everything back to
it's initialisation point
from being drawn
* @param animation_manager		A pointer to a tfx_animation_manager_t that you want to reset
*/
tfxAPI void ResetAnimationManager(tfx_animation_manager_t *animation_manager);

/*
Frees all data from the animation manager including sprite data, metrics and instances.
from being drawn
* @param animation_manager		A pointer to a tfx_animation_manager_t that you want to reset
*/
tfxAPI void FreeAnimationManager(tfx_animation_manager_t *animation_manager);

/*
Get the tfx_animation_buffer_metrics_t from an animation manager. This will contain the info you need to upload the sprite data,
offsets and animation instances to the GPU. Only offsets and animation instances need to be uploaded to the GPU each frame. Sprite
data can be done ahead of time.
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @returns						tfx_animation_buffer_metrics_t containing buffer sizes
*/
tfxAPI inline tfx_animation_buffer_metrics_t GetAnimationBufferMetrics(tfx_animation_manager_t *animation_manager) {
	return animation_manager->buffer_metrics;
}

/*
Get the total number of sprites that need to be drawn by an animation manager this frame. You can use this in your renderer
to draw your sprite instances
* @param animation_manager		A pointer to a tfx_animation_manager_t where the effect animation is being managed
* @returns						tfxU32 of the number of sprites
*/
tfxAPI inline tfxU32 GetTotalSpritesThatNeedDrawing(tfx_animation_manager_t *animation_manager) {
	return animation_manager->buffer_metrics.total_sprites_to_draw;
}

/*
Get the total number of instances being processed by an animation manager. This will not necessarily be the same number as
the instances being rendered if some are being culled in your custom callback if your using one.
* @param animation_manager		A pointer to a tfx_animation_manager_t that you want to clear
* @returns int					The number of instances being updated
*/
tfxAPI inline tfxU32 GetTotalInstancesBeingUpdated(tfx_animation_manager_t *animation_manager) {
	return animation_manager->instances_in_use[animation_manager->current_in_use_buffer].size();
}

/*
Create the image data required for GPU shaders such as animation viewer. The image data will contain data such as uv coordinates
that the shaders can use to create the sprite data. Once you have built the data you can use GetLibraryImageData to get the buffer
and upload it to the gpu.
* @param library				A pointer to some image data where the image data is. You can use GetParticleShapes with a tfx_library_t or tfx_animation_manager_t for this
* @param uv_lookup				A function pointer to a function that you need to set up in order to get the uv coordinates from whatever renderer you're using
*/
tfxAPI tfx_gpu_shapes_t BuildGPUShapeData(tfx_vector_t<tfx_image_data_t> *particle_shapes, tfx_vec4_t(uv_lookup)(void *ptr, tfx_gpu_image_data_t *image_data, int offset));

/*
Get a pointer to the GPU shapes which you can use in a memcpy
* @param particle_shapes		A pointer the tfx_gpu_shapes_t
*/
tfxAPI inline void *GetGPUShapesPointer(tfx_gpu_shapes_t *particle_shapes) {
	return particle_shapes->list.data;
}

/*
Get a pointer to the particle shapes data in the animation manager. This can be used with BuildGPUShapeData when you want to upload the data to the GPU
* @param animation_manager		A pointer the tfx_animation_manager_t
*/
tfxAPI inline tfx_vector_t<tfx_image_data_t> *GetParticleShapes(tfx_animation_manager_t *animation_manager) {
	return &animation_manager->particle_shapes.data;
}

/*
Get a pointer to the particle shapes data in the animation manager. This can be used with BuildGPUShapeData when you want to upload the data to the GPU
* @param animation_manager		A pointer the tfx_animation_manager_t
*/
tfxAPI inline tfx_vector_t<tfx_image_data_t> *GetParticleShapes(tfx_library_t *library) {
	return &library->particle_shapes.data;
}

/*
Get the number of shapes in the GPU Shape Data buffer. Make sure you call BuildGPUShapeData first or they'll be nothing to return
* @param library				A pointer to a tfx_animation_manager_t where the image data will be created.
* @returns tfxU32				The number of shapes in the buffer
*/
tfxAPI inline tfxU32 GetGPUShapeCount(tfx_gpu_shapes_t *particle_shapes) {
	return particle_shapes->list.size();
}

/*
Get the size in bytes of the GPU image data in a tfx_library_t
* @param library				A pointer to a tfx_library_t where the image data exists.
* @returns size_t				The size in bytes of the image data
*/
tfxAPI inline size_t GetGPUShapesSizeInBytes(tfx_gpu_shapes_t *particle_shapes) {
	return particle_shapes->list.size_in_bytes();
}

/*
Get the total number of sprites in an animation manger's sprite data buffer
* @param animation_manager		A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32				The number of sprites in the buffer
*/
tfxAPI inline tfxU32 GetTotalSpriteDataCount(tfx_animation_manager_t *animation_manager) {
	if (animation_manager->flags & tfxAnimationManagerFlags_is_3d) {
		return animation_manager->sprite_data_3d.current_size;
	}
	return animation_manager->sprite_data_2d.current_size;
}

/*
Get the total number of sprites in an animation manger's sprite data buffer
* @param animation_manager		A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32				The number of sprites in the buffer
*/
tfxAPI inline size_t GetSpriteDataSizeInBytes(tfx_animation_manager_t *animation_manager) {
	if (animation_manager->flags & tfxAnimationManagerFlags_is_3d) {
		return animation_manager->sprite_data_3d.size_in_bytes();
	}
	return animation_manager->sprite_data_2d.size_in_bytes();
}

/*
Get the buffer memory address for the sprite data in an animation manager
* @param animation_manager		A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*				A pointer to the sprite data memory
*/
tfxAPI inline void* GetSpriteDataBufferPointer(tfx_animation_manager_t *animation_manager) {
	if (animation_manager->flags & tfxAnimationManagerFlags_is_3d) {
		return animation_manager->sprite_data_3d.data;
	}
	return animation_manager->sprite_data_2d.data;
}

/*
Get the size in bytes of the offsets buffer in an animation manager
* @param animation_manager		A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t				Size in bytes of the offsets buffer
*/
tfxAPI inline size_t GetOffsetsSizeInBytes(tfx_animation_manager_t *animation_manager) {
	return animation_manager->offsets.current_size * sizeof(tfxU32);
}

/*
Get the size in bytes of the render queue of animation instances buffer in an animation manager
* @param animation_manager		A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t				Size in bytes of the instances buffer
*/
tfxAPI inline size_t GetAnimationInstancesSizeInBytes(tfx_animation_manager_t *animation_manager) {
	return animation_manager->render_queue.current_size * sizeof(tfx_animation_instance_t);
}

/*
Get the size in bytes of the animation emitter properties list
* @param animation_manager		A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns size_t				Size in bytes of the properties bufffer
*/
tfxAPI inline size_t GetAnimationEmitterPropertySizeInBytes(tfx_animation_manager_t *animation_manager) {
	return animation_manager->emitter_properties.current_size * sizeof(tfx_animation_emitter_properties_t);
}

/*
Get the number of emitter properties being using by the animation manager
* @param animation_manager		A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns tfxU32				Number of emitter properties
*/
tfxAPI inline tfxU32 GetAnimationEmitterPropertyCount(tfx_animation_manager_t *animation_manager) {
	return animation_manager->emitter_properties.current_size;
}

/*
Get the buffer memory address for the sprite data in an animation manager
* @param animation_manager		A pointer to a tfx_animation_manager_t to get the sprite data from
* @returns void*				A pointer to the sprite data memory
*/
tfxAPI void* GetAnimationEmitterPropertiesBufferPointer(tfx_animation_manager_t *animation_manager);

/*
Reset an effect template and make it empty so you can use it to store another effect.
* @param t						A pointer to a tfx_effect_template_t
*/
tfxAPI void ResetTemplate(tfx_effect_template_t *t);

/*
Get the root effect from the template
* @param t						A pointer to a tfx_effect_template_t
* @returns						A pointer to the root effect
*/
tfxAPI tfx_effect_emitter_t *GetEffectFromTemplate(tfx_effect_template_t *t);

/*
Get an emitter or sub effect from an effect template.
* @param t						A pointer to a tfx_effect_template_t
* @param path					A path to the emitter or sub effect that you want to retrieve. Must be a valid path. Example path might be: "Explosion/Smoke"
* @returns						A pointer to the root effect
*/
tfxAPI tfx_effect_emitter_t *GetEmitterFromTemplate(tfx_effect_template_t *t, tfx_str256_t *path);

/*
Set the user data for any effect or emitter in the effect template. This user data will get passed through to any update callback functions
* @param t						A pointer to a tfx_effect_template_t
* @param path					A path to the effect or emitter in the effect template
* @param data					A pointer to the user data
*/
tfxAPI void SetTemplateUserData(tfx_effect_template_t *t, tfx_str256_t *path, void *data);

/*
Set the user data for the root effect in an effect template
* @param t						A pointer to a tfx_effect_template_t
* @param data					A pointer to the user data
*/
tfxAPI void SetTemplateEffectUserData(tfx_effect_template_t *t, void *data);

/*
Set the same user data for all effects and emitters/sub effects in the effect template
* @param t						A pointer to a tfx_effect_template_t
* @param data					A pointer to the user data that will be set to all effects and emitters in the template
*/
tfxAPI void SetTemplateUserDataAll(tfx_effect_template_t *t, void *data);

/*
Set an update callback for the root effect in the effect template.
* @param t						A pointer to a tfx_effect_template_t
* @param update_callback		A pointer to the call back function
*/
tfxAPI void SetTemplateEffectUpdateCallback(tfx_effect_template_t *t, void(*update_callback)(tfx_particle_manager_t *pm, tfxEffectID effect_index));

/*
Pre-record this effect into a sprite cache so that you can play the effect back without the need to actually caclulate particles in realtime.
	* @param pm			Reference to a pm that will be used to run the particle simulation and record the sprite data
	* @param path		const *char of a path to the emitter in the effect.Must be a valid path, for example: "My Effect/My Emitter"
	* @param camera		Array of 3 floats with the camera position (only needed for 3d effects that are sorted by depth
*/
tfxAPI void RecordTemplateEffect(tfx_effect_template_t *t, tfx_particle_manager_t *pm, float update_frequency, float camera_position[3]);

/*
Disable an emitter within an effect. Disabling an emitter will stop it being added to the particle manager when calling AddEffectToParticleManager
* @param path		const *char of a path to the emitter in the effect. Must be a valid path, for example: "My Effect/My Emitter"
*/
tfxAPI void DisableTemplateEmitter(tfx_effect_template_t *t, const char *path);

/*
Enable an emitter within an effect so that it is added to the particle manager when calling AddEffectToParticleManager. Emitters are enabled by default.
* @param path		const *char of a path to the emitter in the effect. Must be a valid path, for example: "My Effect/My Emitter"
*/
tfxAPI void EnableTemplateEmitter(tfx_effect_template_t *t, const char *path);

/*
Scale all nodes on a global graph graph of the effect
* @param global_type		tfx_graph_type of the global graph that you want to scale. Must be a global graph or an assert will be called
* @param amount				A float of the amount that you want to scale the multiplier by.
*/
tfxAPI void ScaleTemplateGlobalMultiplier(tfx_effect_template_t *t, tfx_graph_type global_type, float amount);

/*
Scale all nodes on an emitter graph
* @param emitter_path		const *char of the emitter path
* @param global_type		tfx_graph_type of the emitter graph that you want to scale. Must be an emitter graph or an assert will be called
* @param amount				A float of the amount that you want to scale the graph by.
*/
tfxAPI void ScaleTemplateEmitterGraph(tfx_effect_template_t *t, const char *emitter_path, tfx_graph_type graph_type, float amount);

/*
Set the single spawn amount for an emitter. Only affects emitters that have the single spawn flag set.
* @param emitter_path		const *char of the emitter path
* @param amount				A float of the amount that you want to set the single spawn amount to.
*/
tfxAPI void SetTemplateSingleSpawnAmount(tfx_effect_template_t *t, const char *emitter_path, tfxU32 amount);

/*
Interpolate between 2 tfxVec3s. You can make use of this in your render function when rendering sprites and interpolating between captured and current positions
* @param tween				The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param world				The current tvxVec3 position
* @param captured			The captured tvxVec3 position
* @returns tfx_vec3_t			The interpolated tfx_vec3_t
*/
tfxAPI inline tfx_vec3_t Tween3d(float tween, const tfx_vec3_t *world, const tfx_vec3_t *captured) {
	tfx_vec3_t tweened;
	tweened = *world * tween + *captured * (1.f - tween);
	return tweened;
}

/*
Interpolate between 2 tfxVec2s. You can make use of this in your render function when rendering sprites and interpolating between captured and current positions
* @param tween		The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param world		The current tvxVec2 position
* @param captured	The captured tvxVec2 position
* @returns tfx_vec2_t	The interpolated tfx_vec2_t
*/
tfxAPI inline tfx_vec2_t Tween2d(float tween, const tfx_vec2_t *world, const tfx_vec2_t *captured) {
	tfx_vec2_t tweened;
	tweened = *world * tween + *captured * (1.f - tween);
	return tweened;
}

/*
Interpolate between 2 float. You can make use of this in your render function when rendering sprites and interpolating between captured and current float values like intensity
* @param tween		The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param world		The current tvxVec2 position
* @param captured	The captured tvxVec2 position
* @returns tfx_vec2_t	The interpolated tfx_vec2_t
*/
tfxAPI inline float TweenFloat(float tween, const float current, const float captured) {
	return current * tween + captured * (1.f - tween);
}

#ifdef tfxINTEL
/*
Interpolate between 2 colors in tfx_rgba8_t format. You can make use of this in your render function when rendering sprites and interpolating between captured and current colors
* @param tween                The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param current            The current tfx_rgba8_t color
* @param captured            The captured tfx_rgba8_t color
* @returns tfx_rgba8_t            The interpolated tfx_rgba8_t
*/
tfxAPI inline tfx_rgba8_t TweenColor(float tween, const tfx_rgba8_t current, const tfx_rgba8_t captured) {
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
    return tfx_rgba8_t(packed.a[0], packed.a[1], packed.a[2], packed.a[3]);
}

/*
Interpolate all sprite transform data in a single function. This will interpolate position, scale and rotation.
* @param tween                The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param current            The current transform struct of the sprite
* @param captured            The captured transform struct of the sprite
* @returns tfx_wide_lerp_transform_result_t            The interpolated transform data in a tfx_wide_lerp_transform_result_t
*/
tfxAPI inline tfx_wide_lerp_transform_result_t InterpolateSpriteTransform(const __m128 *tween, const tfx_sprite_transform3d_t *current, const tfx_sprite_transform3d_t *captured) {
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

#elifdef tfxARM

/*
Interpolate between 2 colors in tfx_rgba8_t format. You can make use of this in your render function when rendering sprites and interpolating between captured and current colors
* @param tween                The interpolation value between 0 and 1. You should pass in the value from your timing function
* @param current            The current tfx_rgba8_t color
* @param captured            The captured tfx_rgba8_t color
* @returns tfx_rgba8_t            The interpolated tfx_rgba8_t
*/
tfxAPI inline tfx_rgba8_t TweenColor(float tween, const tfx_rgba8_t current, const tfx_rgba8_t captured) {
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
tfxAPI inline tfx_wide_lerp_transform_result_t InterpolateSpriteTransform(const tfx128 *tween, const tfx_sprite_transform3d_t *current, const tfx_sprite_transform3d_t *captured) {
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

} //namespace
#endif
