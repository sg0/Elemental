/*
   Copyright (c) 2009-2015, Jack Poulson
                      2013, Jeff Hammond
                      2013, Jed Brown
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
http://opensource.org/licenses/BSD-2-Clause
*/

#include "El.hpp"
#include <cassert>

typedef unsigned char* UCP;

namespace
{
inline void SafeMpi (int mpiError)
{
    DEBUG_ONLY(
        if( mpiError != MPI_SUCCESS )    
        {
            char errorString[MPI_MAX_ERROR_STRING];
            int lengthOfErrorString;
            MPI_Error_string( mpiError, errorString, &lengthOfErrorString );
            El::RuntimeError( std::string(errorString) );
        }
    )
}
}			// anonymous namespace

namespace El
{
namespace mpi
{

bool CommSameSizeAsInteger()
{ return sizeof(MPI_Comm) == sizeof(int); }

bool GroupSameSizeAsInteger()
{ return sizeof(MPI_Group) == sizeof(int); }

// MPI environmental routines
// ==========================

void Initialize (int &argc, char **&argv)
{
    MPI_Init (&argc, &argv);
}

int InitializeThread (int &argc, char **&argv,
                      int required)
{
    int provided;
#ifdef EL_HAVE_MPI_INIT_THREAD
    MPI_Init_thread (&argc, &argv, required, &provided);
#else
    MPI_Init (&argc, &argv);
    provided = 0;	// equivalent to MPI_THREAD_SINGLE
#endif
    return provided;
}

void Finalize ()
{
    MPI_Finalize ();
}

bool Initialized ()
{
    int initialized;

    MPI_Initialized (&initialized);
    return initialized;
}

bool Finalized ()
{
    int finalized;

    MPI_Finalized (&finalized);
    return finalized;
}

int QueryThread ()
{
    int provided;

#ifdef EL_HAVE_MPI_QUERY_THREAD
    MPI_Query_thread (&provided);
#else
    provided = 0;	// equivalent to MPI_THREAD_SINGLE
#endif
    return provided;
}

void Abort (Comm comm, int errCode)
{
    MPI_Abort (comm.comm, errCode);
}

double Time ()
{
    return MPI_Wtime ();
}

void Create (UserFunction * func, bool commutes, Op & op)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Create"))
    SafeMpi (MPI_Op_create (func, commutes, &op.op));
}

void Free (Op & op)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Free"))
    SafeMpi (MPI_Op_free (&op.op));
}

// Communicator manipulation
// =========================

int WorldRank ()
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::WorldRank"))
    return Rank (mpi::COMM_WORLD);
}

int Rank (Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Rank"))
    if (comm != COMM_NULL)
    {
        int rank;

        SafeMpi (MPI_Comm_rank (comm.comm, &rank));
        return rank;
    }
    else
        return mpi::UNDEFINED;
}

int Size (Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Size"))
    if (comm != COMM_NULL)
    {
        int size;

        SafeMpi (MPI_Comm_size (comm.comm, &size));
        return size;
    }
    else
        return mpi::UNDEFINED;
}

void Create (Comm parentComm, Group subsetGroup,
             Comm & subsetComm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Create"))
    SafeMpi (MPI_Comm_create
             (parentComm.comm, subsetGroup.group,
              &subsetComm.comm));
}

void Dup (Comm original, Comm & duplicate)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Dup"))
    SafeMpi (MPI_Comm_dup
             (original.comm, &duplicate.comm));
}

void Split (Comm comm, int color, int key, Comm & newComm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Split"))
    SafeMpi (MPI_Comm_split
             (comm.comm, color, key, &newComm.comm));
}

void Free (Comm & comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Free"))
    SafeMpi (MPI_Comm_free (&comm.comm));
}

bool Congruent (Comm comm1, Comm comm2)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Congruent")) int
    result;
    SafeMpi (MPI_Comm_compare
             (comm1.comm, comm2.comm, &result));
    return (result == MPI_IDENT
            || result == MPI_CONGRUENT);
}

void ErrorHandlerSet (Comm comm,
                      ErrorHandler errorHandler)
{
    DEBUG_ONLY (CallStackEntry
                cse ("mpi::ErrorHandlerSet"))
#ifdef EL_HAVE_MPI_COMM_SET_ERRHANDLER
    SafeMpi (MPI_Comm_set_errhandler
             (comm.comm, errorHandler));
#else
    SafeMpi (MPI_Errhandler_set
             (comm.comm, errorHandler));
#endif
}

// Cartesian communicator routines
// ===============================

void CartCreate
(Comm comm, int numDims, const int *dimensions,
 const int *periods, bool reorder, Comm & cartComm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::CartCreate"))
    SafeMpi
    (MPI_Cart_create
     (comm.comm, numDims,
      const_cast < int *>(dimensions),
      const_cast < int *>(periods), reorder,
      &cartComm.comm));
}

void CartSub (Comm comm, const int *remainingDims,
              Comm & subComm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::CartSub"))
    SafeMpi (MPI_Cart_sub
             (comm.comm,
              const_cast < int *>(remainingDims),
              &subComm.comm));
}

// Group manipulation
// ==================

int Rank (Group group)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Rank")) int
    rank;
    SafeMpi (MPI_Group_rank (group.group, &rank));
    return rank;
}

int Size (Group group)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Size")) int
    size;
    SafeMpi (MPI_Group_size (group.group, &size));
    return size;
}

void CommGroup (Comm comm, Group & group)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::CommGroup"))
    SafeMpi (MPI_Comm_group
             (comm.comm, &group.group));
}

void Dup (Group group, Group & newGroup)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Dup"))
    // For some reason, MPI_Group_dup does not exist
    Excl (group, 0, 0, newGroup);
}

void Union (Group groupA, Group groupB, Group & newGroup)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Union"))
    SafeMpi (MPI_Group_union
             (groupA.group, groupB.group,
              &newGroup.group));
}

void Incl (Group group, int n, const int *ranks,
           Group & subGroup)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Incl"))
    SafeMpi (MPI_Group_incl
             (group.group, n,
              const_cast < int *>(ranks),
              &subGroup.group));
}

void Excl (Group group, int n, const int *ranks,
           Group & subGroup)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Excl"))
    SafeMpi (MPI_Group_excl
             (group.group, n,
              const_cast < int *>(ranks),
              &subGroup.group));
}

void Difference (Group parent, Group subset,
                 Group & complement)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Difference"))
    SafeMpi (MPI_Group_difference
             (parent.group, subset.group,
              &complement.group));
}

void Free (Group & group)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Free"))
    SafeMpi (MPI_Group_free (&group.group));
}

// Rank translations
// =================

int Translate (Group origGroup, int origRank,
               Group newGroup)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Translate")) int
    newRank;
    Translate (origGroup, 1, &origRank, newGroup,
               &newRank);
    return newRank;
}

int Translate (Comm origComm, int origRank,
               Group newGroup)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Translate")) int
    newRank;
    Translate (origComm, 1, &origRank, newGroup,
               &newRank);
    return newRank;
}

int Translate (Group origGroup, int origRank,
               Comm newComm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Translate")) int
    newRank;
    Translate (origGroup, 1, &origRank, newComm,
               &newRank);
    return newRank;
}

int Translate (Comm origComm, int origRank, Comm newComm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Translate")) int
    newRank;
    Translate (origComm, 1, &origRank, newComm, &newRank);
    return newRank;
}

void Translate
(Group origGroup, int size, const int *origRanks,
 Group newGroup, int *newRanks)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Translate"))
    SafeMpi
    (MPI_Group_translate_ranks
     (origGroup.group, size,
      const_cast < int *>(origRanks),
      newGroup.group, newRanks));
}

void Translate
(Comm origComm, int size, const int *origRanks,
 Group newGroup, int *newRanks)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Translate"))
    Group origGroup;

    CommGroup (origComm, origGroup);
    Translate (origGroup, size, origRanks, newGroup,
               newRanks);
    Free (origGroup);
}

void Translate
(Group origGroup, int size, const int *origRanks,
 Comm newComm, int *newRanks)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Translate"))
    Group newGroup;

    CommGroup (newComm, newGroup);
    Translate (origGroup, size, origRanks, newGroup,
               newRanks);
    Free (newGroup);
}

void Translate
(Comm origComm, int size, const int *origRanks,
 Comm newComm, int *newRanks)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Translate"))
    Group origGroup, newGroup;

    CommGroup (origComm, origGroup);
    CommGroup (newComm, newGroup);
    Translate (origGroup, size, origRanks, newGroup,
               newRanks);
    Free (origGroup);
    Free (newGroup);
}

// DERIVED Datatype creation
// =========================
// FIXME these functions for DDT creation are 
// completely untested
#ifdef EL_USE_DERIVED_DATATYPE
void StridedDatatype (El_strided_t* stride_descr,
	Datatype old_type, Datatype* new_type,
	size_t* source_dims)
{
    int old_type_size;
    SafeMpi (MPI_Type_size (old_type, &old_type_size));
    int *dims = NULL, *sizes = NULL;

    // count of blocks must be non-zero
    assert (stride_descr->num > 0);
    // size is NULL
    assert (stride_descr->sizes != NULL);
    // offset is NULL
    assert (stride_descr->offsets != NULL);

    // check for contiguous transfers
    if ((source_dims == NULL) && (stride_descr->num == 1))
    {
	int elem_count = stride_descr->sizes[0] / old_type_size;
	// derived datatype is not a multiple of original type
	assert ((stride_descr->sizes[0] % old_type_size == 0));
	SafeMpi ( MPI_Type_contiguous (elem_count, old_type, new_type) );
	return;
    }
    // offsets should be monotonic increasing
    for (int i = 1; i < stride_descr->num; i++)
	assert (stride_descr->offsets[i] >= stride_descr->offsets[i - 1]);
    /* Notes: 
     * Sayan: This weird hack is because MPI_Type_create_subarray throws an error when
     * stride_descr->sizes and source_dims is passed directly (probably type mismatch?) */
    /* heap */
    dims = new int[stride_descr->num];
    sizes = new int[stride_descr->num];

    for (int i = 0; i < stride_descr->num; i++)
    {
	dims[i] = EL_INT_SAFE_CAST (source_dims[i]);
	sizes[i] = EL_INT_SAFE_CAST (stride_descr->sizes[i]);
    }

    SafeMpi ( MPI_Type_create_subarray (stride_descr->num, reinterpret_cast<const int *>(dims),
		reinterpret_cast<const int *>(sizes),
		reinterpret_cast<int *>(stride_descr->offsets), MPI_ORDER_C,
		old_type, new_type) );

    delete[] dims;
    delete[] sizes;
}

void VectorDatatype (El_iov_t * vect_descr,
	Datatype old_type, Datatype * new_type,
	vector_pattern_t data_pattern)
{
    int old_type_size;
    int stride;
    int fixed_block_fixed_stride = 1,  // MPI_Type_vector
	fixed_block_var_stride = 1;	     // MPI_Type_hindexed_block
    /* defaults:
     * var_block_var_stride=1 - MPI_Type_hindexed
     * var_block_fixed_stride=1 -  MPI_Type_hindexed 
     */
    SafeMpi ( MPI_Type_size (old_type, &old_type_size) );
    // count of blocks must be non-zero
    assert (vect_descr->count > 0);
    // size is NULL
    assert (vect_descr->sizes != NULL);
    // offset is NULL
    assert (vect_descr->offsets != NULL);
    // check for contiguous transfers
    if (vect_descr->count == 1)
    {
	int elem_count = vect_descr->sizes[0] / old_type_size;
	// derived datatype is not a multiple of original type
	assert (vect_descr->sizes[0] % old_type_size == 0);
	SafeMpi ( MPI_Type_contiguous (elem_count, old_type, new_type) );
	return;
    }
    // offsets should be monotonic increasing
    for (int i = 1; i < vect_descr->count; i++)
	assert (vect_descr->offsets[i] >= vect_descr->offsets[i - 1]);

    // identify the pattern of strides, fixed or varying
    if (data_pattern == UNKNOWN_BLOCK_STRIDE)
    {
	stride = (vect_descr->offsets[1] - vect_descr->offsets[0]);
	for (int i = 1; i < vect_descr->count; i++)
	{
	    // check for fixed blocklengths and fixed strides
	    if ((vect_descr->sizes[i] == vect_descr->sizes[i - 1]) &&
		    (stride ==
		     (vect_descr->offsets[i] - vect_descr->offsets[i - 1])))
		fixed_block_fixed_stride++;

	    // check for fixed blocklengths and variable strides
	    if ((vect_descr->sizes[i] == vect_descr->sizes[i - 1]) &&
		    !(stride ==
			(vect_descr->offsets[i] - vect_descr->offsets[i - 1])))
		fixed_block_var_stride++;
	}
    }

    if (data_pattern == FIXED_BLOCK_FIXED_STRIDE)
	fixed_block_fixed_stride = vect_descr->count;

    if (data_pattern == FIXED_BLOCK_VAR_STRIDE)
	fixed_block_var_stride = vect_descr->count;

    // check if constant strides, if yes 
    // then create _type_vector, else 
    // _type_hindexed 
    if (fixed_block_fixed_stride == vect_descr->count)
    {				// _vector
	int stride = ((vect_descr->offsets[1] - vect_descr->offsets[0])
		/ old_type_size);
	int blocklength = vect_descr->sizes[0];
	SafeMpi ( MPI_Type_vector (vect_descr->count, blocklength,
		    stride, old_type, new_type) );
    }
    else if (fixed_block_var_stride == vect_descr->count)	// _hindexed_block
	SafeMpi ( MPI_Type_create_hindexed_block (vect_descr->count, vect_descr->sizes[0],
		    vect_descr->offsets, old_type, new_type) );
    else				// _hindexed
	SafeMpi ( MPI_Type_create_hindexed (vect_descr->count,
		    (const int *) vect_descr->sizes,
		    vect_descr->offsets, old_type, new_type) );
}
#endif // EL_USE_DERIVED_DATATYPE

// MPI-3 RMA functions
// ==================

#if MPI_VERSION>=3 && defined(EL_ENABLE_RMA_AXPY)
long ReadInc (Window & win, Aint offset, long inc, int fop_root)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::ReadInc"))
    long otemp;			
    SafeMpi ( MPI_Fetch_and_op (&inc, &otemp, MPI_LONG, fop_root, offset, MPI_SUM,
	    win) );
    SafeMpi ( MPI_Win_flush_local (fop_root, win) );

    return otemp;
}

void SetWindowProp (Window & window, int prop)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::SetWindowProp"))
    Info info;

    SafeMpi (MPI_Info_create (&info));

    if (prop & (1 << 0))	// strict
        SafeMpi (MPI_Info_set
                 (info, "accumulate_ordering",
                  "rar,raw,war,waw"));


    if (prop & (1 << 1))	// partial
        SafeMpi (MPI_Info_set
                 (info, "accumulate_ordering",
                  "rar,waw"));

    if (prop & (1 << 2))	// none
        SafeMpi (MPI_Info_set
                 (info, "accumulate_ops",
                  "same_op_no_op"));

    SafeMpi (MPI_Win_set_info (window, info));
}

//NOTE assuming MPI_MODE_NOCHECK
void WindowLock (int rank, Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::WindowLock"))
    SafeMpi (MPI_Win_lock
             (MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK,
              window));
}

void WindowLock (Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::WindowLock"))
    SafeMpi (MPI_Win_lock_all
             (MPI_MODE_NOCHECK, window));
}

void WindowUnlock (int rank, Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::WindowUnlock"))
    SafeMpi (MPI_Win_unlock (rank, window));
}

void WindowUnlock (Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::WindowUnlock"))
    SafeMpi (MPI_Win_unlock_all (window));
}

// RMA Utilities
void WindowCreate (void *baseptr, int size, Comm comm, Window & window)
{
    DEBUG_ONLY( CallStackEntry cse ("mpi::WindowCreate") )

    SafeMpi( MPI_Win_create
             ( baseptr, (MPI_Aint) size, 1, MPI_INFO_NULL,
              comm.comm, &window ) );
#ifdef EL_NO_ACC_ORDERING
    SetWindowProp( window, NO_ACC_ORDERING );
#endif
}


void CheckBounds (Window & window, Datatype win_type, Datatype type,
                  size_t count, ptrdiff_t target_offset)
{
    int flag, type_size, win_type_size;
    size_t displ;
    void * dest=NULL;

    SafeMpi (MPI_Type_size (type, &type_size));
    SafeMpi (MPI_Type_size (win_type, &win_type_size));
    Aint lb, extent;

    SafeMpi (MPI_Win_get_attr(window, MPI_WIN_BASE, dest, &flag /* unused */));

    /* Calculate displacement from beginning of the window */
    if (dest == MPI_BOTTOM)
        displ = 0;
    else
        displ = (size_t) ((uint8_t*)((uint8_t*)dest + target_offset * type_size) - (uint8_t*)dest);

    SafeMpi (MPI_Type_get_true_extent(type, &lb, &extent));

    // invalid remote address
    assert (displ >= 0 && displ < win_type_size);
    // transfer out of range
    assert (displ + count*extent <= win_type_size);
}

#ifdef EL_EXPLICIT_PROGRESS
void Progress ( Comm comm )
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Progress"))
    int flag;
    SafeMpi (MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, 
		comm.comm, &flag, MPI_STATUS_IGNORE));
}
#endif

void WindowFree (Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::WindowFree"))
    SafeMpi (MPI_Win_free (&window));
}

void * GetWindowBase (Window & window)
{
    DEBUG_ONLY( CallStackEntry cse ("mpi::GetWindowBase") )

    int flag = 0;
    // TODO check flag
    void * attribute_val = NULL;
    SafeMpi( MPI_Win_get_attr
             ( window, MPI_WIN_BASE, &attribute_val, &flag) );
    return attribute_val;
}

void WindowAllocate (int size, Comm comm, Window & window)
{
    DEBUG_ONLY( CallStackEntry cse ("mpi::WindowAllocate") )
    
    void * base = NULL;
    SafeMpi( MPI_Win_allocate
             ( (MPI_Aint) size, 1, MPI_INFO_NULL,
              comm.comm, &base, &window) );

#ifdef EL_NO_ACC_ORDERING
    SetWindowProp( window, NO_ACC_ORDERING );
#endif
}

// put
template<typename R>
void Iput (const R* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Iput"))
#ifdef EL_ENSURE_PUT_ATOMICITY
    SafeMpi (MPI_Accumulate
             (source, origin_count, TypeMap<R>(),
              target_rank, disp, target_count,
              TypeMap<R>(), MPI_REPLACE, window));
#else
    SafeMpi (MPI_Put
             (source, origin_count, TypeMap<R>(),
              target_rank, disp, target_count,
              TypeMap<R>(), window));
#endif
}

template<typename R>
void Iput (const Complex<R>* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Iput"))
#ifdef EL_ENSURE_PUT_ATOMICITY
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi (MPI_Accumulate
             (source, 2*origin_count, TypeMap<R>(),
              target_rank, disp, 2*target_count,
              TypeMap<R>(), MPI_REPLACE, window));
#else
    SafeMpi (MPI_Accumulate
             (source, origin_count, TypeMap<Complex<R>>(),
              target_rank, disp, target_count,
              TypeMap<Complex<R>>(), MPI_REPLACE, window));
#endif
#else
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi (MPI_Put
             (source, 2*origin_count, TypeMap<R>(),
              target_rank, disp, 2*target_count,
              TypeMap<R>(), window));
#else
    SafeMpi (MPI_Put
             (source, origin_count, TypeMap<Complex<R>>(),
              target_rank, disp, target_count,
              TypeMap<Complex<R>>(), window));
#endif
#endif
}

template<typename R>
void Rput (const R* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window,
           Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Rput"))
#ifdef EL_ENSURE_PUT_ATOMICITY
    SafeMpi (MPI_Raccumulate
             (source, origin_count, TypeMap<R>(),
              target_rank, disp, target_count,
              TypeMap<R>(), MPI_REPLACE, window, &request));
#else
    SafeMpi (MPI_Rput
             (source, origin_count, TypeMap<R>(),
              target_rank, disp, target_count,
              TypeMap<R>(), window, &request));
#endif
}

template<typename R>
void Rput (const Complex<R>* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window,
           Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Rput"))
#ifdef EL_ENSURE_PUT_ATOMICITY
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi (MPI_Raccumulate
             (source, 2*origin_count, TypeMap<R>(),
              target_rank, disp, 2*target_count,
              TypeMap<R>(), MPI_REPLACE, window, &request));
#else
    SafeMpi (MPI_Raccumulate
             (source, origin_count, TypeMap<Complex<R>>(),
              target_rank, disp, target_count,
              TypeMap<Complex<R>>(), MPI_REPLACE, window, &request));
#endif
#else
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi (MPI_Rput
             (source, 2*origin_count, TypeMap<R>(),
              target_rank, disp, 2*target_count,
              TypeMap<R>(), window, &request));
#else
    SafeMpi (MPI_Rput
             (source, origin_count, TypeMap<Complex<R>>(),
              target_rank, disp, target_count,
              TypeMap<Complex<R>>(), window, &request));
#endif
#endif
}
template void Iput (const byte* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iput (const int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iput (const unsigned* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iput (const long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iput (const unsigned long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Iput (const long long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iput (const unsigned long long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
#endif
template void Iput (const float* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iput (const double* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iput (const Complex<double>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iput (const Complex<float>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);

template void Rput (const byte* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rput (const int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rput (const unsigned* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rput (const long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rput (const unsigned long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Rput (const long long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rput (const unsigned long long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
#endif
template void Rput (const float* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rput (const double* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rput (const Complex<double>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rput (const Complex<float>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);

// when source-target size == 1
template<typename T>
void Iput( const T source, int target_rank, Aint disp, Window& window )
{
    Iput ( &source, 1, target_rank, disp, 1, window );
}

template<typename T>
void Rput( const T source, int target_rank, Aint disp,
           Window& window, Request& request )
{
    Rput ( &source, 1, target_rank, disp, 1, window, request );
}

template void Rput (const byte source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rput (const int source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rput (const unsigned source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rput (const long int source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rput (const unsigned long source, int target_rank,
                    Aint disp, Window & window, Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Rput (const long long int source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rput (const unsigned long long source, int target_rank,
                    Aint disp, Window & window, Request & request);
#endif
template void Rput (const float source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rput (const double source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rput (const Complex<double> source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rput (const Complex<float> source, int target_rank,
                    Aint disp, Window & window, Request & request);

template void Iput (const byte source, int target_rank,
                    Aint disp, Window & window);
template void Iput (const int source, int target_rank,
                    Aint disp, Window & window);
template void Iput (const unsigned source, int target_rank,
                    Aint disp, Window & window);
template void Iput (const long int source, int target_rank,
                    Aint disp, Window & window);
template void Iput (const unsigned long source, int target_rank,
                    Aint disp, Window & window);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Iput (const long long int source, int target_rank,
                    Aint disp, Window & window);
template void Iput (const unsigned long long source, int target_rank,
                    Aint disp, Window & window);
#endif
template void Iput (const float source, int target_rank,
                    Aint disp, Window & window);
template void Iput (const double source, int target_rank,
                    Aint disp, Window & window);
template void Iput (const Complex<double> source, int target_rank,
                    Aint disp, Window & window);
template void Iput (const Complex<float> source, int target_rank,
                    Aint disp, Window & window);
// get
template<typename R>
void Iget (R* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Iget"))
#ifdef EL_ENSURE_GET_ATOMICITY
    SafeMpi (MPI_Get_accumulate
             (NULL, 0, TypeMap<R>(), source,
              origin_count, TypeMap<R>(),
              target_rank, disp, target_count,
              TypeMap<R>(), MPI_NO_OP, window));
#else
    SafeMpi (MPI_Get
             (source, origin_count, TypeMap<R>(),
              target_rank, disp, target_count,
              TypeMap<R>(), window));
#endif
}

template<typename R>
void Iget (Complex<R>* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Iget"))
#ifdef EL_ENSURE_GET_ATOMICITY
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi (MPI_Get_accumulate
             (NULL, 0, TypeMap<R>(), source,
              2*origin_count, TypeMap<R>(),
              target_rank, disp, 2*target_count,
              TypeMap<R>(), MPI_NO_OP, window));
#else
    SafeMpi (MPI_Get_accumulate
             (NULL, 0, TypeMap<Complex<R>>(), source,
              origin_count, TypeMap<Complex<R>>(),
              target_rank, disp, target_count,
              TypeMap<Complex<R>>(), MPI_NO_OP, window));
#endif
#else
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi (MPI_Get
             (source, 2*origin_count, TypeMap<R>(),
              target_rank, disp, 2*target_count,
              TypeMap<R>(), window));
#else
    SafeMpi (MPI_Get
             (source, origin_count, TypeMap<Complex<R>>(),
              target_rank, disp, target_count,
              TypeMap<Complex<R>>(), window));
#endif
#endif
}

template<typename R>
void Rget (R* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window,
           Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Rget"))
#ifdef EL_ENSURE_GET_ATOMICITY
    SafeMpi (MPI_Rget_accumulate
             (NULL, 0, TypeMap<R>(), source,
              origin_count, TypeMap<R>(),
              target_rank, disp, target_count,
              TypeMap<R>(), MPI_NO_OP, window,
              &request));
#else
    SafeMpi (MPI_Rget
             (source, origin_count, TypeMap<R>(),
              target_rank, disp, target_count,
              TypeMap<R>(), window, &request));
#endif
}

template<typename R>
void Rget (Complex<R>* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window,
           Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Rget"))
#ifdef EL_ENSURE_GET_ATOMICITY
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi (MPI_Rget_accumulate
             (NULL, 0, TypeMap<R>(), source,
              2*origin_count, TypeMap<R>(),
              target_rank, disp, 2*target_count,
              TypeMap<R>(), MPI_NO_OP, window, &request));
#else
    SafeMpi (MPI_Rget_accumulate
             (NULL, 0, TypeMap<Complex<R>>(), source,
              origin_count, TypeMap<Complex<R>>(),
              target_rank, disp, target_count,
              TypeMap<Complex<R>>(), MPI_NO_OP, window, &request));
#endif
#else
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi (MPI_Rget
             (source, 2*origin_count, TypeMap<R>(),
              target_rank, disp, 2*target_count,
              TypeMap<R>(), window, &request));
#else
    SafeMpi (MPI_Rget
             (source, origin_count, TypeMap<Complex<R>>(),
              target_rank, disp, target_count,
              TypeMap<Complex<R>>(), window, &request));
#endif
#endif
}
template void Iget (byte* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iget (int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iget (unsigned* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iget (long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iget (unsigned long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Iget (long long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iget (unsigned long long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
#endif
template void Iget (float* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iget (double* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iget (Complex<double>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iget (Complex<float>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);

template void Rget (byte* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rget (int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rget (unsigned* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rget (long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rget (unsigned long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Rget (long long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rget (unsigned long long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
#endif
template void Rget (float* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rget (double* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rget (Complex<double>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Rget (Complex<float>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);

// when source-target size == 1
template<typename T>
void Iget( T source, int target_rank, Aint disp, Window& window )
{
    Iget ( &source, 1, target_rank, disp, 1, window );
}

template<typename T>
void Rget( T source, int target_rank, Aint disp,
           Window& window, Request& request )
{
    Rget ( &source, 1, target_rank, disp, 1, window, request );
}

template void Rget (byte source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rget (int source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rget (unsigned source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rget (long int source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rget (unsigned long source, int target_rank,
                    Aint disp, Window & window, Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Rget (long long int source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rget (unsigned long long source, int target_rank,
                    Aint disp, Window & window, Request & request);
#endif
template void Rget (float source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rget (double source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rget (Complex<double> source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Rget (Complex<float> source, int target_rank,
                    Aint disp, Window & window, Request & request);

template void Iget (byte source, int target_rank,
                    Aint disp, Window & window);
template void Iget (int source, int target_rank,
                    Aint disp, Window & window);
template void Iget (unsigned source, int target_rank,
                    Aint disp, Window & window);
template void Iget (long int source, int target_rank,
                    Aint disp, Window & window);
template void Iget (unsigned long source, int target_rank,
                    Aint disp, Window & window);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Iget (long long int source, int target_rank,
                    Aint disp, Window & window);
template void Iget (unsigned long long source, int target_rank,
                    Aint disp, Window & window);
#endif
template void Iget (float source, int target_rank,
                    Aint disp, Window & window);
template void Iget (double source, int target_rank,
                    Aint disp, Window & window);
template void Iget (Complex<double> source, int target_rank,
                    Aint disp, Window & window);
template void Iget (Complex<float> source, int target_rank,
                    Aint disp, Window & window);

// acc
template<typename R>
void Iacc (const R* source, int origin_count, int target_rank,
           Aint disp, int target_count, Op op, Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Iaccumulate"))
	
    SafeMpi (MPI_Accumulate
             (source, origin_count,
              TypeMap<R>(), target_rank, disp,
              target_count, TypeMap<R>(), op.op,
              window));
}

template<typename R>
void Iacc (const Complex<R>* source, int origin_count, int target_rank,
           Aint disp, int target_count, Op op, Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Iaccumulate"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi (MPI_Accumulate
             (source, 2*origin_count,
              TypeMap<R>(), target_rank, disp,
              2*target_count, TypeMap<R>(), op.op,
              window));
#else
    SafeMpi (MPI_Accumulate
             (source, origin_count,
              TypeMap<Complex<R>>(), target_rank, disp,
              target_count, TypeMap<Complex<R>>(), op.op,
              window));
#endif
}

template<typename R>
void Racc (const R* source, int origin_count, int target_rank,
           Aint disp, int target_count, Op op, Window & window,
           Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Raccumulate"))
        SafeMpi (MPI_Raccumulate
                 (source, origin_count,
                  TypeMap<R>(), target_rank, disp,
                  target_count, TypeMap<R>(), op.op,
                  window, &request));
}

template<typename R>
void Racc (const Complex<R>* source, int origin_count, int target_rank,
           Aint disp, int target_count, Op op, Window & window,
           Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Raccumulate"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi (MPI_Raccumulate
             (source, 2*origin_count,
              TypeMap<R>(), target_rank, disp,
              2*target_count, TypeMap<R>(), op.op,
              window, &request));
#else
    SafeMpi (MPI_Raccumulate
             (source, origin_count,
              TypeMap<Complex<R>>(), target_rank, disp,
              target_count, TypeMap<Complex<R>>(), op.op,
              window, &request));
#endif
}
template void Iacc (const byte* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window);
template void Iacc (const int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window);
template void Iacc (const unsigned* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window);
template void Iacc (const long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window);
template void Iacc (const unsigned long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Iacc (const long long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window);
template void Iacc (const unsigned long long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window);
#endif
template void Iacc (const float* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window);
template void Iacc (const double* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window);
template void Iacc (const Complex<double>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window);
template void Iacc (const Complex<float>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window);

template void Racc (const byte* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window, Request & request);
template void Racc (const int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window, Request & request);
template void Racc (const unsigned* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window, Request & request);
template void Racc (const long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window, Request & request);
template void Racc (const unsigned long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window, Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Racc (const long long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window, Request & request);
template void Racc (const unsigned long long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window, Request & request);
#endif
template void Racc (const float* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window, Request & request);
template void Racc (const double* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window, Request & request);
template void Racc (const Complex<double>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window, Request & request);
template void Racc (const Complex<float>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Op op, Window & window, Request & request);

// op = SUM
template<typename R>
void Iacc (const R* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window)
{
    Iacc ( source, origin_count, target_rank, disp, target_count, SUM, window );
}

template<typename R>
void Racc (const R* source, int origin_count, int target_rank,
           Aint disp, int target_count, Window & window,
           Request & request)
{
    Racc ( source, origin_count, target_rank, disp, target_count, SUM, window, request );
}

template void Iacc (const byte* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iacc (const int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iacc (const unsigned* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iacc (const long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iacc (const unsigned long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Iacc (const long long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iacc (const unsigned long long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
#endif
template void Iacc (const float* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iacc (const double* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iacc (const Complex<double>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);
template void Iacc (const Complex<float>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window);

template void Racc (const byte* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Racc (const int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Racc (const unsigned* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Racc (const long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Racc (const unsigned long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Racc (const long long int* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Racc (const unsigned long long* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
#endif
template void Racc (const float* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Racc (const double* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Racc (const Complex<double>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);
template void Racc (const Complex<float>* source, int origin_count, int target_rank,
                    Aint disp, int target_count, Window & window, Request & request);

// when source-target size == 1 and op = SUM
template<typename T>
void Iacc (const T source, int target_rank, Aint disp, Window & window)
{
    Iacc ( &source, 1, target_rank, disp, 1, SUM, window );
}

template<typename T>
void Racc (const T source, int target_rank, Aint disp, Window & window,
           Request & request)
{
    Racc ( &source, 1, target_rank, disp, 1, SUM, window, request );
}

template void Racc (const byte source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Racc (const int source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Racc (const unsigned source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Racc (const long int source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Racc (const unsigned long source, int target_rank,
                    Aint disp, Window & window, Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Racc (const long long int source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Racc (const unsigned long long source, int target_rank,
                    Aint disp, Window & window, Request & request);
#endif
template void Racc (const float source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Racc (const double source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Racc (const Complex<double> source, int target_rank,
                    Aint disp, Window & window, Request & request);
template void Racc (const Complex<float> source, int target_rank,
                    Aint disp, Window & window, Request & request);

template void Iacc (const byte source, int target_rank,
                    Aint disp, Window & window);
template void Iacc (const int source, int target_rank,
                    Aint disp, Window & window);
template void Iacc (const unsigned source, int target_rank,
                    Aint disp, Window & window);
template void Iacc (const long int source, int target_rank,
                    Aint disp, Window & window);
template void Iacc (const unsigned long source, int target_rank,
                    Aint disp, Window & window);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Iacc (const long long int source, int target_rank,
                    Aint disp, Window & window);
template void Iacc (const unsigned long long source, int target_rank,
                    Aint disp, Window & window);
#endif
template void Iacc (const float source, int target_rank,
                    Aint disp, Window & window);
template void Iacc (const double source, int target_rank,
                    Aint disp, Window & window);
template void Iacc (const Complex<double> source, int target_rank,
                    Aint disp, Window & window);
template void Iacc (const Complex<float> source, int target_rank,
                    Aint disp, Window & window);

// Synchronization
// ---------------
void Flush (int target_rank, Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Flush"))
    SafeMpi (MPI_Win_flush (target_rank, window));
}

void Flush (Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Flush"))
    SafeMpi (MPI_Win_flush_all (window));
}

void FlushLocal (int target_rank, Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::FlushLocal"))
    SafeMpi (MPI_Win_flush_local (target_rank, window));
}

void FlushLocal (Window & window)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::FlushLocal"))
    SafeMpi (MPI_Win_flush_local_all (window));
}
#endif // EL_ENABLE_RMA_AXPY

// Various utilities
// =================
// Free request
void RequestFree (Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::RequestFree"))
    SafeMpi (MPI_Request_free (&request));
}

// Wait until every process in comm reaches this statement
void Barrier (Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Barrier"))
    SafeMpi (MPI_Barrier (comm.comm));
}

#if MPI_VERSION>=3 && defined(EL_USE_IBARRIER_FOR_AXPY)
void IBarrier (Comm comm, Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::IBarrier"))
    SafeMpi (MPI_Ibarrier (comm.comm, &request));
}
#endif


// Test for completion
bool Test (Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Test")) 
    
    Status status;
    int flag;

    SafeMpi (MPI_Test (&request, &flag, &status));
    if (flag)
	return true;
    else
	return false;
}

bool Test (Request & request, Status & status)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Test")) 
    int flag;
    SafeMpi (MPI_Test (&request, &flag, &status));

    if (flag)
	return true;
    else
	return false;
}

bool Testany (int count, Request * requests, int &indx,
              Status & status)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Testany")) 
    int flag;
    SafeMpi (MPI_Testany
             (count, requests, &indx, &flag, &status));
     if (flag)
	return true;
    else
	return false;
}

bool Testany (int count, Request * requests, int &indx)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Testany")) 
    
    int flag;
    Status status;

    SafeMpi (MPI_Testany
             (count, requests, &indx, &flag, &status)); 
    if (flag)
	return true;
    else
	return false;
}

bool Testany (int count, Request * requests)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Testany")) 

    int flag, indx;
    Status status;

    SafeMpi (MPI_Testany
             (count, requests, &indx, &flag, &status));
    if (flag)
	return true;
    else
	return false;
}

// Ensure that the request finishes before continuing
void Wait (Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Wait")) 
    SafeMpi (MPI_Wait (&request, MPI_STATUS_IGNORE));
}

// Ensure that the request finishes before continuing
void Wait (Request & request, Status & status)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Wait"))
    SafeMpi (MPI_Wait (&request, &status));
}

// Ensure that several requests finish before continuing
void WaitAll (int numRequests, Request * requests)
{

    DEBUG_ONLY(CallStackEntry cse("mpi::WaitAll"))
    vector<Status> statuses( numRequests );
    SafeMpi( MPI_Waitall( numRequests, requests, statuses.data() ) );
}

// Ensure that several requests finish before continuing
void WaitAll (int numRequests, Request * requests,
              Status * statuses)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::WaitAll"))
    SafeMpi (MPI_Waitall
             (numRequests, requests, statuses));
}

// Ensure that any requests finish before continuing
void WaitAny (int numRequests, Request * requests, Int * index)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::WaitAny"))
    SafeMpi (MPI_Waitany
             (numRequests, requests, index, MPI_STATUS_IGNORE));
}

// Nonblocking test for message completion
bool IProbe (int source, int tag, Comm comm,
             Status & status)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::IProbe")) int
    flag;
    SafeMpi (MPI_Iprobe
             (source, tag, comm.comm, &flag, &status));
    return flag;
}

bool IProbe (int source, Comm comm, Status & status)
{
    return IProbe (source, mpi::ANY_TAG, comm, status);
}

bool IProbe (Comm comm, Status & status)
{
    return IProbe (mpi::ANY_SOURCE, mpi::ANY_TAG, comm, status);
}

void Probe (int source, int tag, Comm comm, Status & status)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Probe"))
    SafeMpi (MPI_Probe(source, tag, comm.comm, &status));
}

void Probe (int source, Comm comm, Status & status)
{
    Probe (source, mpi::ANY_TAG, comm, status);
}

void Probe (Comm comm, Status & status)
{
    Probe (mpi::ANY_SOURCE, mpi::ANY_TAG, comm, status);
}

bool IMprobe (int source, int tag, Comm comm,
             Status & status, Message & message)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::IMprobe")) 
    int flag;
    SafeMpi (MPI_Improbe
             (source, tag, comm.comm, &flag, &message, &status));
    return flag;
}

template < typename T > int GetCount (Status & status)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::GetCount")) int
    count;
    SafeMpi (MPI_Get_count
             (&status, TypeMap < T > (), &count));
    return count;
}
template int GetCount < byte > (Status & status);
template int GetCount < int >(Status & status);
template int GetCount < unsigned >(Status & status);
template int GetCount < long int >(Status & status);
template int GetCount < unsigned long >(Status & status);

#ifdef EL_HAVE_MPI_LONG_LONG
template int GetCount < long long int >(Status & status);
template int GetCount <
unsigned long long >(Status & status);
#endif
template int GetCount < float >(Status & status);
template int GetCount < double >(Status & status);
template int GetCount < Complex <
float >>(Status & status);
template int GetCount < Complex <
double >>(Status & status);

template < typename R >
void TaggedSend (const R * buf, int count, int to,
                 int tag, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Send"))
    SafeMpi (MPI_Send
             (const_cast < R * >(buf), count,
              TypeMap < R > (), to, tag, comm.comm));
}

template < typename R >
void TaggedSend (const Complex < R > *buf, int count,
                 int to, int tag, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Send"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Send
     (const_cast < Complex < R > *>(buf), 2 * count,
      TypeMap < R > (), to, tag, comm.comm));
#else
    SafeMpi
    (MPI_Send
     (const_cast < Complex < R > *>(buf), count,
      TypeMap < Complex < R >> (), to, tag,
      comm.comm));
#endif
}

template void TaggedSend (const byte * buf, int count,
                          int to, int tag, Comm comm);
template void TaggedSend (const int *buf, int count,
                          int to, int tag, Comm comm);
template void TaggedSend (const unsigned *buf, int count,
                          int to, int tag, Comm comm);
template void TaggedSend (const long int *buf, int count,
                          int to, int tag, Comm comm);
template void TaggedSend (const unsigned long *buf,
                          int count, int to, int tag,
                          Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedSend (const long long int *buf,
                          int count, int to, int tag,
                          Comm comm);
template void TaggedSend (const unsigned long long *buf,
                          int count, int to, int tag,
                          Comm comm);
#endif
template void TaggedSend (const float *buf, int count,
                          int to, int tag, Comm comm);
template void TaggedSend (const double *buf, int count,
                          int to, int tag, Comm comm);
template void TaggedSend (const Complex < float >*buf,
                          int count, int to, int tag,
                          Comm comm);
template void TaggedSend (const Complex < double >*buf,
                          int count, int to, int tag,
                          Comm comm);

template < typename T >
void Send (const T * buf, int count, int to,
           Comm comm)
{
    TaggedSend (buf, count, to, 0, comm);
}

template void Send (const byte * buf, int count, int to,
                    Comm comm);
template void Send (const int *buf, int count, int to,
                    Comm comm);
template void Send (const unsigned *buf, int count,
                    int to, Comm comm);
template void Send (const long int *buf, int count,
                    int to, Comm comm);
template void Send (const unsigned long *buf, int count,
                    int to, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Send (const long long int *buf, int count,
                    int to, Comm comm);
template void Send (const unsigned long long *buf,
                    int count, int to, Comm comm);
#endif
template void Send (const float *buf, int count, int to,
                    Comm comm);
template void Send (const double *buf, int count, int to,
                    Comm comm);
template void Send (const Complex < float >*buf,
                    int count, int to, Comm comm);
template void Send (const Complex < double >*buf,
                    int count, int to, Comm comm);

template < typename T >
void TaggedSend (T b, int to, int tag, Comm comm)
{
    TaggedSend (&b, 1, to, tag, comm);
}

template void TaggedSend (byte b, int to, int tag,
                          Comm comm);
template void TaggedSend (int b, int to, int tag,
                          Comm comm);
template void TaggedSend (unsigned b, int to, int tag,
                          Comm comm);
template void TaggedSend (long int b, int to, int tag,
                          Comm comm);
template void TaggedSend (unsigned long b, int to,
                          int tag, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedSend (long long int b, int to,
                          int tag, Comm comm);
template void TaggedSend (unsigned long long b, int to,
                          int tag, Comm comm);
#endif
template void TaggedSend (float b, int to, int tag,
                          Comm comm);
template void TaggedSend (double b, int to, int tag,
                          Comm comm);
template void TaggedSend (Complex < float >b, int to,
                          int tag, Comm comm);
template void TaggedSend (Complex < double >b, int to,
                          int tag, Comm comm);

template < typename T > void Send (T b, int to, Comm comm)
{
    TaggedSend (b, to, 0, comm);
}

template void Send (byte b, int to, Comm comm);
template void Send (int b, int to, Comm comm);
template void Send (unsigned b, int to, Comm comm);
template void Send (long int b, int to, Comm comm);
template void Send (unsigned long b, int to, Comm comm);

#ifdef EL_HAVE_MPI_LONG_LONG
template void Send (long long int b, int to, Comm comm);
template void Send (unsigned long long b, int to,
                    Comm comm);
#endif
template void Send (float b, int to, Comm comm);
template void Send (double b, int to, Comm comm);
template void Send (Complex < float >b, int to,
                    Comm comm);
template void Send (Complex < double >b, int to,
                    Comm comm);

template < typename R >
void TaggedISend
(const R * buf, int count, int to, int tag, Comm comm,
 Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::TaggedISend"))
    SafeMpi
    (MPI_Isend
     (const_cast < R * >(buf), count,
      TypeMap < R > (), to, tag, comm.comm,
      &request));
}

template < typename R >
void TaggedISend
(const Complex < R > *buf, int count, int to, int tag,
 Comm comm, Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::TaggedISend"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Isend
     (const_cast < Complex < R > *>(buf), 2 * count,
      TypeMap < R > (), to, tag, comm.comm,
      &request));
#else
    SafeMpi
    (MPI_Isend
     (const_cast < Complex < R > *>(buf), count,
      TypeMap < Complex < R >> (), to, tag, comm.comm,
      &request));
#endif
}

template void TaggedISend (const byte * buf, int count,
                           int to, int tag, Comm comm,
                           Request & request);
template void TaggedISend (const int *buf, int count,
                           int to, int tag, Comm comm,
                           Request & request);
template void TaggedISend (const unsigned *buf, int count,
                           int to, int tag, Comm comm,
                           Request & request);
template void TaggedISend (const long int *buf, int count,
                           int to, int tag, Comm comm,
                           Request & request);
template void TaggedISend (const unsigned long *buf,
                           int count, int to, int tag,
                           Comm comm, Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedISend (const long long int *buf,
                           int count, int to, int tag,
                           Comm comm, Request & request);
template void TaggedISend (const unsigned long long *buf,
                           int count, int to, int tag,
                           Comm comm, Request & request);
#endif
template void TaggedISend (const float *buf, int count,
                           int to, int tag, Comm comm,
                           Request & request);
template void TaggedISend (const double *buf, int count,
                           int to, int tag, Comm comm,
                           Request & request);
template void TaggedISend (const Complex < float >*buf,
                           int count, int to, int tag,
                           Comm comm, Request & request);
template void TaggedISend (const Complex < double >*buf,
                           int count, int to, int tag,
                           Comm comm, Request & request);

template < typename T >
void ISend
(const T * buf, int count, int to, Comm comm,
 Request & request)
{
    TaggedISend (buf, count, to, 0, comm, request);
}

template void ISend (const byte * buf, int count, int to,
                     Comm comm, Request & request);
template void ISend (const int *buf, int count, int to,
                     Comm comm, Request & request);
template void ISend (const unsigned *buf, int count,
                     int to, Comm comm,
                     Request & request);
template void ISend (const long int *buf, int count,
                     int to, Comm comm,
                     Request & request);
template void ISend (const unsigned long *buf, int count,
                     int to, Comm comm,
                     Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void ISend (const long long int *buf, int count,
                     int to, Comm comm,
                     Request & request);
template void ISend (const unsigned long long *buf,
                     int count, int to, Comm comm,
                     Request & request);
#endif
template void ISend (const float *buf, int count, int to,
                     Comm comm, Request & request);
template void ISend (const double *buf, int count, int to,
                     Comm comm, Request & request);
template void ISend (const Complex < float >*buf,
                     int count, int to, Comm comm,
                     Request & request);
template void ISend (const Complex < double >*buf,
                     int count, int to, Comm comm,
                     Request & request);

template < typename T >
void TaggedISend (T b, int to, int tag, Comm comm,
                  Request & request)
{
    TaggedISend (&b, 1, to, tag, comm, request);
}

template void TaggedISend (byte buf, int to, int tag,
                           Comm comm, Request & request);
template void TaggedISend (int buf, int to, int tag,
                           Comm comm, Request & request);
template void TaggedISend (unsigned buf, int to, int tag,
                           Comm comm, Request & request);
template void TaggedISend (long int buf, int to, int tag,
                           Comm comm, Request & request);
template void TaggedISend (unsigned long buf, int to,
                           int tag, Comm comm,
                           Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedISend (long long int buf, int to,
                           int tag, Comm comm,
                           Request & request);
template void TaggedISend (unsigned long long buf, int to,
                           int tag, Comm comm,
                           Request & request);
#endif
template void TaggedISend (float buf, int to, int tag,
                           Comm comm, Request & request);
template void TaggedISend (double buf, int to, int tag,
                           Comm comm, Request & request);
template void TaggedISend (Complex < float >buf, int to,
                           int tag, Comm comm,
                           Request & request);
template void TaggedISend (Complex < double >buf, int to,
                           int tag, Comm comm,
                           Request & request);

template < typename T >
void ISend (T b, int to, Comm comm, Request & request)
{
    TaggedISend (b, to, 0, comm, request);
}

template void ISend (byte buf, int to, Comm comm,
                     Request & request);
template void ISend (int buf, int to, Comm comm,
                     Request & request);
template void ISend (unsigned buf, int to, Comm comm,
                     Request & request);
template void ISend (long int buf, int to, Comm comm,
                     Request & request);
template void ISend (unsigned long buf, int to, Comm comm,
                     Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void ISend (long long int buf, int to, Comm comm,
                     Request & request);
template void ISend (unsigned long long buf, int to,
                     Comm comm, Request & request);
#endif
template void ISend (float buf, int to, Comm comm,
                     Request & request);
template void ISend (double buf, int to, Comm comm,
                     Request & request);
template void ISend (Complex < float >buf, int to,
                     Comm comm, Request & request);
template void ISend (Complex < double >buf, int to,
                     Comm comm, Request & request);

template < typename R >
void TaggedISSend
(const R * buf, int count, int to, int tag, Comm comm,
 Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::ISSend"))

    SafeMpi
    (MPI_Issend
     (const_cast < R * >(buf), count,
      TypeMap < R > (), to, tag, comm.comm,
      &request));
}

template < typename R >
void TaggedISSend
(const Complex < R > *buf, int count, int to, int tag,
 Comm comm, Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::ISSend"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Issend
     (const_cast < Complex < R > *>(buf), 2 * count,
      TypeMap < R > (), to, tag, comm.comm,
      &request));
#else
    SafeMpi
    (MPI_Issend
     (const_cast < Complex < R > *>(buf), count,
      TypeMap < Complex < R >> (), to, tag, comm.comm,
      &request));
#endif
}

template void TaggedISSend (const byte * buf, int count,
                            int to, int tag, Comm comm,
                            Request & request);
template void TaggedISSend (const int *buf, int count,
                            int to, int tag, Comm comm,
                            Request & request);
template void TaggedISSend (const unsigned *buf,
                            int count, int to, int tag,
                            Comm comm, Request & request);
template void TaggedISSend (const long int *buf,
                            int count, int to, int tag,
                            Comm comm, Request & request);
template void TaggedISSend (const unsigned long *buf,
                            int count, int to, int tag,
                            Comm comm, Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedISSend (const long long int *buf,
                            int count, int to, int tag,
                            Comm comm, Request & request);
template void TaggedISSend (const unsigned long long *buf,
                            int count, int to, int tag,
                            Comm comm, Request & request);
#endif
template void TaggedISSend (const float *buf, int count,
                            int to, int tag, Comm comm,
                            Request & request);
template void TaggedISSend (const double *buf, int count,
                            int to, int tag, Comm comm,
                            Request & request);
template void TaggedISSend (const Complex < float >*buf,
                            int count, int to, int tag,
                            Comm comm, Request & request);
template void TaggedISSend (const Complex < double >*buf,
                            int count, int to, int tag,
                            Comm comm, Request & request);

template < typename T >
void ISSend (const T * buf, int count, int to,
             Comm comm, Request & request)
{
    TaggedISSend (buf, count, to, 0, comm, request);
}

template void ISSend (const byte * buf, int count, int to,
                      Comm comm, Request & request);
template void ISSend (const int *buf, int count, int to,
                      Comm comm, Request & request);
template void ISSend (const unsigned *buf, int count,
                      int to, Comm comm,
                      Request & request);
template void ISSend (const long int *buf, int count,
                      int to, Comm comm,
                      Request & request);
template void ISSend (const unsigned long *buf, int count,
                      int to, Comm comm,
                      Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void ISSend (const long long int *buf, int count,
                      int to, Comm comm,
                      Request & request);
template void ISSend (const unsigned long long *buf,
                      int count, int to, Comm comm,
                      Request & request);
#endif
template void ISSend (const float *buf, int count, int to,
                      Comm comm, Request & request);
template void ISSend (const double *buf, int count,
                      int to, Comm comm,
                      Request & request);
template void ISSend (const Complex < float >*buf,
                      int count, int to, Comm comm,
                      Request & request);
template void ISSend (const Complex < double >*buf,
                      int count, int to, Comm comm,
                      Request & request);

template < typename T >
void TaggedISSend (T b, int to, int tag, Comm comm,
                   Request & request)
{
    TaggedISSend (&b, 1, to, tag, comm, request);
}

template void TaggedISSend (byte b, int to, int tag,
                            Comm comm, Request & request);
template void TaggedISSend (int b, int to, int tag,
                            Comm comm, Request & request);
template void TaggedISSend (unsigned b, int to, int tag,
                            Comm comm, Request & request);
template void TaggedISSend (long int b, int to, int tag,
                            Comm comm, Request & request);
template void TaggedISSend (unsigned long b, int to,
                            int tag, Comm comm,
                            Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedISSend (long long int b, int to,
                            int tag, Comm comm,
                            Request & request);
template void TaggedISSend (unsigned long long b, int to,
                            int tag, Comm comm,
                            Request & request);
#endif
template void TaggedISSend (float b, int to, int tag,
                            Comm comm, Request & request);
template void TaggedISSend (double b, int to, int tag,
                            Comm comm, Request & request);
template void TaggedISSend (Complex < float >b, int to,
                            int tag, Comm comm,
                            Request & request);
template void TaggedISSend (Complex < double >b, int to,
                            int tag, Comm comm,
                            Request & request);

template < typename R >
void TaggedRecv (R * buf, int count, int from,
                 int tag, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Recv")) Status
    status;
    SafeMpi (MPI_Recv
             (buf, count, TypeMap < R > (), from, tag,
              comm.comm, &status));
}

template < typename R >
void TaggedRecv (Complex < R > *buf, int count,
                 int from, int tag, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Recv")) Status
    status;
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Recv
     (buf, 2 * count, TypeMap < R > (), from, tag,
      comm.comm, &status));
#else
    SafeMpi
    (MPI_Recv
     (buf, count, TypeMap < Complex < R >> (), from,
      tag, comm.comm, &status));
#endif
}

template void TaggedRecv (byte * buf, int count, int from,
                          int tag, Comm comm);
template void TaggedRecv (int *buf, int count, int from,
                          int tag, Comm comm);
template void TaggedRecv (unsigned *buf, int count,
                          int from, int tag, Comm comm);
template void TaggedRecv (long int *buf, int count,
                          int from, int tag, Comm comm);
template void TaggedRecv (unsigned long *buf, int count,
                          int from, int tag, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedRecv (long long int *buf, int count,
                          int from, int tag, Comm comm);
template void TaggedRecv (unsigned long long *buf,
                          int count, int from, int tag,
                          Comm comm);
#endif
template void TaggedRecv (float *buf, int count, int from,
                          int tag, Comm comm);
template void TaggedRecv (double *buf, int count,
                          int from, int tag, Comm comm);
template void TaggedRecv (Complex < float >*buf,
                          int count, int from, int tag,
                          Comm comm);
template void TaggedRecv (Complex < double >*buf,
                          int count, int from, int tag,
                          Comm comm);

template < typename R >
void TaggedRecvS (R * buf, int count, int from,
                 int tag, Comm comm, Status & status)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Recv")) 
    SafeMpi (MPI_Recv
             (buf, count, TypeMap < R > (), from, tag,
              comm.comm, &status));
}

template < typename R >
void TaggedRecvS (Complex < R > *buf, int count,
                 int from, int tag, Comm comm, Status & status)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Recv"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Recv
     (buf, 2 * count, TypeMap < R > (), from, tag,
      comm.comm, &status));
#else
    SafeMpi
    (MPI_Recv
     (buf, count, TypeMap < Complex < R >> (), from,
      tag, comm.comm, &status));
#endif
}

template void TaggedRecvS (byte * buf, int count, int from,
                          int tag, Comm comm, Status & status);
template void TaggedRecvS (int *buf, int count, int from,
                          int tag, Comm comm, Status & status);
template void TaggedRecvS (unsigned *buf, int count,
                          int from, int tag, Comm comm, Status & status);
template void TaggedRecvS (long int *buf, int count,
                          int from, int tag, Comm comm, Status & status);
template void TaggedRecvS (unsigned long *buf, int count,
                          int from, int tag, Comm comm, Status & status);
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedRecvS (long long int *buf, int count,
                          int from, int tag, Comm comm, Status & status);
template void TaggedRecvS (unsigned long long *buf,
                          int count, int from, int tag,
                          Comm comm, Status & status);
#endif
template void TaggedRecvS (float *buf, int count, int from,
                          int tag, Comm comm, Status & status);
template void TaggedRecvS (double *buf, int count,
                          int from, int tag, Comm comm, Status & status);
template void TaggedRecvS (Complex < float >*buf,
                          int count, int from, int tag,
                          Comm comm, Status & status);
template void TaggedRecvS (Complex < double >*buf,
                          int count, int from, int tag,
                          Comm comm, Status & status);

// matching recv
template < typename R >
void TaggedMrecv (R * buf, int count, Message & msg)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Mrecv")) 
    Status status;
    SafeMpi (MPI_Mrecv
             (buf, count, TypeMap < R > (), 
	      &msg, &status));
}

template < typename R >
void TaggedMrecv (Complex < R > *buf, int count, Message & msg)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Mrecv")) 
    Status status;
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Mrecv
     (buf, 2 * count, TypeMap < R > (), &msg, &status));
#else
    SafeMpi
    (MPI_Mrecv
     (buf, count, TypeMap < Complex < R >> (), 
      &msg, &status));
#endif
}

template void TaggedMrecv (byte * buf, int count, Message & msg);
template void TaggedMrecv (int *buf, int count, Message & msg);
template void TaggedMrecv (unsigned *buf, int count, Message & msg);
template void TaggedMrecv (long int *buf, int count, Message & msg);
template void TaggedMrecv (unsigned long *buf, int count, Message & msg);
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedMrecv (long long int *buf, int count, Message & msg);
template void TaggedMrecv (unsigned long long *buf,
                          int count, Message & msg);
#endif
template void TaggedMrecv (float *buf, int count, Message & msg);
template void TaggedMrecv (double *buf, int count, Message & msg);
template void TaggedMrecv (Complex < float >*buf,
                          int count, Message & msg);
template void TaggedMrecv (Complex < double >*buf,
                          int count, Message & msg);

template < typename T >
void Recv (T * buf, int count, int from, Comm comm)
{
    TaggedRecv (buf, count, from, mpi::ANY_TAG, comm);
}

template void Recv (byte * buf, int count, int from,
                    Comm comm);
template void Recv (int *buf, int count, int from,
                    Comm comm);
template void Recv (unsigned *buf, int count, int from,
                    Comm comm);
template void Recv (long int *buf, int count, int from,
                    Comm comm);
template void Recv (unsigned long *buf, int count,
                    int from, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Recv (long long int *buf, int count,
                    int from, Comm comm);
template void Recv (unsigned long long *buf, int count,
                    int from, Comm comm);
#endif
template void Recv (float *buf, int count, int from,
                    Comm comm);
template void Recv (double *buf, int count, int from,
                    Comm comm);
template void Recv (Complex < float >*buf, int count,
                    int from, Comm comm);
template void Recv (Complex < double >*buf, int count,
                    int from, Comm comm);

template < typename T > T TaggedRecv (int from, int tag,
                                      Comm comm)
{
    T b;

    TaggedRecv (&b, 1, from, tag, comm);
    return b;
}

template byte TaggedRecv (int from, int tag, Comm comm);
template int TaggedRecv (int from, int tag, Comm comm);
template unsigned TaggedRecv (int from, int tag,
                              Comm comm);
template long int TaggedRecv (int from, int tag,
                              Comm comm);
template unsigned long TaggedRecv (int from, int tag,
                                   Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int TaggedRecv (int from, int tag,
                                   Comm comm);
template unsigned long long TaggedRecv (int from, int tag,
                                        Comm comm);
#endif
template float TaggedRecv (int from, int tag, Comm comm);
template double TaggedRecv (int from, int tag, Comm comm);
template Complex < float >TaggedRecv (int from, int tag,
                                      Comm comm);
template Complex < double >TaggedRecv (int from, int tag,
                                       Comm comm);

template < typename T > T Recv (int from, Comm comm)
{
    return TaggedRecv < T > (from, mpi::ANY_TAG, comm);
}

template byte Recv (int from, Comm comm);
template int Recv (int from, Comm comm);
template unsigned Recv (int from, Comm comm);
template long int Recv (int from, Comm comm);
template unsigned long Recv (int from, Comm comm);

#ifdef EL_HAVE_MPI_LONG_LONG
template long long int Recv (int from, Comm comm);
template unsigned long long Recv (int from, Comm comm);
#endif
template float Recv (int from, Comm comm);
template double Recv (int from, Comm comm);
template Complex < float >Recv (int from, Comm comm);
template Complex < double >Recv (int from, Comm comm);

template < typename R >
void TaggedIRecv
(R * buf, int count, int from, int tag, Comm comm,
 Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::TaggedIRecv"))
    SafeMpi
    (MPI_Irecv
     (buf, count, TypeMap < R > (), from, tag,
      comm.comm, &request));
}

template < typename R >
void TaggedIRecv
(Complex < R > *buf, int count, int from, int tag,
 Comm comm, Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::TaggedIRecv"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Irecv
     (buf, 2 * count, TypeMap < R > (), from, tag,
      comm.comm, &request));
#else
    SafeMpi
    (MPI_Irecv
     (buf, count, TypeMap < Complex < R >> (), from,
      tag, comm.comm, &request));
#endif
}

template void TaggedIRecv (byte * buf, int count,
                           int from, int tag, Comm comm,
                           Request & request);
template void TaggedIRecv (int *buf, int count, int from,
                           int tag, Comm comm,
                           Request & request);
template void TaggedIRecv (unsigned *buf, int count,
                           int from, int tag, Comm comm,
                           Request & request);
template void TaggedIRecv (long int *buf, int count,
                           int from, int tag, Comm comm,
                           Request & request);
template void TaggedIRecv (unsigned long *buf, int count,
                           int from, int tag, Comm comm,
                           Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedIRecv (long long int *buf, int count,
                           int from, int tag, Comm comm,
                           Request & request);
template void TaggedIRecv (unsigned long long *buf,
                           int count, int from, int tag,
                           Comm comm, Request & request);
#endif
template void TaggedIRecv (float *buf, int count,
                           int from, int tag, Comm comm,
                           Request & request);
template void TaggedIRecv (double *buf, int count,
                           int from, int tag, Comm comm,
                           Request & request);
template void TaggedIRecv (Complex < float >*buf,
                           int count, int from, int tag,
                           Comm comm, Request & request);
template void TaggedIRecv (Complex < double >*buf,
                           int count, int from, int tag,
                           Comm comm, Request & request);

template < typename T >
void IRecv (T * buf, int count, int from, Comm comm,
            Request & request)
{
    TaggedIRecv (buf, count, from, mpi::ANY_TAG, comm,
                 request);
}

template void IRecv (byte * buf, int count, int from,
                     Comm comm, Request & request);
template void IRecv (int *buf, int count, int from,
                     Comm comm, Request & request);
template void IRecv (unsigned *buf, int count, int from,
                     Comm comm, Request & request);
template void IRecv (long int *buf, int count, int from,
                     Comm comm, Request & request);
template void IRecv (unsigned long *buf, int count,
                     int from, Comm comm,
                     Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void IRecv (long long int *buf, int count,
                     int from, Comm comm,
                     Request & request);
template void IRecv (unsigned long long *buf, int count,
                     int from, Comm comm,
                     Request & request);
#endif
template void IRecv (float *buf, int count, int from,
                     Comm comm, Request & request);
template void IRecv (double *buf, int count, int from,
                     Comm comm, Request & request);
template void IRecv (Complex < float >*buf, int count,
                     int from, Comm comm,
                     Request & request);
template void IRecv (Complex < double >*buf, int count,
                     int from, Comm comm,
                     Request & request);

template < typename T >
T TaggedIRecv (int from, int tag, Comm comm,
               Request & request)
{
    T b;

    TaggedIRecv (&b, 1, from, tag, comm, request);
    return b;
}

template byte TaggedIRecv (int from, int tag, Comm comm,
                           Request & request);
template int TaggedIRecv (int from, int tag, Comm comm,
                          Request & request);
template unsigned TaggedIRecv (int from, int tag,
                               Comm comm,
                               Request & request);
template long int TaggedIRecv (int from, int tag,
                               Comm comm,
                               Request & request);
template unsigned long TaggedIRecv (int from, int tag,
                                    Comm comm,
                                    Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int TaggedIRecv (int from, int tag,
                                    Comm comm,
                                    Request & request);
template unsigned long long TaggedIRecv (int from,
        int tag,
        Comm comm,
        Request &
        request);
#endif
template float TaggedIRecv (int from, int tag, Comm comm,
                            Request & request);
template double TaggedIRecv (int from, int tag, Comm comm,
                             Request & request);
template Complex < float >TaggedIRecv (int from, int tag,
                                       Comm comm,
                                       Request & request);
template Complex < double >TaggedIRecv (int from, int tag,
                                        Comm comm,
                                        Request &
                                        request);

template < typename T >
T IRecv (int from, Comm comm, Request & request)
{
    return TaggedIRecv < T > (from, mpi::ANY_TAG, comm,
                              request);
}

template byte IRecv (int from, Comm comm,
                     Request & request);
template int IRecv (int from, Comm comm,
                    Request & request);
template unsigned IRecv (int from, Comm comm,
                         Request & request);
template long int IRecv (int from, Comm comm,
                         Request & request);
template unsigned long IRecv (int from, Comm comm,
                              Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int IRecv (int from, Comm comm,
                              Request & request);
template unsigned long long IRecv (int from, Comm comm,
                                   Request & request);
#endif
template float IRecv (int from, Comm comm,
                      Request & request);
template double IRecv (int from, Comm comm,
                       Request & request);
template Complex < float >IRecv (int from, Comm comm,
                                 Request & request);
template Complex < double >IRecv (int from, Comm comm,
                                  Request & request);

template < typename R >
void TaggedSendRecv
(const R * sbuf, int sc, int to, int stag,
 R * rbuf, int rc, int from, int rtag, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::SendRecv"))
    Status status;

    SafeMpi (MPI_Sendrecv
             (const_cast < R * >(sbuf), sc,
              TypeMap < R > (), to, stag, rbuf, rc,
              TypeMap < R > (), from, rtag, comm.comm,
              &status));
}

template < typename R >
void TaggedSendRecv
(const Complex < R > *sbuf, int sc, int to, int stag,
 Complex < R > *rbuf, int rc, int from, int rtag,
 Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::SendRecv"))
    Status status;

#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Sendrecv
     (const_cast < Complex < R > *>(sbuf), 2 * sc,
      TypeMap < R > (), to, stag, rbuf, 2 * rc,
      TypeMap < R > (), from, rtag, comm.comm,
      &status));
#else
    SafeMpi
    (MPI_Sendrecv
     (const_cast < Complex < R > *>(sbuf),
      sc, TypeMap < Complex < R >> (), to, stag,
      rbuf,
      rc, TypeMap < Complex < R >> (), from, rtag,
      comm.comm, &status));
#endif
}

template void TaggedSendRecv
(const byte * sbuf, int sc, int to, int stag,
 byte * rbuf, int rc, int from, int rtag, Comm comm);
template void TaggedSendRecv
(const int *sbuf, int sc, int to, int stag,
 int *rbuf, int rc, int from, int rtag, Comm comm);
template void TaggedSendRecv
(const unsigned *sbuf, int sc, int to, int stag,
 unsigned *rbuf, int rc, int from, int rtag,
 Comm comm);
template void TaggedSendRecv (const long int *sbuf,
                              int sc, int to, int stag,
                              long int *rbuf, int rc,
                              int from, int rtag,
                              Comm comm);
template void TaggedSendRecv (const unsigned long *sbuf,
                              int sc, int to, int stag,
                              unsigned long *rbuf, int rc,
                              int from, int rtag,
                              Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedSendRecv
(const long long int *sbuf, int sc, int to, int stag,
 long long int *rbuf, int rc, int from, int rtag,
 Comm comm);
template void TaggedSendRecv (const unsigned long long
                              *sbuf, int sc, int to,
                              int stag,
                              unsigned long long *rbuf,
                              int rc, int from, int rtag,
                              Comm comm);
#endif
template void TaggedSendRecv
(const float *sbuf, int sc, int to, int stag,
 float *rbuf, int rc, int from, int rtag, Comm comm);
template void TaggedSendRecv
(const double *sbuf, int sc, int to, int stag,
 double *rbuf, int rc, int from, int rtag, Comm comm);
template void TaggedSendRecv
(const Complex < float >*sbuf, int sc, int to,
 int stag, Complex < float >*rbuf, int rc, int from,
 int rtag, Comm comm);
template void TaggedSendRecv (const Complex <
                              double >*sbuf, int sc,
                              int to, int stag,
                              Complex < double >*rbuf,
                              int rc, int from, int rtag,
                              Comm comm);

template < typename T >
void SendRecv
(const T * sbuf, int sc, int to,
 T * rbuf, int rc, int from, Comm comm)
{
    TaggedSendRecv (sbuf, sc, to, 0, rbuf, rc, from,
                    mpi::ANY_TAG, comm);
}

template void SendRecv
(const byte * sbuf, int sc, int to,
 byte * rbuf, int rc, int from, Comm comm);
template void SendRecv
(const int *sbuf, int sc, int to,
 int *rbuf, int rc, int from, Comm comm);
template void SendRecv
(const unsigned *sbuf, int sc, int to,
 unsigned *rbuf, int rc, int from, Comm comm);
template void SendRecv
(const long int *sbuf, int sc, int to,
 long int *rbuf, int rc, int from, Comm comm);
template void SendRecv
(const unsigned long *sbuf, int sc, int to,
 unsigned long *rbuf, int rc, int from, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void SendRecv
(const long long int *sbuf, int sc, int to,
 long long int *rbuf, int rc, int from, Comm comm);
template void SendRecv
(const unsigned long long *sbuf, int sc, int to,
 unsigned long long *rbuf, int rc, int from,
 Comm comm);
#endif
template void SendRecv
(const float *sbuf, int sc, int to,
 float *rbuf, int rc, int from, Comm comm);
template void SendRecv
(const double *sbuf, int sc, int to,
 double *rbuf, int rc, int from, Comm comm);
template void SendRecv
(const Complex < float >*sbuf, int sc, int to,
 Complex < float >*rbuf, int rc, int from, Comm comm);
template void SendRecv
(const Complex < double >*sbuf, int sc, int to,
 Complex < double >*rbuf, int rc, int from,
 Comm comm);

template < typename T >
T TaggedSendRecv (T sb, int to, int stag, int from,
                  int rtag, Comm comm)
{
    T rb;

    TaggedSendRecv (&sb, 1, to, stag, &rb, 1, from, rtag,
                    comm);
    return rb;
}

template byte TaggedSendRecv
(byte sb, int to, int stag, int from, int rtag,
 Comm comm);
template int TaggedSendRecv (int sb, int to, int stag,
                             int from, int rtag,
                             Comm comm);
template unsigned TaggedSendRecv (unsigned sb, int to,
                                  int stag, int from,
                                  int rtag, Comm comm);
template long int TaggedSendRecv (long int sb, int to,
                                  int stag, int from,
                                  int rtag, Comm comm);
template unsigned long TaggedSendRecv (unsigned long sb,
                                       int to, int stag,
                                       int from, int rtag,
                                       Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int TaggedSendRecv
(long long int sb, int to, int stag, int from,
 int rtag, Comm comm);
template unsigned long long TaggedSendRecv (unsigned long
        long sb,
        int to,
        int stag,
        int from,
        int rtag,
        Comm comm);
#endif
template float TaggedSendRecv
(float sb, int to, int stag, int from, int rtag,
 Comm comm);
template double TaggedSendRecv (double sb, int to,
                                int stag, int from,
                                int rtag, Comm comm);
template Complex < float >TaggedSendRecv (Complex <
        float >sb,
        int to,
        int stag,
        int from,
        int rtag,
        Comm comm);
template Complex < double >TaggedSendRecv (Complex <
        double >sb,
        int to,
        int stag,
        int from,
        int rtag,
        Comm comm);

template < typename T >
T SendRecv (T sb, int to, int from, Comm comm)
{
    return TaggedSendRecv (sb, to, 0, from, mpi::ANY_TAG,
                           comm);
}

template byte SendRecv (byte sb, int to, int from,
                        Comm comm);
template int SendRecv (int sb, int to, int from,
                       Comm comm);
template unsigned SendRecv (unsigned sb, int to, int from,
                            Comm comm);
template long int SendRecv (long int sb, int to, int from,
                            Comm comm);
template unsigned long SendRecv (unsigned long sb, int to,
                                 int from, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int SendRecv (long long int sb, int to,
                                 int from, Comm comm);
template unsigned long long SendRecv (unsigned long long
                                      sb, int to,
                                      int from,
                                      Comm comm);
#endif
template float SendRecv (float sb, int to, int from,
                         Comm comm);
template double SendRecv (double sb, int to, int from,
                          Comm comm);
template Complex < float >SendRecv (Complex < float >sb,
                                    int to, int from,
                                    Comm comm);
template Complex < double >SendRecv (Complex < double >sb,
                                     int to, int from,
                                     Comm comm);

template < typename R >
void TaggedSendRecv
(R * buf, int count, int to, int stag, int from,
 int rtag, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::SendRecv"))
    Status status;

    SafeMpi (MPI_Sendrecv_replace
             (buf, count, TypeMap < R > (), to, stag,
              from, rtag, comm.comm, &status));
}

template < typename R >
void TaggedSendRecv
(Complex < R > *buf, int count, int to, int stag,
 int from, int rtag, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::SendRecv"))
    Status status;

#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Sendrecv_replace
     (buf, 2 * count, TypeMap < R > (), to, stag,
      from, rtag, comm.comm, &status));
#else
    SafeMpi
    (MPI_Sendrecv_replace
     (buf, count, TypeMap < Complex < R >> (),
      to, stag, from, rtag, comm.comm, &status));
#endif
}

template void TaggedSendRecv
(byte * buf, int count, int to, int stag, int from,
 int rtag, Comm comm);
template void TaggedSendRecv (int *buf, int count, int to,
                              int stag, int from,
                              int rtag, Comm comm);
template void TaggedSendRecv (unsigned *buf, int count,
                              int to, int stag, int from,
                              int rtag, Comm comm);
template void TaggedSendRecv (long int *buf, int count,
                              int to, int stag, int from,
                              int rtag, Comm comm);
template void TaggedSendRecv (unsigned long *buf,
                              int count, int to, int stag,
                              int from, int rtag,
                              Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedSendRecv
(long long int *buf, int count, int to, int stag,
 int from, int rtag, Comm comm);
template void TaggedSendRecv (unsigned long long *buf,
                              int count, int to, int stag,
                              int from, int rtag,
                              Comm comm);
#endif
template void TaggedSendRecv
(float *buf, int count, int to, int stag, int from,
 int rtag, Comm comm);
template void TaggedSendRecv (double *buf, int count,
                              int to, int stag, int from,
                              int rtag, Comm comm);
template void TaggedSendRecv (Complex < float >*buf,
                              int count, int to, int stag,
                              int from, int rtag,
                              Comm comm);
template void TaggedSendRecv (Complex < double >*buf,
                              int count, int to, int stag,
                              int from, int rtag,
                              Comm comm);

template < typename T >
void SendRecv (T * buf, int count, int to, int from,
               Comm comm)
{
    TaggedSendRecv (buf, count, to, 0, from, mpi::ANY_TAG,
                    comm);
}

template void SendRecv
(byte * buf, int count, int to, int from, Comm comm);
template void SendRecv
(int *buf, int count, int to, int from, Comm comm);
template void SendRecv
(unsigned *buf, int count, int to, int from,
 Comm comm);
template void SendRecv (long int *buf, int count, int to,
                        int from, Comm comm);
template void SendRecv (unsigned long *buf, int count,
                        int to, int from, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void SendRecv
(long long int *buf, int count, int to, int from,
 Comm comm);
template void SendRecv (unsigned long long *buf,
                        int count, int to, int from,
                        Comm comm);
#endif
template void SendRecv
(float *buf, int count, int to, int from, Comm comm);
template void SendRecv
(double *buf, int count, int to, int from, Comm comm);
template void SendRecv
(Complex < float >*buf, int count, int to, int from,
 Comm comm);
template void SendRecv (Complex < double >*buf, int count,
                        int to, int from, Comm comm);

template < typename R >
void Broadcast (R * buf, int count, int root,
                Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Broadcast"))
    SafeMpi (MPI_Bcast
             (buf, count, TypeMap < R > (), root,
              comm.comm));
}

template < typename R >
void Broadcast (Complex < R > *buf, int count,
                int root, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Broadcast"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi (MPI_Bcast
             (buf, 2 * count, TypeMap < R > (), root,
              comm.comm));
#else
    SafeMpi (MPI_Bcast
             (buf, count, TypeMap < Complex < R >> (),
              root, comm.comm));
#endif
}

template void Broadcast (byte * buf, int count, int root,
                         Comm comm);
template void Broadcast (int *buf, int count, int root,
                         Comm comm);
template void Broadcast (unsigned *buf, int count,
                         int root, Comm comm);
template void Broadcast (long int *buf, int count,
                         int root, Comm comm);
template void Broadcast (unsigned long *buf, int count,
                         int root, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Broadcast (long long int *buf, int count,
                         int root, Comm comm);
template void Broadcast (unsigned long long *buf,
                         int count, int root, Comm comm);
#endif
template void Broadcast (float *buf, int count, int root,
                         Comm comm);
template void Broadcast (double *buf, int count, int root,
                         Comm comm);
template void Broadcast (Complex < float >*buf, int count,
                         int root, Comm comm);
template void Broadcast (Complex < double >*buf,
                         int count, int root, Comm comm);

template < typename T > void Broadcast (T & b, int root,
                                        Comm comm)
{
    Broadcast (&b, 1, root, comm);
}

template void Broadcast (byte & b, int root, Comm comm);
template void Broadcast (int &b, int root, Comm comm);
template void Broadcast (unsigned &b, int root,
                         Comm comm);
template void Broadcast (long int &b, int root,
                         Comm comm);
template void Broadcast (unsigned long &b, int root,
                         Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Broadcast (long long int &b, int root,
                         Comm comm);
template void Broadcast (unsigned long long &b, int root,
                         Comm comm);
#endif
template void Broadcast (float &b, int root, Comm comm);
template void Broadcast (double &b, int root, Comm comm);
template void Broadcast (Complex < float >&b, int root,
                         Comm comm);
template void Broadcast (Complex < double >&b, int root,
                         Comm comm);

#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
template < typename R >
void IBroadcast (R * buf, int count, int root,
                 Comm comm, Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::IBroadcast"))
    SafeMpi
    (MPI_Ibcast
     (buf, count, TypeMap < R > (), root, comm.comm,
      &request));
}

template < typename R >
void IBroadcast
(Complex < R > *buf, int count, int root, Comm comm,
 Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::IBroadcast"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Ibcast
     (buf, 2 * count, TypeMap < R > (), root,
      comm.comm, &request));
#else
    SafeMpi
    (MPI_Ibcast
     (buf, count, TypeMap < Complex < R >> (), root,
      comm.comm, &request));
#endif
}

template void IBroadcast (byte * buf, int count, int root,
                          Comm comm, Request & request);
template void IBroadcast (int *buf, int count, int root,
                          Comm comm, Request & request);
template void IBroadcast (unsigned *buf, int count,
                          int root, Comm comm,
                          Request & request);
template void IBroadcast (long int *buf, int count,
                          int root, Comm comm,
                          Request & request);
template void IBroadcast (unsigned long *buf, int count,
                          int root, Comm comm,
                          Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void IBroadcast (long long int *buf, int count,
                          int root, Comm comm,
                          Request & request);
template void IBroadcast (unsigned long long *buf,
                          int count, int root, Comm comm,
                          Request & request);
#endif
template void IBroadcast (float *buf, int count, int root,
                          Comm comm, Request & request);
template void IBroadcast (double *buf, int count,
                          int root, Comm comm,
                          Request & request);
template void IBroadcast (Complex < float >*buf,
                          int count, int root, Comm comm,
                          Request & request);
template void IBroadcast (Complex < double >*buf,
                          int count, int root, Comm comm,
                          Request & request);

template < typename T >
void IBroadcast (T & b, int root, Comm comm,
                 Request & request)
{
    IBroadcast (&b, 1, root, comm, request);
}

template void IBroadcast (byte & b, int root, Comm comm,
                          Request & request);
template void IBroadcast (int &b, int root, Comm comm,
                          Request & request);
template void IBroadcast (unsigned &b, int root,
                          Comm comm, Request & request);
template void IBroadcast (long int &b, int root,
                          Comm comm, Request & request);
template void IBroadcast (unsigned long &b, int root,
                          Comm comm, Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void IBroadcast (long long int &b, int root,
                          Comm comm, Request & request);
template void IBroadcast (unsigned long long &b, int root,
                          Comm comm, Request & request);
#endif
template void IBroadcast (float &b, int root, Comm comm,
                          Request & request);
template void IBroadcast (double &b, int root, Comm comm,
                          Request & request);
template void IBroadcast (Complex < float >&b, int root,
                          Comm comm, Request & request);
template void IBroadcast (Complex < double >&b, int root,
                          Comm comm, Request & request);
#endif // ifdef EL_HAVE_NONBLOCKING_COLLECTIVES

template < typename R >
void Gather
(const R * sbuf, int sc, R * rbuf, int rc, int root,
 Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Gather"))
    SafeMpi
    (MPI_Gather
     (const_cast < R * >(sbuf), sc, TypeMap < R > (),
      rbuf, rc, TypeMap < R > (), root, comm.comm));
}

template < typename R >
void Gather
(const Complex < R > *sbuf, int sc,
 Complex < R > *rbuf, int rc, int root, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Gather"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Gather
     (const_cast < Complex < R > *>(sbuf), 2 * sc,
      TypeMap < R > (), rbuf, 2 * rc,
      TypeMap < R > (), root, comm.comm));
#else
    SafeMpi
    (MPI_Gather
     (const_cast < Complex < R > *>(sbuf), sc,
      TypeMap < Complex < R >> (), rbuf, rc,
      TypeMap < Complex < R >> (), root, comm.comm));
#endif
}

template void Gather (const byte * sbuf, int sc,
                      byte * rbuf, int rc, int root,
                      Comm comm);
template void Gather (const int *sbuf, int sc, int *rbuf,
                      int rc, int root, Comm comm);
template void Gather (const unsigned *sbuf, int sc,
                      unsigned *rbuf, int rc, int root,
                      Comm comm);
template void Gather (const long int *sbuf, int sc,
                      long int *rbuf, int rc, int root,
                      Comm comm);
template void Gather (const unsigned long *sbuf, int sc,
                      unsigned long *rbuf, int rc,
                      int root, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Gather (const long long int *sbuf, int sc,
                      long long int *rbuf, int rc,
                      int root, Comm comm);
template void Gather (const unsigned long long *sbuf,
                      int sc, unsigned long long *rbuf,
                      int rc, int root, Comm comm);
#endif
template void Gather (const float *sbuf, int sc,
                      float *rbuf, int rc, int root,
                      Comm comm);
template void Gather (const double *sbuf, int sc,
                      double *rbuf, int rc, int root,
                      Comm comm);
template void Gather (const Complex < float >*sbuf,
                      int sc, Complex < float >*rbuf,
                      int rc, int root, Comm comm);
template void Gather (const Complex < double >*sbuf,
                      int sc, Complex < double >*rbuf,
                      int rc, int root, Comm comm);

#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
template < typename R >
void IGather
(const R * sbuf, int sc,
 R * rbuf, int rc, int root, Comm comm,
 Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::IGather"))
    SafeMpi
    (MPI_Igather
     (const_cast < R * >(sbuf), sc, TypeMap < R > (),
      rbuf, rc, TypeMap < R > (), root, comm.comm,
      &request));
}

template < typename R >
void IGather
(const Complex < R > *sbuf, int sc,
 Complex < R > *rbuf, int rc, int root, Comm comm,
 Request & request)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::IGather"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Igather
     (const_cast < Complex < R > *>(sbuf), 2 * sc,
      TypeMap < R > (), rbuf, 2 * rc,
      TypeMap < R > (), root, comm.comm, &request));
#else
    SafeMpi
    (MPI_Igather
     (const_cast < Complex < R > *>(sbuf), sc,
      TypeMap < Complex < R >> (), rbuf, rc,
      TypeMap < Complex < R >> (), root, comm.comm,
      &request));
#endif
}

template void IGather
(const byte * sbuf, int sc,
 byte * rbuf, int rc, int root, Comm comm,
 Request & request);
template void IGather (const int *sbuf, int sc, int *rbuf,
                       int rc, int root, Comm comm,
                       Request & request);
template void IGather (const unsigned *sbuf, int sc,
                       unsigned *rbuf, int rc, int root,
                       Comm comm, Request & request);
template void IGather (const long int *sbuf, int sc,
                       long int *rbuf, int rc, int root,
                       Comm comm, Request & request);
template void IGather (const unsigned long *sbuf, int sc,
                       unsigned long *rbuf, int rc,
                       int root, Comm comm,
                       Request & request);
#ifdef EL_HAVE_MPI_LONG_LONG
template void IGather
(const long long int *sbuf, int sc,
 long long int *rbuf, int rc, int root, Comm comm,
 Request & request);
template void IGather (const unsigned long long *sbuf,
                       int sc, unsigned long long *rbuf,
                       int rc, int root, Comm comm,
                       Request & request);
#endif
template void IGather
(const float *sbuf, int sc,
 float *rbuf, int rc, int root, Comm comm,
 Request & request);
template void IGather (const double *sbuf, int sc,
                       double *rbuf, int rc, int root,
                       Comm comm, Request & request);
template void IGather (const Complex < float >*sbuf,
                       int sc, Complex < float >*rbuf,
                       int rc, int root, Comm comm,
                       Request & request);
template void IGather (const Complex < double >*sbuf,
                       int sc, Complex < double >*rbuf,
                       int rc, int root, Comm comm,
                       Request & request);
#endif // ifdef EL_HAVE_NONBLOCKING_COLLECTIVES

template < typename R >
void Gather
(const R * sbuf, int sc,
 R * rbuf, const int *rcs, const int *rds, int root,
 Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Gather"))
    SafeMpi
    (MPI_Gatherv
     (const_cast < R * >(sbuf),
      sc,
      TypeMap < R > (),
      rbuf,
      const_cast < int *>(rcs),
      const_cast < int *>(rds), TypeMap < R > (),
      root, comm.comm));
}

template < typename R >
void Gather
(const Complex < R > *sbuf, int sc,
 Complex < R > *rbuf, const int *rcs, const int *rds,
 int root, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Gather"))
#ifdef EL_AVOID_COMPLEX_MPI
    const int commRank = Rank( comm );
    const int commSize = Size( comm );
    vector<int> rcsDouble, rdsDouble;
    if( commRank == root )
    {
        rcsDouble.resize (commSize);
        rdsDouble.resize (commSize);
        for (int i = 0; i < commSize; ++i)
        {
            rcsDouble[i] = 2 * rcs[i];
            rdsDouble[i] = 2 * rds[i];
        }
    }
    SafeMpi
    (MPI_Gatherv
     (const_cast < Complex < R > *>(sbuf), 2 * sc,
      TypeMap < R > (), rbuf, rcsDouble.data (),
      rdsDouble.data (), TypeMap < R > (), root,
      comm.comm));
#else
    SafeMpi
    (MPI_Gatherv
     (const_cast < Complex < R > *>(sbuf),
      sc,
      TypeMap < Complex < R >> (),
      rbuf,
      const_cast < int *>(rcs),
      const_cast < int *>(rds),
      TypeMap < Complex < R >> (), root,
      comm.comm));
#endif
}

template void Gather
(const byte * sbuf, int sc,
 byte * rbuf, const int *rcs, const int *rds,
 int root, Comm comm);
template void Gather (const int *sbuf, int sc, int *rbuf,
                      const int *rcs, const int *rds,
                      int root, Comm comm);
template void Gather (const unsigned *sbuf, int sc,
                      unsigned *rbuf, const int *rcs,
                      const int *rds, int root,
                      Comm comm);
template void Gather (const long int *sbuf, int sc,
                      long int *rbuf, const int *rcs,
                      const int *rds, int root,
                      Comm comm);
template void Gather (const unsigned long *sbuf, int sc,
                      unsigned long *rbuf, const int *rcs,
                      const int *rds, int root,
                      Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Gather
(const long long int *sbuf, int sc,
 long long int *rbuf, const int *rcs, const int *rds,
 int root, Comm comm);
template void Gather (const unsigned long long *sbuf,
                      int sc, unsigned long long *rbuf,
                      const int *rcs, const int *rds,
                      int root, Comm comm);
#endif
template void Gather
(const float *sbuf, int sc,
 float *rbuf, const int *rcs, const int *rds,
 int root, Comm comm);
template void Gather (const double *sbuf, int sc,
                      double *rbuf, const int *rcs,
                      const int *rds, int root,
                      Comm comm);
template void Gather (const Complex < float >*sbuf,
                      int sc, Complex < float >*rbuf,
                      const int *rcs, const int *rds,
                      int root, Comm comm);
template void Gather (const Complex < double >*sbuf,
                      int sc, Complex < double >*rbuf,
                      const int *rcs, const int *rds,
                      int root, Comm comm);

template < typename R >
void AllGather
(const R * sbuf, int sc, R * rbuf, int rc, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::AllGather"))
#ifdef EL_USE_BYTE_ALLGATHERS
    SafeMpi
    (MPI_Allgather
     ((UCP) const_cast < R * >(sbuf), sizeof (R) * sc,
      MPI_UNSIGNED_CHAR, (UCP) rbuf, sizeof (R) * rc,
      MPI_UNSIGNED_CHAR, comm.comm));
#else
    SafeMpi
    (MPI_Allgather
     (const_cast < R * >(sbuf), sc, TypeMap < R > (),
      rbuf, rc, TypeMap < R > (), comm.comm));
#endif
}

template < typename R >
void AllGather
(const Complex < R > *sbuf, int sc,
 Complex < R > *rbuf, int rc, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::AllGather"))
#ifdef EL_USE_BYTE_ALLGATHERS
    SafeMpi
    (MPI_Allgather
     ((UCP) const_cast < Complex < R > *>(sbuf),
      2 * sizeof (R) * sc, MPI_UNSIGNED_CHAR,
      (UCP) rbuf, 2 * sizeof (R) * rc,
      MPI_UNSIGNED_CHAR, comm.comm));
#else
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Allgather
     (const_cast < Complex < R > *>(sbuf), 2 * sc,
      TypeMap < R > (), rbuf, 2 * rc,
      TypeMap < R > (), comm.comm));
#else
    SafeMpi
    (MPI_Allgather
     (const_cast < Complex < R > *>(sbuf), sc,
      TypeMap < Complex < R >> (), rbuf, rc,
      TypeMap < Complex < R >> (), comm.comm));
#endif
#endif
}

template void AllGather (const byte * sbuf, int sc,
                         byte * rbuf, int rc, Comm comm);
template void AllGather (const int *sbuf, int sc,
                         int *rbuf, int rc, Comm comm);
template void AllGather (const unsigned *sbuf, int sc,
                         unsigned *rbuf, int rc,
                         Comm comm);
template void AllGather (const long int *sbuf, int sc,
                         long int *rbuf, int rc,
                         Comm comm);
template void AllGather (const unsigned long *sbuf,
                         int sc, unsigned long *rbuf,
                         int rc, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllGather (const long long int *sbuf,
                         int sc, long long int *rbuf,
                         int rc, Comm comm);
template void AllGather (const unsigned long long *sbuf,
                         int sc, unsigned long long *rbuf,
                         int rc, Comm comm);
#endif
template void AllGather (const float *sbuf, int sc,
                         float *rbuf, int rc, Comm comm);
template void AllGather (const double *sbuf, int sc,
                         double *rbuf, int rc, Comm comm);
template void AllGather (const Complex < float >*sbuf,
                         int sc, Complex < float >*rbuf,
                         int rc, Comm comm);
template void AllGather (const Complex < double >*sbuf,
                         int sc, Complex < double >*rbuf,
                         int rc, Comm comm);

template < typename R >
void AllGather
(const R * sbuf, int sc,
 R * rbuf, const int *rcs, const int *rds, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::AllGather"))
#ifdef EL_USE_BYTE_ALLGATHERS
    const int commSize = Size( comm );
    vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = sizeof (R) * rcs[i];
        byteRds[i] = sizeof (R) * rds[i];
    }
    SafeMpi
    (MPI_Allgatherv
     ((UCP) const_cast < R * >(sbuf), sizeof (R) * sc,
      MPI_UNSIGNED_CHAR, (UCP) rbuf, byteRcs.data (),
      byteRds.data (), MPI_UNSIGNED_CHAR, comm.comm));
#else
    SafeMpi
    (MPI_Allgatherv
     (const_cast < R * >(sbuf),
      sc,
      TypeMap < R > (),
      rbuf,
      const_cast < int *>(rcs),
      const_cast < int *>(rds), TypeMap < R > (),
      comm.comm));
#endif
}

template < typename R >
void AllGather
(const Complex < R > *sbuf, int sc,
 Complex < R > *rbuf, const int *rcs, const int *rds,
 Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::AllGather"))
#ifdef EL_USE_BYTE_ALLGATHERS
    const int commSize = Size( comm );
    vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = 2 * sizeof (R) * rcs[i];
        byteRds[i] = 2 * sizeof (R) * rds[i];
    }
    SafeMpi
    (MPI_Allgatherv
     ((UCP) const_cast < Complex < R > *>(sbuf),
      2 * sizeof (R) * sc, MPI_UNSIGNED_CHAR,
      (UCP) rbuf, byteRcs.data (), byteRds.data (),
      MPI_UNSIGNED_CHAR, comm.comm));
#else
 #ifdef EL_AVOID_COMPLEX_MPI
    const int commSize = Size( comm );
    vector<int> realRcs( commSize ), realRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        realRcs[i] = 2 * rcs[i];
        realRds[i] = 2 * rds[i];
    }
    SafeMpi
    (MPI_Allgatherv
     (const_cast < Complex < R > *>(sbuf), 2 * sc,
      TypeMap < R > (), rbuf, realRcs.data (),
      realRds.data (), TypeMap < R > (), comm.comm));
#else
    SafeMpi
    (MPI_Allgatherv
     (const_cast < Complex < R > *>(sbuf),
      sc,
      TypeMap < Complex < R >> (),
      rbuf,
      const_cast < int *>(rcs),
      const_cast < int *>(rds),
      TypeMap < Complex < R >> (), comm.comm));
#endif
#endif
}

template void AllGather
(const byte * sbuf, int sc,
 byte * rbuf, const int *rcs, const int *rds,
 Comm comm);
template void AllGather (const int *sbuf, int sc,
                         int *rbuf, const int *rcs,
                         const int *rds, Comm comm);
template void AllGather (const unsigned *sbuf, int sc,
                         unsigned *rbuf, const int *rcs,
                         const int *rds, Comm comm);
template void AllGather (const long int *sbuf, int sc,
                         long int *rbuf, const int *rcs,
                         const int *rds, Comm comm);
template void AllGather (const unsigned long *sbuf,
                         int sc, unsigned long *rbuf,
                         const int *rcs, const int *rds,
                         Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllGather
(const long long int *sbuf, int sc,
 long long int *rbuf, const int *rcs, const int *rds,
 Comm comm);
template void AllGather (const unsigned long long *sbuf,
                         int sc, unsigned long long *rbuf,
                         const int *rcs, const int *rds,
                         Comm comm);
#endif
template void AllGather
(const float *sbuf, int sc,
 float *rbuf, const int *rcs, const int *rds,
 Comm comm);
template void AllGather (const double *sbuf, int sc,
                         double *rbuf, const int *rcs,
                         const int *rds, Comm comm);
template void AllGather (const Complex < float >*sbuf,
                         int sc, Complex < float >*rbuf,
                         const int *rcs, const int *rds,
                         Comm comm);
template void AllGather (const Complex < double >*sbuf,
                         int sc, Complex < double >*rbuf,
                         const int *rcs, const int *rds,
                         Comm comm);

template < typename R >
void Scatter
(const R * sbuf, int sc, R * rbuf, int rc, int root,
 Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Scatter"))
    SafeMpi
    (MPI_Scatter
     (const_cast < R * >(sbuf), sc, TypeMap < R > (),
      rbuf, rc, TypeMap < R > (), root, comm.comm));
}

template < typename R >
void Scatter
(const Complex < R > *sbuf, int sc,
 Complex < R > *rbuf, int rc, int root, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Scatter"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Scatter
     (const_cast < Complex < R > *>(sbuf), 2 * sc,
      TypeMap < R > (), rbuf, 2 * rc,
      TypeMap < R > (), root, comm.comm));
#else
    SafeMpi
    (MPI_Scatter
     (const_cast < Complex < R > *>(sbuf), sc,
      TypeMap < Complex < R >> (), rbuf, rc,
      TypeMap < Complex < R >> (), root, comm.comm));
#endif
}

template void Scatter
(const byte * sbuf, int sc,
 byte * rbuf, int rc, int root, Comm comm);
template void Scatter
(const int *sbuf, int sc, int *rbuf, int rc, int root,
 Comm comm);
template void Scatter (const unsigned *sbuf, int sc,
                       unsigned *rbuf, int rc, int root,
                       Comm comm);
template void Scatter (const long int *sbuf, int sc,
                       long int *rbuf, int rc, int root,
                       Comm comm);
template void Scatter (const unsigned long *sbuf, int sc,
                       unsigned long *rbuf, int rc,
                       int root, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Scatter
(const long long int *sbuf, int sc,
 long long int *rbuf, int rc, int root, Comm comm);
template void Scatter
(const unsigned long long *sbuf, int sc,
 unsigned long long *rbuf, int rc, int root,
 Comm comm);
#endif
template void Scatter
(const float *sbuf, int sc,
 float *rbuf, int rc, int root, Comm comm);
template void Scatter
(const double *sbuf, int sc,
 double *rbuf, int rc, int root, Comm comm);
template void Scatter
(const Complex < float >*sbuf, int sc,
 Complex < float >*rbuf, int rc, int root, Comm comm);
template void Scatter
(const Complex < double >*sbuf, int sc,
 Complex < double >*rbuf, int rc, int root,
 Comm comm);

template < typename R >
void Scatter (R * buf, int sc, int rc, int root,
              Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Scatter"))
    const int commRank = Rank (comm);

    if (commRank == root)
    {
#ifdef EL_HAVE_MPI_IN_PLACE
        SafeMpi
        (MPI_Scatter
         (buf, sc, TypeMap < R > (),
          MPI_IN_PLACE, rc, TypeMap < R > (), root,
          comm.comm));
#else
        const int commSize = Size( comm );
        vector<R> sendBuf( sc*commSize );
        MemCopy( sendBuf.data(), buf, sc*commSize );
        SafeMpi
        (MPI_Scatter
         (sendBuf.data (), sc, TypeMap < R > (),
          buf, rc, TypeMap < R > (), root,
          comm.comm));
#endif
    }
    else
    {
        SafeMpi
        (MPI_Scatter
         (0, sc, TypeMap < R > (),
          buf, rc, TypeMap < R > (), root,
          comm.comm));
    }
}

template < typename R >
void Scatter (Complex < R > *buf, int sc, int rc,
              int root, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Scatter"))
    const int commRank = Rank (comm);

    if (commRank == root)
    {
#ifdef EL_AVOID_COMPLEX_MPI
#ifdef EL_HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Scatter
          ( buf,          2*sc, TypeMap<R>(), 
            MPI_IN_PLACE, 2*rc, TypeMap<R>(), root, comm.comm ) );
# else
        const int commSize = Size( comm );
        vector<Complex<R>> sendBuf( sc*commSize );
        MemCopy( sendBuf.data(), buf, sc*commSize );
        SafeMpi
        (MPI_Scatter
         (sendBuf.data (), 2 * sc, TypeMap < R > (),
          buf, 2 * rc, TypeMap < R > (), root,
          comm.comm));
#endif
#else
#ifdef EL_HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Scatter
          ( buf,          sc, TypeMap<Complex<R>>(), 
            MPI_IN_PLACE, rc, TypeMap<Complex<R>>(), root, comm.comm ) );
# else
        const int commSize = Size( comm );
        vector<Complex<R>> sendBuf( sc*commSize );
        MemCopy( sendBuf.data(), buf, sc*commSize );
        SafeMpi
        (MPI_Scatter
         (sendBuf.data (), sc,
          TypeMap < Complex < R >> (), buf, rc,
          TypeMap < Complex < R >> (), root,
          comm.comm));
#endif
#endif
    }
    else
    {
#ifdef EL_AVOID_COMPLEX_MPI
        SafeMpi
        (MPI_Scatter
         (0, 2 * sc, TypeMap < R > (),
          buf, 2 * rc, TypeMap < R > (), root,
          comm.comm));
#else
        SafeMpi
        (MPI_Scatter
         (0, sc, TypeMap < Complex < R >> (),
          buf, rc, TypeMap < Complex < R >> (),
          root, comm.comm));
#endif
    }
}

template void Scatter (byte * buf, int sc, int rc,
                       int root, Comm comm);
template void Scatter (int *buf, int sc, int rc, int root,
                       Comm comm);
template void Scatter (unsigned *buf, int sc, int rc,
                       int root, Comm comm);
template void Scatter (long int *buf, int sc, int rc,
                       int root, Comm comm);
template void Scatter (unsigned long *buf, int sc, int rc,
                       int root, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Scatter (long long int *buf, int sc, int rc,
                       int root, Comm comm);
template void Scatter (unsigned long long *buf, int sc,
                       int rc, int root, Comm comm);
#endif
template void Scatter (float *buf, int sc, int rc,
                       int root, Comm comm);
template void Scatter (double *buf, int sc, int rc,
                       int root, Comm comm);
template void Scatter (Complex < float >*buf, int sc,
                       int rc, int root, Comm comm);
template void Scatter (Complex < double >*buf, int sc,
                       int rc, int root, Comm comm);

template < typename R >
void AllToAll
(const R * sbuf, int sc, R * rbuf, int rc, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::AllToAll"))
    SafeMpi
    (MPI_Alltoall
     (const_cast < R * >(sbuf), sc, TypeMap < R > (),
      rbuf, rc, TypeMap < R > (), comm.comm));
}

template < typename R >
void AllToAll
(const Complex < R > *sbuf, int sc,
 Complex < R > *rbuf, int rc, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::AllToAll"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Alltoall
     (const_cast < Complex < R > *>(sbuf), 2 * sc,
      TypeMap < R > (), rbuf, 2 * rc,
      TypeMap < R > (), comm.comm));
#else
    SafeMpi
    (MPI_Alltoall
     (const_cast < Complex < R > *>(sbuf), sc,
      TypeMap < Complex < R >> (), rbuf, rc,
      TypeMap < Complex < R >> (), comm.comm));
#endif
}

template void AllToAll
(const byte * sbuf, int sc, byte * rbuf, int rc,
 Comm comm);
template void AllToAll (const int *sbuf, int sc,
                        int *rbuf, int rc, Comm comm);
template void AllToAll (const unsigned *sbuf, int sc,
                        unsigned *rbuf, int rc,
                        Comm comm);
template void AllToAll (const long int *sbuf, int sc,
                        long int *rbuf, int rc,
                        Comm comm);
template void AllToAll (const unsigned long *sbuf, int sc,
                        unsigned long *rbuf, int rc,
                        Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllToAll
(const long long int *sbuf, int sc,
 long long int *rbuf, int rc, Comm comm);
template void AllToAll
(const unsigned long long *sbuf, int sc,
 unsigned long long *rbuf, int rc, Comm comm);
#endif
template void AllToAll
(const float *sbuf, int sc, float *rbuf, int rc,
 Comm comm);
template void AllToAll (const double *sbuf, int sc,
                        double *rbuf, int rc, Comm comm);
template void AllToAll (const Complex < float >*sbuf,
                        int sc, Complex < float >*rbuf,
                        int rc, Comm comm);
template void AllToAll (const Complex < double >*sbuf,
                        int sc, Complex < double >*rbuf,
                        int rc, Comm comm);

template < typename R >
void AllToAll
(const R * sbuf, const int *scs, const int *sds,
 R * rbuf, const int *rcs, const int *rds, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::AllToAll"))
    SafeMpi
    (MPI_Alltoallv
     (const_cast < R * >(sbuf),
      const_cast < int *>(scs),
      const_cast < int *>(sds),
      TypeMap < R > (),
      rbuf,
      const_cast < int *>(rcs),
      const_cast < int *>(rds), TypeMap < R > (),
      comm.comm));
}

template < typename R >
void AllToAll
(const Complex < R > *sbuf, const int *scs,
 const int *sds, Complex < R > *rbuf, const int *rcs,
 const int *rds, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::AllToAll"))
#ifdef EL_AVOID_COMPLEX_MPI
    int p;
    MPI_Comm_size( comm.comm, &p );
    vector<int> scsDoubled(p);
    vector<int> sdsDoubled(p);
    vector<int> rcsDoubled(p);
    vector<int> rdsDoubled(p);
    for( int i=0; i<p; ++i )
        scsDoubled[i] = 2*scs[i];
    for( int i=0; i<p; ++i )
        sdsDoubled[i] = 2*sds[i];
    for( int i=0; i<p; ++i )
        rcsDoubled[i] = 2*rcs[i];
    for( int i=0; i<p; ++i )
        rdsDoubled[i] = 2*rds[i];
    SafeMpi
    ( MPI_Alltoallv
      ( const_cast<Complex<R>*>(sbuf),
              scsDoubled.data(), sdsDoubled.data(), TypeMap<R>(),
        rbuf, rcsDoubled.data(), rdsDoubled.data(), TypeMap<R>(), comm.comm ) );
#else
    SafeMpi
    ( MPI_Alltoallv
      ( const_cast<Complex<R>*>(sbuf), 
        const_cast<int*>(scs), 
        const_cast<int*>(sds), 
        TypeMap<Complex<R>>(),
        rbuf, 
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<Complex<R>>(),
        comm.comm ) );
#endif
}

template void AllToAll
(const byte * sbuf, const int *scs, const int *sds,
 byte * rbuf, const int *rcs, const int *rds,
 Comm comm);
template void AllToAll (const int *sbuf, const int *scs,
                        const int *sds, int *rbuf,
                        const int *rcs, const int *rds,
                        Comm comm);
template void AllToAll (const unsigned *sbuf,
                        const int *scs, const int *sds,
                        unsigned *rbuf, const int *rcs,
                        const int *rds, Comm comm);
template void AllToAll (const long int *sbuf,
                        const int *scs, const int *sds,
                        long int *rbuf, const int *rcs,
                        const int *rds, Comm comm);
template void AllToAll (const unsigned long *sbuf,
                        const int *scs, const int *sds,
                        unsigned long *rbuf,
                        const int *rcs, const int *rds,
                        Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllToAll
(const long long int *sbuf, const int *scs,
 const int *sds, long long int *rbuf, const int *rcs,
 const int *rds, Comm comm);
template void AllToAll (const unsigned long long *sbuf,
                        const int *scs, const int *sds,
                        unsigned long long *rbuf,
                        const int *rcs, const int *rds,
                        Comm comm);
#endif
template void AllToAll
(const float *sbuf, const int *scs, const int *sds,
 float *rbuf, const int *rcs, const int *rds,
 Comm comm);
template void AllToAll (const double *sbuf,
                        const int *scs, const int *sds,
                        double *rbuf, const int *rcs,
                        const int *rds, Comm comm);
template void AllToAll (const Complex < float >*sbuf,
                        const int *scs, const int *sds,
                        Complex < float >*rbuf,
                        const int *rcs, const int *rds,
                        Comm comm);
template void AllToAll (const Complex < double >*sbuf,
                        const int *scs, const int *sds,
                        Complex < double >*rbuf,
                        const int *rcs, const int *rds,
                        Comm comm);

template < typename T >
void Reduce
(const T * sbuf, T * rbuf, int count, Op op, int root,
 Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Reduce"))
    if (count != 0)
    {
        SafeMpi (MPI_Reduce
                 (const_cast < T * >(sbuf), rbuf, count,
                  TypeMap < T > (), op.op, root,
                  comm.comm));
    }
}

template < typename R >
void Reduce
(const Complex < R > *sbuf,
 Complex < R > *rbuf, int count, Op op, int root,
 Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Reduce"))
    if (count != 0)
    {
#ifdef EL_AVOID_COMPLEX_MPI
        if (op == SUM)
        {
            SafeMpi
            (MPI_Reduce
             (const_cast < Complex < R > *>(sbuf),
              rbuf, 2 * count, TypeMap < R > (),
              op.op, root, comm.comm));
        }
        else
        {
            SafeMpi
            (MPI_Reduce
             (const_cast < Complex < R > *>(sbuf),
              rbuf, count,
              TypeMap < Complex < R >> (), op.op,
              root, comm.comm));
        }
#else
        SafeMpi
        (MPI_Reduce
         (const_cast < Complex < R > *>(sbuf),
          rbuf, count, TypeMap < Complex < R >> (),
          op.op, root, comm.comm));
#endif
    }
}

template void Reduce (const byte * sbuf, byte * rbuf,
                      int count, Op op, int root,
                      Comm comm);
template void Reduce (const int *sbuf, int *rbuf,
                      int count, Op op, int root,
                      Comm comm);
template void Reduce (const unsigned *sbuf,
                      unsigned *rbuf, int count, Op op,
                      int root, Comm comm);
template void Reduce (const long int *sbuf,
                      long int *rbuf, int count, Op op,
                      int root, Comm comm);
template void Reduce (const unsigned long *sbuf,
                      unsigned long *rbuf, int count,
                      Op op, int root, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Reduce (const long long int *sbuf,
                      long long int *rbuf, int count,
                      Op op, int root, Comm comm);
template void Reduce (const unsigned long long *sbuf,
                      unsigned long long *rbuf, int count,
                      Op op, int root, Comm comm);
#endif
template void Reduce (const float *sbuf, float *rbuf,
                      int count, Op op, int root,
                      Comm comm);
template void Reduce (const double *sbuf, double *rbuf,
                      int count, Op op, int root,
                      Comm comm);
template void Reduce (const Complex < float >*sbuf,
                      Complex < float >*rbuf, int count,
                      Op op, int root, Comm comm);
template void Reduce (const Complex < double >*sbuf,
                      Complex < double >*rbuf, int count,
                      Op op, int root, Comm comm);
template void Reduce (const ValueInt < Int > *sbuf,
                      ValueInt < Int > *rbuf, int count,
                      Op op, int root, Comm comm);
template void Reduce (const ValueInt < float >*sbuf,
                      ValueInt < float >*rbuf, int count,
                      Op op, int root, Comm comm);
template void Reduce (const ValueInt < double >*sbuf,
                      ValueInt < double >*rbuf, int count,
                      Op op, int root, Comm comm);
template void Reduce (const ValueIntPair < Int > *sbuf,
                      ValueIntPair < Int > *rbuf,
                      int count, Op op, int root,
                      Comm comm);
template void Reduce (const ValueIntPair < float >*sbuf,
                      ValueIntPair < float >*rbuf,
                      int count, Op op, int root,
                      Comm comm);
template void Reduce (const ValueIntPair < double >*sbuf,
                      ValueIntPair < double >*rbuf,
                      int count, Op op, int root,
                      Comm comm);

template < typename T >
void Reduce (const T * sbuf, T * rbuf, int count,
             int root, Comm comm)
{
    Reduce (sbuf, rbuf, count, mpi::SUM, root, comm);
}

template void Reduce (const byte * sbuf, byte * rbuf,
                      int count, int root, Comm comm);
template void Reduce (const int *sbuf, int *rbuf,
                      int count, int root, Comm comm);
template void Reduce (const unsigned *sbuf,
                      unsigned *rbuf, int count, int root,
                      Comm comm);
template void Reduce (const long int *sbuf,
                      long int *rbuf, int count, int root,
                      Comm comm);
template void Reduce (const unsigned long *sbuf,
                      unsigned long *rbuf, int count,
                      int root, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Reduce (const long long int *sbuf,
                      long long int *rbuf, int count,
                      int root, Comm comm);
template void Reduce (const unsigned long long *sbuf,
                      unsigned long long *rbuf, int count,
                      int root, Comm comm);
#endif
template void Reduce (const float *sbuf, float *rbuf,
                      int count, int root, Comm comm);
template void Reduce (const double *sbuf, double *rbuf,
                      int count, int root, Comm comm);
template void Reduce (const Complex < float >*sbuf,
                      Complex < float >*rbuf, int count,
                      int root, Comm comm);
template void Reduce (const Complex < double >*sbuf,
                      Complex < double >*rbuf, int count,
                      int root, Comm comm);
template void Reduce (const ValueInt < Int > *sbuf,
                      ValueInt < Int > *rbuf, int count,
                      int root, Comm comm);
template void Reduce (const ValueInt < float >*sbuf,
                      ValueInt < float >*rbuf, int count,
                      int root, Comm comm);
template void Reduce (const ValueInt < double >*sbuf,
                      ValueInt < double >*rbuf, int count,
                      int root, Comm comm);
template void Reduce (const ValueIntPair < Int > *sbuf,
                      ValueIntPair < Int > *rbuf,
                      int count, int root, Comm comm);
template void Reduce (const ValueIntPair < float >*sbuf,
                      ValueIntPair < float >*rbuf,
                      int count, int root, Comm comm);
template void Reduce (const ValueIntPair < double >*sbuf,
                      ValueIntPair < double >*rbuf,
                      int count, int root, Comm comm);

template < typename T > T Reduce (T sb, Op op, int root,
                                  Comm comm)
{
    T rb;

    Reduce (&sb, &rb, 1, op, root, comm);
    return rb;
}

template byte Reduce (byte sb, Op op, int root,
                      Comm comm);
template int Reduce (int sb, Op op, int root, Comm comm);
template unsigned Reduce (unsigned sb, Op op, int root,
                          Comm comm);
template long int Reduce (long int sb, Op op, int root,
                          Comm comm);
template unsigned long Reduce (unsigned long sb, Op op,
                               int root, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int Reduce (long long int sb, Op op,
                               int root, Comm comm);
template unsigned long long Reduce (unsigned long long sb,
                                    Op op, int root,
                                    Comm comm);
#endif
template float Reduce (float sb, Op op, int root,
                       Comm comm);
template double Reduce (double sb, Op op, int root,
                        Comm comm);
template Complex < float >Reduce (Complex < float >sb,
                                  Op op, int root,
                                  Comm comm);
template Complex < double >Reduce (Complex < double >sb,
                                   Op op, int root,
                                   Comm comm);
template ValueInt < Int > Reduce (ValueInt < Int > sb,
                                  Op op, int root,
                                  Comm comm);
template ValueInt < float >Reduce (ValueInt < float >sb,
                                   Op op, int root,
                                   Comm comm);
template ValueInt < double >Reduce (ValueInt < double >sb,
                                    Op op, int root,
                                    Comm comm);
template ValueIntPair < Int > Reduce (ValueIntPair < Int >
                                      sb, Op op, int root,
                                      Comm comm);
template ValueIntPair < float >Reduce (ValueIntPair <
                                       float >sb, Op op,
                                       int root,
                                       Comm comm);
template ValueIntPair < double >Reduce (ValueIntPair <
                                        double >sb, Op op,
                                        int root,
                                        Comm comm);

template < typename T > T Reduce (T sb, int root,
                                  Comm comm)
{
    T rb;

    Reduce (&sb, &rb, 1, mpi::SUM, root, comm);
    return rb;
}

template byte Reduce (byte sb, int root, Comm comm);
template int Reduce (int sb, int root, Comm comm);
template unsigned Reduce (unsigned sb, int root,
                          Comm comm);
template long int Reduce (long int sb, int root,
                          Comm comm);
template unsigned long Reduce (unsigned long sb, int root,
                               Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int Reduce (long long int sb, int root,
                               Comm comm);
template unsigned long long Reduce (unsigned long long sb,
                                    int root, Comm comm);
#endif
template float Reduce (float sb, int root, Comm comm);
template double Reduce (double sb, int root, Comm comm);
template Complex < float >Reduce (Complex < float >sb,
                                  int root, Comm comm);
template Complex < double >Reduce (Complex < double >sb,
                                   int root, Comm comm);
template ValueInt < Int > Reduce (ValueInt < Int > sb,
                                  int root, Comm comm);
template ValueInt < float >Reduce (ValueInt < float >sb,
                                   int root, Comm comm);
template ValueInt < double >Reduce (ValueInt < double >sb,
                                    int root, Comm comm);
template ValueIntPair < Int > Reduce (ValueIntPair < Int >
                                      sb, int root,
                                      Comm comm);
template ValueIntPair < float >Reduce (ValueIntPair <
                                       float >sb,
                                       int root,
                                       Comm comm);
template ValueIntPair < double >Reduce (ValueIntPair <
                                        double >sb,
                                        int root,
                                        Comm comm);

template < typename T >
void Reduce (T * buf, int count, Op op, int root,
             Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Reduce"))
    if (count != 0)
    {
        const int commRank = Rank (comm);

        if (commRank == root)
        {
#ifdef EL_HAVE_MPI_IN_PLACE
            SafeMpi
            (MPI_Reduce
             (MPI_IN_PLACE, buf, count,
              TypeMap < T > (), op.op, root,
              comm.comm));
#else
            vector<T> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            (MPI_Reduce
             (sendBuf.data (), buf, count,
              TypeMap < T > (), op.op, root,
              comm.comm));
#endif
        }
        else
            SafeMpi
            (MPI_Reduce
             (buf, 0, count, TypeMap < T > (),
              op.op, root, comm.comm));
    }
}

template < typename R >
void Reduce (Complex < R > *buf, int count, Op op,
             int root, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::Reduce"))
    if (count != 0)
    {
        const int commRank = Rank (comm);

#ifdef EL_AVOID_COMPLEX_MPI
        if (op == SUM)
        {
            if (commRank == root)
            {
# ifdef EL_HAVE_MPI_IN_PLACE
                SafeMpi
                ( MPI_Reduce
                  ( MPI_IN_PLACE, buf, 2*count, TypeMap<R>(), op.op, 
                    root, comm.comm ) );
# else
                std::vector < Complex <
                R >> sendBuf (count);
                MemCopy (sendBuf.data (), buf,
                         count);
                SafeMpi (MPI_Reduce
                         (sendBuf.data (), buf,
                          2 * count,
                          TypeMap < R > (), op.op,
                          root, comm.comm));
#endif
            }
            else
                SafeMpi
                (MPI_Reduce
                 (buf, 0, 2 * count,
                  TypeMap < R > (), op.op, root,
                  comm.comm));
        }
        else
        {
            if (commRank == root)
            {
#ifdef EL_HAVE_MPI_IN_PLACE
                SafeMpi
                (MPI_Reduce
                 (MPI_IN_PLACE, buf, count,
                  TypeMap < Complex < R >> (),
                  op.op, root, comm.comm));
#else
                std::vector < Complex <
                R >> sendBuf (count);
                MemCopy (sendBuf.data (), buf,
                         count);
                SafeMpi (MPI_Reduce
                         (sendBuf.data (), buf,
                          count,
                          TypeMap < Complex <
                          R >> (), op.op, root,
                          comm.comm));
#endif
            }
            else
                SafeMpi
                (MPI_Reduce
                 (buf, 0, count,
                  TypeMap < Complex < R >> (),
                  op.op, root, comm.comm));
        }
#else
        if (commRank == root)
        {
# ifdef EL_HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Reduce
              ( MPI_IN_PLACE, buf, count, TypeMap<Complex<R>>(), op.op, 
                root, comm.comm ) );
#else
        std::vector < Complex <
        R >> sendBuf (count);
        MemCopy (sendBuf.data (), buf, count);
        SafeMpi
        (MPI_Reduce
         (sendBuf.data (), buf, count,
          TypeMap < Complex < R >> (), op.op,
          root, comm.comm));
#endif
        }
        else
            SafeMpi
            (MPI_Reduce
             (buf, 0, count,
              TypeMap < Complex < R >> (), op.op,
              root, comm.comm));
#endif
    }
}

template void Reduce (byte * buf, int count, Op op,
                      int root, Comm comm);
template void Reduce (int *buf, int count, Op op,
                      int root, Comm comm);
template void Reduce (unsigned *buf, int count, Op op,
                      int root, Comm comm);
template void Reduce (long int *buf, int count, Op op,
                      int root, Comm comm);
template void Reduce (unsigned long *buf, int count,
                      Op op, int root, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Reduce (long long int *buf, int count,
                      Op op, int root, Comm comm);
template void Reduce (unsigned long long *buf, int count,
                      Op op, int root, Comm comm);
#endif
template void Reduce (float *buf, int count, Op op,
                      int root, Comm comm);
template void Reduce (double *buf, int count, Op op,
                      int root, Comm comm);
template void Reduce (Complex < float >*buf, int count,
                      Op op, int root, Comm comm);
template void Reduce (Complex < double >*buf, int count,
                      Op op, int root, Comm comm);
template void Reduce (ValueInt < Int > *buf, int count,
                      Op op, int root, Comm comm);
template void Reduce (ValueInt < float >*buf, int count,
                      Op op, int root, Comm comm);
template void Reduce (ValueInt < double >*buf, int count,
                      Op op, int root, Comm comm);
template void Reduce (ValueIntPair < Int > *buf,
                      int count, Op op, int root,
                      Comm comm);
template void Reduce (ValueIntPair < float >*buf,
                      int count, Op op, int root,
                      Comm comm);
template void Reduce (ValueIntPair < double >*buf,
                      int count, Op op, int root,
                      Comm comm);

template < typename T >
void Reduce (T * buf, int count, int root, Comm comm)
{
    Reduce (buf, count, mpi::SUM, root, comm);
}

template void Reduce (byte * buf, int count, int root,
                      Comm comm);
template void Reduce (int *buf, int count, int root,
                      Comm comm);
template void Reduce (unsigned *buf, int count, int root,
                      Comm comm);
template void Reduce (long int *buf, int count, int root,
                      Comm comm);
template void Reduce (unsigned long *buf, int count,
                      int root, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void Reduce (long long int *buf, int count,
                      int root, Comm comm);
template void Reduce (unsigned long long *buf, int count,
                      int root, Comm comm);
#endif
template void Reduce (float *buf, int count, int root,
                      Comm comm);
template void Reduce (double *buf, int count, int root,
                      Comm comm);
template void Reduce (Complex < float >*buf, int count,
                      int root, Comm comm);
template void Reduce (Complex < double >*buf, int count,
                      int root, Comm comm);
template void Reduce (ValueInt < Int > *buf, int count,
                      int root, Comm comm);
template void Reduce (ValueInt < float >*buf, int count,
                      int root, Comm comm);
template void Reduce (ValueInt < double >*buf, int count,
                      int root, Comm comm);
template void Reduce (ValueIntPair < Int > *buf,
                      int count, int root, Comm comm);
template void Reduce (ValueIntPair < float >*buf,
                      int count, int root, Comm comm);
template void Reduce (ValueIntPair < double >*buf,
                      int count, int root, Comm comm);

template < typename T >
void AllReduce (const T * sbuf, T * rbuf, int count,
                Op op, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::AllReduce"))
    if (count != 0)
    {
        SafeMpi (MPI_Allreduce
                 (const_cast < T * >(sbuf), rbuf, count,
                  TypeMap < T > (), op.op, comm.comm));
    }
}

template < typename R >
void AllReduce
(const Complex < R > *sbuf, Complex < R > *rbuf,
 int count, Op op, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::AllReduce"))
    if (count != 0)
    {
#ifdef EL_AVOID_COMPLEX_MPI
        if (op == SUM)
        {
            SafeMpi
            (MPI_Allreduce
             (const_cast < Complex < R > *>(sbuf),
              rbuf, 2 * count, TypeMap < R > (),
              op.op, comm.comm));
        }
        else
        {
            SafeMpi
            (MPI_Allreduce
             (const_cast < Complex < R > *>(sbuf),
              rbuf, count,
              TypeMap < Complex < R >> (), op.op,
              comm.comm));
        }
#else
        SafeMpi
        (MPI_Allreduce
         (const_cast < Complex < R > *>(sbuf),
          rbuf, count, TypeMap < Complex < R >> (),
          op.op, comm.comm));
#endif
    }
}

template void AllReduce (const byte * sbuf, byte * rbuf,
                         int count, Op op, Comm comm);
template void AllReduce (const int *sbuf, int *rbuf,
                         int count, Op op, Comm comm);
template void AllReduce (const unsigned *sbuf,
                         unsigned *rbuf, int count, Op op,
                         Comm comm);
template void AllReduce (const long int *sbuf,
                         long int *rbuf, int count, Op op,
                         Comm comm);
template void AllReduce (const unsigned long *sbuf,
                         unsigned long *rbuf, int count,
                         Op op, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllReduce (const long long int *sbuf,
                         long long int *rbuf, int count,
                         Op op, Comm comm);
template void AllReduce (const unsigned long long *sbuf,
                         unsigned long long *rbuf,
                         int count, Op op, Comm comm);
#endif
template void AllReduce (const float *sbuf, float *rbuf,
                         int count, Op op, Comm comm);
template void AllReduce (const double *sbuf, double *rbuf,
                         int count, Op op, Comm comm);
template void AllReduce (const Complex < float >*sbuf,
                         Complex < float >*rbuf,
                         int count, Op op, Comm comm);
template void AllReduce (const Complex < double >*sbuf,
                         Complex < double >*rbuf,
                         int count, Op op, Comm comm);
template void AllReduce (const ValueInt < Int > *sbuf,
                         ValueInt < Int > *rbuf,
                         int count, Op op, Comm comm);
template void AllReduce (const ValueInt < float >*sbuf,
                         ValueInt < float >*rbuf,
                         int count, Op op, Comm comm);
template void AllReduce (const ValueInt < double >*sbuf,
                         ValueInt < double >*rbuf,
                         int count, Op op, Comm comm);
template void AllReduce (const ValueIntPair < Int > *sbuf,
                         ValueIntPair < Int > *rbuf,
                         int count, Op op, Comm comm);
template void AllReduce (const ValueIntPair <
                         float >*sbuf,
                         ValueIntPair < float >*rbuf,
                         int count, Op op, Comm comm);
template void AllReduce (const ValueIntPair <
                         double >*sbuf,
                         ValueIntPair < double >*rbuf,
                         int count, Op op, Comm comm);

template < typename T >
void AllReduce (const T * sbuf, T * rbuf, int count,
                Comm comm)
{
    AllReduce (sbuf, rbuf, count, mpi::SUM, comm);
}

template void AllReduce (const byte * sbuf, byte * rbuf,
                         int count, Comm comm);
template void AllReduce (const int *sbuf, int *rbuf,
                         int count, Comm comm);
template void AllReduce (const unsigned *sbuf,
                         unsigned *rbuf, int count,
                         Comm comm);
template void AllReduce (const long int *sbuf,
                         long int *rbuf, int count,
                         Comm comm);
template void AllReduce (const unsigned long *sbuf,
                         unsigned long *rbuf, int count,
                         Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllReduce (const long long int *sbuf,
                         long long int *rbuf, int count,
                         Comm comm);
template void AllReduce (const unsigned long long *sbuf,
                         unsigned long long *rbuf,
                         int count, Comm comm);
#endif
template void AllReduce (const float *sbuf, float *rbuf,
                         int count, Comm comm);
template void AllReduce (const double *sbuf, double *rbuf,
                         int count, Comm comm);
template void AllReduce (const Complex < float >*sbuf,
                         Complex < float >*rbuf,
                         int count, Comm comm);
template void AllReduce (const Complex < double >*sbuf,
                         Complex < double >*rbuf,
                         int count, Comm comm);
template void AllReduce (const ValueInt < Int > *sbuf,
                         ValueInt < Int > *rbuf,
                         int count, Comm comm);
template void AllReduce (const ValueInt < float >*sbuf,
                         ValueInt < float >*rbuf,
                         int count, Comm comm);
template void AllReduce (const ValueInt < double >*sbuf,
                         ValueInt < double >*rbuf,
                         int count, Comm comm);
template void AllReduce (const ValueIntPair < Int > *sbuf,
                         ValueIntPair < Int > *rbuf,
                         int count, Comm comm);
template void AllReduce (const ValueIntPair <
                         float >*sbuf,
                         ValueIntPair < float >*rbuf,
                         int count, Comm comm);
template void AllReduce (const ValueIntPair <
                         double >*sbuf,
                         ValueIntPair < double >*rbuf,
                         int count, Comm comm);

template < typename T > T AllReduce (T sb, Op op,
                                     Comm comm)
{
    T rb;

    AllReduce (&sb, &rb, 1, op, comm);
    return rb;
}

template byte AllReduce (byte sb, Op op, Comm comm);
template int AllReduce (int sb, Op op, Comm comm);
template unsigned AllReduce (unsigned sb, Op op,
                             Comm comm);
template long int AllReduce (long int sb, Op op,
                             Comm comm);
template unsigned long AllReduce (unsigned long sb, Op op,
                                  Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int AllReduce (long long int sb, Op op,
                                  Comm comm);
template unsigned long long AllReduce (unsigned long long
                                       sb, Op op,
                                       Comm comm);
#endif
template float AllReduce (float sb, Op op, Comm comm);
template double AllReduce (double sb, Op op, Comm comm);
template Complex < float >AllReduce (Complex < float >sb,
                                     Op op, Comm comm);
template Complex < double >AllReduce (Complex <
                                      double >sb, Op op,
                                      Comm comm);
template ValueInt < Int > AllReduce (ValueInt < Int > sb,
                                     Op op, Comm comm);
template ValueInt < float >AllReduce (ValueInt <
                                      float >sb, Op op,
                                      Comm comm);
template ValueInt < double >AllReduce (ValueInt <
                                       double >sb, Op op,
                                       Comm comm);
template ValueIntPair < Int > AllReduce (ValueIntPair <
        Int > sb, Op op,
        Comm comm);
template ValueIntPair < float >AllReduce (ValueIntPair <
        float >sb,
        Op op,
        Comm comm);
template ValueIntPair < double >AllReduce (ValueIntPair <
        double >sb,
        Op op,
        Comm comm);

template < typename T > T AllReduce (T sb, Comm comm)
{
    return AllReduce (sb, mpi::SUM, comm);
}

template byte AllReduce (byte sb, Comm comm);
template int AllReduce (int sb, Comm comm);
template unsigned AllReduce (unsigned sb, Comm comm);
template long int AllReduce (long int sb, Comm comm);
template unsigned long AllReduce (unsigned long sb,
                                  Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int AllReduce (long long int sb,
                                  Comm comm);
template unsigned long long AllReduce (unsigned long long
                                       sb, Comm comm);
#endif
template float AllReduce (float sb, Comm comm);
template double AllReduce (double sb, Comm comm);
template Complex < float >AllReduce (Complex < float >sb,
                                     Comm comm);
template Complex < double >AllReduce (Complex <
                                      double >sb,
                                      Comm comm);
template ValueInt < Int > AllReduce (ValueInt < Int > sb,
                                     Comm comm);
template ValueInt < float >AllReduce (ValueInt <
                                      float >sb,
                                      Comm comm);
template ValueInt < double >AllReduce (ValueInt <
                                       double >sb,
                                       Comm comm);
template ValueIntPair < Int > AllReduce (ValueIntPair <
        Int > sb,
        Comm comm);
template ValueIntPair < float >AllReduce (ValueIntPair <
        float >sb,
        Comm comm);
template ValueIntPair < double >AllReduce (ValueIntPair <
        double >sb,
        Comm comm);

template < typename T >
void AllReduce (T * buf, int count, Op op, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::AllReduce"))
    if (count != 0)
    {
#ifdef EL_HAVE_MPI_IN_PLACE
        SafeMpi
        (MPI_Allreduce
         (MPI_IN_PLACE, buf, count,
          TypeMap < T > (), op.op, comm.comm));
#else
        vector<T> sendBuf( count );
        MemCopy( sendBuf.data(), buf, count );
        SafeMpi
        (MPI_Allreduce
         (sendBuf.data (), buf, count,
          TypeMap < T > (), op.op, comm.comm));
#endif
    }
}

template < typename R >
void AllReduce (Complex < R > *buf, int count, Op op,
                Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::AllReduce"))
    if (count != 0)
    {
#ifdef EL_AVOID_COMPLEX_MPI
        if (op == SUM)
        {
#ifdef EL_HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Allreduce
              ( MPI_IN_PLACE, buf, 2*count, TypeMap<R>(), op.op, comm.comm ) );
# else
            vector<Complex<R>> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            (MPI_Allreduce
             (sendBuf.data (), buf, 2 * count,
              TypeMap < R > (), op.op,
              comm.comm));
#endif
        }
        else
        {
#ifdef EL_HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Allreduce
              ( MPI_IN_PLACE, buf, count, TypeMap<Complex<R>>(), 
                op.op, comm.comm ) );
# else
            vector<Complex<R>> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            (MPI_Allreduce
             (sendBuf.data (), buf, count,
              TypeMap < Complex < R >> (), op.op,
              comm.comm));
#endif
        }
#else
#ifdef EL_HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Allreduce
          ( MPI_IN_PLACE, buf, count, TypeMap<Complex<R>>(), op.op, 
            comm.comm ) );
# else
        vector<Complex<R>> sendBuf( count );
        MemCopy( sendBuf.data(), buf, count );
        SafeMpi
        (MPI_Allreduce
         (sendBuf.data (), buf, count,
          TypeMap < Complex < R >> (), op.op,
          comm.comm));
#endif
#endif
    }
}

template void AllReduce (byte * buf, int count, Op op,
                         Comm comm);
template void AllReduce (int *buf, int count, Op op,
                         Comm comm);
template void AllReduce (unsigned *buf, int count, Op op,
                         Comm comm);
template void AllReduce (long int *buf, int count, Op op,
                         Comm comm);
template void AllReduce (unsigned long *buf, int count,
                         Op op, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllReduce (long long int *buf, int count,
                         Op op, Comm comm);
template void AllReduce (unsigned long long *buf,
                         int count, Op op, Comm comm);
#endif
template void AllReduce (float *buf, int count, Op op,
                         Comm comm);
template void AllReduce (double *buf, int count, Op op,
                         Comm comm);
template void AllReduce (Complex < float >*buf, int count,
                         Op op, Comm comm);
template void AllReduce (Complex < double >*buf,
                         int count, Op op, Comm comm);
template void AllReduce (ValueInt < Int > *buf, int count,
                         Op op, Comm comm);
template void AllReduce (ValueInt < float >*buf,
                         int count, Op op, Comm comm);
template void AllReduce (ValueInt < double >*buf,
                         int count, Op op, Comm comm);
template void AllReduce (ValueIntPair < Int > *buf,
                         int count, Op op, Comm comm);
template void AllReduce (ValueIntPair < float >*buf,
                         int count, Op op, Comm comm);
template void AllReduce (ValueIntPair < double >*buf,
                         int count, Op op, Comm comm);

template < typename T >
void AllReduce (T * buf, int count, Comm comm)
{
    AllReduce (buf, count, mpi::SUM, comm);
}

template void AllReduce (byte * buf, int count,
                         Comm comm);
template void AllReduce (int *buf, int count, Comm comm);
template void AllReduce (unsigned *buf, int count,
                         Comm comm);
template void AllReduce (long int *buf, int count,
                         Comm comm);
template void AllReduce (unsigned long *buf, int count,
                         Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllReduce (long long int *buf, int count,
                         Comm comm);
template void AllReduce (unsigned long long *buf,
                         int count, Comm comm);
#endif
template void AllReduce (float *buf, int count,
                         Comm comm);
template void AllReduce (double *buf, int count,
                         Comm comm);
template void AllReduce (Complex < float >*buf, int count,
                         Comm comm);
template void AllReduce (Complex < double >*buf,
                         int count, Comm comm);
template void AllReduce (ValueInt < Int > *buf, int count,
                         Comm comm);
template void AllReduce (ValueInt < float >*buf,
                         int count, Comm comm);
template void AllReduce (ValueInt < double >*buf,
                         int count, Comm comm);
template void AllReduce (ValueIntPair < Int > *buf,
                         int count, Comm comm);
template void AllReduce (ValueIntPair < float >*buf,
                         int count, Comm comm);
template void AllReduce (ValueIntPair < double >*buf,
                         int count, Comm comm);

template < typename R >
void ReduceScatter (R * sbuf, R * rbuf, int rc, Op op,
                    Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::ReduceScatter"))
#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = Size (comm);
    const int commRank = Rank (comm);

    AllReduce (sbuf, rc * commSize, op, comm);
    MemCopy (rbuf, &sbuf[commRank * rc], rc);
#elif defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
    SafeMpi
    (MPI_Reduce_scatter_block
     (sbuf, rbuf, rc, TypeMap < R > (), op.op,
      comm.comm));
#else
    const int commSize = Size (comm);

    Reduce (sbuf, rc * commSize, op, 0, comm);
    Scatter (sbuf, rc, rbuf, rc, 0, comm);
#endif
}

template < typename R >
void ReduceScatter
(Complex < R > *sbuf, Complex < R > *rbuf, int rc,
 Op op, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::ReduceScatter"))
#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = Size (comm);
    const int commRank = Rank (comm);

    AllReduce (sbuf, rc * commSize, op, comm);
    MemCopy (rbuf, &sbuf[commRank * rc], rc);
#elif defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    (MPI_Reduce_scatter_block
     (sbuf, rbuf, 2 * rc, TypeMap < R > (), op.op,
      comm.comm));
#else
    SafeMpi
    (MPI_Reduce_scatter_block
     (sbuf, rbuf, rc, TypeMap < Complex < R >> (),
      op.op, comm.comm));
#endif
#else
    const int commSize = Size (comm);

    Reduce (sbuf, rc * commSize, op, 0, comm);
    Scatter (sbuf, rc, rbuf, rc, 0, comm);
#endif
}

template void ReduceScatter (byte * sbuf, byte * rbuf,
                             int rc, Op op, Comm comm);
template void ReduceScatter (int *sbuf, int *rbuf, int rc,
                             Op op, Comm comm);
template void ReduceScatter (unsigned *sbuf,
                             unsigned *rbuf, int rc,
                             Op op, Comm comm);
template void ReduceScatter (long int *sbuf,
                             long int *rbuf, int rc,
                             Op op, Comm comm);
template void ReduceScatter (unsigned long *sbuf,
                             unsigned long *rbuf, int rc,
                             Op op, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void ReduceScatter (long long int *sbuf,
                             long long int *rbuf, int rc,
                             Op op, Comm comm);
template void ReduceScatter (unsigned long long *sbuf,
                             unsigned long long *rbuf,
                             int rc, Op op, Comm comm);
#endif
template void ReduceScatter (float *sbuf, float *rbuf,
                             int rc, Op op, Comm comm);
template void ReduceScatter (double *sbuf, double *rbuf,
                             int rc, Op op, Comm comm);
template void ReduceScatter (Complex < float >*sbuf,
                             Complex < float >*rbuf,
                             int rc, Op op, Comm comm);
template void ReduceScatter (Complex < double >*sbuf,
                             Complex < double >*rbuf,
                             int rc, Op op, Comm comm);

template < typename T >
void ReduceScatter (T * sbuf, T * rbuf, int rc,
                    Comm comm)
{
    ReduceScatter (sbuf, rbuf, rc, mpi::SUM, comm);
}

template void ReduceScatter (byte * sbuf, byte * rbuf,
                             int rc, Comm comm);
template void ReduceScatter (int *sbuf, int *rbuf, int rc,
                             Comm comm);
template void ReduceScatter (unsigned *sbuf,
                             unsigned *rbuf, int rc,
                             Comm comm);
template void ReduceScatter (long int *sbuf,
                             long int *rbuf, int rc,
                             Comm comm);
template void ReduceScatter (unsigned long *sbuf,
                             unsigned long *rbuf, int rc,
                             Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void ReduceScatter (long long int *sbuf,
                             long long int *rbuf, int rc,
                             Comm comm);
template void ReduceScatter (unsigned long long *sbuf,
                             unsigned long long *rbuf,
                             int rc, Comm comm);
#endif
template void ReduceScatter (float *sbuf, float *rbuf,
                             int rc, Comm comm);
template void ReduceScatter (double *sbuf, double *rbuf,
                             int rc, Comm comm);
template void ReduceScatter (Complex < float >*sbuf,
                             Complex < float >*rbuf,
                             int rc, Comm comm);
template void ReduceScatter (Complex < double >*sbuf,
                             Complex < double >*rbuf,
                             int rc, Comm comm);

template < typename T > T ReduceScatter (T sb, Op op,
        Comm comm)
{
    T rb;

    ReduceScatter (&sb, &rb, 1, op, comm);
    return rb;
}

template byte ReduceScatter (byte sb, Op op, Comm comm);
template int ReduceScatter (int sb, Op op, Comm comm);
template unsigned ReduceScatter (unsigned sb, Op op,
                                 Comm comm);
template long int ReduceScatter (long int sb, Op op,
                                 Comm comm);
template unsigned long ReduceScatter (unsigned long sb,
                                      Op op, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int ReduceScatter (long long int sb,
                                      Op op, Comm comm);
template unsigned long long ReduceScatter (unsigned long
        long sb, Op op,
        Comm comm);
#endif
template float ReduceScatter (float sb, Op op, Comm comm);
template double ReduceScatter (double sb, Op op,
                               Comm comm);
template Complex < float >ReduceScatter (Complex <
        float >sb, Op op,
        Comm comm);
template Complex < double >ReduceScatter (Complex <
        double >sb,
        Op op,
        Comm comm);

template < typename T > T ReduceScatter (T sb, Comm comm)
{
    return ReduceScatter (sb, mpi::SUM, comm);
}

template byte ReduceScatter (byte sb, Comm comm);
template int ReduceScatter (int sb, Comm comm);
template unsigned ReduceScatter (unsigned sb, Comm comm);
template long int ReduceScatter (long int sb, Comm comm);
template unsigned long ReduceScatter (unsigned long sb,
                                      Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int ReduceScatter (long long int sb,
                                      Comm comm);
template unsigned long long ReduceScatter (unsigned long
        long sb,
        Comm comm);
#endif
template float ReduceScatter (float sb, Comm comm);
template double ReduceScatter (double sb, Comm comm);
template Complex < float >ReduceScatter (Complex <
        float >sb,
        Comm comm);
template Complex < double >ReduceScatter (Complex <
        double >sb,
        Comm comm);

template < typename R >
void ReduceScatter (R * buf, int rc, Op op, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::ReduceScatter"))
#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = Size (comm);
    const int commRank = Rank (comm);

    AllReduce (buf, rc * commSize, op, comm);
    if (commRank != 0)
        MemCopy (buf, &buf[commRank * rc], rc);
#elif defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
#ifdef EL_HAVE_MPI_IN_PLACE
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( MPI_IN_PLACE, buf, rc, TypeMap<R>(), op.op, comm.comm ) );
# else
    const int commSize = Size( comm );
    vector<R> sendBuf( rc*commSize );
    MemCopy( sendBuf.data(), buf, rc*commSize );
    SafeMpi
    (MPI_Reduce_scatter_block
     (sendBuf.data (), buf, rc, TypeMap < R > (),
      op.op, comm.comm));
#endif
#else
    const int commSize = Size (comm);

    Reduce (buf, rc * commSize, op, 0, comm);
    Scatter (buf, rc, rc, 0, comm);
#endif
}

// TODO: Handle case where op is not summation
template < typename R >
void ReduceScatter (Complex < R > *buf, int rc, Op op,
                    Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::ReduceScatter"))
#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = Size (comm);
    const int commRank = Rank (comm);

    AllReduce (buf, rc * commSize, op, comm);
    if (commRank != 0)
        MemCopy (buf, &buf[commRank * rc], rc);
#elif defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
#ifdef EL_AVOID_COMPLEX_MPI
#ifdef EL_HAVE_MPI_IN_PLACE
    SafeMpi
    (MPI_Reduce_scatter_block
     (MPI_IN_PLACE, buf, 2 * rc, TypeMap < R > (),
      op.op, comm.comm));
#else
    const int commSize = Size (comm);
    std::vector < Complex < R >> sendBuf (rc * commSize);
    MemCopy (sendBuf.data (), buf, rc * commSize);
    SafeMpi
    (MPI_Reduce_scatter_block
     (sendBuf.data (), buf, 2 * rc, TypeMap < R > (),
      op.op, comm.comm));
#endif
#else
#ifdef EL_HAVE_MPI_IN_PLACE
    SafeMpi
    (MPI_Reduce_scatter_block
     (MPI_IN_PLACE, buf, rc,
      TypeMap < Complex < R >> (), op.op, comm.comm));
#else
    const int commSize = Size (comm);

    std::vector < Complex < R >> sendBuf (rc * commSize);
    MemCopy (sendBuf.data (), buf, rc * commSize);
    SafeMpi
    (MPI_Reduce_scatter_block
     (sendBuf.data (), buf, rc,
      TypeMap < Complex < R >> (), op.op, comm.comm));
#endif
#endif
#else
    const int commSize = Size (comm);

    Reduce (buf, rc * commSize, op, 0, comm);
    Scatter (buf, rc, rc, 0, comm);
#endif
}

template void ReduceScatter (byte * buf, int rc, Op op,
                             Comm comm);
template void ReduceScatter (int *buf, int rc, Op op,
                             Comm comm);
template void ReduceScatter (unsigned *buf, int rc, Op op,
                             Comm comm);
template void ReduceScatter (long int *buf, int rc, Op op,
                             Comm comm);
template void ReduceScatter (unsigned long *buf, int rc,
                             Op op, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void ReduceScatter (long long int *buf, int rc,
                             Op op, Comm comm);
template void ReduceScatter (unsigned long long *buf,
                             int rc, Op op, Comm comm);
#endif
template void ReduceScatter (float *buf, int rc, Op op,
                             Comm comm);
template void ReduceScatter (double *buf, int rc, Op op,
                             Comm comm);
template void ReduceScatter (Complex < float >*buf,
                             int rc, Op op, Comm comm);
template void ReduceScatter (Complex < double >*buf,
                             int rc, Op op, Comm comm);

template < typename T >
void ReduceScatter (T * buf, int rc, Comm comm)
{
    ReduceScatter (buf, rc, mpi::SUM, comm);
}

template void ReduceScatter (byte * buf, int rc,
                             Comm comm);
template void ReduceScatter (int *buf, int rc, Comm comm);
template void ReduceScatter (unsigned *buf, int rc,
                             Comm comm);
template void ReduceScatter (long int *buf, int rc,
                             Comm comm);
template void ReduceScatter (unsigned long *buf, int rc,
                             Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void ReduceScatter (long long int *buf, int rc,
                             Comm comm);
template void ReduceScatter (unsigned long long *buf,
                             int rc, Comm comm);
#endif
template void ReduceScatter (float *buf, int rc,
                             Comm comm);
template void ReduceScatter (double *buf, int rc,
                             Comm comm);
template void ReduceScatter (Complex < float >*buf,
                             int rc, Comm comm);
template void ReduceScatter (Complex < double >*buf,
                             int rc, Comm comm);

template < typename R >
void ReduceScatter
(const R * sbuf, R * rbuf, const int *rcs, Op op,
 Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::ReduceScatter"))
    SafeMpi
    (MPI_Reduce_scatter
     (const_cast < R * >(sbuf),
      rbuf, const_cast < int *>(rcs),
      TypeMap < R > (), op.op, comm.comm));
}

template < typename R >
void ReduceScatter
(const Complex < R > *sbuf, Complex < R > *rbuf,
 const int *rcs, Op op, Comm comm)
{
    DEBUG_ONLY (CallStackEntry cse ("mpi::ReduceScatter"))
#ifdef EL_AVOID_COMPLEX_MPI
    if (op == SUM)
    {
        int p;

        MPI_Comm_size( comm.comm, &p );
        vector<int> rcsDoubled(p);
        for( int i=0; i<p; ++i )
            rcsDoubled[i] = 2*rcs[i];
        SafeMpi
        (MPI_Reduce_scatter
         (const_cast < Complex < R > *>(sbuf),
          rbuf, rcsDoubled.data (),
          TypeMap < R > (), op.op, comm.comm));
    }
    else
    {
        SafeMpi
        (MPI_Reduce_scatter
         (const_cast < Complex < R > *>(sbuf),
          rbuf, const_cast < int *>(rcs),
          TypeMap < Complex < R >> (), op.op,
          comm.comm));
    }
#else
    SafeMpi
    (MPI_Reduce_scatter
     (const_cast < Complex < R > *>(sbuf),
      rbuf, const_cast < int *>(rcs),
      TypeMap < Complex < R >> (), op.op,
      comm.comm));
#endif
}

template void ReduceScatter (const byte * sbuf,
                             byte * rbuf, const int *rcs,
                             Op op, Comm comm);
template void ReduceScatter (const int *sbuf, int *rbuf,
                             const int *rcs, Op op,
                             Comm comm);
template void ReduceScatter (const unsigned *sbuf,
                             unsigned *rbuf,
                             const int *rcs, Op op,
                             Comm comm);
template void ReduceScatter (const long int *sbuf,
                             long int *rbuf,
                             const int *rcs, Op op,
                             Comm comm);
template void ReduceScatter (const unsigned long *sbuf,
                             unsigned long *rbuf,
                             const int *rcs, Op op,
                             Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void ReduceScatter (const long long int *sbuf,
                             long long int *rbuf,
                             const int *rcs, Op op,
                             Comm comm);
template void ReduceScatter (const unsigned long long
                             *sbuf,
                             unsigned long long *rbuf,
                             const int *rcs, Op op,
                             Comm comm);
#endif
template void ReduceScatter (const float *sbuf,
                             float *rbuf, const int *rcs,
                             Op op, Comm comm);
template void ReduceScatter (const double *sbuf,
                             double *rbuf, const int *rcs,
                             Op op, Comm comm);
template void ReduceScatter (const Complex < float >*sbuf,
                             Complex < float >*rbuf,
                             const int *rcs, Op op,
                             Comm comm);
template void ReduceScatter (const Complex <
                             double >*sbuf,
                             Complex < double >*rbuf,
                             const int *rcs, Op op,
                             Comm comm);

template < typename T >
void ReduceScatter (const T * sbuf, T * rbuf,
                    const int *rcs, Comm comm)
{
    ReduceScatter (sbuf, rbuf, rcs, mpi::SUM, comm);
}

template void ReduceScatter (const byte * sbuf,
                             byte * rbuf, const int *rcs,
                             Comm comm);
template void ReduceScatter (const int *sbuf, int *rbuf,
                             const int *rcs, Comm comm);
template void ReduceScatter (const unsigned *sbuf,
                             unsigned *rbuf,
                             const int *rcs, Comm comm);
template void ReduceScatter (const long int *sbuf,
                             long int *rbuf,
                             const int *rcs, Comm comm);
template void ReduceScatter (const unsigned long *sbuf,
                             unsigned long *rbuf,
                             const int *rcs, Comm comm);
#ifdef EL_HAVE_MPI_LONG_LONG
template void ReduceScatter (const long long int *sbuf,
                             long long int *rbuf,
                             const int *rcs, Comm comm);
template void ReduceScatter (const unsigned long long
                             *sbuf,
                             unsigned long long *rbuf,
                             const int *rcs, Comm comm);
#endif

template void ReduceScatter( const float* sbuf, float* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const double* sbuf, double* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const Complex<float>* sbuf, Complex<float>* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const Complex<double>* sbuf, Complex<double>* rbuf, const int* rcs, Comm comm );

void VerifySendsAndRecvs
( const vector<int>& sendCounts,
  const vector<int>& recvCounts, mpi::Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::VerifySendsAndRecvs"))
    const int commSize = mpi::Size( comm );
    vector<int> actualRecvCounts(commSize);
    mpi::AllToAll
    ( &sendCounts[0],       1,
      &actualRecvCounts[0], 1, comm );
    for( int proc=0; proc<commSize; ++proc )
        if( actualRecvCounts[proc] != recvCounts[proc] )
            LogicError
            ("Expected recv count of ",recvCounts[proc],
             " but recv'd ",actualRecvCounts[proc]," from process ",proc);
}

template<typename T>
void SparseAllToAll
( const vector<T>& sendBuffer,
  const vector<int>& sendCounts, const vector<int>& sendDispls,
        vector<T>& recvBuffer,
  const vector<int>& recvCounts, const vector<int>& recvDispls,
        mpi::Comm comm )
{
#ifdef EL_USE_CUSTOM_ALLTOALLV
    const int commSize = mpi::Size( comm );
    int numSends=0,numRecvs=0;
    for( int proc=0; proc<commSize; ++proc )
    {
        if( sendCounts[proc] != 0 )
            ++numSends;
        if( recvCounts[proc] != 0 )
            ++numRecvs;
    }
    vector<mpi::Status> statuses(numSends+numRecvs);
    vector<mpi::Request> requests(numSends+numRecvs);
    int rCount=0;
    for( int proc=0; proc<commSize; ++proc )
    {
        int count = recvCounts[proc];
        int displ = recvDispls[proc];
        if( count != 0 )
            mpi::IRecv
            ( &recvBuffer[displ], count, proc, comm, requests[rCount++] );
    }
#ifdef EL_BARRIER_IN_ALLTOALLV
    // This should help ensure that recvs are posted before the sends
    mpi::Barrier( comm );
#endif
    for( int proc=0; proc<commSize; ++proc )
    {
        int count = sendCounts[proc];
        int displ = sendDispls[proc];
        if( count != 0 )
            mpi::ISend
            ( &sendBuffer[displ], count, proc, comm, requests[rCount++] );
    }
    mpi::WaitAll( numSends+numRecvs, &requests[0], &statuses[0] );
#else
    mpi::AllToAll
    ( &sendBuffer[0], &sendCounts[0], &sendDispls[0],
      &recvBuffer[0], &recvCounts[0], &recvDispls[0], comm );
#endif
}

#define PROTO(T) \
  template void SparseAllToAll \
  ( const vector<T>& sendBuffer, \
    const vector<int>& sendCounts, const vector<int>& sendDispls, \
          vector<T>& recvBuffer, \
    const vector<int>& recvCounts, const vector<int>& recvDispls, \
          mpi::Comm comm );
#include "El/macros/Instantiate.h"

} // namespace mpi
} // namespace El
