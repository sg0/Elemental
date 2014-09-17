/*
This file is part of Elemental and is under the BSD 2-Clause License,
which can be found in the LICENSE file in the root directory, or at
http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include <cassert>

// TODO Use DDT for put/get/acc when EL_USE_DERIVED_TYPE is defined
// TODO bring back const interfaces
namespace El
{
template<typename T>
AxpyInterface2<T>::AxpyInterface2()
    : GlobalArrayPut_(0), GlobalArrayGet_(0),
    matrices_(0), coords_(0),
    toBeAttachedForGet_(false), toBeAttachedForPut_(false),
    attached_(false), detached_(false)
{ }

template<typename T>
AxpyInterface2<T>::AxpyInterface2( DistMatrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::AxpyInterface2"))

    attached_ 		= false;
    detached_ 		= false;
    toBeAttachedForGet_ = true;
    toBeAttachedForPut_ = true;
    GlobalArrayPut_ 	= &Z;
    GlobalArrayGet_ 	= &Z;

    const Grid& g = Z.Grid();
	const Int p = g.Size ();
    
    if ( matrices_.empty() )
    {
	struct matrix_params_ mp;
	mp.data_.resize(p);
	mp.requests_.resize(p);
	mp.statuses_.resize(p);   
	mp.base_ = NULL;
	// push back new matrix_params created
	// with default constructor
	matrices_.push_back( mp );
    }
   
    if ( coords_.empty() )
    {
	struct coord_params_ cp;
	cp.coord_.resize(p);
	cp.requests_.resize(p);
	cp.statuses_.resize(p);   
	cp.base_ = NULL;
	// push back new matrix_params created
	// with default constructor
	coords_.push_back( cp );
    }
}

template<typename T>
AxpyInterface2<T>::~AxpyInterface2()
{
        if( std::uncaught_exception() )
        {
            std::ostringstream os;
            os << "Uncaught exception detected during AxpyInterface2 destructor "
            "that required a call to Detach. Instead of allowing for the "
            "possibility of Detach throwing another exception and "
            "resulting in a 'terminate', we instead immediately dump the "
            "call stack (if not in RELEASE mode) since the program will "
            "likely hang:" << std::endl;
            std::cerr << os.str();
            DEBUG_ONLY(DumpCallStack())
        }
        else
        {
            Detach();
        }
}

template<typename T>
void AxpyInterface2<T>::Attach( DistMatrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Attach"))
    // attached_ will be only set in Attach
    // and only unset in Detach
    if (!attached_)
        attached_ = true;
    else
        LogicError("Must detach before reattaching.");
	
    // the matrix base_ is not known until 
    // an update operation (put/get/acc)
    // so it is kept blank

    // if DistMatrix is non-const, all one-sided
    // transfers -- put, get and acc are possible
    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
    {
	GlobalArrayPut_ 	= &Z;
	toBeAttachedForPut_ 	= true;
	GlobalArrayGet_ 	= &Z;
	toBeAttachedForGet_ 	= true;

	const Grid& g = Z.Grid();
	const Int p = g.Size ();

	if ( matrices_.empty() )
	{
	    struct matrix_params_ mp;
	    mp.data_.resize(p);
	    mp.requests_.resize(p);
	    mp.statuses_.resize(p);   
	    mp.base_ = NULL;
	    // push back new matrix_params created
	    // with default constructor
	    matrices_.push_back( mp );
	}

	if ( coords_.empty() )
	{
	    struct coord_params_ cp;
	    cp.coord_.resize(p);
	    cp.requests_.resize(p);
	    cp.statuses_.resize(p);   
	    cp.base_ = NULL;
	    // push back new matrix_params created
	    // with default constructor
	    coords_.push_back( cp );
	}
    }
}

template<typename T>
Int AxpyInterface2<T>::NextIndexMatrix (
	Int target,
	Int dataSize, 
	T * base_address,
	Int *mindex)
{
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface2::NextIndexMatrix"))
    
    assert ( base_address != NULL );

    Int matrixIndex = 0;
    DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    const Int p = g.Size();
    const Int numMatrices = matrices_.size();

    // search for matrix base
    for (Int m = 0; m < numMatrices; m++)
    {
	if ( matrices_[m].base_ == base_address )
	{
	    matrixIndex = m;
	    break;
	}
	if ( matrices_[m].base_ == NULL )
	{
	    matrices_[m].base_ = base_address;
	    matrixIndex = m;
	    break;    
	}
	matrixIndex = m+1;
    }
    
    // need to create new object
    if ( matrixIndex == numMatrices)
    {
	struct matrix_params_ mp;
	mp.data_.resize(p);
	mp.requests_.resize(p);
	mp.statuses_.resize(p);   
	mp.base_ = NULL;
	// push back new matrix_params created
	// with default constructor
	matrices_.push_back( mp );
	matrices_[matrixIndex].base_ = base_address;
    }
    // go through the request, data, 
    // status objects
    const Int numCreated = matrices_[matrixIndex].data_[target].size ();
    DEBUG_ONLY (if (numCreated != Int (matrices_[matrixIndex].requests_[target].size ()) ||
		numCreated != Int (matrices_[matrixIndex].statuses_[target].size ()))
	    LogicError ("size mismatch");)

	for (Int i = 0; i < numCreated; ++i)
	{
	    // If this request is still running, 
	    // test to see if it finished.
	    if (matrices_[matrixIndex].statuses_[target][i])
	    {
		const bool finished = mpi::Test (matrices_[matrixIndex].requests_[target][i]);
		matrices_[matrixIndex].statuses_[target][i] = !finished;
	    }

	    if (!matrices_[matrixIndex].statuses_[target][i])
	    {
		matrices_[matrixIndex].statuses_[target][i] = true;
		matrices_[matrixIndex].data_[target][i].resize ( dataSize );
		*mindex = matrixIndex;
		return i;
	    }
	}

    matrices_[matrixIndex].data_[target].resize ( numCreated + 1 );
    matrices_[matrixIndex].data_[target][numCreated].resize ( dataSize );
    matrices_[matrixIndex].requests_[target].push_back ( mpi::REQUEST_NULL );
    matrices_[matrixIndex].statuses_[target].push_back ( true );
    *mindex = matrixIndex;

    return numCreated;
}

template<typename T>
Int AxpyInterface2<T>::NextIndexCoord (
	Int i, Int j,
	Int target,
	T * base_address,
	Int *cindex)
{
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface2::NextIndexCoord"))
    
    assert ( base_address != NULL );

    Int coordIndex = 0;
    DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    const Int p = g.Size();
    const Int numCoords = coords_.size();

    // search for matrix base
    for (Int m = 0; m < numCoords; m++)
    {
	if ( coords_[m].base_ == base_address )
	{
	    coordIndex = m;
	    break;
	}
	if ( coords_[m].base_ == NULL )
	{
	    coords_[m].base_ = base_address;
	    coordIndex = m;
	    break;    
	}
	coordIndex = m+1;
    }
    
    // need to create new object
    if ( coordIndex == numCoords )
    {
	struct coord_params_ cp;
	cp.coord_.resize(p);
	cp.requests_.resize(p);
	cp.statuses_.resize(p);   
	cp.base_ = NULL;
	// push back new matrix_params created
	// with default constructor
	coords_.push_back( cp );
	coords_[coordIndex].base_ = base_address;
    }
    // go through the request, data, 
    // status objects
    const Int numCreated = coords_[coordIndex].coord_[target].size ();
    DEBUG_ONLY (if (numCreated != Int (coords_[coordIndex].requests_[target].size ()) ||
		numCreated != Int (matrices_[coordIndex].statuses_[target].size ()))
	    LogicError ("size mismatch");)

	for (Int i = 0; i < numCreated; ++i)
	{
	    // If this request is still running, 
	    // test to see if it finished.
	    if (coords_[coordIndex].statuses_[target][i])
	    {
		const bool finished = mpi::Test (coords_[coordIndex].requests_[target][i]);
		coords_[coordIndex].statuses_[target][i] = !finished;
	    }

	    if (!coords_[coordIndex].statuses_[target][i])
	    {
		coords_[coordIndex].statuses_[target][i] = true;
		coords_[coordIndex].coord_[target][i][0] = i;
		coords_[coordIndex].coord_[target][i][1] = j;
		*cindex = coordIndex;
		return i;
	    }
	}

    coords_[coordIndex].coord_[target].resize ( numCreated + 1 );
    coords_[coordIndex].coord_[target][numCreated][0] = i;
    coords_[coordIndex].coord_[target][numCreated][1] = j;
    coords_[coordIndex].requests_[target].push_back ( mpi::REQUEST_NULL );
    coords_[coordIndex].statuses_[target].push_back ( true );
    *cindex = coordIndex;

    return numCreated;
}

template<typename T>
void AxpyInterface2<T>::Iput( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Iput"))

    if( i < 0 || j < 0 )
        LogicError("Submatrix offsets must be non-negative");
    if ( !toBeAttachedForPut_ )
        LogicError("Global matrix cannot be updated");

    DistMatrix<T>& Y = *GlobalArrayPut_;
    //do boundary checks
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError("Submatrix out of bounds of global matrix");

    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = (Y.ColAlign() + i) % r;
    const Int rowAlign = (Y.RowAlign() + j) % c;

    const Int XLDim = Z.LDim();
    // local matrix width and height
    const Int height = Z.Height();
    const Int width = Z.Width();

    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;

    const Int YLDim = Y.LDim ();
	
    Int matrix_index, coord_index;

    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        // number of entries in my PE
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight * localWidth;
	
        if( numEntries != 0 )
        {
            const Int destination = receivingRow + r*receivingCol;
            T* XBuffer = Z.Buffer();
	    const Int dindex =
		NextIndexMatrix (destination,
			numEntries, 
			XBuffer,
			&matrix_index);

	    DEBUG_ONLY (if
			(Int (matrices_[matrix_index].data_[destination][dindex].size ()) !=
			 numEntries) LogicError ("Error in NextIndexMatrix");)

	    T *sendBuffer = matrices_[matrix_index].data_[destination][dindex].data ();
            for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }
	    // put request 
	    mpi::TaggedISend (sendBuffer, numEntries, destination, 
			DATA_PUT_TAG, g.VCComm (), 
			matrices_[matrix_index].requests_[destination][dindex]);
		    
	    // send coordinates
	    const Int cindex =
		NextIndexCoord (i, j,
			destination,
			XBuffer,
			&coord_index);

	    Int *coord = reinterpret_cast<Int *>(coords_[coord_index].coord_[destination][cindex].data ());
	    coord[0] = i; coord[1] = j;
	    mpi::TaggedISend (coord, 2, destination, COORD_PUT_TAG, g.VCComm (), 
		    coords_[coord_index].requests_[destination][cindex]);
	}
        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

template<typename T>
void AxpyInterface2<T>::Iget( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Iget"))
    // a call to Attach with a non-const DistMatrix must set
    // toBeAttachedForGet_ also, if not then it is assumed that
    // the DistMatrix isn't attached
    if ( !toBeAttachedForGet_ )
	LogicError ("Cannot perform this operation as matrix is not attached.");
    DistMatrix<T>& X = *GlobalArrayGet_;

    const Int height = Z.Height ();
    const Int width = Z.Width ();
    
    if (i + height > X.Height () || j + width > X.Width ())
	LogicError ("Invalid AxpyGlobalToLocal submatrix");

    const Grid & g = X.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int p = g.Size ();

    std::vector<T> recvVector_;
    Int coord_index;

    T* XBuffer = Z.Buffer();
    // Send out the requests to all processes in the grid
    for (Int rank = 0; rank < p; ++rank)
    {
	// send coordinates, no need to send a separate
	// request object
	const Int cindex =
		NextIndexCoord (i, j,
			rank,
			XBuffer,
			&coord_index);
	Int *coord = reinterpret_cast<Int *>(coords_[coord_index].coord_[rank][cindex].data ());
        coord[0] = i; coord[1] = j;
	mpi::TaggedISend (coord, 2, rank, 
		REQUEST_GET_TAG, g.VCComm (), 
		coords_[coord_index].requests_[rank][cindex]);
    }

    // Receive all of the replies
    Int numReplies = 0;
    while (numReplies < p)
    {
	mpi::Status status;
	HandleGlobalToLocalData ( Z );
	if (mpi::IProbe
		(mpi::ANY_SOURCE, DATA_GET_TAG, g.VCComm (), status))
	{
	    const Int source = status.MPI_SOURCE;
	    // Ensure that we have a recv buffer
	    const Int count = mpi::GetCount <T> (status);
	    recvVector_.resize (count);
	    T *recvBuffer = recvVector_.data ();

	    // Receive the data
	    mpi::TaggedRecv
		(recvBuffer, count, source, DATA_GET_TAG, g.VCComm ());

	    // Compute the local heights and offsets
	    const Int myRow = g.Row ();
	    const Int myCol = g.Col ();
	    const Int colAlign = (X.ColAlign () + i) % r;
	    const Int rowAlign = (X.RowAlign () + j) % c;
	    const Int colShift = Shift (myRow, colAlign, r);
	    const Int rowShift = Shift (myCol, rowAlign, c);
	    const Int localHeight = Length (height, colShift, r);
	    const Int localWidth = Length (width, rowShift, c);

	    // Unpack the local matrix
	    for (Int t = 0; t < localWidth; ++t)
	    {
		T *YCol = X.Buffer (0, rowShift + t * c);
		const T *XCol = &recvBuffer[t * localHeight];
		for (Int s = 0; s < localHeight; ++s)
		    YCol[colShift + s * r] = XCol[s];
	    }
	    ++numReplies;
	    recvVector_.clear();
	}
    }
}

// accumulate = Update Y(i:i+height-1,j:j+width-1) += X,
// where X is height x width
template<typename T>
void AxpyInterface2<T>::Iacc( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Iacc"))

    if( i < 0 || j < 0 )
        LogicError("Submatrix offsets must be non-negative");
    if ( !toBeAttachedForPut_ )
        LogicError("Global matrix cannot be updated");

    DistMatrix<T>& Y = *GlobalArrayPut_;
   
    //do boundary checks
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
        LogicError("Submatrix out of bounds of global matrix");

    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();
    const Int colAlign = (Y.ColAlign() + i) % r;
    const Int rowAlign = (Y.RowAlign() + j) % c;
 
    Int matrix_index, coord_index;
    
    const Int XLDim = Z.LDim();
    // local matrix width and height
    const Int height = Z.Height();
    const Int width = Z.Width();

    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;

    const Int YLDim = Y.LDim();

    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
        // number of entries in my PE
        const Int localHeight = Length( height, colShift, r );
        const Int localWidth = Length( width, rowShift, c );
        const Int numEntries = localHeight * localWidth;
	
        if( numEntries != 0  )
        {
            const Int destination = receivingRow + r*receivingCol;
            T* XBuffer = Z.Buffer();
	    
	    // send data 
     	    const Int dindex =
		NextIndexMatrix (destination,
			numEntries, 
			XBuffer,
			&matrix_index);

	    DEBUG_ONLY (if
			(Int (matrices_[matrix_index].data_[destination][dindex].size ()) !=
			 numEntries) LogicError ("Error in NextIndexMatrix");)
	    
	    T *sendBuffer = matrices_[matrix_index].data_[destination][dindex].data ();
	    for( Int t=0; t<localWidth; ++t )
            {
                T* thisSendCol = &sendBuffer[t*localHeight];
                T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
                for( Int s=0; s<localHeight; ++s )
                    thisSendCol[s] = thisXCol[colShift+s*r];
            }

	    mpi::TaggedISend (sendBuffer, numEntries, destination, 
	    	DATA_ACC_TAG, g.VCComm(), 
	   	matrices_[matrix_index].requests_[destination][dindex]);
	
	   // send coordinates
	    const Int cindex =
		NextIndexCoord (i, j,
			destination,
			XBuffer,
			&coord_index);
	    
	    Int *coord = reinterpret_cast<Int *>(coords_[coord_index].coord_[destination][cindex].data());
	    coord[0] = i; coord[1] = j;
	    mpi::TaggedISend (coord, 2, destination, 
		    COORD_ACC_TAG, g.VCComm(), 
		    coords_[coord_index].requests_[destination][cindex]);
	}
        receivingRow = (receivingRow + 1) % r;
        if( receivingRow == 0 )
            receivingCol = (receivingCol + 1) % c;
    }
}

// wait 
template<typename T>
void AxpyInterface2<T>::WaitMatrix ( Matrix<T>& Z )
{
    DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    const Int p = g.Size();
    const Int numMatrices = matrices_.size();
    Int matrixIndex;
    const T *base_address = Z.LockedBuffer();

    // search for matrix base
    for (Int m = 0; m < numMatrices; m++)
    {
	if ( matrices_[m].base_ == base_address )
	{
	    matrixIndex = m;
	    break;
	}
	matrixIndex = m+1;
    }
    // matrix not found
    if ( matrixIndex == numMatrices)
	return;

    for (int rank = 0; rank < p; ++rank)
    {
	if ( matrices_[matrixIndex].statuses_[rank].size() == 0 )
	    continue;
	const Int numRequests = matrices_[matrixIndex].requests_[rank].size ();
	for (int i = 0; i < numRequests; i++)
	{
	    mpi::Wait ( matrices_[matrixIndex].requests_[rank][i] );
	    matrices_[matrixIndex].statuses_[rank][i] = false;
	}
    }
}

// progress communication for a particular matrix
// progress requests
template<typename T>
bool AxpyInterface2<T>::TestMatrix ( Matrix<T>& Z )
{
    DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    const Int p = g.Size();
    const Int numMatrices = matrices_.size();
    Int matrixIndex;
    const T *base_address = Z.LockedBuffer();

    // search for matrix base
    for (Int m = 0; m < numMatrices; m++)
    {
	if ( matrices_[m].base_ == base_address )
	{
	    matrixIndex = m;
	    break;
	}
	matrixIndex = m+1;
    }

    // matrix not found
    if ( matrixIndex == numMatrices)
	return true;

    for (int rank = 0; rank < p; ++rank)
    {
	if ( matrices_[matrixIndex].statuses_[rank].size() == 0 )
	    continue;
	const Int numStatuses = matrices_[matrixIndex].requests_[rank].size ();
	for (int i = 0; i < numStatuses; i++)
	{
	    matrices_[matrixIndex].statuses_[rank][i] = !mpi::Test ( matrices_[matrixIndex].requests_[rank][i] );
	    if ( matrices_[matrixIndex].statuses_[rank][i] )
		return false;
	}
    }
    return true;
}

template<typename T>
bool AxpyInterface2<T>::TestCoord ( Matrix<T>& Z )
{
    DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    const Int p = g.Size();
    const Int numCoords = coords_.size();
    Int coordIndex;
    const T *base_address = Z.LockedBuffer();

    // search for coord base
    for (Int m = 0; m < numCoords; m++)
    {
	if ( coords_[m].base_ == base_address )
	{
	    coordIndex = m;
	    break;
	}
	coordIndex = m+1;
    }

    // coord not found
    if ( coordIndex == numCoords)
	return true;

    for (int rank = 0; rank < p; ++rank)
    {
	if ( coords_[coordIndex].statuses_[rank].size() == 0 )
	    continue;
	const Int numStatuses = coords_[coordIndex].requests_[rank].size ();
	for (int i = 0; i < numStatuses; i++)
	{
	    coords_[coordIndex].statuses_[rank][i] = !mpi::Test ( coords_[coordIndex].requests_[rank][i] );
	    if ( coords_[coordIndex].statuses_[rank][i] )
		return false;
	}
    }
    return true;
}

// flush ensures local and remote completion
// this interface assumes a send has been issued
// and will post a matching receive and progress
template<typename T>
void AxpyInterface2<T>::Flush( Matrix<T>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Flush"))
    if( !toBeAttachedForPut_ && !toBeAttachedForGet_ )
	LogicError("Must initiate transfer before flushing.");

    DistMatrix<T>& Y = *GlobalArrayPut_;
    const Grid& g = Y.Grid();
    
    mpi::Status status;
    bool DONE = false;
    mpi::Request nb_bar_request;
    bool nb_bar_active = false;
    
    while ( !DONE )
    {
	// similar to HandleXYZ functions in original AxpyInterface
	if ( mpi::IProbe (mpi::ANY_SOURCE, mpi::ANY_TAG, g.VCComm(), status) )
	{
	    switch (status.MPI_TAG)
	    {
		case DATA_PUT_TAG:
		    {
			const Int count = mpi::GetCount <T> (status);
			HandleLocalToGlobalData ( Z, count, status.MPI_SOURCE );
			break;
		    }
		case DATA_ACC_TAG:
		    {
			const Int count = mpi::GetCount <T> (status);
			HandleLocalToGlobalAcc ( Z, count, status.MPI_SOURCE );
			break;
		    }
		case REQUEST_GET_TAG:
		    {
			HandleGlobalToLocalData ( Z );
			break;
		    }
	    }
	}
	
	if ( nb_bar_active )
	{
	    DONE = mpi::Test ( nb_bar_request );
	}
	else
	{
	    // check if all sends (data or request) are 
	    // complete for a particular matrix
	    if ( TestMatrix( Z ) && TestCoord( Z ) )
	    {
		mpi::IBarrier ( g.VCComm(), nb_bar_request );
		nb_bar_active = true;
	    }
	}
    }
}

template < typename T > 
void AxpyInterface2<T>::HandleLocalToGlobalData ( Matrix<T>& Z, Int count, Int source )
{
    DistMatrix<T> &Y = *GlobalArrayPut_;
    const Grid & g = Y.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int myRow = g.Row ();
    const Int myCol = g.Col ();
    int height = Z.Height();
    int width = Z.Width();
    // data vector
    std::vector<T> getVector_;
    getVector_.resize (count);

    DEBUG_ONLY (if (count < Int (sizeof (T)))
	    LogicError ("Count was too small");)
    DEBUG_ONLY (if (Int (getVector_.size ()) != count) 
	    LogicError ("Not enough space allocated");)

    // post receive for coordinates
    Int coord[2];
    mpi::TaggedRecv (coord, 2, source, 
	    COORD_PUT_TAG, g.VCComm());
    Int i = coord[0]; 
    Int j = coord[1];
 
    // post receive for data
    T *getBuffer = getVector_.data();
    mpi::TaggedRecv (getBuffer, count, source, 
	    DATA_PUT_TAG, g.VCComm());
        
    // Update Y
    const T *XBuffer = reinterpret_cast < const T * >(getBuffer);
    const Int colAlign = (Y.ColAlign() + i) % r;
    const Int rowAlign = (Y.RowAlign() + j) % c;
    const Int colShift = Shift (myRow, colAlign, r);
    const Int rowShift = Shift (myCol, rowAlign, c);

    const Int localHeight = Length (height, colShift, r);
    const Int localWidth = Length (width, rowShift, c);
    const Int iLocalOffset = Length (i, Y.ColShift(), r);
    const Int jLocalOffset = Length (j, Y.RowShift(), c);

    for (Int t = 0; t < localWidth; ++t)
    {
	T *YCol = Y.Buffer (iLocalOffset, jLocalOffset + t);
	const T *XCol = &XBuffer[t * localHeight];
	for (Int s = 0; s < localHeight; ++s)
	    YCol[s] = XCol[s];
    }
    // Free the memory
    getVector_.clear();
}

// replica of above function except this accumulates
template < typename T > 
void AxpyInterface2<T>::HandleLocalToGlobalAcc ( Matrix<T>& Z, Int count, Int source )
{
    DistMatrix<T> &Y = *GlobalArrayPut_;
    const Grid & g = Y.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int myRow = g.Row ();
    const Int myCol = g.Col ();
    const int height = Z.Height();
    const int width = Z.Width();

    // data buffer
    std::vector<T> getVector_;
    getVector_.resize (count);

    DEBUG_ONLY (if (count < Int (sizeof (T)))
	    LogicError ("Count was too small");)

    DEBUG_ONLY (if (Int (getVector_.size ()) != count) 
	    LogicError ("Not enough space allocated");)
    
    // post receive for coordinates
    Int coord[2];
    mpi::TaggedRecv (coord, 2, source, 
	    COORD_ACC_TAG, g.VCComm());
    Int i = coord[0]; Int j = coord[1];
    
    // post receive for data
    T *getBuffer = getVector_.data();
    mpi::TaggedRecv (getBuffer, count, source, 
	    DATA_ACC_TAG, g.VCComm());

    // Update Y
    const T *XBuffer = reinterpret_cast < const T * >(getBuffer);
    const Int colAlign = (Y.ColAlign() + i) % r;
    const Int rowAlign = (Y.RowAlign() + j) % c;
    const Int colShift = Shift (myRow, colAlign, r);
    const Int rowShift = Shift (myCol, rowAlign, c);

    const Int localHeight = Length (height, colShift, r);
    const Int localWidth = Length (width, rowShift, c);
    const Int iLocalOffset = Length (i, Y.ColShift(), r);
    const Int jLocalOffset = Length (j, Y.RowShift(), c);

    for (Int t = 0; t < localWidth; ++t)
    {
	T *YCol = Y.Buffer (iLocalOffset, jLocalOffset + t);
	const T *XCol = &XBuffer[t * localHeight];
	for (Int s = 0; s < localHeight; ++s)
	    YCol[s] += XCol[s];
    }
    // Free the memory
    getVector_.clear();
}

// handle request for data, post a matching issend
template < typename T >
void AxpyInterface2<T>::HandleGlobalToLocalData ( Matrix<T>& Z )
{
    DEBUG_ONLY (CallStackEntry cse ("AxpyInterface::HandleGlobalToLocalData"))

    if ( !toBeAttachedForGet_ )
	LogicError("Local matrix cannot be updated");

    DistMatrix<T>& Y = *GlobalArrayGet_;
    const Grid& g = Y.Grid();
    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();
    const Int myRow = g.Row();
    const Int myCol = g.Col();

    Int matrix_index;

    mpi::Status status;

    if (mpi::IProbe (mpi::ANY_SOURCE, REQUEST_GET_TAG, g.VCComm (), status))
    {
	const Int source = status.MPI_SOURCE;

	// post receive for coordinates
	Int coord[2] = {-1, -1};
	mpi::TaggedRecv (coord, 2, source, 
		REQUEST_GET_TAG, g.VCComm());
    	Int i = coord[0]; Int j = coord[1];
 
	const Int colAlign = (Y.ColAlign() + i) % r;
	const Int rowAlign = (Y.RowAlign() + j) % c;

	const Int XLDim = Z.LDim();
	// local matrix width and height
	const Int height = Z.Height();
	const Int width = Z.Width();

	const Int colShift = Shift (myRow, colAlign, r);
	const Int rowShift = Shift (myCol, rowAlign, c);
	const Int localHeight = Length (height, colShift, r);
	const Int localWidth = Length (width, rowShift, c);

	const Int iLocalOffset = Length (i, Y.ColShift (), r);
	const Int jLocalOffset = Length (j, Y.RowShift (), c);

	const Int numEntries = localHeight * localWidth;

	DEBUG_ONLY (if (numEntries < Int (sizeof (T)))
		LogicError ("Count was too small");)

	T* XBuffer = Z.Buffer();
	const Int index =
	NextIndexMatrix (source,
		numEntries, 
		XBuffer,
		&matrix_index);

	DEBUG_ONLY (if
		(Int (matrices_[matrix_index].data_[source][index].size ()) !=
		 numEntries) LogicError ("Error in NextIndexMatrix");)
	
	T *replyBuffer = matrices_[matrix_index].data_[source][index].data ();
	for (Int t = 0; t < localWidth; ++t)
	{
	    T *sendCol = &replyBuffer[t * localHeight];
	    const T *XCol = Y.LockedBuffer (iLocalOffset, jLocalOffset + t);
	    MemCopy (sendCol, XCol, localHeight);
	}

	// Fire off non-blocking send
	mpi::TaggedISend (replyBuffer, numEntries, source, 
		DATA_GET_TAG, g.VCComm (), 
		matrices_[matrix_index].requests_[source][index]);
    }
}

// blocking update routines
template<typename T>
void AxpyInterface2<T>::Put( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Put"))

    if( i < 0 || j < 0 )
	LogicError("Submatrix offsets must be non-negative");
    if ( !toBeAttachedForPut_ )
	LogicError("Global matrix cannot be updated");

    DistMatrix<T>& Y = *GlobalArrayPut_;

    //do boundary checks
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
	LogicError("Submatrix out of bounds of global matrix");

    const Grid& g = Y.Grid();

    const Int XLDim = Z.LDim();

    const Int height = Z.Height();
    const Int width = Z.Width();

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();

    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();

    Int matrix_index, coord_index;

    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
 
    const Int colAlign = (Y.ColAlign() + i) % r;	    
    const Int rowAlign = (Y.RowAlign() + j) % c;

    const Int YLDim = Y.LDim();
	
    for( Int step=0; step<p; ++step )
    {
	const Int colShift = Shift( receivingRow, colAlign, r );
	const Int rowShift = Shift( receivingCol, rowAlign, c );

	const Int localHeight = Length( height, colShift, r );
	const Int localWidth = Length( width, rowShift, c );
	
	const Int numEntries = localHeight * localWidth;

	if( numEntries != 0  )
	{
	    const Int destination = receivingRow + r*receivingCol;

	    T* XBuffer = Z.Buffer();
	    // send data 
	    const Int dindex =
		NextIndexMatrix (destination,
			numEntries, 
			XBuffer,
			&matrix_index);

	    DEBUG_ONLY (if
		    (Int (matrices_[matrix_index].data_[destination][dindex].size ()) !=
		     numEntries) LogicError ("Error in NextIndexMatrix");)

	    T *sendBuffer = matrices_[matrix_index].data_[destination][dindex].data ();
	   
	    for( Int t=0; t<localWidth; ++t )
	    {
		T* thisSendCol = &sendBuffer[t*localHeight];
		T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
		for( Int s=0; s<localHeight; ++s )
		    thisSendCol[s] = thisXCol[colShift+s*r];
	    }

	    mpi::TaggedISend (sendBuffer, numEntries, destination, 
		    DATA_PUT_TAG, g.VCComm(), 
		    matrices_[matrix_index].requests_[destination][dindex]);

	    // send coordinates
	    const Int cindex =
		NextIndexCoord (i, j,
			destination,
			XBuffer,
			&coord_index);

	    Int *coord = reinterpret_cast<Int *>(coords_[coord_index].coord_[destination][cindex].data());
	    coord[0] = i; coord[1] = j;
	    mpi::TaggedISend (coord, 2, destination, 
		    COORD_PUT_TAG, g.VCComm(), 
		    coords_[coord_index].requests_[destination][cindex]);
	}

	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;
    }

    // progress my sends
    TestMatrix ( Z );
    TestCoord ( Z );

    // my receives
    std::vector<T> getVector;
    bool flag = true;
    
    while ( flag )
    {
	mpi::Status status;
        flag = mpi::IProbe( mpi::ANY_SOURCE, DATA_PUT_TAG, g.VCComm (), status);
	
	if ( flag )
	{
	    const Int source = status.MPI_SOURCE;
	    const Int count = mpi::GetCount <T> (status);

	    // post receive for coordinates
	    Int coord[2];
	    mpi::TaggedRecv (coord, 2, source, 
		    COORD_PUT_TAG, g.VCComm());
	    Int i = coord[0]; Int j = coord[1];

	    // post recv for data	
	    getVector.resize (count);
	    T *getBuffer = getVector.data ();

	    mpi::TaggedRecv (getBuffer, count, source, 
		    DATA_PUT_TAG, g.VCComm ());

	    // Update Y
	    const T *XBuffer = reinterpret_cast < const T * >(getBuffer);

	    const Int colAlign = (Y.ColAlign () + i) % r;
	    const Int rowAlign = (Y.RowAlign () + j) % c;

	    const Int colShift = Shift (g.Row(), colAlign, r);
	    const Int rowShift = Shift (g.Col(), rowAlign, c);

	    const Int localHeight = Length (height, colShift, r);
	    const Int localWidth = Length (width, rowShift, c);

	    const Int iLocalOffset = Length (i, Y.ColShift (), r);
	    const Int jLocalOffset = Length (j, Y.RowShift (), c);

	    for (Int t = 0; t < localWidth; ++t)
	    {
		T *YCol = Y.Buffer (iLocalOffset, jLocalOffset + t);
		const T *XCol = &XBuffer[t * localHeight];
		for (Int s = 0; s < localHeight; ++s)
		    YCol[s] = XCol[s];
	    }
	}
    }
   
    // wait for my sends
    WaitMatrix ( Z );
    
    getVector.clear();
}

template<typename T>
void AxpyInterface2<T>::Get( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Get"))
    // a call to Attach with a non-const DistMatrix must set
    // toBeAttachedForGet_ also, if not then it is assumed that
    // the DistMatrix isn't attached
    if ( !toBeAttachedForGet_ )
	LogicError ("Cannot perform this operation as matrix is not attached.");
    DistMatrix<T>& X = *GlobalArrayGet_;

    const Int height = Z.Height ();
    const Int width = Z.Width ();
    
    if (i + height > X.Height () || j + width > X.Width ())
	LogicError ("Invalid AxpyGlobalToLocal submatrix");

    const Grid & g = X.Grid ();
    const Int r = g.Height ();
    const Int c = g.Width ();
    const Int p = g.Size ();

    Int matrix_index;
     	
    const Int XLDim = Z.LDim();	
    const Int colAlign = (X.ColAlign() + i) % r;
    const Int rowAlign = (X.RowAlign() + j) % c;
    const Int iLocalOffset = Length (i, X.ColShift (), r);
    const Int jLocalOffset = Length (j, X.RowShift (), c);

    Int receivingRow = g.Row();
    Int receivingCol = g.Col();
   
    for( Int step=0; step<p; ++step )
    {
        const Int colShift = Shift( receivingRow, colAlign, r );
        const Int rowShift = Shift( receivingCol, rowAlign, c );
	
	const Int localHeight = Length (height, colShift, r);
	const Int localWidth = Length (width, rowShift, c);

	const Int numEntries = localHeight * localWidth;

	DEBUG_ONLY (if (numEntries < Int (sizeof (T)))
		LogicError ("Count was too small");)

	if( numEntries != 0 )
	{
	    const Int source = receivingRow + r*receivingCol;
	    T* XBuffer = Z.Buffer();
	
	    const Int index =
	    NextIndexMatrix (source,
		    numEntries, 
		    XBuffer,
		    &matrix_index);

	    DEBUG_ONLY (if
		(Int (matrices_[matrix_index].data_[source][index].size ()) !=
		 numEntries) LogicError ("Error in NextIndexMatrix");)

            T *replyBuffer = matrices_[matrix_index].data_[source][index].data ();
		
	    for (Int t = 0; t < localWidth; ++t)
	    {
		T *sendCol = &replyBuffer[t * localHeight];
		const T *XCol = X.LockedBuffer (iLocalOffset, jLocalOffset + t);
		MemCopy (sendCol, XCol, localHeight);
	    }
		
	    // Fire off non-blocking send
	    mpi::TaggedISend (replyBuffer, numEntries, source, 
		    DATA_GET_TAG, g.VCComm (), 
		    matrices_[matrix_index].requests_[source][index]);
	}
	
	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;
    }

    // progress my sends
    TestMatrix (Z);
 
    std::vector<T> recvVector_;
    // Receive all of the replies
    while ( 1 )
    {
	mpi::Status status;
	if (mpi::IProbe
		(mpi::ANY_SOURCE, DATA_GET_TAG, g.VCComm (), status))
	{
	    const Int source = status.MPI_SOURCE;
	    // Ensure that we have a recv buffer
	    const Int count = mpi::GetCount <T> (status);
	    recvVector_.resize (count);
	    T *recvBuffer = recvVector_.data ();

	    // Receive the data
	    mpi::TaggedRecv
		(recvBuffer, count, source, DATA_GET_TAG, g.VCComm ());

	    // Compute the local heights and offsets
	    const Int myRow = g.Row ();
	    const Int myCol = g.Col ();
	    
	    const Int colAlign = (X.ColAlign () + i) % r;
	    const Int rowAlign = (X.RowAlign () + j) % c;
	    
	    const Int colShift = Shift (myRow, colAlign, r);
	    const Int rowShift = Shift (myCol, rowAlign, c);

	    const Int localHeight = Length (height, colShift, r);
	    const Int localWidth = Length (width, rowShift, c);

	    // Unpack the local matrix
	    for (Int t = 0; t < localWidth; ++t)
	    {
		T *YCol = X.Buffer (0, rowShift + t * c);
		const T *XCol = &recvBuffer[t * localHeight];
		for (Int s = 0; s < localHeight; ++s)
		    YCol[colShift + s * r] = XCol[s];
	    }
	}
	recvVector_.clear();
    }
}

// accumulate = Update Y(i:i+height-1,j:j+width-1) += X,
// where X is height x width
template<typename T>
void AxpyInterface2<T>::Acc( Matrix<T>& Z, Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Acc"))

    if( i < 0 || j < 0 )
	LogicError("Submatrix offsets must be non-negative");
    if ( !toBeAttachedForPut_ )
	LogicError("Global matrix cannot be updated");

    DistMatrix<T>& Y = *GlobalArrayPut_;

    //do boundary checks
    if( i+Z.Height() > Y.Height() || j+Z.Width() > Y.Width() )
	LogicError("Submatrix out of bounds of global matrix");

    const Grid& g = Y.Grid();

    const Int XLDim = Z.LDim();

    const Int height = Z.Height();
    const Int width = Z.Width();

    const Int r = g.Height();
    const Int c = g.Width();
    const Int p = g.Size();

    const Int myProcessRow = g.Row();
    const Int myProcessCol = g.Col();

    Int matrix_index, coord_index;

    Int receivingRow = myProcessRow;
    Int receivingCol = myProcessCol;
 
    const Int colAlign = (Y.ColAlign() + i) % r;	    
    const Int rowAlign = (Y.RowAlign() + j) % c;

    const Int YLDim = Y.LDim();
	
    for( Int step=0; step<p; ++step )
    {
	const Int colShift = Shift( receivingRow, colAlign, r );
	const Int rowShift = Shift( receivingCol, rowAlign, c );

	const Int localHeight = Length( height, colShift, r );
	const Int localWidth = Length( width, rowShift, c );
	
	const Int numEntries = localHeight * localWidth;

	if( numEntries != 0  )
	{
	    const Int destination = receivingRow + r*receivingCol;

	    T* XBuffer = Z.Buffer();
	    // send data 
	    const Int dindex =
		NextIndexMatrix (destination,
			numEntries, 
			XBuffer,
			&matrix_index);

	    DEBUG_ONLY (if
		    (Int (matrices_[matrix_index].data_[destination][dindex].size ()) !=
		     numEntries) LogicError ("Error in NextIndexMatrix");)

	    T *sendBuffer = matrices_[matrix_index].data_[destination][dindex].data ();
	   
	    for( Int t=0; t<localWidth; ++t )
	    {
		T* thisSendCol = &sendBuffer[t*localHeight];
		T* thisXCol = &XBuffer[(rowShift+t*c)*XLDim];
		for( Int s=0; s<localHeight; ++s )
		    thisSendCol[s] = thisXCol[colShift+s*r];
	    }

	    mpi::TaggedISend (sendBuffer, numEntries, destination, 
		    DATA_ACC_TAG, g.VCComm(), 
		    matrices_[matrix_index].requests_[destination][dindex]);

	    // send coordinates
	    const Int cindex =
		NextIndexCoord (i, j,
			destination,
			XBuffer,
			&coord_index);

	    Int *coord = reinterpret_cast<Int *>(coords_[coord_index].coord_[destination][cindex].data());
	    coord[0] = i; coord[1] = j;
	    mpi::TaggedISend (coord, 2, destination, 
		    COORD_ACC_TAG, g.VCComm(), 
		    coords_[coord_index].requests_[destination][cindex]);
	}

	receivingRow = (receivingRow + 1) % r;
	if( receivingRow == 0 )
	    receivingCol = (receivingCol + 1) % c;
    }

    // progress my sends
    TestMatrix ( Z );
    TestCoord ( Z );

    // my receives
    std::vector<T> getVector;
    bool flag = true;
    
    while ( flag )
    {
	mpi::Status status;
        flag = mpi::IProbe( mpi::ANY_SOURCE, DATA_ACC_TAG, g.VCComm (), status);
	
	if ( flag )
	{
	    const Int source = status.MPI_SOURCE;
	    const Int count = mpi::GetCount <T> (status);

	    // post receive for coordinates
	    Int coord[2];
	    mpi::TaggedRecv (coord, 2, source, 
		    COORD_ACC_TAG, g.VCComm());
	    Int i = coord[0]; Int j = coord[1];

	    // post recv for data	
	    getVector.resize (count);
	    T *getBuffer = getVector.data ();

	    mpi::TaggedRecv (getBuffer, count, source, 
		    DATA_ACC_TAG, g.VCComm ());

	    // Update Y
	    const T *XBuffer = reinterpret_cast < const T * >(getBuffer);

	    const Int colAlign = (Y.ColAlign () + i) % r;
	    const Int rowAlign = (Y.RowAlign () + j) % c;

	    const Int colShift = Shift (g.Row(), colAlign, r);
	    const Int rowShift = Shift (g.Col(), rowAlign, c);

	    const Int localHeight = Length (height, colShift, r);
	    const Int localWidth = Length (width, rowShift, c);

	    const Int iLocalOffset = Length (i, Y.ColShift (), r);
	    const Int jLocalOffset = Length (j, Y.RowShift (), c);

	    for (Int t = 0; t < localWidth; ++t)
	    {
		T *YCol = Y.Buffer (iLocalOffset, jLocalOffset + t);
		const T *XCol = &XBuffer[t * localHeight];
		for (Int s = 0; s < localHeight; ++s)
		    YCol[s] += XCol[s];
	    }
	}
    }
   
    // wait for my sends
    WaitMatrix ( Z );
    
    getVector.clear();
}

// detach collectively 
template<typename T>
void AxpyInterface2<T>::Detach()
{
    DEBUG_ONLY(CallStackEntry cse("AxpyInterface2::Detach"))
    // destructor will call detach again...
    if (detached_)
	return;
    if( !attached_ )
	LogicError("Must attach before detaching.");

    const Grid& g = ( toBeAttachedForPut_ ?
	    GlobalArrayPut_->Grid() :
	    GlobalArrayGet_->Grid() );

    mpi::Barrier( g.VCComm() );

    attached_ 		= false;
    detached_ 		= true;
    toBeAttachedForPut_ = false;
    toBeAttachedForGet_ = false;

    GlobalArrayPut_ 	= 0;
    GlobalArrayGet_ 	= 0;

    matrices_.clear();
    coords_.clear();
}

template class AxpyInterface2<Int>;
template class AxpyInterface2<float>;
template class AxpyInterface2<double>;
template class AxpyInterface2<Complex<float>>;
template class AxpyInterface2<Complex<double>>;

} // namespace El
