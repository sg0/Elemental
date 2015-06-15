/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Public section
// ##############

// Constructors and destructors
// ============================

template<typename T>
AbstractDistMatrix<T>::AbstractDistMatrix( const El::Grid& grid, int root )
: viewType_(OWNER),
  height_(0), width_(0),
  matrix_(0,0,true),
  colConstrained_(false), rowConstrained_(false), rootConstrained_(false),
  colAlign_(0), rowAlign_(0),
  root_(root), grid_(&grid)
{ }

template<typename T>
AbstractDistMatrix<T>::AbstractDistMatrix( AbstractDistMatrix<T>&& A ) 
EL_NOEXCEPT
: viewType_(A.viewType_),
  height_(A.height_), width_(A.width_), 
  colConstrained_(A.colConstrained_), rowConstrained_(A.rowConstrained_),
  rootConstrained_(A.rootConstrained_),
  colAlign_(A.colAlign_), rowAlign_(A.rowAlign_),
  colShift_(A.colShift_), rowShift_(A.rowShift_), 
  root_(A.root_),
  grid_(A.grid_)
{ matrix_.ShallowSwap( A.matrix_ ); }

// Optional to override
// --------------------

template<typename T>
AbstractDistMatrix<T>::~AbstractDistMatrix() { }

// Assignment and reconfiguration
// ==============================

template<typename T>
AbstractDistMatrix<T>& 
AbstractDistMatrix<T>::operator=( AbstractDistMatrix<T>&& A )
{
    if( Viewing() || A.Viewing() )
    {
        Copy( A, *this );
    }
    else
    {
        matrix_.ShallowSwap( A.matrix_ );
        viewType_ = A.viewType_;
        height_ = A.height_;
        width_ = A.width_;
        colConstrained_ = A.colConstrained_;
        rowConstrained_ = A.rowConstrained_;
        rootConstrained_ = A.rootConstrained_;
        colAlign_ = A.colAlign_;
        rowAlign_ = A.rowAlign_;
        colShift_ = A.colShift_;
        rowShift_ = A.rowShift_;
        root_ = A.root_;
        grid_ = A.grid_;
    }
    return *this;
}

template<typename T>
const AbstractDistMatrix<T>&
AbstractDistMatrix<T>::operator=( const DistMultiVec<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::operator="))
    Copy( A, *this );
    return *this;
}

template<typename T>
void
AbstractDistMatrix<T>::Empty()
{
    matrix_.Empty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
    colAlign_ = 0;
    rowAlign_ = 0;
    colConstrained_ = false;
    rowConstrained_ = false;
    rootConstrained_ = false;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::EmptyData()
{
    matrix_.Empty_();
    viewType_ = OWNER;
    height_ = 0;
    width_ = 0;
}

template<typename T>
void
AbstractDistMatrix<T>::SetGrid( const El::Grid& grid )
{
    if( grid_ != &grid )
    {
        grid_ = &grid; 
        Empty();
    }
}

#if defined(EL_USE_WIN_ALLOC_FOR_RMA) && \
	!defined(EL_USE_WIN_CREATE_FOR_RMA)
template<typename T>
void
AbstractDistMatrix<T>::SetDim( Int height, Int width )
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::SetDim");
        AssertNotLocked();
        if( Viewing() && (height > height_ || width > width_) )
            LogicError("Tried to increase the size of a view");
    )
    height_ = height; 
    width_ = width;    
    
    matrix_.SetDim_ ( Length(height,ColShift(),ColStride()),
    	              Length(width,RowShift(),RowStride()) );
    //std::cout << "ADM:: height = " << Length(height,ColShift(),ColStride()) 
    //	<< " width = " << Length(width,RowShift(),RowStride()) << "\n";
}

template<typename T>
void
AbstractDistMatrix<T>::SetWindowBase( T* baseptr )
{
    DEBUG_ONLY( CallStackEntry cse("ADM::SetWindowBase"); )
    matrix_.SetWindowBase_( baseptr );
}
#endif

template<typename T>
void
AbstractDistMatrix<T>::Resize( Int height, Int width )
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::Resize");
        AssertNotLocked();
        if( Viewing() && (height > height_ || width > width_) )
            LogicError("Tried to increase the size of a view");
    )
    height_ = height; 
    width_ = width;
    if( Participating() )
        matrix_.Resize_
        ( Length(height,ColShift(),ColStride()),
          Length(width,RowShift(),RowStride()) );
}

template<typename T>
void
AbstractDistMatrix<T>::Resize( Int height, Int width, Int ldim )
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::Resize");
        AssertNotLocked();
        if( Viewing() && 
            (height > height_ || width > width_ || ldim > matrix_.LDim()) )
            LogicError("Tried to increase the size of a view");
    )
    height_ = height; 
    width_ = width;
    if( Participating() )
        matrix_.Resize_
        ( Length(height,ColShift(),ColStride()),
          Length(width,RowShift(),RowStride()), ldim );
}

template<typename T>
void
AbstractDistMatrix<T>::MakeConsistent( bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeConsistent"))

    const Int msgLength = 9;
    Int message[msgLength];
    if( CrossRank() == Root() )
    {
        message[0] = viewType_;
        message[1] = height_;
        message[2] = width_;
        message[3] = colConstrained_;
        message[4] = rowConstrained_;
        message[5] = rootConstrained_;
        message[6] = colAlign_;
        message[7] = rowAlign_;
        message[8] = root_;
    }

    const El::Grid& g = *grid_;
    if( !g.InGrid() && !includingViewers )
        LogicError("Non-participating process called MakeConsistent");
    if( g.InGrid() )
    {
        // TODO: Ensure roots are consistent within each cross communicator
        mpi::Broadcast( message, msgLength, Root(), CrossComm() );
    }
    if( includingViewers )
    {
        const Int vcRoot = g.VCToViewingMap(0);
        mpi::Broadcast( message, msgLength, vcRoot, g.ViewingComm() );
    }
    const ViewType newViewType    = static_cast<ViewType>(message[0]);
    const Int newHeight           = message[1]; 
    const Int newWidth            = message[2];
    const bool newConstrainedCol  = message[3];
    const bool newConstrainedRow  = message[4];
    const bool newConstrainedRoot = message[5];
    const Int newColAlign         = message[6];
    const Int newRowAlign         = message[7];
    const int root                = message[8];

    root_            = root;
    viewType_        = newViewType;
    colConstrained_  = newConstrainedCol;
    rowConstrained_  = newConstrainedRow;
    rootConstrained_ = newConstrainedRoot;
    colAlign_        = newColAlign;
    rowAlign_        = newRowAlign;

    SetShifts();
    Resize( newHeight, newWidth );
}

template<typename T>
void
AbstractDistMatrix<T>::MakeSizeConsistent( bool includingViewers )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeSizeConsistent"))

    const Int msgSize = 2;
    Int message[msgSize];
    if( CrossRank() == Root() )
    {
        message[0] = height_;
        message[1] = width_;
    }

    const El::Grid& g = *grid_;
    if( !g.InGrid() && !includingViewers )
        LogicError("Non-participating process called MakeSizeConsistent");
    if( g.InGrid() )
        mpi::Broadcast( message, msgSize, Root(), CrossComm() );
    if( includingViewers )
    {
        const Int vcRoot = g.VCToViewingMap(0);
        mpi::Broadcast( message, msgSize, vcRoot, g.ViewingComm() );
    }
    const Int newHeight = message[0]; 
    const Int newWidth  = message[1];
    Resize( newHeight, newWidth );
}

// Realignment
// -----------

template<typename T>
void
AbstractDistMatrix<T>::Align( int colAlign, int rowAlign, bool constrain )
{ 
    DEBUG_ONLY(CallStackEntry cse("ADM::Align"))
    const bool requireChange = colAlign_ != colAlign || rowAlign_ != rowAlign;
    DEBUG_ONLY(
        if( Viewing() && requireChange )
            LogicError("Tried to realign a view");
    )
    if( requireChange )
        Empty();
    if( constrain )
    {
        colConstrained_ = true;
        rowConstrained_ = true;
    }
    colAlign_ = colAlign;
    rowAlign_ = rowAlign;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::AlignCols( int colAlign, bool constrain )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ADM::AlignCols");
        if( Viewing() && colAlign_ != colAlign )
            LogicError("Tried to realign a view");
    )
    if( colAlign_ != colAlign )
        EmptyData();
    if( constrain )
        colConstrained_ = true;
    colAlign_ = colAlign;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::AlignRows( int rowAlign, bool constrain )
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ADM::AlignRows");
        if( Viewing() && rowAlign_ != rowAlign )
            LogicError("Tried to realign a view");
    )
    if( rowAlign_ != rowAlign )
        EmptyData();
    if( constrain )
        rowConstrained_ = true;
    rowAlign_ = rowAlign;
    SetShifts();
}

template<typename T>
void
AbstractDistMatrix<T>::FreeAlignments() 
{ 
    if( !Viewing() )
    {
        colConstrained_ = false;
        rowConstrained_ = false;
        rootConstrained_ = false;
    }
    else
        LogicError("Cannot free alignments of views");
}

template<typename T>
void
AbstractDistMatrix<T>::SetRoot( int root, bool constrain )
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::SetRoot");
        if( root < 0 || root >= mpi::Size(CrossComm()) )
            LogicError("Invalid root");
    )
    if( root != root_ )
        Empty();
    root_ = root;
    if( constrain )
        rootConstrained_ = true;
}

template<typename T>
void
AbstractDistMatrix<T>::AlignWith
( const El::DistData& data, bool constrain, bool allowMismatch )
{ 
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignWith"))
    AlignColsWith( data, constrain, allowMismatch );
    AlignRowsWith( data, constrain, allowMismatch );
}

template<typename T>
void
AbstractDistMatrix<T>::AlignColsWith
( const El::DistData& data, bool constrain, bool allowMismatch )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignColsWith"))
    SetGrid( *data.grid );
    SetRoot( data.root );
    if(      data.colDist == ColDist() || data.colDist == PartialColDist() )
        AlignCols( data.colAlign, constrain );
    else if( data.rowDist == ColDist() || data.rowDist == PartialColDist() )
        AlignCols( data.rowAlign, constrain );
    else if( data.colDist == PartialUnionColDist() )
        AlignCols( data.colAlign % ColStride(), constrain );
    else if( data.rowDist == PartialUnionColDist() )
        AlignCols( data.rowAlign % ColStride(), constrain );
    else if( ColDist()    != CollectedColDist() && 
             data.colDist != CollectedColDist() && 
             data.rowDist != CollectedColDist() && !allowMismatch )
        LogicError("Nonsensical alignment");
}

template<typename T>
void AbstractDistMatrix<T>::AlignRowsWith
( const El::DistData& data, bool constrain, bool allowMismatch )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignRowsWith"))
    SetGrid( *data.grid );
    SetRoot( data.root );
    if(      data.colDist == RowDist() || data.colDist == PartialRowDist() )
        AlignRows( data.colAlign, constrain );
    else if( data.rowDist == RowDist() || data.rowDist == PartialRowDist() )
        AlignRows( data.rowAlign, constrain );
    else if( data.colDist == PartialUnionRowDist() )
        AlignRows( data.colAlign % RowStride(), constrain );
    else if( data.rowDist == PartialUnionRowDist() )
        AlignRows( data.rowAlign % RowStride(), constrain );
    else if( RowDist()    != CollectedRowDist() && 
             data.colDist != CollectedRowDist() && 
             data.rowDist != CollectedRowDist() && !allowMismatch )
        LogicError("Nonsensical alignment");
}

template<typename T>
void
AbstractDistMatrix<T>::AlignAndResize
( int colAlign, int rowAlign, Int height, Int width, 
  bool force, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignAndResize"))
    if( !Viewing() )
    {
        if( force || !ColConstrained() )
        {
            colAlign_ = colAlign;
            SetColShift(); 
        }
        if( force || !RowConstrained() )
        {
            rowAlign_ = rowAlign;
            SetRowShift();
        }
    }
    if( constrain )
    {
        colConstrained_ = true;
        rowConstrained_ = true;
    }
    if( force && (colAlign_ != colAlign || rowAlign_ != rowAlign) )
        LogicError("Could not set alignments"); 
    Resize( height, width );
}

template<typename T>
void
AbstractDistMatrix<T>::AlignColsAndResize
( int colAlign, Int height, Int width, bool force, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignColsAndResize"))
    if( !Viewing() && (force || !ColConstrained()) )
    {
        colAlign_ = colAlign;
        SetColShift(); 
    }
    if( constrain )
        colConstrained_ = true;
    if( force && colAlign_ != colAlign )
        LogicError("Could not set col alignment");
    Resize( height, width );
}

template<typename T>
void
AbstractDistMatrix<T>::AlignRowsAndResize
( int rowAlign, Int height, Int width, bool force, bool constrain )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::AlignRowsAndResize"))
    if( !Viewing() && (force || !RowConstrained()) )
    {
        rowAlign_ = rowAlign;
        SetRowShift(); 
    }
    if( constrain )
        rowConstrained_ = true;
    if( force && rowAlign_ != rowAlign )
        LogicError("Could not set row alignment");
    Resize( height, width );
}

// Buffer attachment
// -----------------

template<typename T>
void
AbstractDistMatrix<T>::Attach
( Int height, Int width, const El::Grid& g, 
  int colAlign, int rowAlign, T* buffer, Int ldim, int root )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Attach"))
    Empty();

    grid_ = &g;
    root_ = root;
    height_ = height;
    width_ = width;
    colAlign_ = colAlign;
    rowAlign_ = rowAlign;
    colConstrained_ = true;
    rowConstrained_ = true;
    rootConstrained_ = true;
    viewType_ = VIEW;
    SetShifts();
    if( Participating() )
    {
        Int localHeight = Length(height,colShift_,ColStride());
        Int localWidth = Length(width,rowShift_,RowStride());
        matrix_.Attach_( localHeight, localWidth, buffer, ldim );
    }
}

template<typename T>
void
AbstractDistMatrix<T>::Attach
( Int height, Int width, const El::Grid& g,
  int colAlign, int rowAlign, El::Matrix<T>& A, int root )
{
    // TODO: Assert that the local dimensions are correct
    Attach( height, width, g, colAlign, rowAlign, A.Buffer(), A.LDim(), root );
}

template<typename T>
void
AbstractDistMatrix<T>::Attach( const El::Grid& g, El::Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Attach"))
    if( g.Size() != 1 )
        LogicError("Assumed a grid size of one");
    Attach( A.Height(), A.Width(), g, 0, 0, A.Buffer(), A.LDim() );
}

template<typename T>
void
AbstractDistMatrix<T>::LockedAttach
( Int height, Int width, const El::Grid& g, 
  int colAlign, int rowAlign, const T* buffer, Int ldim, int root )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::LockedAttach"))
    Empty();

    grid_ = &g;
    root_ = root;
    height_ = height;
    width_ = width;
    colAlign_ = colAlign;
    rowAlign_ = rowAlign;
    colConstrained_ = true;
    rowConstrained_ = true;
    rootConstrained_ = true;
    viewType_ = LOCKED_VIEW;
    SetShifts();
    if( Participating() )
    {
        Int localHeight = Length(height,colShift_,ColStride());
        Int localWidth = Length(width,rowShift_,RowStride());
        matrix_.LockedAttach_( localHeight, localWidth, buffer, ldim );
    }
}

template<typename T>
void
AbstractDistMatrix<T>::LockedAttach
( Int height, Int width, const El::Grid& g, 
  int colAlign, int rowAlign, const El::Matrix<T>& A, int root )
{
    // TODO: Assert that the local dimensions are correct
    LockedAttach
    ( height, width, g, colAlign, rowAlign, A.LockedBuffer(), A.LDim(), root );
}

template<typename T>
void
AbstractDistMatrix<T>::LockedAttach( const El::Grid& g, const El::Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::LockedAttach"))
    if( g.Size() != 1 )
        LogicError("Assumed a grid size of one");
    LockedAttach( A.Height(), A.Width(), g, 0, 0, A.LockedBuffer(), A.LDim() );
}

// Basic queries
// =============

// Global matrix information
// -------------------------

template<typename T>
Int AbstractDistMatrix<T>::Height() const { return height_; }
template<typename T>
Int AbstractDistMatrix<T>::Width() const { return width_; }

template<typename T>
Int AbstractDistMatrix<T>::DiagonalLength( Int offset ) const
{ return El::DiagonalLength(height_,width_,offset); }

template<typename T>
bool AbstractDistMatrix<T>::Viewing() const { return IsViewing( viewType_ ); }
template<typename T>
bool AbstractDistMatrix<T>::Locked() const { return IsLocked( viewType_ ); }

// Local matrix information
// ------------------------

template<typename T>
Int AbstractDistMatrix<T>::LocalHeight() const { return matrix_.Height(); }
template<typename T>
Int AbstractDistMatrix<T>::LocalWidth() const { return matrix_.Width(); }
template<typename T>
Int AbstractDistMatrix<T>::LDim() const { return matrix_.LDim(); }

template<typename T>
El::Matrix<T>& 
AbstractDistMatrix<T>::Matrix() { return matrix_; }
template<typename T>
const El::Matrix<T>& 
AbstractDistMatrix<T>::LockedMatrix() const { return matrix_; }

template<typename T>
size_t
AbstractDistMatrix<T>::AllocatedMemory() const { return matrix_.MemorySize(); }

template<typename T>
T*
AbstractDistMatrix<T>::Buffer() { return matrix_.Buffer(); }

template<typename T>
T*
AbstractDistMatrix<T>::Buffer( Int iLoc, Int jLoc )
{ return matrix_.Buffer(iLoc,jLoc); }

template<typename T>
const T*
AbstractDistMatrix<T>::LockedBuffer() const
{ return matrix_.LockedBuffer(); }

template<typename T>
const T*
AbstractDistMatrix<T>::LockedBuffer( Int iLoc, Int jLoc ) const
{ return matrix_.LockedBuffer(iLoc,jLoc); }

// Distribution information
// ------------------------

template<typename T>
const El::Grid& AbstractDistMatrix<T>::Grid() const { return *grid_; }

template<typename T>
int AbstractDistMatrix<T>::ColAlign() const { return colAlign_; }
template<typename T>
int AbstractDistMatrix<T>::RowAlign() const { return rowAlign_; }

template<typename T>
int AbstractDistMatrix<T>::ColShift() const { return colShift_; }
template<typename T>
int AbstractDistMatrix<T>::RowShift() const { return rowShift_; }

template<typename T>
bool AbstractDistMatrix<T>::ColConstrained() const { return colConstrained_; }
template<typename T>
bool AbstractDistMatrix<T>::RowConstrained() const { return rowConstrained_; }
template<typename T>
bool AbstractDistMatrix<T>::RootConstrained() const { return rootConstrained_; }

template<typename T>
bool AbstractDistMatrix<T>::Participating() const
{ return grid_->InGrid() && (CrossRank()==root_); }

template<typename T>
int AbstractDistMatrix<T>::RowOwner( Int i ) const
{ return int((i+ColAlign()) % ColStride()); }
template<typename T>
int AbstractDistMatrix<T>::ColOwner( Int j ) const
{ return int((j+RowAlign()) % RowStride()); }
template<typename T>
int AbstractDistMatrix<T>::Owner( Int i, Int j ) const
{ return RowOwner(i)+ColOwner(j)*ColStride(); }

template<typename T>
Int AbstractDistMatrix<T>::LocalRow( Int i ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("ADM::LocalRow");
        if( !IsLocalRow(i) )
            LogicError("Requested local index of non-local row");
    )
    return LocalRowOffset(i);
}

template<typename T>
Int AbstractDistMatrix<T>::LocalCol( Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::LocalCol");
        if( !IsLocalCol(j) )
            LogicError("Requested local index of non-local column");
    )
    return LocalColOffset(j);
}

template<typename T>
Int AbstractDistMatrix<T>::LocalRowOffset( Int i ) const
{ return Length_(i,ColShift(),ColStride()); }
template<typename T>
Int AbstractDistMatrix<T>::LocalColOffset( Int j ) const
{ return Length_(j,RowShift(),RowStride()); }

template<typename T>
Int AbstractDistMatrix<T>::GlobalRow( Int iLoc ) const
{ return ColShift() + iLoc*ColStride(); }
template<typename T>
Int AbstractDistMatrix<T>::GlobalCol( Int jLoc ) const
{ return RowShift() + jLoc*RowStride(); }

template<typename T>
bool AbstractDistMatrix<T>::IsLocalRow( Int i ) const
{ return Participating() && RowOwner(i) == ColRank(); }
template<typename T>
bool AbstractDistMatrix<T>::IsLocalCol( Int j ) const
{ return Participating() && ColOwner(j) == RowRank(); }
template<typename T>
bool AbstractDistMatrix<T>::IsLocal( Int i, Int j ) const
{ return IsLocalRow(i) && IsLocalCol(j); }

template<typename T>
mpi::Comm AbstractDistMatrix<T>::PartialColComm() const { return ColComm(); }
template<typename T>
mpi::Comm AbstractDistMatrix<T>::PartialRowComm() const { return RowComm(); }

template<typename T>
mpi::Comm AbstractDistMatrix<T>::PartialUnionColComm() const
{ return mpi::COMM_SELF; }
template<typename T>
mpi::Comm AbstractDistMatrix<T>::PartialUnionRowComm() const
{ return mpi::COMM_SELF; }

template<typename T>
int AbstractDistMatrix<T>::PartialColStride() const { return ColStride(); }
template<typename T>
int AbstractDistMatrix<T>::PartialRowStride() const { return RowStride(); }

template<typename T>
int AbstractDistMatrix<T>::PartialUnionColStride() const { return 1; }
template<typename T>
int AbstractDistMatrix<T>::PartialUnionRowStride() const { return 1; }

template<typename T>
int AbstractDistMatrix<T>::ColRank() const { return mpi::Rank(ColComm()); }
template<typename T>
int AbstractDistMatrix<T>::RowRank() const { return mpi::Rank(RowComm()); }

template<typename T>
int AbstractDistMatrix<T>::PartialColRank() const
{ return mpi::Rank(PartialColComm()); }
template<typename T>
int AbstractDistMatrix<T>::PartialRowRank() const
{ return mpi::Rank(PartialRowComm()); }

template<typename T>
int AbstractDistMatrix<T>::PartialUnionColRank() const
{ return mpi::Rank(PartialUnionColComm()); }
template<typename T>
int AbstractDistMatrix<T>::PartialUnionRowRank() const
{ return mpi::Rank(PartialUnionRowComm()); }

template<typename T>
int AbstractDistMatrix<T>::DistRank() const
{ return mpi::Rank(DistComm()); }
template<typename T>
int AbstractDistMatrix<T>::CrossRank() const
{ return mpi::Rank(CrossComm()); }
template<typename T>
int AbstractDistMatrix<T>::RedundantRank() const
{ return mpi::Rank(RedundantComm()); }

template<typename T>
int AbstractDistMatrix<T>::Root() const { return root_; }

// Single-entry manipulation
// =========================

// Global entry manipulation
// -------------------------

template<typename T>
T
AbstractDistMatrix<T>::Get( Int i, Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::Get");
        if( !grid_->InGrid() )
            LogicError("Get should only be called in-grid");
    )
    T value;
    if( CrossRank() == Root() )
    {
        const int owner = Owner( i, j );
        if( owner == DistRank() )
            value = GetLocal( LocalRow(i), LocalCol(j) );
        mpi::Broadcast( value, owner, DistComm() );
    }
    mpi::Broadcast( value, Root(), CrossComm() ); 
    return value;
}

template<typename T>
Base<T>
AbstractDistMatrix<T>::GetRealPart( Int i, Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::GetRealPart");
        if( !grid_->InGrid() )
            LogicError("Get should only be called in-grid");
    )
    Base<T> value;
    if( CrossRank() == Root() )
    {
        const int owner = Owner( i, j );
        if( owner == DistRank() )
            value = GetLocalRealPart( LocalRow(i), LocalCol(j) );
        mpi::Broadcast( value, owner, DistComm() );
    }
    mpi::Broadcast( value, Root(), CrossComm() );
    return value;
}

template<typename T>
Base<T>
AbstractDistMatrix<T>::GetImagPart( Int i, Int j ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("ADM::GetImagPart");
        if( !grid_->InGrid() )
            LogicError("Get should only be called in-grid");
    )
    Base<T> value;
    if( IsComplex<T>::val )
    {
        if( CrossRank() == Root() )
        {
            const int owner = Owner( i, j );
            if( owner == DistRank() )
                value = GetLocalRealPart( LocalRow(i), LocalCol(j) );
            mpi::Broadcast( value, owner, DistComm() );
        }
        mpi::Broadcast( value, Root(), CrossComm() );
    }
    else
        value = 0;
    return value;
}

template<typename T>
void
AbstractDistMatrix<T>::Set( Int i, Int j, T value )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Set"))
    if( IsLocal(i,j) )
        SetLocal( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::SetRealPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetRealPart"))
    if( IsLocal(i,j) )
        SetLocalRealPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::SetImagPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::SetImagPart"))
    if( IsLocal(i,j) )
        SetLocalImagPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::Update( Int i, Int j, T value )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Update"))
    if( IsLocal(i,j) )
        UpdateLocal( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::UpdateRealPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateRealPart"))
    if( IsLocal(i,j) )
        UpdateLocalRealPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::UpdateImagPart( Int i, Int j, Base<T> value )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::UpdateImagPart"))
    if( IsLocal(i,j) )
        UpdateLocalImagPart( LocalRow(i), LocalCol(j), value );
}

template<typename T>
void
AbstractDistMatrix<T>::MakeReal( Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::MakeReal"))
    if( IsLocal(i,j) )
        MakeLocalReal( LocalRow(i), LocalCol(j) );
}

template<typename T>
void
AbstractDistMatrix<T>::Conjugate( Int i, Int j )
{
    DEBUG_ONLY(CallStackEntry cse("ADM::Conjugate"))
    if( IsLocal(i,j) )
        ConjugateLocal( LocalRow(i), LocalCol(j) );
}

// Local entry manipulation
// ------------------------

template<typename T>
T
AbstractDistMatrix<T>::GetLocal( Int i, Int j ) const
{ return matrix_.Get(i,j); }

template<typename T>
Base<T>
AbstractDistMatrix<T>::GetLocalRealPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetRealPart(iLoc,jLoc); }

template<typename T>
Base<T>
AbstractDistMatrix<T>::GetLocalImagPart( Int iLoc, Int jLoc ) const
{ return matrix_.GetImagPart(iLoc,jLoc); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Set(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocalRealPart( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.SetRealPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::SetLocalImagPart( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.SetImagPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocal( Int iLoc, Int jLoc, T alpha )
{ matrix_.Update(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocalRealPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.UpdateRealPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::UpdateLocalImagPart
( Int iLoc, Int jLoc, Base<T> alpha )
{ matrix_.UpdateImagPart(iLoc,jLoc,alpha); }

template<typename T>
void
AbstractDistMatrix<T>::MakeLocalReal( Int iLoc, Int jLoc )
{ matrix_.MakeReal( iLoc, jLoc ); }

template<typename T>
void
AbstractDistMatrix<T>::ConjugateLocal( Int iLoc, Int jLoc )
{ matrix_.Conjugate( iLoc, jLoc ); }

// Diagonal manipulation
// =====================
template<typename T>
bool AbstractDistMatrix<T>::DiagonalAlignedWith
( const El::DistData& d, Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::DiagonalAlignedWith"))
    if( Grid() != *d.grid )
        return false;

    const Int diagRoot = DiagonalRoot(offset);
    if( diagRoot != d.root )
        return false;

    const int diagAlign = DiagonalAlign(offset);
    const Dist UDiag = DiagCol( ColDist(), RowDist() ); 
    const Dist VDiag = DiagRow( ColDist(), RowDist() );
    if( d.colDist == UDiag && d.rowDist == VDiag )
        return d.colAlign == diagAlign;
    else if( d.colDist == VDiag && d.rowDist == UDiag )
        return d.rowAlign == diagAlign;
    else
        return false;
}

template<typename T>
int AbstractDistMatrix<T>::DiagonalRoot( Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::DiagonalRoot"))
    const El::Grid& grid = Grid();

    if( ColDist() == MC && RowDist() == MR )
    {
        // Result is an [MD,* ] or [* ,MD]
        int owner;
        if( offset >= 0 )
        {
            const int procRow = ColAlign();
            const int procCol = (RowAlign()+offset) % RowStride();
            owner = procRow + ColStride()*procCol;
        }
        else
        {
            const int procRow = (ColAlign()-offset) % ColStride();
            const int procCol = RowAlign();
            owner = procRow + ColStride()*procCol;
        }
        return grid.DiagPath(owner);
    }
    else if( ColDist() == MR && RowDist() == MC )
    {
        // Result is an [MD,* ] or [* ,MD]
        int owner;
        if( offset >= 0 )
        {
            const int procCol = ColAlign();
            const int procRow = (RowAlign()+offset) % RowStride();
            owner = procRow + ColStride()*procCol;
        }
        else
        {
            const int procCol = (ColAlign()-offset) % ColStride();
            const int procRow = RowAlign();
            owner = procRow + ColStride()*procCol;
        }
        return grid.DiagPath(owner);
    }
    else
        return Root();
}

template<typename T>
int AbstractDistMatrix<T>::DiagonalAlign( Int offset ) const
{
    DEBUG_ONLY(CallStackEntry cse("ADM::DiagonalAlign"))
    const El::Grid& grid = Grid();

    if( ColDist() == MC && RowDist() == MR )
    {
        // Result is an [MD,* ] or [* ,MD]
        int owner;
        if( offset >= 0 )
        {
            const int procRow = ColAlign();
            const int procCol = (RowAlign()+offset) % RowStride();
            owner = procRow + ColStride()*procCol;
        }
        else
        {
            const int procRow = (ColAlign()-offset) % ColStride();
            const int procCol = RowAlign();
            owner = procRow + ColStride()*procCol;
        }
        return grid.DiagPathRank(owner);
    }
    else if( ColDist() == MR && RowDist() == MC )
    {
        // Result is an [MD,* ] or [* ,MD]
        int owner;
        if( offset >= 0 )
        {
            const int procCol = ColAlign();
            const int procRow = (RowAlign()+offset) % RowStride();
            owner = procRow + ColStride()*procCol;
        }
        else
        {
            const int procCol = (ColAlign()-offset) % ColStride();
            const int procRow = RowAlign();
            owner = procRow + ColStride()*procCol;
        }
        return grid.DiagPathRank(owner);
    }
    else if( ColDist() == STAR )
    {
        // Result is a [V,* ] or [* ,V]
        if( offset >= 0 )
            return (RowAlign()+offset) % RowStride();
        else
            return RowAlign();
    }
    else
    {
        // Result is [U,V] or [V,U], where V is either STAR or CIRC
        if( offset >= 0 )
            return ColAlign();
        else
            return (ColAlign()-offset) % ColStride();
    }
}

// Assertions
// ==========

template<typename T>
void
AbstractDistMatrix<T>::ComplainIfReal() const
{ 
    if( !IsComplex<T>::val )
        LogicError("Called complex-only routine with real data");
}

template<typename T>
void
AbstractDistMatrix<T>::AssertNotLocked() const
{
    if( Locked() )
        LogicError("Assertion that matrix not be a locked view failed");
}

template<typename T>
void
AbstractDistMatrix<T>::AssertNotStoringData() const
{
    if( matrix_.MemorySize() > 0 )
        LogicError("Assertion that matrix not be storing data failed");
}

template<typename T>
void
AbstractDistMatrix<T>::AssertValidEntry( Int i, Int j ) const
{
    if( i < 0 || i >= Height() || j < 0 || j >= Width() )
        LogicError
        ("Entry (",i,",",j,") is out of bounds of ",Height(),
         " x ",Width()," matrix");
}

template<typename T>
void
AbstractDistMatrix<T>::AssertValidSubmatrix
( Int i, Int j, Int height, Int width ) const
{
    if( i < 0 || j < 0 )
        LogicError("Indices of submatrix were negative");
    if( height < 0 || width < 0 )
        LogicError("Dimensions of submatrix were negative");
    if( (i+height) > Height() || (j+width) > Width() )
        LogicError
        ("Submatrix is out of bounds: accessing up to (",i+height-1,
         ",",j+width-1,") of ",Height()," x ",Width()," matrix");
}

template<typename T> 
void
AbstractDistMatrix<T>::AssertSameSize( Int height, Int width ) const
{
    if( Height() != height || Width() != width )
        LogicError("Assertion that matrices be the same size failed");
}

// Private section
// ###############

// Exchange metadata with another matrix
// =====================================

template<typename T>
void 
AbstractDistMatrix<T>::ShallowSwap( AbstractDistMatrix<T>& A )
{
    matrix_.ShallowSwap( A.matrix_ );
    std::swap( viewType_, A.viewType_ );
    std::swap( height_ , A.height_ );
    std::swap( width_, A.width_ );
    std::swap( colConstrained_, A.colConstrained_ );
    std::swap( rowConstrained_, A.rowConstrained_ );
    std::swap( rootConstrained_, A.rootConstrained_ );
    std::swap( colAlign_, A.colAlign_ );
    std::swap( rowAlign_, A.rowAlign_ );
    std::swap( colShift_, A.colShift_ );
    std::swap( rowShift_, A.rowShift_ );
    std::swap( root_, A.root_ );
    std::swap( grid_, A.grid_ );
}

// Modify the distribution metadata
// ================================

template<typename T>
void
AbstractDistMatrix<T>::SetShifts()
{
    if( Participating() )
    {
        colShift_ = Shift(ColRank(),colAlign_,ColStride());
        rowShift_ = Shift(RowRank(),rowAlign_,RowStride());
    }
    else
    {
        colShift_ = 0;
        rowShift_ = 0;
    }
}

template<typename T>
void
AbstractDistMatrix<T>::SetColShift()
{
    if( Participating() )
        colShift_ = Shift(ColRank(),colAlign_,ColStride());
    else
        colShift_ = 0;
}

template<typename T>
void
AbstractDistMatrix<T>::SetRowShift()
{
    if( Participating() )
        rowShift_ = Shift(RowRank(),rowAlign_,RowStride());
    else
        rowShift_ = 0;
}

// Outside of class
// ----------------

template<typename T> 
void
AssertConforming1x2
( const AbstractDistMatrix<T>& AL, const AbstractDistMatrix<T>& AR )
{
    if( AL.Height() != AR.Height() )    
        LogicError
        ("1x2 not conformant:\n",
         DimsString(AL,"Left"),"\n",DimsString(AR,"Right"));
    if( AL.ColAlign() != AR.ColAlign() )
        LogicError("1x2 is misaligned");
}

template<typename T> 
void
AssertConforming2x1
( const AbstractDistMatrix<T>& AT, const AbstractDistMatrix<T>& AB )
{
    if( AT.Width() != AB.Width() )
        LogicError
        ("2x1 is not conformant:\n",
         DimsString(AT,"Top"),"\n",DimsString(AB,"Bottom"));
    if( AT.RowAlign() != AB.RowAlign() )
        LogicError("2x1 is not aligned");
}

template<typename T> 
void
AssertConforming2x2
( const AbstractDistMatrix<T>& ATL, const AbstractDistMatrix<T>& ATR, 
  const AbstractDistMatrix<T>& ABL, const AbstractDistMatrix<T>& ABR ) 
{
    if( ATL.Width() != ABL.Width() || ATR.Width() != ABR.Width() ||
        ATL.Height() != ATR.Height() || ABL.Height() != ABR.Height() )
        LogicError
        ("2x2 is not conformant:\n",
         DimsString(ATL,"TL"),"\n",DimsString(ATR,"TR"),"\n",
         DimsString(ABL,"BL"),"\n",DimsString(ABR,"BR"));
    if( ATL.ColAlign() != ATR.ColAlign() ||
        ABL.ColAlign() != ABR.ColAlign() ||
        ATL.RowAlign() != ABL.RowAlign() ||
        ATR.RowAlign() != ABR.RowAlign() )
        LogicError("2x2 set of matrices must aligned to combine");
}

// Instantiations for {Int,Real,Complex<Real>} for each Real in {float,double}
// ###########################################################################

#ifndef EL_RELEASE
 #define PROTO(T) \
  template class AbstractDistMatrix<T>;\
  template void AssertConforming1x2\
  ( const AbstractDistMatrix<T>& AL, const AbstractDistMatrix<T>& AR );\
  template void AssertConforming2x1\
  ( const AbstractDistMatrix<T>& AT, const AbstractDistMatrix<T>& AB );\
  template void AssertConforming2x2\
  ( const AbstractDistMatrix<T>& ATL, const AbstractDistMatrix<T>& ATR,\
    const AbstractDistMatrix<T>& ABL, const AbstractDistMatrix<T>& ABR );
#else
 #define PROTO(T) template class AbstractDistMatrix<T>;
#endif

#include "El/macros/Instantiate.h"

} // namespace El
