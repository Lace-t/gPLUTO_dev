////////////////////////////////////////////////////////////////////////
/*!
  \file
  \brief Array class definition

  \author A. Mignone (andrea.mignone@unito.it)
  \date   March 18, 2026
*/
////////////////////////////////////////////////////////////////////////

// --------------------------------------------------------
// 1. Template class for 1D arrays 
// --------------------------------------------------------

#ifndef CHECK_ARY_MEM_BOUNDS
  #define LOC_CHECK_ARY_MEM_BOUNDS
  #define CHECK_ARY_MEM_BOUNDS  NO
#endif

template<typename T, bool MANAGED=false>
class Array1DTemplated {
public:
  
  // Default constructor initializes an empty array
  #pragma omp declare target
  #pragma acc routine seq
  Array1DTemplated() : stride0_(-1), data_(nullptr) {}
  #pragma omp end declare target
  
  // Constructor that allocates an array of size n0 (only if MANAGED is false)
  #pragma omp declare target
  #pragma acc routine seq
  Array1DTemplated(int n0) {
    static_assert(MANAGED==false);
    create(n0);
  }
  #pragma omp end declare target

  // Constructor that assigns an existing data pointer (only if MANAGED is true)
  #pragma omp declare target
  #pragma acc routine seq
  Array1DTemplated(T* data, int n0) {
    static_assert(MANAGED==true);
    create(data,n0);
  }
  #pragma omp end declare target
  // Deleted assignment operator to prevent accidental copying
  #pragma omp declare target  
  #pragma acc routine seq
  Array1DTemplated operator=(const Array1DTemplated&) = delete;
  #pragma omp end declare target
  // Deleted copy constructor to prevent accidental copying
  #pragma omp declare target
  #pragma acc routine seq
  Array1DTemplated(const Array1DTemplated &) = delete;
  #pragma omp end declare target
  // Destructor deallocates memory if not managed
  #pragma omp declare target
  #pragma acc routine seq
  ~Array1DTemplated() {
    if constexpr (!MANAGED)
      destroy();
  }
  #pragma omp end declare target
  // Checks if the array is allocated
  bool isAllocated() {
    return data_ != nullptr;
  }
  
  // Allocates memory for the array and initializes metadata
  void create(int n0) {
    n0_=n0;
    stride0_ = 1;
    total_size_ = n0_;
    data_ = new T[total_size_];
    g_usedMemory += n0_*sizeof(T);
#if DEBUG == YES
printf("1D memory: +%f = (number of elements %d * size %d)/1.e+6\n"
       ,((float)(n0_*sizeof(T))/1.e6), n0_*n1_*n2_, sizeof(T));
#endif
  }

  // Assigns an existing data pointer without allocating new memory
  #pragma omp declare target
  #pragma acc routine seq
  void create(T* data, int n0) {
    n0_=n0;
    stride0_ = 1;
    total_size_ = n0_;
    data_ = data;
  }
  #pragma omp end declare target
  // Overloaded function call operator for accessing elements with optional bounds checking
  #pragma omp declare target
  #pragma acc routine seq
  T &operator () (int i) const {
    int indx = i*stride0_;
    // printf("indx: %d\n", indx);
    #if CHECK_ARY_MEM_BOUNDS
    if( i >= n0_ || i < 0 ){
      printf("! Array1DTemplated. First index out of bound\n");
      printf("  Accessing index %d for n0_ size %d\n", i, n0_);
      printf("  Abort()\n");
      abort();
    }
    if (indx < 0 || indx >= total_size_ ){
      printf("! Array1DTemplated. Index out of bound.\n");
      printf("  indx: %d\n", indx);
      printf("  total_size_: %d\n", total_size_);
      printf("  stride0_: %d\n", stride0_);
      printf("  Error accessing array at %d\n",i);
      printf("  Abort()\n");
      abort();
    }
    #endif
    return data_[indx];
  }
  #pragma omp end declare target
  // Deallocates memory and resets metadata fields
  #pragma omp declare target
  #pragma acc routine seq
  void destroy() {
    // #pragma acc routine seq
    delete [] data_;
    data_=nullptr;
    stride0_=-1;
    n0_=0;
  }
  #pragma omp end declare target
protected:
  T  *data_; // Pointer to array data
  int n0_; // Number of elements
  int stride0_; // Stride for indexing
  int total_size_; // Total allocated size
  char fill_[12]; // Manual padding for memory alignment
};

// Type aliases for commonly used array types
using Ary1D=Array1DTemplated<double,true>;
using Ary1DInt=Array1DTemplated<int,true>;
using Ary1Duint32_t=Array1DTemplated<uint32_t,true>;
using AryNotManaged1D=Array1DTemplated<double,false>;
using Ary1DNotManagedInt=Array1DTemplated<int,false>;

// --------------------------------------------------------
// 2. Template class for 2D arrays 
// --------------------------------------------------------

template<typename T, bool MANAGED=false>
class Array2DTemplated {
public:
  
  // Default constructor initializes an empty array
  #pragma omp declare target
  #pragma acc routine seq
  Array2DTemplated() : stride0_(-1), stride1_(-1), data_(nullptr) {}
  #pragma omp end declare target
  // Constructor that allocates an array of size n0*n1 (only if MANAGED is false)
  #pragma omp declare target
  #pragma acc routine seq
  Array2DTemplated(int n0, int n1) {
    static_assert(MANAGED==false);
    create(n0,n1);
  }
  #pragma omp end declare target
  // Constructor that assigns an existing data pointer (only if MANAGED is true)
  #pragma omp declare target
  #pragma acc routine seq
  Array2DTemplated(T* data, int n0, int n1) {
    static_assert(MANAGED==true);
    create(data,n0,n1);
  }
  #pragma omp end declare target
  // Deleted assignment operator to prevent accidental copying
  #pragma omp declare target
  #pragma acc routine seq
  Array2DTemplated operator=(const Array2DTemplated&) = delete;
  #pragma omp end declare target
  // Deleted copy constructor to prevent accidental copying
  #pragma omp declare target
  #pragma acc routine seq
  Array2DTemplated(const Array2DTemplated &) = delete;
  #pragma omp end declare target
  // Destructor deallocates memory if not managed
  #pragma omp declare target
  #pragma acc routine seq
  ~Array2DTemplated() {
    if constexpr (!MANAGED)
      destroy();
  }
  
  // Checks if the array is allocated
  bool isAllocated() {
    return data_ != nullptr;
  }

  #pragma omp end declare target
  // Allocates memory for the array and initializes metadata
  void create(int n0, int n1) {
    n0_=n0;n1_=n1;
    stride0_ = 1;
    stride1_ = n0_;
    total_size_ = n0_*n1_;
    data_ = new T[total_size_];
    g_usedMemory += n0_*n1_*sizeof(T);
#if DEBUG == YES
printf("2D memory: +%f = (number of elements %d * size %d)/1.e+6\n"
       ,((float)(n0_*n1_*sizeof(T))/1.e6), n0_*n1_*n2_, sizeof(T));
#endif
  }

  // Assigns an existing data pointer without allocating new memory
  #pragma omp declare target
  #pragma acc routine seq
  void create(T* data, int n0, int n1) {
    n0_=n0;n1_=n1;

    // #ifdef _OPENACC
    stride0_ = 1;
    stride1_ = n0_;
    // #else
    // stride0_ = n1_;
    // stride1_ = 1;
    // #endif

    total_size_ = n0_*n1_;
    data_ = data;

  }
  #pragma omp end declare target
  // Overloaded function call operator for accessing elements with optional bounds checking
  #pragma omp declare target  
  #pragma acc routine seq
  T &operator () (int i, int nv) const {
    int indx = i*stride0_+nv*stride1_;
    #if CHECK_ARY_MEM_BOUNDS
    if( i >= n0_ || i < 0 ){
      printf("! Array2DTemplated. First index out of bound\n");
      printf("  Accessing index %d for n0_ size %d\n", i, n0_);
      printf("  Abort()\n");
      abort();
    }
    if( nv >= n1_ || nv < 0 ){
      printf("! Array2DTemplated. Second index out of bound\n");
      printf("  Accessing index %d for n1_ size %d\n", nv, n1_);
      printf("  Abort()\n");
      abort();
    }
    if (indx < 0 || indx >= total_size_ ){
      printf("! Array2DTemplated. Index out of bound.\n");
      printf("  indx: %d\n", indx);
      printf("  total_size_: %d\n", total_size_);
      printf("  stride0_: %d\n", stride0_);
      printf("  stride1_: %d\n", stride1_);
      printf("  Error accessing array at %d, %d\n",i,nv);
      printf("  Abort()\n");
      abort();
    }
    #endif
    return data_[indx];
  }
  #pragma omp end declare target
  // Deallocates memory and resets metadata fields
  #pragma omp declare target
  #pragma acc routine seq
  void destroy() {
    // #pragma acc routine seq
    delete [] data_;
    data_=nullptr;
    stride0_=stride1_=-1;
    n0_=n1_=0;
  }
  #pragma omp end declare target
  protected:
    T  *data_; // Pointer to array data
    int n0_,n1_; // Number of elements
    int stride0_,stride1_; // Strides for indexing
    int total_size_; // Total allocated size
};

// Type aliases for commonly used array types
using Ary2D=Array2DTemplated<double,true>;
using Ary2DInt=Array2DTemplated<int,true>;
// using AryNotManaged2D=Array2DTemplated<double,false>;
// using Ary2DNotManagedInt=Array2DTemplated<int,false>;

// --------------------------------------------------------
// 3. Template class for 3D arrays 
// --------------------------------------------------------

template<typename T>
class Array3DTemplated {
public:
  // Default constructor initializes an empty array
//  Array3DTemplated() : stride0_(-1), stride1_(-1), stride2_(-1), data_(nullptr) {}
  Array3DTemplated() : stride0_(-1), stride1_(-1), stride2_(-1),
                       data_(nullptr), location_(-1) {}

//  Constructor that allocates an array of size n0*n1*n2
  // Array3DTemplated(int n0, int n1, int n2) {
  //   create(n0,n1,n2);
  // }

  Array3DTemplated(int n0, int n1, int n2) : location_(-1) {
    create(n0,n1,n2);
  }

  // The location inside the cell
  int location_; 

  // Deleted assignment operator to prevent accidental copying
  Array3DTemplated operator=(const Array3DTemplated&) = delete;
  // Deleted copy constructor to prevent accidental copying
  Array3DTemplated(const Array3DTemplated &) = delete;

  // Destructor deallocates memory
  ~Array3DTemplated() {
    destroy();
  }

  // Checks if the array is allocated
  bool isAllocated() {
    return data_ != nullptr;
  }

  // void deepCopy(Array3DTemplated &dest) {
  //   dest.create(n0_,n1_,n2_);
  //   #pragma acc parallel loop
  //   for (int i=0;i<total_size_;++i)
  //     dest.data_[i] = data_[i];
  // }

  // Allocates memory for the array and initializes metadata
  void create(int n0, int n1, int n2) {

    n0_=n0;n1_=n1;n2_=n2;

    stride0_ = 1;
    stride1_ = n0_;
    stride2_ = n1_*n0_;

    total_size_ = n0_*n1_*n2_;
    data_ = new T[total_size_];
    g_usedMemory += n0_*n1_*n2_*sizeof(T);
#if DEBUG == YES
printf("3D memory: +%f = (number of elements %d * size %d)/1.e+6\n"
       ,((float)(n0_*n1_*n2_*sizeof(T))/1.e6), n0_*n1_*n2_, sizeof(T));
#endif
  }

  // Overloaded function call operator for accessing elements with optional bounds checking
  #pragma omp declare target
  #pragma acc routine seq
  T &operator () (int i, int j, int k) const {
    int indx;
    indx = i*stride0_+j*stride1_+k*stride2_;
    #if CHECK_ARY_MEM_BOUNDS
    if( i >= n0_ || i < 0 ){
      printf("! Array3DTemplated. First index out of bound\n");
      printf("  Accessing index %d for n0_ size %d\n", i, n0_);
      printf("  Abort()\n");
      abort();
    }
    if( j >= n1_ || j < 0){
      printf("! Array3DTemplated. Second index out of bound\n");
      printf("  Accessing index %d for n1_ size %d\n", j, n1_);
      printf("  Abort()\n");
      abort();
    }
    if( k >= n2_ || k < 0 ){
      printf("! Array3DTemplated. Third index out of bound\n");
      printf("  Accessing index %d for n2_ size %d\n", k, n2_);
      printf("  Abort()\n");
      abort();
    }
    if (indx < 0 || indx >= total_size_ ){
      printf("! Array3DTemplated. Index out of bound.\n");
      printf("  indx: %d\n", indx);
      printf("  total_size_: %d\n", total_size_);
      printf("  stride0_: %d\n", stride0_);
      printf("  stride1_: %d\n", stride1_);
      printf("  stride2_: %d\n", stride2_);
      printf("  Error accessing array at %d, %d, %d\n",i,j,k);
      printf("  Abort()\n");
      abort();
    }
    #endif
    return data_[indx];
  }
  #pragma omp end declare target
  void setValue(T value = 0) {
    #pragma acc parallel loop
    for (int i=0;i<total_size_;++i)
      data_[i] = value;
  }

  // Deallocates memory and resets metadata fields
  void destroy() {
    delete [] data_;
    data_=nullptr;
    stride0_=stride1_=stride2_=-1;
    n0_=n1_=n2_=0;
  }

  void printAry3D() const {
    if (data_==nullptr) {
        printf("Array is not allocated.\n");
        return;
    }
    printf("Dimensions: %d x %d x %d\n", n0_, n1_, n2_);
    printf("Strides: %d, %d, %d\n", stride0_, stride1_, stride2_);
    printf("Total elements: %d\n", total_size_);
    printf("Data pointer address: %p\n", static_cast<void*>(data_));
  }

protected:
  T  *data_; // Pointer to array data
  int n0_,n1_,n2_; // Number of elements
  int stride0_,stride1_,stride2_; // Strides for indexing
  int total_size_; // Total allocated size
  char fill_[20]; // Manual padding (round to 64 bytes) 
};

// Type aliases for commonly used array types
using Ary3D=Array3DTemplated<double>;
using Ary3DInt=Array3DTemplated<int>;
using Ary3Dbool=Array3DTemplated<bool>;
using Ary3Duint32_t=Array3DTemplated<uint32_t>;
// using Ary3Dchar=Array3DTemplated<char>;
using Ary3Dschar=Array3DTemplated<signed char>;

// --------------------------------------------------------
// 4. Template class for 4D arrays 
// --------------------------------------------------------

template<typename T>
class Array4DTemplated {
public:
  // Default constructor initializes an empty array
  Array4DTemplated() : stride0_(-1), stride1_(-1), stride2_(-1), stride3_(-1), data_(nullptr) {}
  
  // Constructor that allocates an array of size n0*n1*n2*n3
  Array4DTemplated(int n0, int n1, int n2, int n3) {
    create(n0,n1,n2,n3);
  }

  // Deleted assignment operator to prevent accidental copying
  Array4DTemplated operator=(const Array4DTemplated&) = delete;
  // Deleted copy constructor to prevent accidental copying
  Array4DTemplated(const Array4DTemplated &) = delete;

  // Destructor deallocates memory
  ~Array4DTemplated() {
    destroy();
  }

  // Checks if the array is allocated
  bool isAllocated() {
    return data_ != nullptr;
  }

  // void deepCopy(Array4DTemplated &dest) {
  //   dest.create(n0_,n1_,n2_, n3_);
  //   #pragma acc parallel loop
  //   for (int i=0;i<total_size_;++i)
  //     dest.data_[i] = data_[i];
  // }

  // Allocates memory for the array and initializes metadata
  #pragma omp declare target
  void create(int n0, int n1, int n2, int n3) {

    n0_=n0;n1_=n1;n2_=n2;n3_=n3;

// #ifdef _OPENACC
    stride0_ = 1;
    stride1_ = n0_;
    stride2_ = n1_*n0_;
    stride3_ = n2_*n1_*n0_;
// #else
//     stride0_ = n3_;
//     stride1_ = n3_*n0_;
//     stride2_ = n3_*n0_*n1_;
//     stride3_ = 1;
// #endif
    total_size_ = n0_*n1_*n2_*n3_;
    data_ = new T[total_size_];
    g_usedMemory += n0_*n1_*n2_*n3_*sizeof(T);
#if DEBUG == YES
printf("4D memory: +%f = (number of elements %d * size %d)/1.e+6\n"
       ,((float)(n0_*n1_*n2_*n3_*sizeof(T))/1.e6), n0_*n1_*n2_*n3_, sizeof(T));
#endif

  }
  #pragma omp end declare target

  
  // Overloaded function call operator for accessing elements with optional bounds checking
  #pragma omp declare target
  #pragma acc routine seq
  T &operator () (int i, int j, int k, int nv) const {
    int indx = i*stride0_+j*stride1_+k*stride2_+nv*stride3_;
    // printf("indx: %d\n", indx);

    #if CHECK_ARY_MEM_BOUNDS
    if( i >= n0_ || i < 0 ){
      printf("! Array4DTemplated. First index out of bound\n");
      printf("  Accessing index %d for n0_ size %d\n", i, n0_);
      printf("  Abort()\n");
      abort();
    }
    if( j >= n1_ || j < 0){
      printf("! Array4DTemplated. Second index out of bound\n");
      printf("  Accessing index %d for n1_ size %d\n", j, n1_);
      printf("  Abort()\n");
      abort();
    }
    if( k >= n2_ || k < 0 ){
      printf("! Array4DTemplated. Third index out of bound\n");
      printf("  Accessing index %d for n2_ size %d\n", k, n2_);
      printf("  Abort()\n");
      abort();
    }
    if( nv >= n3_ || nv < 0 ){
      printf("! Array4DTemplated. Fourth index out of bound\n");
      printf("  Accessing index %d for n3_ size %d\n", k, n3_);
      printf("  Abort()\n");
      abort();
    }
    if (indx < 0 || indx >= total_size_ ){
      printf("! Array4DTemplated. Index out of bound.\n");
      printf("  indx: %d\n", indx);
      printf("  total_size_: %d\n", total_size_);
      printf("  stride0_: %d\n", stride0_);
      printf("  stride1_: %d\n", stride1_);
      printf("  stride2_: %d\n", stride2_);
      printf("  stride3_: %d\n", stride3_);
      printf("  Error accessing array at %d, %d, %d, %d\n",i,j,k, nv);
      printf("  Abort()\n");
      abort();
    }
    #endif
    return data_[indx];
  }
  #pragma omp end declare target

  // void setValue(T value = 0) {
  //   #pragma acc parallel loop
  //   for (int i=0;i<total_size_;++i)
  //     data_[i] = value;
  // }

  // Deallocates memory and resets metadata fields
  void destroy() {
    delete [] data_;
    data_=nullptr;
    stride0_=stride1_=stride2_=stride3_=-1;
    n0_=n1_=n2_=n3_=0;
  }

protected:
  T  *data_; // Pointer to array data
  int n0_,n1_,n2_,n3_; // Number of elements
  int stride0_,stride1_,stride2_,stride3_; // Strides for indexing
  int total_size_; // Total allocated size
  char fill_[20]; // Manual padding (round to 64 bytes)
  
};

// Type aliases for commonly used array types
using Ary4D=Array4DTemplated<double>;
// using Ary4DInt=Array4DTemplated<int>;

  // --------------------------------------------------------
  // ACC_MALLOC Template class for 4D arrays 
  // --------------------------------------------------------

template<typename T>
class Array4DTemplated_acc {
public:
  // Default constructor initializes an empty array
  Array4DTemplated_acc() : stride0_(-1), stride1_(-1), stride2_(-1), stride3_(-1), data_(nullptr) {}
  // Constructor that allocates an array of size n0*n1*n2
  #pragma acc routine seq
  Array4DTemplated_acc(int n0, int n1, int n2,int n3) {
    create(n0,n1,n2,n3);
  }

  // Deleted assignment operator to prevent accidental copying
  Array4DTemplated_acc operator=(const Array4DTemplated_acc&) = delete;
  // Deleted copy constructor to prevent accidental copying
  Array4DTemplated_acc(const Array4DTemplated_acc &) = delete;

  // Destructor deallocates memory
  #pragma acc routine seq
  ~Array4DTemplated_acc() {
    destroy();
  }

  // Checks if the array is allocated
  #pragma acc routine seq
  bool isAllocated() {
    return data_ != nullptr;
  }

  // void deepCopy(Array4DTemplated_acc &dest) {
  //   dest.create(n0_,n1_,n2_);
  //   #pragma acc parallel loop
  //   for (int i=0;i<total_size_;++i)
  //     dest.data_[i] = data_[i];
  // }

  // Allocates memory for the array and initializes metadata
  void create(int n0, int n1, int n2, int n3) {

    n0_=n0;n1_=n1;n2_=n2;n3_=n3;

    stride0_ = 1;
    stride1_ = n0_;
    stride2_ = n1_*n0_;
    stride3_ = n2_*n1_*n0_;

    total_size_ = n0_*n1_*n2_*n3_;
    // data_ = new T[total_size_];
#ifdef _OPENACC
data_ = static_cast<T*>(acc_malloc(sizeof(T) * total_size_));
#else
    data_ = new T[total_size_];
#endif

    g_usedMemory += n0_*n1_*n2_*n3_*sizeof(T);
#ifdef DEBUG
printf("4D memory: +%f = (number of elements %d * size %d)/1.e+6\n"
       ,((float)(n0_*n1_*n2_*n3_*sizeof(T))/1.e6), n0_*n1_*n2_*n3_, sizeof(T));
#endif
  }

  // Overloaded function call operator for accessing elements with optional bounds checking
  #pragma acc routine seq
  T &operator () (int i, int j, int k, int nv) const {
    int indx = i*stride0_+j*stride1_+k*stride2_+nv*stride3_;
    #if CHECK_ARY_MEM_BOUNDS
    if( i >= n0_ || i < 0 ){
      printf("! Array4DTemplated_acc. First index out of bound\n");
      printf("  Accessing index %d for n0_ size %d\n", i, n0_);
      printf("  Abort()\n");
      abort();
    }
    if( j >= n1_ || j < 0){
      printf("! Array4DTemplated_acc. Second index out of bound\n");
      printf("  Accessing index %d for n1_ size %d\n", j, n1_);
      printf("  Abort()\n");
      abort();
    }
    if( k >= n2_ || k < 0 ){
      printf("! Array4DTemplated_acc. Third index out of bound\n");
      printf("  Accessing index %d for n2_ size %d\n", k, n2_);
      printf("  Abort()\n");
      abort();
    }
    if( nv >= n3_ || nv < 0 ){
      printf("! Array4DTemplated. Fourth index out of bound\n");
      printf("  Accessing index %d for n3_ size %d\n", k, n3_);
      printf("  Abort()\n");
      abort();
    }
    if (indx < 0 || indx >= total_size_ ){
      printf("! Array4DTemplated_acc. Index out of bound.\n");
      printf("  indx: %d\n", indx);
      printf("  total_size_: %d\n", total_size_);
      printf("  stride0_: %d\n", stride0_);
      printf("  stride1_: %d\n", stride1_);
      printf("  stride2_: %d\n", stride2_);
      printf("  Error accessing array at %d, %d, %d\n",i,j,k);
      printf("  Abort()\n");
      abort();
    }
    #endif
    return data_[indx];
  }

  // void setValue(T value = 0) {
  //   #pragma acc parallel loop
  //   for (int i=0;i<total_size_;++i)
  //     data_[i] = value;
  // }

  // Deallocates memory and resets metadata fields
  void destroy() {
    // delete [] data_;
    // data_=nullptr;
    #ifdef _OPENACC
    if (data_) {
      acc_free(data_);  // NOT delete[]
      data_ = nullptr;
    }
    #else
    delete [] data_;
    data_=nullptr;
    #endif
    stride0_=stride1_=stride2_=stride3_=-1;
    n0_=n1_=n2_=n3_=0;
  }

  // void acc_create() {
  //   #pragma acc enter data copyin(this[:1]) 
  //   #pragma acc enter data create(data_[0:total_size_]) 
  //   #pragma acc update device(data_) 
  // }
  
  // void acc_copyin() {
  //   #pragma acc enter data copyin(this[:1]) 
  //   #pragma acc enter data copyin(data_[0:total_size_]) 
  //   #pragma acc update device(data_) 
  // }

  // void acc_copyout() {
  //   #pragma acc exit data copyout(data_[0:total_size_])
  //   #pragma acc exit data delete(this[:1])
  // }
  
  // void acc_delete() {
  //   #pragma acc exit data delete(data_[0:total_size_])
  //   #pragma acc exit data delete(this[:1])
  // }

protected:
  T  *data_; // Pointer to array data
  int n0_,n1_,n2_,n3_; // Number of elements
  int stride0_,stride1_,stride2_,stride3_; // Strides for indexing
  int total_size_; // Total allocated size
  char fill_[12]; // Manual padding (round to 64 bytes) 
public:
  T* devicePtr() const { return data_; }
};

// Type aliases for commonly used array types
using Ary4D_acc=Array4DTemplated_acc<double>;
// using Ary4DInt_acc=Array4DTemplated_acc<int>;
// using Ary4Duint32_t_acc=Array4DTemplated_acc<uint32_t>;
// // using Ary4Dchar=Array4DTemplated_acc<char>;
// using Ary4Dschar_acc=Array4DTemplated_acc<signed char>;


#ifdef LOC_CHECK_ARY_MEM_BOUNDS
#undef LOC_CHECK_ARY_MEM_BOUNDS
#undef CHECK_ARY_MEM_BOUNDS
#endif
