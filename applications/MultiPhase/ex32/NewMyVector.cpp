/*=========================================================================

 Program: FEMuS
 Module: MyVector
 Authors: Eugenio Aulisa

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "MyVector.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <boost/mpi/datatype.hpp>

namespace newfemus {

  // ******************
  // ctor / dtor / mpi init
  // ******************

  template <class Type>
  MyVector<Type>::MyVector()
    : _begin(0),
      _end(0),
      _size(0),
      _iproc(0),
      _nprocs(1),
      _serial(true),
      _vecIsAllocated(false),
      _status("UNINITIALIZED"),
      _lproc(0),
      _MY_MPI_DATATYPE(MPI_DATATYPE_NULL) {

    init_mpi();
    _lproc = _iproc;
    _MY_MPI_DATATYPE = boost::mpi::get_mpi_datatype(Type());
  }

  template <class Type>
  MyVector<Type>::MyVector(const unsigned &size, const Type &value)
    : MyVector() {
    resize(size, value);
  }

  template <class Type>
  MyVector<Type>::~MyVector() = default;

  template <class Type>
  void MyVector<Type>::init_mpi() {
    int ip, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &ip);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    _iproc  = static_cast<unsigned>(ip);
    _nprocs = static_cast<unsigned>(np);
  }

  // ******************
  // memory / layout
  // ******************

  template <class Type>
  void MyVector<Type>::resize(const unsigned &size, const Type value) {
    _vec.resize(size, value);
    _vecIsAllocated = true;
    _serial = true;

    _begin = 0;
    _end   = size;
    _size  = size;

    _offset.clear();
    _offset.push_back(0);
    _offset.push_back(size);
  }

  template <class Type>
  void MyVector<Type>::resize(const std::vector<unsigned> &offset,
                              const Type value) {
    // distributed resize according to offset
    _offset = offset;
    _nprocs = static_cast<unsigned>(_offset.size() - 1);

    _begin = _offset[_iproc];
    _end   = _offset[_iproc + 1];
    _size  = _end - _begin;

    _vec.resize(_size, value);
    _vecIsAllocated = true;
    _serial = false;
  }

  template <class Type>
  void MyVector<Type>::clear() {
    _vec.clear();
    _vec2.clear();
    _offset.clear();

    _begin = _end = _size = 0;
    _vecIsAllocated = false;
    _serial = true;
    _status = "UNINITIALIZED";
  }

  // ******************
  // indexing / size
  // ******************

  template <class Type>
  unsigned MyVector<Type>::size() const {
    return _size;
  }

  template <class Type>
  unsigned MyVector<Type>::begin() const {
    return _begin;
  }

  template <class Type>
  unsigned MyVector<Type>::end() const {
    return _end;
  }

  template <class Type>
  Type &MyVector<Type>::operator[](const unsigned &i) {
    // assuming caller uses local range [begin,end)
    assert(i >= _begin && i < _end);
    return _vec[i - _begin];
  }

  template <class Type>
  const Type &MyVector<Type>::operator[](const unsigned &i) const {
    assert(i >= _begin && i < _end);
    return _vec[i - _begin];
  }

  // ******************
  // parallel layout
  // ******************

  template <class Type>
  void MyVector<Type>::buildOffset() {
    // compute global size
    unsigned gsize = 0;
    if (_serial) {
      gsize = static_cast<unsigned>(_vec.size());
    } else if (!_offset.empty()) {
      gsize = _offset.back();
    }

    _offset.assign(_nprocs + 1, 0);
    unsigned loc = gsize / _nprocs;
    unsigned rem = gsize % _nprocs;

    _offset[0] = 0;
    for (unsigned p = 1; p < _nprocs; ++p) {
      _offset[p] = _offset[p - 1] + loc;
      if (p <= rem) _offset[p] += 1;
    }
    _offset[_nprocs] = gsize;

    _begin = _offset[_iproc];
    _end   = _offset[_iproc + 1];
    _size  = _end - _begin;
  }

  template <class Type>
  std::vector<unsigned> MyVector<Type>::getOffset() const {
    return _offset;
  }

  template <class Type>
  void MyVector<Type>::scatter() {
    buildOffset();
    scatter(_offset);
  }

  template <class Type>
  void MyVector<Type>::scatter(const std::vector<unsigned> &offset) {
    if (!_serial) {
      std::cout << "Error in MyVector::scatter(offset): vector is in "
                << status() << " status" << std::endl;
      abort();
    }

    if (offset.size() != _nprocs + 1) {
      std::cout << "Error in MyVector::scatter(offset): wrong offset size"
                << std::endl;
      abort();
    }

    _offset = offset;
    _begin  = _offset[_iproc];
    _end    = _offset[_iproc + 1];
    _size   = _end - _begin;

    _vec2.swap(_vec);
    _vec.resize(_size);

    if (!_vec2.empty()) {
      std::copy(_vec2.begin() + _begin,
                _vec2.begin() + _end,
                _vec.begin());
    }

    _serial = false;
  }

  template <class Type>
  void MyVector<Type>::stack() {
    if (!_serial) {
      // local sizes → offsets
      unsigned locsize = static_cast<unsigned>(_vec.size());
      std::vector<unsigned> sizes(_nprocs, 0);

      MPI_Allgather(&locsize, 1, MPI_UNSIGNED,
                    sizes.data(), 1, MPI_UNSIGNED,
                    MPI_COMM_WORLD);

      _offset.assign(_nprocs + 1, 0);
      for (unsigned p = 0; p < _nprocs; ++p) {
        _offset[p + 1] = _offset[p] + sizes[p];
      }

      _begin = _offset[_iproc];
      _end   = _offset[_iproc + 1];
      _size  = _end - _begin;
    } else {
      _offset.clear();
      _offset.push_back(0);
      _offset.push_back(static_cast<unsigned>(_vec.size()));
      _begin = 0;
      _end   = static_cast<unsigned>(_vec.size());
      _size  = _end;
    }
  }

  // ******************
  // collectives / status
  // ******************

  template <class Type>
  void MyVector<Type>::broadcast(const unsigned &lproc) {
    if (!_vecIsAllocated) {
      std::cout << "Error in MyVector::broadcast: vector is uninitialized"
                << std::endl;
      abort();
    }

    unsigned bsize = static_cast<unsigned>(_vec.size());
    MPI_Bcast(&bsize, 1, MPI_UNSIGNED, lproc, MPI_COMM_WORLD);

    if (_iproc != lproc) {
      _vec2.swap(_vec);
      _vec.resize(bsize);
    }

    MPI_Bcast(_vec.empty() ? nullptr : _vec.data(),
              static_cast<int>(bsize),
              _MY_MPI_DATATYPE,
              lproc,
              MPI_COMM_WORLD);

    _begin = 0;
    _end   = bsize;
    _size  = bsize;

    _offset.clear();
    _offset.push_back(0);
    _offset.push_back(bsize);

    _lproc = lproc;
  }

  template <class Type>
  void MyVector<Type>::clearBroadcast() {
    if (_lproc != _iproc) {
      _vec.swap(_vec2);
      std::vector<Type>().swap(_vec2);
    }
    _lproc = _iproc;
  }

  template <class Type>
  const std::string &MyVector<Type>::status() {
    if (!_vecIsAllocated)      _status = "UNINITIALIZED";
    else if (_serial)          _status = "SERIAL";
    else                       _status = "PARALLEL";
    return _status;
  }

  template <class Type>
  void MyVector<Type>::localize(std::vector<Type> &v_local) const {
    if (_offset.empty()) {
      std::cout << "Error in MyVector::localize: offset not initialized"
                << std::endl;
      abort();
    }

    // ensure v_local large enough
    unsigned global_size = _offset.back();
    v_local.resize(global_size);

    for (unsigned kp = 0; kp < _nprocs; ++kp) {
      unsigned size_kp = _offset[kp + 1] - _offset[kp];
      MPI_Bcast(&v_local[_offset[kp]],
                static_cast<int>(size_kp),
                _MY_MPI_DATATYPE,
                kp,
                MPI_COMM_WORLD);
    }
  }

  // ******************
  // explicit instantiation
  // ******************

  template class MyVector<double>;
  template class MyVector<float>;
  template class MyVector<int>;
  template class MyVector<long int>;
  template class MyVector<short unsigned int>;
  template class MyVector<unsigned int>;
  template class MyVector<long unsigned int>;
  template class MyVector<char>;

} // end namespace newfemus
