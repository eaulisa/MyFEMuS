
/*=========================================================================

 Program: FEMuS
 Module: MyMatrix
 Authors: Eugenio Aulisa

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <string>

#include <iostream>
#include <vector>
#include <stdlib.h>

#include <mpi.h>
#include <boost/mpi/datatype.hpp>

#include "MyMatrix.hpp"

namespace femus {

  // ******************
  template <class Type> MyMatrix<Type>::MyMatrix() {
    init();
  }

  // ******************
  template <class Type> MyMatrix<Type>::MyMatrix(const unsigned &rsize, const unsigned &csize, const Type value) {
    init();
    resize(rsize, csize, value);
  }

  // ******************
  template <class Type> MyMatrix<Type>::MyMatrix(const std::vector < unsigned > &offset, const unsigned &csize, const Type value) {
    init();
    resize(offset, csize, value);
  }

  // ******************
  template <class Type> MyMatrix<Type>::MyMatrix(const MyVector < unsigned > &rowSize, const Type value) {

    init();

    _rowSize = rowSize;
    _begin = _rowSize.begin();
    _end = _rowSize.end();
    _size = _rowSize.size();

    if((_rowSize.status()).compare("UNINITIALIZED") == 0) {
      std::cout << "In function MyMatrix(const MyVector < unsigned > &rowSize, const Type value), rowSize is UNINITIALIZED" << std::endl;
      abort();
    }
    else  if((_rowSize.status()).compare("PARALLEL") == 0) {
      _offset = _rowSize.getOffset();
      _rowOffset.resize(_offset);
      _matSize.resize(_nprocs);
      _matSize.scatter();
      _serial = false;
    }
    else {
      _rowOffset.resize(_size);
    }

    unsigned matsize = 0;
    for(unsigned i = _rowOffset.begin(); i < _rowOffset.end(); i++) {
      _rowOffset[i] = matsize ;
      matsize += _rowSize[i];
    }

    _mat.resize(matsize, value);

    if(!_serial)
      _matSize[_iproc] = matsize;

    _matIsAllocated = true;
  }

  //*******************
  template <class Type> void MyMatrix<Type>::init() {
    int iproc, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    _iproc = static_cast < unsigned >(iproc);
    _nprocs = static_cast < unsigned >(nprocs);

    Type dummy = 0;
    _MY_MPI_DATATYPE = boost::mpi::get_mpi_datatype(dummy);

    _matIsAllocated = false;
    _serial = true;

    _begin = 0;
    _end = 0;
    _size = 0;

    _lproc = _iproc;
  }

  // ******************
  template <class Type> MyMatrix<Type>::~MyMatrix() {
    clear();
  }

  // ******************
  template <class Type> void MyMatrix<Type>::resize(const unsigned &rsize, const unsigned &csize, const Type value) {

    _begin = 0;
    _end = rsize;
    _size = _end - _begin;

    _rowSize.resize(_size, csize);
    _rowOffset.resize(_size);

    unsigned matsize = 0;
    for(unsigned i = _rowOffset.begin(); i < _rowOffset.end(); i++) {
      _rowSize[i] = csize;
      _rowOffset[i] = matsize ;
      matsize += _rowSize[i];
    }

    _matSize.resize(_nprocs);
    _matSize[_iproc] = _size * csize;
    _mat.resize(_matSize[_iproc], value);

    _matIsAllocated = true;
    _serial = true;

  }

  // ******************
  template <class Type> void MyMatrix<Type>::resize(const std::vector < unsigned > &offset, const unsigned &csize, const Type value) {

    _offset = offset;
    if(_nprocs != _offset.size() - 1) {
      std::cout << "Error in MyMatrix.Resize(...), offset.size()" << _offset.size()
                << "!= from nprocs+1=" << _nprocs + 1 << std::endl;
      abort();
    }

    _begin = _offset[_iproc];
    _end = _offset[_iproc + 1];
    _size = _end - _begin;

    _rowSize.resize(offset, csize);
    _rowOffset.resize(offset);

    unsigned matsize = 0;
    for(unsigned i = _rowOffset.begin(); i < _rowOffset.end(); i++) {
      _rowSize[i] = csize;
      _rowOffset[i] = matsize ;
      matsize += _rowSize[i];
    }

    _matSize.resize(_nprocs);
    _matSize.scatter();
    _matSize[_iproc] = _size * csize;
    _mat.resize(_matSize[_iproc], value);

    _matIsAllocated = true;
    _serial = false;

  }

  // ******************
  template <class Type> void MyMatrix<Type>::shrinkToFit(Type remove) {
    MyVector <unsigned> rowSize2 = _rowSize;
    MyVector <unsigned> rowOffset2 = _rowOffset;
    unsigned counter = 0;
    for(unsigned i = begin(); i < end(); i++) {
      rowOffset2[i] = counter;
      unsigned j = 0;
      while(j < end(i) && _mat[_rowOffset[i] + j] != remove) {
        counter++;
        j++;
      }
      rowSize2[i] = j;
    }

    _mat2.resize(counter);
    if(!_serial)
      _matSize[_iproc] = counter;

    for(unsigned i = begin(); i < end(); i++) {
      for(unsigned j = 0; j < rowSize2[i]; j++) {
        _mat2[rowOffset2[i] + j] = _mat[_rowOffset[i] + j];
      }
    }

    _rowSize = rowSize2;
    _rowOffset = rowOffset2;

    _mat.swap(_mat2);
    std::vector<Type>().swap(_mat2);

  }

  // ******************
  template <class Type> void MyMatrix<Type>::shrinkToFit(MyVector<unsigned> &rowSize2) {

    MyVector <unsigned> rowOffset2 = _rowOffset;
    unsigned counter = 0;
    for(unsigned i = begin(); i < end(); i++) {
      rowOffset2[i] = counter;
      counter += rowSize2[i];
    }

    _mat2.resize(counter);
    if(!_serial)
      _matSize[_iproc] = counter;

    for(unsigned i = begin(); i < end(); i++) {
      for(unsigned j = 0; j < rowSize2[i]; j++) {
        _mat2[rowOffset2[i] + j] = _mat[_rowOffset[i] + j];
      }
    }

    _rowSize = rowSize2;
    _rowOffset = rowOffset2;

    _mat.swap(_mat2);
    std::vector<Type>().swap(_mat2);

  }

  // ******************
  template <class Type> void MyMatrix<Type>::clear() {
    std::vector<Type>().swap(_mat);
    std::vector<Type>().swap(_mat2);
    std::vector<unsigned>().swap(_offset);
    _rowSize.clear();
    _rowOffset.clear();
    _matSize.clear();
    _begin = 0;
    _end = 0;
    _size = 0;
    _lproc = _iproc;
    _matIsAllocated = false;
    _serial = true;
  }

  template <class Type> Type* MyMatrix<Type>::operator[](const unsigned &i) {
    return (_rowSize[i] == 0) ? nullptr : _mat.data() + _rowOffset[i];
  }

  template <class Type> Type& MyMatrix<Type>::operator()(const unsigned &i, const unsigned &j) {
    return _mat[ _rowOffset[i] + j];
  }

  // ******************
  template <class Type> unsigned MyMatrix<Type>::size() {
    return _size;
  }

  template <class Type> unsigned MyMatrix<Type>::size(const unsigned &i) {
    return (_matIsAllocated) ? _rowSize[i] : 0;
  }

  // ******************
  template <class Type> unsigned MyMatrix<Type>::begin() {
    return _begin;
  }

  template <class Type> unsigned MyMatrix<Type>::begin(const unsigned &i) {
    return 0;
  }

  // ******************
  template <class Type> unsigned MyMatrix<Type>::end() {
    return _end;
  }

  template <class Type> unsigned MyMatrix<Type>::end(const unsigned &i) {
    return (_matIsAllocated) ? _rowSize[i] : 0;
  }

  // ******************
  template <class Type> void MyMatrix<Type>::scatter(const std::vector<unsigned> &offset) {
    _offset = offset;

    if (!_serial) {
      std::cout << "Error in MyMatrix.scatter(), matrix is in " << status() << " status" << std::endl;
      abort();
    }

    if (_nprocs != _offset.size() - 1) {
      std::cout << "Error in MyMatrix.resize(...), offset.size() != from nprocs" << std::endl;
      abort();
    }

    _begin = _offset[_iproc];
    _end = _offset[_iproc + 1];
    _size = _end - _begin;

    MyVector<unsigned> rowOffset2 = _rowOffset;

    _rowSize.scatter(_offset);
    _rowOffset.scatter(_offset);

    unsigned matsize = 0;
    for (unsigned i = _begin; i < _end; i++) {
      _rowOffset[i] = matsize;
      if (_rowSize[i] > UINT_MAX - matsize) MPI_Abort(MPI_COMM_WORLD, 1);
      matsize += _rowSize[i];
    }

    _matSize.resize(_nprocs);
    _matSize.scatter();
    _matSize[_iproc] = matsize;

    _mat.swap(_mat2);
    _mat.resize(_matSize[_iproc]);

    for (unsigned i = _offset[_iproc]; i < _offset[_iproc + 1]; i++) {
      for (unsigned j = begin(i); j < end(i); j++) {
        _mat[_rowOffset[i] + j] = _mat2[rowOffset2[i] + j];
      }
    }

    std::vector<Type>().swap(_mat2);

    _serial = false;
    _matIsAllocated = true;
  }

  template <class Type> void MyMatrix<Type>::scatter(const std::vector<unsigned> &offset, const unsigned &lproc) {

    if(!_serial) {
      std::cout << "Error in MyMatrix.scatter(), matrix is in " << status() << " status" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if(_iproc == lproc) {
      _offset = offset;
      if(_offset.size() != _nprocs + 1) {
        std::cout << "Error in MyMatrix.scatter(), offset.size() != nprocs + 1" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    else {
      _offset.resize(_nprocs + 1);
    }

    MPI_Bcast(_offset.data(), static_cast<int>(_nprocs + 1), MPI_UNSIGNED, lproc, MPI_COMM_WORLD);

    _begin = _offset[_iproc];
    _end = _offset[_iproc + 1];
    _size = _end - _begin;

    std::vector<int> matSize_send;
    std::vector<int> matOffset_send;

    if(_iproc == lproc) {
      matSize_send.resize(_nprocs, 0);
      matOffset_send.resize(_nprocs, 0);
      for(unsigned p = 0; p < _nprocs; p++) {
        unsigned matsize = 0;
        if(_offset[p] < _offset[p + 1]) {
          matOffset_send[p] = static_cast<int>(_rowOffset[_offset[p]]);
        }
        for(unsigned i = _offset[p]; i < _offset[p + 1]; i++) {
          matsize += _rowSize[i];
        }
        matSize_send[p] = static_cast<int>(matsize);
      }
    }

    _rowSize.setSerial(true);
    _rowSize.scatter(_offset, lproc);

    _rowOffset.setSerial(true);
    _rowOffset.scatter(_offset, lproc);

    unsigned matsize = 0;
    for(unsigned i = _begin; i < _end; i++) {
      _rowOffset[i] = matsize;
      matsize += _rowSize[i];
    }

    _matSize.resize(_nprocs);
    _matSize.scatter();
    _matSize[_iproc] = matsize;

    _mat.swap(_mat2);
    _mat.resize(_matSize[_iproc]);

    MPI_Scatterv(_mat2.data(), matSize_send.data(), matOffset_send.data(), _MY_MPI_DATATYPE, _mat.data(), static_cast<int>(_matSize[_iproc]), _MY_MPI_DATATYPE, lproc, MPI_COMM_WORLD);

    std::vector<Type>().swap(_mat2);

    _serial = false;
    _matIsAllocated = true;

  }

  // ******************
  template <class Type> void MyMatrix<Type>::scatter() {
    buildOffset();
    scatter(_offset);
  }

  // ******************
  template <class Type> void MyMatrix<Type>::stack(const unsigned lproc) {
    if (_serial || !_matIsAllocated) {
      std::cout << "Error in MyMatrix.stack(), matrix must be PARALLEL" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (_size > static_cast<unsigned>(INT_MAX) || _mat.size() > static_cast<unsigned>(INT_MAX)) MPI_Abort(MPI_COMM_WORLD, 1);

    int count_send = static_cast<int>(_mat.size());
    std::vector<int> counts_recv((_iproc == lproc) ? _nprocs : 0, 0);
    std::vector<int> offsets_recv((_iproc == lproc) ? _nprocs : 0, 0);

    MPI_Gather(&count_send, 1, MPI_INT, counts_recv.data(), 1, MPI_INT, lproc, MPI_COMM_WORLD);

    std::vector<unsigned> rowSizes_recv;
    std::vector<Type> values_recv;
    std::vector<MPI_Request> communication((_iproc == lproc) ? 2 * _nprocs : 2, MPI_REQUEST_NULL);
    int total_recv = 0;

    if (_iproc == lproc) {
      if (_offset[_nprocs] > static_cast<unsigned>(INT_MAX)) MPI_Abort(MPI_COMM_WORLD, 1);

      for (unsigned jproc = 0; jproc < _nprocs; jproc++) {
        if (counts_recv[jproc] > INT_MAX - total_recv) MPI_Abort(MPI_COMM_WORLD, 1);
        offsets_recv[jproc] = total_recv;
        total_recv += counts_recv[jproc];
      }

      rowSizes_recv.resize(_offset[_nprocs]);
      values_recv.resize(total_recv);

      for (unsigned i = _begin; i < _end; i++) rowSizes_recv[i] = _rowSize[i];
      std::copy(_mat.begin(), _mat.end(), values_recv.begin() + offsets_recv[lproc]);

      for (unsigned jproc = 0; jproc < lproc; jproc++) {
        unsigned rows_recv = _offset[jproc + 1] - _offset[jproc];

        if (rows_recv > 0) MPI_Irecv(rowSizes_recv.data() + _offset[jproc], static_cast<int>(rows_recv), MPI_UNSIGNED, static_cast<int>(jproc), 230, MPI_COMM_WORLD, &communication[2 * jproc]);
        if (counts_recv[jproc] > 0) MPI_Irecv(values_recv.data() + offsets_recv[jproc], counts_recv[jproc], _MY_MPI_DATATYPE, static_cast<int>(jproc), 231, MPI_COMM_WORLD, &communication[2 * jproc + 1]);
      }

      for (unsigned jproc = lproc + 1; jproc < _nprocs; jproc++) {
        unsigned rows_recv = _offset[jproc + 1] - _offset[jproc];

        if (rows_recv > 0) MPI_Irecv(rowSizes_recv.data() + _offset[jproc], static_cast<int>(rows_recv), MPI_UNSIGNED, static_cast<int>(jproc), 230, MPI_COMM_WORLD, &communication[2 * jproc]);
        if (counts_recv[jproc] > 0) MPI_Irecv(values_recv.data() + offsets_recv[jproc], counts_recv[jproc], _MY_MPI_DATATYPE, static_cast<int>(jproc), 231, MPI_COMM_WORLD, &communication[2 * jproc + 1]);
      }
    }
    else {
      if (_size > 0) MPI_Isend(&_rowSize[_begin], static_cast<int>(_size), MPI_UNSIGNED, lproc, 230, MPI_COMM_WORLD, &communication[0]);
      if (count_send > 0) MPI_Isend(_mat.data(), count_send, _MY_MPI_DATATYPE, lproc, 231, MPI_COMM_WORLD, &communication[1]);
    }

    MPI_Waitall(static_cast<int>(communication.size()), communication.data(), MPI_STATUSES_IGNORE);

    _mat.swap(values_recv);
    std::vector<Type>().swap(values_recv);
    std::vector<Type>().swap(_mat2);

    unsigned globalSize = _offset[_nprocs];

    for (unsigned i = 0; i <= lproc; i++) _offset[i] = 0;
    for (unsigned i = lproc + 1; i < _offset.size(); i++) _offset[i] = globalSize;

    _rowSize.clear();
    _rowOffset.clear();
    _rowSize.resize(_offset);
    _rowOffset.resize(_offset);

    _begin = _offset[_iproc];
    _end = _offset[_iproc + 1];
    _size = _end - _begin;

    unsigned offset = 0;
    for (unsigned i = _begin; i < _end; i++) {
      _rowSize[i] = rowSizes_recv[i];
      _rowOffset[i] = offset;
      offset += rowSizes_recv[i];
    }

    _matSize[_iproc] = static_cast<unsigned>(total_recv);
    _lproc = lproc;
    _serial = false;
    _matIsAllocated = true;
  }

  // ******************
  template <class Type> void MyMatrix<Type>::buildOffset() {
    unsigned locsize = _size / _nprocs;
    unsigned reminder = _size % _nprocs;
    _offset.resize(_nprocs + 1);
    _offset[0] = 0;
    for(unsigned i = 1; i < _nprocs; i++) {
      _offset[i] = _offset[i - 1] + locsize;
      if(i <= reminder) _offset[i] += 1;
    }
    _offset[_nprocs] = _size;
  }

  // ******************
  template <class Type> std::vector<unsigned>  MyMatrix<Type>::getOffset() {
    return _offset;
  }

  // ******************

  template <class Type> void MyMatrix<Type>::broadcast(const unsigned &lproc) {
    if (_serial) {
      std::cout << "Error in MyMatrix.broadcast(), matrix is in " << status() << " status" << std::endl;
      abort();
    }

    _matSize.broadcast(lproc);
    if (_matSize[lproc] > static_cast<unsigned>(INT_MAX)) MPI_Abort(MPI_COMM_WORLD, 1);

    _rowSize.broadcast(lproc);
    _rowOffset.broadcast(lproc);

    if (_iproc != lproc) {
      _mat.swap(_mat2);
      _mat.resize(_matSize[lproc]);
    }

    MPI_Bcast(_mat.data(), static_cast<int>(_matSize[lproc]), _MY_MPI_DATATYPE, lproc, MPI_COMM_WORLD);

    _begin = _offset[lproc];
    _end = _offset[lproc + 1];
    _size = _end - _begin;
    _lproc = lproc;
  }

  // ******************
  template <class Type> void MyMatrix<Type>::clearBroadcast() {

    _matSize.clearBroadcast();
    _rowSize.clearBroadcast();
    _rowOffset.clearBroadcast();

    if(_lproc != _iproc) {
      _mat.swap(_mat2);
      std::vector<Type>().swap(_mat2);
      _begin = _offset[_iproc];
      _end = _offset[_iproc + 1];
      _size = _end - _begin;
    }
  }

  // ****************
  template <class Type> const std::string & MyMatrix<Type>::status() {

    if(!_matIsAllocated)
      _status = "UNINITIALIZED";
    else if(_serial)
      _status = "SERIAL";
    else
      _status = "PARALLEL";

    return _status;
  }

  // ******************

  // Explicit template instantiation
  template class MyMatrix<float>;
  template class MyMatrix<double>;
  template class MyMatrix<long double>;
  template class MyMatrix<int>;
  template class MyMatrix<short int>;
  template class MyMatrix<long int>;
  template class MyMatrix<short unsigned int>;
  template class MyMatrix<unsigned int>;
  template class MyMatrix<long unsigned int>;
  template class MyMatrix<char>;

} //end namespace femus
