/** MEXWRAPPER.H
 * Author: Manuel Lopez
 * Email: jmlopez.rod@gmail.com
 *	  jmlopez@math.uh.edu
 *
 * Department of Mathematics
 * University of Houston
 *
 * License: http://creativecommons.org/licenses/by-sa/3.0/
 *
 * URL: http://jmlopez-rod.github.com/mexWrapper
 *
 * History:
 *
 * @version %0.0.1%, %06/19/2010% 	Simple matrix wrapper
 * @version %0.0.2%, %07/11/2010% 	Added support for cells
 * @version %0.0.3%, %11/12/2010% 	Only for MATLAB
 * @version %0.0.5%, %12/06/2010% 	Added mlCPPMatrix(CPP only) and 
 * 					mlPrMatrix(Changes values by pointing to an mxArray) 
 * @version %0.1.0%, %04/07/2011% 	Indexing changed to match matlab indexing.
 * 					Function init has been replaced by setSize.
 * @version %0.5.0%, %10/05/2011% 	New implementation with templates. Allows
 *					the creation of multidimensional arrays.
 *					Everything is an "array": Requires the
 *					the type and dimension.
 * @version %1.0.0%, %06/23/2012% 	Made available to everyone via github.
 */

#ifndef MEXWRAPPER
#define MEXWRAPPER

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
using std::vector;

/*******************************************************************************
 * MACROS AND TYPES:
 ******************************************************************************/
#define IX int
#define CIX const IX
#define CCHAR const char*
typedef enum mxType {IN=0, PR, OUT, CPP, ML, CELL};
#ifndef N_INPUTS
#define N_INPUTS(n, error_msg) \
if (nrhs != n) mexErrMsgTxt(error_msg)
#endif
#ifndef N_OUTPUTS
#define N_OUTPUTS(n, error_msg) \
if (nlhs != n) mexErrMsgTxt(error_msg)
#endif
#ifndef GET_SCALAR
#define GET_SCALAR(ptr, type, var, err) \
if (!mxIsDouble(ptr) || mxIsComplex(ptr) || mxGetN(ptr)*mxGetM(ptr)!=1) mexErrMsgTxt(err); \
const type var =  mxGetScalar(ptr)
#endif
// THESE MACROS ARE PRIVATE: DO NOT USE!!!!
#define __WRITE_ARRAY_OPERATOR__ \
	double& operator[] (CIX& i) {return _re[i-1];} \
	const double& operator[] (CIX& i) const {return _re[i-1];} \
	double& im(CIX& i) {return _im[i-1];} \
	const double& im(CIX& i) const {return _im[i-1];}

#define __WRITE_ACCESS_FUNCTIONS__ \
	CIX* size() {return _dim;} \
	CIX& size(CIX& i) {return _dim[i-1];} \
	double* rePr() {return _re;} \
	double* imPr() {return _im;} \
	mxComplexity complexity() {return (_im == NULL) ? mxREAL : mxCOMPLEX;}

#define __WRITE_CELL_ACCESS_FUNCTIONS__ \
	CIX* size() {return _dim;} \
	CIX& size(CIX& i) {return _dim[i-1];} 
/*******************************************************************************
 * BASE CLASSES:
 ******************************************************************************/
template<int ndim>
class mxBase {
protected:
	IX _dim[ndim];
	IX _Dim[ndim];
	double* _re;
	double* _im;
public:
	mxBase() : _re(NULL), _im(NULL) {}
	mxBase(CIX* d, CIX& nd, CCHAR err): _re(NULL), _im(NULL) {
		if (nd != ndim) mexErrMsgTxt(err);
		_dim[0] = d[0]; 
		_Dim[0] = 1;
		for (IX i=1; i < ndim; ++i) {
			_dim[i] = d[i];
			_Dim[i] = _Dim[i-1]*d[i-1];
		}
	}
	virtual ~mxBase() = 0;
	__WRITE_ARRAY_OPERATOR__
	double& operator() (IX* i) {
		IX offset = i[0]-1;
		for (IX n = 1; n < ndim; ++n) offset += (i[n]-1)*_Dim[n];
		return _re[offset];
	}
	const double& operator() (IX* i) const {
		IX offset = i[0]-1;
		for (IX n = 1; n < ndim; ++n) offset += (i[n]-1)*_Dim[n];
		return _re[offset];
	}
	double& im(IX* i) {
		IX offset = i[0]-1;
		for (IX n = 1; n < ndim; ++n) offset += (i[n]-1)*_Dim[n];
		return _im[offset];
	}
	const double& im(IX* i) const {
		IX offset = i[0]-1;
		for (IX n = 1; n < ndim; ++n) offset += (i[n]-1)*_Dim[n];
		return _im[offset];
	}
	void init(CIX* d, CIX& nd, CCHAR err) {
		if (nd != ndim) mexErrMsgTxt(err);
		_dim[0] = d[0]; 
		_Dim[0] = 1;
		for (IX i=1; i < ndim; ++i) {
			_dim[i] = d[i];
			_Dim[i] = _Dim[i-1]*d[i-1];
		}
	}
	CIX dim() {return ndim;}
	__WRITE_ACCESS_FUNCTIONS__
};
template<int ndim> inline mxBase<ndim>::~mxBase() {}
/** Specialization for Matrices **/
template<>
class mxBase<2> {
protected:
	IX _dim[2];
	double* _re;
	double* _im;
public:
	mxBase(): _re(NULL), _im(NULL) {}
	mxBase(CIX* d, CIX& nd, CCHAR err): _re(NULL), _im(NULL) {
		if (nd != 2) mexErrMsgTxt(err);
		_dim[0] = d[0]; _dim[1] = d[1];
	}
	virtual ~mxBase() = 0;
	__WRITE_ARRAY_OPERATOR__
	double& operator() (CIX& i1, CIX& i2) {return _re[(i1-1)+(i2-1)*_dim[0]];}
	const double& operator() (CIX& i1, CIX& i2) const {return _re[(i1-1)+(i2-1)*_dim[0]];}
	double& im(CIX& i1, CIX& i2) {return _im[(i1-1)+(i2-1)*_dim[0]];}
	const double& im(CIX& i1, CIX& i2) const {return _im[(i1-1)+(i2-1)*_dim[0]];}
	void init(CIX* d, CIX& nd, CCHAR err) {
		if (nd != 2) mexErrMsgTxt(err);
		_dim[0] = d[0]; _dim[1] = d[1];
	}
	CIX dim() {return 2;}
	__WRITE_ACCESS_FUNCTIONS__
};
mxBase<2>::~mxBase() {}
/** Specialization for Arrays of Matrices **/
template<>
class mxBase<3> {
protected:
	IX _dim[3];
	IX _Dim[3];
	double* _re;
	double* _im;
public:
	mxBase(): _re(NULL), _im(NULL) {}
	mxBase(CIX* d, CIX& nd, CCHAR err): _re(NULL), _im(NULL) {
		if (nd != 3) mexErrMsgTxt(err);
		_dim[0] = d[0]; _dim[1] = d[1]; _dim[2] = d[2];
		_Dim[0] =    1; _Dim[1] = d[0]; _Dim[2] = d[0]*d[1];
	}
	virtual ~mxBase() = 0;
	__WRITE_ARRAY_OPERATOR__
	double& operator() (CIX& i1, CIX& i2, CIX& i3) {return _re[(i1-1)+(i2-1)*_Dim[1]+(i3-1)*_Dim[2]];}
	const double& operator() (CIX& i1, CIX& i2, CIX& i3) const {return _re[(i1-1)+(i2-1)*_Dim[1]+(i3-1)*_Dim[2]];}
	double& im(CIX& i1, CIX& i2, CIX& i3) {return _im[(i1-1)+(i2-1)*_Dim[1]+(i3-1)*_Dim[2]];}
	const double& im(CIX& i1, CIX& i2, CIX& i3) const {return _im[(i1-1)+(i2-1)*_Dim[1]+(i3-1)*_Dim[2]];}
	void init(CIX* d, CIX& nd, CCHAR err) {
		if (nd != 3) mexErrMsgTxt(err);
		_dim[0] = d[0]; _dim[1] = d[1]; _dim[2] = d[2];
		_Dim[0] =    1; _Dim[1] = d[0]; _Dim[2] = d[0]*d[1];
	}
	CIX dim() {return 3;}
	__WRITE_ACCESS_FUNCTIONS__
};
mxBase<3>::~mxBase() {}
/*******************************************************************************
 * ARRAYS:
 ******************************************************************************/
template<mxType type, int ndim>
class array : public mxBase<ndim>  {
public:
	array() {}
	virtual ~array() {}
};
/** IN **/
template<int ndim>
class array<IN, ndim> : public mxBase<ndim> {
	const mxArray*& mlArray;
public:
	array(const mxArray*& in, CCHAR err): 
	mxBase<ndim>::mxBase(mxGetDimensions(in), mxGetNumberOfDimensions(in), err), mlArray(in) {
		if (!mxIsDouble(in) || mxIsSparse(in)) mexErrMsgTxt(err);
		this->_re = mxGetPr(in);
		this->_im = mxGetPi(in);
	}
	mxArray* getPr(void) {return (mxArray*)mlArray;}
};
/** PR **/
template<int ndim>
class array<PR, ndim> : public mxBase<ndim> {
public:
	array(const mxArray*& in, CCHAR err): 
	mxBase<ndim>::mxBase(mxGetDimensions(in), mxGetNumberOfDimensions(in), err) {
		if (!mxIsDouble(in) || mxIsSparse(in)) mexErrMsgTxt(err);
		this->_re = mxGetPr(in);
		this->_im = mxGetPi(in);
	}
	void setPr(const mxArray*& in, const char* err) {
		this->init(mxGetDimensions(in), mxGetNumberOfDimensions(in), err);
		if (!mxIsDouble(in) || mxIsSparse(in)) mexErrMsgTxt(err);
		this->_re = mxGetPr(in);
		this->_im = mxGetPi(in);
	}
};
/** OUT **/
template<int ndim>
class array<OUT, ndim> : public mxBase<ndim> {
	mxArray*& mlArray;
public:
	array(mxArray*& out): mxBase<ndim>::mxBase(), mlArray(out) {
		mlArray = NULL;
	}
	void setSize(CIX* d, mxComplexity ComplexFlag) {
		if (mlArray != NULL) mxDestroyArray(mlArray);
		this->init(d, ndim, "ERROR: array<OUT, ndim>::setSize dimension mismatch.");
		mlArray = mxCreateNumericArray(ndim, d, mxDOUBLE_CLASS, ComplexFlag);
		this->_re = mxGetPr(mlArray);
		this->_im = mxGetPi(mlArray);
	}
	mxArray*& getPr(void) {return mlArray;}
	void setPr(mxArray*& pr) {
		if (mlArray != NULL) mxDestroyArray(mlArray);
		this->init(mxGetDimensions(pr), mxGetNumberOfDimensions(pr), "ERROR: array<OUT, ndim>::setPr dimension mismatch.");
		mlArray = pr;
		this->_re = mxGetPr(mlArray);
		this->_im = mxGetPi(mlArray);
	}
	void copy(array<IN, ndim>& in) {
		if (mlArray != NULL) mxDestroyArray(mlArray);
		this->init(in._dim, ndim, "ERROR: array<OUT, ndim>::copy(array<IN, ndim>) dimension mismatch.");
		mlArray = mxDuplicateArray(in.getPr());
		this->_re = mxGetPr(mlArray);
		this->_im = mxGetPi(mlArray);
	}
	void copy(double* re_data, double* im_data, CIX* d) {
		if (mlArray != NULL) mxDestroyArray(mlArray);
		this->init(d, ndim, "ERROR: array<OUT, ndim>::copy(re, im, dim) dimension mismatch.");
		if (im_data == NULL) {
			mlArray = mxCreateNumericArray(ndim, d, mxDOUBLE_CLASS, mxREAL);
		} else {
			mlArray = mxCreateNumericArray(ndim, d, mxDOUBLE_CLASS, mxCOMPLEX);
			this->_im = mxGetPi(mlArray);
			memcpy(this->_im, im_data, mxGetNumberOfElements(mlArray)*sizeof(double));
		}
		this->_re = mxGetPr(mlArray);
		memcpy(this->_re, re_data, mxGetNumberOfElements(mlArray)*sizeof(double));
	}
};
/** CPP **/
template<int ndim>
class array<CPP, ndim> : public mxBase<ndim> {
private:
	void clear() {
		if (this->_re != NULL) {
			delete [] this->_re;
			this->_re = NULL;
		}
		if (this->_im != NULL) {
			delete [] this->_im;
			this->_im = NULL;
		}
	}
public:
	array(): mxBase<ndim>::mxBase() {}
	~array() {clear();}
	void setSize(CIX* d, mxComplexity ComplexFlag) {
		clear();
		this->init(d, ndim, "ERROR: array<CPP, ndim>::setSize dimension mismatch.");
		int i=0, s=1; while(i<ndim) s*=d[i++];
		this->_re = new double[s];
		if (ComplexFlag == mxCOMPLEX) {
			this->_im = new double[s];
		}
	}
	void copy(double* re_data, double* im_data, CIX* d) {
		clear();
		this->init(d, ndim, "ERROR: array<CPP, ndim>::copy(re, im, dim) dimension mismatch.");
		int i=0, s=1; while(i<ndim) s*=d[i++];
		if (im_data == NULL) {
			this->_im = NULL;
		} else {
			this->_im = new double[s];
			memcpy(this->_im, im_data, s*sizeof(double));
		}
		this->_re = new double[s];
		memcpy(this->_re, re_data, s*sizeof(double));
	}
};
/** ML **/
template<int ndim>
class array<ML, ndim> : public mxBase<ndim> {
private:
	mxArray* mlArray;
	void clear() {
		if (mlArray != NULL) {
			mxDestroyArray(mlArray);
			mlArray = NULL;
		}
	}
public:
	array(): mxBase<ndim>::mxBase(), mlArray(NULL) {}
	~array() {clear();}
	void setSize(CIX* d, mxComplexity ComplexFlag) {
		clear();
		this->init(d, ndim, "ERROR: array<ML, ndim>::setSize dimension mismatch.");
		mlArray = mxCreateNumericArray(ndim, d, mxDOUBLE_CLASS, ComplexFlag);
		this->_re = mxGetPr(mlArray);
		this->_im = mxGetPi(mlArray);
	}
	void copy(double* re_data, double* im_data, CIX* d) {
		clear();
		this->init(d, ndim, "ERROR: array<ML, ndim>::copy(re, im, dim) dimension mismatch.");
		if (im_data == NULL) {
			mlArray = mxCreateNumericArray(ndim, d, mxDOUBLE_CLASS, mxREAL);
		} else {
			mlArray = mxCreateNumericArray(ndim, d, mxDOUBLE_CLASS, mxCOMPLEX);
			this->_im = mxGetPi(mlArray);
			memcpy(this->_im, im_data, mxGetNumberOfElements(mlArray)*sizeof(double));
		}
		this->_re = mxGetPr(mlArray);
		memcpy(this->_re, re_data, mxGetNumberOfElements(mlArray)*sizeof(double));
	}
	void copy(mxArray* m) {
		clear();
		this->init(mxGetDimensions(m), mxGetNumberOfDimensions(m), "ERROR: array<ML, ndim>::copy(mxArray*) dimension mismatch.");
		mlArray = mxDuplicateArray(m);
		this->_re = mxGetPr(mlArray);
		this->_im = mxGetPi(mlArray);
	}
	mxArray* getPr(void) {return mlArray;}
	void setPr(mxArray*& m) {
		clear();
		this->init(mxGetDimensions(m), mxGetNumberOfDimensions(m), "ERROR: array<ML, ndim>::setPr(mxArray*) dimension mismatch.");
		mlArray = m;
		this->_re = mxGetPr(mlArray);
		this->_im = mxGetPi(mlArray);
	}
};
/** CELL: (Do not use this in code, this is made for a cell object) **/
template<int ndim>
class array<CELL, ndim> : public mxBase<ndim> {
private:
	mxArray* mlArray;
	bool isRef;
public:
	array(): mxBase<ndim>::mxBase(), mlArray(NULL), isRef(false) {}
	array(mxArray* in, bool ref = false): isRef(ref),
	mxBase<ndim>::mxBase(mxGetDimensions(in), mxGetNumberOfDimensions(in), "ERROR: array<CELL, ndim>::array(mxArray*) dimension mismatch.") {
		if (ref) {
			if (!mxIsDouble(in) || mxIsSparse(in)) mexErrMsgTxt("array<CELL, ndim>::array(mxArray*): Invalid mxArray.");
			mlArray = in;
		} else mlArray = mxDuplicateArray(in);
		this->_re = mxGetPr(mlArray);
		this->_im = mxGetPi(mlArray);
	}
	void free() {
		if (!isRef && mlArray != NULL) {
			mxDestroyArray(mlArray);
			mlArray = NULL;
		}
	}
	void setSize(CIX* d, mxComplexity ComplexFlag) {
		free();
		this->init(d, ndim, "ERROR: array<CELL, ndim>::setSize dimension mismatch.");
		mlArray = mxCreateNumericArray(ndim, d, mxDOUBLE_CLASS, ComplexFlag);
		this->_re = mxGetPr(mlArray);
		this->_im = mxGetPi(mlArray);
	}
	void copy(double* re_data, double* im_data, CIX* d) {
		free();
		this->init(d, ndim, "ERROR: array<CELL, ndim>::copy(re, im, dim) dimension mismatch.");
		if (im_data == NULL) {
			mlArray = mxCreateNumericArray(ndim, d, mxDOUBLE_CLASS, mxREAL);
		} else {
			mlArray = mxCreateNumericArray(ndim, d, mxDOUBLE_CLASS, mxCOMPLEX);
			this->_im = mxGetPi(mlArray);
			memcpy(this->_im, im_data, mxGetNumberOfElements(mlArray)*sizeof(double));
		}
		this->_re = mxGetPr(mlArray);
		memcpy(this->_re, re_data, mxGetNumberOfElements(mlArray)*sizeof(double));
	}
	void copy(mxArray* m) {
		free();
		this->init(mxGetDimensions(m), mxGetNumberOfDimensions(m), "ERROR: array<ML, ndim>::copy(mxArray*) dimension mismatch.");
		mlArray = mxDuplicateArray(m);
		this->_re = mxGetPr(mlArray);
		this->_im = mxGetPi(mlArray);
	}
	mxArray* getPr(void) {return mlArray;}
	void setPr(mxArray* m) {
		free();
		this->init(mxGetDimensions(m), mxGetNumberOfDimensions(m), "ERROR: array<ML, ndim>::setPr(mxArray*) dimension mismatch.");
		mlArray = m;
		this->_re = mxGetPr(mlArray);
		this->_im = mxGetPi(mlArray);
	}
};
/** setSize functions **/
template<mxType type>
void setSize(array<type, 1>& m, CIX& d1, mxComplexity ComplexFlag) {
	int dim[1] = {d1};
	m.setSize(dim, ComplexFlag);
}
template<mxType type>
void setSize(array<type, 2>& m, CIX& d1, CIX& d2, mxComplexity ComplexFlag) {
	int dim[2] = {d1, d2};
	m.setSize(dim, ComplexFlag);
}
template<mxType type>
void setSize(array<type, 3>& m, CIX& d1, CIX& d2, CIX& d3, mxComplexity ComplexFlag) {
	int dim[3] = {d1, d2, d3};
	m.setSize(dim, ComplexFlag);
}
/*******************************************************************************
 * BASE CLASSES FOR CELLS:
 ******************************************************************************/
template<int ndim, int mxDim>
class cellBase {
protected:
	IX _dim[ndim];
	IX _Dim[ndim];
	array<CELL, mxDim> **_ix;
public:
	cellBase() : _ix(NULL) { int i=0; while(i<ndim) _dim[i++] = 0; }
	cellBase(CIX* d, CIX& nd, CCHAR err): _ix(NULL) {
		if (nd != ndim) mexErrMsgTxt(err);
		_dim[0] = d[0]; 
		_Dim[0] = 1;
		for (IX i=1; i < ndim; ++i) {
			_dim[i] = d[i];
			_Dim[i] = _Dim[i-1]*d[i-1];
		}
	}
	virtual ~cellBase() = 0;
	
	array<CELL, mxDim>& operator[] (CIX& i) {return *_ix[i-1];}
	const array<CELL, mxDim>& operator[] (CIX& i) const {return *_ix[i-1];}
	
	array<CELL, mxDim>& operator() (IX* i) {
		IX offset = i[0]-1;
		for (IX n = 1; n < ndim; ++n) offset += (i[n]-1)*_Dim[n];
		return *_ix[offset];
	}
    const array<CELL, mxDim>& operator() (IX* i) const {
		IX offset = i[0]-1;
		for (IX n = 1; n < ndim; ++n) offset += (i[n]-1)*_Dim[n];
		return *_ix[offset];
	}

	void init(CIX* d, CIX& nd, CCHAR err) {
		if (nd != ndim) mexErrMsgTxt(err);
		_dim[0] = d[0]; 
		_Dim[0] = 1;
		for (IX i=1; i < ndim; ++i) {
			_dim[i] = d[i];
			_Dim[i] = _Dim[i-1]*d[i-1];
		}
	}
	CIX dim() {return ndim;}
	__WRITE_CELL_ACCESS_FUNCTIONS__
};
template<int ndim, int mxDim> inline cellBase<ndim, mxDim>::~cellBase() {}
/** Specialization for Cell Matrices **/
template<int mxDim>
class cellBase<2, mxDim> {
protected:
	IX _dim[2];
	array<CELL, mxDim> **_ix;
public:
	cellBase(): _ix(NULL) { int i=0; while(i<2) _dim[i++] = 0; }
	cellBase(CIX* d, CIX& nd, CCHAR err): _ix(NULL) {
		if (nd != 2) mexErrMsgTxt(err);
		_dim[0] = d[0]; _dim[1] = d[1];
	}
	virtual ~cellBase() = 0;
	
	array<CELL, mxDim>& operator[] (CIX& i) {return *_ix[i-1];}
	const array<CELL, mxDim>& operator[] (CIX& i) const {return *_ix[i-1];}
	
	array<CELL, mxDim>& operator() (CIX& i1, CIX& i2) {return *_ix[(i1-1)+(i2-1)*_dim[0]];}
	const array<CELL, mxDim>& operator() (CIX& i1, CIX& i2) const {return *_ix[(i1-1)+(i2-1)*_dim[0]];}
	void init(CIX* d, CIX& nd, CCHAR err) {
		if (nd != 2) mexErrMsgTxt(err);
		_dim[0] = d[0]; _dim[1] = d[1];
	}
	CIX dim() {return 2;}
	__WRITE_CELL_ACCESS_FUNCTIONS__
};
template<int mxDim> inline cellBase<2, mxDim>::~cellBase() {}
/** Specialization for Arrays of Matrices **/
template<int mxDim>
class cellBase<3, mxDim> {
protected:
	IX _dim[3];
	IX _Dim[3];
	array<CELL, mxDim> **_ix;
public:
	cellBase(): _ix(NULL) { int i=0; while(i<3) _dim[i++] = 0; }
	cellBase(CIX* d, CIX& nd, CCHAR err): _ix(NULL) {
		if (nd != 3) mexErrMsgTxt(err);
		_dim[0] = d[0]; _dim[1] = d[1]; _dim[2] = d[2];
		_Dim[0] =    1; _Dim[1] = d[0]; _Dim[2] = d[0]*d[1];
	}
	virtual ~cellBase() = 0;

	array<CELL, mxDim>& operator[] (CIX& i) {return *_ix[i-1];}
	const array<CELL, mxDim>& operator[] (CIX& i) const {return *_ix[i-1];}
	
	array<CELL, mxDim>& operator() (CIX& i1, CIX& i2, CIX& i3) {return *_ix[(i1-1)+(i2-1)*_Dim[1]+(i3-1)*_Dim[2]];}
	const array<CELL, mxDim>& operator() (CIX& i1, CIX& i2, CIX& i3) const {return *_ix[(i1-1)+(i2-1)*_Dim[1]+(i3-1)*_Dim[2]];}
	void init(CIX* d, CIX& nd, CCHAR err) {
		if (nd != 3) mexErrMsgTxt(err);
		_dim[0] = d[0]; _dim[1] = d[1]; _dim[2] = d[2];
		_Dim[0] =    1; _Dim[1] = d[0]; _Dim[2] = d[0]*d[1];
	}
	CIX dim() {return 3;}
	__WRITE_CELL_ACCESS_FUNCTIONS__
};
template<int mxDim> inline cellBase<3, mxDim>::~cellBase() {}
/*******************************************************************************
 * CELLS:
 ******************************************************************************/
template<mxType type, int ndim, int mxDim>
class cell : public cellBase<ndim, mxDim>  {
public:
	cell() {}
	virtual ~cell() {}
};
/** IN **/
template<int ndim, int mxDim>
class cell<IN, ndim, mxDim> : public cellBase<ndim, mxDim> {
	const mxArray*& mlArray;
public:
	cell(const mxArray*& in, CCHAR err): 
	cellBase<ndim, mxDim>::cellBase(mxGetDimensions(in), mxGetNumberOfDimensions(in), err), mlArray(in) {
		if (!mxIsCell(in)) mexErrMsgTxt(err);
		int total = 1, i=0; while(i<ndim) total*=this->_dim[i++];
		this->_ix = new array<CELL, mxDim>*[total];
		for (int j=0; j < total; ++j) {
			if (mxGetCell(in, j) == NULL) mexErrMsgTxt("ERROR: cell<IN, ndim, mxDim>::cell(const mxArray*): Cell contents do not have the correct dimensions.");
			this->_ix[j] = new array<CELL, mxDim>(mxGetCell(in, j), true);
		}
	}
	~cell(){
		int total = 1, i=0; while(i<ndim) total*=this->_dim[i++];
		for (int j=0; j < total; ++j) {
			delete this->_ix[j];
		}
		delete[] this->_ix;
	}
	mxArray* getPr(void) {return (mxArray*)mlArray;}
};
/** OUT **/
template<int ndim, int mxDim>
class cell<OUT, ndim, mxDim> : public cellBase<ndim, mxDim> {
	mxArray*& mlArray;
public:
	cell(mxArray*& out): mlArray(out) {}
	~cell(){
		mlArray = mxCreateCellArray(ndim, this->_dim);
		int total = 1, i=0; while(i<ndim) total*=this->_dim[i++];
		for (int j=0; j < total; ++j) {
			if ((this->_ix[j])->getPr() != NULL){
				mxSetCell(mlArray, j, (this->_ix[j])->getPr());
			}
			delete this->_ix[j];
		}
		delete[] this->_ix;
	}
	void setSize(CIX* d) {
		int total = 1, i=0; while(i < ndim) total*=this->_dim[i++];
		for (int j=0; j < total; ++j) {
			if (this->_ix[j] != NULL) delete this->_ix[j];
		}
		delete[] this->_ix;
		
		this->init(d, ndim, "ERROR: cell<OUT, ndim, mxDim>::setSize dimension mismatch.");
		total = 1, i=0; while(i<ndim) total*=this->_dim[i++];
		this->_ix = new array<CELL, mxDim>*[total];
		for (int j=0; j < total; ++j)
			this->_ix[j] = new array<CELL, mxDim>();
	}
	mxArray* getPr(void) {return mlArray;}
};
/** setSize functions **/
template<mxType type, int mxDim>
void setSize(cell<type, 1, mxDim>& m, CIX& d1) {
	int dim[1] = {d1};
	m.setSize(dim);
}
template<mxType type, int mxDim>
void setSize(cell<type, 2, mxDim>& m, CIX& d1, CIX& d2) {
	int dim[2] = {d1, d2};
	m.setSize(dim);
}
template<mxType type, int mxDim>
void setSize(cell<type, 3, mxDim>& m, CIX& d1, CIX& d2, CIX& d3) {
	int dim[3] = {d1, d2, d3};
	m.setSize(dim);
}
#endif
