//
//  lisk.ipp
//  LiSK
//
//  Created by Sebastian on 08/03/16.
//  Copyright © 2016 Sebastian Kirchner. All rights reserved.
//
//	This file is part of LiSK.
//
//	LiSK is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	LiSK is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with LiSK.  If not, see <http://www.gnu.org/licenses/>.

#include "lisk.hpp"
#include <string>
#include <cmath>
#include <cassert>
#include <type_traits>


/************************************************************************
*                                                                       *
*                   Implementation of LiSK                              *
*                                                                       *
*************************************************************************/
/*******************************************
 *  C'stors                                *
 *******************************************/
template <typename numtype>
LiSK<numtype>::LiSK(const size_t n, const size_t prec) :
        _prec{prec<17 ? 17:prec},
        _constants_clR{n},
        _constants{n},
        _iep_d{0.,std::pow(10,-(17-_offset))},
        _iep_clN{cln::complex(0,cln::cl_F(("1.e-"+std::to_string(prec-_offset)+"_"+std::to_string(prec)).c_str()))}
{
    
        // Check for supported argument type
    if (!std::is_same<numtype, cp_d>::value && !std::is_same<numtype, cln::cl_N>::value) {
        throw std::runtime_error("Momentarily LiSK only supports 'std::complex<double>' and 'cln::cl_N' as types!");
    }
    
        // Initilise constants by respecting user choice and the following default minimum
	std::unordered_map<std::string, size_t> default_max_values = {
            {"nHarmNum",3},{"nPosZeta",4},{"nNegZeta",100},{"nFactorials",100},{"nBernNum",100},{"nLiCij",100},{"nLiBn_eq59",50}};
	
	for (auto x : _user_max_values) {
		default_max_values[x.first] = std::max(default_max_values[x.first], x.second);
	}
	
    _init(default_max_values);
	
}

template <typename numtype>
template <typename T>
LiSK<numtype>::constants<T>::constants(const size_t Li_n) {
    
    max_values["nweight"] = Li_n;
    
        // Adapt vector size of expected weights n for Li_n
        // LiBn_eq59 is specialised to n<=4, where n=1 does not
        // contribute
    LiCij.resize(Li_n);
    LiBn_eq59.resize(3);
}


/************************************************************************
 *                   Public Wrappers									*
 ************************************************************************/
/********************************************
 * Li(n,x)                                  *
 ********************************************/
template <typename numtype>
numtype LiSK<numtype>::Li(const size_t n, const numtype x) {
	
    assert(n>0);
    switch (n) {
        case 1:
            return _Li1(_add_iep(x));
            
        case 2:
            return _Li2(_add_iep(x));
            
        case 3:
            return _Li3(_add_iep(x));
            
        case 4:
            return _Li4(_add_iep(x));
            
        default:
            return _Lin_basis(n, _add_iep(x));
    }
}

template <typename numtype> inline
numtype LiSK<numtype>::Li1(const numtype x) {
    return _Li1(_add_iep(x));
}

template <typename numtype> inline
numtype LiSK<numtype>::Li2(const numtype x) {
    return _Li2(_add_iep(x));
}

template <typename numtype> inline
numtype LiSK<numtype>::Li3(const numtype x) {
    return _Li3(_add_iep(x));
}

template <typename numtype> inline
numtype LiSK<numtype>::Li4(const numtype x) {
    return _Li4(_add_iep(x));
}

/********************************************
 * Li22(x,y)                                *
 ********************************************/
template <typename numtype>
numtype LiSK<numtype>::Li22(const numtype x, const numtype y) {
	
	const auto one = Re((numtype)1);
	const numtype tx = _add_iep(x), ty = _add_iep(y);
	
		//Special cases from (A.1)-(A.3) of 1601.02649
	if (is_zero(x) || is_zero(y)) return 0;
	
	if (is_zero(x-one) && is_zero(y-one)) return (((numtype)3)*_constants.PosZeta[4])/((numtype)4);
	if (is_zero(x+one) && is_zero(y+one)) return -(((numtype)3)*_constants.PosZeta[4])/((numtype)16);
	
	if (is_zero(one/x-y)) {
		
			// Check if inside of radius of convergence or if stuffle
			// relation is needed
		if (abs(x)<=1) {
			const numtype lmty = log(-ty), Li2 = _Li2(ty), &zeta2 = _constants.PosZeta[2];
			return ((numtype)3)*_Li4(ty) - Li2*(Li2 + lmty*lmty + ((numtype)6)*zeta2)/((numtype)2) - (((numtype)2)*zeta2*zeta2)/((numtype)5);
			
		} else {
				// Li_22(y,x)
			const numtype lmtx = log(-tx), Li2x = _Li2(tx), &zeta2 = _constants.PosZeta[2];
			const numtype Li22_yx = ((numtype)3)*_Li4(tx) - Li2x*(Li2x + lmtx*lmtx + ((numtype)6)*zeta2)/((numtype)2)
                    - (((numtype)2)*zeta2*zeta2)/((numtype)5);
			
			return -Li22_yx - _Li4(tx*ty) + Li2x*_Li2(ty);
		}
	}
	
	return _Li22(tx, ty);
}


/************************************************************************
 *                   Initialisation										*
 ************************************************************************/
/************************************************
 * Auxillary functions for Init                 *
 ************************************************/
    // Max value from vector
template <typename numtype>
template <typename T>
T LiSK<numtype>::_vec_max(const std::vector<T> v) const {
    
    T tmp = (v.size()==0 ? T() : v[0]);
    for(auto i : v) tmp = std::max(tmp, i);
    return tmp;
}

    // Harmonic numbers
template <typename numtype>
void LiSK<numtype>::_HarmNum_extend(){
    
        // nstart and nend define the start and end index of H_n
    auto &v = _constants_clR.HarmNum;
    const size_t nstart = v.size();
    const size_t nend   = _vec_max<size_t>({_constants_clR.max_values["nHarmNum"]+1,nstart});
    
    for(size_t n=nstart; n<nend; ++n){
        cln::cl_RA tmp = 0;
        if(n>0){
            for(size_t k=1;k<=n; ++k){
                tmp += cln::cl_R("1")/k;
            }
        }
        v.push_back(tmp);
    }
}

    // Positve zeta values
template <typename numtype>
void LiSK<numtype>::_PosZeta_extend() {
    
        // nstart and nend define the start and end index of zeta(n)
    auto &v = _constants_clR.PosZeta;
    const size_t nstart = v.size();
    const size_t nend   = _vec_max<size_t>({_constants_clR.max_values["nPosZeta"]+1,nstart});
    
    for (size_t n=nstart; n<nend; ++n) {
		if(n==0) v.push_back(-cln::cl_R("1/2"));
        else if (n==1) v.push_back(0);
        else v.push_back(cln::zeta((int)n, cln::float_format(_prec)));
    }
}

    // Factorials
template <typename numtype>
void LiSK<numtype>::_Factorials_extend() {
    
        // nstart and nend define the start and end index of n!
    auto &v = _constants_clR.Factorials;
    const size_t nstart = v.size();
    const size_t nend   = _vec_max<size_t>({_constants_clR.max_values["nFactorials"]+1,nstart});
    
    for (size_t n=nstart; n<nend; ++n) {
        v.push_back(cln::factorial((uintL)n));
    }
}

    // Bernoulli numbers. Also zeros are stored for now.
template <typename numtype>
void LiSK<numtype>::_BernNum_extend() {
    
        // nstart and nend define the start and end index of B_n
    auto &v = _constants_clR.BernNum;
    const size_t nstart = v.size();
    const size_t nend   = _vec_max<size_t>({_constants_clR.max_values["nBernNum"]+1,nstart});
    
    for (size_t n=nstart; n<nend; ++n) {
		if (n==0) v.push_back(cln::cl_I(1));
        else{
            cln::cl_R tmp = 0;
            for (size_t k=0; k<n; ++k) {
                tmp += cln::binomial((uintL)n+1, (uintL)k)*v[k];
            }
            v.push_back(-tmp/(n+1));
        }
    }
}

    // Negative zeta values
template <typename numtype>
void LiSK<numtype>::_NegZeta_extend() {
    
    auto &v = _constants_clR.NegZeta;
    const size_t nmax = _constants_clR.max_values["nNegZeta"];
    if (nmax==0) throw std::runtime_error("nmax==0 not allowed in NegZeta. Set nmax>0!");
    
        // nstart and nend define the start and end index of zeta(-n)
    const size_t nstart = v.size()+1;
    const size_t nend   = _vec_max<size_t>({nmax+1,nstart});
    
    for (size_t n=nstart; n<nend; ++n) {
        v.push_back(-_constants_clR.BernNum[n+1]/(n+1));
    }
}

    // Li constants Cij
template <typename numtype>
void LiSK<numtype>::_LiCij_extend() {
    
        // If requested weight is higher than previous one -> adapt size
    auto &v = _constants_clR.LiCij;
    const size_t nweight = _constants_clR.max_values["nweight"];
    if (nweight>v.size()){
        _constants_clR.LiCij.resize(nweight);
        _constants.LiCij.resize(nweight);
    }
    
        // i is vector index, i.e. i=n-1 for Li_n
    for (size_t i=0; i<nweight; ++i) {
        
            // nstart and nend define the start and end index of j in C_ij
        const size_t nstart = v[i].size();
        const size_t nend   = _vec_max<size_t>({_constants_clR.max_values["nLiCij"]+1,nstart});
        
        for (size_t j=nstart; j<nend; ++j) {
            if (!i) (j==0 ? v[0].push_back(cln::cl_I(1)) : v[0].push_back(0));
            else {
                
                cln::cl_R tmp = 0;
                for (size_t k=0; k<=j; ++k) {
                        // recursion formula from eq.(5.11). Multiply again by (k+1)! since (j+1)! is included in LiCij vector definition
                    tmp += cln::binomial((uintL)j, (uintL)k)*_constants_clR.BernNum[j-k]*v[i-1][k]/(k+1)*_constants_clR.Factorials[k+1];
                }
                v[i].push_back(tmp/_constants_clR.Factorials[j+1]);
            }
        }
    }
}

    // Li constants LiBn_eq59. Specialised for 2<=n<=4
template <typename numtype>
void LiSK<numtype>::_LiBn_eq59_extend() {
    
        // m = n-2 for each Li_n, n>=2
    auto &v = _constants_clR.LiBn_eq59;
    const size_t nmax = _constants_clR.max_values["nLiBn_eq59"];
    for (size_t m=0; m<3; ++m) {
            // nstart and nend define start and end index n in LiBn_eq59
        const size_t nstart = v[m].size()+1;
        const size_t nend   = _vec_max<size_t>({nmax+1,nstart});
        
        for (size_t n=nstart; n<nend ; ++n) {
            v[m].push_back(_constants_clR.BernNum[2*n]/(2*n)/_constants_clR.Factorials[2*n+m+1]);
        }
    }
}

/********************************************
 *  Init                                    *
 ********************************************/
template <typename numtype>
void LiSK<numtype>::_init(const std::unordered_map<std::string, size_t> max_values) {
	
        // Set input values and choose the correct limits
    std::unordered_map<std::string, size_t> &mv = _constants_clR.max_values;
    for (auto x : max_values) {
        if (mv.count(x.first)>0) mv[x.first] = x.second;
        else mv.insert(x);
    }
    
    mv["nBernNum"]      = _vec_max<size_t>({mv["nBernNum"],mv["nLiCij"],mv["nNegZeta"]+1,2*mv["nLiBn_eq59"]});
    mv["nFactorials"]   = _vec_max<size_t>({mv["nFactorials"],mv["nLiCij"]+1,2*mv["nLiBn_eq59"]+3,mv["nweight"]});
    mv["nPosZeta"]      = _vec_max<size_t>({mv["nPosZeta"],mv["nweight"]});
    mv["nHarmNum"]      = _vec_max<size_t>({mv["nHarmNum"],mv["nweight"]-1});
        
        // Extend constants vectors
    _HarmNum_extend(); _PosZeta_extend(); _Factorials_extend(); _BernNum_extend();
    _NegZeta_extend(); _LiCij_extend(); _LiBn_eq59_extend();
    
		// Convert CLN constants to requested and used type
	_constants.max_values = _constants_clR.max_values;
	
	for (size_t i=_constants.HarmNum.size(); i<_constants_clR.HarmNum.size(); ++i) {
		_constants.HarmNum.push_back(convert(_constants_clR.HarmNum[i]));
	}
	for (size_t i=_constants.PosZeta.size(); i<_constants_clR.PosZeta.size(); ++i) {
		_constants.PosZeta.push_back(convert(_constants_clR.PosZeta[i]));
	}
	for (size_t i=_constants.Factorials.size(); i<_constants_clR.Factorials.size(); ++i) {
		_constants.Factorials.push_back(convert(_constants_clR.Factorials[i]));
	}
	for (size_t i=_constants.BernNum.size(); i<_constants_clR.BernNum.size(); ++i) {
		_constants.BernNum.push_back(convert(_constants_clR.BernNum[i]));
	}
	for (size_t i=_constants.NegZeta.size(); i<_constants_clR.NegZeta.size(); ++i) {
		_constants.NegZeta.push_back(convert(_constants_clR.NegZeta[i]));
	}
	for (size_t i=0; i<_constants_clR.LiCij.size(); ++i) {
		for (size_t j=_constants.LiCij[i].size(); j<_constants_clR.LiCij[i].size(); ++j) {
			_constants.LiCij[i].push_back(convert(_constants_clR.LiCij[i][j]));
		}
	}
	for (size_t i=0; i<_constants_clR.LiBn_eq59.size(); ++i) {
		for (size_t j=_constants.LiBn_eq59[i].size(); j<_constants_clR.LiBn_eq59[i].size(); ++j) {
			_constants.LiBn_eq59[i].push_back(convert(_constants_clR.LiBn_eq59[i][j]));
		}
	}
}


/************************************************************************
 *                   Implementation Li(n,x)								*
 ************************************************************************/
/********************************************
 *  Li(n,x)                                 *
 ********************************************/
template <typename numtype> inline
numtype LiSK<numtype>::_Li1(const numtype& x){
	return -log(((numtype)1)-x);
}

template <typename numtype>
numtype LiSK<numtype>::_Li2(const numtype& x){
	
	const auto one = Re((numtype)1);
	
	if (Re(x)<=one/2 && abs(x)<=one) {
		const numtype al = -log(one-x), res = al - al*al/((numtype)4);
		return _Li_sumBn(res, al, 2,true);
	}
	else if (Re(x)>one/2 && abs(x-one)<=one){
		const numtype al = -log(x), res = _constants.PosZeta[2]-al-al*al/((numtype)4)+al*log(al);
		return _Li_sumBn(res, al, 2);
	}
	else{
			// If this case is reached use inversion relation. Li(1/t) is then
			// always evaluated by the first if-condition
		return -Li2(one/x) - log(-x)*log(-x)/((numtype)2)-_constants.PosZeta[2];
	}
}

template <typename numtype>
numtype LiSK<numtype>::_Li3(const numtype& x){
	
	const auto one = Re((numtype)1);
	const numtype two = 2, three = 3, four = 4;
	
	if (Re(x)<=one/2 && abs(x)<=one) {
		return _Lin_1mexpal(3, x);
	}
	else if (Re(x)>one/2 && abs(x-one)<=one){
		const numtype al  = -log(x), alsq(al*al);
		const numtype res = _constants.PosZeta[3]-_constants.PosZeta[2]*al+(three*alsq)/four+alsq*al/(three*four)-alsq/two*log(al);
		return _Li_sumBn(res, al, 3,false,true);
	}
	else{
			// If this case is reached use inversion relation. Li(1/t) is then
			// always evaluated by the first if-condition
		const numtype lx = log(-x);
		return Li3(one/x) - Pow(lx, 3)/(two*three) - _constants.PosZeta[2]*lx;
	}
}

template <typename numtype>
numtype LiSK<numtype>::_Li4(const numtype& x){
	
	const auto one = Re((numtype)1);
	const numtype two = 2, four = 4, six = 6, eight = 8, eleven = 11;
	
	if (Re(x)<=one/2 && abs(x)<=one) {
		return _Lin_1mexpal(4, x);
	}
	else if (Re(x)>one/2 && abs(x-one)<=one){
		const numtype al  = -log(x), alsq(al*al);
		const numtype res = _constants.PosZeta[4] - _constants.PosZeta[3]*al + _constants.PosZeta[2]/two*alsq - (eleven*alsq*al)/(six*six)
                - alsq*alsq/(six*eight) + alsq*al/six*log(al);
		return _Li_sumBn(res, al, 4);
	}
	else{
			// If this case is reached use inversion relation. Li(1/t) is then
			// always evaluated by the first if-condition
		const numtype lxsq = log(-x)*log(-x), &zeta2 = _constants.PosZeta[2];
		return -Li4(one/x) - lxsq*lxsq/(four*six) - zeta2/two*lxsq - (((numtype)7)*zeta2*zeta2)/((numtype)10);
	}
}

template <typename numtype>
numtype LiSK<numtype>::_Lin_basis(const size_t n, const numtype& x){
	
	const auto one = Re((numtype)1);
	
	if (Re(x)<=one/2 && abs(x)<=one) {
		return _Lin_1mexpal(n, x);
	}
	else if (Re(x)>one/2 && abs(x-one)<=one){
		return _Lin_expal(n, x);
	}
	else{
			// If this case is reached use inversion relation. Li(1/t) is then
			// always evaluated by the first if-condition
		return _Lin_inverse(n, x);
	}
}

/********************************************
 *  Li(n,x=1-exp(-alpha)), eq.(5.10)        *
 ********************************************/
template <typename numtype>
numtype LiSK<numtype>::_Lin_1mexpal(const size_t n, const numtype& x){
	
	assert(n>0); // Li_n -> n>0
	const auto one = Re((numtype)1);
	const numtype al = -log(one-x);
	numtype alpow(al), res = 0, tmp = 0, prev_res = 0;
	size_t j = 0;
	
		// Adapt constants if necessary
	if (n>_constants_clR.max_values["nweight"]) {
		_constants_clR.max_values["nweight"] = n;
		_init(_constants_clR.max_values);
	}
	
	do {
			// Compute more constants if necessary
		if (j>_constants_clR.max_values["nLiCij"]){
			_constants_clR.max_values["nLiCij"] = j + _step;
			_init(_constants_clR.max_values);
		}
		prev_res = res;
		tmp      = _constants.LiCij[n-1][j]*alpow;
		res      += tmp;
		alpow    *= al;
		j++;
	} while ( _ErrorEstimate(res, prev_res, tmp) );
	
	return res;
}

/********************************************
 *  Li(n,x=exp(-alpha)), eq.(5.8)           *
 ********************************************/
template <typename numtype>
numtype LiSK<numtype>::_Lin_expal(const size_t n, const numtype& x){
	
	assert(n>0); // Li_n -> n>0
	const numtype al = -log(x);
	numtype alpow = 1, tmp = 0, res = 0, prev_res = 0;
	numtype sgn = 1; size_t m = 0;
	auto &mv = _constants_clR.max_values;
	
		// Adapt constants if necessary
	if (n-1>mv["nHarmNum"] || n-1>mv["nFactorials"] || n>mv["nPosZeta"]) {
		mv["nHarmNum"] = mv["nFactorials"] = n-1;
		mv["nPosZeta"] = n;
		_init(mv);
	}
	
	res = (Pow(-((numtype)1), (int)n-1)/_constants.Factorials[n-1])*Pow(al, (int)n-1)*(_constants.HarmNum[n-1]-log(al));
	
	do {
			// Adapt constants if necessary
		if (m>mv["nFactorials"] || std::abs((int)(n-m))>mv["nNegZeta"]){
			mv["nFactorials"]   = m + _step;
			mv["nNegZeta"]      = std::abs((int)(n-m)) + _step;
			_init(mv);
		}
		
		prev_res = res;
		const long ndiff = n-m;
		if (m==n-1) {
			tmp = 0;
		} else {
			tmp = (ndiff>=0 ? _constants.PosZeta[ndiff] : _constants.NegZeta[-ndiff-1])/_constants.Factorials[m]*sgn*alpow;
		}
		res     += tmp;
		alpow   *= al;
		sgn     *= -1;
		m++;
	} while ( _ErrorEstimate(res, prev_res, tmp) );
	
	return res;
}

/********************************************************
 *  Inversion formula for Li(n,x), eq.(5.6)				*
 ********************************************************/
template <typename numtype>
numtype LiSK<numtype>::_Lin_inverse(const size_t n, const numtype &x) {
	
	assert(n>=2); // n=1 is trivial
	const auto one = Re((numtype)1), two = Re((numtype)2);
	const numtype sgn = ((n-1)%2) ? -1 : 1;
	const numtype lx = log(-x);
	numtype res = 0;
	
		// Adapt constants if necessary
	if (n>_constants_clR.max_values["nweight"]) {
		_constants_clR.max_values["nweight"] = n;
		_init(_constants_clR.max_values);
	}
	
	for (int r=1; r<=n/2; ++r) {
		res += ((Pow(((numtype)2), 1-2*r)-one)*_constants.PosZeta[2*r]/_constants.Factorials[n-2*r])*Pow(lx, (int)n-2*r);
	}
	
		// Eq.(5.10) is always the correct one to call
	return sgn*_Lin_1mexpal(n, one/x) - Pow(lx, (int)n)/_constants.Factorials[n] + two*res;
}

/********************************************************
 *  sum B_{2m}/(2m(2m+n-1)!)*al^{2m*n-1} from eq.(5.9)  *
 ********************************************************/
template <typename numtype>
numtype LiSK<numtype>::_Li_sumBn(const numtype &precal, const numtype &al, const size_t n, const bool fac_2n, const bool sign){
	
	assert(n>=2 && n<=4); // only for Li_n, n=2,3,4
	const auto one = Re((numtype)1), two = Re((numtype)2);
	const numtype alsq(al*al), alfac(Pow(al,(int)n-1));
	numtype alpow(alsq), tmp = 0, res(precal), prev_res = 0;
	size_t m = 1;
	
	do {
			// Compute more constants if necessary
		if (m>_constants_clR.max_values["nLiBn_eq59"]){
			_constants_clR.max_values["nLiBn_eq59"] = m + _step;
			_init(_constants_clR.max_values);
		}
		
		prev_res = res;
		tmp      = _constants.LiBn_eq59[n-2][m-1];
		if (fac_2n) tmp *= two*m;
		if (sign)   tmp *= -one;
		tmp     *= alpow*alfac;
		res     += tmp;
		alpow   *= alsq;
		m++;
	} while (_ErrorEstimate(res, prev_res, tmp));
	
	return res;
}


/************************************************************************
 *                   Implementation Li22(x,y)							*
 ************************************************************************/
/********************************************
 *  Li22(x,y)								*
 ********************************************/
template <typename numtype>
numtype LiSK<numtype>::_Li22(const numtype& x, const numtype& y) {
	
	const auto one = Re((numtype)1);
	const auto ax = abs(x), axy = abs(x*y);
	
	if (ax<one && axy<one) {
		return _Li22_orig(x, y);
	}
	else if (ax>=one && axy>=one){
		return _Li22_inversion(x, y);
	}
	else{
		return _Li22_stuffle(x, y);
	}
}

/********************************************
 *  Original definition of Li22(x,y)        *
 *  from eq.(6.1)                           *
 ********************************************/
template <typename numtype>
numtype LiSK<numtype>::_Li22_orig(const numtype &x, const numtype &y) const {
	
    // Initilize by performing first step (i=2)
    numtype xpow = x*x, ypow = y, res = xpow*ypow/((numtype)4), tmp = 0, prev_res = 0;
    // log the sum over j
    numtype prev_y = ypow;
    size_t i = 3;
    
    do {
        xpow *= x;
        ypow *= y;
        prev_y += ypow/((numtype)((i-1)*(i-1)));
        
        // tmp and prev_res are only needed for the precision check
        tmp  = xpow*prev_y/((numtype)(i*i));
        prev_res = res;
        res  += tmp;
		i++;
        
        } while (_ErrorEstimate(res, prev_res, tmp));
    
    return res;
}

/********************************************
 *  Stuffle relation for Li22(x,y)          *
 *  from eq.(6.3)                           *
 ********************************************/
template <typename numtype> inline
numtype LiSK<numtype>::_Li22_stuffle(const numtype &x, const numtype &y) {
    return -_Li22(y, x) - _Li4(x*y) + _Li2(x)*_Li2(y);
}

/********************************************
 *  Inversion relation for Li22(x,y)        *
 *  from eq.(6.4)                           *
 ********************************************/
template <typename numtype>
numtype LiSK<numtype>::_Li22_inversion(const numtype &x, const numtype &y) {
	
	const auto one = Re((numtype)1);
	const numtype &Pisq_6 = _constants.PosZeta[2], mxy = -x*y, invx = one/x, lmxy = log(mxy), lmxy2 = lmxy*lmxy, lx2 = Pow(log(-x),2);
	
	return _Li22(invx, one/y) - _Li4(-mxy) + ((numtype)3)*(_Li4(invx)+_Li4(y)) + ((numtype)2)*(_Li3(invx)-_Li3(y))*lmxy
            + _Li2(invx)*(Pisq_6+lmxy2/((numtype)2)) + _Li2(y)*(lmxy2-lx2)/((numtype)2);
}


/************************************************************************
 *                   Specialised auxillary functions					*
 ************************************************************************/
/*************************************************************************
 *                                                                       *
 *            Double Precision Specific Implementation                   *
 *                                                                       *
 *************************************************************************/
/********************************************
 *  Convert to double                       *
 ********************************************/
template <> inline
cp_d LiSK<cp_d>::convert(const cln::cl_R &x) const {
    return {cln::double_approx(x),0};
}
    
/********************************************
 *  Add imaginary part                      *
 ********************************************/
template <> inline
cp_d LiSK<cp_d>::_add_iep(const cp_d &x) const {
    return x - _iep_d;
}


/********************************************
 *  ErrorEstimate                           *
 ********************************************/
template <> inline
bool LiSK<cp_d>::_ErrorEstimate(const cp_d& res, const cp_d& prev_res, const cp_d& current) const{
    return res!=prev_res || cln::zerop((cln::cl_R)abs(current));
}
	
/********************************************
 *  Power									*
 ********************************************/
template <> inline
cp_d LiSK<cp_d>::Pow(const cp_d &x, const int n) const{
	return std::pow(x,n);
}
	
/********************************************
 *  is_zero									*
 ********************************************/
template <> inline
bool LiSK<cp_d>::is_zero(const cp_d &x) const{
	return cln::zerop(cln::complex(x.real(), x.imag()));
}

/*************************************************************************
 *                                                                       *
 *            Arbitrary Precision Specific Implementation                *
 *                                                                       *
 *************************************************************************/
/********************************************
 *  Convert to CLN							*
 ********************************************/
template <> inline
cln::cl_N LiSK<cln::cl_N>::convert(const cln::cl_R &x) const {
	return _cln_tofloat ? cln::cl_float(x, cln::float_format(_prec)) : x;
}
	
/********************************************
 *  Adjust digits                           *
 ********************************************/
template<> inline
cln::cl_N LiSK<cln::cl_N>::_adjustDigits(const cln::cl_N &x) const{
    
    const auto prec   = cln::float_format(_prec);
    const cln::cl_F r = cln::cl_float(cln::realpart(x), prec);
    const cln::cl_F i = cln::cl_float(cln::imagpart(x), prec);
    
    return cln::complex(r, i);
}

/********************************************
 *  Add imaginary part                      *
 ********************************************/
template <> inline
cln::cl_N LiSK<cln::cl_N>::_add_iep(const cln::cl_N &x) const {
    return _adjustDigits(x) - _iep_clN;
}
    
/********************************************
 *  is_zero									*
 ********************************************/
template <> inline
bool LiSK<cln::cl_N>::is_zero(const cln::cl_N &x) const{
    return cln::zerop(x);
}
    
/********************************************
 *  ErrorEstimate                           *
 ********************************************/    
template <> inline
bool LiSK<cln::cl_N>::_ErrorEstimate(const cln::cl_N& res, const cln::cl_N& prev_res, const cln::cl_N& current) const{
    return res!=prev_res || is_zero(current);
}
	
/********************************************
 *  Power									*
 ********************************************/
template <> inline
cln::cl_N LiSK<cln::cl_N>::Pow(const cln::cl_N &x, const int n) const{
	return cln::expt(x, (cln::cl_I)n);
}
