//
/// @file lisk.hpp
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

#ifndef lisk_hpp
#define lisk_hpp

#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include <unordered_map>
#include <cln/cln.h>

/** LiSK namespace */
namespace LiSK {
    
    /** 
     * Type defs
     * @note Define subsequently used abbreviations
     */
    typedef std::complex<double> cp_d;
    
    /**
     * LiSK - Library for Evaluating Classical Polylogarithms and Li_22
     * @note A lightweight library for evaluating classical polylogarithms
     * Li_n and the special function Li_{22} based on arXiv:1601.02649.
     * Currently a double precision as well as an arbitrary precision
     * implementation is available; the latter through the CLN library by 
     * Bruno Haible ( http://www.ginac.de/CLN/ ).
     * CLN is an essential part of this library and is also used in the
     * initialisation phase of the double prescision implementation.
     * @note Although this library is header-only the actual implementation
     * is stored in lisk.cpp which is included at the end of lisk.hpp.
     *
     */
    template <typename numtype>
    class LiSK {
		
    private:
        
            // The user can change the setup as follows:
            //
            // 1.)  Set the step size _step for extending the vectors of constants.
            //      A small value reduces the initialisation time but requires the
            //      initialisation more often. Larger values vice versa. (default = 5)
            //
            // 2.)  Vary _offset to determine the order of magnitude of the 'small'
            //      imaginary part, i.e. iep = i*10^(offset-prec). With prec=17
            //      for double approximation. Otherwise prec is set during
            //      construction of the LiSK object. (default = 2)
            //
            // 3.)  Enable/disable the conversion to floating point numbers of the
            //      stored constants, i.e. the floating point approximation of rational
            //      numbers are stored in the look-up tables instead of the exact value.
            //      IMPORTANT: Although enabling (true) the conversion leads to
            //      performance improvements, it does not work well for all arguments
            //      in the complex plane. Especially real input parameters seem to pose
            //      a sincere problem. Handle with care! (default = false)
			//
			// 4.)	Influence LiSKs initialisation phase by raising the number of pre-
			//		calculated constants in _user_max_values. Higher values lead to
			//		a longer initialisation time but less re-computing during the main
			//		run. Default minimum values are set during the construction nevertheless.
        
        /** Step size
         * @note Set step size used during computation of constants
         * in case an extension is required.
         */
        const unsigned _step = 5;
        
        /** Small positive imaginary part
         * @note Small positive imaginary part which will be assigned as x->x-iep.
         * For the double precision implementation we set _iep_d= i*10^(offset-17)
         * and in the arbitrary precision case as _iep_clN = i*10^(offset-prec),
         * with prec set in constructor. The offset can be variied, but the user
         * has to know what he's doing (default = 2).
         * @{
         */
        const int _offset = 2;
        const cp_d _iep_d;
        const cln::cl_N _iep_clN;
        /** @} */
        
        /** Enable/disable floating point conversion of constants
         *  @note If this flag is set true all constants in the look-up
         *  tables are given by their floating point approximation,
         *  respecting the required precision. However, problems can occur
         *  when using real input variables. Be careful! (default = false)
         */
        const bool _cln_tofloat = false;
		
		/** User defined upper limits for various constants. 
		 *  @note The limits, which determine how many constants are pre-calculated
		 *  can be set here. The keys should not be changed. Note: A minimum will be
		 *  assumed during construction and the user choice may be overwritten. But
		 *  higher values will be respected.
		 */
		const std::unordered_map<std::string, size_t> _user_max_values {
			{"nHarmNum",3},			// Harmonic numbers H(n)
			{"nPosZeta",4},			// Zeta(n) for n>=0
			{"nNegZeta",100},		// Zeta(n) for n<0
			{"nFactorials",100},	// Factorials n!
			{"nBernNum",100},		// Bernoulli numbers B(n)
			{"nLiCij",200},			// Cij/(j+1)! from eq.(5.10)
			{"nLiBn_eq59",100}		// B_{2n}/(2n*(2n+m)!) from eq.(5.9)
		};
		
		
            // ### Do NOT change anything from here ###
		
        /** Precision value
         *  @note Only relevant in arbitrary precision case.
         *  Minimum is double precision with _prec=17
         */
        const size_t _prec;
                
        /** ErrorEstimate
         *  @note Try to estimate how far we have to expand to get the required
         *  precision. Use practical approach as also used in GiNaC, i.e.
         *  check if new result differs from previous one. Due to included
         *  zeros in the constants we need to check for them.
         *  @param res Latest computed result
         *  @param prev_res In previous step computed result
         *  @param current Check current if zero
         *  @return True, if res!=prev_res || tmp==0. False otherwise
         */
        bool _ErrorEstimate(const numtype& res, const numtype& prev_res, const numtype& current) const;
        
        /** Constants
         *  @note Collection of constants which are precomputed during 
         *  initialisation and stored. If additional constants are 
         *  required during runtime the existing set will be extended.
         *  Initial maximum values can be set in _init.
         *  @{
         */
        template <typename T>
        struct constants{
            /** Harmonic numbers H(n) for n=1,2,3,... */
            std::vector<T> HarmNum;
            /** zeta(n) for n=0,(1),2,3,4. n=1 is kept but set to 0 */
            std::vector<T> PosZeta;
            /** zeta(n) for n<0 and n odd. All negative even n give zero */
            std::vector<T> NegZeta;
            /** Factorials n! (n>170 exceeds double precision). */
            std::vector<T> Factorials;
            /** Bernoulli numbers of first kind, i.e. B_1 = -1/2. The vector is
             *  {B_0,B_1,B_2n} with n=1,2,3,...*/
            std::vector<T> BernNum;
            /** Li constants Cij/(j+1)! from eq.(5.10) */
            std::vector<std::vector<T>> LiCij;
            /** Li constants B_{2n}/(2n*(2n+m)!) from eq.(5.9) */
            std::vector<std::vector<T>> LiBn_eq59;
            /** Max values for available constants. The keys are named following the
             *  convention "n"+"vector name", e.g. nHarmNum. In addition "nweight" is set
             *  and determines the highest weigth for which the constants are computed.
             */
            std::unordered_map<std::string, size_t> max_values;
            
            /** C'stor
             *  @note Pre-compute constants needed for the evaluation of Li_n and Li_22.
             *  @param Li_n Set expected weight of classical polylogarithms
             */
            constants(const size_t Li_n);
        };
		constants<numtype> _constants;
        constants<cln::cl_R> _constants_clR;	// Real numbers because this respects the
												// different nature of rationals and floats
        /** @} */
        
        /** Initialise constants
         *  @param max_values Map setting the maximum for each constant indicated by its key value. See constants::max_value
         */
		void _init(const std::unordered_map<std::string, size_t> max_values);
		
		/** Convert CLN to desired type
		 *  @param x Value to be converted
		 *  @return x in desired type
		 */
		numtype convert(const cln::cl_R& x) const;
		
        /** Maximum of vector entries
         *  @param v Vector from which maximal value is determined, e.g. {x1,x2}
         *  @return Maximal value of v, e.g. x1 if x1>x2
         */
        template<typename T>
        T _vec_max(std::vector<T> v) const;
        
        /** Extend Harmonic numbers
         *  @note Extends HarmNum vector up to nHarmNum, starting with H_0=0.
         */
        void _HarmNum_extend();
        
        /** Extend positve Zeta values
         *  @note Extends PosZeta vector with zeta(n), n=0,1,2,... and zeta(1)=0.
         */
        void _PosZeta_extend();
        
        /** Extend factorials
         *  @note Extends Factorials vector n! with n=0,1,2,...
         */
        void _Factorials_extend();
        
        /** Extend Bernoulli numbers
         *  @note Extends BernNum vector B_n with n=0,1,2,3,... . B_n=0 for odd n>1. 
         *  Zeros are kept for now. If that poses a performance or memory problem 
         *  later on it could be changed.
         */
        void _BernNum_extend();
        
        /** Negative zeta values
         *  @note Extends NegZeta vector with zeta(-n), n=1,2,3.... zeta(-n)=0 for even n>0. 
         *  Zeros are kept for now.
         */
        void _NegZeta_extend();
        
        /** Extend Li constants LiCij
         *  @note Compute the coefficients tCij = LiCij/(j+1)! from eq.(5.10). Extends
         *  the multi-dim vector LiCij containing tCij with i=1,2,3,4 and j=0,1,2,...
         */
        void _LiCij_extend();
        
        /** Extend Li constants LiBn_eq59
         *  @note Compute the coefficients tBnm = B_{2n}/(2n*(2n+m)!) for n>0 and m=1,2,3
         */
        void _LiBn_eq59_extend();
        
        /** Li_n(1-exp(-alpha))
         *  @note Implementation of eq.(5.10) up to desired accuracy.
         *  @param n Weight of Li_n
         *  @param x Point of evaluation, i.e. x=1-exp(-alpha)
         *  @return Li_n(x) up to desired accuracy.
         */
        numtype _Lin_1mexpal(const size_t n, const numtype& x);
        
        /** Li sums
         *  @note Compute Li_n(~al) = precal + sign * sum_{m=1} B_{2m}/(2m*(2m+n-1)!)*al^{2m+n-1} as needed
         *  for the classical polylogarithms (eq.(5.9)).
         *  @param precal Terms to be added to sum
         *  @param al Accordingly mapped argument of Li_n(x). See eq.(5.9) and eq.(5.10)
         *  @param n Weight specifier Li_n
         *  @param fac_2n  If true, each summand will be multiplied by 2m
         *  @param sign If true, the entire sum is multiplied by -1
         *  @return Li_n(~al)
         */
        numtype _Li_sumBn(const numtype &precal,const numtype &al, const size_t n, const bool fac_2n=false, const bool sign=false);
        
        /** General Li_n(x)
         *  @param x Point of evaluation
         *  @param n Weight of Li_n
         *  @return Li_n(x)
         */
        numtype _Lin_basis(const size_t n, const numtype& x);
        
        /** Inversion formula for Li(n,x)
         *  @note Inversion formula from eq.(5.6). Used for Li_n with n>4.
         *  @param x Point of evaluation
         *  @param n Weight of Li_n
         *  @return Li_n(x) where |x|>1
         */
        numtype _Lin_inverse(const size_t n, const numtype& x);
        
        /** All order Li_n(exp(-alpha))
         *  @note Implementation of eq.(5.8) for all n>0.
         *  @param n Weight specifier Li_n
         *  @param x Point of evaluation, i.e. x=exp(-alpha)
         *  @return Li_n(x) up to desired accuracy
         */
        numtype _Lin_expal(const size_t n, const numtype& x);
		
		/** Power function
		 *  @note Power function for exponents of integer type.
		 *  Nothing more is needed here.
		 *  @param x Basis
		 *  @param n Exponent
		 *  @return x^n
		 */
		numtype Pow(const numtype& x, const int n) const;
		
		/** Real part
		 *  @param x Complex number
		 *  @return Re(x)
		 *  @{
		 */
		inline
		cln::cl_R Re(const cln::cl_N& x) const{
			return cln::realpart(x);
		};
		
		inline
		double Re(const std::complex<double>& x) const{
			return x.real();
		};
		/** @} */
		
		/** Is_zero
		 *  @note Using CLNs zerop() to determine if number
		 *  is zero.
		 *  @param x Complex number
		 *  @return True, if x=0. False, otherwise
		 */
		bool is_zero(const numtype& x) const;		
		
        /** Add imaginary part
         *  @note Add small imaginary part depending on precision.
         *  We follow the widely used convention for mathematical
         *  software, e.g. "The C Standard", Mathematica, GiNaC, etc.
         *  @param x Complex value
         *  @return x-iep
         */
        numtype _add_iep(const numtype& x) const;
        
        /** Adjust digits
         *  @note Adjust digits of input variables such that they fullfil
         *  desired precision. Only relevant for arbitrary precision part.
         *  @param x Complex value
         *  @return x adjusted to prec digits
         */
        inline
        cln::cl_N _adjustDigits(const cln::cl_N& x) const;
        
        /** Li(n,x)
         *  @note Classical polylogarithm of weight n at point x but no
         *  imaginary part is assigned
         *  @param n Weight of classical polylogarithm
         *  @param x Point of evaluation
         *  @return Li(n,x)
         *  @{
         */
        numtype _Li1(const numtype& x);
        numtype _Li2(const numtype& x);
        numtype _Li3(const numtype& x);
        numtype _Li4(const numtype& x);
        /** @} */
        
        /** Li22(x,y)
         *  @note No imaginary parts will be assigned.
         *  @param x Point of evaluation
         *  @param y Point of evaluation
         *  @return Li22(x,y)
         */
        numtype _Li22(const numtype& x, const numtype& y);
        
        /** Original Li22(x,y)
         *  @note Original definition of Li22(x,y) from (6.1) if
         *  |x|<1 and |xy|<1.
         *  @param x Point of evaluation
         *  @param y Point of evaluation
         *  @return Li22(x,y)
         */
        numtype _Li22_orig(const numtype& x, const numtype& y) const;
        
        /** Stuffle relation for Li22(x,y)
         *  @param x Point of evaluation
         *  @param y Point of evaluation
         *  @return -Li22(y,x) - Li4(x*y) + Li2(x)*Li2(y)
         */
        inline
        numtype _Li22_stuffle(const numtype& x, const numtype& y);
        
        /** Inversion relation for Li22(x,y)
         *  @param x Point of evaluation
         *  @param y Point of evaluation
         *  @return Inversion relation from eq.(6.4)
         */
        numtype _Li22_inversion(const numtype& x, const numtype& y);
                
    public:
        
        /** C'stor
         *  @note During initialisation (construction) constants for polylogarithms up to
         *  weight n are computed. Occurences of polylogarithms with higher weights during
         *  runtime are possible but might lead to longer evaluation times due to adaption 
         *  of the constants. If weights n>4 are expected to occur n should be specified 
         *  accordingly to reduce evaluation time during the actual computation. 
         *  If CLN is used for arbitrary precision then the desired precision is set with prec.
         *  @param n Prepare for polylogs of weight n
         *  @param prec Set precision when using CLN
         */
        LiSK(const size_t n=4, const size_t prec=34);
        
        /** Li(n,x)
         *  @note Classical polylogarithm of weight n at point x. Public wrapper
         *  for the actual private implementation. A small imaginary part is assigned,
         *  i.e. x->x+iep. For the n<=4 the abbreviations Lin(x)=Li(n,x) are supported.
         *  @param n Weight of classical polylogarithm
         *  @param x Point of evaluation
         *  @return Li(n,x)
         *  @{
         */
        numtype Li(const size_t n, const numtype x);
        numtype Li1(const numtype x);
        numtype Li2(const numtype x);
        numtype Li3(const numtype x);
        numtype Li4(const numtype x);
        /** @} */
        
        /** Li22(x,y)
         *  @note Public wrapper for Li22(x,y). A small imaginary part is assigned to x 
         *  and y, i.e. x->x+iep and y->y+iep.
         *  @param x Point of evaluation
         *  @param y Point of evaluation
         *  @return Li22(x,y)
         */
        numtype Li22(const numtype x, const numtype y);
                
    };
    
/*************************************************************************
*                   Implementation of LiSK                              *
*************************************************************************/
#include "lisk.ipp"
    
} // End LiSK namespace

#endif /* lisk_hpp */
