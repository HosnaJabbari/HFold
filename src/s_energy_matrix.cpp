/***************************************************************************
                          s_energy_matrix.cpp  -  description
                             -------------------
    begin                : Fri Apr 12 2002
    copyright            : (C) 2002 by Mirela Andronescu
    email                : andrones@cs.ubc.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

 // This is the V matrix

#include "constants.h"
#include "h_struct.h"
#include "h_common.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <string>

#include "s_energy_matrix.h"




s_energy_matrix::s_energy_matrix (std::string seq, cand_pos_t length, vrna_param_t *params)
// The constructor
{
    params_ = params;
    make_pair_matrix();
    S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);

    n = length;
    seq_= seq;
    

    // an vector with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
	index.resize(n);
    cand_pos_t total_length = (n *(n+1))/2;
    index[0] = 0;
    for (cand_pos_t i=1; i < n; i++)
        index[i] = index[i-1]+n-i+1;

	WM.resize(total_length,INF);
    // this array holds V(i,j), and what (i,j) encloses: hairpin loop, stack pair, internal loop or multi-loop
	nodes.resize(total_length);
}


s_energy_matrix::~s_energy_matrix ()
// The destructor
{
	free(S_);
	free(S1_);
}

/**
 * @brief Gives the WM(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
 * I am adding +1 to all S as I haven't shifted the variables over to 1->n instead of 0->n-1
 * 
 * @param vij The V(i,j) energy
 * @param vi1j The V(i+1,j) energy
 * @param vij1 The V(i,j-1) energy
 * @param vi1j1 The V(i+1,j-1) energy
*/
energy_t s_energy_matrix::E_MLStem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, paramT* params,cand_pos_t i, cand_pos_t j, const  cand_pos_t& n, std::vector<Node> &tree){

	energy_t e = INF,en=INF;

	pair_type type = pair[S[i]][S[j]];

	
	if ((tree[i].pair < -1 && tree[j].pair < -1) || (tree[i].pair == j)) {
		en = vij; // i j
		if (en != INF) {
			if (params->model_details.dangles == 2){
				base_type mm5 = i>1 ? S[i-1] : -1;
            	base_type mm3 = j<n ? S[j+1] : -1;
				en += E_MLstem(type, mm5, mm3, params);
			}
			else{
				en += E_MLstem(type, -1, -1, params);
			}
			e = MIN2(e, en);
		}
	}
	if(params->model_details.dangles == 1){
		const base_type mm5 = S[i], mm3 = S[j];

		if (((tree[i+1].pair < -1 && tree[j].pair < -1) || (tree[i+1].pair == j)) && tree[i].pair < 0) {
      		en = (j-i-1 >TURN) ? vi1j : INF; // i+1 j
      		if (en != INF) {
        		en += params->MLbase;

            	type = pair[S[i+1]][S[j]];
            	en += E_MLstem(type, mm5, -1, params);

        		e = MIN2(e, en);
      		}
    	}

		if (((tree[i].pair < -1 && tree[j-1].pair < -1) || (tree[i].pair == j-1)) && tree[j].pair < 0) {
      		en = (j-1-i>TURN) ? vij1 : INF; // i j-1
      		if (en != INF) {
       			en += params->MLbase;

            	type = pair[S[i]][S[j-1]];
            	en += E_MLstem(type, -1, mm3, params);
 
        		e = MIN2(e, en);
      		}
    	}
    	if (((tree[i+1].pair < -1 && tree[j-1].pair < -1) || (tree[i+1].pair == j-1)) && tree[i].pair < 0 && tree[j].pair<0) {
      		en = (j-1-i-1>TURN) ? vi1j1 : INF; // i+1 j-1
      		if (en != INF) {
        		en += 2 * params->MLbase;

        		type = pair[S[i+1]][S[j-1]];
        		en += E_MLstem(type, mm5, mm3, params);
        
				e = MIN2(e, en);
      		}
    	} 
		
	}


    return e;
}

/**
* @brief Computes the multiloop V contribution. This gives back essentially VM(i,j).
* 
* Added plus 1 to all S's as I haven't changed it over to 1->n from 0->n-1
* 
* @param dmli1 Row of WM2 from one iteration ago
* @param dmli2 Row of WM2 from two iterations ago 
*/
energy_t s_energy_matrix::E_MbLoop(const energy_t WM2ij, const energy_t WM2ip1j, const energy_t WM2ijm1, const energy_t WM2ip1jm1, const short* S, paramT* params, cand_pos_t i, cand_pos_t j, std::vector<Node> &tree){

	energy_t e = INF,en = INF;
  	pair_type tt  = pair[S[j]][S[i]];
	bool pairable = (tree[i].pair <-1 && tree[j].pair <-1) || (tree[i].pair == j);
	
	/* double dangles */
	switch(params->model_details.dangles){
		case 2:
			if (pairable) {
			e = WM2ij;

			if (e != INF) {

				base_type si1 = S[i+1];
				base_type sj1 = S[j-1];

				e += E_MLstem(tt, sj1, si1, params) + params->MLclosing;
			}

			}
			break;

		case 1:
			/**
			* ML pair D0
			*  new closing pair (i,j) with mb part [i+1,j-1]  
			*/
			
			if (pairable) {
        		e = WM2ij;

        		if (e != INF) {

          			e += E_MLstem(tt, -1, -1, params) + params->MLclosing;

        		}
      		}
      		/** 
			* ML pair 5
			* new closing pair (i,j) with mb part [i+2,j-1] 
			*/

      		if (pairable && tree[i+1].pair < 0) {
        		en = WM2ip1j;

        		if (en != INF) {

          			base_type si1 =  S[i+1];

          			en += E_MLstem(tt, -1, si1, params) + params->MLclosing + params->MLbase;
      
        		}
      		}
      		e   = MIN2(e, en);
			
			/** 
			* ML pair 3
			* new closing pair (i,j) with mb part [i+1, j-2] 
			*/
			if (pairable && tree[j-1].pair < 0) {
				en = WM2ijm1;

				if (en != INF) {
					base_type sj1 = S[j-1];

					en += E_MLstem(tt, sj1, -1, params) + params->MLclosing + params->MLbase; 
				}
			}
			e   = MIN2(e, en);
			/** 
			* ML pair 53
			* new closing pair (i,j) with mb part [i+2.j-2]
			*/
			if (pairable && tree[i+1].pair < 0 && tree[j-1].pair <0) {
				en = WM2ip1jm1;			

				if (en != INF) {

					base_type si1 = S[i+1];
					base_type sj1 = S[j-1];

					en += E_MLstem(tt, sj1, si1, params) + params->MLclosing + 2 * params->MLbase;
				}
			}
			e   = MIN2(e, en);
      		break;
		case 0:
			if (pairable) {
				e = WM2ij;

				if (e != INF) {
					e += E_MLstem(tt, -1, -1, params) + params->MLclosing;
				}
			}
			break; 
	}


	return e;
}

void s_energy_matrix::compute_energy_WM_restricted (cand_pos_t i, cand_pos_t j, std::vector<Node> &tree)
// compute de MFE of a partial multi-loop closed at (i,j), the restricted case
{
    if(j-i+1<4) return;
	energy_t tmp;
    // ++j;
	cand_pos_t ij = index[(i-1)]+(j-1)-(i-1);
	cand_pos_t iplus1j = index[(i-1)+1]+(j-1)-(i-1)-1;
	cand_pos_t ijminus1 = index[(i-1)]+(j-1)-1-(i-1);

	WM[ij] = E_MLStem(get_energy(i-1,j-1),get_energy(i+1-1,j-1),get_energy(i-1,j-1-1),get_energy(i+1-1,j-1-1),S_,params_,i,j,n,tree);

	if (tree[i].pair <= -1)
	{
		tmp = WM[iplus1j] + params_->MLbase;
		if (tmp < WM[ij])
		{
			WM[ij] = tmp;
		}
	}

	if (tree[j].pair <= -1)
	{
		tmp = WM[ijminus1] + params_->MLbase;
		if (tmp < WM[ij])
		{
			WM[ij] = tmp;
		}
	}

	for (cand_pos_t k=i; k < j; k++)
	{
		cand_pos_t ik = index[(i-1)]+(k-1)-(i-1);
		cand_pos_t kplus1j = index[(k-1)+1]+(j-1)-(k-1)-1;
		tmp = WM[ik] + WM[kplus1j];
		if (tmp < WM[ij])
		{
			WM[ij] = tmp;
		}
	}
    
}

energy_t s_energy_matrix::compute_energy_VM_restricted (cand_pos_t i, cand_pos_t j, std::vector<Node> &tree)
// compute the MFE of a multi-loop closed at (i,j), the restricted case
{
    energy_t min = INF;
	// i--;
	// j--;
    for (cand_pos_t k = i+2; k <= j-3; ++k)
    {
		// Original versions (i will be subracting 1 from the indices for WM2 so that i can switch 0->n-1 to 1->n)
		// energy_t WM2ij = get_energy_WM(i+1,k) + get_energy_WM(k+1,j-1);
        // energy_t WM2ip1j = get_energy_WM(i+2,k) + get_energy_WM(k+1,j-1);
        // energy_t WM2ijm1 = get_energy_WM(i+1,k) + get_energy_WM(k+1,j-2);
        // energy_t WM2ip1jm1 = get_energy_WM(i+2,k) + get_energy_WM(k+1,j-2);

        energy_t WM2ij = get_energy_WM(i+1-1,k-1) + get_energy_WM(k+1-1,j-1-1);
        energy_t WM2ip1j = get_energy_WM(i+2-1,k-1) + get_energy_WM(k+1-1,j-1-1);
        energy_t WM2ijm1 = get_energy_WM(i+1-1,k-1) + get_energy_WM(k+1-1,j-2-1);
        energy_t WM2ip1jm1 = get_energy_WM(i+2-1,k-1) + get_energy_WM(k+1-1,j-2-1);
        min = std::min(min,E_MbLoop(WM2ij,WM2ip1j,WM2ijm1,WM2ip1jm1,S_,params_,i,j,tree));
    }
    return min;
}


/**
 * @brief This code returns the hairpin energy for a given base pair.
 * @param i The left index in the base pair
 * @param j The right index in the base pair
*/
PARAMTYPE s_energy_matrix::HairpinE(const std::string& seq, const short* S, const short* S1,  const paramT* params, cand_pos_t i, cand_pos_t j) {
	
	const int ptype_closing = pair[S[i]][S[j]];

	if (ptype_closing==0) return INF;

	return E_Hairpin(j-i-1,ptype_closing,S1[i+1],S1[j-1],&seq.c_str()[i-1], const_cast<paramT *>(params));
}

/**
 * @brief restricted version
*/
energy_t s_energy_matrix::compute_internal_restricted(cand_pos_t i, cand_pos_t j, const paramT *params, std::vector<int> &up){
	energy_t v_iloop = INF;
	cand_pos_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
	const int ptype_closing = pair[S_[i]][S_[j]];
	for ( cand_pos_t k=i+1; k<=max_k; ++k) {
		
		cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
        if((up[k-1]>=(k-i-1))){
            for (int l=j-1; l>=min_l; --l) {
                if(up[j-1]>=(j-l-1)){
                    energy_t v_iloop_kl = E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<paramT *>(params)) + get_energy(k-1,l-1);
                    v_iloop = std::min(v_iloop,v_iloop_kl);
                }
            }
        }
	}
	return v_iloop;
}

energy_t s_energy_matrix::compute_stack(cand_pos_t i, cand_pos_t j, const paramT *params){

	const int ptype_closing = pair[S_[i]][S_[j]];
	cand_pos_t k = i+1;
    cand_pos_t l = j-1;
    return E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<paramT *>(params)) + get_energy(k-1,l-1);
}

PARAMTYPE s_energy_matrix::compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, const paramT *params){

	const int ptype_closing = pair[S_[i]][S_[j]];
    return E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<paramT *>(params)) + get_energy(k-1,l-1);
}

void s_energy_matrix::compute_energy_restricted (cand_pos_t i, cand_pos_t j, sparse_tree &tree)
// compute the V(i,j) value, if the structure must be restricted
{
    energy_t min, min_en[3];
    cand_pos_t k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;


	const pair_type ptype_closing = pair[S_[i]][S_[j]];
    const bool unpaired = (tree.tree[i].pair<-1 && tree.tree[j].pair<-1);
	const bool paired = (tree.tree[i].pair == j && tree.tree[j].pair == i);
    if (paired || unpaired)    // if i and j can pair
    {
        bool canH = !(tree.up[j-1]<(j-i-1));
        if(canH) min_en[0] = HairpinE(seq_,S_,S1_,params_,i,j);

        min_en[1] = compute_internal_restricted(i,j,params_,tree.up);
        min_en[2] = compute_energy_VM_restricted(i,j,tree.tree);
    }

    for (k=0; k<3; k++)
    {
        if (min_en[k] < min)
        {
            min = min_en[k];
            min_rank = k;
        }
    }

    switch (min_rank)
    {
        case  0: type = HAIRP; break;
        case  1: type = INTER; break;
        case  2: type = MULTI; break;
        default: type = NONE;
    }

    // if (min < INF/2) {
    //     int ij = index[i]+j-i;
    //     nodes[ij].energy = min;
    //     nodes[ij].type = type;
    // }

	 if (min < INF/2) {
        int ij = index[(i-1)]+(j-1)-(i-1);
        nodes[ij].energy = min;
        nodes[ij].type = type;
    }
}

//Mateo 13 Sept 2023
void s_energy_matrix::compute_hotspot_energy (cand_pos_t i, cand_pos_t j, bool is_stack)
{
    //printf("in compute_hotspot_energy i:%d j:%d\n",i,j);
    energy_t energy = 0;
    if(is_stack){
        energy = compute_stack(i+1,j+1,params_);
        // printf("stack: %d\n",energy);
    }else{
        energy = HairpinE(seq_,S_,S1_,params_,i+1,j+1);
        // printf("hairpin: %d\n",energy);
    }
        
    cand_pos_t ij = index[i]+j-i;
    nodes[ij].energy = energy;
    return;
}
