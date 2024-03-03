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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <string>
// Hosna, March 5, 2012
// malloc.h is not needed in my mac as stdlib.h does the same
//#include <malloc.h>

#include "constants.h"
#include "h_struct.h"
#include "externs.h"
#include "h_common.h"
#include "s_energy_matrix.h"




s_energy_matrix::s_energy_matrix (std::string seq, cand_pos_t length, vrna_param_t *params)
// The constructor
{
    this->VM = NULL;
    params_ = params;
    make_pair_matrix();
    S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);

    seqlen = length;
    seq_= seq;
    


    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    index = new cand_pos_t [length];
    cand_pos_t total_length = (length *(length+1))/2;
    index[0] = 0;
    for (cand_pos_t i=1; i < length; i++)
        index[i] = index[i-1]+length-i+1;

    // this array holds V(i,j), and what (i,j) encloses: hairpin loop, stack pair, internal loop or multi-loop
    nodes = new free_energy_node [total_length];
    if (nodes == NULL) giveup ("Cannot allocate memory", "s_energy_matrix");
}


s_energy_matrix::~s_energy_matrix ()
// The destructor
{
    delete [] index;
    delete [] nodes;
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
* @brief Computes the internal loop value for V
* 
* @param V V array
* @param i row index
* @param j column index
*/
energy_t s_energy_matrix::compute_internal(cand_pos_t i, cand_pos_t j, const paramT *params){
	energy_t v_iloop = INF;
	cand_pos_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
	const int ptype_closing = pair[S_[i]][S_[j]];
	for ( int k=i+1; k<=max_k; ++k) {
		
		cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
        for (cand_pos_t l=j-1; l>=min_l; --l) {
            energy_t v_iloop_kl = E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<paramT *>(params)) + get_energy(k-1,l-1);
            v_iloop = std::min(v_iloop,v_iloop_kl);
            
        }
        
	}
	return v_iloop;
}

/**
 * @brief restricted version
*/
energy_t s_energy_matrix::compute_internal_restricted(cand_pos_t i, cand_pos_t j, const paramT *params, str_features *fres){
	energy_t v_iloop = INF;
	cand_pos_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
	const int ptype_closing = pair[S_[i]][S_[j]];
	for ( cand_pos_t k=i+1; k<=max_k; ++k) {
		
		cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
        if (!exists_restricted (i-1,k-1,fres)){
            for (int l=j-1; l>=min_l; --l) {
                if (!exists_restricted (l-1,j-1,fres)){
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


void s_energy_matrix::compute_energy (cand_pos_t i, cand_pos_t j)
// compute the V(i,j) value
{
    PARAMTYPE min, min_en[4];
    int k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;
    
    int ptype_closing = pair[S_[i+1]][S_[j+1]];
    if (ptype_closing>0)    // if i and j can pair
    {
        // compute free energy of hairpin loop, stack pair, internal loop and multi-loop
        min_en[0] = HairpinE(seq_,S_,S1_,params_,i+1,j+1);
        if (i<=j-TURN-1)
        {
            // TODO: uncomment
            if (!ignore_internal)
                min_en[1] = compute_internal(i+1,j+1,params_);
                // min_en[1] = compute_internal(i,j,params_);
            if (!ignore_multi)
                min_en[2] = VM->compute_energy (i, j);
        }
    }

    // see which of them is the minimum
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

    if (min_rank > -1)
    {
        if (debug)
            printf ("V(%d,%d) type %c energy %d\n", i, j, type, min);
    }

    if (min < INF/2)
    {
        int ij = index[i]+j-i;
        nodes[ij].energy = min;
        nodes[ij].type = type;
    }
}


void s_energy_matrix::compute_energy_restricted (cand_pos_t i, cand_pos_t j, str_features *fres)
// compute the V(i,j) value, if the structure must be restricted
{
    PARAMTYPE min, min_en[3];
    int k, min_rank;
    char type;

    min_rank = -1;
    min = INF/2;
    min_en[0] = INF;
    min_en[1] = INF;
    min_en[2] = INF;

    // printf("random thing is %d",params_);
	// Hosna, March 26, 2012
	// if the restricted base pairs are non-canonical then checking for can_pair only will cause missing those base pairs
	int ptype_closing = pair[S_[i+1]][S_[j+1]];
    if (ptype_closing>0 || fres[i].pair == j)    // if i and j can pair
    {

        if (!exists_restricted (i, j, fres)) min_en[0] = HairpinE(seq_,S_,S1_,params_,i+1,j+1);

        min_en[1] = compute_internal_restricted(i+1,j+1,params_,fres);
        min_en[2] = VM->compute_energy_restricted (i, j, fres);
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

    if (min_rank > -1 && debug) {
    	printf ("V(%d,%d) type %c energy %d\n", i, j, type, min);
   	}

    if (min < INF/2) {
        int ij = index[i]+j-i;
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
        
    //printf ("V(%d,%d) is_stack: %d energy %d\n", i, j, is_stack, energy);
    cand_pos_t ij = index[i]+j-i;
    nodes[ij].energy = energy;
    return;
}
