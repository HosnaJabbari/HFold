/***************************************************************************
                          s_energy_matrix.h  -  description
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

// the V matrix

#ifndef ENERGY_MATRIX_H
#define ENERGY_MATRIX_H

#include "s_internal_loop.h"
#include "s_multi_loop.h"
#include "base_types.hh"
#include <string>


class s_energy_matrix
{
    public:

        friend class s_multi_loop;

        s_energy_matrix (std::string seq, cand_pos_t length, vrna_param_t *params);
        // The constructor

        ~s_energy_matrix ();
        // The destructor

        vrna_param_t *params_;

        void set_loops (s_multi_loop *VM)
        // Set the local loops to the given values
        {
            this->VM = VM;
        }

        short *S_;
        short *S1_;
        // VM_sub should be NULL if you don't want suboptimals

        void compute_energy (int i, int j);
        // compute the V(i,j) value

        void compute_energy_restricted (cand_pos_t i, cand_pos_t j, str_features *fres);


        free_energy_node* get_node (cand_pos_t i, cand_pos_t j) { cand_pos_t ij = index[i]+j-i; return &nodes[ij]; }
        // return the node at (i,j)

        // May 15, 2007. Added "if (i>=j) return INF;"  below. It was miscalculating the backtracked structure.
        PARAMTYPE get_energy (cand_pos_t i, cand_pos_t j) { if (i>=j) return INF; cand_pos_t ij = index[i]+j-i; return nodes[ij].energy; }
        // return the value at V(i,j)

        char get_type (cand_pos_t i, cand_pos_t j) { cand_pos_t ij = index[i]+j-i; return nodes[ij].type; }
        // return the type at V(i,j)
         //Mateo 13 Sept 2023
        void compute_hotspot_energy (cand_pos_t i, cand_pos_t j, bool is_stack);

        energy_t HairpinE(const std::string& seq, const short* S, const short* S1,  const paramT* params, cand_pos_t i, cand_pos_t j);
        energy_t compute_stack(cand_pos_t i, cand_pos_t j, const paramT *params);
        energy_t compute_internal(cand_pos_t i, cand_pos_t j, const paramT *params);
        energy_t compute_internal_restricted(cand_pos_t i, cand_pos_t j, const paramT *params, str_features *fres);
        energy_t compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, const paramT *params);


    // better to have protected variable rather than private, it's necessary for Hfold
    protected:
    //private:
        s_multi_loop *VM;

       
        std::string seq_;
        int seqlen;                // sequence length
        int *index;                // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
        free_energy_node *nodes;   // the free energy and type (i.e. base pair closing a hairpin loops, stacked pair etc), for each i and j
};



#endif
