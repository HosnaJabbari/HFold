#ifndef PSEUDO_LOOP_H_
#define PSEUDO_LOOP_H_
#include "base_types.hh"
#include "h_struct.h"
#include "h_common.h"
#include <stdio.h>
#include <string.h>
#include "V_final.h"
#include "VM_final.h"

class VM_final;
class V_final;
class pseudo_loop{

public:
	// constructor
	pseudo_loop(std::string seq, char* restricted, V_final *V, VM_final *VM, vrna_param_t *params);

	// destructor
	~pseudo_loop();

    void compute_energies(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

    // energy_t get_energy(cand_pos_t i, cand_pos_t j);
	// in order to be able to check the border values consistantly
	// I am adding these get functions
	energy_t get_WI(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

	energy_t get_VP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	energy_t get_WMB(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	energy_t get_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, sparse_tree &tree);
	energy_t get_WIP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

	energy_t get_VPP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

	// based on discussion with Anne, we changed WMB to case 2 and WMBP(containing the rest of the recurrences)
	energy_t get_WMBP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

    int is_weakly_closed(cand_pos_t i, cand_pos_t j);
    int is_empty_region(cand_pos_t i, cand_pos_t j);

    void back_track(char *structure, minimum_fold *f, seq_interval *cur_interval, sparse_tree &tree);

	// Hosna, May 1st, 2012
	// We need a specific back track function for pkonly case
	void back_track_pkonly(char *structure, minimum_fold *f, seq_interval *cur_interval);

    void set_stack_interval(seq_interval *stack_interval);
    seq_interval *get_stack_interval(){return stack_interval;}
    char *get_structure(){return structure;}
    minimum_fold *get_minimum_fold(){return f;}

private:

	cand_pos_t n;
	char *restricted;
	std::string seq;

    VM_final *VM;	        // multi loop object
    V_final *V;		        // the V object

	h_str_features *fres;
	seq_interval *stack_interval;
	char *structure;
	minimum_fold *f;
	vrna_param_t *params_;


	//Hosna
    energy_t *WI;				// the loop inside a pseudoknot (in general it looks like a W but is inside a pseudoknot)
    energy_t *VP;				// the loop corresponding to the pseudoknotted region of WMB
    energy_t *WMB;				// the main loop for pseudoloops and bands
	energy_t *WMBP; 				// the main loop to calculate WMB
	energy_t *WIP;				// the loop corresponding to WI'
    energy_t *VPP;				// the loop corresponding to VP'
    energy_t *BE;				// the loop corresponding to BE
    cand_pos_t *index;				// the array to keep the index of two dimensional arrays like WI and weakly_closed

	short *S_;
	short *S1_;

    // function to allocate space for the arrays
    void allocate_space();

    void compute_WI(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	// Hosna: This function is supposed to fill in the WI array

	void compute_VP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	// Hosna: this function is supposed to fill the VP array

	void compute_WMB(cand_pos_t i,cand_pos_t j, sparse_tree &tree);
	// Hosna: this function is supposed to fill the WMB array

	// based on discussion with Anne, we changed WMB to case 2 and WMBP(containing the rest of the recurrences)
	void compute_WMBP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	// this is the helper recurrence to fill the WMB array

	void compute_WIP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	// Hosna: this function is supposed to fill the WIP array

	void compute_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, sparse_tree &tree);
	// Hosna: this function is supposed to fill the BE array

	void compute_VPP(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	// Hosna: this function is supposed to fill the VPP array

	// Hosna Feb 8th, 2007:
	// I have to calculate the e_stP in a separate function
	energy_t get_e_stP(cand_pos_t i, cand_pos_t j);
	energy_t get_e_intP(cand_pos_t i,cand_pos_t ip, cand_pos_t jp, cand_pos_t j);
	energy_t compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, const paramT *params);

  	// Hosna: Feb 19th 2007
  	// used for backtracking
  	void insert_node (cand_pos_t i, cand_pos_t j, char type);//, seq_interval *stack_interval);

};
#endif /*PSEUDO_LOOP_H_*/
