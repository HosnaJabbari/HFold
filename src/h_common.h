#ifndef H_COMMON_H_
#define H_COMMON_H_

#include "h_struct.h"

#define P_WMB			'R'
#define P_VP			'D'
#define P_VPP			'E'
#define P_WI			'G'
#define P_BE			'J'
#define P_WIP			'L'
#define P_WMBP			'T'
#define P_V				'A' // This case is only for the cases that we have some pairings in input structure that are not valid in terms of simfold restrictions


#define NOT_COVERED		-1
#define STACK_EMPTY -1 // originally this value is 0, which I think is wrong! Hosna, March 8, 2012
#define RESTRICTED_UNPAIR -1
#define FREE_TO_PAIR	-2

void detect_original_pairs_arcs(char *structure, int *p_table, int *arc_table);
void detect_original_PKed_pairs(char *structure, int *p_table);
double compute_h_sensitivity (char *ref_structure, char *pred_structure);
double compute_h_ppv (char *ref_structure, char *pred_structure);

// Hosna: helper function to fill in the weakly_closed array
void detect_weakly_closed(h_str_features *fres, int *weakly_closed, int nb_nucleotides, int *index);

// Hosna: helper function to fill in the not_paired_all array
void detect_not_paired_all(h_str_features *fres, int *not_paired_all, int nb_nucleotides, int *index);

// Hosna: this function fills the bs table which keeps track of
// bs and Bs for each i and l
void detect_border_bs(h_str_features *fres, int** border_bs, int nb_nucleotides);

// Hosna: this function filld the bps table which keeps track of
// b' and B' for each l and j
void detect_border_bps(h_str_features *fres, int** border_bps, int nb_nucleotides);

void h_init (stack_ds *st);
void h_push (stack_ds *st, int el);
int h_pop (stack_ds *st);

void detect_h_structure_features (char *structure, h_str_features *f);

#define MIN(A, B)      ((A) < (B) ? (A) : (B))
#define MAX(A, B)      ((A) > (B) ? (A) : (B))

double compute_accuracy (char *ref_structure, char *pred_structure);
double compute_sensitivity (char *ref_structure, char *pred_structure);
double compute_ppv (char *ref_structure, char *pred_structure);

void giveup (const char *string1,const char *string2);
// to add: variable nb of parameters, as in scanf, printf

int nuc_to_int (char nucleotide);
// PRE:  nucleotide is 'A', 'C', 'G' or 'T'
// POST: Return 0 for A, 1 for C, 2 for G, 3 for T

void detect_original_pairs(char *structure, int *p_table);
// PRE:  structure contains the desired structure
// POST: p_table will contain the index of each base pair
//               or -1 if it does not pair
// Feb 28, 2008: structure can also have:
//  - angles: < or >, which denote the ends of a pseudoknot. In that case, p_table would still be filled in the same way.
//      The assumption is that the <> pairs are always nested within parentheses.
//      That is, a structure like this (<)> is not possible.
//  - x, which denotes that I should ignore that part. p_table would be -3 in that case/

int valid_structure (int i, int j, char *structure);
// returns 1 if this structure is valid (i.e. complete), 0 if it's partial

void detect_structure_features (char *structure, str_features *f);
// PRE:  None
// POST: The variable f is filled with structure features, i.e. what type of elementary structure
//       this base is closing (such as stacked pair, hairpin loop etc.)

int exists_restricted (int i, int j, str_features *fres);

int is_structured (int i, int j, char *structure);
// return 1 if structure has some parentheses between i and j inclusive
// return 0 otherwise

#endif /*H_COMMON_H_*/
