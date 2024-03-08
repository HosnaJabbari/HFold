#ifndef H_COMMON_H_
#define H_COMMON_H_

#include "h_struct.h"


//Hosna: June 22, 2007

// #define PS_penalty = 960 		//exterior pseudoloop initiation penalty (9.6 Kcal/mol)
// #define PSM_penalty 		1500		//penalty for introducing pseudoknot inside a multiloop (15 Kcal/mol)
// #define PSP_penalty 		1500		//penalty for introducing pseudoknot inside a pseudoloop (15 Kcal/mol)
// #define PB_penalty 			20			//band penalty (0.2 Kcal/mol)
// #define PUP_penalty			10			//penalty for an un-paired base in a pseudoloop or a band (0.1 Kcal/mol)
// #define PPS_penalty 		10			//penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band(0.1 Kcal/mol)

// #define a_penalty			340			//penalty for introducing a multiloop (3.4 Kcal/mol)
// #define b_penalty			40			//penalty for base pair in a multiloop (0.4 Kcal/mol)
// #define c_penalty			0			//penalty for un-paired base in a multi-loop

// #define e_stP_penalty		0.83		// e_stP = 0.83 * e_s
// #define e_intP_penalty		0.83		// e_intP = 0.83 * e_int

// #define ap_penalty			340			//penalty for introducing a multiloop that spans a band (3.4 Kcal/mol)
// #define bp_penalty			40			//base pair penalty for a multiloop that spans a band (0.4 Kcal/mol)
// #define cp_penalty          0			//penalty for unpaired base in a multiloop that spans a band



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

// Hosna, March 19, 2012
// the original value of the matrices should be set to -INF and then be changed to their correct value
#define MINUS_INF             -1600000      // a very small value (minus infinity)

void detect_original_pairs_arcs(char *structure, int *p_table, int *arc_table);
void detect_original_PKed_pairs(char *structure, int *p_table);
double compute_h_sensitivity (char *ref_structure, char *pred_structure);
double compute_h_ppv (char *ref_structure, char *pred_structure);
// void h_fill_data_structures_with_new_parameters (char *filename);
// void h_fill_data_structures_with_new_parameters (double *param);
// int h_create_string_params();

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

// Hosna June 26, 2007
// I need functions to convert str_features to h_str_features and the other way around

h_str_features *convert_str_features_to_h_str_features(str_features *f);
str_features *convert_h_str_features_to_str_features(h_str_features *f);

void detect_h_structure_features (char *structure, h_str_features *f);



#define isY(i)  (i==U || i==C)
#define isR(i)  (i==A || i==G)

#define MIN(A, B)      ((A) < (B) ? (A) : (B))
#define MAX(A, B)      ((A) > (B) ? (A) : (B))

// I precomputed these values in s_partition_function.cpp
//#define EXPA   (exp (misc.multi_offset * oneoverRT))
//#define EXPB(X)   (exp (((X)*misc.multi_helix_penalty) * oneoverRT))
//#define EXPC(X)   (exp (((X)*misc.multi_free_base_penalty) * oneoverRT))

#define IFD    if (ignore_dangles)

#define AU_penalty(X,Y)  ((((X) != C || (Y) != G) && ((X) != G || (Y) != C))?misc.terminal_AU_penalty:0)
#define has_AU_penalty(X,Y)  ((((X) != C || (Y) != G) && ((X) != G || (Y) != C))?1:0)
#define AU_penalty_enthalpy(X,Y)  ((((X) != C || (Y) != G) && ((X) != G || (Y) != C))?enthalpy_misc.terminal_AU_penalty:0)

//#define asymmetry_penalty(size1, size2) (MIN (misc.asymmetry_penalty_max_correction, abs (size1-size2) * misc.asymmetry_penalty_array [MIN (2, MIN ((size1), (size2)))-1]))

// PARAMTYPE asymmetry_penalty (int size1, int size2);

// #define asymmetry_penalty_enthalpy(size1, size2) (MIN (enthalpy_misc.asymmetry_penalty_max_correction, abs (size1-size2) * enthalpy_misc.asymmetry_penalty_array [MIN (2, MIN ((size1), (size2)))-1]))

#define IGINF(x) (((x) == INF)?0:(x))
// ignore infinite values

// does not modify numbers


double compute_accuracy (char *ref_structure, char *pred_structure);
double compute_sensitivity (char *ref_structure, char *pred_structure);
double compute_ppv (char *ref_structure, char *pred_structure);

// double compute_pf_sensitivity (char *ref_structure, s_partition_function *part, double threshold);
// compute the sensitivity obtained after thresholding the base pair probabilities
// part is the partition function object, which contains base pair probabilities
// returns -1 if undefined (denominator is 0)

// double compute_pf_ppv (char *ref_structure, s_partition_function *part, double threshold);
// compute the positive predictive value obtained after thresholding the base pair probabilities
// part is the partition function object, which contains base pair probabilities
// returns -1 if undefined (denominator is 0)

void giveup (const char *string1,const char *string2);
// to add: variable nb of parameters, as in scanf, printf

void giveup2 (const char *string1,const char *string2, FILE *file);
// to add: variable nb of parameters, as in scanf, printf

void create_random_sequence (int length, char *sequence);
// function to create uniformly random sequences - for demonstration purposes

void create_random_restricted (char *sequence, char *restricted);
// sequence is an input argument
// restricted is the output argument

void remove_space (char *structure);
// PRE: none
// POST: remove the space(s) from structure, if any; modifies structure

void empty_string (char * str);

int can_pair (int base1, int base2);
// PRE:  base1 and base2 are nucleotides over the alphabet {A, C, G, T, U}
// POST: return 1 if they can pair, 0 otherwise

int watson_crick (int base1, int base2);
// PRE:  base1 and base2 are nucleotides over the alphabet {A, C, G, T, U}
// POST: return 1 if they are watson crick pair, 0 otherwise

int nuc_to_int (char nucleotide);
// PRE:  nucleotide is 'A', 'C', 'G' or 'T'
// POST: Return 0 for A, 1 for C, 2 for G, 3 for T


char int_to_nuc (int inuc);


int is_nucleotide (char base);
// PRE:  base is a character
// POST: return true if base is a nucleotide (A, C, G or T)
//       return false otherwise


void check_sequence (char *sequence);
// check sequence for length and alphabet


// PARAMTYPE penalty_by_size (int size, char type);
// PRE:  size is the size of the loop
//       type is HAIRP or INTER or BULGE
// POST: return the penalty by size of the loop

//PARAMTYPE IL_penalty_by_size_2D (int size1, int size2);

// PARAMTYPE penalty_by_size_enthalpy (int size, char type);

void substr (char *source, int begin, int end, char *dest);
// PRE:  begin and end are smaller than strlen(source)
// POST: Put in dest what is in source between position begin and position end

void replace_str_piece (char *sequence, int position, char *seq);
// PRE:  begin + strlen(seq) < strlen (sequence)
// POST: In sequence, at position, replace what is was by seq

void reverse_complement_of_seq (const char *seq, char *complem);
// PRE:  seq is a sequence
// POST: complement and reverse sequence and put the result into compl

void insert_space (char *structure, int place);
// PRE:  None
// POST: insert a space at the specified place, in structure

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

int complementary_bases (char b1, char b2);
// returns 1 if b1 and b2 are complementary bases

int self_complementary (char *sequence);
// return 1 if this sequence is self-complementary
// self_complementary means the first half is the reverse complement of the second half
// if length (sequence) is an odd number, the middle base does not matter

int exists_restricted (int i, int j, str_features *fres);
int exists_restricted_ptable (int i, int j, int *ptable);

int is_structured (int i, int j, char *structure);
// return 1 if structure has some parentheses between i and j inclusive
// return 0 otherwise

// void print_stacking_energies();
// void print_tstacki_dangling_energies();
// void print_stack_dangling_energies();
// void print_stack_equation_dangling();
// void print_int22_tstacki();

// void read_parsi_options_from_file (char *filename);
// the file should contain values (0 or 1) for each of the parsi options. For example the following is all-parsimonious options
//     parsi_tstackh = 1;
//     parsi_tstacki = 1;
//     parsi_asymmetry = 1;
//     parsi_int11 = 1;
//     parsi_int21 = 1;
//     parsi_int22 = 1;
//     parsi_bulge1 = 1;
//     parsi_dangles = 1;
//     parsi_others = 1;
//     parsi_length = 1;
//     parsi_special = 1;


int loss (int first, int last);
// known_pairings contains the pairings from 0 to n-1, of the reference structure
// pred_pairings contains the pairings of a potential structure on the region first-last inclusive
//      the other regions don't matter
// Returns the "Hamming" distance between known_pairings and pred_pairings on the region first-last inclusive
// This function is used for the loss-augmented prediction
// Written on August 9, 2008
// Note: Maybe this measure is better than the Hamming distance:
//      (# correctly predicted bp - # incorrectly predicted bp) / # true bp.
//      This will be in (-inf,1], but it only includes the base pairs,
//      whereas the Hamming distance measure also includes the unpaired bases.

// Added on Sep 3, 2008, for loss-augmented prediction
double compute_distance (char *ref_structure, char *pred_structure);
// It has to be the same mathematical function as the one implemented in "loss"

#endif /*H_COMMON_H_*/
