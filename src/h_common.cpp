#include "h_common.h"
#include "constants.h"
#include "h_externs.h"
#include "h_struct.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


/*
 * This function is just the same as detect_original_pairs
 * but it also calculates the arcs for each base and saves them in arc_table
 *
 * The algorithm for finding the arcs is as follows:
 * When checking each base in the input structure,
 *
 * case 1) if (structure[i] == '.' or ' ' or '_')
 * 		1-a) if stack is not empty, then put arc_table[i] = top element on the stack
 * 		1-b) otherwise arc_table[i] = -1
 *
 * case 2) if (structure[i] == '(')
 * 		2-a) if stack is not empty, then put arc_table[i] = top element on the stack and push i in the stack
 * 		2-b) otherwise arc_table[i] = -1 and push i in the stack
 *
 * case 3) if (strcuture[i] == ')'), pop the stack and match it with j
 * 		3-a) if stack is not empty, then put arc_table[i] = top element on the stack
 * 		3-b) otherwise put arc_table[i] = -1
 *
 */

void detect_original_pairs_arcs(char *structure, int *p_table, int *arc_table)
// PRE:  structure contains the desired structure
// POST: pairs will contain the index of each base pair
//               or -1 if it does not pair
//
{
		//printf("structure: %s\n", structure);
        int i, j, struct_len;
        stack_ds st;
        h_init (&st);
		//printf("The given structure is: \n %s\n", structure);
        struct_len = strlen (structure);
	// Hosna March 8, 2012
	// since index i starts at 0 and stack top is also set to 0 to show stack is empty, if i=0 is paired then we have incorrect arc values!
	// So I am introducing STACK_EMPTY = -1 to h_common.h and change h_init and h_pop accordingly

        for (i=0; i < struct_len; i++)
          {
			  // Hosna March 8, 2012
			  // changing nested ifs to switch for optimality
			  switch (structure[i])
				{
					case 'x':
					{
					  p_table[i] = RESTRICTED_UNPAIR;
					  if (st.top > STACK_EMPTY){//0){
						arc_table[i] = st.elem[st.top];
					  }else{
						arc_table[i] = NOT_COVERED;
					  }
					}
						break;
					case '.':
					case '_':
					{
					  p_table[i] = FREE_TO_PAIR;
					  if (st.top > STACK_EMPTY){//0){;
						arc_table[i] = st.elem[st.top];
					  }else{
						arc_table[i] = NOT_COVERED;
					  }
					}
						break;
					case '(':
						{
							if (st.top > STACK_EMPTY){//0){
								arc_table[i] = st.elem[st.top];
							  }else{
								arc_table[i] = NOT_COVERED;
							  }
							  h_push (&st, i);
						}
						break;
					case ')':
					  {
						j = h_pop (&st);
						p_table[i] = j;
						p_table[j] = i;
						  if (st.top > STACK_EMPTY){//0){
							arc_table[i] = st.elem[st.top];
							}else{
								arc_table[i] = NOT_COVERED;
							}
					  }
						break;
			  }
          }
        if (st.top != STACK_EMPTY)//0)
        {
            fprintf(stderr,"The given structure is not valid: %d more left parentheses than right parentheses\n", st.top);
            exit (1);
        }

}

/**
 * This function is just like the above function except the case that
 * it can handle the pseudoknotted structures of at most density 2
 *
 * the density 2 structures can be presented with ( and [ in dot parenthesis format
 * so we need at most 2 stacks to keep track of the base pairings
 * we call the stacks st and st_brack for ( and [ respectively.
 *
 *
 */

void detect_original_PKed_pairs(char *structure, int *p_table)
// PRE:  structure contains the desired structure
// POST: pairs will contain the index of each base pair
//               or -1 if it does not pair
{
	int i, j, struct_len;
	stack_ds st; //stach used for (
    stack_ds st_brack; // stack used for [
	h_init (&st);
	h_init (&st_brack);
	struct_len = strlen (structure);
	for (i=0; i < struct_len; i++) {
		// Hosna March 8, 2012
		// changing nested ifs to switch for optimality
		switch (structure[i]) {
			case '.':
				p_table[i] = -1;
				break;

			case ' ':
			case '_':
				p_table[i] = -2;
				break;

			case '(':
				h_push (&st, i);
				break;

			case '[':
				h_push (&st_brack, i);
				break;

			case ')':
				j = h_pop (&st);
				p_table[i] = j;
				p_table[j] = i;
				break;

			case ']':
				j = h_pop (&st_brack);
				p_table[i] = j;
				p_table[j] = i;
				break;
		}
	}

	if (st.top != STACK_EMPTY || st_brack.top != STACK_EMPTY) { //0 || st_brack.top != 0)
    	fprintf(stderr,"The given structure is not valid: %d more left parenthesis than right parentheses\n", st.top);
        exit (1);
    }
}

/*
 * The algorithm for finding the weakly closed regions for a given pseudoknot free structure is as follows
 *
 * initialization:
 * 				open = -1;
 * 				weakly_closed array initialized to 0
 *
 *
 * case 1) if the region is of the form [i,i] and i is not paired,
 * then it is considered as a weakly closed region and we put weakly_closed[ii] = 1
 * otherwise we put open = i
 *
 * case 2) if we are considering region [i,j] where i!=j and we know region [i, j-1] is weakly closed
 *		2-a) if j is not paired, then [i,j] is also weakly closed and we put weakly_closed[ij] = 1
 * 		2-b) if j is paired and we have i <= bp(j) < j then [i,j] is weakly_closed and we put weakly_closed[ij] = 1
 * 		THIS CASE CAN NEVER HAPPERN ==> REMOVE FROM CODE
 * 		2-c) otherwise [i,j] is NOT weakly closed
 * 			if (open == -1) then open = j
 *
 * case 3) if we are considering region [i,j] where i!=j and we know region [i,j-1] is NOT weakly closed
 * 		3-a) if j is paired and bp(j) == open then [i,j] is weakly closed and we put weakly_closed[ij] = 1 and open = -1
 * 		3-b) otherwise the region is still NOT weakly closed
 *
 */

void detect_weakly_closed(h_str_features *fres, int *weakly_closed, int nb_nucleotides, int *index){
	int i,j;
	int open = -1;

	for (i = 0; i < nb_nucleotides; i++){
		open = -1;
		for (j = i; j < nb_nucleotides; j++){
			int ij = index[i]+j-i; // index[i]+j-i gives the index ij
			if (i == j ){
				if (fres[j].pair < 0){ // j is not paired
					weakly_closed[ij]= 1;
				}else{
					open = j;
				}
			}else if (weakly_closed[ij-1]){
				if (fres[j].pair < 0){ // j is not paired
					weakly_closed[ij] = 1;
				}else if (open == -1){
					open = j;
				}
			}else if(fres[j].pair >= 0 && open == fres[j].pair){
				weakly_closed[ij] = 1;
				open = -1;
			}
		}
	}
}
/*
 * Hosna: Feb 12, 2007
 * algorithm for finding if region [i,j] is an empty region
 * We define an empty region as follows:
 * region [i,j] is an empty region if for all k i<k<j, k is unpaired
 *
 * in our program we only need to know empty regions based on the original structure
 *
 * for every base i and j, such that i <= j
 * 1) if (i == j && i is not paired in G)
 * then region [i,i] is an empty region and
 * we put not_paired_all[i,j] = 1
 *
 * 2) else if j is not paired in G and we know that region[i,j-1] is an empty region
 * then region [i,j] is an empty region and we put
 * not_paired_all[i,j] = 1
 *
 * 3) otherwise region[i,j] is not an empty region and we put
 * not_paired_all[i,j] = 0
 *
 */

void detect_not_paired_all(h_str_features *fres, int *not_paired_all, int nb_nucleotides, int *index){
	int i, j;
	for(i = 0; i < nb_nucleotides; i++){
		for(j = i; j < nb_nucleotides; j++){
			int ij = index[i]+j-i;
			if (i == j){
				if (fres[i].pair < 0){
					not_paired_all[ij]=1;
				}else{
					not_paired_all[ij] = 0;
				}
			}else if (not_paired_all[ij-1] == 1 && fres[j].pair < 0){
				not_paired_all[ij] = 1;
			}
			else{
				not_paired_all[ij] = 0;
			}
		}
	}
}



/* Hosna:
 * The algorithm for finding b(i,l) and B(l,j) is as follows:
 *
 *  for every base l, if l is in an arc and l is not paired (i.e. f[l].arc != -1)
 *  and every i (0 < i < nb_nucleotides)
 *
 *  1) if (i <= arc(l)) then we are finding b(i,l):
 * 		1-a) temp = arc(arc(l))
 * 		1-b) if (temp == -1 || temp <= i) then we put
 * 			b(i,l) = b'(i,l) (i.e. border_bs[l][i] = arc(l))
 * 		1-c) otherwise
 * 				1-c-i) while (arc(temp) ! = -1 && arc(temp) > i)
 * 							temp = arc(temp)
 * 				1-c-ii) b(i,l) = temp and we put border_bs[l][i] = temp
 *
 *  2) if (i >= pair(arc(l))) then we are finding B(l,j):
 * 		1-a) temp = arc(arc(l))
 * 		1-b) if (temp == -1 || pair(temp) >= i) then we put
 * 			B(l,i) = B'(l,i) (i.e. border_bs[l][i] = pair(arc(l)))
 * 		1-c) otherwise
 * 				1-c-i) while (arc(temp) ! = -1 && pair(arc(temp)) < i)
 * 							temp = arc(temp)
 * 				1-c-ii) B(l,i) = pair(temp) and we put border_bs[l][i] = pair(temp)
 *
 *  3) if (i > arc(l) && i < pair(arc(l))) then border_bs[l][i] = -1
 */

void detect_border_bs(h_str_features *fres, int** border_bs, int nb_nucleotides){

	int l,i;
	for (l = 0; l < nb_nucleotides; l++){
		for (i = 0; i < nb_nucleotides; i++){
			int cover_l = fres[l].arc, pair_l=fres[l].pair; // Hosna March 8, 2012, using local varibales for optimality
			if (cover_l == -1 || pair_l >= 0){
				border_bs[l][i] = -2;
			}else{
				if (i <= cover_l){
					int temp = fres[cover_l].arc;
					if (temp == -1 || temp < i){ // Hosna: Jan 31, 2007: temp <= i changed to < to include i itself too
						border_bs[l][i] = cover_l;
					}else{
						while(fres[temp].arc != -1 && fres[temp].arc >= i){ // Hosna: Jan 31, 2007: fres[temp].arc > i changed to >= to include i itself too
							temp = fres[temp].arc;
						}
						border_bs[l][i] = temp;
					}
				}
				if (i >= fres[cover_l].pair){
					int temp = fres[cover_l].arc;
					if (temp == -1 || fres[temp].pair > i){ //Hosna: Jan 31, 2007: fres[temp].pair >= i changed to > to include i itself too
						border_bs[l][i] = fres[cover_l].pair;
					}else{
						while(fres[temp].arc != -1 && fres[fres[temp].arc].pair <= i){ //Hosna: Jan 31, 2007: fres[fres[temp].arc].pair < i changed to <= to include i itself too
							temp = fres[temp].arc;
						}
						border_bs[l][i] = fres[temp].pair;
					}
				}
				if (i > cover_l && i < fres[cover_l].pair){
					border_bs[l][i] = -1;
				}
			}
		}
	}
}

/* Hosna:
 * The algorithm for finding b'(i,l) and B'(l,j) is as follows:
 *
 *  for every base l, if l is in an arc and l is not paired (i.e. f[l].arc != -1)
 *  and every i (0 < i < nb_nucleotides)
 *
 * 	1) if (i < arc(l)) then b'(i,l) = arc(l)
 *  and we put the corresponding value in border_bps[l][i]
 *
 *  2) if (i > pair(arc(l))) then B'(l,i) = pair(arc(l))
 *  and we put the corresponding value in border_bps[l][i]
 *
 *  3) if (i >= arc(l) && i<= pair(arc(l))) then border[l][i] = -1
 *
 */

void detect_border_bps(h_str_features *fres, int** border_bps, int nb_nucleotides){

	int l, i;
	for (l = 0; l < nb_nucleotides; l++){
		for (i = 0; i < nb_nucleotides ; i++){
			int cover_l = fres[l].arc, pair_l=fres[l].pair; // Hosna March 8, 2012, using local varibales for optimality
			if (cover_l  == -1 || pair_l >= 0){
				border_bps[l][i] = -2;//INF;//-2;
			}else{
				if (i <= cover_l ){		//Hosna: Jan 31, 2007: < changed to <= to include i itself too
					border_bps[l][i] = cover_l ;
				}
				if (i >= fres[cover_l ].pair){ //Hosna: Jan 31, 2007: > changed to >= to include i itself too
					border_bps[l][i] = fres[fres[l].arc].pair;
				}
				if ( i > cover_l  && i < fres[cover_l].pair){ //Hosna: Jan 31, 2007: >= and <= changed to > and < to include i itself too
					border_bps[l][i] = -1;
				}
			}
		}
	}
}

void h_init (stack_ds *st)
// PRE:  None
// POST: Initialize the stack st
{
	st->top = STACK_EMPTY;//0;
}

void h_push (stack_ds *st, int el)
// PRE:  st is a valid stack
// POST: Push an element to the stack
{
	st->top = st->top +1;
    st->elem[st->top] = el;
}

int h_pop (stack_ds *st)
// PRE:  st is a valid stack, that is not empty
// POST: pop an element from the stack and return it
{
    if (st->top <= STACK_EMPTY)//0)
    {
        fprintf(stderr,"The given structure is not valid: more right parentheses than left parentheses\n");
        exit (1);
    }
    int result = st->elem[st->top];
    st->top = st->top -1 ;
    return result;
}

void detect_h_structure_features (char *structure, h_str_features *f)
// PRE:  None
// POST: The variable f is filled with structure features, i.e. what type of elementary structure
//       this base is closing (such as stacked pair, hairpin loop etc.)
{
    int num_branches, i, j;
    int p_table[MAXSLEN];
    int arc_table[MAXSLEN];
    int bri[MAX_BRANCHES];
    int nb_nucleotides;
    nb_nucleotides = strlen(structure);
    detect_original_pairs_arcs (structure, p_table, arc_table);
    for (i=0; i < nb_nucleotides; i++)
    {
		// Hosna, March 8, 2012
		// use local variables instead of getting array values all the time
		int i_pair = p_table[i];
		f[i].pair = i_pair;//p_table[i];
        f[i].arc = arc_table[i];


        if (i_pair>i)//p_table[i] > i)
        {
            f[i_pair].pair = i;
            // check if it is stacked pair
			int i_pair_plus1 = p_table[i+1];
            if (i_pair_plus1 == i_pair-1 && i_pair_plus1 > i+1)//(p_table[i+1] == p_table[i]-1 && p_table[i+1] > i+1)
            {
                f[i].type = INTER;
                f[i_pair].type = INTER;//f[p_table[i]].type = STACK;
                continue;
            }
            // check if it is hairpin, internal loop or multi-loop
            num_branches = 0;
            for (j=i+1; j<i_pair; j++)//(j=i+1; j < p_table[i]; j++)
            {
                if (p_table[j] > j)
                {
                    bri[num_branches] = j;
                    num_branches++;
                    j = p_table[j];
                }
            }
            if (num_branches == 0)  // hairpin
            {
                f[i].type = HAIRP;
                f[i_pair].type = HAIRP; //f[p_table[i]].type = HAIRP;
            }
            else if (num_branches == 1) // internal loop
            {
                f[i].type = INTER;
                f[i_pair].type = INTER; //f[p_table[i]].type = INTER;
                f[i].num_branches = 1;
                f[i].bri[0] = bri[0];
            }
            else    // multi loop
            {
                f[i].type = MULTI;
                f[i_pair].type = MULTI; //f[p_table[i]].type = MULTI;
                f[i].num_branches = num_branches;
                for (j=0; j < num_branches; j++)
                    f[i].bri[j] = bri[j];
            }
        }

    }
}

/*
 * Hosna: January 10, 2008
 * The following two functions are modified versions of
 * the functions found in simfold/src/common/common.cpp
 * the modifications are to make them work for density-2 structures
 *
 */

double compute_h_sensitivity (char *ref_structure, char *pred_structure)
// returns 0 if undefined (denominator is 0)
{
    int ptable_ref[MAXSLEN];
    int ptable_pred[MAXSLEN];
    int distance;
    int len, i;
    double sens;
    int num_correct_bp;
    int num_true_bp;

    len = strlen(ref_structure);
    detect_original_PKed_pairs (ref_structure, ptable_ref);
    detect_original_PKed_pairs (pred_structure, ptable_pred);
    num_correct_bp = 0;
    num_true_bp = 0;
    for (i=0; i < len; i++)
    {
        if (ptable_ref[i] > -1)    // paired base
        {
            num_true_bp++;
            if (ptable_pred[i] == ptable_ref[i])
                num_correct_bp++;
        }
    }
    if (num_true_bp == 0)
        return -1.0;
    sens = num_correct_bp*1.0/num_true_bp;
    return sens;
}


double compute_h_ppv (char *ref_structure, char *pred_structure)
// returns 0 if undefined (denominator is 0)
{
    int ptable_ref[MAXSLEN];
    int ptable_pred[MAXSLEN];
    int distance;
    int len, i;
    double ppv;
    int num_correct_bp;
    int num_pred_bp;

    len = strlen(ref_structure);
    detect_original_PKed_pairs (ref_structure, ptable_ref);
    detect_original_PKed_pairs (pred_structure, ptable_pred);
    num_correct_bp = 0;
    num_pred_bp = 0;
    for (i=0; i < len; i++)
    {
        if (ptable_ref[i] > -1 && ptable_pred[i] == ptable_ref[i])    // paired base
            num_correct_bp++;
        if (ptable_pred[i] > -1)    // paired base
            num_pred_bp++;
    }
    if (num_pred_bp == 0)
        return -1.0;
    ppv = num_correct_bp*1.0/num_pred_bp;
    return ppv;
}

int is_structured (int i, int j, char *structure)
// return 1 if structure has some parentheses between i and j inclusive
// return 0 otherwise
{
    int k;
    for (k=i; k <= j; k++)
    {
        if (structure[k] != '.')
        {
            return 1;
        }
    }
    return 0;
}

double compute_accuracy (char *ref_structure, char *pred_structure)
{
    int ptable_ref[MAXSLEN];
    int ptable_pred[MAXSLEN];
    int distance;
    int len, i;
    double accuracy;

    len = strlen(ref_structure);
    detect_original_pairs (ref_structure, ptable_ref);
    detect_original_pairs (pred_structure, ptable_pred);
    distance = 0;
    for (i=0; i < len; i++)
    {
        if (ptable_pred[i] != ptable_ref[i])
            distance ++;
    }
    accuracy = 1.0-(double)distance/len;
    return accuracy;
}

double compute_sensitivity (char *ref_structure, char *pred_structure)
// returns 0 if undefined (denominator is 0)
{
    int ptable_ref[MAXSLEN];
    int ptable_pred[MAXSLEN];
    int distance;
    int len, i;
    double sens;
    int num_correct_bp;
    int num_true_bp;

    len = strlen(ref_structure);
    detect_original_pairs (ref_structure, ptable_ref);
    detect_original_pairs (pred_structure, ptable_pred);
    num_correct_bp = 0;
    num_true_bp = 0;
    for (i=0; i < len; i++)
    {
        if (ptable_ref[i] > -1)    // paired base
        {
            num_true_bp++;
            if (ptable_pred[i] == ptable_ref[i])
                num_correct_bp++;
        }
    }
    if (num_true_bp == 0)
        return 0.0;
    sens = num_correct_bp*1.0/num_true_bp;
    return sens;
}


double compute_ppv (char *ref_structure, char *pred_structure)
// returns 0 if undefined (denominator is 0)
{
    int ptable_ref[MAXSLEN];
    int ptable_pred[MAXSLEN];
    int distance;
    int len, i;
    double ppv;
    int num_correct_bp;
    int num_pred_bp;

    len = strlen(ref_structure);
    detect_original_pairs (ref_structure, ptable_ref);
    detect_original_pairs (pred_structure, ptable_pred);
    num_correct_bp = 0;
    num_pred_bp = 0;
    for (i=0; i < len; i++)
    {
        if (ptable_ref[i] > -1 && ptable_pred[i] == ptable_ref[i])    // paired base
            num_correct_bp++;
        if (ptable_pred[i] > -1)    // paired base
            num_pred_bp++;
    }
    if (num_pred_bp == 0)
        return 0.0;
    ppv = num_correct_bp*1.0/num_pred_bp;
    return ppv;
}

void giveup (const char *string1,const char *string2)
// to add: variable nb of parameters, as in scanf, printf
{
    char temp[100];
    sprintf (temp, "%s %s", string1, string2);
    perror (temp);
    exit(1);
}

int nuc_to_int (char nucleotide)
// PRE:  nucleotide is 'A', 'C', 'G' or 'T'
// POST: Return 0 for A, 1 for C, 2 for G, 3 for T
{
    switch(nucleotide)
    {
        case 'a':
        case 'A': return A;
        // In the thermodynamic set, the experiments in Chen_Turner_2005 and Chen_Turner_2006b contain P. For the MODEL = SIMPLE, we consider this is an A.
        case 'p':
        case 'P': return A;
        case 'c':
        case 'C': return C;
        case 'g':
        case 'G': return G;
		case 'u':
        case 'U': return U;
		case 'x':
        case 'X': return X;
        default : return U;
    }
}

void init (stack_ds *st)
// PRE:  None
// POST: Initialize the stack st
{
        st->top = 0;
}

void push (stack_ds *st, int el)
// PRE:  st is a valid stack
// POST: Push an element to the stack
{
        st->elem[st->top++] = el;
}

int pop (stack_ds *st)
// PRE:  st is a valid stack, that is not empty
// POST: pop an element from the stack and return it
{
    if (st->top <= 0)
    {
        fprintf(stderr,"The given structure is not valid: more right parentheses than left parentheses\n");
        exit (1);
    }
        return st->elem[--st->top];
}



void detect_original_pairs(char *structure, int *p_table)
// PRE:  structure contains the desired structure
// POST: p_table will contain the index of each base pair
//               or -1 if it does not pair
// Feb 28, 2008: structure can also have:
//  - angles: < or >, which denote the ends of a pseudoknot. In that case, p_table would still be filled in the same way.
//      The assumption is that the <> pairs are always nested within parentheses.
//      That is, a structure like this (<)> is not possible.
//  - x, which denotes that I should ignore that part. p_table would be -3 in that case/
{
        int i, j, struct_len;
        stack_ds st;
        init (&st);
        struct_len = strlen (structure);
        for (i=0; i < struct_len; i++)
          {
            if (structure[i] == 'x')
              p_table[i] = -1;
            else if (structure[i] == '.' || structure[i] == '_')
              p_table[i] = -2;
            else if ((structure[i] == 'x') || (structure[i] == 'X'))
              p_table[i] = -3;
            else if (structure[i] == '(' || structure[i] == '<')
              push (&st, i);
            else if (structure[i] == ')' || structure[i] == '>')
              {
                j = pop (&st);
                p_table[i] = j;
                p_table[j] = i;
              }
          }
        if (st.top != 0)
        {
            fprintf(stderr,"The given structure is not valid: %d more left parentheses than right parentheses: %s\n", st.top, structure);
            exit (1);
        }
}


int valid_structure (int i, int j, char *structure)
// returns 1 if this structure is valid (i.e. complete), 0 if it's partial
{
    int k;
    stack_ds st;
    init (&st);
    for (k=i; k <= j; k++)
    {
        if (structure[k] == '(')
            push (&st, k);
        else if (structure[k] == ')')
        {
            if (st.top == 0)
                return 0;    // more right parentheses
            pop (&st);
        }
    }
    if (st.top != 0)
        return 0;    // more left parentheses
    return 1;
}


void detect_structure_features (char *structure, str_features *f)
// PRE:  None
// POST: The variable f is filled with structure features, i.e. what type of elementary structure
//       this base is closing (such as stacked pair, hairpin loop etc.)
// Modified on Feb 28, 2008
//  Added the case when structure can also have angles and x's.
{
    int num_branches, i, j;
    int p_table[MAXSLEN];
    int bri[MAX_BRANCHES];
    int nb_nucleotides;

    nb_nucleotides = strlen(structure);
    detect_original_pairs (structure, p_table);
    for (i=0; i < nb_nucleotides; i++)
    {
        f[i].pair = p_table[i];
        if (p_table[i] > i)
        {
            f[p_table[i]].pair = i;
            // this base pair might be angle brackets or parentheses.
            if (structure[i] == '<')
            {
                // first make sure the pair is also angle bracket
                if (structure[p_table[i]] != '>')
                {
                   fprintf(stderr,"ERROR! structure is not valid, position %d should be > and is %c\n%s\n", p_table[i], structure[p_table[i]], structure);
                    exit(1);
                }
                f[i].type = NONE;
                f[p_table[i]].type = NONE;

                continue;
            }
            // if we got here, it means the base pair was ()
            // just make sure the partner is )
            if (structure[p_table[i]] != ')')
            {
               fprintf(stderr,"ERROR! structure is not valid, position %d should be ) and is %c\n%s\n", p_table[i], structure[p_table[i]], structure);
                exit(1);
            }
            // check if it is hairpin, internal loop or multi-loop
            num_branches = 0;
            // one of the branches can be <xxx> and this is a multi-loop no matter what
            int is_multi_loop = 0;
            for (j=i+1; j < p_table[i]; j++)
            {
                if (p_table[j] > j)
                {
                    if (structure[j] == '<')
                        is_multi_loop = 1;
                    bri[num_branches] = j;
                    num_branches++;
                    j = p_table[j];
                }
            }
            if (num_branches == 0)  // hairpin
            {
                int ignore_hairpin = 0;
                // check whether this hairpin loop should be ignored
                for (j=i+1; j < p_table[i]; j++)
                {
                    if (structure[j] == 'X' || structure[j] == 'x')
                    {
                        ignore_hairpin = 1;
                        break;
                    }
                }
                if (ignore_hairpin)
                {
                    f[i].type = NONE;
                    f[p_table[i]].type = NONE;
                }
                else
                {
                    f[i].type = HAIRP;
                    f[p_table[i]].type = HAIRP;
                }
            }
            else if (num_branches == 1 && !is_multi_loop) // internal loop
            {
                // TODO: test if we have x's inside
                f[i].type = INTER;
                f[p_table[i]].type = INTER;
                f[i].num_branches = 1;
                f[i].bri[0] = bri[0];
            }
            else    // multi loop
            {
                // TODO: test if we have x's inside
                f[i].type = MULTI;
                f[p_table[i]].type = MULTI;
                f[i].num_branches = num_branches;
                for (j=0; j < num_branches; j++)
                    f[i].bri[j] = bri[j];
            }
        }
    }

    // Ian Wark August 8 2017
    // exists_restricted is a very time consuming function because it is called so many times.
    // We will then pre-compute it and store it to save time.

    // Allocate the vectors
    f->exists_restricted_arr.resize(nb_nucleotides);
    for (int i = 0; i < nb_nucleotides; ++i) {
        f->exists_restricted_arr[i].resize(nb_nucleotides);

        // compute whether there are any restricted base pairs between i and j
        for (int j = i; j < nb_nucleotides; ++j) {
            // default 0
            f->exists_restricted_arr[i][j] = 0;

            for (int k = i+1; k < j; ++k)
            {
                if (f[k].pair > -1) {
                    // if restricted base pair, set to 1
                    f->exists_restricted_arr[i][j] = 1;
                }
            }

        }
    }
}

int exists_restricted (int i, int j, str_features *fres)
{
    return fres->exists_restricted_arr[i][j];
}