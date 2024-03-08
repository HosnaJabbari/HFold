#ifndef H_STRUCT_H_
#define H_STRUCT_H_

#include "constants.h"
#include <vector>

// the data structure stored in the V array
typedef struct minimum_fold
{
    short int pair;
    char type;                  // type can be 'H', 'S', 'I', 'M'
    char filled;                // I think this is not used any more
    minimum_fold()
    {
        pair = -1;
        type = NONE;
        filled = 'N';
    }
} minimum_fold;

// another way to represent a structure. Used to measure free energy of a give structure etc.
typedef struct str_features
{
    short int pair;
    char type;                   // type can be 'H', 'S', 'I', 'M' etc
    short int num_branches;
    int bri[MAX_BRANCHES];      // the i of each branch

    // Ian Wark August 8 2017
    // exists restricted is a very common function
    // that is ultimately very time consuming
    // precompute the results and save in this array
    std::vector< std::vector<int> > exists_restricted_arr;

    str_features()
    {
        pair = -1;
        type = NONE;
        num_branches = 0;
    }
} str_features;



typedef struct
{
        int top;
        int elem[MAXSLEN];
} stack_ds;



// This node is used to keep the intervals that need to be further backtracked
struct seq_interval
{
  int i;
  int j;
  PARAMTYPE energy;                        // it is used
  char type;
  seq_interval* next;

  void copy (seq_interval *other)
  {
    other->i = i;
    other->j = j;
    other->energy = energy;
    other->type = type;
  }
};


struct struct_node
{
    minimum_fold* f;                    // an array
    seq_interval* intervals;            // M: a linked list
    PARAMTYPE bot_en;                         // not used?
    PARAMTYPE energy;                         // M: min energy of any structure starting with the partial structure so far
    char* structure;
    struct_node* previous;              // M: made doubly linked list, to be able to keep the size < limit
    struct_node* next;

    struct_node()
    {
        f = NULL;
        intervals = NULL;
        structure = NULL;
        previous = NULL;
        next = NULL;
    }
};


typedef struct seq_node
//class seq_node
{
    char* structure;
    seq_node* next;
}seq_node;


struct free_energy_node
{
    PARAMTYPE energy;
    char type;          // type may be: N (NONE), H (HAIRPIN), S (STACKED), I (INTERNAL), M (MULTI)
    free_energy_node()
    {
        energy = 10000; // INF
        type = NONE;
    }
};
// another way to represent a structure. Used to measure free energy of a give structure etc.
typedef struct h_str_features
{
    short int pair;
    char type;                   // type can be 'H', 'S', 'I', 'M' etc
    short int num_branches;
    int bri[MAX_BRANCHES];      // the i of each branch
    int arc;					// keeps the left base pair of the arc
    h_str_features()
    {
        pair = -1;
        type = NONE;
        num_branches = 0;
        arc = -1;
    }
} h_str_features;

#endif /*H_STRUCT_H_*/
