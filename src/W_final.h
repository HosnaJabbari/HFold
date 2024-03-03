#ifndef W_FINAL_H_
#define W_FINAL_H_

#include "hotspot.hh"
#include "VM_final.h"
#include "V_final.h"
#include "pseudo_loop.h"
#include "base_types.hh"
// #include "s_min_folding.h"
#include "s_energy_matrix.h"
#include "s_multi_loop.h"
#include "h_common.h"
#include <string>

void get_hotspots(std::string seq,std::vector<Hotspot> &hotspot_list, int max_hotspot, vrna_param_s *params);
int distance(int left, int right);
void expand_hotspot(s_energy_matrix *V, Hotspot &hotspot, int n);
//Mateo 2024
//comparison function for hotspot so we can use it when sorting
bool compare_hotspot_ptr(Hotspot &a, Hotspot &b);
 

class W_final{
	public:
		W_final(std::string seq, char *cseq, char *res, bool pk_free);
        // constructor for the restricted mfe case

        ~W_final ();
        // The destructor

        double hfold ();

        vrna_param_t *params_;
        std::string structure;        // MFE structure
        // PRE:  the init_data function has been called;
        //       the space for structure has been allocate
        // POST: fold sequence, return the MFE structure in structure, and return the MFE

		// PRE:  the init_data function has been called;
		//       the space for structure has been allocate
		// POST: fold sequence, return the MFE structure in structure, and return the MFE



    protected:
    	// Hosna: June 18th, 2007:
        // this pointer is the main part of the Hierarchical fold program
        // and corresponds to WMB recurrence
        pseudo_loop *WMB;
        // pointer to the final V matrix
        V_final *v;

        s_multi_loop *VM;       // multi loop object
        s_energy_matrix *V;     // the V object
        PARAMTYPE *W;                 // the W exterior loop array
        int n;     // sequence length (number of nucleotides)
        int* int_sequence;      // sequence in integer representation (faster)  
        seq_interval *stack_interval;  // used for backtracking
        minimum_fold *f;        // the minimum folding, see structs.h
        std::string seq_;
        short *S_;
	    short *S1_;
        char *restricted;    // restricted structure given as input - restricts base pairs eg (________) 
        char* sequence;
        bool pk_free = false;
        

        // pointer to the final VM matrix
        VM_final *vm;

        void insert_node (int i, int j, char type);

        void space_allocation();

        // allocate the necessary memory
        double fold_sequence_restricted ();

        void backtrack_restricted (seq_interval *cur_interval, str_features *fres);
        // backtrack, the restricted case

		// backtrack, the restricted case with pk only base pairs

        void compute_W_restricted (int j, str_features *fres);
        // fill the W array, the restricted case

		// fill the W array, with addition of just pseudoknotted base pairs to the original structure

        int compute_W_br2_restricted (int j, str_features *fres, int &must_choose_this_branch);


        int compute_W_br3_restricted (int j, str_features *fres);

        energy_t E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, paramT* params, const cand_pos_t i,const cand_pos_t j, cand_pos_t n, str_features *fres);

};

#endif /*W_FINAL_H_*/
