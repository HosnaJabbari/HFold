#ifndef VM_FINAL_H_
#define VM_FINAL_H_

#include <stdio.h>
#include "h_common.h"
#include "h_struct.h"
#include "pseudo_loop.h"
#include "V_final.h"
#include <string>

class pseudo_loop;
class V_final; 

class VM_final{
public:
	VM_final(std::string seq, cand_pos_t n,vrna_param_t *params);
	~VM_final();
	//void set_WMB_matrix(pseudo_loop *WMB) {this->WMB = WMB; }
	void set_V_matrix (V_final *Vf) { 
		this->v = Vf;

//		printf("VM: set_V_matrix successful!\n"); 
	}
	void set_WMB_matrix(pseudo_loop *wmb){
		this->wmb = wmb;
//		printf("VM: set_WMB_matrix successful!\n");
	}
	void compute_energy(cand_pos_t i, cand_pos_t j, str_features *fres, sparse_tree &tree);
	energy_t get_energy(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	//int get_energy_pk_only(int i, int j, str_features *fres); //April 4, 2012
	void WM_compute_energy(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	//void WM_compute_energy_pkonly(int i, int j);
	void set_WM_matrix(int *m){this->WM = m;}
	energy_t get_energy_WM(cand_pos_t i, cand_pos_t j, sparse_tree &tree);
	//int get_energy_WM_pkonly(int i, int j);

	vrna_param_t *params_;
	short *S_;
	short *S1_;
	
protected:

    int *sequence;                 // the entire sequence for which we compute the energy. 
                                       //     Each base is converted into integer, because it's faster.
    cand_pos_t n;                    // sequence length

	
    //s_energy_matrix *V;            // a pointer to the free energy matrix V
    V_final *v;
        
    pseudo_loop *wmb;				// a pointer to the pseudo_loop matrix
    
    std::vector<cand_pos_t> index;    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    int *WM;      // WM - 2D array (actually n*(n-1)/2 long 1D array)
    int *VM;
};

#endif /*VM_FINAL_H_*/
