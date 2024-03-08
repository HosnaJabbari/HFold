#include "VM_final.h"
#include "h_externs.h"
#include <iostream>

VM_final::VM_final(std::string seq, cand_pos_t n, vrna_param_t *params)
{
	this->n = n;
	this->v = NULL;
	this->wmb = NULL;
    params_ = params;
    make_pair_matrix();
    S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);

    index.resize(n);    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    cand_pos_t total_length = (this->n *(this->n+1))/2;
    index[0] = 0;
    for (cand_pos_t i=1; i < this->n; i++)
        index[i] = index[i-1]+this->n-i+1;

    WM = new int [total_length];
    if (WM == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (cand_pos_t i=0; i < total_length; i++) WM[i] = INF;

    VM = new int [total_length];
    if (VM == NULL) giveup ("Cannot allocate memory", "VM_final");
    for (cand_pos_t i=0; i < total_length; i++) VM[i] = INF;

//    printf("an object of VM_final was successfully created! \n");
}

VM_final::~VM_final()
{
    delete [] WM;
    delete [] VM;
    free(S_);
    free(S1_);
}

void VM_final::compute_energy(cand_pos_t i, cand_pos_t j, str_features *fres, sparse_tree &tree){
	// Hosna June 26, 2007
	// I have to figure out how to calculate the energy here
    if(j-i+1<4) return;
	// here comes the copied part from simfold with all dangling energies
	int min = INF, tmp, k;
    int iplus1k;
    int kplus1jminus1;
    int iplus2k;
    int kplus1jminus2;
    for (k = i+2; k <= j-3; k++)
    {
        iplus1k = index[i+1] + k -i-1;
        kplus1jminus1 = index[k+1] + j-1 -k-1;
        iplus2k = index[i+2] + k -i-2;
        kplus1jminus2 = index[k+1] + j-2 -k-1;

        energy_t WM2ij = WM[iplus2k] + WM[kplus1jminus1];
        energy_t WM2ip1j = WM[iplus1k] + WM[kplus1jminus2];
        energy_t WM2ijm1 = WM[iplus1k] + WM[kplus1jminus2];
        energy_t WM2ip1jm1 = WM[iplus2k] + WM[kplus1jminus2];


        min = std::min(min,v->v->E_MbLoop(WM2ij,WM2ip1j,WM2ijm1,WM2ip1jm1,S_,params_,i+1,j+1,tree.tree));
    }

    min += params_->MLclosing;

	energy_t wmb_energy = this->wmb->get_WMB(i,j,tree) + a_penalty + PSM_penalty;
	cand_pos_t ij = index[i]+j-i;
	VM[ij] = std::min(min,wmb_energy);
	//printf("VM[%d,%d] = %d \n",i,j,VM[ij]);

}

energy_t VM_final::get_energy(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	cand_pos_t ij = index[i]+j-i;
	if (i >= j){
		return INF;
	}
    // An empty region is weakly closed by definition so removing
	if (tree.weakly_closed(i+1,j+1)){
		return VM[ij];
	}
	return INF;
}
 
/**
 *  PRE: simfold's WM matrix has been filled for i and j
 *  and now we need to fill in the WM matrix that hfold needs
 *
 */
void VM_final::WM_compute_energy(cand_pos_t i, cand_pos_t j, sparse_tree &tree){

	// int s_wm = s_vm->get_energy_WM(i,j);
    energy_t s_wm = v->v->get_energy_WM(i,j);

	// Hosna: July 5th, 2007
	// add a b_penalty to this case to match the previous cases
	energy_t wmb_energy = wmb->get_WMB(i,j,tree)+PSM_penalty+b_penalty;
	energy_t min = std::min(s_wm, wmb_energy);
	cand_pos_t ij = index[i]+j-i;
	this->WM[ij] = min;
//	printf("hfold's WM min = %d \n",min);
}



energy_t VM_final::get_energy_WM(cand_pos_t i, cand_pos_t j , sparse_tree &tree){
	if (i >= j || !tree.weakly_closed(i+1,j+1) ){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
//	printf("hfold's WM(%d,%d) = %d \n", i,j,WM[ij]);
	return this->WM[ij];

}
