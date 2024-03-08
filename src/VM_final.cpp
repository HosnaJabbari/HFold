#include "VM_final.h"
#include "externs.h"
#include "h_externs.h"
#include <iostream>

VM_final::VM_final(int *seq, cand_pos_t n)
{
	this->n = n;
	sequence = seq;
	this->v = NULL;
	this->wmb = NULL;

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
}

void VM_final::compute_energy(cand_pos_t i, cand_pos_t j, str_features *fres, sparse_tree &tree){
	// Hosna June 26, 2007
	// I have to figure out how to calculate the energy here

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


        tmp = WM[iplus1k] + WM[kplus1jminus1];
        if (tmp < min)
            min = tmp;

        if (fres[i+1].pair <= -1)
        {
            tmp = WM[iplus2k] + WM[kplus1jminus1] +
                dangle_top [sequence [i]]
                [sequence [j]]
                [sequence [i+1]] +
                misc.multi_free_base_penalty;
            if (tmp < min)
                min = tmp;
        }
        if (fres[j-1].pair <= -1)
        {
            tmp = WM[iplus1k] + WM[kplus1jminus2] +
                dangle_bot [sequence[i]]
                [sequence[j]]
                [sequence[j-1]] +
                misc.multi_free_base_penalty;
            if (tmp < min)
                min = tmp;
        }
        if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
        {
            tmp = WM[iplus2k] + WM[kplus1jminus2] +
                dangle_top [sequence [i]]
                [sequence [j]]
                [sequence [i+1]] +
                dangle_bot [sequence[i]]
                [sequence[j]]
                [sequence[j-1]] +
                2 * misc.multi_free_base_penalty;
            if (tmp < min)
            min = tmp;
        }
    }

    min += misc.multi_helix_penalty + misc.multi_offset +
           AU_penalty (sequence[i], sequence[j]);

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
