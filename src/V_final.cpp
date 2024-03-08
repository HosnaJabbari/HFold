#include "h_struct.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "V_final.h"

V_final::V_final(cand_pos_t n){
	index.resize(n);
    cand_pos_t total_length = (n *(n+1))/2;
    index[0] = 0;
    for (cand_pos_t i=1; i < n; i++)
        index[i] = index[i-1]+n-i+1;

    type = new int[total_length];
    if (type == NULL) giveup ("Cannot allocate memory", "V_final");
    for (cand_pos_t i = 0; i < total_length; i++) type[i] = -1;
//	printf("an object of V_final was successfully created! \n");


}

V_final::~V_final(){
	delete [] type;
}

void V_final::setloops(s_energy_matrix *v, VM_final *vm){
	this->v = v;
	this->vm = vm;

//	printf("V_final loops were successfully set! \n");
}

energy_t V_final::get_energy(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	// Hosna: June 28th, 2007 -- I believe this is unnessecary as we are just collecting the energies -- If any of these cases occur, the energies would already be INF -- Mateo 2024
	// if (i >= j || (fres[i].pair > -1 && fres[i].pair != j) || (fres[j].pair > -1 && fres[j].pair != i)){
	// 	return INF;
	// }

	energy_t v_energy = v->get_energy(i,j);
	/*if (i>=13 && j<=39 && v_energy<INF){
	printf("V_final: v_energy(%d,%d) = %d and type = %c\n", i,j,v_energy, get_type(i,j));
	}*/
	energy_t vm_energy = vm->get_energy(i,j,tree);
//	printf("V_final: vm_energy(%d,%d) = %d \n", i,j,vm_energy);
	cand_pos_t ij = index[i]+j-i;
	if (v_energy < vm_energy){
		type[ij] = 0;
	}else{
		type[ij] = 1;
	}
	return MIN(v_energy,vm_energy);
}

char V_final::get_type(cand_pos_t i, cand_pos_t j){
	cand_pos_t ij = index[i]+j-i;
	if (type[ij] == 0) // comes from v
	{
		return v->get_type(i,j);
	}
	return MULTI;
}