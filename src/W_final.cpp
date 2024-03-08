
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "pseudo_loop.h"
#include "V_final.h"
#include "W_final.h"
#include "constants.h"
#include "h_struct.h"
#include "externs.h"
#include "h_externs.h"
#include "h_common.h"


// Hosna June 20th, 2007
// calls the constructor for s_min_folding
// to create all the matrixes required for simfold
// and then calls allocate_space in here to allocate
// space for WMB and V_final
W_final::W_final(std::string seq,std::string res,char *cseq, char *restricted, bool pk_free) : params_(scale_parameters())
{
	seq_ = seq;
	this->sequence = cseq;
	this->restricted = restricted;
	this->res = res;
	this->n = seq.length();
	make_pair_matrix();
	params_->model_details.dangles = 1;
    S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);
	this->pk_free = pk_free;
	W.resize(n,0);
	space_allocation();
}


W_final::~W_final()
{
	delete vm;
	delete v;
	delete WMB;
	delete V;
	delete [] f;
    delete [] int_sequence;
	free(params_);
	free(S_);
	free(S1_);

	// Ian Wark July 21 2017
	// we don't need this, this is done in s_min_folding
	//delete [] int_sequence;
}

// Hosna June 20th, 2007
// allocates space for WMB object and V_final
void W_final::space_allocation(){

	// From simfold
	f = new minimum_fold [n];
    if (f == NULL) giveup ("Cannot allocate memory", "energy");
	int_sequence = new int[n];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (int i=0; i < n; i++) int_sequence[i] = nuc_to_int(sequence[i]);

    V = new s_energy_matrix (seq_, n,params_);
    if (V == NULL) giveup ("Cannot allocate memory", "energy");
	structure = std::string (n,'.');

		// Hosna June 20th, 2007
	vm = new VM_final(this->int_sequence,this->n);
	if (vm == NULL) giveup ("Cannot allocate memory", "W_final");
	// Hosna June 20th, 2007
	v = new V_final(n);
	if (v == NULL) giveup ("Cannot allocate memory", "W_final");
	v->setloops(this->V,vm);

	// Hosna: June 20th 2007
    WMB = new pseudo_loop (seq_,sequence,restricted,v,vm,params_);
    if (WMB == NULL) giveup ("Cannot allocate memory", "W_final");

    // Hosna: June 20th 2007
    vm->set_V_matrix(v);
    vm->set_WMB_matrix(WMB);


}




double W_final::hfold(sparse_tree &tree){

	h_str_features *h_fres;
    if ((h_fres = new h_str_features[n]) == NULL) giveup ("Cannot allocate memory", "h_str_features");
    // detect the structure features
    detect_h_structure_features (restricted, h_fres);
    WMB->set_features(h_fres);
    WMB->initialize();

    str_features *fres;
    if ((fres = new str_features[n]) == NULL) giveup ("Cannot allocate memory", "str_features");
    // detect the structure features
    detect_structure_features (restricted, fres);

	// for (int j=1; j <= n; ++j)
    // {
    //     for (int i=1; i<j; ++i)
    //     {
	// 		const bool evaluate = tree.weakly_closed(i,j);
	// 		const pair_type ptype_closing = pair[S_[i]][S_[j]];
	// 		const bool restricted = tree.tree[i].pair == -1 || tree.tree[j].pair == -1;
	// 		// const bool pkonly = (!pk_only || paired);

	// 		if(ptype_closing> 0 && evaluate && !restricted)
    //         V->compute_energy_restricted (i,j,tree);


    //     }
    //     // if I put this before V calculation, WM(i,j) cannot be calculated, because it returns infinity
	// 	V->compute_energy_WM_restricted(j,tree.tree);
    // }

		for (int i = n; i >=1; --i)
		{	
			for (int j =i; j<=n; ++j)//for (i=0; i<=j; i++)
			{
				const bool evaluate = tree.weakly_closed(i,j);
				const pair_type ptype_closing = pair[S_[i]][S_[j]];
				const bool restricted = tree.tree[i].pair == -1 || tree.tree[j].pair == -1;
				// const bool pkonly = (!pk_only || paired);

				if(ptype_closing> 0 && evaluate && !restricted)
				V->compute_energy_restricted (i,j,tree);
				
				V->compute_energy_WM_restricted(i,j,tree.tree);

				WMB->compute_energies(i-1,j-1,tree);

				vm->WM_compute_energy(i-1,j-1,tree);		
			}

		}

	for (cand_pos_t j= TURN+1; j <= n; j++){
		energy_t m1 = INF;
		energy_t m2 = INF;
		energy_t m3 = INF;
		if(tree.tree[j-1].pair < 0) m1 = W[j-1-1];
		
		
		for (cand_pos_t i=1; i<=j-TURN-1; i++){
		 	// m2 = compute_W_br2_restricted (j, fres, must_choose_this_branch);
			energy_t acc = (i>1) ? W[i-1-1]: 0;
			m2 = std::min(m2,acc + E_ext_Stem(v->get_energy(i-1,j-1,tree),v->get_energy(i+1-1,j-1,tree),v->get_energy(i-1,j-1-1,tree),v->get_energy(i+1-1,j-1-1,tree),S_,params_,i,j,n,tree.tree));
			if(tree.weakly_closed(i,j)) m3 = std::min(m3,acc + WMB->get_WMB(i-1,j-1,tree) + PS_penalty);
			}
		W[j-1] = std::min({m1,m2,m3});
	}


	// I can't do this at the time at the moment as, unlike with candidates, there is no way for me get a partitionable index that is already calculated if I make indices like the above i.e. W[i] is not calculated yet and I can't add it to W[j]
	// for (cand_pos_t i = n-1; i>=0; --i){
	// 	energy_t acc = (i-1>0) ? W[i-1]: 0;
		
	// 	for (cand_pos_t j=i; j<n; ++j){
	// 		energy_t m1 = INF;
	// 		energy_t m2 = INF;
	// 		energy_t m3 = INF;
	// 	 	// m2 = compute_W_br2_restricted (j, fres, must_choose_this_branch);
	// 		if(fres[j].pair < 0) m1 = W[j-1];
	// 		m2 = std::min(m2,E_ext_Stem(v->get_energy(i,j),v->get_energy(i+1,j),v->get_energy(i,j-1),v->get_energy(i+1,j-1),S_,params_,i,j,n,fres) + acc);
	// 		m3 = std::min(m3,WMB->get_energy(i,j) + acc + PS_penalty);
	// 		W[j] = std::min(m1,std::min(m2,m3));
	// 	}
	// }
    double energy = W[n-1]/100.0;






    // // backtrack
    // // first add (0,n-1) on the stack
    stack_interval = new seq_interval;
    stack_interval->i = 0;
    stack_interval->j = n - 1;
    stack_interval->energy = W[n-1];
    stack_interval->type = FREE;
    stack_interval->next = NULL;

    seq_interval *cur_interval = stack_interval;

    while ( cur_interval != NULL)
    {
        stack_interval = stack_interval->next;
        backtrack_restricted (cur_interval, fres,tree);
        delete cur_interval;    // this should make up for the new in the insert_node
        cur_interval = stack_interval;
    }
    delete [] h_fres;
    delete [] fres;
    return energy;

}

/**
 * @brief Gives the W(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
 * Until the changes to fres, I am adding +1 to the ptype closing and Si and Sj's to make them match - Mateo 2024
 * 
 * @param vij The V(i,j) energy
 * @param vi1j The V(i+1,j) energy
 * @param vij1 The V(i,j-1) energy
 * @param vi1j1 The V(i+1,j-1) energy
*/
energy_t W_final::E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, paramT* params, const cand_pos_t i,const cand_pos_t j, cand_pos_t n, std::vector<Node> &tree){

	energy_t e = INF,en = INF;
  	pair_type tt  = pair[S[i]][S[j]];
	
    if ((tree[i].pair <-1 && tree[j].pair <-1) || (tree[i].pair == j && tree[j].pair == i)) {
				en = vij; // i j

				if (en != INF) {
					if (params->model_details.dangles == 2){
						base_type si1 = i>1 ? S[i-1] : -1;
                		base_type sj1 = j<n ? S[j+1] : -1;
                        en += vrna_E_ext_stem(tt, si1, sj1, params);
					}
                    else{
                        en += vrna_E_ext_stem(tt, -1, -1, params);
					}

                    e = MIN2(e, en);
					
				}

	}

	if(params->model_details.dangles  == 1){
        tt  = pair[S[i+1]][S[j]];
        if (((tree[i+1].pair <-1 && tree[j].pair <-1) || (tree[i+1].pair == j)) && tree[i].pair<0) {
            en = (j-i-1>TURN) ? vi1j : INF; //i+1 j

            if (en != INF) {

                base_type si1 = S[i];
                en += vrna_E_ext_stem(tt, si1, -1, params);
            }

            e = MIN2(e,en);

        }
        tt  = pair[S[i]][S[j-1]];
        if (((tree[i].pair <-1 && tree[j-1].pair <-1) || (tree[i].pair == j-1)) && tree[j].pair<0) {
            en = (j-1-i>TURN) ? vij1 : INF; // i j-1
            if (en != INF) {

                base_type sj1 = S[j];

                en += vrna_E_ext_stem(tt, -1, sj1, params);
            }
            e = MIN2(e,en);

        }
        tt  = pair[S[i+1]][S[j-1]];
        if (((tree[i+1].pair <-1 && tree[j-1].pair <-1) || (tree[i+1].pair == j-1)) && tree[i].pair < 0 && tree[j].pair<0) {
            en = (j-1-i-1>TURN) ? vi1j1 : INF; // i+1 j-1

            if (en != INF) {

                base_type si1 = S[i];
                base_type sj1 = S[j];

                en += vrna_E_ext_stem(tt, si1, sj1, params);
            }
            e = MIN2(e,en);
        }
	}
	return e;
}

void W_final::backtrack_restricted(seq_interval *cur_interval, str_features *fres, sparse_tree &tree){
    char type;


	// printf("type is %c and i is %d and j is %d\n",cur_interval->type,cur_interval->i,cur_interval->j);
	//Hosna, March 8, 2012
	// changing nested if to switch for optimality
	switch (cur_interval->type){
		case LOOP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (i >= j)
				return;
			f[i].pair = j;
			f[j].pair = i;

			// Hosna Jun. 28 2007
			// if the pairing is part of original structure, put '(' and ')' in the structure
			// otherwise make it '[' and ']' -- changed to () if pseudoknot-free and [] if pseudoknotted -Mateo
			structure[i] = '(';
			structure[j] = ')';
		

			type = v->get_type (i,j);
			if(i==4 && j==32) printf(" i is %d and j is %d and type is %c\n",i,j,type);
			
			if (debug)
				printf ("\t(%d,%d) LOOP - type %c\n", i,j,type);
			// Hosna, March 8, 2012
			// changing nested ifs to switch for optimality
			switch (type){
				case HAIRP:
			//else if (type == HAIRP)
				{
					f[i].type = HAIRP;
					f[j].type = HAIRP;
				}
					break;
				case INTER:
			//else if (type == INTER)
				{
					f[i].type = INTER;
					f[j].type = INTER;
					// detect the other closing pair
					int ip, jp, best_ip, best_jp, minq;
					int tmp, min = INF;
					int max_ip = std::min(j-TURN-2,i+MAXLOOP+1);
					for (ip = i+1; ip <= max_ip; ip++)
					{
						if (!exists_restricted (i,ip,fres)){
							int min_l=std::max(ip+TURN+1 + MAXLOOP+2, ip+j-i) - MAXLOOP-2;
							for (int jp = j-1; jp >= min_l; --jp)
							{
								if(!exists_restricted (jp,j,fres)){
							
									tmp = V->compute_int(i+1,j+1,ip+1,jp+1,params_);
									if (tmp < min)
									{
										min = tmp;
										best_ip = ip;
										best_jp = jp;
									}
								}
							}
						}
					}
					if (best_ip < best_jp)
						insert_node (best_ip, best_jp, LOOP);
					else
					{
						fprintf(stderr,"NOT GOOD RESTR INTER, i=%d, j=%d, best_ip=%d, best_jp=%d\n", i, j, best_ip, best_jp);
						exit (0);
					}
				}
					break;
				case MULTI:
			//else if (type == MULTI)
				{
					f[i].type = MULTI;
					f[j].type = MULTI;
					int k, best_k = -1, best_row = -1;
					int tmp= INF, min = INF;
					for (k = i+1; k <= j-1; k++)
					  {
						tmp = vm->get_energy_WM (i+1,k,tree) + vm->get_energy_WM (k+1, j-1,tree);
						if (tmp < min)
						  {
							min = tmp;
							best_k = k;
							best_row = 1;
						  }
						  // TODO:
						  // Hosna, May 1st, 2012
						  // do I need to check for non-canonical base pairings here as well so the dangle values not be INF??
						if (fres[i+1].pair <= -1)
						{
							tmp = vm->get_energy_WM (i+2,k,tree) + vm->get_energy_WM (k+1, j-1,tree) +
							dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
							misc.multi_free_base_penalty;
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 2;
							}
						}
						if (fres[j-1].pair <= -1)
						{
							tmp = vm->get_energy_WM (i+1,k,tree) + vm->get_energy_WM (k+1, j-2,tree) +
							dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
							misc.multi_free_base_penalty;
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 3;
							}
						}
						if (fres[i+1].pair <= -1 && fres[j-1].pair <= -1)
						{
							tmp = vm->get_energy_WM (i+2,k,tree) + vm->get_energy_WM (k+1, j-2,tree) +
							dangle_top [int_sequence[i]][int_sequence[j]][int_sequence[i+1]] +
							dangle_bot [int_sequence[i]][int_sequence[j]][int_sequence[j-1]] +
							2*misc.multi_free_base_penalty;
							if (tmp < min)
							{
								min = tmp;
								best_k = k;
								best_row = 4;
							}
						}

						// Hosna: June 28, 2007
						// the last branch of VM, which is WMB_(i+1),(j-1)
						// Hosna: July 5th, 2007
						// add a b_penalty to this case to match the previous cases
						tmp = WMB->get_WMB(i+1,j-1,tree)+ a_penalty + PSM_penalty+ b_penalty;
						if (tmp < min)
						{
							min = tmp;
							best_row = 5;
						}

					  }
					switch (best_row)
					  {
					  case 1:
		//              	printf("M_WM(%d,%d) branch 1: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-1);
						insert_node (i+1, best_k, M_WM);
						insert_node (best_k+1, j-1, M_WM);
						break;
					  case 2:
		//              	printf("M_WM(%d,%d) branch 2: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-1);
						insert_node (i+2, best_k, M_WM);
						insert_node (best_k+1, j-1, M_WM);
						break;
					  case 3:
		//              	printf("M_WM(%d,%d) branch 3: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+1,best_k,best_k+1,j-2);
						insert_node (i+1, best_k, M_WM);
						insert_node (best_k+1, j-2, M_WM);
						break;
					  case 4:
		//              	printf("M_WM(%d,%d) branch 4: pushing M_WM(%d,%d) and M_WM(%d,%d) \n", i,j,i+2,best_k,best_k+1,j-2);
						insert_node (i+2, best_k, M_WM);
						insert_node (best_k+1, j-2, M_WM);
						break;
					  // Hosna: June 28, 2007
					  // the last branch of VM, which is WMB_(i+1),(j-1)
					  case 5:
		//              	printf("M_WM(%d,%d) branch 3: pushing P_WMB(%d,%d)\n", i,j,i+1,j-1);
						insert_node(i+1,j-1, P_WMB);
						break;
					  }
				}
					break;
			}
		}
			break;
		case FREE:
		{
			int j = cur_interval->j;

			if (j==0) return;

			int min = INF, tmp, best_row, i, best_i, acc, energy_ij;
			if (debug)
				printf ("\t(0,%d) FREE\n", j);

			// this case is for j unpaired, so I have to check that.
			if (fres[j].pair <= -1)
			{
				tmp = W[j-1];
				if (tmp < min)
				{
					min = tmp;
					best_row = 0;
				}
			}
			for (i=0; i<=j-1; i++)    // no TURN
			{

				// Don't need to make sure i and j don't have to pair with something else
				//  it's INF, done in fold_sequence_restricted
				acc = (i-1>0) ? W[i-1] : 0;
				energy_ij = v->get_energy(i,j,tree);

				if (energy_ij < INF)
				{
					tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j]) + acc;
					if (tmp < min)
					{
					min = tmp;
					best_i = i;
					best_row = 1;
					}
					
				}

				if (fres[i].pair <= -1)
				{
					energy_ij = v->get_energy(i+1,j,tree);
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j]) + acc;
						tmp += dangle_bot [int_sequence[j]]
							[int_sequence[i+1]]
							[int_sequence[i]];
						
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 2;
							if(i==4 && j==242) printf("here");
						}
						
					}
				}
				if (fres[j].pair <= -1)
				{
					energy_ij = v->get_energy(i,j-1,tree);
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i],int_sequence[j-1]) + acc;
						tmp += dangle_top [int_sequence[j-1]]
							[int_sequence[i]]
							[int_sequence[j]];
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 3;
						}
					}
				}
				if (fres[i].pair <= -1 && fres[j].pair <= -1)
				{
					energy_ij = v->get_energy(i+1,j-1,tree);
					if (energy_ij < INF)
					{
						tmp = energy_ij + AU_penalty (int_sequence[i+1],int_sequence[j-1]) + acc;
						tmp += dangle_bot [int_sequence[j-1]]
							[int_sequence[i+1]]
							[int_sequence[i]];
						tmp += dangle_top [int_sequence[j-1]]
							[int_sequence[i+1]]
							[int_sequence[j]];
						if (tmp < min)
						{
							min = tmp;
							best_i = i;
							best_row = 4;
						}
					}
				}
			}

			// Hosna: June 28, 2007
			// the last branch of W, which is WMB_i,j
	//        energy_ij = WMB->get_energy(0,j);
	//        if (energy_ij < INF){
	//          	tmp = energy_ij + PS_penalty;
	//           	if (tmp < min){
	//           		min = tmp;
	//           		best_row = 5;
	//           	}
	//        }
		// Hosna June 30, 2007
		// The following would not take care of when
		// we have some unpaired bases before the start of the WMB
		for (i=0; i<=j-1; i++)
		{
			// Hosna: July 9, 2007
			// We only chop W to W + WMB when the bases before WMB are free
			if (i == 0 || (WMB->is_weakly_closed(0,i-1) && WMB->is_weakly_closed(i,j))){

				acc = (i-1>0) ? W[i-1]: 0;

				energy_ij = WMB->get_WMB(i,j,tree);

				if (energy_ij < INF)
				{
					tmp = energy_ij + PS_penalty + acc;

					if (tmp < min)
					{
						min = tmp;
						best_row = 5;
						best_i = i;
					}
				}

				// I have to condition on  fres[i].pair <= -1 to make sure that i can be unpaired
				if (fres[i].pair <= -1 && i+1 < j)
				{
					energy_ij = WMB->get_WMB(i+1,j,tree);
					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min)
						{
							min = tmp;
							best_row = 6;
							best_i = i;
						}
					}
				}

				// I have to condition on  fres[j].pair <= -1 to make sure that j can be unpaired
				if (fres[j].pair <= -1 && i < j-1)
				{
					energy_ij = WMB->get_WMB(i,j-1,tree);
					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min)
						{
							min = tmp;
							best_row = 7;
							best_i = i;
						}
					}
				}

				if (fres[i].pair <= -1 && fres[j].pair <= -1 && i+1 < j-1)
				{
					energy_ij = WMB->get_WMB(i+1,j-1,tree);
					if (energy_ij < INF)
					{
						tmp = energy_ij + PS_penalty + acc;
						if (tmp < min)
						{
							min = tmp;
							best_row = 8;
							best_i = i;
						}
					}
				}
			}
		}
			switch (best_row)
			{
				case 0:
					//printf("W(%d) case 0: inserting Free (0,%d)\n",j,j-1);
					insert_node (0, j-1, FREE); break;
				case 1:
					//printf("W(%d) case 1: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, LOOP);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 2:
					//printf("W(%d) case 2: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, LOOP);
					if (best_i >= 0)// Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 3:
					//printf("W(%d) case 3: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, LOOP);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 4:
					//printf("W(%d) case 4: inserting Loop(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, LOOP);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				// Hosna: June 28, 2007
				// the last branch of W, which is WMB_i,j
				case 5:
					//printf("W(%d) case 5: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j,best_i-1);
					insert_node (best_i, j, P_WMB);
					if (best_i-1 > 0)     // it was TURN instead of 0  - not sure if TURN shouldn't be here
						insert_node (0, best_i-1, FREE);
					break;
				case 6:
					//printf("W(%d) case 6: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j,best_i);
					insert_node (best_i+1, j, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
				case 7:
					//printf("W(%d) case 7: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i,j-1,best_i-1);
					insert_node (best_i, j-1, P_WMB);
					if (best_i-1 > 0)
						insert_node (0, best_i-1, FREE);
					break;
				case 8:
					//printf("W(%d) case 8: inserting WMB(%d,%d) and Free (0,%d)\n",j,best_i+1,j-1,best_i);
					insert_node (best_i+1, j-1, P_WMB);
					if (best_i >= 0) // Hosna, March 26, 2012, was best_i-1 instead of best_i
						insert_node (0, best_i, FREE);
					break;
			}
		}
			break;
		case M_WM:
//  else if(cur_interval->type == M_WM)
		{
			  int i = cur_interval->i;
			  int j = cur_interval->j;
			  int tmp, min = INF;
			  int best_k, best_row;

			  if (debug)
				printf ("\t (%d,%d) M_WM\n", i,j);

			  tmp = v->get_energy(i,j,tree) +
				AU_penalty (int_sequence[i], int_sequence[j]) +
				misc.multi_helix_penalty;
			  if (tmp < min)
				{
				  min = tmp;
				  best_row = 1;
				}
			  if (fres[i].pair <= -1)
			  {
				  tmp = v->get_energy(i+1,j,tree) +
						AU_penalty (int_sequence[i+1], int_sequence[j]) +
						dangle_bot [int_sequence[j]]
						[int_sequence[i+1]]
						[int_sequence[i]] +
						misc.multi_helix_penalty +
						misc.multi_free_base_penalty;
				  if (tmp < min)
				  {
					  min = tmp;
					  best_row = 2;
				  }
			  }
			  if (fres[j].pair <= -1)
			  {
				  tmp = v->get_energy(i,j-1,tree) +
						AU_penalty (int_sequence[i], int_sequence[j-1]) +
						dangle_top [int_sequence[j-1]]
									[int_sequence[i]]
									[int_sequence[j]] +
						misc.multi_helix_penalty +
						misc.multi_free_base_penalty;
				  if (tmp < min)
				  {
					  min = tmp;
					  best_row = 3;
				  }
			  }
			  if (fres[i].pair <= -1 && fres[j].pair <= -1)
			  {
				  tmp = v->get_energy(i+1,j-1,tree) +
						AU_penalty (int_sequence[i+1], int_sequence[j-1]) +
						dangle_bot [int_sequence[j-1]]
									[int_sequence[i+1]]
									[int_sequence[i]] +
						dangle_top [int_sequence[j-1]]
									[int_sequence[i+1]]
									[int_sequence[j]] +
						misc.multi_helix_penalty +
						2*misc.multi_free_base_penalty;
				  if (tmp < min)
				  {
					  min = tmp;
					  best_row = 4;
				  }
			  }
			  if (fres[i].pair <= -1)
			  {
				  tmp = vm->get_energy_WM (i+1,j,tree) + misc.multi_free_base_penalty;
				  if (tmp < min)
				  {
					  min = tmp;
					  best_row = 5;
				  }
			  }
			  if (fres[j].pair <= -1)
			  {
				  tmp = vm->get_energy_WM (i,j-1,tree) + misc.multi_free_base_penalty;
				  if (tmp < min)
				  {
					  min = tmp;
					  best_row = 6;
				  }
			  }

			  for (int k=i; k < j; k++)
				{
					tmp = vm->get_energy_WM (i, k,tree) + vm->get_energy_WM (k+1, j,tree);
					if (tmp < min)
					  {
						min = tmp;
						best_k = k;
						best_row = 7;
					  }
				}
			  // Hosna: June 28, 2007
			  // the last branch of WW, which is WMB_i,j
			  tmp = WMB->get_WMB(i,j,tree)+PSM_penalty;
			  if (tmp < min){
				min = tmp;
				best_row = 8;
			  }

			  switch (best_row)
				{
				  case 1: insert_node (i, j, LOOP); break;
				  case 2: insert_node (i+1, j, LOOP); break;
				  case 3: insert_node (i, j-1, LOOP); break;
				  case 4: insert_node (i+1, j-1, LOOP); break;
				  case 5:
					if (j-i-1 > 0)
					  insert_node (i+1, j, M_WM);
					break;
				  case 6:
					if (j-1-i > 0)
					  insert_node (i, j-1, M_WM);
					break;
				  case 7:
					if (best_k-i > 0)
					  insert_node (i, best_k, M_WM);
					if (j-best_k-1 > 0)
					  insert_node (best_k+1, j, M_WM);
					break;
				  // Hosna: June 28, 2007
				  // the last branch of W, which is WMB_i,j
				  case 8:
					insert_node(i,j,P_WMB);
					break;
				  }
			}
			break;
    // Hosna: Feb 19th 2007
		case P_WMB:
   // else if(cur_interval->type == P_WMB)
		{
			char struc[n+1];
			strcpy(struc,structure.c_str());
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(struc,f,cur_interval,tree);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WMBP:
    // Hosna: April 18th, 2007
    // changed WMB to case 2 and WMBP
   // else if(cur_interval->type == P_WMBP)
		{
			char struc[n+1];
			strcpy(struc,structure.c_str());
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(struc,f,cur_interval,tree);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VP:
    //else if(cur_interval->type == P_VP)
		{
			char struc[n+1];
			strcpy(struc,structure.c_str());
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(struc,f,cur_interval,tree);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_VPP:
    //else if(cur_interval->type == P_VPP)
		{
			char struc[n+1];
			strcpy(struc,structure.c_str());
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(struc,f,cur_interval,tree);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WI:
    //else if(cur_interval->type == P_WI)
		{
			char struc[n+1];
			strcpy(struc,structure.c_str());
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(struc,f,cur_interval,tree);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_BE:
    //else if(cur_interval->type == P_BE)
		{	
			char struc[n+1];
			strcpy(struc,structure.c_str());
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(struc,f,cur_interval,tree);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		case P_WIP:
    //else if(cur_interval->type == P_WIP)
		{
			char struc[n+1];
			strcpy(struc,structure.c_str());
			WMB->set_stack_interval(stack_interval);
			WMB->back_track(struc,f,cur_interval,tree);
			stack_interval = WMB->get_stack_interval();
			structure = WMB->get_structure();
			f = WMB->get_minimum_fold();
		}
			break;
		default:
			printf("Should not be here!\n");
	}


}

void W_final::insert_node (int i, int j, char type)
  // insert at the beginning
{
    seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;
}


//Mateo 13 Sept 2023
//return number of bases in between the two inclusive index
int distance(int left, int right){
    return (right-left-1);
}

//Mateo 13 Sept 2023
//given a initial hotspot which is a hairpin loop, keep trying to add a arc to form a larger stack
void expand_hotspot(s_energy_matrix *V, Hotspot &hotspot, int n){
    //printf("\nexpanding hotspot: i: %d j: %d\n",hotspot->get_left_inner_index(),hotspot->get_right_inner_index());
    double energy = 0;
    // int non_gc_penalty = 0;
    // int dangle_top_penalty = 0;
    // int dangle_bot_penalty = 0;
	int dangle_penalty = 0;

    //calculation for the hairpin that is already in there
    V->compute_hotspot_energy(hotspot.get_left_outer_index(),hotspot.get_right_outer_index(),0);


    //try to expand by adding a arc right beside the current out most arc
    while(hotspot.get_left_outer_index()-1 >= 0 && hotspot.get_right_outer_index()+1 <= n-1){
		int ptype_closing = pair[V->S_[hotspot.get_left_outer_index()-1+1]][V->S_[hotspot.get_right_outer_index()+1+1]];
        if(ptype_closing>0){
            hotspot.move_left_outer_index();
            hotspot.move_right_outer_index();
            hotspot.increment_size();
            V->compute_hotspot_energy(hotspot.get_left_outer_index(),hotspot.get_right_outer_index(),1);
            //printf("AU(i:%d,j:%d) = %d\n",hotspot->get_left_outer_index(),hotspot->get_right_outer_index(), AU_penalty (int_sequence[hotspot->get_left_outer_index()],int_sequence[hotspot->get_right_outer_index()]));
        }else{
            break;
        }
    }

    // non_gc_penalty += AU_penalty (int_sequence[hotspot.get_left_outer_index()],int_sequence[hotspot.get_right_outer_index()]);

    //if current out left-1 >= 0  (aka still have spot on left side of curent left out)
    //if current out right+1 <= nb_nuc-1 (aka still have spot on right side of curent right out)
    // if(hotspot.get_left_outer_index() - 1 >= 0 && hotspot.get_right_outer_index() + 1 <= n-1){
    //     int i = hotspot.get_left_outer_index()-1;
    //     int j = hotspot.get_right_outer_index()+1;
	// 	int tt = pair[V->S_[i+1]][V->S_[j+1]];
	// 	int si1 = i>0 ? V->S_[i-1] : -1;
	// 	int sj1 = j<n-1 ? V->S_[j+1] : -1;
	// 	dangle_penalty = vrna_E_ext_stem(tt, si1, sj1, V->params_);
    //     // dangle_bot_penalty = dangle_bot [int_sequence[j-1]][int_sequence[i+1]][int_sequence[i]];
    //     // dangle_top_penalty = dangle_top [int_sequence[j-1]][int_sequence[i+1]][int_sequence[j]];
    //     //printf("i: %d j: %d, dangle_bot: %d dangle_top: %d\n",i,j,dangle_bot_penalty,dangle_top_penalty);
    // }else if(hotspot.get_left_outer_index() - 1 >= 0){
    //     int i = hotspot.get_left_outer_index()-1;
    //     int j = nb_nucleotides-1;
    //     // dangle_bot_penalty = dangle_bot [int_sequence[j]][int_sequence[i+1]][int_sequence[i]];
    //     //printf("i: %d j: %d, dangle_bot: %d \n",i,j,dangle_bot_penalty);
    // }else if(hotspot.get_right_outer_index() + 1 <= nb_nucleotides-1){
    //     int i = 0;
    //     int j = hotspot.get_right_outer_index()+1;
    //     // dangle_top_penalty = dangle_top [int_sequence [j-1]][int_sequence [i]][int_sequence [j]];
    //     //printf("i: %d j: %d, angle_top: %d\n",i,j,dangle_top_penalty);
    // }
	 int i = hotspot.get_left_outer_index()-1;
	int j = hotspot.get_right_outer_index()+1;
	int tt = pair[V->S_[i+1]][V->S_[j+1]];
	int si1 = i>0 ? V->S_[i-1] : -1;
	int sj1 = j<n-1 ? V->S_[j+1] : -1;
	dangle_penalty = vrna_E_ext_stem(tt, si1, sj1, V->params_);


    energy = V->get_energy(hotspot.get_left_outer_index(),hotspot.get_right_outer_index());

    // printf("here and %d\n",energy);
    //printf("energy: %lf, AU_total: %d, dangle_top_total: %d, dangle_bot_total: %d\n",energy,non_gc_penalty,dangle_top_penalty,dangle_bot_penalty);

    energy = (energy + dangle_penalty) / 100;

    hotspot.set_energy(energy);
    //printf("done: %d %d %d %d\n",hotspot->get_left_outer_index(),hotspot->get_left_inner_index(),hotspot->get_right_inner_index(),hotspot->get_right_outer_index());
    return;
}

//Mateo 13 Sept 2023
//look for every possible hairpin loop, and try to add a arc to form a larger stack with at least min_stack_size bases
void get_hotspots(std::string seq,std::vector<Hotspot> &hotspot_list,int max_hotspot, vrna_param_s *params){
    
	int n = seq.length();
	s_energy_matrix *V;
	V = new s_energy_matrix (seq,n,params);
	make_pair_matrix();
    int min_bp_distance = 3;
    int min_stack_size = 3; //the hotspot must be a stack of size >= 3
    // Hotspot current_hotspot;
    //start at min_stack_size-1 and go outward to try to add more arcs to form bigger stack because we cannot expand more than min_stack_size from there anyway
    for(int i = min_stack_size-1; i < n; i++){
        for(int j = i; j < n; j++){
			int ptype_closing = pair[V->S_[i+1]][V->S_[j+1]];
            if(ptype_closing>0 && distance(i,j) >= min_bp_distance){
                // current_hotspot = new Hotspot(i,j,nb_nucleotides);
                Hotspot current_hotspot(i,j,n);

                expand_hotspot(V,current_hotspot,n);


                if(current_hotspot.get_size() < min_stack_size || current_hotspot.is_invalid_energy()){

                }else{
                    
                    current_hotspot.set_structure();
                    hotspot_list.push_back(current_hotspot);

                }
            }
        }
    }

    //make sure we only keep top 20 hotspot with lowest energy
    std::sort(hotspot_list.begin(), hotspot_list.end(),compare_hotspot_ptr);
    while(hotspot_list.size() > max_hotspot){
        // delete hotspot_list.back();
        hotspot_list.pop_back();
    }

    //if no hotspot found, add all _ as restricted
    if(hotspot_list.size() == 0){
        // Hotspot* hotspot = new Hotspot(0,nb_nucleotides-1,nb_nucleotides);
        Hotspot hotspot(0,n-1,n);
        hotspot.set_default_structure();
        hotspot_list.push_back(hotspot);
    }
	delete V;

    return;
}

bool compare_hotspot_ptr(Hotspot &a, Hotspot &b) { 
    return (a.get_energy() < b.get_energy()); 
}