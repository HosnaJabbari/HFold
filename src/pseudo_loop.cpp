#include "pseudo_loop.h"
#include "constants.h"
#include "h_struct.h"
#include "h_externs.h"
#include "h_common.h"

#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>

#include "VM_final.h"

// set to -1 because >= 0 means there is already a base pair there,
// and -1 means restricted struture says there is no base pair there.
#define FRES_RESTRICTED_MIN -1

pseudo_loop::pseudo_loop(std::string seq, char* restricted, s_energy_matrix *V, VM_final *VM, short *S, short *S1, vrna_param_t *params)
{
	this->seq = seq;
	this->restricted = restricted;
	this->V = V;
	this->VM = VM;
	S_ = S;
	S1_ = S1;
	params_ = params;
	make_pair_matrix();
    allocate_space();
}

void pseudo_loop::allocate_space()
{
    n = seq.length();

    index.resize(n);
    cand_pos_t total_length = (n *(n+1))/2;
    index[0] = 0;
    for (cand_pos_t i=1; i < n; i++)
        index[i] = index[i-1]+n-i+1;

    WI.resize(total_length,0);

    VP.resize(total_length,INF);

    WMB.resize(total_length,INF);

    WMBP.resize(total_length,INF);

    WIP.resize(total_length,INF);

    VPP.resize(total_length,INF);

    BE.resize(total_length,0);

}

pseudo_loop::~pseudo_loop()
{
}

void pseudo_loop::compute_energies(cand_pos_t i, cand_pos_t j, sparse_tree &tree)
{

	compute_VP(i,j,tree); // Hosna, March 14, 2012, changed the positionof computing VP from after BE to befor WMBP

	compute_WMBP(i,j,tree);

    compute_WMB(i,j,tree);

    compute_WI(i,j,tree);

    compute_WIP(i,j,tree);

    compute_VPP(i,j,tree);

	compute_BE(tree.tree[j+1].pair-1,j,tree.tree[i+1].pair-1,i,tree);

}
// Added +1 to fres/tree indices as they are 1 ahead at the moment
void pseudo_loop::compute_WI(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	energy_t min = INF, m1 = INF, m2= INF, m3= INF;
	cand_pos_t ij = index[i]+j-i;
	if (WI[ij] != 0){ //calculated before
		return;
	}

	//base cases
	// if [i,j] is not weakly closed then WI[i,j] = INF
	if (!tree.weakly_closed(i+1,j+1)){
		WI[ij] = INF;
		return;
	}

	// branch 4, one base
	if (i == j){
		WI[ij] = PUP_penalty;
		return;
	}

	// Hosna: Feb 16, 2007:
	// we don't need to check to see if i and j are inside an arc
	// because they are not in an arc in G but they will be in an arc in G'
	if (tree.tree[i+1].parent->index != tree.tree[j+1].parent->index){
		WI[ij] = INF;
		return;
	}

// Hosna: July 2nd, 2007
// in branch 1 of WI, we can have a case like
// ((..))((...))
// such that both i and j are paired but we can chop them

	for (cand_pos_t t = i; t< j; t++){
		int wi_1 = get_WI(i,t,tree);
		int wi_2 = get_WI(t+1,j,tree);
		int energy = wi_1 + wi_2;
		m1 = (m1 > energy)? energy : m1;
	}
	// branch 2:

	if ((tree.tree[i+1].pair == j+1 && tree.tree[j+1].pair == i+1) ||(tree.tree[i+1].pair < FRES_RESTRICTED_MIN && tree.tree[j+1].pair < FRES_RESTRICTED_MIN)){
		// Hosna, April 16th, 2007
		// changed back to see if it will work fine
		// Hosna: April 19th, 2007
		// I think we should call the restricted version

		energy_t v_ener = (i>j)? INF: V->get_energy(i,j);
		m2 = v_ener + PPS_penalty;
	}

	m3 = get_WMB(i,j,tree) + PSP_penalty + PPS_penalty;

	min = std::min(m1,std::min(m2,m3));
	WI[ij] = min;
}


void pseudo_loop::compute_VP(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	cand_pos_t ij = index[i]+j-i;
	if (VP[ij] != INF){//has been calculated before
		return;
	}
	const pair_type ptype_closing = pair[S_[i+1]][S_[j+1]];
	
	// base cases:
	// a) i == j => VP[ij] = INF
	// b) [i,j] is a weakly_closed region => VP[ij] = INF
	// c) i or j is paired in original structure => VP[ij] = INF

	if (i == j || j-i<4 || tree.weakly_closed(i+1,j-1) || tree.tree[i+1].pair >= FRES_RESTRICTED_MIN || tree.tree[j+1].pair >= FRES_RESTRICTED_MIN || ptype_closing == 0)	{
		VP[ij] = INF;
		return;
	}
	else{
		energy_t m1 = INF, m2 = INF, m3 = INF, m4= INF, m5 = INF, m6 = INF, m7 = INF; //different branches
		
		// Borders -- added one to i and j to make it fit current bounds but also subtracted 1 from answer as the tree bounds are shifted as well
		cand_pos_t Bp_ij = tree.Bp(i+1,j+1)-1;
		cand_pos_t B_ij = tree.B(i+1,j+1)-1;
		cand_pos_t b_ij = tree.b(i+1,j+1)-1;
		cand_pos_t bp_ij = tree.bp(i+1,j+1)-1;
		pair_type ptype_closing = pair[S_[i+1]][S_[j+1]];
		
		//branchs:
		// 1) inArc(i) and NOT_inArc(j)
		// WI(i+1)(B'(i,j)-1)+WI(B(i,j)+1)(j-1)

		// Hosna April 9th, 2007
		// need to check the borders as they may be negative
		if((tree.tree[i+1].parent->index-1) > 0 && (tree.tree[j+1].parent->index-1) < (tree.tree[i+1].parent->index-1) && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0){
			energy_t WI_ipus1_BPminus = get_WI(i+1,Bp_ij - 1,tree) ;
			energy_t WI_Bplus_jminus = get_WI(B_ij + 1,j-1,tree);
			m1 =   WI_ipus1_BPminus + WI_Bplus_jminus;
		}

		// 2) NOT_inArc(i) and inArc(j)
		// WI(i+1)(b(i,j)-1)+WI(b'(i,j)+1)(j-1)

		// Hosna April 9th, 2007
		// checking the borders as they may be negative
		if ((tree.tree[i+1].parent->index-1) < (tree.tree[j+1].parent->index-1) && (tree.tree[j+1].parent->index-1) > 0 && b_ij>= 0 && bp_ij >= 0 && Bp_ij < 0){
			energy_t WI_i_plus_b_minus = get_WI(i+1,b_ij - 1,tree);
			energy_t WI_bp_plus_j_minus = get_WI(bp_ij + 1,j-1,tree);
			m2 = WI_i_plus_b_minus + WI_bp_plus_j_minus;
		}

		// 3) inArc(i) and inArc(j)
		// WI(i+1)(B'(i,j)-1)+WI(B(i,j)+1)(b(i,j)-1)+WI(b'(i,j)+1)(j-1)

		// Hosna April 9th, 2007
		// checking the borders as they may be negative
		if((tree.tree[i+1].parent->index-1) > 0 && (tree.tree[j+1].parent->index-1) > 0 && Bp_ij >= 0 && B_ij >= 0  && b_ij >= 0 && bp_ij>= 0){
			energy_t WI_i_plus_Bp_minus = get_WI(i+1,Bp_ij - 1,tree);
			energy_t WI_B_plus_b_minus = get_WI(B_ij + 1,b_ij - 1,tree);
			energy_t WI_bp_plus_j_minus = get_WI(bp_ij +1,j - 1,tree);
			m3 = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
		}

		// 4) NOT_paired(i+1) and NOT_paired(j-1) and they can pair together
		// e_stP(i,i+1,j-1,j) + VP(i+1)(j-1)
		pair_type ptype_closingip1jm1 = pair[S_[i+1+1]][S_[j-1+1]];
		if((tree.tree[i+1+1].pair) < FRES_RESTRICTED_MIN && (tree.tree[j-1+1].pair) < FRES_RESTRICTED_MIN && ptype_closingip1jm1>0){
			m4 = get_e_stP(i,j)+ get_VP(i+1,j-1,tree);
		}

		// 5) NOT_paired(r) and NOT_paired(rp)
		//  VP(i,j) = e_intP(i,ip,jp,j) + VP(ip,jp)
		// Hosna, April 6th, 2007
		// whenever we use get_borders we have to check for the correct values
		cand_pos_t min_borders = std::min((cand_pos_tu) Bp_ij, (cand_pos_tu) b_ij);
		cand_pos_t edge_i = std::min(i+MAXLOOP+1,j-TURN-1);
		min_borders = std::min({min_borders,edge_i});
//		printf("B'(%d,%d) = %d, b(%d,%d) = %d, min_borders = %d\n",i,j,get_Bp(i,j),i,j,get_b(i,j), min_borders);
		for (cand_pos_t k = i+1; k < min_borders; ++k){
			// Hosna: April 20, 2007
			// i and ip and j and jp should be in the same arc
			// also it should be the case that [i+1,ip-1] && [jp+1,j-1] are empty regions

			if (tree.tree[k+1].pair < FRES_RESTRICTED_MIN && (tree.tree[i+1].parent->index == tree.tree[k+1].parent->index) && (tree.up[(k+1)-1] > ((k+1)-(i+1)-1))){
				// Hosna, April 6th, 2007
				// whenever we use get_borders we have to check for the correct values
				cand_pos_t max_borders = std::max(bp_ij,B_ij)+1;
				cand_pos_t edge_j = k+j-i-MAXLOOP-2;
				max_borders = std::max({max_borders,edge_j});
//				printf("b'(%d,%d) = %d, B(%d,%d) = %d, max_borders = %d\n",i,j,get_bp(i,j),i,j,get_B(i,j), max_borders);
				for (cand_pos_t l = j-1; l > max_borders ; --l){

					pair_type ptype_closingkj = pair[S_[k+1]][S_[l+1]];
					if (tree.tree[l+1].pair < FRES_RESTRICTED_MIN && ptype_closingkj>0 && (tree.up[(j+1)-1] > ((j+1)-(l+1)-1))){
						// Hosna: April 20, 2007
						// i and ip and j and jp should be in the same arc
						if (tree.tree[j+1].parent->index == tree.tree[l+1].parent->index){
							energy_t tmp = get_e_intP(i,k,l,j) + get_VP(k,l,tree);

							m5 = std::min(m5,tmp);
						}
					}
				}
			}
		}

		// 6) VP(i,j) = WIP(i+1,r-1) + VPP(r,j-1)
		// Hosna April 9th, 2007
		// checking the borders as they may be negative numbers
		cand_pos_t min_Bp_j = j;
		if (Bp_ij > 0 && Bp_ij < n && Bp_ij < min_Bp_j){
			min_Bp_j = Bp_ij;
		}
		for (cand_pos_t r = i+1; r < min_Bp_j ; ++r){


			if (tree.tree[r+1].pair < FRES_RESTRICTED_MIN){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				int tmp = get_WIP(i+1,r-1,tree) + get_VPP(r,j-1,tree) + ap_penalty + 2*bp_penalty;
				m6 = std::min(m6,tmp);
			}
		}


		// 7) VP(i,j) = VPP(i+1,r) + WIP(r+1,j-1)
		// Hosna April 9th, 2007
		// checking the borders as they may be negative numbers
		cand_pos_t max_i_bp = i;
		if (bp_ij > 0 && bp_ij < n && bp_ij > max_i_bp){
			max_i_bp = bp_ij;
		}
		for (cand_pos_t r = max_i_bp+1; r < j ; ++r){

			if (tree.tree[r+1].pair < FRES_RESTRICTED_MIN){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				energy_t tmp = get_VPP(i+1,r,tree) + get_WIP(r+1,j-1,tree)+ ap_penalty + 2* bp_penalty;
				m7 = std::min(m7,tmp);
			}
		}

		//finding the min energy
		energy_t vp_h = std::min({m1,m2,m3});
		energy_t vp_iloop = std::min({m4,m5});
		energy_t vp_split = std::min({m6,m7});
		energy_t min = std::min({vp_h,vp_iloop,vp_split});
		

		VP[ij] = min;

	}
}

void pseudo_loop::compute_WMBP(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	cand_pos_t ij = index[i]+j-i;
	if (WMBP[ij] != INF){
		return;
	}
	//base case
	if (i == j){
		WMBP[ij] = INF;
		return;
	}
	// Hosna: July 6th, 2007
	// added impossible cases

	if ((tree.tree[i+1].pair >= FRES_RESTRICTED_MIN && tree.tree[i+1].pair > j+1)
	||  (tree.tree[j+1].pair >= FRES_RESTRICTED_MIN && tree.tree[j+1].pair < i+1)
	||  (tree.tree[i+1].pair >= FRES_RESTRICTED_MIN && tree.tree[i+1].pair < i+1 )
	||  (tree.tree[j+1].pair >= FRES_RESTRICTED_MIN && j+1 < tree.tree[j+1].pair)){
		WMB[ij] = INF;
		return;
	}
	else{
		energy_t m1 = INF, m3 = INF, m4 = INF, m5 = INF;
		// if not paired(j) and paired(i) then
		// WMBP(i,j) = 2*Pb + min_{i<l<bp(i)}(BE(i,bp(i),b'(i,l),bp(b'(i,l)))+WI(b'+1,l-1)+VP(l,j))
		if(tree.tree[j+1].pair < 0 && tree.tree[i+1].pair >= 0){
			energy_t tmp = INF;
			// Hosna: June 29, 2007
			// if j is inside i's arc then the l should be
			// less than j not bp(i)
			// check with Anne
			for (cand_pos_t l = i+1; l < j; l++){
				// Hosna, March 14, 2007
				// fixing the for loop

				// Hosna, April 9th, 2007
				// checking the borders as they may be negative
				cand_pos_t bp_il = tree.bp(i+1,l+1)-1;
				if(bp_il >= 0 && bp_il < n && l+TURN <= j){
					energy_t BE_energy = get_BE(i,tree.tree[i+1].pair-1,bp_il,tree.tree[bp_il+1].pair-1,tree);
					energy_t WI_energy = get_WI(bp_il +1,l-1,tree);
					energy_t VP_energy = get_VP(l,j,tree);
					energy_t sum = BE_energy + WI_energy + VP_energy;
					tmp = std::min(tmp,sum);
				}
			}
			m1 = 2*PB_penalty + tmp;

		}
			

		// 3)
		cand_pos_t b_ij = tree.b(i+1,j+1)-1;
		if (tree.tree[j+1].pair < 0){
			energy_t tmp = INF;
			for (cand_pos_t l = i+1; l<j ; l++)	{
				// Hosna, April 6th, 2007
				// whenever we use get_borders we have to check for the correct values
				cand_pos_t B_lj = tree.B(l+1,j+1)-1;
				cand_pos_t Bp_lj = tree.Bp(l+1,j+1)-1;
				if (tree.tree[l+1].parent->index > -1 && B_lj >= 0 && B_lj < n && Bp_lj >= 0 && Bp_lj<n){
					// Hosna: April 19th, 2007
					// the chosen l should be less than border_b(i,j)
					if (b_ij >= 0 && b_ij < n && l < b_ij){

						// Hosna: June 29 2007
						// after going over the program with Cristina, we noticed that
						// l should be < B'(i,j)
						//if (l < get_Bp(i,j) && l+TURN <= j){

						// Hosna: July 5th, 2007:
						// as long as we have i <= arc(l)< j we are fine
						if (i <= tree.tree[l+1].parent->index-1 && tree.tree[l+1].parent->index-1 < j && l+TURN <=j){
							energy_t sum = get_BE(tree.tree[B_lj+1].pair-1,B_lj,tree.tree[Bp_lj+1].pair-1,Bp_lj,tree)+ get_WMBP(i,l-1,tree)+ get_VP(l,j,tree);
							tmp = std::min(tmp,sum);
						}
					}
				}
				// Hosna: April 5th
				// after going over the WMB recurrence with Anne, we think we should add another p_b penalty
				// to the 3rd case ==> 2*P_b
				m3 = 2*PB_penalty + tmp;
			}
		}

		// 4) WMB(i,j) = VP(i,j) + P_b
		energy_t tmp = get_VP(i,j,tree) + PB_penalty;
		if (tmp < m4){
			m4 = tmp;
		}
		// 5) WMB(i,j) = min_{i<l<j}{WMB(i,l)+WI(l+1,j)} if bp(j)<j
		// Hosna: Feb 5, 2007
		if(tree.tree[j+1].pair < j+1){
			for(cand_pos_t l = i+1; l<j; l++){
				// Hosna: March 14th, 2007
				// I think l cannot be paired

				// Hosna: April 18th, 2007
				// l and j should be in the same arc
				if (tree.tree[l+1].pair < 0 && tree.tree[l+1].parent->index > -1 && tree.tree[j+1].parent->index > -1 && tree.tree[l+1].parent->index == tree.tree[l+1].parent->index){
					tmp = get_WMBP(i,l,tree) + get_WI(l+1,j,tree);
					m5 = std::min(m5,tmp);

				}
			}
		}

		// get the min for WMB
		WMBP[ij] = std::min(std::min(m1,m3),std::min(m4,m5));
	}


}

void pseudo_loop::compute_WMB(cand_pos_t  i, cand_pos_t  j, sparse_tree &tree){
	cand_pos_t ij = index[i]+j-i;
	if (WMB[ij] != INF){
		return;
	}
	//base case
	if (i == j){
		WMB[ij] = INF;
		return;
	}
	// Hosna: July 6th, 2007
	// added impossible cases

	if ((tree.tree[i+1].pair >= FRES_RESTRICTED_MIN && tree.tree[i+1].pair > j+1)
	||  (tree.tree[j+1].pair >= FRES_RESTRICTED_MIN && tree.tree[j+1].pair < i+1)
	||  (tree.tree[i+1].pair >= FRES_RESTRICTED_MIN && tree.tree[i+1].pair < i+1 )
	||  (tree.tree[j+1].pair >= FRES_RESTRICTED_MIN && j+1 < tree.tree[j+1].pair)){
		WMB[ij] = INF;
		return;
	}
	else{
		energy_t m2 = INF, mWMBP = INF;
		// 2)
		if (tree.tree[j+1].pair >= 0 && j+1 > tree.tree[j+1].pair){
			cand_pos_t bp_j = tree.tree[j+1].pair-1;
			energy_t temp = INF;
			for (cand_pos_t l = (bp_j +1); (l < j); l++){
				// Hosna: April 24, 2007
				// correct case 2 such that a multi-pseudoknotted
				// loop would not be treated as case 2

				cand_pos_t Bp_lj = tree.Bp(l+1,j+1)-1;
				if (Bp_lj >= 0 && Bp_lj<n){
					energy_t sum = get_BE(bp_j,j,tree.tree[Bp_lj+1].pair-1,Bp_lj,tree) + get_WMBP(i,l,tree) + get_WI(l+1,Bp_lj-1,tree);
					temp = std::min(temp,sum);
				}

			}
			m2 = PB_penalty + temp;
		}
		// check the WMBP value
		mWMBP =  get_WMBP(i,j,tree);

		// get the min for WMB
		WMB[ij] = MIN(m2,mWMBP);
	}
}

void pseudo_loop::compute_WIP(cand_pos_t  i, cand_pos_t  j, sparse_tree &tree){
	cand_pos_t ij = index[i]+j-i;
	if (WIP[ij] < INF/2){ // was calculated before
		return;
	}
	if (tree.tree[i+1].parent->index != tree.tree[j+1].parent->index || i == j || !tree.weakly_closed(i+1,j+1)){
		WIP[ij] = INF;
		return;
	}
	energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;

    // Ian Wark July 19 2017
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs

	// branch 1:
	if (tree.tree[i+1].pair < FRES_RESTRICTED_MIN){
		m1 = get_WIP(i+1,j,tree) + cp_penalty;
	}
	// branch 2:
	if (tree.tree[j+1].pair < FRES_RESTRICTED_MIN){
		m2 = get_WIP(i,j-1,tree) + cp_penalty;
	}
	//branch 3:
	for (cand_pos_t t = i; t <j; t++){
		int tmp = get_WIP(i,t,tree) + get_WIP(t+1,j,tree);
		m3 = std::min(m3,tmp);
	}

	pair_type ptype_closing = pair[S_[i+1]][S_[j+1]];
	// branch 4:
	if (tree.tree[i+1].pair == j+1 || (tree.tree[i+1].pair < FRES_RESTRICTED_MIN && tree.tree[j+1].pair < FRES_RESTRICTED_MIN && ptype_closing >0)){

		m4 = V->get_energy(i,j) + bp_penalty;

	}

	// branch 5:
	m5 = get_WMB(i,j,tree) + PSM_penalty + bp_penalty;

	WIP[ij] = MIN(MIN(m1,MIN(m2,m3)),MIN(m4,m5));

}

void pseudo_loop::compute_VPP(cand_pos_t i, cand_pos_t  j, sparse_tree &tree){
	cand_pos_t ij = index[i]+j-i;
	if (VPP[ij] != INF){ // computed before
		return;
	}
	if (i == j  || tree.weakly_closed(i+1,j+1)){
		VPP[ij] = INF;
		return;
	}
	energy_t m1 = INF, m2 = INF,m3 = INF, m4 = INF;

	//branch 1:
	// Hosna April 9th, 2007
	// checking the borders as they may be negative numbers
	cand_pos_t max_i_bp = i;
	cand_pos_t bp_ij = tree.bp(i+1,j+1)-1;
	if (bp_ij > 0 && bp_ij < n && bp_ij > max_i_bp){
		max_i_bp = bp_ij;
	}
	for (cand_pos_t r = max_i_bp+1; r < j; r++ ){
        // Ian Wark July 19 2017
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (tree.tree[r+1].pair < FRES_RESTRICTED_MIN){
			energy_t tmp = get_VP(i,r,tree) + get_WIP(r+1,j,tree);
			m1 = std::min(m1,tmp);
		}
	}

	//branch 2:
	// Hosna April 9th, 2007
	// checking the borders as they may be negative numbers
	cand_pos_t min_Bp_j = j;
	cand_pos_t Bp_ij = tree.Bp(i+1,j+1)-1;
	if (Bp_ij > 0 && Bp_ij < n && bp_ij < min_Bp_j){
		min_Bp_j = Bp_ij;
	}
	for (cand_pos_t r = i+1; r < min_Bp_j; r++){
        // Ian Wark July 19 2017
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (tree.tree[r+1].pair < FRES_RESTRICTED_MIN){
			energy_t tmp = get_WIP(i,r-1,tree) + get_VP(r,j,tree);
			m2 = std::min(m2,tmp);
		}
	}

	// Branch 3:
	for (cand_pos_t r = max_i_bp+1; r < j; r++ ){
        // Ian Wark July 19 2017
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (tree.tree[r+1].pair < FRES_RESTRICTED_MIN && tree.up[j+1-1] >= (j-(r+1)-1)){
			energy_t tmp = get_VP(i,r,tree) + (cp_penalty *(j-r)); // check the (j-r) part
			m3 = std::min(m3,tmp);
		}
	}

	// Branch 4:

	for (cand_pos_t r = i+1; r < min_Bp_j; r++){
        // Ian Wark July 19 2017
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (tree.tree[r+1].pair < FRES_RESTRICTED_MIN && tree.up[r-1+1-1] >= ((r-1)-(i)-1)){
			energy_t tmp = (cp_penalty * (r-i)) + get_VP(r,j,tree);
			m4 = std::min(m4,tmp);
		}
	}

	VPP[ij] = std::min({m1,m2,m3,m4}); //MIN(MIN(m1,m2),MIN(m3,m4));
}

void pseudo_loop::compute_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, sparse_tree &tree){


    // Ian Wark July 19 2017
    // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if (!( i >= 0 && i <= ip && ip < jp && jp <= j && j < n && tree.tree[i+1].pair >= FRES_RESTRICTED_MIN && tree.tree[j+1].pair >= FRES_RESTRICTED_MIN && tree.tree[ip+1].pair >= FRES_RESTRICTED_MIN && tree.tree[jp+1].pair >= FRES_RESTRICTED_MIN && tree.tree[i+1].pair == j+1 && tree.tree[j+1].pair == i+1 && tree.tree[ip+1].pair == jp+1 && tree.tree[jp+1].pair == ip+1)){ //impossible cases
		return;
	}
	cand_pos_t iip = index[i]+ip-i;
	if (BE[iip] != 0){ // computed before
		return;
	}
	// base case: i.j and ip.jp must be in G
	if (tree.tree[i+1].pair != j+1 || tree.tree[ip+1].pair != jp+1){
		BE[iip] = INF;
		return;
	}

	// base case:
	if(i == ip && j == jp && i<j){
		BE[iip] = 0;
		return;
	}

	energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;
	// 1) bp(i+1) == j-1
	if (tree.tree[i+1+1].pair == j-1+1){
		m1 = get_e_stP(i,j) + get_BE(i+1,j-1,ip,jp,tree);

	}

	// cases 2-5 are all need an l s.t. i<l<=ip and jp<=bp(l)<j
	for (cand_pos_t l = i+1; l<= ip ; l++){

		// Hosna: March 14th, 2007
		if (tree.tree[l+1].pair >= FRES_RESTRICTED_MIN && jp <= tree.tree[l+1].pair-1 && tree.tree[l+1].pair-1 < j){
			// Hosna, March 15, 2007
			// since not_paired_all[i,l] includes i and l themselves
			// and in BE energy calculation we are looking for the oepn region (i,l)
			// we have to look at not_paired_all[i+1,l-1]
			cand_pos_t lp = tree.tree[l].pair-1;
			cand_pos_t il = index[i]+l-i;
			cand_pos_t lpj = index[lp]+j-lp;
			// 2)
			// Hosna June 29, 2007
			// when we pass a stacked pair instead of an internal loop to e_int, it returns underflow,
			// so I am checking explicitely that we won't have stems instead of internal loop
			bool empty_region_il = (tree.up[(l+1)-1] >= l-i-1); //empty between i+1 and lp-1
			bool empty_region_lpj = (tree.up[(j+1)-1] >= j-lp-1); // empty between l+1 and ip-1
			bool weakly_closed_il = tree.weakly_closed(i+1+1,l-1+1); // weakly closed between i+1 and lp-1
			bool weakly_closed_lpj = tree.weakly_closed(lp+1+1,j-1+1); // weakly closed between l+1 and ip-1


			if (empty_region_il && empty_region_lpj){//&& !(ip == (i+1) && jp==(j-1)) && !(l == (i+1) && lp == (j-1))){
				energy_t tmp = get_e_intP(i,l,lp,j)+ get_BE(l,lp,ip,jp,tree);
				m2 = std::min(m2,tmp);
			}

			// 3)
			if (weakly_closed_il && weakly_closed_lpj){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				energy_t tmp = get_WIP(i+1,l-1,tree) + get_BE(l,lp,ip,jp,tree) + get_WIP(lp+1,j-1,tree)+ ap_penalty + 2* bp_penalty;
				m3 = std::min(m3,tmp);
			}

			// 4)
			if (weakly_closed_il && empty_region_lpj){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				energy_t tmp = get_WIP(i+1,l-1,tree) + get_BE(l,lp,ip,jp,tree) + cp_penalty * (j-lp+1) + ap_penalty + 2*bp_penalty;
				m4 = std::min(m4,tmp);
			}

			// 5)
			if (empty_region_il && weakly_closed_lpj){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				energy_t tmp = ap_penalty + 2*bp_penalty + (cp_penalty * (l-i+1)) + get_BE(l,lp,ip,jp,tree) + get_WIP(lp+1,j-1,tree);
				m5 = std::min(m5,tmp);
			}
		}
	}

	// finding the min and putting it in BE[iip]
	BE[iip] = std::min({m1,m2,m3,m4,m5});
}

energy_t pseudo_loop::get_WI(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	if (i>j){
		return 0;
	}
	cand_pos_t ij = index[i]+j-i;
	//printf("get_WI(%d,%d), after computation its value = %d!\n",i,j, WI[ij]);
	return WI[ij];


}

energy_t pseudo_loop::get_VP(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	// Hosna, March 16, 2012
	// two bases should be at least 3 bases apart

	if (j-i < TURN || i >= j || tree.tree[i+1].pair >= 0 || tree.tree[j+1].pair >= 0 || tree.weakly_closed(i+1,j+1)){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;

	return VP[ij];

}
energy_t pseudo_loop::get_WMB(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	// Hosna: July 6th, 2007
	// added impossible cases
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN ||(tree.tree[i+1].pair >= 0 && tree.tree[i+1].pair > j+1) || (tree.tree[j+1].pair >= 0 && tree.tree[j+1].pair < i+1) || (tree.tree[i+1].pair >= 0 && tree.tree[i+1].pair < i+1 ) || (tree.tree[j+1].pair >= 0 && j+1 < tree.tree[j+1].pair)){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;

	return WMB[ij];
}

// Hosna: April 18th, 2007
// changed WMB to case 2 and WMBP
energy_t pseudo_loop::get_WMBP(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	// Hosna: July 6th, 2007
	// added impossible cases
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i< TURN || (tree.tree[i+1].pair >= 0 && tree.tree[i+1].pair > j+1) || (tree.tree[j+1].pair >= 0 && tree.tree[j+1].pair < i+1) || (tree.tree[i+1].pair >= 0 && tree.tree[i+1].pair < i+1 ) || (tree.tree[j+1].pair >= 0 && j+1 < tree.tree[j+1].pair)){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;
	return WMBP[ij];
}

energy_t pseudo_loop::get_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, sparse_tree &tree){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i>= TURN && i >= 0 && i <= ip && ip < jp && jp <= j && j < n && tree.tree[i+1].pair >=0 && tree.tree[j+1].pair >= 0 && tree.tree[ip+1].pair >= 0 && tree.tree[jp+1].pair >= 0 && tree.tree[i+1].pair == j+1 && tree.tree[j+1].pair == i+1 && tree.tree[ip+1].pair == jp+1 && tree.tree[jp+1].pair == ip+1){
		if(i == ip && j == jp && i<j){
			return 0;
		}
		cand_pos_t iip = index[i]+ip-i;

		return BE[iip];
	}else{
		return INF;
	}
}

energy_t pseudo_loop::get_WIP(cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN || i >= j || !tree.weakly_closed(i+1,j+1)){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;

	return WIP[ij];
}

energy_t pseudo_loop::get_VPP(cand_pos_t i, cand_pos_t j,sparse_tree &tree){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN || i >= j || tree.weakly_closed(i+1,j+1)){
		return INF;
	}
	cand_pos_t ij = index[i]+j-i;

	return VPP[ij];

}

energy_t pseudo_loop::compute_int(cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l, const paramT *params){

	const pair_type ptype_closing = pair[S_[i]][S_[j]];
    return E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<paramT *>(params));
}

energy_t pseudo_loop::get_e_stP(cand_pos_t i, cand_pos_t j){
	if (i+1 == j-1){ // TODO: do I need something like that or stack is taking care of this?
		return INF;
	}
	energy_t ss = compute_int(i+1,j+1,i+1+1,j-1+1,params_);
	return lrint(e_stP_penalty * ss);
}

energy_t pseudo_loop::get_e_intP(cand_pos_t i, cand_pos_t ip, cand_pos_t jp, cand_pos_t j){
	// Hosna Feb 12th, 2007:
	// this function is only being called in branch 5 of VP
	// and branch 2 of BE
	// in both cases regions [i,ip] and [jp,j] are closed regions
	energy_t e_int = compute_int(i+1,j+1,ip+1,jp+1,params_);

	// Hosna April 3rd, 2007
	// based on the discussion with Anne, we decided to have
	// e_intP = 0.83 * e_int
	energy_t energy = lrint(e_intP_penalty * e_int);
	return energy;
}

void pseudo_loop::back_track(std::string structure, minimum_fold *f, seq_interval *cur_interval, sparse_tree &tree)
{
	this->structure = structure;
	this->f = f;
	// Hosna March 8, 2012
	// changing the nested if structure to switch for optimality
	switch (cur_interval->type)
	{
			case P_WMB:
			{
				int i = cur_interval->i;
				int j = cur_interval->j;
				if (i >= j)
					return;
				int tmp = INF, best_l = -1, best_row = -1, min = INF;

				// case 1
				if (tree.tree[j+1].pair >= 0 && j+1 > tree.tree[j+1].pair){
					int l, acc = INF;
					int bp_j = tree.tree[j+1].pair-1;
					for (l = bp_j +1; l < j; l++){
						// Hosna: April 24, 2007
						// correct case 2 such that a multi-pseudoknotted
						// loop would not be treated as case 2
						cand_pos_t Bp_lj = tree.Bp(l+1,j+1)-1;

						if (Bp_lj >= 0 && Bp_lj<n){
							int sum = get_BE(bp_j,j,tree.tree[Bp_lj+1].pair-1,Bp_lj,tree) + get_WMBP(i,l,tree) + get_WI(l+1,Bp_lj-1,tree);
							if (acc > sum){
								acc = sum;
								best_l = l;
							}
						}
					}
					tmp = PB_penalty + acc;
					if (tmp < min){
						min = tmp;
						best_row = 1;
					}

				}
				// case WMBP
				tmp = get_WMBP(i,j,tree);
				if (tmp < min){
					min = tmp;
					best_row = 2;
				}

				switch (best_row)
				{
					case 1:
						if (best_l > -1){
							insert_node(i,best_l,P_WMBP);
							insert_node(best_l +1,tree.Bp(best_l+1,j+1)-1-1,P_WI);
							insert_node(tree.tree[j+1].pair-1,tree.tree[tree.Bp(best_l+1,j+1)-1+1].pair-1, P_BE);
						}
						break;
					case 2:
						insert_node(i,j,P_WMBP);
						break;
				}
			}
				break;
			case P_WMBP:
			{
				int i = cur_interval->i;
				int j = cur_interval->j;
				if (i >= j)
					return;
				int tmp = INF, best_l = -1, best_row = -1, min = INF;

				// case 1
				if(tree.tree[j+1].pair < 0 && tree.tree[i+1].pair >= 0){
					int l, l1 = -1;
					// Hosna: June 29, 2007
					// if j is inside i's arc then the l should be
					// less than j not bp(i)
					// check with Anne
					// Hosna: July 5th, 2007:
					// if we have bp(i)> j then we should not have come to the WMBP
					for (l = i+1; l < j; l++){
						// Hosna, April 9th, 2007
						// checking the borders as they may be negative
						// Hosna: July 5th, 2007:
						cand_pos_t bp_il = tree.bp(i+1,l+1)-1;
						// removed bp(l)<0 as VP should handle that
						if( bp_il >= 0 &&  bp_il < n && l+TURN <= j){
							// Hosna: April 19th, 2007
							// the chosen l should be less than border_b(i,j)
							int BE_energy = get_BE(i,tree.tree[i+1].pair-1,bp_il,tree.tree[bp_il+1].pair-1,tree);
							int WI_energy = get_WI(bp_il +1,l-1,tree);
							int VP_energy = get_VP(l,j,tree);
							int sum = BE_energy + WI_energy + VP_energy;
							if (tmp > sum){
								tmp = sum;
								l1 = l;
							}
						}
					}
					tmp = 2*PB_penalty + tmp;
					if (tmp < min){
						min = tmp;
						best_row = 1;
						best_l = l1;
					}
				}
				// case 3
				cand_pos_t b_ij = tree.b(i+1,j+1)-1;
				if (tree.tree[j+1].pair < 0){
					int l, acc = INF;
					int l3 = -1;
					for (l = i+1; l<j ; l++)	{
						// Hosna, April 6th, 2007
						// whenever we use get_borders we have to check for the correct values
						cand_pos_t B_lj = tree.B(l+1,j+1)-1;
						cand_pos_t Bp_lj = tree.Bp(l+1,j+1)-1;
						if (tree.tree[l+1].parent->index > -1 && B_lj >= 0 && B_lj < n && Bp_lj >= 0 && Bp_lj<n){
							// Hosna: April 19th, 2007
							// the chosen l should be less than border_b(i,j)
							if (b_ij >= 0 && b_ij < n && l < b_ij){
								// Hosna: June 29 2007
								// after going over the program with Cristina, we noticed that
								// l should be < B'(i,j)
								// Hosna: July 5th, 2007:
								// as long as we have i <= arc(l)< j we are fine
								if (i <= tree.tree[l+1].parent->index-1 && tree.tree[l+1].parent->index-1 < j && l+TURN <=j){
									int sum = get_BE(tree.tree[B_lj+1].pair-1,B_lj,tree.tree[Bp_lj+1].pair-1,Bp_lj,tree)+ get_WMBP(i,l-1,tree)+ get_VP(l,j,tree);
									if (acc > sum){
										acc = sum;
										l3 = l;
									}
								}
							}
						}
					}
					// Hosna: April 5th
					// after going over the WMB recurrence with Anne, we think we should add another p_b penalty
					// to the 3rd case ==> 2*P_b
					tmp = 2 *PB_penalty + acc;
					if (tmp < min){
						min = tmp;
						best_row = 3;
						best_l = l3;
					}
				}
				// case 4
				tmp = get_VP(i,j,tree) + PB_penalty;
				if (tmp < min){
					min = tmp;
					best_row = 4;
				}

				// case 5
				if(tree.tree[j+1].pair < j+1){
					int l, acc = INF;
					for(l = i+1; l<j; l++){
						// Hosna: April 18th, 2007
						// l and j should be in the same arc
						if (tree.tree[l+1].pair < 0 && tree.tree[l+1].parent->index > -1 && tree.tree[j+1].parent->index > -1 && tree.tree[l+1].parent->index == tree.tree[l+1].parent->index){
						tmp = get_WMBP(i,l,tree) + get_WI(l+1,j,tree);
							if (tmp < min){
								min = tmp;
								best_l = l;
								best_row = 5;
							}
						}
					}
				}

				switch (best_row)
				{
					case 1:
						if (best_l > -1){
							insert_node(i,tree.bp(i+1,best_l+1)-1,P_BE);
							insert_node(tree.bp(i+1,best_l+1)-1+1,best_l-1,P_WI);
							insert_node(best_l,j,P_VP);
						}
						break;
					case 3:
						if (best_l > -1){
							insert_node(i,best_l -1,P_WMBP);
							insert_node(best_l,j,P_VP);
							insert_node(tree.tree[tree.B(best_l+1,j+1)].pair-1,tree.tree[tree.Bp(best_l+1,j+1)].pair-1,P_BE);
						}
						break;
					case 4:
						insert_node(i,j,P_VP);
						break;
					case 5:
						if (best_l > -1){
							insert_node(i,best_l,P_WMBP);
							insert_node(best_l +1,j,P_WI);
						}
						break;
				}

			}
				break;
			case P_VP:
			{
				cand_pos_t i = cur_interval->i;
				cand_pos_t j = cur_interval->j;
				if (i>=j){
					return;
				}
				f[i].pair = j;
				f[j].pair = i;
				this->structure[i] = '[';
				this->structure[j] = ']';
				//printf("----> original VP: adding (%d,%d) <-------\n",i,j);
				f[i].type = P_VP;
				f[j].type = P_VP;

				int min = INF, tmp = INF, best_ip = INF, best_jp = INF, best_row = -1, best_r = INF;
				cand_pos_t Bp_ij = tree.Bp(i+1,j+1)-1;
				cand_pos_t B_ij = tree.B(i+1,j+1)-1;
				cand_pos_t b_ij = tree.b(i+1,j+1)-1;
				cand_pos_t bp_ij = tree.bp(i+1,j+1)-1;
				//case 1
				// Hosna April 9th, 2007
				// need to check the borders as they may be negative
				if((tree.tree[i+1].parent->index-1) > 0 && (tree.tree[j+1].parent->index-1) < (tree.tree[i+1].parent->index-1) && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0){
					energy_t WI_ipus1_BPminus = get_WI(i+1,Bp_ij - 1,tree) ;
					energy_t WI_Bplus_jminus = get_WI(B_ij + 1,j-1,tree);
					tmp =   WI_ipus1_BPminus + WI_Bplus_jminus;
					if (tmp < min){
						min = tmp;
						best_row = 1;
					}
				}
				//case 2
				// Hosna April 9th, 2007
				// checking the borders as they may be negative
				if ((tree.tree[i+1].parent->index-1) < (tree.tree[j+1].parent->index-1) && (tree.tree[j+1].parent->index-1) > 0 && b_ij>= 0 && bp_ij >= 0 && Bp_ij < 0){
					energy_t WI_i_plus_b_minus = get_WI(i+1,b_ij - 1,tree);
					energy_t WI_bp_plus_j_minus = get_WI(bp_ij + 1,j-1,tree);
					tmp = WI_i_plus_b_minus + WI_bp_plus_j_minus;
					if (tmp < min){
						min = tmp;
						best_row = 2;
					}
				}
				//case 3
				// Hosna April 9th, 2007
				// checking the borders as they may be negative
				if((tree.tree[i+1].parent->index-1) > 0 && (tree.tree[j+1].parent->index-1) > 0 && Bp_ij >= 0 && B_ij >= 0  && b_ij >= 0 && bp_ij>= 0){
					energy_t WI_i_plus_Bp_minus = get_WI(i+1,Bp_ij - 1,tree);
					energy_t WI_B_plus_b_minus = get_WI(B_ij + 1,b_ij - 1,tree);
					energy_t WI_bp_plus_j_minus = get_WI(bp_ij +1,j - 1,tree);
					tmp = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
					if (tmp < min){
						min = tmp;
						best_row = 3;
					}
				}
				//case 4
				pair_type ptype_closingip1jm1 = pair[S_[i+1+1]][S_[j+1+1]];
				if(tree.tree[i+1+1].pair < 0 && tree.tree[j-1+1].pair < 0 && ptype_closingip1jm1){
					tmp = get_e_stP(i,j)+ get_VP(i+1,j-1,tree);
					if (tmp < min){
						min = tmp;
						best_row = 4;
					}
				}
				
				cand_pos_t min_borders = std::min((cand_pos_tu) Bp_ij, (cand_pos_tu) b_ij);
				cand_pos_t edge_i = std::min(i+MAXLOOP+1,j-TURN-1);
				min_borders = std::min({min_borders,edge_i});
				for (cand_pos_t k = i+1; k < min_borders; ++k){
					// Hosna: April 20, 2007
					// i and ip and j and jp should be in the same arc
					// it should also be the case that [i+1,ip-1] && [jp+1,j-1] are empty regions
					if (tree.tree[k+1].pair < FRES_RESTRICTED_MIN && (tree.tree[i+1].parent->index == tree.tree[k+1].parent->index) && (tree.up[(k+1)-1] > ((k+1)-(i+1)-1))){
						// Hosna, April 9th, 2007
						// whenever we use get_borders we have to check for the correct values
						cand_pos_t max_borders = std::max(bp_ij,B_ij)+1;
						cand_pos_t edge_j = k+j-i-MAXLOOP-2;
						max_borders = std::max({max_borders,edge_j});
						for (cand_pos_t l = j-1; l > max_borders ; --l){
							pair_type ptype_closingkj = pair[S_[k+1]][S_[l+1]];
							if (tree.tree[l+1].pair < FRES_RESTRICTED_MIN && ptype_closingkj>0 && (tree.up[(j+1)-1] > ((j+1)-(l+1)-1))){
								// Hosna: April 20, 2007
								// i and ip and j and jp should be in the same arc
								if (tree.tree[j+1].parent->index == tree.tree[l+1].parent->index){
									tmp = get_e_intP(i,k,l,j) + get_VP(k,l,tree);
									if (tmp < min){
										min = tmp;
										best_row = 5;
										best_ip = k;
										best_jp = l;
									}
								}
							}
						}
					}
				}
				//case 6
				// Hosna April 9th, 2007
				// checking the borders as they may be negative numbers
				cand_pos_t min_Bp_j = j;
				if (Bp_ij > 0 && Bp_ij < n && Bp_ij < min_Bp_j){
					min_Bp_j = Bp_ij;
				}
				for (cand_pos_t r = i+1; r < min_Bp_j ; r++){
					if (tree.tree[r+1].pair < 0){
						// Hosna: July 5th, 2007
						// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
						// in a multiloop that spans a band
						tmp = get_WIP(i+1,r-1,tree) + get_VPP(r,j-1,tree) + ap_penalty + 2* bp_penalty;
						if (tmp < min){
							min = tmp;
							best_row = 6;
							best_r = r;
						}
					}
				}
				//case 7
				// Hosna April 9th, 2007
				// checking the borders as they may be negative numbers
				cand_pos_t max_i_bp = i;
				if (bp_ij > 0 && bp_ij < n && bp_ij > max_i_bp){
					max_i_bp = bp_ij;
				}
				for (cand_pos_t r = max_i_bp+1; r < j ; r++){
					if (tree.tree[r+1].pair < 0){
						// Hosna: July 5th, 2007
						// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
						// in a multiloop that spans a band
						tmp = get_VPP(i+1,r,tree) + get_WIP(r+1,j-1,tree) + ap_penalty + 2* bp_penalty;
						if (tmp < min){
							min = tmp;
							best_row = 7;
							best_r = r;
						}
					}
				}

				switch (best_row)
				{
					case 1:
						if (i+1 <= Bp_ij-1){
							insert_node(i+1,Bp_ij-1,P_WI);
						}
						if ( B_ij+1 <= j-1){
							insert_node( B_ij+1,j-1,P_WI);
						}
						break;
					case 2:
						if (i+1 <= b_ij-1){
							insert_node(i+1,b_ij-1,P_WI);
						}
						if (bp_ij+1 <= j-1){
							insert_node(bp_ij+1,j-1,P_WI);
						}
						break;
					case 3:
						if (i+1 <= Bp_ij-1){
							insert_node(i+1,Bp_ij-1,P_WI);
						}
						if (B_ij+1 <= b_ij-1){
							insert_node(B_ij+1,b_ij-1,P_WI);
						}
						if (bp_ij+1 <= j-1){
							insert_node(bp_ij+1,j-1,P_WI);
						}
						break;
					case 4:
						if (i+1 <= j-1){
							insert_node(i+1,j-1,P_VP);
						}
						break;
					case 5:
						if (best_ip <= best_jp){
							insert_node(best_ip,best_jp,P_VP);
						}
						break;
					case 6:
						if (i+1 <= best_r-1){
							insert_node(i+1,best_r-1,P_WIP);
						}
						if (best_r <= j-1){
							insert_node(best_r,j-1,P_VPP);
						}
						break;
					case 7:
						if (i+1 <= best_r){
							insert_node(i+1,best_r,P_VPP);
						}
						if (best_r+1 <= j-1){
							insert_node(best_r+1,j-1,P_WIP);
						}
						break;
				}
			}
				break;
			case P_VPP:
			{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (i >= j){
				return;
			}
			int min = INF, tmp = INF, best_r = INF, best_row = -1;
			cand_pos_t Bp_ij = tree.Bp(i+1,j+1)-1;
			cand_pos_t bp_ij = tree.bp(i+1,j+1)-1;

			//case 1
			// Hosna April 9th, 2007
			// checking the borders as they may be negative numbers
			cand_pos_t max_i_bp = i;
			if (bp_ij > 0 && bp_ij < n && bp_ij > max_i_bp){
				max_i_bp = bp_ij;
			}
			for (cand_pos_t r = max_i_bp+1; r < j; r++ ){
				if (tree.tree[r+1].pair < 0){
					tmp = get_VP(i,r,tree) + get_WIP(r+1,j,tree);
					if (tmp < min){
						min = tmp;
						best_row = 1;
						best_r = r;
					}
				}
			}
			//case 2:
			// Hosna April 9th, 2007
			// checking the borders as they may be negative numbers
			cand_pos_t min_Bp_j = j;
			if (Bp_ij > 0 && Bp_ij < n && Bp_ij < min_Bp_j){
				min_Bp_j = Bp_ij;
			}
			for (cand_pos_t r = i+1; r < min_Bp_j; r++){
				if (tree.tree[r+1].pair < 0){
					tmp = get_WIP(i,r-1,tree) + get_VP(r,j,tree);
					if (tmp < min){
						min = tmp;
						best_row = 2;
						best_r = r;
					}
				}
			}

			// Branch 3:
			// will have to recheck this empty code 
			for (cand_pos_t r = max_i_bp+1; r < j; r++ ){
				if (tree.tree[r+1].pair < 0 && tree.up[j+1] > ((j+1)-(r+1+1)+1)){
					tmp = get_VP(i,r,tree) + cp_penalty *(j-r); // check the (j-r) part
					if (tmp < min){
						min = tmp;
						best_row = 3;
						best_r = r;
					}
				}
			}

			// Branch 4:
			for (cand_pos_t r = i+1; r < min_Bp_j; r++){
				if (tree.tree[r+1].pair < 0 && tree.up[r-1+1] > ((r-1+1)-(i+1)+1)){
					tmp = cp_penalty * (r-i) + get_VP(r,j,tree);
					if (tmp < min){
						min = tmp;
						best_row = 4;
						best_r = r;
					}
				}
			}
			switch(best_row)
			{
				case 1:
					if (best_r != INF){
						if (i <= best_r){
							insert_node(i,best_r,P_VP);
						}
						if (best_r+1 <= j){
							insert_node(best_r +1,j,P_WIP);
						}
					}
					break;
				case 2:
					if (best_r != INF){
						if (i <= best_r-1){
							insert_node(i,best_r-1,P_WIP);
						}
						if (best_r <= j){
							insert_node(best_r,j,P_VP);
						}
					}
					break;
				case 3:
					if (best_r != INF){
						if (i <= best_r){
							insert_node(i,best_r,P_VP);
						}
					}
					break;
				case 4:
					if (best_r != INF){
						if (best_r <= j){
							insert_node(best_r,j,P_VP);
						}
					}
					break;
			}
		}
			break;
		case P_WI:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (i >= j){
				return;
			}
			int min = INF, tmp = INF, best_row = -1, best_t= -1;

				// Hosna: July 2nd, 2007
				// in branch 1 of WI, we can have a case like
				// ((..))((...))
				// such that both i and j are paired but we can chop them
				for (cand_pos_t t = i; t< j; t++){
					int wi_1 = get_WI(i,t,tree);
					int wi_2 = get_WI(t+1,j,tree);
					tmp = wi_1 + wi_2;
					if(tmp < min){
						min = tmp;
						best_row = 1;
						best_t = t;
					}
				}
	//		}
			//case 2
			if ((tree.tree[i+1].pair == j+1 && tree.tree[j+1].pair == i+1)||(tree.tree[i+1].pair < 0 && tree.tree[j+1].pair < 0)){
				// Hosna, April 16th, 2007
				// changed back to see if it will work fine
				// Hosna: April 19th, 2007
				// I think we should call the restricted version

				tmp = V->get_energy(i,j) + PPS_penalty;
				if(tmp < min){
					min = tmp;
					best_row = 2;
				}
			}
			//case 3
			// Hosna: April 20, 2007
			// removed the penalty of PPS

			// Hosna: July 5th, 2007
			// Anne said we should put PPS back
			// change PSM to PSP
			tmp = get_WMB(i,j,tree) + PSP_penalty + PPS_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 3;
			}
			switch (best_row)
			{
				case 1:
					if (best_t != -1){
						if (i <= best_t){
							insert_node(i,best_t,P_WI);
						}
						if (best_t+1 <= j){
							insert_node(best_t+1,j,P_WI);
						}
					}
					break;
				case 2:
					if (i < j){
						insert_node(i,j,LOOP);
					}
					break;
				case 3:
					if (i < j){
						insert_node(i,j,P_WMB);
					}
					break;
			}
		}
			break;
		case P_BE:
		{
			cand_pos_t i = cur_interval->i;
			cand_pos_t j = tree.tree[i+1].pair-1;
			cand_pos_t ip = cur_interval->j;
			cand_pos_t jp = tree.tree[ip+1].pair-1;
			if (i > ip || i > j || ip > jp || jp > j){
				return;
			}

			f[i].pair = j;
			f[j].pair = i;
			this->structure[i] = '(';
			this->structure[j] = ')';
			f[i].type = P_BE;
			f[j].type = P_BE;
			f[ip].pair = jp;
			f[jp].pair = ip;
			this->structure[ip] = '(';
			this->structure[jp] = ')';
			f[ip].type = P_BE;
			f[jp].type = P_BE;

			energy_t min = INF, tmp = INF;
			cand_pos_t best_row = -1, best_l = INF;
			//case 1
			if (tree.tree[i+1+1].pair == j-1+1){
				tmp = get_e_stP(i,j) + get_BE(i+1,j-1,ip,jp,tree);
				if(tmp < min){
					min = tmp;
					best_row = 1;
				}
			}
			for (cand_pos_t l = i+1; l<= ip ; l++){
				if (tree.tree[l+1].pair >= 0 && jp+1 <= tree.tree[l+1].pair && tree.tree[l+1].pair < j+1){
				int lp = tree.tree[l+1].pair-1;
				int il = index[i]+l-i;
				int lpj = index[lp]+j-lp;

				bool empty_region_il = (tree.up[(l+1)-1] >= l-i-1); //empty between i+1 and lp-1
				bool empty_region_lpj = (tree.up[(j+1)-1] >= j-lp-1); // empty between l+1 and ip-1
				bool weakly_closed_il = tree.weakly_closed(i+1+1,l-1+1); // weakly closed between i+1 and lp-1
				bool weakly_closed_lpj = tree.weakly_closed(lp+1+1,j-1+1); // weakly closed between l+1 and ip-1

				if (empty_region_il && empty_region_lpj){
					tmp = get_e_intP(i,l,lp,j)+ get_BE(l,lp,ip,jp,tree);
					if (min > tmp){
						min = tmp;
						best_row = 2;
						best_l = l;
					}
				}

				// case 3
				if (weakly_closed_il && weakly_closed_lpj){
					tmp = get_WIP(i+1,l-1,tree) + get_BE(l,lp,ip,jp,tree) + get_WIP(lp+1,j-1,tree);
					if (min > tmp){
						min = tmp;
						best_row = 3;
						best_l = l;
					}
				}

				// case 4
				if (weakly_closed_il && empty_region_lpj){
					// Hosna: July 5th, 2007
					// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
					// in a multiloop that spans a band
					tmp = get_WIP(i+1,l-1,tree) + get_BE(l,lp,ip,jp,tree) + c_penalty * (j-lp+1) + ap_penalty + 2* bp_penalty;
					if (min > tmp){
						min = tmp;
						best_row = 4;
						best_l = l;
					}
				}

				// case 5
				if (empty_region_il && weakly_closed_lpj){
					// Hosna: July 5th, 2007
					// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
					// in a multiloop that spans a band
					tmp = ap_penalty + 2* bp_penalty+ c_penalty * (l-i+1) + get_BE(l,lp,ip,jp,tree) + get_WIP(lp+1,j-1,tree);
					if (min > tmp){
						min = tmp;
						best_row = 5;
						best_l = l;
					}
				}
				}
			}
			switch(best_row)
			{
				case 1:
					if (i+1 <= ip){
						insert_node(i+1,ip,P_BE);
					}
					break;
				case 2:
					if (best_l != INF && best_l <= ip){
						insert_node(best_l,ip,P_BE);
					}
					break;
				case 3:
					if (best_l != INF){
						if (i+1 <= best_l-1){
							insert_node(i+1,best_l-1,P_WIP);
						}
						if (best_l <= ip){
							insert_node(best_l,ip,P_BE);
						}
						if (tree.tree[best_l+1].pair +1 <= j-1+1){
							insert_node(tree.tree[best_l+1].pair-1 +1,j-1,P_WIP);
						}
					}
					break;
				case 4:
					if (best_l != INF){
						if (i+1 <= best_l-1){
							insert_node(i+1,best_l-1,P_WIP);
						}
						if (best_l <= ip){
							insert_node(best_l,ip,P_BE);
						}
					}
					break;
				case 5:
					if (best_l != INF){
						if (best_l <= ip){
							insert_node(best_l,ip,P_BE);
						}
						if (tree.tree[best_l+1].pair +1 <= j-1+1){
							insert_node(tree.tree[best_l+1].pair -1 +1,j-1,P_WIP);
						}
					}
					break;
			}

		}
			break;
		case P_WIP:
		{
			int i = cur_interval->i;
			int j = cur_interval->j;
			if (i == j){
				return;
			}
			int min = INF, tmp = INF, best_row = -1, best_t = INF;
			//case 1
			if (tree.tree[i+1].pair < 0){
				tmp = get_WIP(i+1,j,tree) + cp_penalty;
				if(tmp < min){
					min = tmp;
					best_row = 1;
				}
			}
			//case 2
			if (tree.tree[j+1].pair < 0){
				tmp = get_WIP(i,j-1,tree) + cp_penalty;
				if (tmp < min){
					min = tmp;
					best_row = 2;
				}
			}
			//case 3
			for (cand_pos_t t = i; t <j; t++){
				tmp = get_WIP(i,t,tree) + get_WIP(t+1,j,tree);
				if (tmp < min){
					min = tmp;
					best_row = 3;
					best_t = t;
				}
			}
			//case 4
			pair_type ptype_closing = pair[S_[i+1]][S_[j+1]];
			if (tree.tree[i+1].pair == j+1 || (tree.tree[i+1].pair < 0 && tree.tree[j+1].pair < 0 && ptype_closing)){
				tmp = V->get_energy(i,j)+ bp_penalty;
				if (tmp < min){
					min = tmp;
					best_row = 4;
				}
			}
			//case 5
			tmp = get_WMB(i,j,tree) + PSM_penalty + bp_penalty;
			if (tmp < min){
				min = tmp;
				best_row = 5;
			}
			switch(best_row)
			{
				case 1:
					if (i+1 <= j){
						insert_node(i+1,j,P_WIP);
					}
					break;
				case 2:
					if (i <= j-1){
						insert_node(i,j-1,P_WIP);
					}
					break;
				case 3:
					if (best_t != INF){
						if (i <= best_t){
							insert_node(i,best_t,P_WIP);
						}
						if (best_t+1 <= j){
							insert_node(best_t +1,j,P_WIP);
						}
					}
					break;
				case 4:
					if (i < j){
						insert_node(i,j,LOOP);
					}
					break;
				case 5:
					if (i < j){
						insert_node(i,j,P_WMB);
					}
					break;
			}
		}
			break;

	//default:
		//printf("Should not happen!!!");
	}
}

void pseudo_loop::insert_node(int i, int j, char type)
{
	seq_interval *tmp;
    tmp = new seq_interval;
    tmp->i = i;
    tmp->j = j;
    tmp->type = type;
    tmp->next = stack_interval;
    stack_interval = tmp;

}

void pseudo_loop::set_stack_interval(seq_interval *stack_interval){
	this->stack_interval = stack_interval;
}
