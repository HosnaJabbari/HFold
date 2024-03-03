#include "pseudo_loop.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "constants.h"
#include "h_struct.h"
#include "h_externs.h"
#include "h_common.h"
#include "VM_final.h"
#include "V_final.h"

// Ian Wark July 19 2017
// constant that defines what fres[i].pair will be compared against (>=) for impossible cases
// set to -1 because >= 0 means there is already a base pair there,
// and -1 means restricted struture says there is no base pair there.
#define FRES_RESTRICTED_MIN -1

pseudo_loop::pseudo_loop(std::string seq,char *sequence, char* restricted, V_final *V, VM_final *VM, vrna_param_t *params)
{
	this->sequence = sequence;
	this->seq = seq;
	this->restricted = restricted;
	this->V = V;
	this->VM = VM;
	params_ = params;
	S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);
	make_pair_matrix();
    allocate_space();
}

void pseudo_loop::allocate_space()
{
    int i;
    nb_nucleotides = strlen(sequence);
    needs_computation = 0; // Hosna, March 14, 2012 I need to remove this variable!! to make everything a lot faster

    index = new int [nb_nucleotides];
    int total_length = (nb_nucleotides *(nb_nucleotides+1))/2;
    index[0] = 0;
    for (int i=1; i < nb_nucleotides; i++)
        index[i] = index[i-1]+nb_nucleotides-i+1;

    WI = new int [total_length];
    if (WI == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < total_length; i++) WI[i] = 0; // if i == j -> p_up

    weakly_closed = new int[total_length];
    if (weakly_closed == NULL) giveup ("Cannot allocate memory", "weakly_closed");
    for (i=0; i < total_length; i++) weakly_closed[i] = 0;

    not_paired_all = new int[total_length];
    if (not_paired_all == NULL) giveup ("Cannot allocate memory", "not_paired_all");
    for (i=0; i < total_length; i++) not_paired_all[i] = 0;

    VP = new int[total_length];
    if (VP == NULL) giveup ("Cannot allocate memory", "VP");
    for (i=0; i < total_length; i++) VP[i] = INF;

    WMB = new int[total_length];
    if (WMB == NULL) giveup ("Cannot allocate memory", "WMB");
    for (i=0; i < total_length; i++) WMB[i] = INF;

    WMBP = new int[total_length];
	if (WMBP == NULL) giveup("Cannot allocate memory","WMBP");
	for (i=0; i < total_length; i++) WMBP[i] = INF;

    WIP = new int[total_length];
    if (WIP == NULL) giveup ("Cannot allocate memory", "WIP");
    for (i=0; i < total_length; i++) WIP[i] = INF;


    VPP = new int[total_length];
    if (VPP == NULL) giveup ("Cannot allocate memory", "VPP");
    for (i=0; i < total_length; i++) VPP[i] = INF;

    BE = new int[total_length];
    if (BE == NULL) giveup ("Cannot allocate memory", "BE");
    for (i=0; i < total_length; i++) BE[i] = 0; //check

    border_bs = new int*[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) border_bs[i] = new int[nb_nucleotides];

    border_bps = new int*[nb_nucleotides];
    for(i = 0; i < nb_nucleotides; i++) border_bps[i] = new int[nb_nucleotides];


    int_sequence = new int[nb_nucleotides];
    if (int_sequence == NULL) giveup ("Cannot allocate memory", "energy");
    for (i=0; i < nb_nucleotides; i++) int_sequence[i] = nuc_to_int(sequence[i]);

}

pseudo_loop::~pseudo_loop()
{
    delete [] WI;
    delete [] WIP;
    delete [] VP;
    delete [] VPP;
    delete [] WMB;
    delete [] WMBP;

    delete [] BE;
    delete [] weakly_closed;
    delete [] not_paired_all;


    // Ian Wark July 21 2017
    // border_bs is array of arrays
    // need to delete sub arrays as well
    for(int i = 0; i < nb_nucleotides; i++) {
        delete [] border_bs[i];
        delete [] border_bps[i];
    }

    delete [] border_bs;
    delete [] border_bps;

    delete [] index;
    delete [] int_sequence;
}

void pseudo_loop::set_features(h_str_features *f){
	fres = f;
}

int pseudo_loop::is_weakly_closed(int i, int j){
	// base case: if i > j then the region is weakly closed
	if (i>j){
		return 1;
	}
	int ij = index[i]+j-i;
	if (weakly_closed[ij] == 1)
		return 1;
	return 0;
}


int pseudo_loop::is_empty_region(int i, int j){
	//base case: if i> j then the region is empty
	if (i>j){
		return 1;
	}
	int ij = index[i]+j-i;
	if (not_paired_all[ij] == 1){
		return 1;
	}
	return 0;
}

void pseudo_loop::initialize(){

	int i, j;

    //Hosna: before going further, we should fill up the weakly closed array
    detect_weakly_closed(fres, weakly_closed, nb_nucleotides, index);
    detect_not_paired_all(fres, not_paired_all, nb_nucleotides, index);
    detect_border_bs(fres,border_bs, nb_nucleotides);
    detect_border_bps(fres,border_bps, nb_nucleotides);

}

void pseudo_loop::compute_energies(int i, int j)
{

	// Hosna, April 18th, 2007
	// based on discussion with Anne, we changed WMB to case 2 and WMBP(containing the rest of the recurrences)

	//	if(debug){
	//		printf("calculating VP(%d,%d) \n",i,j);
	//	}
	compute_VP(i,j,fres); // Hosna, March 14, 2012, changed the positionof computing VP from after BE to befor WMBP


//	if(debug){
//		printf("calculating WMBP(%d,%d) \n",i,j);
//	}::compute_WI
	compute_WMBP(i,j,fres);
//	if(debug){
//		printf("calculating WMB(%d,%d) \n",i,j);
//	}
    compute_WMB(i,j,fres);
//	if(debug){
//		printf("calculating WI(%d,%d) \n",i,j);
//	}
    compute_WI(i,j,fres);
//    if(debug){
//		printf("calculating WIP(%d,%d) \n",i,j);
//	}
    compute_WIP(i,j,fres);
//    if(debug){
//		printf("calculating VPP(%d,%d) \n",i,j);
//	}
    compute_VPP(i,j,fres);
//	if(debug){
//		printf("calculating BE(%d,%d) \n",i,j);
//	}

	compute_BE(fres[j].pair,j,fres[i].pair,i,fres);

//    if (debug){
//    	printf("WI(%d,%d) = %d \n", i,j, get_WI(i,j));
//    	printf("WIP(%d,%d) = %d \n", i,j, get_WIP(i,j));
//    	printf("VPP(%d,%d) = %d \n", i,j, get_VPP(i,j));
//    	printf("BE(%d,%d,%d,%d) = %d \n", i,fres[i].pair,j,fres[j].pair, get_BE(i,fres[i].pair,j,fres[j].pair));
//    	printf("VP(%d,%d) = %d \n", i,j, get_VP(i,j));
//    	printf("WMBP(%d,%d) = %d \n", i,j, get_WMBP(i,j));
//    	printf("WMB(%d,%d) = %d \n", i,j, get_WMB(i,j));
//    }
}

void pseudo_loop::compute_WI(int i, int j , h_str_features *fres){
	int min = INF, m1 = INF, m2= INF, m3= INF;
	int ij = index[i]+j-i;
	if (WI[ij] != 0){ //calculated before
		return;
	}

	//base cases
	// if [i,j] is not weakly closed then WI[i,j] = INF
	if (is_weakly_closed(i,j) == 0){
		WI[ij] = INF;
		return;
	}

	// branch 4, one base
	if (i == j){
		WI[ij] = PUP_penalty;
		return;
	}
	// Hosna: Feb 12, 2007
	// changed this part to see if it works better

	// Hosna: Feb 16, 2007:
	// we don't need to check to see if i and j are inside an arc
	// because they are not in an arc in G but they will be in an arc in G'
	if (fres[i].arc != fres[j].arc){
		WI[ij] = INF;
		return;
	}

// Hosna: July 2nd, 2007
// in branch 1 of WI, we can have a case like
// ((..))((...))
// such that both i and j are paired but we can chop them

	// branch 1:
//	if (fres[i].pair < 0 && fres[j].pair < 0)
//	{
	int t;
	for (t = i; t< j; t++){
		int wi_1 = get_WI(i,t);
		int wi_2 = get_WI(t+1,j);
		int energy = wi_1 + wi_2;
		m1 = (m1 > energy)? energy : m1;
	}
//	}
	// branch 2:

	// Ian Wark July 19 2017
	// fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if ((fres[i].pair == j && fres[j].pair == i)
	||(fres[i].pair < FRES_RESTRICTED_MIN && fres[j].pair < FRES_RESTRICTED_MIN)){
		// Hosna, April 16th, 2007
		// changed back to see if it will work fine
		// Hosna: April 19th, 2007
		// I think we should call the restricted version

		int v_ener = (i>j)? INF: V->get_energy(i,j);
		m2 = v_ener + PPS_penalty;
	}
	// branch 3:
	// Hosna: April 20, 2007
	// removed the penalty of PPS

	// Hosna: July 5th, 2007
	// Anne said we should put PPS back
	// change PSM to PSP

	m3 = get_WMB(i,j) + PSP_penalty + PPS_penalty;


	min = MIN(m1,MIN(m2,m3));
	WI[ij] = min;
}


void pseudo_loop::compute_WI_pkonly(int i, int j , h_str_features *fres){
	int min = INF, m1 = INF, m2= INF, m3= INF;
	int ij = index[i]+j-i;
	if (WI[ij] != 0){ //calculated before
		return;
	}

	//base cases
	// if [i,j] is not weakly closed then WI[i,j] = INF
	if (is_weakly_closed(i,j) == 0){
		WI[ij] = INF;
		return;
	}

	// branch 4, one base
	if (i == j){
		WI[ij] = PUP_penalty;
		return;
	}
	if (fres[i].arc != fres[j].arc){
		WI[ij] = INF;
		return;
	}


	// branch 1:
	int t;
	for (t = i; t< j; t++){
		int wi_1 = get_WI(i,t);
		int wi_2 = get_WI(t+1,j);
		int energy = wi_1 + wi_2;
		m1 = (m1 > energy)? energy : m1;
	}
	// branch 2:
	if (fres[i].pair == j && fres[j].pair == i){

		int v_ener = (i>j)? INF: V->get_energy(i,j);
		m2 = v_ener + PPS_penalty;

	}
	// branch 3:
	m3 = get_WMB(i,j) + PSP_penalty + PPS_penalty;


	min = MIN(m1,MIN(m2,m3));
	WI[ij] = min;
}


void pseudo_loop::compute_VP(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (VP[ij] != INF){//has been calculated before
//		if (debug)
//		{
//			printf("VP(%d,%d) was calculated before ==> VP(%d,%d) = %d \n",i,j,i,j,VP[ij]);
//		}
		return;
	}
	// base cases:
	// a) i == j => VP[ij] = INF
	// b) [i,j] is a weakly_closed region => VP[ij] = INF
	// c) i or j is paired in original structure => VP[ij] = INF

	if (i == j || j-i<4 || weakly_closed[ij] == 1 || fres[i].pair >= FRES_RESTRICTED_MIN || fres[j].pair >= FRES_RESTRICTED_MIN || can_pair(int_sequence[i],int_sequence[j]) != 1)	{
		VP[ij] = INF;
		return;
	}
	else{
		int m1 = INF, m2 = INF, m3 = INF, m4= INF, m5 = INF, m6 = INF, m7 = INF, m8 = INF; //different branches
		//branchs:
		// 1) inArc(i) and NOT_inArc(j)
		// WI(i+1)(B'(i,j)-1)+WI(B(i,j)+1)(j-1)

		// Hosna April 9th, 2007
		// need to check the borders as they may be negative
		if(fres[i].arc > -1 && fres[j].arc == -1 && get_Bp(i,j) >= 0 && get_Bp(i,j)< nb_nucleotides && get_B(i,j) >= 0 && get_B(i,j) < nb_nucleotides){
			int Bp_i = get_Bp(i,j);
			int B_i = get_B(i,j);
			int WI_ipus1_BPminus = get_WI(i+1,Bp_i - 1) ;
			int WI_Bplus_jminus = get_WI(B_i + 1,j-1);
			m1 =   WI_ipus1_BPminus + WI_Bplus_jminus;
		}

		// 2) NOT_inArc(i) and inArc(j)
		// WI(i+1)(b(i,j)-1)+WI(b'(i,j)+1)(j-1)

		// Hosna April 9th, 2007
		// checking the borders as they may be negative
		if (fres[i].arc == -1 && fres[j].arc > -1 && get_b(i,j)>= 0 && get_b(i,j) < nb_nucleotides && get_bp(i,j) >= 0 && get_bp(i,j) < nb_nucleotides){
			int b_i = get_b(i,j);
			int bp_i = get_bp(i,j);
			int WI_i_plus_b_minus = get_WI(i+1,b_i - 1);
			int WI_bp_plus_j_minus = get_WI(bp_i + 1,j-1);
			m2 = WI_i_plus_b_minus + WI_bp_plus_j_minus;
		}

		// 3) inArc(i) and inArc(j)
		// WI(i+1)(B'(i,j)-1)+WI(B(i,j)+1)(b(i,j)-1)+WI(b'(i,j)+1)(j-1)

		// Hosna April 9th, 2007
		// checking the borders as they may be negative
		if(fres[i].arc > -1 && fres[j].arc > -1 && get_Bp(i,j) >= 0 && get_Bp(i,j) < nb_nucleotides && get_B(i,j) >= 0 && get_B(i,j) < nb_nucleotides && get_b(i,j) >= 0 && get_b(i,j) < nb_nucleotides && get_bp(i,j)>= 0 && get_bp(i,j) < nb_nucleotides){
			int Bp_i = get_Bp(i,j);
			int B_i = get_B(i,j);
			int b_i = get_b(i,j);
			int bp_i = get_bp(i,j);
			int WI_i_plus_Bp_minus = get_WI(i+1,Bp_i - 1);
			int WI_B_plus_b_minus = get_WI(B_i + 1,b_i - 1);
			int WI_bp_plus_j_minus = get_WI(bp_i +1,j - 1);
			m3 = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
		}

		// 4) NOT_paired(i+1) and NOT_paired(j-1) and they can pair together
		// e_stP(i,i+1,j-1,j) + VP(i+1)(j-1)
		if(fres[i+1].pair < FRES_RESTRICTED_MIN && fres[j-1].pair < FRES_RESTRICTED_MIN && can_pair(int_sequence[i+1],int_sequence[j-1])){
			m4 = get_e_stP(i,j)+ get_VP(i+1,j-1);
		}

		// 5) NOT_paired(r) and NOT_paired(rp)
		//  VP(i,j) = e_intP(i,ip,jp,j) + VP(ip,jp)
		int ip, jp;
		int max_borders;
		// Hosna, April 6th, 2007
		// whenever we use get_borders we have to check for the correct values
		int min_borders = 0; // what if both are negative
		if (get_Bp(i,j)> 0 && get_Bp(i,j) < nb_nucleotides && get_b(i,j) >0 && get_b(i,j) < nb_nucleotides){
			min_borders = MIN(get_Bp(i,j),get_b(i,j));
		}else if (get_b(i,j) > 0 && get_b(i,j) < nb_nucleotides && (get_Bp(i,j) < 0 || get_Bp(i,j) > nb_nucleotides)){
			min_borders = get_b(i,j);
		}else if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && (get_b(i,j) < 0 || get_b(i,j) > nb_nucleotides)){
			min_borders = get_Bp(i,j);
		}
//		printf("B'(%d,%d) = %d, b(%d,%d) = %d, min_borders = %d\n",i,j,get_Bp(i,j),i,j,get_b(i,j), min_borders);
		for (ip = i+1; ip < min_borders; ip++){
			// Hosna: April 20, 2007
			// i and ip and j and jp should be in the same arc
			// also it should be the case that [i+1,ip-1] && [jp+1,j-1] are empty regions

			// Ian Wark July 19 2017
            // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
            // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
			if (fres[ip].pair < FRES_RESTRICTED_MIN && (fres[i].arc == fres[ip].arc) && is_empty_region(i+1,ip-1) == 1){
				// Hosna, April 6th, 2007
				// whenever we use get_borders we have to check for the correct values
				max_borders= 0;
				if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides){
					max_borders = MAX(get_bp(i,j),get_B(i,j));
				}else if (get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides && (get_bp(i,j) < 0 || get_bp(i,j) > nb_nucleotides)){
					max_borders = get_B(i,j);
				}else if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && (get_B(i,j) < 0 || get_B(i,j) > nb_nucleotides)){
					max_borders = get_bp(i,j);
				}
//				printf("b'(%d,%d) = %d, B(%d,%d) = %d, max_borders = %d\n",i,j,get_bp(i,j),i,j,get_B(i,j), max_borders);
				for (jp = max_borders+1; jp < j ; jp++){
                    // Ian Wark July 19 2017
                    // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
                    // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
					if (fres[jp].pair < FRES_RESTRICTED_MIN && can_pair(int_sequence[ip],int_sequence[jp]) && is_empty_region(jp+1,j-1) == 1){
						// Hosna: April 20, 2007
						// i and ip and j and jp should be in the same arc
						if (fres[j].arc == fres[jp].arc ){
							int tmp = get_e_intP(i,ip,jp,j) + get_VP(ip,jp);

							m5 = std::min(m5,tmp);
						}
					}
				}
			}
		}

		// 6) VP(i,j) = WIP(i+1,r-1) + VPP(r,j-1)
		int r;
		// Hosna April 9th, 2007
		// checking the borders as they may be negative numbers
		int min_Bp_j = j;
		if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_Bp(i,j) < min_Bp_j){
			min_Bp_j = get_Bp(i,j);
		}
		for (r = i+1; r < min_Bp_j ; r++){
            // Ian Wark July 19 2017
            // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
            // otherwise it will create pairs in spots where the restricted structure says there should be no pairs

			if (fres[r].pair < FRES_RESTRICTED_MIN){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				int tmp = get_WIP(i+1,r-1) + get_VPP(r,j-1) + ap_penalty + 2*bp_penalty;
				m6 = std::min(m6,tmp);
			}
		}


		// 7) VP(i,j) = VPP(i+1,r) + WIP(r+1,j-1)
		// Hosna April 9th, 2007
		// checking the borders as they may be negative numbers
		int max_i_bp = i;
		if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
			max_i_bp = get_bp(i,j);
		}
		for (r = max_i_bp+1; r < j ; r++){

			if (fres[r].pair < FRES_RESTRICTED_MIN){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				int tmp = get_VPP(i+1,r) + get_WIP(r+1,j-1)+ ap_penalty + 2* bp_penalty;
				m7 = std::min(m7,tmp);
			}
		}

		//finding the min energy
		int min = MIN(MIN(m1,m8),MIN(m2,m3));
		min = (min > MIN(m4,m5))? MIN(m4,m5) : min;
		min = (min > MIN(m6,m7))? MIN(m6,m7) : min;

		VP[ij] = min;

	}
}

void pseudo_loop::compute_WMBP(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
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

	// Ian Wark July 19 2017
	// fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if ((fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair > j)
	||  (fres[j].pair >= FRES_RESTRICTED_MIN && fres[j].pair < i)
	||  (fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair < i )
	||  (fres[j].pair >= FRES_RESTRICTED_MIN && j < fres[j].pair)){
		WMB[ij] = INF;
		return;
	}
	else{
		int m1 = INF, m3 = INF, m4 = INF, m5 = INF;
		// if not paired(j) and paired(i) then
		// WMBP(i,j) = 2*Pb + min_{i<l<bp(i)}(BE(i,bp(i),b'(i,l),bp(b'(i,l)))+WI(b'+1,l-1)+VP(l,j))
		if(fres[j].pair < 0 && fres[i].pair >= 0){
			int tmp = INF, l, l_min=-1;
			// Hosna: June 29, 2007
			// if j is inside i's arc then the l should be
			// less than j not bp(i)
			// check with Anne
//			for (l = i+1; l < MIN(fres[i].pair,j); l++){
			// Hosna: July 5th, 2007:
			// if we have bp(i)> j then we should not have come to the WMBP
			for (l = i+1; l < j; l++){
				// Hosna, March 14, 2007
				// fixing the for loop

				// Hosna, April 9th, 2007
				// checking the borders as they may be negative
//				if(fres[l].pair < 0 && get_bp(i,l) >= 0 && get_bp(i,l) < nb_nucleotides && l+TURN <= j){
				// Hosna: July 5th, 2007:
				// removed bp(l)<0 as VP should handle that
				if(get_bp(i,l) >= 0 && get_bp(i,l) < nb_nucleotides && l+TURN <= j){
					int bp_i_l = get_bp(i,l);
					int BE_energy = get_BE(i,fres[i].pair,bp_i_l,fres[bp_i_l].pair);
					int WI_energy = get_WI(bp_i_l +1,l-1);
					int VP_energy = get_VP(l,j);
					int sum = BE_energy + WI_energy + VP_energy;
					if (tmp > sum){
						tmp = sum;
						l_min = l;
					}
				}
			}
			m1 = 2*PB_penalty + tmp;

		}

		// 3)
		if (fres[j].pair < 0){
			int l, temp = INF, l_min=-1;
			for (l = i+1; l<j ; l++)	{
				// Hosna, April 6th, 2007
				// whenever we use get_borders we have to check for the correct values

				if (fres[l].arc > -1 && get_B(l,j) >= 0 && get_B(l,j) < nb_nucleotides && get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
					// Hosna: April 19th, 2007
					// the chosen l should be less than border_b(i,j)
					if (get_b(i,j) >= 0 && get_b(i,j) < nb_nucleotides && l < get_b(i,j)){

						// Hosna: June 29 2007
						// after going over the program with Cristina, we noticed that
						// l should be < B'(i,j)
	//					if (l < get_Bp(i,j) && l+TURN <= j){

						// Hosna: July 5th, 2007:
						// as long as we have i <= arc(l)< j we are fine
						if (i <= fres[l].arc && fres[l].arc < j && l+TURN <=j){
							int sum = get_BE(fres[get_B(l,j)].pair,get_B(l,j),fres[get_Bp(l,j)].pair,get_Bp(l,j))+ get_WMBP(i,l-1)+ get_VP(l,j);
							if (temp > sum){
								temp = sum;
								l_min = l;
							}
						}
					}
				}
				// Hosna: April 5th
				// after going over the WMB recurrence with Anne, we think we should add another p_b penalty
				// to the 3rd case ==> 2*P_b
				m3 = 2*PB_penalty + temp;
			}
		}

		// 4) WMB(i,j) = VP(i,j) + P_b
		int temp = get_VP(i,j) + PB_penalty;
		if (temp < m4){
			m4 = temp;
		}
		// 5) WMB(i,j) = min_{i<l<j}{WMB(i,l)+WI(l+1,j)} if bp(j)<j
		// Hosna: Feb 5, 2007
		if(fres[j].pair < j){
			int l,l_min =-1;
			for(l = i+1; l<j; l++){
				// Hosna: March 14th, 2007
				// I think l cannot be paired

				// Hosna: April 18th, 2007
				// l and j should be in the same arc
				if (fres[l].pair < 0 && fres[l].arc > -1 && fres[j].arc > -1 && fres[j].arc == fres[l].arc){
					int temp = get_WMBP(i,l) + get_WI(l+1,j);
					if (temp < m5){
						m5 = temp;
						l_min = l;
					}

				}
			}
		}

		// get the min for WMB
		WMBP[ij] = MIN(MIN(m1,m3),MIN(m4,m5));
	}


}

void pseudo_loop::compute_WMB(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
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

    // Ian Wark July 19 2017
	// fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if ((fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair > j)
	 || (fres[j].pair >= FRES_RESTRICTED_MIN && fres[j].pair < i)
	 || (fres[i].pair >= FRES_RESTRICTED_MIN && fres[i].pair < i )
	 || (fres[j].pair >= FRES_RESTRICTED_MIN && j < fres[j].pair)){
		WMB[ij] = INF;
		return;
	}
	else{
		int m2 = INF, mWMBP = INF;
		// 2)
		if (fres[j].pair >= 0 && j > fres[j].pair){
			int l, l_min=-1;
			int bp_j = fres[j].pair;
			int temp = INF;
			for (l = (bp_j +1); (l < j); l++){
				// Hosna: April 24, 2007
				// correct case 2 such that a multi-pseudoknotted
				// loop would not be treated as case 2

				// Hosna: July 5th, 2007
				// this restriction was removed as it is not needed here
//				if (l > fres[i].pair){
	//				printf("l = %d and bp_l = %d \n",l,fres[l].pair);
					// Hosna April 9th,
					// checking the borders as they may be negative numbers

					// Hosna: July 5th, 2007:
					// we don't need to check that l is unpaired here
//					if (fres[l].pair < 0 && get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
					if (get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
						int sum = get_BE(bp_j,j,fres[get_Bp(l,j)].pair,get_Bp(l,j)) + get_WMBP(i,l) + get_WI(l+1,get_Bp(l,j)-1);
						if (l == 600 && i == 522 && j == 615) {
							int t = 0;
						}
						if (temp > sum){
							temp = sum;
							l_min = l;
						}
					}
//				}

			}
			m2 = PB_penalty + temp;
		}
		// check the WMBP value
		mWMBP =  get_WMBP(i,j);

		// get the min for WMB
		WMB[ij] = MIN(m2,mWMBP);
	}
}

void pseudo_loop::compute_WIP(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (WIP[ij] < INF/2){ // was calculated before
		return;
	}
	if (fres[i].arc != fres[j].arc || i == j || weakly_closed[ij]== 0){
		WIP[ij] = INF;
		return;
	}
	int m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;

    // Ian Wark July 19 2017
	// fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
	// otherwise it will create pairs in spots where the restricted structure says there should be no pairs

	// branch 1:
	if (fres[i].pair < FRES_RESTRICTED_MIN){
		m1 = get_WIP(i+1,j) + cp_penalty;
	}
	// branch 2:
	if (fres[j].pair < FRES_RESTRICTED_MIN){
		m2 = get_WIP(i,j-1) + cp_penalty;
	}
	//branch 3:
	int t;
	for (t = i; t <j; t++){
		int tmp = get_WIP(i,t) + get_WIP(t+1,j);
		if (tmp < m3){
			m3 = tmp;
		}
	}

	// branch 4:
	if (fres[i].pair == j
	|| (fres[i].pair < FRES_RESTRICTED_MIN && fres[j].pair < FRES_RESTRICTED_MIN && can_pair(int_sequence[i],int_sequence[j]))){

		m4 = V->get_energy(i,j)	+ bp_penalty;

	}

	// branch 5:
	m5 = get_WMB(i,j) + PSM_penalty + bp_penalty;

	WIP[ij] = MIN(MIN(m1,MIN(m2,m3)),MIN(m4,m5));

}

void pseudo_loop::compute_WIP_pkonly(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (WIP[ij] < INF/2){ // was calculated before
		return;
	}
	if (fres[i].arc != fres[j].arc || i == j || weakly_closed[ij]== 0){
		WIP[ij] = INF;
		return;
	}
	int m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;

    // Ian Wark July 19 2017
    // fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
    // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	// branch 1:
	if (fres[i].pair < FRES_RESTRICTED_MIN){
		m1 = get_WIP(i+1,j) + cp_penalty;
	}
	// branch 2:
	if (fres[j].pair < FRES_RESTRICTED_MIN){
		m2 = get_WIP(i,j-1) + cp_penalty;
	}
	//branch 3:
	int t;
	for (t = i; t <j; t++){
		int tmp = get_WIP(i,t) + get_WIP(t+1,j);
		if (tmp < m3){
			m3 = tmp;
		}
	}

	// branch 4:
	if (fres[i].pair == j && fres[j].pair == i){

		m4 = V->get_energy(i,j)	+ bp_penalty;

	}

	// branch 5:
	m5 = get_WMB(i,j) + PSM_penalty + bp_penalty;

	WIP[ij] = MIN(MIN(m1,MIN(m2,m3)),MIN(m4,m5));

}





void pseudo_loop::compute_VPP(int i, int j, h_str_features *fres){
	int ij = index[i]+j-i;
	if (VPP[ij] != INF){ // computed before
		return;
	}
	if (i == j  || this->is_weakly_closed(i,j)){
		VPP[ij] = INF;
		return;
	}
	int m1 = INF, m2 = INF;

	// Hosna: July 4th, 2007
	// After discussion with Anne, we figured out that we need to add
	// two more cases to VPP so that it can handle cases that in only one side
	// we have some structure and the other side is empty
	int m3 = INF, m4 = INF;

	//branch 1:
	int r=-1 ;
	// Hosna April 9th, 2007
	// checking the borders as they may be negative numbers
	int max_i_bp = i;
	if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
		max_i_bp = get_bp(i,j);
	}
	for (r = max_i_bp+1; r < j; r++ ){
        // Ian Wark July 19 2017
        // fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (fres[r].pair < FRES_RESTRICTED_MIN){
			int tmp = get_VP(i,r) + get_WIP(r+1,j);
			if (tmp < m1){
				m1 = tmp;
			}
		}
	}

	//branch 2:
	// Hosna April 9th, 2007
	// checking the borders as they may be negative numbers
	int min_Bp_j = j;
	if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_bp(i,j) < min_Bp_j){
		min_Bp_j = get_Bp(i,j);
	}
	for (r = i+1; r < min_Bp_j; r++){
        // Ian Wark July 19 2017
        // fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (fres[r].pair < FRES_RESTRICTED_MIN){
			int tmp = get_WIP(i,r-1) + get_VP(r,j);
			if (tmp < m2){
				m2 = tmp;
			}
		}
	}

	// Branch 3:
	for (r = max_i_bp+1; r < j; r++ ){
        // Ian Wark July 19 2017
        // fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (fres[r].pair < FRES_RESTRICTED_MIN && this->is_empty_region(r+1,j)){
			int tmp = get_VP(i,r) + (cp_penalty *(j-r)); // check the (j-r) part
			if (tmp < m3){
				m3 = tmp;
			}
		}
	}

	// Branch 4:

	for (r = i+1; r < min_Bp_j; r++){
        // Ian Wark July 19 2017
        // fres[i].pair < 0 changed to fres[i].pair < FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
		if (fres[r].pair < FRES_RESTRICTED_MIN && this->is_empty_region(i,r-1)){
			int tmp = (cp_penalty * (r-i)) + get_VP(r,j);
			if (tmp < m4){
				m4 = tmp;
			}
		}
	}
	int min_branches = m1;
	if (m2 < min_branches){
		min_branches = m2;
	}
	if (m3 < min_branches){
		min_branches = m3;
	}
	if (m4 < min_branches){
		min_branches = m4;
	}
	VPP[ij] = min_branches; //MIN(MIN(m1,m2),MIN(m3,m4));
}

void pseudo_loop::compute_BE(int i, int j, int ip, int jp, h_str_features * fres){


    // Ian Wark July 19 2017
    // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
    // otherwise it will create pairs in spots where the restricted structure says there should be no pairs
	if (!( i >= 0 && i <= ip && ip < jp && jp <= j && j < nb_nucleotides && fres[i].pair >= FRES_RESTRICTED_MIN && fres[j].pair >= FRES_RESTRICTED_MIN && fres[ip].pair >= FRES_RESTRICTED_MIN && fres[jp].pair >= FRES_RESTRICTED_MIN && fres[i].pair == j && fres[j].pair == i && fres[ip].pair == jp && fres[jp].pair == ip)){ //impossible cases
		return;
	}
	int iip = index[i]+ip-i;
	if (BE[iip] != 0){ // computed before
		return;
	}
	// base case: i.j and ip.jp must be in G
	if (fres[i].pair != j || fres[ip].pair != jp){
		BE[iip] = INF;
		return;
	}

	// base case:
	if(i == ip && j == jp && i<j){
		BE[iip] = 0;
		return;
	}

	int m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;
	// 1) bp(i+1) == j-1
	if (fres[i+1].pair == j-1){
		m1 = get_e_stP(i,j) + get_BE(i+1,j-1,ip,jp);

	}

	// cases 2-5 are all need an l s.t. i<l<=ip and jp<=bp(l)<j
	int l;
	for (l = i+1; l<= ip ; l++){
        // Ian Wark July 19 2017
        // fres[i].pair >= 0 changed to fres[i].pair >= FRES_RESTRICTED_MIN (which equals -1 at time of writing)
        // otherwise it will create pairs in spots where the restricted structure says there should be no pairs

		// Hosna: March 14th, 2007
		if (fres[l].pair >= FRES_RESTRICTED_MIN && jp <= fres[l].pair && fres[l].pair < j){
			// Hosna, March 15, 2007
			// since not_paired_all[i,l] includes i and l themselves
			// and in BE energy calculation we are looking for the oepn region (i,l)
			// we have to look at not_paired_all[i+1,l-1]
			int lp = fres[l].pair;
			int il = index[i]+l-i;
			int lpj = index[lp]+j-lp;
			// 2)
			// Hosna June 29, 2007
			// when we pass a stacked pair instead of an internal loop to e_int, it returns underflow,
			// so I am checking explicitely that we won't have stems instead of internal loop
			if (is_empty_region(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1 ){//&& !(ip == (i+1) && jp==(j-1)) && !(l == (i+1) && lp == (j-1))){
				int temp = get_e_intP(i,l,lp,j)+ get_BE(l,lp,ip,jp);
				if (m2 > temp){
					m2 = temp;
				}
			}

			// 3)
			if (is_weakly_closed(i+1,l-1) == 1 && is_weakly_closed(lp+1,j-1) == 1){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				int temp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp) + get_WIP(lp+1,j-1)+ ap_penalty + 2* bp_penalty;
				if (m3 > temp){
					m3 = temp;
				}
			}

			// 4)
			if (is_weakly_closed(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				int temp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp) + cp_penalty * (j-lp+1) + ap_penalty + 2*bp_penalty;
				if (m4 > temp){
					m4 = temp;
				}
			}

			// 5)
			if (is_empty_region(i+1,l-1) == 1 && is_weakly_closed(lp+1,j-1) == 1){
				// Hosna: July 5th, 2007
				// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
				// in a multiloop that spans a band
				int temp = ap_penalty + 2*bp_penalty + (cp_penalty * (l-i+1)) + get_BE(l,lp,ip,jp) + get_WIP(lp+1,j-1);
				if (m5 > temp){
					m5 = temp;
				}
			}
		}
	}

	// finding the min and putting it in BE[iip]
	BE[iip] = MIN(m1,MIN(MIN(m2,m3),MIN(m4,m5)));
}

int pseudo_loop::get_WI(int i, int j){
	if (i>j){
		return 0;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && WI[ij] == 0){
		//printf("get_WI(%d,%d), and we need to compute WI (i.e it's 0)!\n",i,j);
		compute_WI(i,j,fres);
	}
	 */
	//printf("get_WI(%d,%d), after computation its value = %d!\n",i,j, WI[ij]);
	return WI[ij];


}

// Hosna, May 1st, 2012
// I don't think we need specific getter function for pkonly case
/*
int pseudo_loop::get_WI_pkonly(int i, int j){
	if (i>j){
		return 0;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	//if (needs_computation == 1 && WI[ij] == 0){
		//printf("get_WI(%d,%d), and we need to compute WI (i.e it's 0)!\n",i,j);
	//	compute_WI_pkonly(i,j,fres);
	//}

	//printf("get_WI(%d,%d), after computation its value = %d!\n",i,j, WI[ij]);
	return WI[ij];


}
*/

int pseudo_loop::get_VP(int i, int j){
	// Hosna, March 16, 2012
	// two bases should be at least 3 bases apart

	if (j-i < TURN || i >= j || fres[i].pair >= 0 || fres[j].pair >= 0 || this->is_weakly_closed(i,j) == 1){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && VP[ij] == INF){
		//printf("get_VP(%d,%d), and we need to compute VP (i.e. it's INF)!\n",i,j);
		compute_VP(i,j,fres);
	}
	 */
	//printf("get_VP(%d,%d), after computation its value = %d!\n",i,j, VP[ij]);
	return VP[ij];

}
int pseudo_loop::get_WMB(int i, int j){
	// Hosna: July 6th, 2007
	// added impossible cases
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN ||(fres[i].pair >= 0 && fres[i].pair > j) || (fres[j].pair >= 0 && fres[j].pair < i) || (fres[i].pair >= 0 && fres[i].pair < i ) || (fres[j].pair >= 0 && j < fres[j].pair)){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && WMB[ij] == INF){
		//printf("get_WMB(%d,%d), and we need to compute WMB (i.e. it's INF)!\n",i,j);
		compute_WMB(i,j,fres);
	}
	 */
	//printf("get_WMB(%d,%d), after computation its value = %d!\n",i,j, WMB[ij]);
	return WMB[ij];
}

// Hosna: April 18th, 2007
// changed WMB to case 2 and WMBP
int pseudo_loop::get_WMBP(int i, int j){
	// Hosna: July 6th, 2007
	// added impossible cases
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i< TURN || (fres[i].pair >= 0 && fres[i].pair > j) || (fres[j].pair >= 0 && fres[j].pair < i) || (fres[i].pair >= 0 && fres[i].pair < i ) || (fres[j].pair >= 0 && j < fres[j].pair)){
		return INF;
	}
	int ij = index[i]+j-i;
	// Hosna, May 1st , 2012
	// these parts are not needed any more
	/*
	if (needs_computation == 1 && WMBP[ij] == INF){
		//printf("get_WMBP(%d,%d), and we need to compute WMBP (i.e. it's INF)!\n",i,j);
		compute_WMBP(i,j,fres);
	}
	 */
	//printf("get_WMBP(%d,%d), after computation its value = %d!\n",i,j, WMBP[ij]);
	return WMBP[ij];
}

int pseudo_loop::get_BE(int i, int j, int ip, int jp){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i>= TURN && i >= 0 && i <= ip && ip < jp && jp <= j && j < nb_nucleotides && fres[i].pair >=0 && fres[j].pair >= 0 && fres[ip].pair >= 0 && fres[jp].pair >= 0 && fres[i].pair == j && fres[j].pair == i && fres[ip].pair == jp && fres[jp].pair == ip){
		if(i == ip && j == jp && i<j){
			return 0;
		}
		int iip = index[i]+ip-i;

		return BE[iip];
	}else{
		return INF;
	}
}

int pseudo_loop::get_WIP(int i, int j){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN || i >= j || this->is_weakly_closed(i,j) != 1){
		return INF;
	}
	int ij = index[i]+j-i;

	return WIP[ij];
}

int pseudo_loop::get_VPP(int i, int j){
	// Hosna, March 16, 2012,
	// i and j should be at least 3 bases apart
	if (j-i < TURN || i >= j || this->is_weakly_closed(i,j) == 1){
		return INF;
	}
	int ij = index[i]+j-i;

	return VPP[ij];

}

// PRE: i< j
int pseudo_loop::get_b(int i,int j){
	// Hosna, April 5th, 2007
	if (i > j){
		return INF;
	}
	int border = MIN(border_bs[j][i],INF);
	return border;
}

// PRE: i<j
int pseudo_loop::get_bp(int i,int j){
	// Hosna, April 5th, 2007
	if (i > j ){
		return -1;
	}
	int border = MAX(border_bps[j][i],-1);
	return border;
}
//PRE: i<j
int pseudo_loop::get_B(int i,int j){
	// Hosna, April 5th, 2007
	if (i > j ){
		return -1;
	}
	int border = MAX(border_bs[i][j],-1);
	return border;
}
//PRE: i<j
int pseudo_loop::get_Bp(int i,int j){
	// Hosna, April 5th, 2007
	if (i > j ){
		return INF;
	}
	int border = MIN(border_bps[i][j],INF);
	return border;
}

PARAMTYPE pseudo_loop::compute_int(int i, int j, int k, int l, const paramT *params){

	const int ptype_closing = pair[S_[i]][S_[j]];
    return E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S_[k]][S_[l]]],S1_[i+1],S1_[j-1],S1_[k-1],S1_[l+1],const_cast<paramT *>(params));
}

int pseudo_loop::
get_e_stP(int i, int j){
	if (i+1 == j-1){ // TODO: do I need something like that or stack is taking care of this?
		return INF;
	}
	int ss = compute_int(i+1,j+1,i+1+1,j-1+1,params_);
	return (int)round(e_stP_penalty * (double)ss);
}

int pseudo_loop::get_e_intP(int i, int ip, int jp, int j){
	// Hosna Feb 12th, 2007:
	// this function is only being called in branch 5 of VP
	// and branch 2 of BE
	// in both cases regions [i,ip] and [jp,j] are closed regions
	int e_int = compute_int(i+1,j+1,ip+1,jp+1,params_);

	// Hosna April 3rd, 2007
	// based on the discussion with Anne, we decided to have
	// e_intP = 0.83 * e_int
	int energy = (int)round(e_intP_penalty * (double)e_int);
	return energy;
}

int pseudo_loop::get_energy(int i, int j){
	return get_WMB(i,j);
}

void pseudo_loop::back_track(char *structure, minimum_fold *f, seq_interval *cur_interval)
{
	this->structure = structure;
	this->f = f;
	needs_computation = 0;
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
				if (fres[j].pair >= 0 && j > fres[j].pair){
					int l, acc = INF;
					int bp_j = fres[j].pair;
					for (l = bp_j +1; l < j; l++){
						// Hosna: April 24, 2007
						// correct case 2 such that a multi-pseudoknotted
						// loop would not be treated as case 2

						// Hosna: July 5th, 2007:
						// We don't need this restriction
		//				if (l > fres[i].pair){
							// Hosna April 9th,
							// checking the borders as they may be negative numbers

							// Hosna: July 5th, 2007:
							// We don't need to check for l not being paired
		//					if (fres[l].pair < 0 && get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
							if (get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
								int sum = get_BE(bp_j,j,fres[get_Bp(l,j)].pair,get_Bp(l,j)) + get_WMBP(i,l) + get_WI(l+1,get_Bp(l,j)-1);
								if (acc > sum){
									acc = sum;
									best_l = l;
								}
							}
		//				}
					}
					tmp = PB_penalty + acc;
					if (tmp < min){
						min = tmp;
						best_row = 1;
					}

				}
				// case WMBP
				tmp = get_WMBP(i,j);
				if (tmp < min){
					min = tmp;
					best_row = 2;
				}

				switch (best_row)
				{
					case 1:
						if (best_l > -1){
							insert_node(i,best_l,P_WMBP);
							insert_node(best_l +1,get_Bp(best_l,j)-1,P_WI);
							insert_node(fres[j].pair,fres[get_Bp(best_l,j)].pair, P_BE);
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
				if(fres[j].pair < 0 && fres[i].pair >= 0){
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
						// removed bp(l)<0 as VP should handle that
						if(get_bp(i,l) >= 0 && get_bp(i,l) < nb_nucleotides && l+TURN <= j){
							// Hosna: April 19th, 2007
							// the chosen l should be less than border_b(i,j)
							int bp_i_l = get_bp(i,l);
							int BE_energy = get_BE(i,fres[i].pair,bp_i_l,fres[bp_i_l].pair);
							int WI_energy = get_WI(bp_i_l +1,l-1);
							int VP_energy = get_VP(l,j);
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
				if (fres[j].pair < 0){
					int l, acc = INF;
					int l3 = -1;
					for (l = i+1; l<j ; l++)	{
						// Hosna, April 6th, 2007
						// whenever we use get_borders we have to check for the correct values
						if (fres[l].arc > -1 && get_B(l,j) >= 0 && get_B(l,j) < nb_nucleotides && get_Bp(l,j) >= 0 && get_Bp(l,j)<nb_nucleotides){
							// Hosna: April 19th, 2007
							// the chosen l should be less than border_b(i,j)
							if (get_b(i,j) >= 0 && get_b(i,j)<nb_nucleotides && l < get_b(i,j)){
								// Hosna: June 29 2007
								// after going over the program with Cristina, we noticed that
								// l should be < B'(i,j)
								// Hosna: July 5th, 2007:
								// as long as we have i <= arc(l)< j we are fine
								if (i <= fres[l].arc && fres[l].arc < j && l+TURN <=j){
									int sum = get_BE(fres[get_B(l,j)].pair,get_B(l,j),fres[get_Bp(l,j)].pair,get_Bp(l,j))+ get_WMBP(i,l-1)+ get_VP(l,j);
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
				tmp = get_VP(i,j) + PB_penalty;
				if (tmp < min){
					min = tmp;
					best_row = 4;
				}

				// case 5
				if(fres[j].pair < j){
					int l, acc = INF;
					for(l = i+1; l<j; l++){
						// Hosna: April 18th, 2007
						// l and j should be in the same arc
						if (fres[l].pair < 0 && fres[l].arc > -1 && fres[j].arc > -1 && fres[l].arc == fres[j].arc){
						tmp = get_WMBP(i,l) + get_WI(l+1,j);
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
							insert_node(i,get_bp(i,best_l),P_BE);
							insert_node(get_bp(i,best_l)+1,best_l-1,P_WI);
							insert_node(best_l,j,P_VP);
						}
						break;
					case 3:
						if (best_l > -1){
							insert_node(i,best_l -1,P_WMBP);
							insert_node(best_l,j,P_VP);
							insert_node(fres[get_B(best_l,j)].pair,fres[get_Bp(best_l,j)].pair,P_BE);
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
				int i = cur_interval->i;
				int j = cur_interval->j;
				if (i>=j){
					return;
				}
				f[i].pair = j;
				f[j].pair = i;
				structure[i] = '[';
				structure[j] = ']';
				//printf("----> original VP: adding (%d,%d) <-------\n",i,j);
				f[i].type = P_VP;
				f[j].type = P_VP;

				int min = INF, tmp = INF, best_ip = INF, best_jp = INF, best_row = -1, best_r = INF;

				//case 1
				// Hosna April 9th, 2007
				// need to check the borders as they may be negative
				if(fres[i].arc > -1 && fres[j].arc == -1 && get_Bp(i,j) >= 0 && get_Bp(i,j)< nb_nucleotides && get_B(i,j) >= 0 && get_B(i,j) < nb_nucleotides){
					int Bp_i = get_Bp(i,j);
					int B_i = get_B(i,j);
					int WI_ipus1_BPminus = get_WI(i+1,Bp_i - 1) ;
					int WI_Bplus_jminus = get_WI(B_i + 1,j-1);
					tmp =   WI_ipus1_BPminus + WI_Bplus_jminus;
					if (tmp < min){
						min = tmp;
						best_row = 1;
					}
				}
				//case 2
				// Hosna April 9th, 2007
				// checking the borders as they may be negative
				if (fres[i].arc == -1 && fres[j].arc > -1 && get_b(i,j)>= 0 && get_b(i,j) < nb_nucleotides && get_bp(i,j) >= 0 && get_bp(i,j) < nb_nucleotides){
					int b_i = get_b(i,j);
					int bp_i = get_bp(i,j);
					int WI_i_plus_b_minus = get_WI(i+1,b_i - 1);
					int WI_bp_plus_j_minus = get_WI(bp_i + 1,j-1);
					tmp = WI_i_plus_b_minus + WI_bp_plus_j_minus;
					if (tmp < min){
						min = tmp;
						best_row = 2;
					}
				}
				//case 3
				// Hosna April 9th, 2007
				// checking the borders as they may be negative
				if(fres[i].arc > -1 && fres[j].arc > -1 && get_Bp(i,j) >= 0 && get_Bp(i,j) < nb_nucleotides && get_B(i,j) >= 0 && get_B(i,j) < nb_nucleotides && get_b(i,j) >= 0 && get_b(i,j) < nb_nucleotides && get_bp(i,j)>= 0 && get_bp(i,j) < nb_nucleotides){
					int Bp_i = get_Bp(i,j);
					int B_i = get_B(i,j);
					int b_i = get_b(i,j);
					int bp_i = get_bp(i,j);
					int WI_i_plus_Bp_minus = get_WI(i+1,Bp_i - 1);
					int WI_B_plus_b_minus = get_WI(B_i + 1,b_i - 1);
					int WI_bp_plus_j_minus = get_WI(bp_i +1,j - 1);
					tmp = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
					if (tmp < min){
						min = tmp;
						best_row = 3;
					}
				}
				//case 4
				if(fres[i+1].pair < 0 && fres[j-1].pair < 0 && can_pair(int_sequence[i+1],int_sequence[j-1])){
					tmp = get_e_stP(i,j)+ get_VP(i+1,j-1);
					if (tmp < min){
						min = tmp;
						best_row = 4;
					}
				}
				
				//case 5
				int ip, jp;
				// Hosna, April 9th, 2007
				// whenever we use get_borders we have to check for the correct values
				int min_borders = 0; // what if both are negative
				if (get_Bp(i,j)> 0 && get_Bp(i,j) < nb_nucleotides && get_b(i,j) >0 && get_b(i,j) < nb_nucleotides){
					min_borders = MIN(get_Bp(i,j),get_b(i,j));
				}else if (get_b(i,j) > 0 && get_b(i,j) < nb_nucleotides && (get_Bp(i,j) < 0 || get_Bp(i,j) > nb_nucleotides)){
					min_borders = get_b(i,j);
				}else if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && (get_b(i,j) < 0 || get_b(i,j) > nb_nucleotides)){
					min_borders = get_Bp(i,j);
				}
				for (ip = i+1; ip < min_borders; ip++){
					// Hosna: April 20, 2007
					// i and ip and j and jp should be in the same arc
					// it should also be the case that [i+1,ip-1] && [jp+1,j-1] are empty regions
					if (fres[ip].pair < 0 && fres[i].arc == fres[ip].arc && is_empty_region(i+1,ip-1)==1){
						// Hosna, April 9th, 2007
						// whenever we use get_borders we have to check for the correct values
						int max_borders= 0;
						if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides){
							max_borders = MAX(get_bp(i,j),get_B(i,j));
						}else if (get_B(i,j) > 0 && get_B(i,j) < nb_nucleotides && (get_bp(i,j) < 0 || get_bp(i,j) > nb_nucleotides)){
							max_borders = get_B(i,j);
						}else if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && (get_B(i,j) < 0 || get_B(i,j) > nb_nucleotides)){
							max_borders = get_bp(i,j);
						}
						for (jp = max_borders+1 ; jp < j ; jp++){
							if (fres[jp].pair < 0 && can_pair(int_sequence[ip],int_sequence[jp]) && is_empty_region(jp+1,j-1)==1){
								// Hosna: April 20, 2007
								// i and ip and j and jp should be in the same arc
								if (fres[j].arc == fres[jp].arc){
									tmp = get_e_intP(i,ip,jp,j) + get_VP(ip,jp);
									if (tmp < min){
										min = tmp;
										best_row = 5;
										best_ip = ip;
										best_jp = jp;
									}
								}
							}
						}
					}
				}
				//case 6
				int r;
				// Hosna April 9th, 2007
				// checking the borders as they may be negative numbers
				int min_Bp_j = j;
				if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_Bp(i,j) < min_Bp_j){
					min_Bp_j = get_Bp(i,j);
				}
				for (r = i+1; r < min_Bp_j ; r++){
					if (fres[r].pair < 0){
						// Hosna: July 5th, 2007
						// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
						// in a multiloop that spans a band
						tmp = get_WIP(i+1,r-1) + get_VPP(r,j-1) + ap_penalty + 2* bp_penalty;
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
				int max_i_bp = i;
				if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
					max_i_bp = get_bp(i,j);
				}
				for (r = max_i_bp+1; r < j ; r++){
					if (fres[r].pair < 0){
						// Hosna: July 5th, 2007
						// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
						// in a multiloop that spans a band
						tmp = get_VPP(i+1,r) + get_WIP(r+1,j-1) + ap_penalty + 2* bp_penalty;
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
						if (i+1 <= get_Bp(i,j)-1){
							insert_node(i+1,get_Bp(i,j)-1,P_WI);
						}
						if (get_B(i,j)+1 <= j-1){
							insert_node(get_B(i,j)+1,j-1,P_WI);
						}
						break;
					case 2:
						if (i+1 <= get_b(i,j)-1){
							insert_node(i+1,get_b(i,j)-1,P_WI);
						}
						if (get_bp(i,j)+1 <= j-1){
							insert_node(get_bp(i,j)+1,j-1,P_WI);
						}
						break;
					case 3:
						if (i+1 <= get_Bp(i,j)-1){
							insert_node(i+1,get_Bp(i,j)-1,P_WI);
						}
						if (get_B(i,j)+1 <= get_b(i,j)-1){
							insert_node(get_B(i,j)+1,get_b(i,j)-1,P_WI);
						}
						if (get_bp(i,j)+1 <= j-1){
							insert_node(get_bp(i,j)+1,j-1,P_WI);
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


			//case 1
			int r ;
			// Hosna April 9th, 2007
			// checking the borders as they may be negative numbers
			int max_i_bp = i;
			if (get_bp(i,j) > 0 && get_bp(i,j) < nb_nucleotides && get_bp(i,j) > max_i_bp){
				max_i_bp = get_bp(i,j);
			}
			for (r = max_i_bp+1; r < j; r++ ){
				if (fres[r].pair < 0){
					tmp = get_VP(i,r) + get_WIP(r+1,j);
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
			int min_Bp_j = j;
			if (get_Bp(i,j) > 0 && get_Bp(i,j) < nb_nucleotides && get_bp(i,j) < min_Bp_j){
				min_Bp_j = get_Bp(i,j);
			}
			for (r = i+1; r < min_Bp_j; r++){
				if (fres[r].pair < 0){
					tmp = get_WIP(i,r-1) + get_VP(r,j);
					if (tmp < min){
						min = tmp;
						best_row = 2;
						best_r = r;
					}
				}
			}

			// Hosna: July 4th, 2007
			// After discussion with Anne, we figured out that we need to add
			// two more cases to VPP so that it can handle cases that in only one side
			// we have some structure and the other side is empty

			// Branch 3:
			for (r = max_i_bp+1; r < j; r++ ){
				if (fres[r].pair < 0 && this->is_empty_region(r+1,j)){
					tmp = get_VP(i,r) + cp_penalty *(j-r); // check the (j-r) part
					if (tmp < min){
						min = tmp;
						best_row = 3;
						best_r = r;
					}
				}
			}

			// Branch 4:

			for (r = i+1; r < min_Bp_j; r++){
				if (fres[r].pair < 0 && this->is_empty_region(i,r-1)){
					tmp = cp_penalty * (r-i) + get_VP(r,j);
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
			//case 1
	//    	if (fres[i].pair < 0 && fres[j].pair < 0)
	//		{
				int t;
				for (t = i; t< j; t++){
					int wi_1 = get_WI(i,t);
					int wi_2 = get_WI(t+1,j);
					tmp = wi_1 + wi_2;
					if(tmp < min){
						min = tmp;
						best_row = 1;
						best_t = t;
					}
				}
	//		}
			//case 2
			if ((fres[i].pair == j && fres[j].pair == i)||(fres[i].pair < 0 && fres[j].pair < 0)){
				// Hosna, April 16th, 2007
				// changed back to see if it will work fine
				// Hosna: April 19th, 2007
				// I think we should call the restricted version
	//			tmp = V->get_energy_restricted(i,j,fres) + PPS_penalty;

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
			tmp = get_WMB(i,j) + PSP_penalty + PPS_penalty;
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
			int i = cur_interval->i;
			int j = fres[i].pair;
			int ip = cur_interval->j;
			int jp = fres[ip].pair;
			if (i > ip || i > j || ip > jp || jp > j){
				return;
			}

			f[i].pair = j;
			f[j].pair = i;
			structure[i] = '(';
			structure[j] = ')';
			f[i].type = P_BE;
			f[j].type = P_BE;
			f[ip].pair = jp;
			f[jp].pair = ip;
			structure[ip] = '(';
			structure[jp] = ')';
			f[ip].type = P_BE;
			f[jp].type = P_BE;

			int min = INF, tmp = INF, best_row = -1, best_l = INF;
			//case 1
			if (fres[i+1].pair == j-1){
				tmp = get_e_stP(i,j) + get_BE(i+1,j-1,ip,jp);
				if(tmp < min){
					min = tmp;
					best_row = 1;
				}
			}
			int l;
			for (l = i+1; l<= ip ; l++){
				if (fres[l].pair >= 0 && jp <= fres[l].pair && fres[l].pair < j){
				int lp = fres[l].pair;
				int il = index[i]+l-i;
				int lpj = index[lp]+j-lp;


				// case 2
	//			if (is_empty_region(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1){
				// Hosna June 29, 2007
				// when we pass a stacked pair instead of an internal loop to e_int, it returns underflow,
				// so I am checking explicitely that we won't have stems instead of internal loop
				if (is_empty_region(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1 ){
					tmp = get_e_intP(i,l,lp,j)+ get_BE(l,lp,ip,jp);
					if (min > tmp){
						min = tmp;
						best_row = 2;
						best_l = l;
					}
				}

				// case 3
				if (is_weakly_closed(i+1,l-1) == 1 && is_weakly_closed(lp+1,j-1) == 1){
					tmp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp) + get_WIP(lp+1,j-1);
					if (min > tmp){
						min = tmp;
						best_row = 3;
						best_l = l;
					}
				}

				// case 4
				if (is_weakly_closed(i+1,l-1) == 1 && is_empty_region(lp+1,j-1) == 1){
					// Hosna: July 5th, 2007
					// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
					// in a multiloop that spans a band
					tmp = get_WIP(i+1,l-1) + get_BE(l,lp,ip,jp) + c_penalty * (j-lp+1) + ap_penalty + 2* bp_penalty;
					if (min > tmp){
						min = tmp;
						best_row = 4;
						best_l = l;
					}
				}

				// case 5
				if (is_empty_region(i+1,l-1) == 1 && is_weakly_closed(lp+1,j-1) == 1){
					// Hosna: July 5th, 2007
					// After meeting with Anne and Cristina --> ap should have 2* bp to consider the biggest and the one that crosses
					// in a multiloop that spans a band
					tmp = ap_penalty + 2* bp_penalty+ c_penalty * (l-i+1) + get_BE(l,lp,ip,jp) + get_WIP(lp+1,j-1);
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
						if (fres[best_l].pair +1 <= j-1){
							insert_node(fres[best_l].pair +1,j-1,P_WIP);
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
						if (fres[best_l].pair +1 <= j-1){
							insert_node(fres[best_l].pair +1,j-1,P_WIP);
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
			if (fres[i].pair < 0){
				tmp = get_WIP(i+1,j) + cp_penalty;
				if(tmp < min){
					min = tmp;
					best_row = 1;
				}
			}
			//case 2
			if (fres[j].pair < 0){
				tmp = get_WIP(i,j-1) + cp_penalty;
				if (tmp < min){
					min = tmp;
					best_row = 2;
				}
			}
			//case 3
			int t;
			for (t = i; t <j; t++){
				tmp = get_WIP(i,t) + get_WIP(t+1,j);
				if (tmp < min){
					min = tmp;
					best_row = 3;
					best_t = t;
				}
			}
			//case 4
			if (fres[i].pair == j || (fres[i].pair < 0 && fres[j].pair < 0 && can_pair(int_sequence[i],int_sequence[j]))){
				tmp = V->get_energy(i,j)+ bp_penalty;
				if (tmp < min){
					min = tmp;
					best_row = 4;
				}
			}
			//case 5
			tmp = get_WMB(i,j) + PSM_penalty + bp_penalty;
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
