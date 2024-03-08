#ifndef V_FINAL_H_
#define V_FINAL_H_

#include "h_struct.h"
#include "s_energy_matrix.h"
#include "VM_final.h"
#include "base_types.hh"
#include <vector>

class VM_final;
class V_final{
	public:
	// constructor
	V_final(cand_pos_t nb_nucleotides);
	~V_final();
	void setloops(s_energy_matrix *v, VM_final *vm);
	energy_t get_energy(cand_pos_t i, cand_pos_t j, sparse_tree &tree);

	char get_type (cand_pos_t i, cand_pos_t j);
    // return the type at v_final(i,j)

	s_energy_matrix *v;
	VM_final *vm;
	protected:
	std::vector<cand_pos_t> index;
	int *type;


};
#endif /*V_FINAL_H_*/
