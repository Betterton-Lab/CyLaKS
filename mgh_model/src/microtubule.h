#ifndef _MICROTUBULE_H
#define _MICROTUBULE_H
#include "tubulin.h"
#include <vector>

#ifndef _PARAMETERS_H
typedef struct system_parameters system_parameters;
#endif
#ifndef _SYSTEM_PROPERTIES_H
typedef struct system_properties system_properties;
#endif

class Microtubule{
	private:

	public:
		int index_;			// Index of MT in mt_list
		int polarity_;		// 0 for plus end on right; 1 for left
		int plus_end_;		// Index of plus_end in lattice
		int minus_end_;		// Index of minus_end in lattice
		int delta_x_;		// Direction motors step (+/- 1)

		int mt_index_adj_;	// Index of adjacent microtubule

		int n_sites_;		// Number of tubulin sites on MT
		int coord_;			// Absolute coordinate of left-most edge of MT

		double kbT_;
		double radius_;
		double height_;
		double eta_inverse_;	
		double site_size_; 
		double big_l_;				// in nm, overall length of MT
		double gamma_;				// in units of pN / nm ???


		Microtubule *neighbor_ = nullptr; //FIXME for 1+ neighbors 
	
		std::vector<Tubulin> lattice_;	// All tubulin sites 

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		Microtubule();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, int i_mt);

		void SetParameters();
		void GenerateLattice();

		void UpdateExtensions();
};
#endif
