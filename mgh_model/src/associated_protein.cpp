#include "master_header.h"
#include "associated_protein.h"

AssociatedProtein::AssociatedProtein(){
}

void AssociatedProtein::Initialize(system_parameters *parameters, 
	system_properties *properties, int ID) {

	ID_ = ID;
	parameters_ = parameters;
	properties_ = properties;
	neighbor_sites_.resize(2*dist_cutoff_ + 1);
	PopulateBindingLookupTable();
}

void AssociatedProtein::PopulateBindingLookupTable(){

	double r_y = 35;		//dist between MTs in nm; static as of now
	binding_weight_lookup_.resize(dist_cutoff_ + 1);
	for(int i_dist = 0; i_dist <= dist_cutoff_; i_dist++){
		double r_x = (double)i_dist*site_size_;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0_;
		double weight = exp(-dr*dr*k_spring_/(2*kbT_));
		binding_weight_lookup_[i_dist] = weight;
//		printf("r: %g, dr: %g,  weight: %g\n", r, dr, weight);
	}
}

void AssociatedProtein::UpdateNeighborSites(){
	
	n_neighbor_sites_ = 0;
	int mt_length = parameters_->length_of_microtubule;
	Tubulin *site = GetActiveHeadSite();
	int i_site = site->index_;
	Microtubule *mt = site->mt_;
	Microtubule *adj_mt = mt->neighbor_;
	int offset = adj_mt->coord_ - mt->coord_;
	// Scan through all potential neighbor sites; only add unoccupied to list 
	int i_entry = 0;
	for(int distance = -dist_cutoff_; distance <= dist_cutoff_; distance++){
		int i_neighbor = i_site + distance - offset;
		// Start index at first bulk site (1) if i_neighbor is 0 or  negative
		if(i_neighbor <= 0){
			distance -= i_neighbor;
		}
		// End scan once last bulk site (mt_length - 2) has been checked
		else if(i_neighbor >= mt_length - 1){
			break;
		} 
		else{
			Tubulin *neighbor = &adj_mt->lattice_[i_neighbor];
			if(neighbor->occupied_ == false){
//				printf("site: %i, neighb: %i   ", i_site, i_neighbor);
//				printf("boop %i\n", i_entry);
				n_neighbor_sites_++;
				neighbor_sites_[i_entry] = neighbor;
				i_entry++;
			}
		}
	}
}

bool AssociatedProtein::NeighborExists(int x_dist){

	Tubulin *neighbor_site = GetNeighborSite(x_dist);
	if(neighbor_site == nullptr)
		return false;
	else
		return true;
}

void AssociatedProtein::UpdateExtension(){

	double r_y = 35;		// static as of now
	double r_x = site_size_*x_dist_;
	double r = sqrt(r_y*r_y + r_x*r_x);
	double extension = r - r_0_; 
	extension_ = extension;
}

void AssociatedProtein::ForceUnbind(){

	double ran = properties_->gsl.GetRanProb(); 
	if(ran < 0.5){
		// Remove xlink head from site
		site_one_->xlink_ = nullptr;
	   	site_one_->occupied_ = false;
		// Update xlink details
		site_one_ = nullptr;
	}
	else{
		// Remove xlink head from site
		site_two_->xlink_ = nullptr;
		site_two_->occupied_ = false;
		// Update xlink details	
		site_two_ = nullptr;	
	}
	extension_ = 0;
	heads_active_--;
}

double AssociatedProtein::GetAnchorCoordinate(){

	// If single bound, use that head; assume other's diffusion avgs out
	if(heads_active_ == 1){
		Tubulin *site = GetActiveHeadSite();
		int index = site->index_;
		int mt_coord = site->mt_->coord_;
		double site_coord = (double)(mt_coord + index);
		return site_coord;
	}
	// If double bound, use avg of head's indices 
	else if(heads_active_ == 2){
		int index_one = site_one_->index_;
		int mt_coord_one = site_one_->mt_->coord_;
		double coord_one = (double)(mt_coord_one + index_one);
		int index_two = site_two_->index_;
		int mt_coord_two = site_two_->mt_->coord_;
		double coord_two = (double)(mt_coord_two + index_two);
		double avg_coord = (coord_one + coord_two)/2;
//		printf("%g O_O\n", avg_coord);
		return avg_coord;
	}
	else{
		printf("not NOT cool bro ... cant get anchor index\n");
		exit(1);
	}
}

double AssociatedProtein::GetBindingWeight(Tubulin *neighbor){

	Tubulin *site = GetActiveHeadSite();
	Microtubule *mt = site->mt_;
	Microtubule *adj_mt = neighbor->mt_;
	if(adj_mt != mt->neighbor_){
		printf("why the microtubules tho (in assiociated protein GBW)\n");
		exit(1);
	}
	int offset = adj_mt->coord_ - mt->coord_;
	// Calculate distance (in x-dir.) between site and neighbor in # of sites 
	int i_site = site->index_;
	int i_neighbor = neighbor->index_;
	int x_dist = abs(i_neighbor - i_site + offset);
	// Get binding weight that corresponds to this x-distance
	double weight = binding_weight_lookup_[x_dist];
	return weight;	
}

Tubulin* AssociatedProtein::GetActiveHeadSite(){

	if(heads_active_ == 1){
		if(site_one_ != nullptr)
			return site_one_;
		else if(site_two_ != nullptr)
			return site_two_;
		else
			printf("what uVX\n");
	}
	else{
		printf("not cool bro...not single bound \n");
		exit(1);
	}
}

Tubulin* AssociatedProtein::GetNeighborSite(int x_dist){
	
	Tubulin *bound_site = GetActiveHeadSite();
	for(int i_entry = 0; i_entry < n_neighbor_sites_; i_entry++){
		Tubulin *neighbor_site = neighbor_sites_[i_entry];
		int offset = neighbor_site->mt_->coord_ - bound_site->mt_->coord_;
		int i_site = bound_site->index_;
		int i_neighb = neighbor_site->index_;
		int neighb_dist = abs(i_neighb - i_site + offset);
		if(neighb_dist == x_dist)
			return neighbor_site;
	}
	return nullptr; 
}
