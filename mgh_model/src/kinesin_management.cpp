#include "master_header.h"
#include "kinesin_management.h"

KinesinManagement::KinesinManagement(){
}

void KinesinManagement::Initialize(system_parameters *parameters, 
        system_properties *properties){

    parameters_ = parameters;
    properties_ = properties;

    GenerateMotors();	
    SetParameters();	
	InitiateLists();
}

void KinesinManagement::GenerateMotors(){

    int n_mts = parameters_->microtubules.count;
    int n_sites = parameters_->microtubules.length;
    // Since only one head has to be bound, the most that will ever
    // be needed (all single-bound) is the total number of sites 
    n_motors_ = n_mts*n_sites;
    motors_.resize(n_motors_);
    for(int ID = 0; ID < n_motors_; ID++){
        motors_[ID].Initialize(parameters_, properties_, ID);
    }
}

void KinesinManagement::SetParameters(){

    double delta_t = parameters_->delta_t;
	double site_size = parameters_->microtubules.site_size;
	/* Statistics for diffusion */
	double D_const = parameters_->motors.diffusion_const;
	double x_squared = (site_size/1000)*(site_size/1000); // convert to um^2
	tau_ = x_squared / (2 * D_const);
	p_diffuse_forward_ = delta_t / tau_;
	p_diffuse_backward_ = delta_t / tau_; 
	// Generate stepping rates based on extension of tether: rates are 
	// increased if stepping towards rest length, and reduced if stepping
	// away from rest length (which increases tether extension)
	rest_dist_ = motors_[0].rest_dist_;
	comp_cutoff_ = motors_[0].comp_cutoff_;
	dist_cutoff_ = motors_[0].dist_cutoff_;
	p_diffuse_to_teth_rest_.resize(2*dist_cutoff_ + 1);
	p_diffuse_from_teth_rest_.resize(2*dist_cutoff_ + 1);
	double kbT = parameters_->kbT;
	double r_0 = motors_[0].r_0_;
	double k_spring = motors_[0].k_spring_;
	double k_eff_slack = motors_[0].k_slack_;
	double r_y = parameters_->microtubules.y_dist / 2;
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		// Calculate tether length for this x_dist as well as +/- 1 it
		double r_x = x_dist_dub * site_size / 2;
		double r_x_fwd = (x_dist_dub + 1) *  site_size / 2;
		double r_x_bck = (x_dist_dub - 1) * site_size / 2;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double r_fwd = sqrt(r_y*r_y + r_x_fwd*r_x_fwd);
		double r_bck = sqrt(r_y*r_y + r_x_bck*r_x_bck);
		// Calculate extension of tether for given x_dist_dub
		double dr = r - r_0; 
		// Calculate extension if motor diffuses toward/away from rest
		double dr_toward; 
		double dr_from; 
		// If extended (pos extension), r_bck is towards rest
		if(dr >= 0){
			dr_toward = r_bck - r_0;
			dr_from = r_fwd - r_0;
		}
		// If compressed (neg extension), r_fwd is towards rest
		else{
			dr_toward = r_fwd - r_0;
			dr_from = r_bck - r_0; 
		}
		double dU_from, 
			   dU_to, 
			   weight_to, 
			   weight_from;
		if(x_dist_dub == 2*rest_dist_){
			// Weights according to Lanksy et al. 
			dU_from = (k_eff_slack/2)*(dr_from*dr_from - dr*dr);
			dU_to = (1/2)*(k_spring*dr_toward*dr_toward - k_eff_slack*dr*dr);
			weight_to = exp(-dU_to/(2*kbT));
		   	weight_from = exp(-dU_from/(2*kbT));
        }
		// use k_spring if extension is positive
		else if(dr > 0){
			dU_to = (k_spring/2)*(dr_toward*dr_toward - dr*dr);
			dU_from = (k_spring/2)*(dr_from*dr_from - dr*dr);
			weight_to = exp(-dU_to/(2*kbT));
			if(x_dist_dub >= 2*dist_cutoff_)
				weight_from = 0;
			else
		   		weight_from = exp(-dU_from/(2*kbT));
		}
		// otherwise, use k_eff to model 'slack' in the tether
		else{
			dU_to = (k_eff_slack/2)*(dr_toward*dr_toward - dr*dr);
			dU_from	= (k_eff_slack/2)*(dr_from*dr_from - dr*dr);
			if(x_dist_dub < 2*comp_cutoff_)
				weight_to = 0;
			else
		   		weight_to = exp(-dU_to/(2*kbT));
			if(x_dist_dub <= 2*comp_cutoff_)
				weight_from = 0;
			else
		   		weight_from = exp(-dU_from/(2*kbT));
		}
		double p_to = weight_to * delta_t / tau_;
		double p_from = weight_from * delta_t / tau_;
		p_diffuse_from_teth_rest_[x_dist_dub] = p_from;
		p_diffuse_to_teth_rest_[x_dist_dub] = p_to;
	}
	/* Statistics for KMC */
    double k_on_i = parameters_->motors.k_on_i;
    double c_motor = parameters_->motors.concentration;
	p_bind_i_ = k_on_i * c_motor * delta_t;
	double c_eff_teth = parameters_->motors.conc_eff_tether;
	p_bind_i_tethered_ = k_on_i * c_eff_teth * delta_t;
	double k_on_ii = parameters_->motors.k_on_ii;
	double c_eff_motor_bind = parameters_->motors.conc_eff_bind;
	p_bind_ii_ = k_on_ii * c_eff_motor_bind * delta_t;
	double k_off_i = parameters_->motors.k_off_i; 
	p_unbind_i_ = k_off_i * delta_t;
    double k_off = parameters_->motors.k_off_ii;
	p_unbind_ii_stepable_ = k_off * delta_t;
	double unbind_ratio = parameters_->motors.k_off_ratio;
	p_unbind_ii_stalled_ = k_off * delta_t / unbind_ratio;
	p_bind_ii_to_teth_.resize(2*dist_cutoff_ + 1);
	p_bind_ii_from_teth_.resize(2*dist_cutoff_ + 1);
	p_unbind_i_tethered_.resize(2*dist_cutoff_ + 1);
	p_unbind_ii_to_teth_.resize(2*dist_cutoff_ + 1);
	p_unbind_ii_from_teth_.resize(2*dist_cutoff_ + 1);
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		// Calculate tether length for this x_dist as well as +/- 1 it
		double r_x = x_dist_dub * site_size / 2;
		double r_x_fwd = (x_dist_dub + 1) *  site_size / 2;
		double r_x_bck = (x_dist_dub - 1) * site_size / 2;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double r_fwd = sqrt(r_y*r_y + r_x_fwd*r_x_fwd);
		double r_bck = sqrt(r_y*r_y + r_x_bck*r_x_bck);
		// Calculate extension of tether for given x_dist_dub
		double dr = r - r_0; 
		// Calculate extension if motor diffuses toward/away from rest
		double dr_toward; 
		double dr_from; 
		// If extended (pos extension), r_bck is towards rest
		if(dr >= 0){
			dr_toward = r_bck - r_0;
			dr_from = r_fwd - r_0;
		}
		// If compressed (neg extension), r_fwd is towards rest
		else{
			dr_toward = r_fwd - r_0;
			dr_from = r_bck - r_0; 
		}
		double U_tot,
			   dU_to,
			   dU_from;
		if(x_dist_dub == 2*rest_dist_){
			U_tot = (k_eff_slack/2)*dr*dr;
			dU_from = (k_eff_slack/2)*(dr_from*dr_from - dr*dr);
			dU_to = (1/2)*(k_spring*dr_toward*dr_toward - k_eff_slack*dr*dr);
        }
		// use k_spring if extension is positive
		else if(dr > 0){
			U_tot = (k_spring/2)*dr*dr;
			dU_to = (k_spring/2)*(dr_toward*dr_toward - dr*dr);
			dU_from = (k_spring/2)*(dr_from*dr_from - dr*dr);
		}
		// otherwise, use k_eff to model 'slack' in the tether
		else{
			U_tot = (k_eff_slack/2)*dr*dr;
			dU_to = (k_eff_slack/2)*(dr_toward*dr_toward - dr*dr);
			dU_from	= (k_eff_slack/2)*(dr_from*dr_from - dr*dr);
		}
		double weight_at = exp(-U_tot/(2*kbT));
		double weight_to = exp(-dU_to/(2*kbT));
		double weight_from = exp(-dU_from/(2*kbT));
		if(x_dist_dub >= 2*dist_cutoff_){
			weight_from = 0;
		}
		if(x_dist_dub <= 2*comp_cutoff_){
			weight_from = 0;
		}
		p_bind_ii_to_teth_[x_dist_dub] = weight_to * p_bind_ii_; 
		p_bind_ii_from_teth_[x_dist_dub] = weight_from * p_bind_ii_; 

		p_unbind_i_tethered_[x_dist_dub] = weight_at * p_unbind_i_;

		p_unbind_ii_to_teth_[x_dist_dub] = weight_to * p_unbind_ii_stalled_;
		p_unbind_ii_from_teth_[x_dist_dub] = 
			weight_from * p_unbind_ii_stalled_; 
	}
	double k_tether_free = parameters_->motors.k_tether_free;
	p_tether_free_ = k_tether_free * c_motor * delta_t;
	p_tether_bound_ = k_tether_free * c_eff_teth * delta_t;
	double k_untether_free = parameters_->motors.k_untether_free;
	p_untether_free_ = k_untether_free * delta_t;
    double motor_speed = parameters_->motors.velocity;
	p_step_ = motor_speed * delta_t / site_size;
	double k_failstep = parameters_->motors.failstep_rate; 
	p_failstep_ = k_failstep * delta_t;
	// Generate untethering and stepping rates for all tether extensions	
	// Everything is 2*dist_cutoff to allow for half-integer distances, 
	// so the 3rd entry will correspond to a distance of 1.5, etc. 
	double k_unteth = parameters_->motors.k_untether;
	double stall_force = parameters_->motors.stall_force;
	p_untether_bound_.resize(2*dist_cutoff_ + 1);
	p_step_to_teth_rest_.resize(2*dist_cutoff_ + 1);
	p_step_from_teth_rest_.resize(2*dist_cutoff_ + 1);
	p_failstep_to_teth_rest_.resize(2*dist_cutoff_ + 1); 
	p_failstep_from_teth_rest_.resize(2*dist_cutoff_ + 1); 
	double fail_coeff = p_failstep_ / p_step_; 
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		double r_x = x_dist_dub * site_size / 2;
		double r = sqrt(r_y*r_y + r_x*r_x);
		double dr = r - r_0;
		double cosine = r_x / r;
		// If extension is positive, treat tether like a spring
		if(dr >= 0){
			// Calculate spring potential energy (PE) for this x_dist
			double U_teth = (k_spring/2)*dr*dr;
			// Get untethering weight for this PE
			double unteth_weight = exp(U_teth/(2*kbT));
			p_untether_bound_[x_dist_dub] = k_unteth*unteth_weight*delta_t;
			// Stepping probability has force-dependent relation
			// (see sci advi suppl)
			double force = dr*k_spring;
			double coeff_to = 1 + cosine * (force / stall_force);
			double coeff_from = 1 - cosine * (force / stall_force); 
			double p_to = coeff_to * p_step_;
			double p_from = coeff_from * p_step_; 
			if(x_dist_dub >= 2*dist_cutoff_){
				p_step_to_teth_rest_[x_dist_dub] = p_to;
				p_step_from_teth_rest_[x_dist_dub] = 0;
				p_failstep_to_teth_rest_[x_dist_dub] = p_to*fail_coeff;
				p_failstep_from_teth_rest_[x_dist_dub] = 0;
			}
			else if(force < stall_force){
				p_step_to_teth_rest_[x_dist_dub] = p_to;
				p_step_from_teth_rest_[x_dist_dub] = p_from;
				p_failstep_to_teth_rest_[x_dist_dub] = p_to*fail_coeff;
				p_failstep_from_teth_rest_[x_dist_dub] = p_from*fail_coeff;
			}
			else{
				p_step_to_teth_rest_[x_dist_dub] = p_to; 
				p_step_from_teth_rest_[x_dist_dub] = 0;
				p_failstep_to_teth_rest_[x_dist_dub] = p_to*fail_coeff;
				p_failstep_from_teth_rest_[x_dist_dub] = 0;
			}
		}
		// Otherwise, use k_eff for slack
		else{
			double U_teth = (k_eff_slack/2)*dr*dr;
			double unteth_weight = exp(U_teth/(2*kbT));
			p_untether_bound_[x_dist_dub] = k_unteth*unteth_weight*delta_t;
			double force = -1 * dr * k_eff_slack;
			double coeff_to = 1 + cosine * (force / stall_force);
			double coeff_from = 1 - cosine * (force / stall_force);
			double p_to = coeff_to * p_step_;
			double p_from = coeff_from * p_step_;
			if(x_dist_dub < 2*comp_cutoff_){
				p_step_to_teth_rest_[x_dist_dub] = 0;
				p_step_from_teth_rest_[x_dist_dub] = 0;
				p_failstep_to_teth_rest_[x_dist_dub] = 0;
				p_failstep_from_teth_rest_[x_dist_dub] = 0;
			}
			else if(x_dist_dub == 2*comp_cutoff_){
				p_step_to_teth_rest_[x_dist_dub] = p_to;
				p_step_from_teth_rest_[x_dist_dub] = 0;
				p_failstep_to_teth_rest_[x_dist_dub] = p_to*fail_coeff;
				p_failstep_from_teth_rest_[x_dist_dub] = 0;
			}
			else if(force < stall_force){
				p_step_to_teth_rest_[x_dist_dub] = p_to;
				p_step_from_teth_rest_[x_dist_dub] = p_from;
				p_failstep_to_teth_rest_[x_dist_dub] = p_to*fail_coeff;
				p_failstep_from_teth_rest_[x_dist_dub] = p_from*fail_coeff;
			}
			else{
				p_step_to_teth_rest_[x_dist_dub] = 0;
				p_step_from_teth_rest_[x_dist_dub] = 0;
				p_failstep_to_teth_rest_[x_dist_dub] = 0;
				p_failstep_from_teth_rest_[x_dist_dub] = 0;
			}
		}
	}
	printf("For motors:\n");
	printf("  rest_dist is %g\n", rest_dist_);
	printf("  comp_cutoff is %i\n", comp_cutoff_);
	printf("  dist_cutoff is %i\n", dist_cutoff_);
	printf("  fail coeff is %g\n", fail_coeff); 	
}

void KinesinManagement::InitiateLists(){

	// One dimensional stuff
	free_tethered_list_.resize(n_motors_);
	bound_i_list_.resize(n_motors_);
	bound_i_bindable_list_.resize(n_motors_); 
	bound_ii_list_.resize(n_motors_);
	stepable_list_.resize(n_motors_);
	stalled_list_.resize(n_motors_);
	bound_ii_tethered_list_.resize(n_motors_);
	stepable_list_.resize(n_motors_);
	// Two dimensional stuff
	n_bound_ii_tethered_.resize(2*dist_cutoff_ + 1);
	n_stepable_to_teth_rest_.resize(2*dist_cutoff_ + 1);
	n_stepable_from_teth_rest_.resize(2*dist_cutoff_ + 1);
	n_stalled_to_teth_rest_.resize(2*dist_cutoff_ + 1);
	n_stalled_from_teth_rest_.resize(2*dist_cutoff_ + 1); 
	bound_ii_tethered_table_.resize(2*dist_cutoff_ + 1);
	stepable_to_rest_table_.resize(2*dist_cutoff_ + 1);
	stepable_from_rest_table_.resize(2*dist_cutoff_ + 1);
	stalled_to_rest_table_.resize(2*dist_cutoff_ + 1);
	stalled_from_rest_table_.resize(2*dist_cutoff_ + 1); 
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		bound_ii_tethered_table_[x_dist_dub].resize(n_motors_);
		stepable_to_rest_table_[x_dist_dub].resize(n_motors_);
		stepable_from_rest_table_[x_dist_dub].resize(n_motors_);
		stalled_to_rest_table_[x_dist_dub].resize(n_motors_);
		stalled_from_rest_table_[x_dist_dub].resize(n_motors_);
		n_bound_ii_tethered_[x_dist_dub] = 0;
		n_stepable_to_teth_rest_[x_dist_dub] = 0;
		n_stepable_from_teth_rest_[x_dist_dub] = 0;
		n_stalled_to_teth_rest_[x_dist_dub] = 0;
		n_stalled_from_teth_rest_[x_dist_dub] = 0;
	}
}

void KinesinManagement::UnboundCheck(Kinesin *motor){

    // Check motor
    if(motor->front_site_ != nullptr || motor->rear_site_ != nullptr){
        printf("Error with motor #%i: classified as unbound, ", motor->ID_);
        printf("but at least one of its heads is attached to tubulin\n");
        exit(1);
    }
    // Check motor_list
    int ID = motor->ID_;
    if(motors_[ID].front_site_ != nullptr
            || motors_[ID].rear_site_ != nullptr){
        printf("Error: motor #%i is out of sync with motor_list\n", ID);
        exit(1);
    }
}

void KinesinManagement::BoundCheck(Kinesin *motor){

    // Check motor
    if(motor->front_site_ == nullptr || motor->rear_site_ == nullptr){
        printf("Error with motor #%i: classified as bound, " , motor->ID_);
        printf("but at least one of its heads is not attached to tubulin\n");
        exit(1);
    }
    // Check motor_list
    int ID = motor->ID_;
    if(motors_[ID].front_site_ == nullptr
            || motors_[ID].rear_site_ == nullptr){
        printf("Error: motor #%i is out of sync with motor_list\n", ID);
        exit(1);
    }
}

int KinesinManagement::GetNumBoundUntethered(){

	int n_untethered = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin* motor = &motors_[i_motor];
		if(motor->heads_active_ > 0
		&& motor->tethered_ == false){
			n_untethered++;
		}
	}
	return n_untethered;
}

Kinesin* KinesinManagement::GetFreeMotor(){

	// Randomly pick a motor from the reservoir
	int i_motor = properties_->gsl.GetRanInt(n_motors_);
	Kinesin *motor = &motors_[i_motor];
	int attempts = 0;
	while(motor->heads_active_ > 0 
	|| motor->tethered_ == true){
		i_motor++;
		if(i_motor == n_motors_)
			i_motor = 0;
		motor = &motors_[i_motor];
		attempts++;
		if(attempts > n_motors_){
			printf("error in get free motor\n");
			exit(1);
		}
	}
	UnboundCheck(motor);
	return motor;
}

Kinesin* KinesinManagement::GetBoundUntetheredMotor(){

	Kinesin* untethered_list[n_motors_];
	int i_entry = 0;
	int n_untethered = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ > 0
		&& motor->tethered_ == false){
			untethered_list[i_entry] = motor;
			i_entry++;
			n_untethered++; 
		}
	}
	int i_motor = properties_->gsl.GetRanInt(n_untethered); 
	Kinesin *motor = untethered_list[i_motor];
	return motor;
}

void KinesinManagement::UpdateFreeTetheredList(){

	int i_entry = 0;
	int n_entries = 0; 
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 0
		&& motor->tethered_ == true){
			free_tethered_list_[i_entry] = motor;
			i_entry++;
			n_entries++;
		}
	}
	if(n_entries != n_free_tethered_){
		printf("something wrong in update_free_tethered.\n");
		exit(1);
	}
}

void KinesinManagement::UpdateBoundIList(){
	
	int i_entry = 0;
	int n_entries = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 1){
			bound_i_list_[i_entry] = motor;
			i_entry++;
			n_entries++;
		}
	}
	if(n_entries != n_bound_i_){
		printf("something wrong in update_pseudo_bound (motors)\n");
		exit(1);
	}
}

void KinesinManagement::UpdateBoundIBindableList(){

	n_bound_i_bindable_ = 0; 
	int i_entry = 0; 
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor]; 
		if(motor->heads_active_ == 1){
			Microtubule *mt = motor->mt_;
			Tubulin *bound_site = motor->GetActiveHeadSite();
			int i_site = bound_site->index_;
			int i_plus = mt->plus_end_;
			int i_minus = mt->minus_end_;
			int dx = mt->delta_x_;
			// Don't access lattice sites that don't exist
			if(i_site == i_plus){
				if(mt->lattice_[i_site - dx].occupied_ == false){
					bound_i_bindable_list_[i_entry] = motor;
					i_entry++;
					n_bound_i_bindable_++;
				}
			}
			else if(i_site == i_minus){
				if(mt->lattice_[i_site + dx].occupied_ == false){
					bound_i_bindable_list_[i_entry] = motor;
					i_entry++;
					n_bound_i_bindable_++;
				}
			}
			else{
				// Add motor if it has an unoccupied site to either side
				if(mt->lattice_[i_site + dx].occupied_ == false
				|| mt->lattice_[i_site - dx].occupied_ == false){
					bound_i_bindable_list_[i_entry] = motor;
					i_entry++;
					n_bound_i_bindable_++; 
				}
			}
		}
	}
}

void KinesinManagement::UpdateBoundIIList(){

	int i_entry = 0;
	int n_entries = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 2){
			if(motor->tethered_ == false){
				bound_ii_list_[i_entry] = motor;
				i_entry++;
				n_entries++;
			}
			else if(motor->xlink_->heads_active_ == 0){
				bound_ii_list_[i_entry] = motor;
				i_entry++;
				n_entries++;
			}
		}
	}
	if(n_entries != n_bound_ii_){
		printf("something wrong in update_unteth_list (motors 1D)");
		printf(" %i in statistics, %i entries tho\n", n_bound_ii_, 
				n_entries);
		exit(1);
	}
}

void KinesinManagement::UpdateStepableList(){

	n_stepable_ = 0; 
	int i_entry = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
        if(motor->heads_active_ == 2){
			if(motor->tethered_ == false){
				Microtubule *mt = motor->mt_;
				int plus_end = mt->plus_end_;
				int dx = mt->delta_x_;
				int i_front_site = motor->front_site_->index_;
				// Exclude plus_end (can't step off of MTs)
				if(i_front_site != plus_end){
					if(mt->lattice_[i_front_site + dx].occupied_ == false){
						stepable_list_[i_entry] = motor;
						i_entry++; 
						n_stepable_++;
					}
				}
			}
			else if(motor->xlink_->heads_active_ == 0){
				Microtubule *mt = motor->mt_;
				int plus_end = mt->plus_end_;
				int dx = mt->delta_x_;
				int i_front_site = motor->front_site_->index_;
				// Exclude plus_end (can't step off of MTs)
				if(i_front_site != plus_end){
					if(mt->lattice_[i_front_site + dx].occupied_ == false){
						stepable_list_[i_entry] = motor;
						i_entry++; 
						n_stepable_++;
					}
				}
			}
		}
	}
}

void KinesinManagement::UpdateStalledList(){

	n_stalled_ = 0; 
	int i_entry = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 2){
			if(motor->tethered_ == false){
				Microtubule *mt = motor->mt_;
				int plus_end = mt->plus_end_;
				int dx = mt->delta_x_;
				int i_front_site = motor->front_site_->index_;
				// Motors on plus-end are automatically considered stalled
				if(i_front_site == plus_end){
					stalled_list_[i_entry] = motor;
					i_entry++;
					n_stalled_++;
				}
				else if(mt->lattice_[i_front_site + dx].occupied_ == true){
					stalled_list_[i_entry] = motor;
					i_entry++; 
					n_stalled_++;
				}
			}
			else if(motor->xlink_->heads_active_ == 0){
				Microtubule *mt = motor->mt_;
				int plus_end = mt->plus_end_;
				int dx = mt->delta_x_;
				int i_front_site = motor->front_site_->index_;
				// Motors on plus-end are automatically considered stalled
				if(i_front_site == plus_end){
					stalled_list_[i_entry] = motor;
					i_entry++;
					n_stalled_++;
				}
				else if(mt->lattice_[i_front_site + dx].occupied_ == true){
					stalled_list_[i_entry] = motor;
					i_entry++; 
					n_stalled_++;
				}

			}
		}
	}
}

void KinesinManagement::UpdateBoundIITetheredList(){

	int i_entry = 0;
	int n_entries = 0;
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				bound_ii_tethered_list_[i_entry] = motor;
				i_entry++;
				n_entries++;
			}
		}
	}
	if(n_entries != n_bound_ii_tethered_tot_){
		printf("something wrong in update_teth_list (motors 1D)");
		printf(" %in statistics, %i entries tho\n", n_bound_ii_tethered_tot_, 
				n_entries);
		exit(1);
	}
}

void KinesinManagement::UpdateBoundIITetheredTable(){

	int i_entry[2*dist_cutoff_ + 1];
	int n_entries[2*dist_cutoff_ + 1]; 
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		i_entry[x_dist_dub] = 0;
		n_entries[x_dist_dub] = 0;
	}
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				motor->UpdateExtension();
				if(motor->tethered_ == true){
					int x_dist_dub = motor->x_dist_doubled_; 
					int index = i_entry[x_dist_dub]; 
					bound_ii_tethered_table_[x_dist_dub][index] = motor;
					i_entry[x_dist_dub]++;
					n_entries[x_dist_dub]++;
				}
			}
		}
	}
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		if(n_entries[x_dist_dub] != n_bound_ii_tethered_[x_dist_dub]){
			printf("something wrong in update_bound_tethered (motors)");
			printf("for ext %i, \n%i in stats but %i entries counted\n", 
					x_dist_dub, n_bound_ii_tethered_[x_dist_dub], 
					n_entries[x_dist_dub]);
			exit(1);
		}
	}
}

void KinesinManagement::UpdateStepableTetheredTables(){

	int i_entry_to[2*dist_cutoff_ + 1];
	int i_entry_from[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		i_entry_to[x_dist_dub] = 0;
		i_entry_from[x_dist_dub] = 0;
		n_stepable_to_teth_rest_[x_dist_dub] = 0; 	
		n_stepable_from_teth_rest_[x_dist_dub] = 0;
	}
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				motor->UpdateExtension();
				Microtubule *mt = motor->mt_;
				int plus_end = mt->plus_end_;
				int delta_x = mt->delta_x_;
				int dx_rest = motor->GetDirectionTowardRest();
				int i_front = motor->front_site_->index_;
				// end pausing
				if(i_front != plus_end
				&& motor->tethered_ == true){
					if(mt->lattice_[i_front + delta_x].occupied_ == false
					&& motor->AtCutoff() == false){
						int x_dist_dub = motor->x_dist_doubled_;
						// if MT's dx is towards rest, add to to_rest list
						if(delta_x == dx_rest){
							int index = i_entry_to[x_dist_dub];
							stepable_to_rest_table_[x_dist_dub][index] 
								= motor;	
							i_entry_to[x_dist_dub]++;
							n_stepable_to_teth_rest_[x_dist_dub]++;
						}
						// otherwise, add to from_rest list
						else if(delta_x == -1 * dx_rest){
							int index = i_entry_from[x_dist_dub];
							stepable_from_rest_table_[x_dist_dub][index] 
								= motor;
							i_entry_from[x_dist_dub]++;
							n_stepable_from_teth_rest_[x_dist_dub]++;
						}
						else{
							printf("hmmm??? teth tables\n");
							exit(1);
						}
					}
				}
			}
		}
	}
}

void KinesinManagement::UpdateStalledTetheredTables(){

	int i_entry_to[2*dist_cutoff_ + 1];
	int i_entry_from[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		i_entry_to[x_dist_dub] = 0;
		i_entry_from[x_dist_dub] = 0;
		n_stalled_to_teth_rest_[x_dist_dub] = 0; 	
		n_stalled_from_teth_rest_[x_dist_dub] = 0;
	}
	for(int i_motor = 0; i_motor < n_motors_; i_motor++){
		Kinesin *motor = &motors_[i_motor];
		if(motor->heads_active_ == 2
		&& motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				motor->UpdateExtension();
				Microtubule *mt = motor->mt_;
				int plus_end = mt->plus_end_;
				int delta_x = mt->delta_x_;
				int dx_rest = motor->GetDirectionTowardRest();
				int i_front = motor->front_site_->index_;
				if(i_front == plus_end){
					int x_dist_dub = motor->x_dist_doubled_;
					// if MT's dx is towards rest, add to to_rest list
					if(delta_x == dx_rest){
						int index = i_entry_to[x_dist_dub];
						stalled_to_rest_table_[x_dist_dub][index] = motor;	
						i_entry_to[x_dist_dub]++;
						n_stalled_to_teth_rest_[x_dist_dub]++;
					}
					// otherwise, add to from_rest list
					else if(delta_x == -1 * dx_rest){
						int index = i_entry_from[x_dist_dub];
						stalled_from_rest_table_[x_dist_dub][index] = motor;
						i_entry_from[x_dist_dub]++;
						n_stalled_from_teth_rest_[x_dist_dub]++;
					}
				}
				else if(motor->AtCutoff() == false){
					if(mt->lattice_[i_front + delta_x].occupied_ == true){
						int x_dist_dub = motor->x_dist_doubled_;
						// if MT's dx is towards rest, add to to_rest list
						if(delta_x == dx_rest){
							int index = i_entry_to[x_dist_dub];
							stalled_to_rest_table_[x_dist_dub][index] 
								= motor;	
							i_entry_to[x_dist_dub]++;
							n_stalled_to_teth_rest_[x_dist_dub]++;
						}
						// otherwise, add to from_rest list
						else if(delta_x == -1 * dx_rest){
							int index = i_entry_from[x_dist_dub];
							stalled_from_rest_table_[x_dist_dub][index] 
								= motor;
							i_entry_from[x_dist_dub]++;
							n_stalled_from_teth_rest_[x_dist_dub]++;
						}
						else{
							printf("hmmm??? teth tables\n");
							exit(1);
						}
					}
				}
			}
		}
	}
}

void KinesinManagement::GenerateDiffusionList(){

	int n_events = 0; 
	// Untethered statistics
	UpdateBoundIIList();
	int n_fwd_unteth = GetNumToStepForward();
	int n_bck_unteth = GetNumToStepBackward();
	while(n_fwd_unteth + n_bck_unteth > n_bound_ii_){
		double ran = properties_->gsl.GetRanProb();
		if(ran < 0.5
		&& n_fwd_unteth > 0)
			n_fwd_unteth--;
		else if(n_bck_unteth > 0)
			n_bck_unteth--;
	}
	n_events += n_fwd_unteth;
	n_events += n_bck_unteth;
	// Tethered statistics
	UpdateBoundIITetheredTable();
	int n_toward_rest[2*dist_cutoff_ + 1];
	int n_from_rest[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		n_toward_rest[x_dist_dub] = GetNumToStepToTethRest(x_dist_dub);
		n_from_rest[x_dist_dub] = GetNumToStepFromTethRest(x_dist_dub);
		int n_to = n_toward_rest[x_dist_dub];
		int n_from = n_from_rest[x_dist_dub];
		int n_tethered = n_bound_ii_tethered_[x_dist_dub];
		while(n_to + n_from > n_tethered){
			double ran = properties_->gsl.GetRanProb();
			double p_to = p_diffuse_to_teth_rest_[x_dist_dub];
			double p_from = p_diffuse_from_teth_rest_[x_dist_dub];
			double tot_prob = p_to + p_from;  // XXX is this correct? XXX
			if(ran < p_to/tot_prob
			&& n_to > 0){
				n_toward_rest[x_dist_dub]--;
				n_to = n_toward_rest[x_dist_dub];
			}
			else if(n_from > 0){
				n_from_rest[x_dist_dub]--;
				n_to = n_from_rest[x_dist_dub];
			}
		}
		n_events += n_toward_rest[x_dist_dub];
		n_events += n_from_rest[x_dist_dub];
	}
	// Only generate list if we have more than zero diffusion events
	if(n_events > 0){
		int pre_list[n_events];
		int diff_index = 0;
		for(int i_event = 0; i_event < n_fwd_unteth; i_event++){
			pre_list[diff_index] = 10;
			diff_index++;
		}
		for(int i_event = 0; i_event < n_bck_unteth; i_event++){
			pre_list[diff_index] = 20;
			diff_index++;
		}
		for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
			int n_step_to = n_toward_rest[x_dist_dub];
			for(int i_event = 0; i_event < n_step_to; i_event++){
				pre_list[diff_index] = 300 + x_dist_dub;
				diff_index++;
			}
			int n_step_from = n_from_rest[x_dist_dub];
			for(int i_event = 0; i_event < n_step_from; i_event++){
				pre_list[diff_index] = 400 + x_dist_dub;
				diff_index++;
			}
		}
		RandomNumberManagement *gsl = &properties_->gsl;
		gsl_ran_shuffle(gsl->rng, pre_list, n_events, sizeof(int));
		diffusion_list_.resize(n_events);
		for(int i_event = 0; i_event < n_events; i_event++){
			diffusion_list_[i_event] = pre_list[i_event];
		}
	}
	else
		diffusion_list_.clear();
}

int KinesinManagement::GetNumToStepForward(){

	int n_bound = n_bound_ii_;
	double p_step = p_diffuse_forward_;
	if(n_bound > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToStepBackward(){

	int n_bound = n_bound_ii_;
	double p_step = p_diffuse_backward_;
	if(n_bound > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToStepToTethRest(int x_dist_doubled){

	int n_bound = n_bound_ii_tethered_[x_dist_doubled];
	double p_step = p_diffuse_to_teth_rest_[x_dist_doubled];
	if(n_bound > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToStepFromTethRest(int x_dist_doubled){

	int n_bound = n_bound_ii_tethered_[x_dist_doubled];
	double p_step = p_diffuse_from_teth_rest_[x_dist_doubled];
	if(n_bound > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_bound);
		return n_to_step;
	}
	else{
		return 0;
	}
}

void KinesinManagement::RunDiffusion(){

//	printf("Start of kinesin diffusion cycle\n");
	GenerateDiffusionList();
	if(diffusion_list_.empty() == false){
		int n_events = diffusion_list_.size();
//		printf("%i KINESIN DIFFUSION EVENTS\n", n_events);
		//printf("%i DIFFUSION EVENTS\n", n_events);
		int x_dist_dub;
		for(int i_event = 0; i_event < n_events; i_event++){
			int diff_event = diffusion_list_[i_event];
			if(diff_event >= 300 && diff_event < 400){
				x_dist_dub = diff_event	% 100;
				diff_event = 30;
			}
			if(diff_event >= 400 && diff_event < 500){
				x_dist_dub = diff_event % 100;
				diff_event = 40;
			}
//			properties_->wallace.PrintMicrotubules(0.000);
			switch(diff_event){
				case 10:
//					printf("unteth step fwd\n");
					RunDiffusion_Forward();
					break;	
				case 20:
//					printf("unteth step bck\n");
					RunDiffusion_Backward();
					break;
				case 30:
//					printf("teth step to (%i)[%i avail]\n", x_dist_dub, 
//							n_bound_ii_tethered_[x_dist_dub]);
					RunDiffusion_To_Teth_Rest(x_dist_dub);
					break;
				case 40:
//					printf("teth step from (%i)[%i avail]\n", x_dist_dub, 
//							n_bound_ii_tethered_[x_dist_dub]);
					RunDiffusion_From_Teth_Rest(x_dist_dub);
					break;
			}
		}
	}
}

void KinesinManagement::RunDiffusion_Forward(){

	UpdateBoundIIList();
	int n_bound = n_bound_ii_;
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_ii_list_[i_entry];
		Microtubule *mt = motor->mt_;
		int i_front = motor->front_site_->index_; 
		int dx = mt->delta_x_;
		int i_plus = mt->plus_end_;
		if(i_front != i_plus){
			if(mt->lattice_[i_front + dx].occupied_ == false){
				// Get new sites
				Tubulin *new_front = &mt->lattice_[i_front + dx];
				Tubulin *new_rear = motor->front_site_;
				Tubulin *old_rear = motor->rear_site_;
				// Update new site
				new_front->motor_ = motor;
				new_front->occupied_ = true;
				// Update old site
				old_rear->motor_ = nullptr;
				old_rear->occupied_ = false;
				// Update motor
				motor->front_site_ = new_front;
				motor->rear_site_ = new_rear;

			}
			else{
	//			printf("oh well fwd\n");
			}
		}
		else{
	//		printf("cant diffuse outta this one brotha\n");
		}
	}
	else{
		printf("ya blew it. we failed to step untethered motor fwd\n");
		exit(1);
	}

}

void KinesinManagement::RunDiffusion_Backward(){

	UpdateBoundIIList();
	int n_bound = n_bound_ii_;
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_ii_list_[i_entry];
		Microtubule *mt = motor->mt_;
		int i_rear = motor->rear_site_->index_; 
		int dx = mt->delta_x_;
		int i_minus = mt->minus_end_;
		if(i_rear != i_minus){
			if(mt->lattice_[i_rear - dx].occupied_ == false){
				// Get new sites
				Tubulin *new_rear = &mt->lattice_[i_rear - dx];
				Tubulin *new_front = motor->rear_site_;
				Tubulin *old_front = motor->front_site_;
				// Update new site
				new_rear->motor_ = motor;
				new_rear->occupied_ = true;
				// Update old site
				old_front->motor_ = nullptr;
				old_front->occupied_ = false;
				// Update motor
				motor->front_site_ = new_front;
				motor->rear_site_ = new_rear;

			}
			else{
	//			printf("oh well bck\n");
			}
		}
		else{
	//		printf("cant diffuse outta this one brotha\n");
		}

	}
	else{
		printf("ya blew it. we failed to step untethered motor bck\n");
		exit(1);
	}
}

void KinesinManagement::RunDiffusion_To_Teth_Rest(int x_dist_doubled){

	int mt_length = parameters_->microtubules.length;
	int mt_array_length = mt_length - 1;	
	UpdateBoundIITetheredTable();
	int n_bound = n_bound_ii_tethered_[x_dist_doubled];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_ii_tethered_table_[x_dist_doubled][i_entry];
		Microtubule *mt = motor->mt_;
		Tubulin *near_site = motor->GetSiteCloserToRest();
		int dx = motor->GetDirectionTowardRest();
		int i_near = near_site->index_;
		if(!(i_near == mt_array_length && dx == 1)
		&& !(i_near == 0 && dx == -1)){
			if(mt->lattice_[i_near + dx].occupied_ == false){
				Tubulin *new_front, *new_rear, 
						*old_front, *old_rear;
				if(near_site == motor->front_site_){
					old_front = near_site;
					new_front = &mt->lattice_[i_near + dx];
					old_rear = motor->rear_site_;
					new_rear = motor->front_site_;
					// Update new site
					new_front->motor_ = motor;
					new_front->occupied_ = true;
					// Update old site
					old_rear->motor_ = nullptr;
					old_rear->occupied_ = false;
				}
				else if(near_site == motor->rear_site_){
					old_rear = near_site;
					new_rear = &mt->lattice_[i_near + dx];
					old_front = motor->front_site_;
					new_front = motor->rear_site_;
					// Update new site
					new_rear->motor_ = motor;
					new_rear->occupied_ = true;
					// Update old site
					old_front->motor_ = nullptr;
					old_front->occupied_ = false;
				}
				else{
					printf("woah woah woah. why is she the WRONG site?");
					printf(" (motors diffuse toward tether)\n");
					exit(1);
				}
				int x_dub_pre = motor->x_dist_doubled_;
				// Update motor
				motor->front_site_ = new_front;
				motor->rear_site_ = new_rear;
				motor->UpdateExtension();
				// Update statistics
				// Make sure an untether event wasn't forced
				if(motor->tethered_ == true){
					int x_dub_post = motor->x_dist_doubled_;
					n_bound_ii_tethered_[x_dub_pre]--;
					n_bound_ii_tethered_[x_dub_post]++;
					// Update prc1 site statistics
					AssociatedProtein *xlink = motor->xlink_;
					AssociatedProteinManagement *prc1 = &properties_->prc1;
					if(xlink->heads_active_ == 1){
						prc1->n_sites_i_tethered_[x_dub_pre]--;
						prc1->n_sites_i_tethered_[x_dub_post]++;
					}
					else if(xlink->heads_active_ == 2){
						int x_dist = xlink->x_dist_;
						prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
						prc1->n_sites_ii_tethered_[x_dub_post][x_dist] += 2;
					}
					else{
						printf("wat in diff_teth_to\n");
						exit(1);
					}
				}
			}
			else{
	//			printf("aww darn teth toward (%i ext)\n", x_dist_doubled);
			}
		}	
	}
	else{
		printf("ya blew it. we failed to step tethered motor towards\n");
		exit(1);
	}
}

void KinesinManagement::RunDiffusion_From_Teth_Rest(int x_dist_doubled){

	int mt_length = parameters_->microtubules.length;
	int mt_array_length = mt_length - 1;	
	UpdateBoundIITetheredTable();
	int n_bound = n_bound_ii_tethered_[x_dist_doubled];
	if(n_bound > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound);
		Kinesin *motor = bound_ii_tethered_table_[x_dist_doubled][i_entry];
		Microtubule *mt = motor->mt_;
		Tubulin *far_site = motor->GetSiteFartherFromRest();
		int dx = motor->GetDirectionTowardRest();
		int i_far = far_site->index_;
//		printf("i: %i, dx: %i\n", i_far, dx);
		if(!(i_far == mt_array_length && dx == -1)
		&& !(i_far == 0 && dx == 1)){
			if(mt->lattice_[i_far - dx].occupied_ == false){
				Tubulin *new_front, *new_rear, 
						*old_front, *old_rear;
				if(far_site == motor->front_site_){
					old_front = far_site;
					new_front = &mt->lattice_[i_far - dx];
					old_rear = motor->rear_site_;
					new_rear = motor->front_site_;
					// Update new site
					new_front->motor_ = motor;
					new_front->occupied_ = true;
					// Update old site
					old_rear->motor_ = nullptr;
					old_rear->occupied_ = false;
				}
				else if(far_site == motor->rear_site_){
					old_rear = far_site;
					new_rear = &mt->lattice_[i_far - dx];
					old_front = motor->front_site_;
					new_front = motor->rear_site_;
					// Update new site  
					new_rear->motor_ = motor;
					new_rear->occupied_ = true;
					// Update old site
					old_front->motor_ = nullptr;
					old_front->occupied_ = false; 
				}
				else{
					printf("woah woah woah. why is she the WRONG site?");
					printf(" (motors diffuse away from tether)\n");
					exit(1);
				}
				int x_dub_pre = motor->x_dist_doubled_;
				// Update motor
				motor->front_site_ = new_front;
				motor->rear_site_ = new_rear;
				motor->UpdateExtension();
				// Update statistics
				// Make sure an untether event wasn't forced
				if(motor->tethered_ == true){
					int x_dub_post = motor->x_dist_doubled_;
					n_bound_ii_tethered_[x_dub_pre]--;
					n_bound_ii_tethered_[x_dub_post]++;
					// Update prc1 site statistics
					AssociatedProtein *xlink = motor->xlink_;
					AssociatedProteinManagement *prc1 = &properties_->prc1;
					if(xlink->heads_active_ == 1){
						prc1->n_sites_i_tethered_[x_dub_pre]--;
						prc1->n_sites_i_tethered_[x_dub_post]++;
					}
					else if(xlink->heads_active_ == 2){
						int x_dist = xlink->x_dist_;
						prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
						prc1->n_sites_ii_tethered_[x_dub_post][x_dist] += 2;
					}
					else{
						printf("wat in diff_teth_from\n");
						exit(1);
					}
				}
			}
			else{
	//			printf("aw darn teth away from (%i ext)\n", x_dist_doubled);
			}
		}
	}
	else{
		printf("ya blew it. we failed to step tethered motor away from\n");
		exit(1);
	}
}

void KinesinManagement::GenerateKMCList(){

	/* see xlinks analog for a much simpler case and explanation */
	int n_events = 0;
	int n_bind_i_free = GetNumToBind_I();
	n_events += n_bind_i_free;
	UpdateFreeTetheredList();
	int n_bind_i_teth = GetNumToBind_I_Tethered();
	n_events += n_bind_i_teth;
	UpdateBoundIBindableList();
	int n_bind_ii = GetNumToBind_II();
	n_events += n_bind_ii;
	UpdateStepableList();
	int n_unbind_stepable = GetNumToUnbind_II_Stepable();
	n_events += n_unbind_stepable;
	int n_step = GetNumToStep(); 
	n_events += n_step;
	while(n_unbind_stepable + n_step > n_stepable_){
		if(n_step > 0){
			n_step--;
			n_events--;
		}
		else if(n_unbind_stepable > 0){
			n_unbind_stepable--;
			n_events--;
		}
	}
	UpdateStalledList();
	int n_unbind_stalled = GetNumToUnbind_II_Stalled();
	n_events += n_unbind_stalled;
	int n_failstep_ut = GetNumToFailstep();
	n_events += n_failstep_ut;
	while(n_unbind_stalled + n_failstep_ut > n_stalled_){
		if(n_failstep_ut > 0){
			n_failstep_ut--;
			n_events--;
		}
		else if(n_unbind_stalled > 0){
			n_unbind_stalled--;
			n_events--;
		}
	}
	UpdateBoundIITetheredList();
	UpdateBoundIList();
	int n_unbind_i = GetNumToUnbind_I();
	n_events += n_unbind_i;
	while(n_unbind_i + n_bind_ii > n_bound_i_){
		double p_bind = p_bind_ii_; 
		double p_unbind = p_unbind_i_;
		double p_tot = p_bind + p_unbind;
		double ran = properties_->gsl.GetRanProb();
		if(n_bind_ii > 0
		&& ran < p_bind / p_tot){
			n_bind_ii--;
			n_events--;
		}
		else if(n_unbind_i > 0){
			n_unbind_i--;
			n_events--;
		}
		else if(n_bind_ii > 0){
			n_bind_ii--;
			n_events--;
		}
	}
	int n_tether_free = GetNumToTether_Free();
	n_events += n_tether_free;
	UpdateBoundIIList();
	int n_tether_bound = GetNumToTether_Bound();
	n_events += n_tether_bound;
	int n_untether_free = GetNumToUntether_Free();
	n_events += n_untether_free;
	while(n_untether_free + n_bind_i_teth > n_free_tethered_){
			double p_unteth = p_untether_free_;
			double p_bind = p_bind_i_tethered_;
			double p_tot = p_unteth + p_bind; 
			double ran = properties_->gsl.GetRanProb();
			if(n_bind_i_teth > 0
			&& ran < p_unteth / p_tot){
				n_bind_i_teth--;
				n_events--;
			}
			else if(n_untether_free > 0){
				n_untether_free--;
				n_events--;
			}
			else if(n_bind_i_teth > 0){
				n_bind_i_teth--;
				n_events--;
			}
	}
	// Handle the statistics of differnt tether extensions separately 
	UpdateStepableTetheredTables();
	UpdateBoundIITetheredTable();
	int n_untether_bound[2*dist_cutoff_ + 1];
	int n_step_to_rest[2*dist_cutoff_ + 1];
	int n_step_from_rest[2*dist_cutoff_ + 1];
	int n_failstep_to_rest[2*dist_cutoff_ + 1];
	int n_failstep_from_rest[2*dist_cutoff_ + 1];
	for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
		n_untether_bound[x_dist_dub] = GetNumToUntether_Bound(x_dist_dub);
		n_step_to_rest[x_dist_dub] = GetNumToStep_ToTethRest(x_dist_dub);
		n_step_from_rest[x_dist_dub] = GetNumToStep_FromTethRest(x_dist_dub);
		n_failstep_to_rest[x_dist_dub] = 
				GetNumToFailstep_ToTethRest(x_dist_dub);
		n_failstep_from_rest[x_dist_dub] =
				GetNumToFailstep_FromTethRest(x_dist_dub);
		int n_unteth = n_untether_bound[x_dist_dub];
		int n_step_to = n_step_to_rest[x_dist_dub];
		int n_step_from = n_step_from_rest[x_dist_dub];
		int n_failstep_to = n_failstep_to_rest[x_dist_dub];
		int n_failstep_from = n_failstep_from_rest[x_dist_dub];
		int n_tethered = n_bound_ii_tethered_[x_dist_dub];
		// Prevent too many steps from occuring (imagine if 2 motors
		// are bound and we roll for 1 unbind but 2 steps)
		while(n_unteth + n_step_to + n_step_from > n_tethered){
			double p_to = p_step_to_teth_rest_[x_dist_dub];
			double p_from = p_step_from_teth_rest_[x_dist_dub];
			double p_tot = p_to + p_from;
			double ran = properties_->gsl.GetRanProb();
			if(ran < p_to / p_tot
			&& n_step_to > 0){
				n_step_to_rest[x_dist_dub]--;
				n_step_to = n_step_to_rest[x_dist_dub];
			}
			else if(n_step_from > 0){
				n_step_from_rest[x_dist_dub]--;
				n_step_from = n_step_from_rest[x_dist_dub];
			}
			else if(n_step_to > 0){
				n_step_to_rest[x_dist_dub]--;
				n_step_to = n_step_to_rest[x_dist_dub];
			}
			else if (n_unteth > 0){
				n_untether_bound[x_dist_dub]--;
				n_unteth = n_untether_bound[x_dist_dub];
			}
		} 
		while(n_unteth + n_failstep_to + n_failstep_from > n_tethered){
			double p_to = p_failstep_to_teth_rest_[x_dist_dub];
			double p_from = p_failstep_from_teth_rest_[x_dist_dub];
			double p_tot = p_to + p_from;
			double ran = properties_->gsl.GetRanProb();
			if(ran < p_to / p_tot
			&& n_failstep_to > 0){
				n_failstep_to_rest[x_dist_dub]--;
				n_failstep_to = n_failstep_to_rest[x_dist_dub];
			}
			else if(n_failstep_from > 0){
				n_failstep_from_rest[x_dist_dub]--;
				n_failstep_from = n_failstep_from_rest[x_dist_dub];
			}
			else if(n_failstep_to > 0){
				n_failstep_to_rest[x_dist_dub]--;
				n_failstep_to = n_failstep_to_rest[x_dist_dub];
			}
			else if (n_unteth > 0){
				n_untether_bound[x_dist_dub]--;
				n_unteth = n_untether_bound[x_dist_dub];
			}
		} 
		n_events += n_untether_bound[x_dist_dub];
		n_events += n_step_to_rest[x_dist_dub];
		n_events += n_step_from_rest[x_dist_dub];
		n_events += n_failstep_to_rest[x_dist_dub];
		n_events += n_failstep_from_rest[x_dist_dub];
	}
    if(n_events > 0){
//		printf("n_events: %i \n", n_events);
        int pre_list[n_events];
		int kmc_index = 0;
		for(int i_event = 0; i_event < n_bind_i_free; i_event++){
//			printf("10\n");
			pre_list[kmc_index] = 10;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_bind_i_teth; i_event++){
//			printf("11\n");
			pre_list[kmc_index] = 11;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_bind_ii; i_event++){
//			printf("12\n");
			pre_list[kmc_index] = 12;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_i; i_event++){
//			printf("20\n");
			pre_list[kmc_index] = 20;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_stepable; i_event++){
//			printf("21\n");
			pre_list[kmc_index] = 21;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_unbind_stalled; i_event++){
//			printf("22\n");
			pre_list[kmc_index] = 22;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_tether_free; i_event++){
//			printf("30\n");
			pre_list[kmc_index] = 30;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_tether_bound; i_event++){
//			printf("31\n");
			pre_list[kmc_index] = 31;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_untether_free; i_event++){
//			printf("50\n");
			pre_list[kmc_index] = 50;
			kmc_index++;
		}

		for(int i_event = 0; i_event < n_step; i_event++){
//			printf("60\n");
			pre_list[kmc_index] = 60;
			kmc_index++;
		}
		for(int i_event = 0; i_event < n_failstep_ut; i_event++){
//			printf("61\n");
			pre_list[kmc_index] = 61;
			kmc_index++;
		}
		for(int x_dist_dub = 0; x_dist_dub <= 2*dist_cutoff_; x_dist_dub++){
			int n_unteth_b = n_untether_bound[x_dist_dub];
			for(int i_event = 0; i_event < n_unteth_b; i_event++){
				pre_list[kmc_index] = 500 + x_dist_dub; 	
//				printf("%i\n", 500 + x_dist_dub);
				kmc_index++;
			}
			int n_step_to = n_step_to_rest[x_dist_dub];
			for(int i_event = 0; i_event < n_step_to; i_event++){
				pre_list[kmc_index] = 600 + x_dist_dub;
//				printf("%i\n", 600 + x_dist_dub);
				kmc_index++;
			}
			int n_step_from = n_step_from_rest[x_dist_dub];
			for(int i_event = 0; i_event < n_step_from; i_event++){
				pre_list[kmc_index] = 700 + x_dist_dub;
//				printf("%i\n", 700 + x_dist_dub);
				kmc_index++;
			}
			int n_failstep_to = n_failstep_to_rest[x_dist_dub];
			for(int i_event = 0; i_event < n_failstep_to; i_event++){
				pre_list[kmc_index] = 800 + x_dist_dub;
//				printf("%i\n", 800 + x_dist_dub);
				kmc_index++;
			}
			int n_failstep_from = n_failstep_from_rest[x_dist_dub]; 
			for(int i_event = 0; i_event < n_failstep_from; i_event++){
				pre_list[kmc_index] = 900 + x_dist_dub;
//				printf("%i\n", 900 + x_dist_dub);
				kmc_index++;
			}
		}
		
		RandomNumberManagement *gsl = &properties_->gsl;
        gsl_ran_shuffle(gsl->rng, pre_list, n_events, sizeof(int));
        kmc_list_.resize(n_events);
        for(int i_event = 0; i_event < n_events; i_event++){
            kmc_list_[i_event] = pre_list[i_event];
        }
    }
    else{
        kmc_list_.clear();
    }
}

int KinesinManagement::GetNumToBind_I(){

    properties_->microtubules.UpdateUnoccupiedList();
    int n_unocc = properties_->microtubules.n_unoccupied_;
	double p_bind = p_bind_i_; 
//	double p_avg = p_bind * n_unocc;
//	int n_to_bind = properties_->gsl.SamplePoissonDist(p_avg);
	if(n_unocc > 0){
		int n_to_bind = properties_->gsl.SampleBinomialDist(p_bind, n_unocc);
		return n_to_bind;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToBind_I_Tethered(){

    double weights_summed = 0;
    // Sum over all tethered but unbound motors
    for(int i_motor = 0; i_motor < n_free_tethered_; i_motor++){
        Kinesin *motor = free_tethered_list_[i_motor];
		motor->UpdateNeighborSites();
        // Get weight of all neighbor sites
        int n_neighbs = motor->n_neighbor_sites_;
        for(int i_neighb = 0; i_neighb < n_neighbs; i_neighb++){
            Tubulin *site = motor->neighbor_sites_[i_neighb];
           	double weight = motor->GetBindingWeight(site);
            weights_summed += weight;
        }
    }
	double p_bind = p_bind_i_tethered_;
	double n_avg = p_bind * weights_summed;
	if(n_avg > 0){
		int n_to_bind = properties_->gsl.SamplePoissonDist(n_avg);
		return n_to_bind;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToBind_II(){

	double p_bind = p_bind_ii_; 
	int n_able = n_bound_i_bindable_;
	if(n_able > 0){
		int n_to_bind = properties_->gsl.SampleBinomialDist(p_bind, n_able);
		return n_to_bind;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToUnbind_I(){

	int n_bound = n_bound_i_; 
	double p_unbind = p_unbind_i_;
	if(n_bound > 0){
		int n_to_unbind = 
			properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
		return n_to_unbind;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToUnbind_II_Stepable(){

	int n_bound = n_stepable_;
	double p_unbind = p_unbind_ii_stepable_;
	if(n_bound > 0){
		int n_to_unbind = 
			properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
		return n_to_unbind;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToUnbind_II_Stalled(){

	int n_bound = n_stalled_;
	double p_unbind = p_unbind_ii_stalled_;
	if(n_bound > 0){
		int n_to_unbind = 
			properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
		return n_to_unbind;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToTether_Free(){

    int n_unteth = properties_->prc1.n_untethered_;
	// Calculate how many free motors tether within delta_t on avg
	double p_teth = p_tether_free_;
	if(n_unteth > 0){
		int n_to_tether = 
			properties_->gsl.SampleBinomialDist(p_teth, n_unteth);
		return n_to_tether;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToTether_Bound(){

	double weights_summed = 0;
	for(int i_motor = 0; i_motor < n_bound_ii_; i_motor++){
		Kinesin *motor = bound_ii_list_[i_motor];
		// Motors bound to free xlinks aren't technically tethered, but
		// already have an xlink partner, so make sure to not count them
		if(motor->tethered_ == false){
			motor->UpdateNeighborXlinks();
			int n_neighbs = motor->n_neighbor_xlinks_; 
			for(int i_neighb = 0; i_neighb < n_neighbs; i_neighb++){
				AssociatedProtein *xlink = motor->neighbor_xlinks_[i_neighb];
				double weight = motor->GetTetheringWeight(xlink);
				weights_summed += weight;
			}
		}
	}
	double n_avg = p_tether_bound_ * weights_summed;
	if(n_avg > 0){
		int n_to_teth = properties_->gsl.SamplePoissonDist(n_avg);
		return n_to_teth;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToUntether_Bound(int x_dist_doubled){

	int n_teth = n_bound_ii_tethered_[x_dist_doubled]; 
	double p_unteth = p_untether_bound_[x_dist_doubled];
	if(n_teth > 0){
		int n_to_unteth = 
			properties_->gsl.SampleBinomialDist(p_unteth, n_teth);
		return n_to_unteth;
	}
	else{
		return 0;
	}
}


int KinesinManagement::GetNumToUntether_Free(){

	int n_teth = n_free_tethered_;
	double p_unteth = p_untether_free_;
	if(n_teth > 0){
		int n_to_unteth = 
			properties_->gsl.SampleBinomialDist(p_unteth, n_teth);
		return n_to_unteth; 
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToStep(){

	int n_stepable = n_stepable_;
	double p_step = p_step_;
	if(n_stepable > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_stepable);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToFailstep(){

	int n_bound = n_stalled_;
	double p_unbind = p_failstep_;
	if(n_bound > 0){
		int n_to_unbind = 
			properties_->gsl.SampleBinomialDist(p_unbind, n_bound);
		return n_to_unbind;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToStep_ToTethRest(int x_dist_doubled){

	int n_stepable = n_stepable_to_teth_rest_[x_dist_doubled];
	double p_step = p_step_to_teth_rest_[x_dist_doubled];
	if(n_stepable > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_stepable);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToStep_FromTethRest(int x_dist_doubled){

	int n_stepable = n_stepable_from_teth_rest_[x_dist_doubled];
	double p_step = p_step_from_teth_rest_[x_dist_doubled];
	if(n_stepable > 0){
		int n_to_step = 
			properties_->gsl.SampleBinomialDist(p_step, n_stepable);
		return n_to_step;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToFailstep_ToTethRest(int x_dist_doubled){

	int n_avail = n_stalled_to_teth_rest_[x_dist_doubled];
	double p_step = p_failstep_to_teth_rest_[x_dist_doubled];
	if(n_avail > 0){
		int n_failstep = 
			properties_->gsl.SampleBinomialDist(p_step, n_avail);
		return n_failstep;
	}
	else{
		return 0;
	}
}

int KinesinManagement::GetNumToFailstep_FromTethRest(int x_dist_doubled){

	int n_avail = n_stalled_from_teth_rest_[x_dist_doubled];
	double p_step = p_failstep_from_teth_rest_[x_dist_doubled];
	if(n_avail > 0){
		int n_failstep = 
			properties_->gsl.SampleBinomialDist(p_step, n_avail);
		return n_failstep;
	}
	else{
		return 0;
	}
}

void KinesinManagement::RunKMC(){

//	printf("Start of Kinesin KMC cycle\n");
    GenerateKMCList();
    if(kmc_list_.empty() == false){
        int n_events = kmc_list_.size();
//		printf("%i MOTOR KMC EVENTS\n", n_events);
		int x_dist_doubled = 0;
        for(int i_event = 0; i_event < n_events; i_event++){
            int kmc_event = kmc_list_[i_event];
			if(kmc_event >= 500 && kmc_event < 600){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 51;
			}
			if(kmc_event >= 600 && kmc_event < 700){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 62; 
			}	
			if(kmc_event >= 700 && kmc_event < 800){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 63;
			}
			if(kmc_event >= 800 && kmc_event < 900){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 64;
			}
			if(kmc_event >= 900 && kmc_event < 1000){
				x_dist_doubled = kmc_event % 100;
				kmc_event = 65;
			}
	//		properties_->wallace.PrintMicrotubules(0.000);
            switch(kmc_event){
				case 10:
	//					printf("free motor pseudo-bound\n");
						KMC_Bind_I();
                        break;
                case 11:
	//					printf("tethered motor pseudo-bound\n");
						KMC_Bind_I_Tethered();
                        break;
				case 12:
	//					printf("pseudo motor bound\n");
						KMC_Bind_II();
						break;
				case 20:
	//					printf("pseudo-bound motor unbound\n");
						KMC_Unbind_I();
						break;
				case 21:
	//					printf("untethered stepable motor unbound\n");
						KMC_Unbind_II_Stepable();
                        break;
				case 22:
	//					printf("untethered stalled motor unbound\n");
						KMC_Unbind_II_Stalled();
						break;
                case 30:
	//					printf("free motor tethered\n");
						KMC_Tether_Free();
                        break;
				case 31:
	//					printf("bound motor tethered\n");
						KMC_Tether_Bound();
						break;
				case 50:
	//					printf("free motor untethered\n");
						KMC_Untether_Free();
						break;
				case 51:
	//					printf("bound motor (ext %i) untethered\n", 
	//							x_dist_doubled);
						KMC_Untether_Bound(x_dist_doubled);
						break;
				case 60:
	//					printf("untethered motor stepped\n");
						KMC_Step();
						break;
						
				case 61:
	//					printf("motor failstepped\n");
						KMC_Failstep();
						break;
						
				case 62:
	//					printf("tethered motor (ext %i) stepped\n", 
	//							x_dist_doubled);
						KMC_Step_ToTethRest(x_dist_doubled);
						break;
				case 63:
	//					printf("tethered motor (ext %i) stepped\n", 
	//							x_dist_doubled);
						KMC_Step_FromTethRest(x_dist_doubled);
						break;
				case 64:
	//					printf("tethered motor (ext %i) failstepped\n", 
	//							x_dist_doubled);
						KMC_Failstep_ToTethRest(x_dist_doubled);
						break;
				case 65:
	//					printf("tethered motor (ext %i) failstepped\n", 
	//							x_dist_doubled);
						KMC_Failstep_FromTethRest(x_dist_doubled);
						break;
            }
        }
    }
}

void KinesinManagement::KMC_Bind_I(){

    // Make sure that at least one unbound motor exists
	properties_->microtubules.UpdateUnoccupiedList();
	if(properties_->microtubules.n_unoccupied_ > 0){	
        Kinesin *motor = GetFreeMotor();
		MicrotubuleManagement *mts = &properties_->microtubules;
		Tubulin *site = mts->GetUnoccupiedSite();
		Microtubule *mt = site->mt_;
		// Update site details
		site->motor_ = motor;
		site->occupied_ = true;
		// Update motor details
		motor->mt_ = mt;
		motor->front_site_ = site;
		motor->heads_active_ = 1;
		// Update statistics
		n_bound_i_++;
	}
	else{
		printf("Error in Bind_I: no unoccupied sites.\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Bind_I_Tethered(){
	
	// Make sure at least one unoccupied site exists
	properties_->microtubules.UpdateUnoccupiedList();
	UpdateFreeTetheredList();
	if(properties_->microtubules.n_unoccupied_ > 0
	&& n_free_tethered_ > 0){
		// Pick a random tethered free motor to bind
		int i_motor = properties_->gsl.GetRanInt(n_free_tethered_);
		Kinesin *motor = free_tethered_list_[i_motor];
		Tubulin *site = motor->GetWeightedNeighborSite();
		int attempts = 0;
		while(site == nullptr){
			if(attempts > 10*n_free_tethered_){
				break;
			}
			i_motor = properties_->gsl.GetRanInt(n_free_tethered_);
			motor = free_tethered_list_[i_motor];
			site = motor->GetWeightedNeighborSite();
			attempts++;	
		}
		if(site != nullptr){
			Microtubule *mt = site->mt_; 
			// Update site details
			site->motor_ = motor;
			site->occupied_ = true;
			// Update motor detail
			motor->mt_ = mt;
			motor->front_site_ = site;
			motor->heads_active_ = 1;
			// Update statistics;
			n_free_tethered_--;
			n_bound_i_++; 
			// update sites for prc1 management
			motor->UpdateExtension();
			if(motor->tethered_ == true){
				if(motor->xlink_->heads_active_ > 0){
					AssociatedProtein *xlink = motor->xlink_;
					AssociatedProteinManagement *prc1 = &properties_->prc1;
					int x_dist_dub = motor->x_dist_doubled_;
					if(xlink->heads_active_ == 1){
						prc1->n_sites_i_tethered_[x_dist_dub]++;	
						prc1->n_sites_i_untethered_--;
					}
					else if(xlink->heads_active_ == 2){
						int x_dist = xlink->x_dist_;
						prc1->n_sites_ii_tethered_[x_dist_dub][x_dist] += 2;
						prc1->n_sites_ii_untethered_[x_dist] -= 2;
					}
					else{
						printf("GOD DAMNIT bind_i_teth motor\n");
						exit(1);
				}
			}
		}
			else{
				printf("NOPE in bind_i_teth\n");
				exit(1);
			}
		}
		else{
//			printf("Failed to bind_i_teth\n");
//			properties_->wallace.PrintMicrotubules(2);
		}
	}
	else if(properties_->microtubules.n_unoccupied_ == 0){
		printf("Error in Bind_I_Tethered: no unoccupied sites. ");
//		properties_->wallace.PrintMicrotubules(0.5);
//		exit(1);
	}
	else{
		printf("Error in Bind_I_Tethered: no tethered free motors ");
//		properties_->wallace.PrintMicrotubules(0.5);
	}
}

void KinesinManagement::KMC_Bind_II(){

	UpdateBoundIBindableList();
	if(n_bound_i_bindable_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound_i_bindable_);
		Kinesin *motor = bound_i_bindable_list_[i_entry];
		Microtubule *mt = motor->mt_;
		Tubulin *bound_site = motor->GetActiveHeadSite();
		int i_site = bound_site->index_;
		int dx = mt->delta_x_;
		int i_plus = mt->plus_end_;
		int i_minus = mt->minus_end_;
		Tubulin *front_site,
				*rear_site; 
		if(motor->heads_active_ != 1){
			printf("nope in pseudos\n");
			exit(1);
		}
		// Choose where to bind second head
		if(i_site == i_plus){
			if(mt->lattice_[i_site - dx].occupied_ == false){
				front_site = bound_site;
				rear_site = &mt->lattice_[i_site - dx];	
			}
			else{
				printf("something wrong in bind i list\n");
				exit(1);
			}
		}
		else if(i_site == i_minus){
			if(mt->lattice_[i_site + dx].occupied_ == false){
				front_site = &mt->lattice_[i_site + dx];
				rear_site = bound_site;
			}
			else{
				printf("something wrong in bind i list\n");
				exit(1);
			}
		}
		else if(mt->lattice_[i_site + dx].occupied_ == false
			 && mt->lattice_[i_site - dx].occupied_ == false){
			double ran = properties_->gsl.GetRanProb();
			if(ran < 0.5){
				front_site = &mt->lattice_[i_site + dx];
				rear_site = bound_site;
			}
			else{
				front_site = bound_site;
				rear_site = &mt->lattice_[i_site - dx];
			}
		}
		else if(mt->lattice_[i_site + dx].occupied_ == false){
			front_site = &mt->lattice_[i_site + dx];
			rear_site = bound_site;
		}
		else if(mt->lattice_[i_site - dx].occupied_ == false){
			front_site = bound_site;
			rear_site = &mt->lattice_[i_site - dx];
		}
		else{
			printf("Error with bind i bound motors list \n");
			exit(1); 
		}
		// Update site details
		front_site->motor_ = motor;
		front_site->occupied_ = true;
		rear_site->motor_ = motor;
		rear_site->occupied_ = true;
		// Update motor details
		motor->front_site_ = front_site;
		motor->rear_site_ = rear_site;
		motor->heads_active_ = 2;	
		// Update statistics
		n_bound_i_--;
		n_bound_i_bindable_--;
		// If motor is tethered, we gotta do a whole lot of bullshit
		if(motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				AssociatedProtein* xlink = motor->xlink_;
				int x_dub_pre = motor->x_dist_doubled_;
				motor->UpdateExtension();
				// weird exception if we triggered a force untether
				if(motor->tethered_ == false){
					// Cancel out the subtraction in the force untether since
					// this particular xlink never contributed to stats, 
					// the subtraction just fucks with stuff
					n_bound_ii_tethered_tot_++;
					n_bound_ii_tethered_[x_dub_pre]++;
				}
				// Otherwise, routine statistic stuff
				else{
					// KMC stuff
					int x_dist_dub = motor->x_dist_doubled_; 
					n_bound_ii_tethered_[x_dist_dub]++;
					n_bound_ii_tethered_tot_++;
					// PRC1 sites for diffusion
					AssociatedProteinManagement *prc1 = &properties_->prc1;
					if(xlink->heads_active_ == 1){
						prc1->n_sites_i_tethered_[x_dub_pre]--;
						prc1->n_sites_i_tethered_[x_dist_dub]++;
					}
					else if(xlink->heads_active_ == 2){
						int x_dist = xlink->x_dist_;
						prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
						prc1->n_sites_ii_tethered_[x_dist_dub][x_dist] += 2;
					}
				}
			}
			// Count motors tethered to free xlinks as untethered
			else{
				n_bound_ii_++;
			}
		}
		// Otherwise, simply add to n_motors that are bound and untethered
		else{
			n_bound_ii_++;
		}
	}
	else{
//		printf("Error in Bind_II: no singly-bound motors. \n");
//		exit(1);
	}
}	

void KinesinManagement::KMC_Unbind_I(){

	UpdateBoundIList();
	if(n_bound_i_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound_i_);
		Kinesin *motor = bound_i_list_[i_entry];
		// Update site details
		Tubulin *site = motor->GetActiveHeadSite();
		site->motor_ = nullptr;
		site->occupied_ = false;
		// Update motor details
		motor->heads_active_ = 0;
		motor->mt_ = nullptr;
		motor->front_site_ = nullptr;
		motor->rear_site_ = nullptr;
		// Update statistics
		if(motor->tethered_ == true){
			if(motor->xlink_->heads_active_ > 0){
				n_free_tethered_++;
				int x_dist_dub = motor->x_dist_doubled_;
				// Update sites for prc1 management
				AssociatedProtein *xlink = motor->xlink_;
				AssociatedProteinManagement *prc1 = &properties_->prc1;	
				if(xlink->heads_active_ == 1){
					prc1->n_sites_i_tethered_[x_dist_dub]--;
					prc1->n_sites_i_untethered_++;
				}
				else{
					int x_dist = xlink->x_dist_;
					prc1->n_sites_ii_tethered_[x_dist_dub][x_dist] -= 2;
					prc1->n_sites_ii_untethered_[x_dist] += 2;
				}
			}
			// If tethered to a free xlink, untether upon unbinding
			else{
				AssociatedProtein* xlink = motor->xlink_;
				motor->tethered_ = false;
				motor->xlink_ = nullptr;
				xlink->tethered_ = false;
				xlink->motor_ = nullptr;
				properties_->prc1.n_free_tethered_--;
			}
		}
		n_bound_i_--;
	}
	else{
		printf("Error in Unbind_I: no pseudo bound motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Unbind_II_Stepable(){

    // Make sure that at least one bound motor exists
	UpdateStepableList();
    if(n_stepable_ > 0){
        // Randomly pick a bound motor (that's not on the boundary)
		int i_entry = 0;
		i_entry = properties_->gsl.GetRanInt(n_stepable_);
		Kinesin *motor = stepable_list_[i_entry];
		// Roll and randomly pick a head to unbind
		double ran = properties_->gsl.GetRanProb();
		// Update site details
		if(ran < 0.5){
			motor->front_site_->motor_ = nullptr;
			motor->front_site_->occupied_ = false;
			motor->front_site_ = nullptr;
		}
		else{
			motor->rear_site_->motor_ = nullptr;
			motor->rear_site_->occupied_ = false;
			motor->rear_site_ = nullptr;
		}
		// Update motor details
		motor->heads_active_--;
		// Update statistics
		n_bound_ii_--; 
		n_bound_i_++; 
	}
	else{
		printf("Error in Unbind_Untethered MOBILE: ");
		printf("no bound untethered motors!\n");
//      exit(1);
    }
}

void KinesinManagement::KMC_Unbind_II_Stalled(){

    // Make sure that at least one bound motor exists
	UpdateStalledList();
    if(n_stalled_ > 0){
        // Randomly pick a bound motor (that's not on the boundary)
		int i_entry = 0;
		i_entry = properties_->gsl.GetRanInt(n_stalled_);
		Kinesin *motor = stalled_list_[i_entry];
		// Roll and randomly pick a head to unbind
		double ran = properties_->gsl.GetRanProb();
		// Update site details
		if(ran < 0.5){
			motor->front_site_->motor_ = nullptr;
			motor->front_site_->occupied_ = false;
			motor->front_site_ = nullptr;
		}
		else{
			motor->rear_site_->motor_ = nullptr;
			motor->rear_site_->occupied_ = false;
			motor->rear_site_ = nullptr;
		}
		// Update motor details
		motor->heads_active_--;
		// Update statistics
		n_bound_ii_--; 
		n_bound_i_++; 
	}
	else{
		printf("Error in Unbind_Untethered STALLED: ");
		printf("no bound untethered motors!\n");
//      exit(1);
    }
}

void KinesinManagement::KMC_Tether_Free(){

	// Make sure there is at least one unbound motor
	Kinesin *motor = GetFreeMotor();
	if(properties_->prc1.n_untethered_ > 0){
		AssociatedProtein *xlink = properties_->prc1.GetUntetheredXlink();
		// Update motor and xlink details
		motor->xlink_ = xlink;
		motor->tethered_ = true;
		xlink->tethered_ = true;
		xlink->motor_ = motor; 	
		// Update statistics 
		n_free_tethered_++;
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		prc1->n_untethered_--;
	}
	else{
		printf("Error in Tether_Free: no bound untethered xlinks\n");
	}
}

void KinesinManagement::KMC_Tether_Bound(){

	UpdateBoundIIList();
	if(n_bound_ii_ > 0){
		// Pick a random free tethered motor
		int i_motor = properties_->gsl.GetRanInt(n_bound_ii_);
		Kinesin *motor = bound_ii_list_[i_motor];
		AssociatedProtein* xlink = nullptr;
		if(motor->tethered_ == false){
			xlink = motor->GetWeightedNeighborXlink();
		}
		int attempts = 0; 
		while(xlink == nullptr){
			if(attempts > 10*n_bound_ii_){
				break;
			}
			i_motor = properties_->gsl.GetRanInt(n_bound_ii_);
			motor = bound_ii_list_[i_motor];
			if(motor->tethered_ == false){
				xlink = motor->GetWeightedNeighborXlink();
			}
			attempts++;
		}
		if(xlink != nullptr){
			// Update motor and xlink details
			motor->xlink_ = xlink;
			motor->tethered_ = true;
			xlink->tethered_ = true;
			xlink->motor_ = motor;
			motor->UpdateExtension();
			if(motor->tethered_ == false){
				printf("what\n");
				exit(1);
			}
			// Update local statistics
			int x_dist_dub = motor->x_dist_doubled_; 
			n_bound_ii_tethered_[x_dist_dub]++;
			n_bound_ii_tethered_tot_++;
			n_bound_ii_--;
			// Update prc1 statistics 
			AssociatedProteinManagement *prc1 = &properties_->prc1;
			prc1->n_untethered_--;
			if(xlink->heads_active_ == 1){
				prc1->n_sites_i_untethered_--;
				prc1->n_sites_i_tethered_[x_dist_dub]++;
			}
			else if(xlink->heads_active_ == 2){
				int x_dist = xlink->x_dist_;
				xlink->UpdateExtension();
				int x_dist_post = xlink->x_dist_;
				if(x_dist != x_dist_post){
					printf("error in tetherbound 2....\n");
					exit(1);
				}
				prc1->n_sites_ii_untethered_[x_dist] -= 2;
				prc1->n_sites_ii_tethered_[x_dist_dub][x_dist] += 2;
			}
			else{
				printf("error in tether_bound...\n");
				exit(1);
			}
		}
		else{
//			printf("Failed to tether bound motor\n");
		}
	}
	else{
		printf("Error in Tether_Bound: no bound untethered motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Untether_Bound(int x_dist_doubled){

	UpdateBoundIITetheredTable();
	int n_bound_tethered = n_bound_ii_tethered_[x_dist_doubled];
	if(n_bound_tethered  > 0){
		int i_entry = properties_->gsl.GetRanInt(n_bound_tethered);
		Kinesin *motor = bound_ii_tethered_table_[x_dist_doubled][i_entry];
		motor->UpdateExtension();
		int x_dub_pre = motor->x_dist_doubled_;
		if(x_dub_pre != x_dist_doubled){
			printf("error in Untether_Bound (motor): ");
			printf("%i received; %i in motor\n", x_dist_doubled, x_dub_pre);
			exit(1);
		}
		AssociatedProtein *xlink = motor->xlink_;
		// Update motor and xlink details
		xlink->motor_ = nullptr;
		xlink->tethered_ = false;
		motor->xlink_ = nullptr;
		motor->tethered_ = false; 
		// Update statistics
		motor->UpdateExtension();
		n_bound_ii_tethered_[x_dub_pre]--;
		n_bound_ii_tethered_tot_--;
		n_bound_ii_++;
		properties_->prc1.n_untethered_++;
		// Update sites for prc1_management
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		if(xlink->heads_active_ == 1){
			prc1->n_sites_i_tethered_[x_dub_pre]--;
			prc1->n_sites_i_untethered_++;
		}
		else if(xlink->heads_active_ ==2){
			int x_dist = xlink->x_dist_;
			prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
			prc1->n_sites_ii_untethered_[x_dist] += 2;
		}
		else{
			printf("god damnit in kmc_untether_bound\n");
			exit(1);
		}
	}
	else{
		printf("Error in Untether_Bound: no bound tethered motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Untether_Free(){

	UpdateFreeTetheredList();
	if(n_free_tethered_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_free_tethered_);
		Kinesin *motor = free_tethered_list_[i_entry];
		AssociatedProtein *xlink = motor->xlink_;
		// Update motor and xlink detail
		xlink->motor_ = nullptr; 
		xlink->tethered_ = false; 
		motor->xlink_ = nullptr;
		motor->tethered_ = false;
		// Update statistics
		n_free_tethered_--;
		properties_->prc1.n_untethered_++;
	}
	else{
		printf("Error in Untether_Free: no free tethered motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Step(){

    // Make sure there is at least one stepable untethered motor
	UpdateStepableList();
	if(n_stepable_ > 0){
		int i_entry = properties_->gsl.GetRanInt(n_stepable_);
		Kinesin *motor = stepable_list_[i_entry];
		Microtubule *mt = motor->mt_; 
		Tubulin *old_front_site = motor->front_site_; 
		Tubulin *old_rear_site = motor->rear_site_;
		int i_old_front = old_front_site->index_;
		int dx = mt->delta_x_;
		Tubulin *new_front_site = &mt->lattice_[i_old_front + dx];
		Tubulin *new_rear_site = old_front_site; 
		// Update front head of motor
		motor->front_site_ = new_front_site;
		new_front_site->motor_ = motor;
		new_front_site->occupied_ = true;
		// Update rear head of motor 
		motor->rear_site_ = new_rear_site;
		old_rear_site->motor_ = nullptr;
		old_rear_site->occupied_ = false;
	}
	else{
//		printf("Error in Step_Untethered: no stepable untethered motors\n");
    }
}

void KinesinManagement::KMC_Failstep(){

	/*
    // Make sure that at least one bound motor exists
	UpdateStalledList();
    if(n_stalled_ > 0){
        // Randomly pick a bound motor (that's not on the boundary)
		int i_entry = properties_->gsl.GetRanInt(n_stalled_);
		Kinesin *motor = stalled_list_[i_entry];
		// Update site details
		motor->rear_site_->motor_ = nullptr;
		motor->rear_site_->occupied_ = false;
		// Update motor details
		motor->heads_active_--;
		motor->rear_site_ = nullptr;
		// Update statistics
		n_bound_ii_--; 
		n_bound_i_++;
		
	}
	else{
		printf("Error in Failstep_UT: no stalled untethered motors!\n");
//      exit(1);
    }
	*/
}

void KinesinManagement::KMC_Step_ToTethRest(int x_dist_doubled){

	UpdateStepableTetheredTables();
	int n_stepable = n_stepable_to_teth_rest_[x_dist_doubled]; 
	if(n_stepable > 0){
		int i_entry = properties_->gsl.GetRanInt(n_stepable);
		Kinesin *motor = stepable_to_rest_table_[x_dist_doubled][i_entry];
		int x_dub_pre = motor->x_dist_doubled_;
		if(x_dub_pre != x_dist_doubled){
			printf("error in Step_Tethered (motor): ");
			printf("%i received; %i in motor\n", x_dist_doubled, x_dub_pre);
			exit(1);
		}
		Microtubule *mt = motor->mt_;
		Tubulin *old_front_site = motor->front_site_;
		Tubulin *old_rear_site = motor->rear_site_;
		int i_old_front = old_front_site->index_;
		int dx = mt->delta_x_;
		int dx_rest = motor->GetDirectionTowardRest();
		if(dx != dx_rest){
			printf("wat in step to rest motor KMC\n");
			exit(1);
		}
//		printf("coord is %g, ext is %i, dx is %i\n", 
//				motor->GetStalkCoordinate(), x_dist_doubled, dx);
		Tubulin *new_front_site = &mt->lattice_[i_old_front + dx];
		Tubulin *new_rear_site = old_front_site; 
		// Update front head of motor
		motor->front_site_ = new_front_site;
		new_front_site->motor_ = motor;
		new_front_site->occupied_ = true;
		// Update rear head of motor 
		motor->rear_site_ = new_rear_site;
		old_rear_site->motor_ = nullptr;
		old_rear_site->occupied_ = false;
		// Update statistics
		motor->UpdateExtension();
		int x_dub_post = motor->x_dist_doubled_;
		int plus_end = motor->mt_->plus_end_; 
		int i_front = motor->front_site_->index_;
		if(x_dub_pre == x_dub_post
		&& i_front != plus_end){
			printf("error in step_teth (motor)\n");
			exit(1);
		}
		n_bound_ii_tethered_[x_dub_pre]--;
		n_bound_ii_tethered_[x_dub_post]++;
		// Update site statistics for prc1
		AssociatedProtein *xlink = motor->xlink_;
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		if(xlink->heads_active_ == 1){
			prc1->n_sites_i_tethered_[x_dub_pre]--;
			prc1->n_sites_i_tethered_[x_dub_post]++;
		}
		else if(xlink->heads_active_ == 2){
			int x_dist = xlink->x_dist_;
			prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
			prc1->n_sites_ii_tethered_[x_dub_post][x_dist] += 2;
		}
		else{
			printf("nope. tethered stepe motor.\n");
			exit(1);
		}
	}
	else{
//		printf("Error in Step_Tethered: no stepable motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Step_FromTethRest(int x_dist_doubled){

	UpdateStepableTetheredTables();
	int n_stepable = n_stepable_from_teth_rest_[x_dist_doubled]; 
	if(n_stepable > 0){
		int i_entry = properties_->gsl.GetRanInt(n_stepable);
		Kinesin *motor = stepable_from_rest_table_[x_dist_doubled][i_entry];
		int x_dub_pre = motor->x_dist_doubled_;
		if(x_dub_pre != x_dist_doubled){
			printf("error in Step_Tethered (motor): ");
			printf("%i received; %i in motor\n", x_dist_doubled, x_dub_pre);
			exit(1);
		}
		Microtubule *mt = motor->mt_;
		Tubulin *old_front_site = motor->front_site_;
		Tubulin *old_rear_site = motor->rear_site_;
		int i_old_front = old_front_site->index_;
		int dx = mt->delta_x_;
		int dx_rest = motor->GetDirectionTowardRest();
		if(dx != -1 * dx_rest){
			printf("wat in step to rest motor KMC\n");
			exit(1);
		}
//		printf("coord is %g, ext is %i, dx is %i\n", 
//				motor->GetStalkCoordinate(), x_dist_doubled, dx);
		Tubulin *new_front_site = &mt->lattice_[i_old_front + dx];
		Tubulin *new_rear_site = old_front_site; 
		// Update front head of motor
		motor->front_site_ = new_front_site;
		new_front_site->motor_ = motor;
		new_front_site->occupied_ = true;
		// Update rear head of motor 
		motor->rear_site_ = new_rear_site;
		old_rear_site->motor_ = nullptr;
		old_rear_site->occupied_ = false;
		// Update statistics
		motor->UpdateExtension();
		int x_dub_post = motor->x_dist_doubled_;
		int plus_end = motor->mt_->plus_end_; 
		int i_front = motor->front_site_->index_;
		if(x_dub_pre == x_dub_post
		&& i_front != plus_end){
			printf("error in step_teth (motor)\n");
			exit(1);
		}
		n_bound_ii_tethered_[x_dub_pre]--;
		n_bound_ii_tethered_[x_dub_post]++;
		// Update site statistics for prc1
		AssociatedProtein *xlink = motor->xlink_;
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		if(xlink->heads_active_ == 1){
			prc1->n_sites_i_tethered_[x_dub_pre]--;
			prc1->n_sites_i_tethered_[x_dub_post]++;
		}
		else if(xlink->heads_active_ == 2){
			int x_dist = xlink->x_dist_;
			prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
			prc1->n_sites_ii_tethered_[x_dub_post][x_dist] += 2;
		}
		else{
			printf("nope. tethered stepe motor.\n");
			exit(1);
		}
	}
	else{
//		printf("Error in Step_Tethered: no stepable motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Failstep_ToTethRest(int x_dist_doubled){

	UpdateStalledTetheredTables();
	int n_stalled = n_stalled_to_teth_rest_[x_dist_doubled]; 
	if(n_stalled > 0){
		int i_entry = properties_->gsl.GetRanInt(n_stalled);
		Kinesin *motor = stalled_to_rest_table_[x_dist_doubled][i_entry];
		int x_dub_pre = motor->x_dist_doubled_;
		if(x_dub_pre != x_dist_doubled){
			printf("error in Step_Tethered (motor): ");
			printf("%i received; %i in motor\n", x_dist_doubled, x_dub_pre);
			exit(1);
		}
		Microtubule *mt = motor->mt_;
		int dx = mt->delta_x_;
		int dx_rest = motor->GetDirectionTowardRest();
		if(dx != dx_rest){
			printf("wat in step to rest motor KMC\n");
			exit(1);
		}
		// Update site details
		motor->rear_site_->motor_ = nullptr;
		motor->rear_site_->occupied_ = false;
		// Update motor details
		motor->heads_active_--;
		motor->rear_site_ = nullptr;
		// Update statistics
		motor->UpdateExtension();
		int x_dub_post = motor->x_dist_doubled_;
		int plus_end = motor->mt_->plus_end_; 
		int i_front = motor->front_site_->index_;
		if(x_dub_pre == x_dub_post
		&& i_front != plus_end){
			printf("error in failstep_teth (motor)\n");
			//exit(1);
		}
		n_bound_ii_tethered_[x_dub_pre]--;
		n_bound_i_++;
		// Update site statistics for prc1
		AssociatedProtein *xlink = motor->xlink_;
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		if(xlink->heads_active_ == 1){
			prc1->n_sites_i_tethered_[x_dub_pre]--;
			prc1->n_sites_i_tethered_[x_dub_post]++;
		}
		else if(xlink->heads_active_ == 2){
			int x_dist = xlink->x_dist_;
			prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
			prc1->n_sites_ii_tethered_[x_dub_post][x_dist] += 2;
		}
		else{
			printf("nope. tethered stall motor.\n");
			exit(1);
		}
	}
	else{
		printf("Error in Failstep_Tethered: no stalled motors!\n");
//		exit(1);
	}
}

void KinesinManagement::KMC_Failstep_FromTethRest(int x_dist_doubled){

	UpdateStalledTetheredTables();
	int n_stalled = n_stalled_from_teth_rest_[x_dist_doubled]; 
	if(n_stalled > 0){
		int i_entry = properties_->gsl.GetRanInt(n_stalled);
		Kinesin *motor = stalled_from_rest_table_[x_dist_doubled][i_entry];
		int x_dub_pre = motor->x_dist_doubled_;
		if(x_dub_pre != x_dist_doubled){
			printf("error in Step_Tethered (motor): ");
			printf("%i received; %i in motor\n", x_dist_doubled, x_dub_pre);
			exit(1);
		}
		Microtubule *mt = motor->mt_;
		int dx = mt->delta_x_;
		int dx_rest = motor->GetDirectionTowardRest();
		if(dx == dx_rest){
			printf("wat in step to rest motor KMC\n");
			exit(1);
		}
		// Update site details
		motor->rear_site_->motor_ = nullptr;
		motor->rear_site_->occupied_ = false;
		// Update motor details
		motor->heads_active_--;
		motor->rear_site_ = nullptr;
		// Update statistics
		motor->UpdateExtension();
		int x_dub_post = motor->x_dist_doubled_;
		int plus_end = motor->mt_->plus_end_; 
		int i_front = motor->front_site_->index_;
		if(x_dub_pre == x_dub_post
		&& i_front != plus_end){
			printf("error in failstep_teth (motor)\n");
			//exit(1);
		}
		n_bound_ii_tethered_[x_dub_pre]--;
		n_bound_i_++;
		// Update site statistics for prc1
		AssociatedProtein *xlink = motor->xlink_;
		AssociatedProteinManagement *prc1 = &properties_->prc1;
		if(xlink->heads_active_ == 1){
			prc1->n_sites_i_tethered_[x_dub_pre]--;
			prc1->n_sites_i_tethered_[x_dub_post]++;
		}
		else if(xlink->heads_active_ == 2){
			int x_dist = xlink->x_dist_;
			prc1->n_sites_ii_tethered_[x_dub_pre][x_dist] -= 2;
			prc1->n_sites_ii_tethered_[x_dub_post][x_dist] += 2;
		}
		else{
			printf("nope. tethered stall motor.\n");
			exit(1);
		}
	}
	else{
		printf("Error in Failstep_Tethered: no stalled motors!\n");
//		exit(1);
	}
}
