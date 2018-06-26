#include "master_header.h"
#include "curator.h"

Curator::Curator(){
}

void Curator::InitializeSimulation(system_parameters *parameters, 
								   system_properties *properties){

	parameters_ = parameters;
	properties_ = properties;
	start_ = clock();
	SetParameters();
	SetExperimentalStage();
	OutputSimDetails();
}

void Curator::SetParameters(){
	
	int n_steps = parameters_->n_steps;
	int n_datapoints = parameters_->n_datapoints;
	data_threshold_ = parameters_->data_threshold;
	range_of_data_ = n_steps - data_threshold_;
	n_pickup_ = range_of_data_/n_datapoints;
	equil_milestone_ = data_threshold_/10;
	data_milestone_ = range_of_data_/10;
}

void Curator::SetExperimentalStage(){

	// Initialize the general science library (gsl) class; 
	// just an easy way of sampling distributions and referencing the RNG
	properties_->gsl.Initialize(parameters_->seed);
	// Initialize microtubules, kinesin4, and prc1 classes 
	properties_->microtubules.Initialize(parameters_, properties_);
	properties_->kinesin4.Initialize(parameters_, properties_); 
	properties_->prc1.Initialize(parameters_, properties_);
}

void Curator::OutputSimDetails(){

	int n_steps = parameters_->n_steps;
	double delta_t = parameters_->delta_t;
    printf("Total simulation duration: %g seconds\n", delta_t*n_steps);
    printf("Timestep duration: %g seconds\n\n", delta_t);
    fflush(stdout);
}

void Curator::OpenFiles(char* sim_name){

	char occupancy_file[160], 
		 motor_ID_file[160], xlink_ID_file[160], 
		 tether_coord_file[160], mt_coord_file[160], 
		 motor_extension_file[160], xlink_extension_file[160],
		 motor_force_file[160], xlink_force_file[160], total_force_file[160];
	// Generate names of output files based on the input simulation name
	sprintf(occupancy_file, "%s_occupancy.file", sim_name);
	sprintf(motor_ID_file, "%s_motorID.file", sim_name);
	sprintf(xlink_ID_file, "%s_xlinkID.file", sim_name);	
	sprintf(tether_coord_file, "%s_tether_coord.file", sim_name);
	sprintf(mt_coord_file, "%s_mt_coord.file", sim_name);
	sprintf(motor_extension_file, "%s_motor_extension.file", sim_name);
	sprintf(xlink_extension_file, "%s_xlink_extension.file", sim_name);
	sprintf(motor_force_file, "%s_motor_force.file", sim_name);
	sprintf(xlink_force_file, "%s_xlink_force.file", sim_name);
	sprintf(total_force_file, "%s_total_force.file", sim_name);
	// Check to see if sim files already exist
	if (FileExists(occupancy_file)){

		printf("Simulation file with this name already exists!\n");
		printf("Do you wish to overwrite? y/n\n");
		std::string response; 
		bool response_unacceptable = true;
		while(response_unacceptable){
			std::getline(std::cin, response);
			if(response == "n"){
				printf("Simulation terminated.\n");
				exit(1); 
			}
			else if(response == "y"){
				printf("Very well.");
			   	printf(" Overwriting data for sim '%s'\n\n", sim_name);
				response_unacceptable = false; 
			}
			else{
				printf("bro I said y or n. try again plz\n");
			}
		}

	}
	// Open occupancy file, which stores the species ID of each occupant 
	// (or -1 for none) for all MT sites during data collection (DC)
	properties_->occupancy_file_ = gfopen(occupancy_file, "w");
	// Open motor ID file, which stores the unique ID of all bound motors 
	// (unbound not tracked) and their respective site indices during DC
	properties_->motor_ID_file_ = gfopen(motor_ID_file, "w");
	// Open xlink ID file, which does the same 
	// as the motor ID file but for xlinks
	properties_->xlink_ID_file_ = gfopen(xlink_ID_file, "w");
	// Open tether coord file, which stores the coordinates 
	// of the anchor points of tethered motors
	properties_->tether_coord_file_ = gfopen(tether_coord_file, "w");
	// Open mt coord file, which stores the coordinates 
	// of the left-most edge of each microtubule during DC
	properties_->mt_coord_file_ = gfopen(mt_coord_file, "w");
	// Open motor extension file, which stores the number of motors 
	// with a certain tether extension for all possible extensions
	properties_->motor_extension_file_ = gfopen(motor_extension_file, "w");
	// Open xlink extension file, which stores the number of stage-2 
	// xlinks at a certain extension for all possible extensions
	properties_->xlink_extension_file_ = gfopen(xlink_extension_file, "w");
	// Open motor force file, which stores the sum 
	// of forces coming from motor tether extensions
	properties_->motor_force_file_ = gfopen(motor_force_file, "w");
	// Open xlink force file, which stores the sum
	// of forces coming from xlink extensions
	properties_->xlink_force_file_ = gfopen(xlink_force_file, "w");
	// Open total force file, which stores the sum of ALL 
	// forces coming from xlink and motor tether extensions
	properties_->total_force_file_ = gfopen(total_force_file, "w");
}

void Curator::PrintMicrotubules(){

    int n_mts = parameters_->n_microtubules;
    int mt_length = parameters_->length_of_microtubule;
	// Figure out which MT is the farthest left 
	int leftmost_coord = 0;
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		int mt_coord = properties_->microtubules.mt_list_[i_mt].coord_;
		if(i_mt == 0)
			leftmost_coord = mt_coord;	
		else if(mt_coord < leftmost_coord)
			leftmost_coord = mt_coord;
	}
	// Print out MTs
    for(int i_mt = n_mts - 1; i_mt >= 0; i_mt--){
		Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
		int mt_coord = mt->coord_;
		int delta = mt_coord - leftmost_coord;
		if(delta < 0){
			printf("wat\n");
			exit(1);
		}
		for(int i_entry = 0; i_entry < delta; i_entry++){
			printf(" ");
		}
        for(int i_site = 0; i_site < mt_length; i_site++){
			Tubulin *site = &mt->lattice_[i_site];
            if(site->occupied_ == false)
                printf("=");
			else if(site->xlink_ != nullptr){
				AssociatedProtein *xlink = site->xlink_;
				if(xlink->heads_active_ == 1){
					if(xlink->tethered_ == false)
						printf("i");
					else
						printf("I");
				}
				else if(xlink->heads_active_ == 2)
					if(xlink->tethered_ == false){
						printf("x");			
//						printf("%i", mt->lattice_[i_site].xlink_->x_dist_);
					}	
					else
						printf("X");
				else{
					printf("no sunny. look in wallace's print\n");
					exit(1);
				}
			}
			else if(site->motor_ != nullptr){
				Kinesin *motor = site->motor_; 
				if(motor->heads_active_ == 1){
					if(motor->tethered_ == false)
						printf("m");
					else
						printf("M");
				}
				else if(motor->heads_active_ == 2){
					int i_front = motor->front_site_->index_;
					int i_rear = motor->rear_site_->index_;
					if(i_front > i_rear){
						if(i_site == i_rear){
							if(motor->tethered_ == false)
								printf("(");
							else
								printf("[");
						}
						if(i_site == i_front){
							if(motor->tethered_ == false)
								printf("%i)", motor->ID_);
							else
								printf("%i]", motor->ID_);
						}
					} 
					else if(i_front < i_rear){
						if(i_site == i_front){	
							if(motor->tethered_ == false)
								printf("(");
							else
								printf("[");
						}
						if(i_site == i_rear){
							if(motor->tethered_ == false)
								printf("%i)", motor->ID_);
							else
								printf("%i]", motor->ID_);
						}
					}
					else{
						printf("error in print MT\n");
						exit(1);
					}
				}
			}
		}
		printf(" %i\n", mt->polarity_);
	}   

/* site coordinate printout below */
	/*
	int mt1_coord = properties_->microtubules.mt_list_[0].coord_;
	int mt2_coord = properties_->microtubules.mt_list_[1].coord_;
	int greater_coord = 0;
	if(mt1_coord > mt2_coord)
		greater_coord = mt1_coord;
	else
		greater_coord = mt2_coord;
	int extra_digits = 0;
	for(int i_site = 0; i_site < mt_length + greater_coord; i_site++){
		if(extra_digits > 0)
			extra_digits--;
		else if(i_site%5 == 0){
			printf("%i", i_site);
			if(i_site < 10)
				extra_digits = 0;
			else if(i_site < 100)
				extra_digits = 1;
			else if(i_site < 1000)
				extra_digits = 2;
			else if(i_site < 10000)
				extra_digits = 3;
			else{
				printf("what the fuck are u doing bro. why do you need more than 10,000 sites??\n");
				exit(1);
			}
		}
		else if(i_site == (mt_length + greater_coord) - 1)
			printf("%i", i_site);
		else
			printf(" ");
	}
	*/
	printf("\n");
}

void Curator::PrintMicrotubules(double pause_duration){

	PrintMicrotubules();
	PauseSim(pause_duration);

}

void Curator::OutputData(){

	int n_mts = parameters_->n_microtubules;
	int mt_length = parameters_->length_of_microtubule;
	// Get file pointers from system properties
	FILE *occupancy_file = properties_->occupancy_file_;
	FILE *motor_ID_file = properties_->motor_ID_file_;
	FILE *xlink_ID_file = properties_->xlink_ID_file_;
	FILE *tether_coord_file = properties_->tether_coord_file_;
	FILE *mt_coord_file = properties_->mt_coord_file_;
	FILE *motor_extension_file = properties_->motor_extension_file_;
	FILE *xlink_extension_file = properties_->xlink_extension_file_;
	FILE *motor_force_file = properties_->motor_force_file_;
	FILE *xlink_force_file = properties_->xlink_force_file_;
	FILE *total_force_file = properties_->total_force_file_;
	// Create arrays to store data at each timestep; ptrs to write it to file 
	double mt_coord_array[n_mts];
	double *mt_coord_ptr = mt_coord_array;
	// For extension statistics, data is on a per-extension basis
	int motor_ext_cutoff = properties_->kinesin4.dist_cutoff_;
	int motor_extension_array[2*motor_ext_cutoff + 1];
	int *motor_extension_ptr = motor_extension_array; 
	int xlink_ext_cutoff = properties_->prc1.dist_cutoff_; 
	int xlink_extension_array[xlink_ext_cutoff + 1];
	int *xlink_extension_ptr = xlink_extension_array;
	// Back to normal per-MT array format
	double motor_force_array[n_mts];
	double *motor_force_ptr = motor_force_array;
	double xlink_force_array[n_mts];
	double *xlink_force_ptr = xlink_force_array;
	double total_force_array[n_mts];
	double *total_force_ptr = total_force_array;	
	// Run through all MTs and get data for each
	for(int i_mt = 0; i_mt < n_mts; i_mt++){
		Microtubule *mt = &properties_->microtubules.mt_list_[i_mt];
		// Create arrays & ptrs for intraMT data 
		int motor_ID_array[mt_length],
			xlink_ID_array[mt_length], 
			occupancy_array[mt_length];
		int	*motor_ID_ptr = motor_ID_array, 
			*xlink_ID_ptr = xlink_ID_array, 
			*occupancy_ptr = occupancy_array;
		double tether_coord_array[mt_length];
		double *teth_coord_ptr = tether_coord_array;
		// Run through all sites on this particular MT
		for(int i_site = 0; i_site < mt_length; i_site++){
			Tubulin *site = &mt->lattice_[i_site];
			// If unoccupied, store the speciesID of tubulin to occupancy file
			// and an ID of -1 (null) to motor/xlink ID files 
			if(site->occupied_ == false){
				occupancy_array[i_site] = site->speciesID_;
				motor_ID_array[i_site] = -1;
				xlink_ID_array[i_site] = -1;
				tether_coord_array[i_site] = -1;
			}
			// If occupied by xlink, store its species ID to occupancy_file,
			// its unique ID to the xlink ID file, and -1 to motor ID file
			else if(site->xlink_ != nullptr){
				occupancy_array[i_site] = site->xlink_->speciesID_;
				motor_ID_array[i_site] = -1;
				xlink_ID_array[i_site] = site->xlink_->ID_;
				tether_coord_array[i_site] = -1;
			}
			// If occupied by motor, store its species ID to occupancy_file, 
			// its unique ID to the motor ID file, and -1 to xlink ID file
			else if(site->motor_ != nullptr){
				occupancy_array[i_site] = site->motor_->speciesID_;
				motor_ID_array[i_site] = site->motor_->ID_;
				xlink_ID_array[i_site] = -1;
				if(site->motor_->tethered_ == true){
					AssociatedProtein* xlink = site->motor_->xlink_;
					double anchor_coord = xlink->GetAnchorCoordinate();
					tether_coord_array[i_site] = anchor_coord; 
				}
				else{
					tether_coord_array[i_site] = -1;
				}
			}
		}
		mt_coord_array[i_mt] = mt->coord_; 
		motor_force_array[i_mt] = mt->GetNetForce_Motors();
		xlink_force_array[i_mt] = mt->GetNetForce_Xlinks();
		total_force_array[i_mt] = mt->GetNetForce();
		// Write the data to respective files one microtubule at a time
		fwrite(occupancy_ptr, sizeof(int), mt_length, occupancy_file);
		fwrite(motor_ID_ptr, sizeof(int), mt_length, motor_ID_file);
		fwrite(xlink_ID_ptr, sizeof(int), mt_length, xlink_ID_file);
		fwrite(teth_coord_ptr, sizeof(double), mt_length, tether_coord_file); 
	}	
	// Scan through kinesin4/prc1 statistics to get extension occupancies 
	for(int i_ext = 0; i_ext <= 2*motor_ext_cutoff; i_ext++){
		KinesinManagement *kinesin4 = &properties_->kinesin4; 
		motor_extension_array[i_ext] = kinesin4->n_bound_tethered_[i_ext];
	}
	for(int i_ext = 0; i_ext <= xlink_ext_cutoff; i_ext++){
		AssociatedProteinManagement *prc1 = &properties_->prc1; 
		xlink_extension_array[i_ext] = prc1->n_double_bound_[i_ext]; 
	}
	// Write the data to respective files one timestep at a time 
	fwrite(mt_coord_ptr, sizeof(double), n_mts, mt_coord_file);
	fwrite(motor_force_ptr, sizeof(double), n_mts, motor_force_file);
	fwrite(xlink_force_ptr, sizeof(double), n_mts, xlink_force_file);
	fwrite(total_force_ptr, sizeof(double), n_mts, total_force_file);
	fwrite(motor_extension_ptr, sizeof(int), 2*motor_ext_cutoff + 1, 
			motor_extension_file);
	fwrite(xlink_extension_ptr, sizeof(int), xlink_ext_cutoff + 1, 
			xlink_extension_file);
	fwrite(motor_force_ptr, sizeof(double), n_mts, motor_force_file);
	fwrite(xlink_force_ptr, sizeof(double), n_mts, xlink_force_file);
	fwrite(total_force_ptr, sizeof(double), n_mts, total_force_file);
}

void Curator::UpdateTimestep(int i_step){

	properties_->current_step_ = i_step;
	// Give updates on equilibrium process (every 10%)
	if(i_step < data_threshold_ && i_step%equil_milestone_ == 0){
		printf("Equilibration is %i percent complete (step # %i)\n", 
		(int)(i_step/equil_milestone_)*10, i_step);
	}
	// Start data collection at appropriate step threshold
	else if(i_step >= data_threshold_){
		int delta = i_step - data_threshold_;
		// Collect data every n_pickup timesteps
		if(delta%n_pickup_ == 0){
			OutputData();
		}
		if(delta%data_milestone_ == 0){
			printf("Data collection is %i percent complete (step # %i)\n",
					(int)(delta/data_milestone_)*10, i_step);
		}
		else if(delta == range_of_data_ - 1){
			printf("Done!");
		}
	}
}

void Curator::PauseSim(double duration){
	
	// Duration should be input in seconds
	pause_dur_.tv_sec = (int) duration;	 
	pause_dur_.tv_nsec = (duration - (int)duration)*100000000;
	nanosleep(&pause_dur_, NULL);
}

void Curator::OutputSimDuration(){

	finish_ = clock();
	sim_duration_ = (double)(finish_ - start_)/CLOCKS_PER_SEC;
/*	stream_ = fopen("sim_duration.dat", "w");
	fprintf(stream_, "Time to execute sim: %f seconds.\n", sim_duration_);
	fclose(stream_);
*/	
	printf(" Time to execute: %f seconds.\n\n", sim_duration_);
}

void Curator::CleanUp(){

	fclose(properties_->occupancy_file_);
	fclose(properties_->motor_ID_file_);
	fclose(properties_->xlink_ID_file_);
	fclose(properties_->tether_coord_file_);
	fclose(properties_->mt_coord_file_);
	fclose(properties_->motor_extension_file_);
	fclose(properties_->xlink_extension_file_);
	fclose(properties_->motor_force_file_);
	fclose(properties_->xlink_force_file_);
	fclose(properties_->total_force_file_);
}

bool Curator::FileExists(std::string file_name){

	struct stat buffer;
	return (stat(file_name.c_str(), &buffer) != -1);
}