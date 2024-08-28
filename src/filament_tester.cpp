#include "cylaks/filament_tester.hpp"
#include "cylaks/protein_tester.hpp"

void FilamentTester::Initialize(ProteinTester *proteins) {

  recorded_force_.resize(Sys::n_steps_run_);
  proteins_ = proteins;
  FilamentManager::proteins_ = dynamic_cast<ProteinManager *>(proteins_);
  FilamentManager::SetParameters();
  FilamentManager::GenerateFilaments();
}

void FilamentTester::UpdateForces() {

  FilamentManager::UpdateForces();
  if (Sys::test_mode_ != "filament_forced_slide") {
    return;
  }
  // Add necessary force to get desired sliding velocity
  // (only applies to x dimension, so i = 0 index)
  double f_required{Sys::slide_velocity_ * protofilaments_[1].gamma_[0]};
  if (Sys::constant_velocity_) {
    double f_applied{f_required - protofilaments_[1].force_[0]};
    if (Sys::i_step_ < Sys::i_pause_ or Sys::i_step_ >= Sys::i_resume_) {
      protofilaments_[1].force_[0] = -f_required;
      // Record applied force
      if (recorded_force_[Sys::i_step_] == 0) {
        recorded_force_[Sys::i_step_] = -f_applied;
      }
    }
  } else {
    /* FORCE-DEP UNBINDING */
    //printf("Petets model \n");
    double f_stall{6.0}; // pN

    //Petes model
    double slide_v = Sys::slide_velocity_;
    double a = slide_v*1 +10;
    double b = slide_v*0.24 +38;
    double c = 7;
    double f_bot{protofilaments_[0].force_[0]};
    double f_per_mot_bot{-f_bot / n_motors_bot};
    double new_velocityb{a-b*f_per_mot_bot+c*f_per_mot_bot*f_per_mot_bot};
    if (new_velocityb>Sys::slide_velocity_){
      new_velocityb=Sys::slide_velocity_;
    }
    if (new_velocityb<0){
      new_velocityb=0;
    }
    //if (new_velocityb<0){
    //  new_velocityb=0;
    //}
    //printf("velocity b is %f \n",new_velocityb);
    protofilaments_[0].force_[0] = new_velocityb*protofilaments_[0].gamma_[0];
    //printf("force is %f \n", protofilaments_[0].force_[0]);
    double f_top{protofilaments_[1].force_[0]};
    double f_per_mot_top{f_top / n_motors_top};
    double new_velocityt{a-b*f_per_mot_top+c*f_per_mot_top*f_per_mot_top};
    if (new_velocityt>Sys::slide_velocity2_){
      new_velocityt=Sys::slide_velocity2_;

    }
    if (new_velocityt<0){
      new_velocityt=0;
    }
    //if (new_velocityt<0){
    //  new_velocityt=0;
    //}
    protofilaments_[1].force_[0] = -new_velocityt*protofilaments_[1].gamma_[0];
    //printf("velocity t is %f \n",new_velocityt);
  }
  //If the overlap has ended, end simulation
  double differance = protofilaments_[0].pos_[0]-protofilaments_[1].pos_[0];
  if (differance>5000) {
    overlap_ended = true;
  }
}
