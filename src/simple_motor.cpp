#include "cylaks/simple_motor.hpp"
#include "cylaks/binding_site.hpp"
#include "cylaks/protofilament.hpp"

bool SimpleMotor::Bind(BindingSite *site, SimpleMotorHead *head) {

  if (site->occupant_ != nullptr) {
    return false;
  }
  site->occupant_ = head;
  head->site_ = site;
  n_heads_active_++;
  if (Sys::test_mode_.empty()) {
    return true;
  }
  return true;
}

bool SimpleMotor::Step(SimpleMotorHead *head) {

  if (head->site_ == nullptr) {
    printf("should DEFINITELY not happen (simpleMotors)\n");
    exit(1);
  }
  BindingSite *old_site{head->site_};
  int dx{head->site_->filament_->dx_};
  BindingSite *new_site(head->site_->GetNeighbor(dx));
  if (new_site->IsOccupied()) {
    return false;
  }
  new_site->occupant_ = head;
  head->site_ = new_site;
  old_site->occupant_ = nullptr;
  return true;
}