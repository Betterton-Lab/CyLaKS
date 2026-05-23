#include "cylaks/simple_motor_head.hpp"
#include "cylaks/simple_motor.hpp"

void SimpleMotorHead::Initialize(size_t sid, size_t id, double radius,
                                 SimpleMotor *parent_ptr,
                                 SimpleMotorHead *other_head_ptr) {
  BindingHead::Initialize(sid, id, radius);
  parent_ = parent_ptr;
  BindingHead::parent_ = dynamic_cast<Protein *>(parent_);
  other_head_ = other_head_ptr;
  BindingHead::other_head_ = dynamic_cast<BindingHead *>(other_head_);
}

bool SimpleMotorHead::Unbind() { return parent_->Unbind(this); }

bool SimpleMotorHead::Step() { return parent_->Step(this); }