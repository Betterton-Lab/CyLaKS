#ifndef _CYLAKS_SIMPLE_MOTOR_HEAD_HPP_
#define _CYLAKS_SIMPLE_MOTOR_HEAD_HPP_
#include "cylaks/binding_head.hpp"
class BindingSite;
class SimpleMotor;

class SimpleMotorHead : public BindingHead {

public:
  SimpleMotor *parent_{nullptr};

private:
public:
  SimpleMotorHead() {}
  void Initialize(size_t sid, size_t id, double radius, SimpleMotor *parent_ptr,
                  SimpleMotorHead *other_head_ptr);
  bool Unbind();
  bool Step();
};
#endif