#ifndef _CYLAKS_SIMPLE_MOTOR_HPP_
#define _CYLAKS_SIMPLE_MOTOR_HPP_
#include "protein.hpp"
#include "simple_motor_head.hpp"

class SimpleMotorHead;

class SimpleMotor : public Protein {

protected:
  double ran_{0.0};

public:
  SimpleMotorHead head_one_; //, head_two_;

protected:
public:
  SimpleMotor() {}
  void Initialize(size_t sid, size_t id) {
    Object::Initialize(sid, id);
    head_one_.Initialize(sid, id, _r_motor_head, this, nullptr);
    // head_two_.Initialize(sid, id, _r_motor_head, this, &head_one_);
  }
  SimpleMotorHead *GetActiveHead() {
    if (n_heads_active_ != 1) {
      Sys::ErrorExit("Motor::GetActiveHead [0]");
      return nullptr;
    }
    if (head_one_.site_ != nullptr) {
      return &head_one_;
    } else {
      Sys::ErrorExit("SimpleMotor::GetActiveHead [1]");
      return nullptr;
    }
  }
  SimpleMotorHead *GetHeadOne() { return &head_one_; }
  SimpleMotorHead *GetHeadTwo() { return nullptr; } // &head_two_; }
  bool Bind(BindingSite *site, SimpleMotorHead *head);
  bool Step(SimpleMotorHead *head);
};
#endif