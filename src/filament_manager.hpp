#ifndef _CYLAKS_FILAMENT_MANAGER_HPP_
#define _CYLAKS_FILAMENT_MANAGER_HPP_
#include "population.hpp"
#include "protofilament.hpp"
#include "system_namespace.hpp"
#include "system_parameters.hpp"
#include "system_rng.hpp"

class ProteinManager;

class FilamentManager {
private:
  bool up_to_date_{false};

  Vec<double> weight_neighbs_bind_;
  Vec<double> weight_neighbs_unbind_;

  size_t n_bd_iterations_{0};
  double dt_eff_{0.0};

  ProteinManager *proteins_{nullptr};

public:
  bool mobile_{false};
  Vec<Protofilament> proto_;
  Vec<BindingSite *> sites_;

  Map<Str, Population<Object>> unoccupied_;

private:
  void SetParameters();
  void GenerateFilaments();

  bool AllFilamentsImmobile();

  void UpdateForces();
  void UpdateLattice();

public:
  FilamentManager() {}
  void Initialize(ProteinManager *proteins) {
    proteins_ = proteins;
    SetParameters();
    GenerateFilaments();
  }
  void AddPop(Str name, Fn<Vec<Object *>(Object *)> sort) {
    unoccupied_.emplace(name, Population<Object>(name, sort, sites_.size()));
  }
  void AddPop(Str name, Fn<Vec<Object *>(Object *)> sort, Vec<size_t> i_size,
              Vec<int> i_min, Fn<Vec<int>(Object *)> get_i) {
    Vec<size_t> sz{i_size[0], i_size[1], i_size[2], sites_.size()};
    unoccupied_.emplace(name, Population<Object>(name, sort, sz, i_min, get_i));
  }
  void FlagForUpdate() { up_to_date_ = false; }
  void UpdateUnoccupied() {
    if (up_to_date_) {
      return;
    }
    up_to_date_ = true;
    for (auto &&pop : unoccupied_) {
      pop.second.ZeroOut();
    }
    // Add sites to unoccupied_ and update weights
    for (auto &&site : sites_) {
      for (auto &&pop : unoccupied_) {
        pop.second.Sort(site);
      }
      int n_neighbs{site->GetNumNeighborsOccupied()};
      site->SetWeight_Bind(weight_neighbs_bind_[n_neighbs]);
      site->SetWeight_Unbind(weight_neighbs_unbind_[n_neighbs]);
    }
    UpdateLattice();
  }
  void RunBD() {
    if (AllFilamentsImmobile()) {
      return;
    }
    UpdateForces();
    for (int i_itr{0}; i_itr < n_bd_iterations_; i_itr++) {
      for (auto &&filament : proto_) {
        filament.UpdatePosition();
      }
      UpdateForces();
    }
  }
};
#endif