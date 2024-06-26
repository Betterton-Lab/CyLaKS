#include "cylaks/protein_tester.hpp"
#include "cylaks/filament_tester.hpp"

void ProteinTester::Initialize(FilamentTester *filaments) {

  filaments_ = filaments;
  ProteinManager::filaments_ = dynamic_cast<FilamentManager *>(filaments_);
  SetTestMode();
}

void ProteinTester::UpdateFilaments() {

  filaments_->UpdateUnoccupied();
  if (Sys::test_mode_ != "filament_ablation") {
    return;
  }
  if (Sys::i_step_ == Sys::ablation_step_) {
    filaments_->protofilaments_[1].pos_[0] += 200.0;
    filaments_->protofilaments_[1].ForceUpdate();
  }
}

void ProteinTester::ReportTestStatistics() {

  for (auto const &entry : test_stats_) {
    Sys::Log("For event %s:\n", entry.first.c_str());
    for (int index{0}; index < entry.second.size(); index++) {
      auto stats = entry.second[index];
      double p{double(stats.first) / stats.second};
      double ref{test_ref_.at(entry.first)[index]};
      Sys::Log("  p[%i] = %.3g (%.3g expected) [%zu / %zu events]\n", index, p,
               ref, stats.first, stats.second);
    }
  }
}

void ProteinTester::SetTestMode() {

  Sys::Log("\n");
  Sys::Log("Initializing test '%s'\n", Sys::test_mode_.c_str());
  if (Sys::test_mode_ == "filament_ablation") {
    InitializeTest_Filament_Ablation();
  } else if (Sys::test_mode_ == "filament_separation") {
    InitializeTest_Filament_Separation();
  } else if (Sys::test_mode_ == "hetero_tubulin") {
    InitializeTest_Filament_HeteroTubulin();
  } else if (Sys::test_mode_ == "kinesin_mutant") {
    InitializeTest_Motor_Heterodimer();
  } else if (Sys::test_mode_ == "motor_lattice_step") {
    InitializeTest_Motor_LatticeStep();
  } else if (Sys::test_mode_ == "motor_lattice_bind") {
    InitializeTest_Motor_LatticeBind();
  } else if (Sys::test_mode_ == "xlink_diffusion") {
    InitializeTest_Xlink_Diffusion();
  } else if (Sys::test_mode_ == "xlink_bind_ii") {
    InitializeTest_Xlink_Bind_II();
  } else {
    Sys::ErrorExit("ProteinTester::SetTestMode()");
  }
}

void ProteinTester::InitializeTest_Filament_Ablation() {

  using namespace Params;
  Filaments::count = 2;
  Sys::Log("  COUNT = 2\n");
  Filaments::n_sites[0] = Filaments::n_sites[1] = 875;
  Filaments::polarity[0] = Filaments::polarity[1] = 0;
  Sys::Log("  N_SITES[0] = %i\n", Filaments::n_sites[0]);
  Sys::Log("  N_SITES[1] = %i\n", Filaments::n_sites[1]);
  Filaments::translation_enabled[0] = false;
  Filaments::translation_enabled[1] = false;
  Filaments::rotation_enabled = false;
  Filaments::x_initial[0] = 0.0;
  Filaments::x_initial[1] = (Filaments::n_sites[0] - 1) * Filaments::site_size;
  Filaments::y_initial[0] = Filaments::y_initial[1] = 0.0;
  Motors::endpausing_active = false;
  printf("Enter ablation time: ");
  Str response;
  std::getline(std::cin, response);
  double t_ablate{(double)std::stod(response)};
  Sys::ablation_step_ = size_t(std::round(t_ablate / dt));
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  filaments_->Initialize(this);
  InitializeEvents();
}

void ProteinTester::InitializeTest_Filament_Separation() {

  using namespace Params;
  // Initialize filament environment
  if (Filaments::n_sites[0] != Filaments::n_sites[1]) {
    printf("\nError! Filaments must be the same length.\n");
    exit(1);
  }
  Motors::c_bulk = 0.0;
  Xlinks::c_bulk = 1.0;
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  Filaments::immobile_until[0] = 0.0;
  Filaments::immobile_until[1] = 0.0;
  Filaments::translation_enabled[0] = false;
  Filaments::translation_enabled[1] = true;
  Filaments::rotation_enabled = false; // true;
  filaments_->Initialize(this);
  int n_xlinks{Sys::n_xlinks_};
  if (n_xlinks == -1) {
    Str response;
    printf("Microtubules are %zu sites in length.\n", Filaments::n_sites[0]);
    // for (int i_pf{0}; i_pf < filaments_->protofilaments_.size(); i_pf++) {
    //   printf("%zu\n", filaments_->protofilaments_[i_pf].sites_.size());
    // }
    printf("Enter number of crosslinkers to insert: ");
    std::getline(std::cin, response);
    n_xlinks = (int)std::stoi(response);
  }
  Sys::Log("%i crosslinkers initialized.\n", n_xlinks);
  int n_places{(int)filaments_->sites_.size() / 2};
  if (n_xlinks > n_places) {
    printf("\nError! Too many crosslinkers for filament length used.\n");
    exit(1);
  }
  // Randomly place crosslinkers on filaments w/ x = 0
  int site_indices[n_places];
  for (int index{0}; index < n_places; index++) {
    site_indices[index] = index;
  }
  SysRNG::Shuffle(site_indices, n_places, sizeof(int));
  for (int i_xlink{0}; i_xlink < n_xlinks; i_xlink++) {
    Protein *xlink{xlinks_.GetFreeEntry()};
    int i_site{site_indices[i_xlink]};
    BindingSite *site_one{&filaments_->protofilaments_[0].sites_[i_site]};
    BindingSite *site_two{&filaments_->protofilaments_[1].sites_[i_site]};
    bool exe_one{xlink->Bind(site_one, &xlink->head_one_)};
    bool exe_two{xlink->Bind(site_two, &xlink->head_two_)};
    if (exe_one and exe_two) {
      bool still_attached{xlink->UpdateExtension()};
      if (still_attached) {
        xlinks_.AddToActive(xlink);
        filaments_->FlagForUpdate();
      } else {
        Sys::ErrorExit("ProteinTester::InitializeTestEnvironment() [2]");
      }
    } else {
      Sys::ErrorExit("ProteinTester::InitializeTestEnvironment() [1]");
    }
  }
  // Poisson distribution; sampled to predict events w/ variable probabilities
  auto poisson = [&](double p, int n) {
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  auto is_doubly_bound = [&](Object *protein) -> Vec<Object *> {
    if (protein->GetNumHeadsActive() == 2) {
      return {protein->GetHeadOne(), protein->GetHeadTwo()};
    }
    return {};
  };
  xlinks_.AddPop("diffuse_ii_to_rest", is_doubly_bound);
  xlinks_.AddPop("diffuse_ii_fr_rest", is_doubly_bound);
  auto exe_diffuse_fwd = [&](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    bool executed{head->Diffuse(1)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      if (!still_attached) {
      }
      // FIXME had to move this from if statement above -- why ?
    }
    filaments_->FlagForUpdate();
    xlinks_.FlagForUpdate();
  };
  auto exe_diffuse_bck = [&](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    bool executed{head->Diffuse(-1)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      if (!still_attached) {
      }
      // FIXME had to move this from if statement above -- why ?
    }
    filaments_->FlagForUpdate();
    xlinks_.FlagForUpdate();
  };
  auto get_weight_diff_ii_to = [](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    return head->GetWeight_Diffuse(1);
  };
  auto get_weight_diff_ii_fr = [](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    return head->GetWeight_Diffuse(-1);
  };
  kmc_.events_.emplace_back("diffuse_ii_to_rest",
                            xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal(),
                            &xlinks_.sorted_.at("diffuse_ii_to_rest").size_,
                            &xlinks_.sorted_.at("diffuse_ii_to_rest").entries_,
                            poisson, get_weight_diff_ii_to, exe_diffuse_fwd);
  kmc_.events_.emplace_back("diffuse_ii_fr_rest",
                            xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal(),
                            &xlinks_.sorted_.at("diffuse_ii_fr_rest").size_,
                            &xlinks_.sorted_.at("diffuse_ii_fr_rest").entries_,
                            poisson, get_weight_diff_ii_fr, exe_diffuse_bck);
}

void ProteinTester::InitializeTest_Filament_HeteroTubulin() {

  double p_hetero{Sys::p_mutant_};
  if (p_hetero == -1.0) {
    printf("Enter fraction of heterogenous tubulin: ");
    Str response_one;
    std::getline(std::cin, response_one);
    p_hetero = (double)std::stod(response_one);
  }
  if (p_hetero < 0.0 or p_hetero > 1.0) {
    printf("Error. Invalid fraction!\n");
    exit(1);
  }
  double bind_aff{Sys::binding_affinity_};
  if (bind_aff == -1.0) {
    printf("Enter decrease in binding affinity ");
    printf("(e.g., 2 will cut p_bind in half): ");
    Str response_two;
    std::getline(std::cin, response_two);
    bind_aff = (double)std::stod(response_two);
  }
  if (bind_aff <= 0.0) {
    printf("Error. Fractional change must be positive!\n");
    exit(1);
  }
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  filaments_->Initialize(this);
  int n_sites{(int)filaments_->sites_.size()};
  int n_hetero{(int)std::round(n_sites * p_hetero)};
  // Randomly place heterogeneous sites on lattice
  int site_indices[n_sites];
  for (int index{0}; index < n_sites; index++) {
    site_indices[index] = index;
  }
  SysRNG::Shuffle(site_indices, n_sites, sizeof(int));
  for (int i_hetero{0}; i_hetero < n_hetero; i_hetero++) {
    int i_site{site_indices[i_hetero]};
    filaments_->sites_[i_site]->SetBindingAffinity(bind_aff);
  }
  InitializeEvents();
}

void ProteinTester::InitializeTest_Motor_Heterodimer() {

  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  filaments_->Initialize(this);
  //  Binomial probabilitiy distribution; sampled to predict most events
  auto binomial = [&](double p, int n) {
    if (n > 0) {
      return SysRNG::SampleBinomial(p, n);
    } else {
      return 0;
    }
  };
  // Poisson distribution; sampled to predict events w/ variable probabilities
  auto poisson = [&](double p, int n) {
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  // head_one_ is catalytic; head_two_ is passive
  auto exe_bind_i = [&](auto *site, auto *pop, auto *fil) {
    if (Sys::i_step_ < pop->step_active_) {
      return;
    }
    auto entry{pop->GetFreeEntry()};
    Sys::Log("bound motor %zu\n", entry->GetID());
    // always bind catalytic head first
    bool executed{entry->Bind(site, &entry->head_one_)};
    if (executed) {
      pop->AddToActive(entry);
      fil->FlagForUpdate();
    }
  };
  auto weight_bind_i = [](auto *site) { return site->GetWeight_Bind(); };
  auto is_unocc = [](Object *site) -> Vec<Object *> {
    if (!site->IsOccupied()) {
      return {site};
    }
    return {};
  };
  filaments_->AddPop("motors", is_unocc);
  kmc_.events_.emplace_back(
      "bind_i", motors_.p_event_.at("bind_i").GetVal(),
      &filaments_->unoccupied_.at("motors").size_,
      &filaments_->unoccupied_.at("motors").entries_, poisson,
      [&](Object *base) {
        return weight_bind_i(dynamic_cast<BindingSite *>(base));
      },
      [&](Object *base) {
        exe_bind_i(dynamic_cast<BindingSite *>(base), &motors_, filaments_);
      });
  // Bind_II
  auto exe_bind_ii = [](auto *bound_head, auto *pop, auto *fil) {
    auto head{bound_head->GetOtherHead()};
    auto site{head->parent_->GetNeighbor_Bind_II()};
    auto executed{head->parent_->Bind(site, head)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    }
  };
  auto weight_bind_ii = [](auto *head) {
    return head->parent_->GetWeight_Bind_II();
  };
  auto is_docked = [](auto *motor) -> Vec<Object *> {
    auto *docked_head{motor->GetDockedHead()};
    if (docked_head != nullptr) {
      return {docked_head->GetOtherHead()};
    }
    return {};
  };
  motors_.AddPop("bind_ii", [&](Object *base) {
    return is_docked(dynamic_cast<Motor *>(base));
  });
  kmc_.events_.emplace_back(
      "bind_ii", motors_.p_event_.at("bind_ii").GetVal(),
      &motors_.sorted_.at("bind_ii").size_,
      &motors_.sorted_.at("bind_ii").entries_, poisson,
      [&](Object *base) {
        return weight_bind_ii(dynamic_cast<CatalyticHead *>(base));
      },
      [&](Object *base) {
        exe_bind_ii(dynamic_cast<CatalyticHead *>(base), &motors_, filaments_);
      });
  // Unbind_II
  auto exe_unbind_ii = [](auto *head, auto *pop, auto *fil) {
    bool executed{head->Unbind()};
    if (executed) {
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    }
  };
  auto weight_unbind_ii = [](auto *head) {
    return head->GetWeight_Unbind_II();
  };
  auto is_ADPP_ii_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 2) {
      // Always unbind active head first if both are ADPP bound
      if (motor->head_one_.ligand_ == CatalyticHead::Ligand::ADPP and
          motor->head_two_.ligand_ == CatalyticHead::Ligand::ADPP) {
        return {&motor->head_one_};
      }
      bool found_head{false};
      CatalyticHead *chosen_head{nullptr};
      if (motor->head_one_.ligand_ == CatalyticHead::Ligand::ADPP) {
        chosen_head = &motor->head_one_;
      }
      if (motor->head_two_.ligand_ == CatalyticHead::Ligand::ADPP) {
        if (chosen_head != nullptr) {
          Sys::ErrorExit("Protein_MGR::is_ADPP_ii_bound()");
        }
        chosen_head = &motor->head_two_;
      }
      if (chosen_head != nullptr) {
        return {chosen_head};
      }
    }
    return {};
  };
  motors_.AddPop("unbind_ii", [&](Object *base) {
    return is_ADPP_ii_bound(dynamic_cast<Motor *>(base));
  });
  kmc_.events_.emplace_back(
      "unbind_ii", motors_.p_event_.at("unbind_ii").GetVal(),
      &motors_.sorted_.at("unbind_ii").size_,
      &motors_.sorted_.at("unbind_ii").entries_, poisson,
      [&](Object *base) {
        return weight_unbind_ii(dynamic_cast<CatalyticHead *>(base));
      },
      [&](Object *base) {
        exe_unbind_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                      filaments_);
      });
  // Unbind_I: Unbind first (singly bound) head of a protein
  auto exe_unbind_i = [](auto *head, auto *pop, auto *fil) {
    bool executed{head->Unbind()};
    if (executed) {
      head->UntetherSatellite();
      pop->RemoveFromActive(head->parent_);
      fil->FlagForUpdate();
    }
  };
  auto weight_unbind_i = [](auto *head) {
    return head->parent_->GetWeight_Unbind_I();
  };
  auto is_ADPP_i_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 1) {
      if (motor->GetActiveHead()->ligand_ == CatalyticHead::Ligand::ADPP) {
        return {motor->GetActiveHead()};
      }
    }
    return {};
  };
  motors_.AddPop("bound_i_ADPP", [&](Object *base) {
    return is_ADPP_i_bound(dynamic_cast<Motor *>(base));
  });
  kmc_.events_.emplace_back(
      "unbind_i", motors_.p_event_.at("unbind_i").GetVal(),
      &motors_.sorted_.at("bound_i_ADPP").size_,
      &motors_.sorted_.at("bound_i_ADPP").entries_, poisson,
      [&](Object *base) {
        return weight_unbind_i(dynamic_cast<CatalyticHead *>(base));
      },
      [&](Object *base) {
        exe_unbind_i(dynamic_cast<CatalyticHead *>(base), &motors_, filaments_);
      });
  // Bind_ATP
  auto exe_bind_ATP = [](auto *head, auto *pop) {
    // printf("boop\n");
    bool executed{head->parent_->Bind_ATP(head)};
    if (executed) {
      pop->FlagForUpdate();
    }
  };
  auto is_NULL_i_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 1) {
      if (motor->GetActiveHead()->ligand_ == CatalyticHead::Ligand::NONE) {
        return {motor->GetActiveHead()};
      }
    }
    return {};
  };
  motors_.AddPop("bound_i_NULL", [&](Object *base) {
    return is_NULL_i_bound(dynamic_cast<Motor *>(base));
  });
  kmc_.events_.emplace_back(
      "bind_ATP_i", motors_.p_event_.at("bind_ATP_i").GetVal(),
      &motors_.sorted_.at("bound_i_NULL").size_,
      &motors_.sorted_.at("bound_i_NULL").entries_, binomial,
      [&](Object *base) {
        exe_bind_ATP(dynamic_cast<CatalyticHead *>(base), &motors_);
      });
  // Hydrolyze_ATP
  if (motors_.active_) {
    auto exe_hydrolyze = [](auto *head, auto *pop) {
      bool executed{head->parent_->Hydrolyze(head)};
      if (executed) {
        pop->FlagForUpdate();
      }
    };
    auto is_ATP_i_bound = [](auto *motor) -> Vec<Object *> {
      if (motor->n_heads_active_ == 1) {
        if (motor->GetActiveHead()->ligand_ == CatalyticHead::Ligand::ATP) {
          return {motor->GetActiveHead()};
        }
      }
      return {};
    };
    motors_.AddPop("bound_i_ATP", [&](Object *base) {
      return is_ATP_i_bound(dynamic_cast<Motor *>(base));
    });
    kmc_.events_.emplace_back(
        "hydrolyze", motors_.p_event_.at("hydrolyze").GetVal(),
        &motors_.sorted_.at("bound_i_ATP").size_,
        &motors_.sorted_.at("bound_i_ATP").entries_, binomial,
        [&](Object *base) {
          exe_hydrolyze(dynamic_cast<CatalyticHead *>(base), &motors_);
        });
  }
  // Diffusion
  auto exe_diff = [](auto *head, auto *pop, auto *fil, int dir) {
    bool executed{head->Diffuse(dir)};
    if (executed) {
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    }
  };
  auto is_singly_bound = [](auto *protein) -> Vec<Object *> {
    if (protein->n_heads_active_ == 1) {
      // only head_two can diffuse
      if (protein->GetActiveHead() == &protein->head_two_) {
        return {&protein->head_two_};
      }
    }
    return {};
  };
  Vec<int> i_min{0, 0, 0};
  Vec<size_t> dim_size{1, 1, _n_neighbs_max + 1};
  auto get_n_neighbs = [](auto *entry) {
    Vec<int> indices_vec{entry->GetNumNeighborsOccupied()};
    return indices_vec;
  };
  motors_.AddPop(
      "bound_i",
      [&](Object *base) {
        return is_singly_bound(dynamic_cast<Motor *>(base));
      },
      dim_size, i_min, get_n_neighbs);
  for (int n_neighbs{0}; n_neighbs < _n_neighbs_max; n_neighbs++) {
    kmc_.events_.emplace_back(
        "diffuse_i_fwd", xlinks_.p_event_.at("diffuse_i_fwd").GetVal(n_neighbs),
        &motors_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
        &motors_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs], binomial,
        [&](Object *base) {
          exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_, 1);
        });
    kmc_.events_.emplace_back(
        "diffuse_i_bck", xlinks_.p_event_.at("diffuse_i_bck").GetVal(n_neighbs),
        &motors_.sorted_.at("bound_i").bin_size_[0][0][n_neighbs],
        &motors_.sorted_.at("bound_i").bin_entries_[0][0][n_neighbs], binomial,
        [&](Object *base) {
          exe_diff(dynamic_cast<BindingHead *>(base), &xlinks_, filaments_, -1);
        });
  }
}

// ! FIXME inconsistency
void ProteinTester::InitializeTest_Motor_LatticeStep() {

  using namespace Params;
  Motors::c_bulk = 1.0;
  Motors::t_active = 0.0;
  Motors::n_runs_to_exit = 1000000;
  // Initialize sim objects
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  // Initialize filaments
  Filaments::count = 1;
  Filaments::n_sites[0] = 1000;
  Filaments::translation_enabled[0] = false;
  Filaments::translation_enabled[1] = false;
  Filaments::rotation_enabled = false;
  Sys::Log("  N_SITES[0] = %i\n", Filaments::n_sites[0]);
  filaments_->Initialize(this);
  printf("Enter test delta (-1 to check against self-coop): ");
  Str response;
  std::getline(std::cin, response);
  int test_delta{(int)std::stoi(response)};
  double test_weight_bind{1.0};
  double test_weight_unbind{1.0};
  if (test_delta > 0) {
    test_weight_bind = Sys::weight_lattice_bind_[test_delta];
    test_weight_unbind = Sys::weight_lattice_unbind_[test_delta];
    for (int delta{0}; delta <= Motors::gaussian_range; delta++) {
      Sys::weight_lattice_bind_[delta] = test_weight_bind;
      Sys::weight_lattice_unbind_[delta] = test_weight_unbind;
    }
  }
  for (auto &&site : filaments_->sites_) {
    site->SetWeight_Bind(test_weight_bind);
    site->SetWeight_Unbind(test_weight_unbind);
  }
  // Initialize statistic trackers
  Vec<Pair<size_t, size_t>> zeros(1, {0, 0});
  Vec<double> p_theory_bind_ii(1, motors_.p_event_.at("bind_ii").GetVal() *
                                      test_weight_bind);
  test_stats_.emplace("bind_ii", zeros);
  test_ref_.emplace("bind_ii", p_theory_bind_ii);
  Vec<double> p_theory_unbind_ii(1, motors_.p_event_.at("unbind_ii").GetVal() *
                                        Square(test_weight_unbind));
  test_stats_.emplace("unbind_ii", zeros);
  test_ref_.emplace("unbind_ii", p_theory_unbind_ii);
  Vec<double> p_theory_unbind_i(1, motors_.p_event_.at("unbind_i").GetVal() *
                                       test_weight_unbind);
  test_stats_.emplace("unbind_i", zeros);
  test_ref_.emplace("unbind_i", p_theory_unbind_i);
  // Place motor head on minus end of microtubule
  BindingSite *site{filaments_->protofilaments_[0].minus_end_};
  Motor *motor{motors_.GetFreeEntry()};
  bool executed{motor->Bind(site, &motor->head_one_)};
  if (executed) {
    motors_.AddToActive(motor);
    filaments_->FlagForUpdate();
  } else {
    Sys::ErrorExit("ProteinTester::InitializeTestEnvironment()");
  }
  auto binomial = [&](double p, int n) {
    if (n > 0) {
      return SysRNG::SampleBinomial(p, n);
    } else {
      return 0;
    }
  };
  // Bind_ATP_I
  auto exe_bind_ATP = [](auto *head, auto *pop) {
    // printf("boop\n");
    bool executed{head->parent_->Bind_ATP(head)};
    if (executed) {
      pop->FlagForUpdate();
    }
  };
  auto is_NULL_i_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 1) {
      if (motor->GetActiveHead()->ligand_ == CatalyticHead::Ligand::NONE) {
        return {motor->GetActiveHead()};
      }
    }
    return {};
  };
  motors_.AddPop("bound_i_NULL", [&](Object *base) {
    return is_NULL_i_bound(dynamic_cast<Motor *>(base));
  });
  kmc_.events_.emplace_back(
      "bind_ATP_i", motors_.p_event_.at("bind_ATP_i").GetVal(),
      &motors_.sorted_.at("bound_i_NULL").size_,
      &motors_.sorted_.at("bound_i_NULL").entries_, binomial,
      [&](Object *base) {
        exe_bind_ATP(dynamic_cast<CatalyticHead *>(base), &motors_);
      });
  // Bind_ATP_II
  auto poisson_ATP = [&](double p, int n) {
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  auto exe_bind_ATP_ii = [](auto *front_head, auto *pop, auto *fil) {
    auto *rear_head{front_head->GetOtherHead()};
    if (front_head->trailing_) {
      return;
    }
    bool unbound{rear_head->Unbind()};
    bool executed{front_head->parent_->Bind_ATP(front_head)};
    if (executed) {
      pop->FlagForUpdate();
      fil->FlagForUpdate();
    }
  };
  auto weight_bind_ATP_ii = [](auto *head) {
    return head->parent_->GetWeight_BindATP_II(head);
  };
  auto is_NULL_ii_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 2) {
      bool found_head{false};
      CatalyticHead *chosen_head{nullptr};
      if (motor->head_one_.ligand_ == CatalyticHead::Ligand::NONE) {
        chosen_head = &motor->head_one_;
      }
      if (motor->head_two_.ligand_ == CatalyticHead::Ligand::NONE) {
        if (chosen_head != nullptr) {
          Sys::ErrorExit("Protein_MGR::is_NULL_ii_bound()");
        }
        chosen_head = &motor->head_two_;
      }
      if (chosen_head != nullptr) {
        return {chosen_head};
      }
    }
    return {};
  };
  motors_.AddPop("bound_ii_NULL", [&](Object *base) {
    return is_NULL_ii_bound(dynamic_cast<Motor *>(base));
  });
  kmc_.events_.emplace_back(
      "bind_ATP_ii", motors_.p_event_.at("bind_ATP_ii").GetVal(),
      &motors_.sorted_.at("bound_ii_NULL").size_,
      &motors_.sorted_.at("bound_ii_NULL").entries_, poisson_ATP,
      [&](Object *base) {
        return weight_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base));
      },
      [&](Object *base) {
        exe_bind_ATP_ii(dynamic_cast<CatalyticHead *>(base), &motors_,
                        filaments_);
      });
  // Hydrolyze
  auto exe_hydrolyze = [](auto *head, auto *pop) {
    bool executed{head->parent_->Hydrolyze(head)};
    if (executed) {
      pop->FlagForUpdate();
    }
  };
  auto is_ATP_i_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 1) {
      if (motor->GetActiveHead()->ligand_ == CatalyticHead::Ligand::ATP) {
        return {motor->GetActiveHead()};
      }
    }
    return {};
  };
  motors_.AddPop("bound_i_ATP", [&](Object *base) {
    return is_ATP_i_bound(dynamic_cast<Motor *>(base));
  });
  kmc_.events_.emplace_back(
      "hydrolyze", motors_.p_event_.at("hydrolyze").GetVal(),
      &motors_.sorted_.at("bound_i_ATP").size_,
      &motors_.sorted_.at("bound_i_ATP").entries_, binomial, [&](Object *base) {
        exe_hydrolyze(dynamic_cast<CatalyticHead *>(base), &motors_);
      });
  // Bind_II
  auto exe_bind_ii = [&](Object *base) {
    auto bound_head{dynamic_cast<CatalyticHead *>(base)};
    auto head{bound_head->GetOtherHead()};
    auto site{head->parent_->GetNeighbor_Bind_II()};
    /*
    // If dock site is plus end, unbind motor and place it on minus end
    if (site == site->filament_->plus_end_) {
      bool exe1{bound_head->Unbind()};
      auto new_site{site->filament_->minus_end_};
      bool exe2{bound_head->parent_->Bind(new_site, bound_head)};
      bool exe3{bound_head->parent_->Bind_ATP(bound_head)};
      bool exe4{bound_head->parent_->Hydrolyze(bound_head)};
      site = bound_head->parent_->GetNeighbor_Bind_II();
    }
    */
    auto executed{head->parent_->Bind(site, head)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      motors_.FlagForUpdate();
      filaments_->FlagForUpdate();
      test_stats_.at("bind_ii")[0].first++;
    }
  };
  auto weight_bind_ii = [](auto *head) {
    return head->parent_->GetWeight_Bind_II();
  };
  auto poisson_bind_ii = [&](double p, int n) {
    test_stats_.at("bind_ii")[0].second += motors_.sorted_.at("bind_ii").size_;
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  auto is_docked = [](auto *motor) -> Vec<Object *> {
    auto *docked_head{motor->GetDockedHead()};
    if (docked_head != nullptr) {
      return {docked_head->GetOtherHead()};
    }
    return {};
  };
  motors_.AddPop("bind_ii", [&](Object *base) {
    return is_docked(dynamic_cast<Motor *>(base));
  });
  kmc_.events_.emplace_back(
      "bind_ii", motors_.p_event_.at("bind_ii").GetVal(),
      &motors_.sorted_.at("bind_ii").size_,
      &motors_.sorted_.at("bind_ii").entries_, poisson_bind_ii,
      [&](Object *base) {
        return weight_bind_ii(dynamic_cast<CatalyticHead *>(base));
      },
      exe_bind_ii);
  // Unbind_II
  auto exe_unbind_ii = [&](Object *base) {
    auto head{dynamic_cast<CatalyticHead *>(base)};
    bool executed{head->Unbind()};
    if (executed) {
      motors_.FlagForUpdate();
      filaments_->FlagForUpdate();
      test_stats_.at("unbind_ii")[0].first++;
    }
  };
  auto poisson_unbind_ii = [&](double p, int n) {
    test_stats_.at("unbind_ii")[0].second +=
        motors_.sorted_.at("unbind_ii").size_;
    // printf("sz = %zu\n", motors_.sorted_.at("unbind_ii").size_);
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  auto weight_unbind_ii = [](auto *head) {
    return head->GetWeight_Unbind_II();
  };
  auto is_ADPP_ii_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 2) {
      bool found_head{false};
      CatalyticHead *chosen_head{nullptr};
      if (motor->head_one_.ligand_ == CatalyticHead::Ligand::ADPP) {
        chosen_head = &motor->head_one_;
      }
      if (motor->head_two_.ligand_ == CatalyticHead::Ligand::ADPP) {
        if (chosen_head != nullptr) {
          Sys::ErrorExit("Protein_MGR::is_ADPP_ii_bound()");
        }
        chosen_head = &motor->head_two_;
      }
      if (chosen_head != nullptr) {
        return {chosen_head};
      }
    }
    return {};
  };
  motors_.AddPop("unbind_ii", [&](Object *base) {
    return is_ADPP_ii_bound(dynamic_cast<Motor *>(base));
  });
  kmc_.events_.emplace_back(
      "unbind_ii", motors_.p_event_.at("unbind_ii").GetVal(),
      &motors_.sorted_.at("unbind_ii").size_,
      &motors_.sorted_.at("unbind_ii").entries_, poisson_unbind_ii,
      [&](Object *base) {
        return weight_unbind_ii(dynamic_cast<CatalyticHead *>(base));
      },
      exe_unbind_ii);
  // Unbind_I
  auto exe_unbind_i = [&](Object *base) {
    // Count stats for unbind_i but do not actually execute it
    test_stats_.at("unbind_i")[0].first++;
    filaments_->FlagForUpdate();
  };
  auto poisson_unbind_i = [&](double p, int n) {
    test_stats_.at("unbind_i")[0].second +=
        motors_.sorted_.at("bound_i_ADPP").size_;
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  auto weight_unbind_i = [](auto *head) {
    return head->parent_->GetWeight_Unbind_I();
  };
  auto is_ADPP_i_bound = [](auto *motor) -> Vec<Object *> {
    if (motor->n_heads_active_ == 1) {
      if (motor->GetActiveHead()->ligand_ == CatalyticHead::Ligand::ADPP) {
        return {motor->GetActiveHead()};
      }
    }
    return {};
  };
  motors_.AddPop("bound_i_ADPP", [&](Object *base) {
    return is_ADPP_i_bound(dynamic_cast<Motor *>(base));
  });
  kmc_.events_.emplace_back(
      "unbind_i", motors_.p_event_.at("unbind_i").GetVal(),
      &motors_.sorted_.at("bound_i_ADPP").size_,
      &motors_.sorted_.at("bound_i_ADPP").entries_, poisson_unbind_i,
      [&](Object *base) {
        return weight_unbind_i(dynamic_cast<CatalyticHead *>(base));
      },
      exe_unbind_i);
}

void ProteinTester::InitializeTest_Motor_LatticeBind() {

  using namespace Params;
  size_t cutoff{Motors::gaussian_range};
  // Set parameters
  Xlinks::c_bulk = 0.0;
  Motors::k_on = 1.0;
  Motors::c_bulk = 10.0;
  Motors::neighb_neighb_energy = 0.0;
  Xlinks::neighb_neighb_energy = 0.0;
  Filaments::count = 1;
  Filaments::n_sites[0] = 2 * cutoff + 1;
  // Initialize sim objects
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  filaments_->Initialize(this);
  // Initialize statistic trackers
  Vec<double> p_theory(cutoff + 1, motors_.p_event_.at("bind_i").GetVal());
  for (int delta{0}; delta <= Motors::gaussian_range; delta++) {
    p_theory[delta] *= Sys::weight_lattice_bind_[delta];
  }
  test_ref_.emplace("bind", p_theory);
  Vec<Pair<size_t, size_t>> zeros(cutoff + 1, {0, 0});
  test_stats_.emplace("bind", zeros);
  // Bind motor head to middle site
  int i_site{(int)Motors::gaussian_range};
  BindingSite *site{&filaments_->protofilaments_[0].sites_[i_site]};
  Motor *motor{motors_.GetFreeEntry()};
  bool executed{motor->Bind(site, &motor->head_one_)};
  if (executed) {
    motors_.AddToActive(motor);
    filaments_->FlagForUpdate();
  } else {
    Sys::ErrorExit("ProteinTester::InitializeTestEnvironment()");
  }
  auto poisson = [&](double p, int n) {
    // Each delta distance has 2 sites available to it each timestep
    // Do not count delta = 0, where 'main' motor is permanently bound to
    for (int delta{1}; delta <= Motors::gaussian_range; delta++) {
      test_stats_.at("bind")[delta].second += 2;
    }
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  auto exe_bind_i = [&](Object *base) {
    BindingSite *site{dynamic_cast<BindingSite *>(base)};
    int i_site{(int)site->index_};
    // 'main' kinesin motor will always be at index_ = lattice_cutoff_
    int delta{abs(i_site - (int)Motors::gaussian_range)};
    test_stats_.at("bind")[delta].first++;
    filaments_->FlagForUpdate();
  };
  auto weight_bind_i = [](auto *site) { return site->GetWeight_Bind(); };
  auto is_unocc = [](Object *site) -> Vec<Object *> {
    if (!site->IsOccupied()) {
      return {site};
    }
    return {};
  };
  filaments_->AddPop("motors", is_unocc);
  kmc_.events_.emplace_back(
      "bind_i", motors_.p_event_.at("bind_i").GetVal(),
      &filaments_->unoccupied_.at("motors").size_,
      &filaments_->unoccupied_.at("motors").entries_, poisson,
      [&](Object *base) {
        return weight_bind_i(dynamic_cast<BindingSite *>(base));
      },
      exe_bind_i);
}

void ProteinTester::InitializeTest_Xlink_Diffusion() {

  using namespace Params;
  // Initialize sim objects
  Xlinks::c_bulk = 1.0;
  Xlinks::t_active = 0.0;
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  // Initialize filaments
  Filaments::count = 2;
  Filaments::n_sites[0] = Filaments::n_sites[1] = 1000;
  Sys::Log("  N_SITES[0] = %i\n", Filaments::n_sites[0]);
  Sys::Log("  N_SITES[1] = %i\n", Filaments::n_sites[1]);
  filaments_->Initialize(this);
  // Initialize stat trackers
  double r_y{std::fabs(Filaments::y_initial[0] - Filaments::y_initial[1])};
  double r_max{xlinks_.r_max_};
  printf("r_max = %g\n", r_max);
  double r_x_max{sqrt(Square(r_max) - Square(r_y))};
  printf("r_x_max = %g\n", r_x_max);
  size_t x_max((size_t)std::ceil(r_x_max / Filaments::site_size));
  printf("x_max = %zu\n", x_max);
  Vec<double> p_theory_to(x_max + 1,
                          xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal());
  Vec<double> p_theory_fr(x_max + 1,
                          xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal());
  Vec<double> wt_bind(x_max + 1, 0.0);
  Vec<double> wt_unbind(x_max + 1, 0.0);
  for (int x{0}; x <= x_max; x++) {
    double r_x{x * Filaments::site_size};
    double r{sqrt(Square(r_x) + Square(r_y))};
    if (r < xlinks_.r_min_ or r > xlinks_.r_max_) {
      wt_bind[x] = 0.0;
      wt_unbind[x] = 0.0;
      continue;
    }
    double dr{r - Params::Xlinks::r_0};
    double dE{0.5 * Params::Xlinks::k_spring * Square(dr)};
    wt_bind[x] = exp(-(1.0 - _lambda_spring) * dE / Params::kbT);
    wt_unbind[x] = exp(_lambda_spring * dE / Params::kbT);
  }
  for (int x{0}; x <= x_max; x++) {
    // Ensure index isn't negative to avoid seg-faults
    int x_to{x > 0 ? x - 1 : 0};
    // Diffusing towards rest is considered an unbinding-type event in
    // regards to Boltzmann factors, since both events let the spring relax
    // (dividing Boltzmann factors yields E(x) - E(x_to) in the exponential)
    double wt_spring_to{wt_unbind[x] / wt_unbind[x_to]};
    // Ensure index isn't out of range to avoid seg-faults
    int x_fr{x < x_max ? x + 1 : 0};
    // Diffusing away from rest is considered a binding-type event in
    // regards to Boltzmann factors, since both events stretch the spring
    // (dividing Boltzmann factors yields E(x_fr) - E(x) in the exponential)
    double wt_spring_fr{wt_bind[x_fr] / wt_bind[x]};
    p_theory_to[x] *= wt_spring_to;
    p_theory_fr[x] *= wt_spring_fr;
  }
  test_ref_.emplace("to_rest", p_theory_to);
  test_ref_.emplace("fr_rest", p_theory_fr);
  Vec<Pair<size_t, size_t>> zeros(x_max + 1, {0, 0});
  test_stats_.emplace("to_rest", zeros);
  test_stats_.emplace("fr_rest", zeros);
  // /*
  for (auto &&entry : test_stats_) {
    printf("For '%s':\n", entry.first.c_str());
    for (int x{0}; x < entry.second.size(); x++) {
      printf("  [%i] = {%zu, %zu}\n", x, entry.second[x].first,
             entry.second[x].second);
    }
  }
  // */
  Protein *xlink{xlinks_.GetFreeEntry()};
  int i_site{(int)std::round(Filaments::n_sites[0] / 2)};
  BindingSite *site_one{&filaments_->protofilaments_[0].sites_[i_site]};
  BindingSite *site_two{&filaments_->protofilaments_[1].sites_[i_site]};
  bool exe_one{xlink->Bind(site_one, &xlink->head_one_)};
  bool exe_two{xlink->Bind(site_two, &xlink->head_two_)};
  if (exe_one and exe_two) {
    bool still_attached{xlink->UpdateExtension()};
    if (still_attached) {
      xlinks_.AddToActive(xlink);
      filaments_->FlagForUpdate();
    } else {
      Sys::ErrorExit("ProteinTester::InitializeTestEnvironment() [2]");
    }
  } else {
    Sys::ErrorExit("ProteinTester::InitializeTestEnvironment() [1]");
  }
  auto is_doubly_bound = [](Object *protein) -> Vec<Object *> {
    if (protein->GetNumHeadsActive() == 2) {
      return {protein->GetHeadOne(), protein->GetHeadTwo()};
    }
    return {};
  };
  xlinks_.AddPop("diffuse_ii_to_rest", is_doubly_bound);
  xlinks_.AddPop("diffuse_ii_fr_rest", is_doubly_bound);
  auto exe_diff_to = [&](Object *base) {
    BindingHead *head{dynamic_cast<BindingHead *>(base)};
    double r_x{head->pos_[0] - head->GetOtherHead()->pos_[0]};
    size_t x{(size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
    bool executed{head->Diffuse(1)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      xlinks_.FlagForUpdate();
      filaments_->FlagForUpdate();
      test_stats_.at("to_rest")[x].first++;
    }
  };
  auto exe_diff_fr = [&](Object *base) {
    BindingHead *head{dynamic_cast<BindingHead *>(base)};
    double r_x{head->pos_[0] - head->GetOtherHead()->pos_[0]};
    size_t x{(size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
    bool executed{head->Diffuse(-1)};
    if (executed) {
      bool still_attached{head->parent_->UpdateExtension()};
      xlinks_.FlagForUpdate();
      filaments_->FlagForUpdate();
      test_stats_.at("fr_rest")[x].first++;
    }
  };
  auto weight_diff_ii = [](auto *head, int dir) {
    return head->GetWeight_Diffuse(dir);
  };
  auto poisson_to = [&](double p, int n) {
    Protein *xlink{xlinks_.active_entries_[0]};
    double r_x{xlink->head_one_.pos_[0] - xlink->head_two_.pos_[0]};
    size_t x{(size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
    if (x != 0) {
      test_stats_.at("to_rest")[x].second += 2;
    }
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  auto poisson_fr = [&](double p, int n) {
    Protein *xlink{xlinks_.active_entries_[0]};
    double r_x{xlink->head_one_.pos_[0] - xlink->head_two_.pos_[0]};
    size_t x{(size_t)std::abs(std::round(r_x / Params::Filaments::site_size))};
    if (x != test_stats_.at("fr_rest").size() - 1) {
      BindingSite *site_one{xlink->head_one_.site_};
      BindingSite *site_two{xlink->head_two_.site_};
      if (site_one != site_one->filament_->plus_end_ and
          site_one != site_one->filament_->minus_end_) {
        test_stats_.at("fr_rest")[x].second += 1;
      }
      if (site_two != site_two->filament_->plus_end_ and
          site_two != site_two->filament_->minus_end_) {
        test_stats_.at("fr_rest")[x].second += 1;
      }
    }
    if (p > 0.0) {
      return SysRNG::SamplePoisson(p);
    } else {
      return 0;
    }
  };
  kmc_.events_.emplace_back(
      "diffuse_ii_to_rest", xlinks_.p_event_.at("diffuse_ii_to_rest").GetVal(),
      &xlinks_.sorted_.at("diffuse_ii_to_rest").size_,
      &xlinks_.sorted_.at("diffuse_ii_to_rest").entries_, poisson_to,
      [&](Object *base) {
        return weight_diff_ii(dynamic_cast<BindingHead *>(base), 1);
      },
      exe_diff_to);
  kmc_.events_.emplace_back(
      "diffuse_ii_fr_rest", xlinks_.p_event_.at("diffuse_ii_fr_rest").GetVal(),
      &xlinks_.sorted_.at("diffuse_ii_fr_rest").size_,
      &xlinks_.sorted_.at("diffuse_ii_fr_rest").entries_, poisson_fr,
      [&](Object *base) {
        return weight_diff_ii(dynamic_cast<BindingHead *>(base), -1);
      },
      exe_diff_fr);
}

void ProteinTester::InitializeTest_Xlink_Bind_II() {

  using namespace Params;
  Xlinks::k_off_ii = 143;
  GenerateReservoirs();
  InitializeWeights();
  SetParameters();
  double r_y{std::fabs(Filaments::y_initial[0] - Filaments::y_initial[1])};
  double r_max{xlinks_.r_max_};
  printf("r_max = %g\n", r_max);
  double r_x_max{sqrt(Square(r_max) - Square(r_y))};
  printf("r_x_max = %g\n", r_x_max);
  int x_max((int)std::ceil(r_x_max / Filaments::site_size));
  printf("x_max = %i\n", x_max);
  // Initialize filament environment
  Filaments::count = 2;
  Sys::Log("  COUNT = 2\n");
  Filaments::n_sites[0] = Filaments::n_sites[1] = 2 * x_max + 1;
  Sys::Log("  N_SITES[0] = %i\n", Filaments::n_sites[0]);
  Sys::Log("  N_SITES[1] = %i\n", Filaments::n_sites[1]);
  Filaments::translation_enabled[0] = false;
  Filaments::translation_enabled[1] = false;
  Filaments::rotation_enabled = false;
  filaments_->Initialize(this);
  Vec<double> p_bind(2 * x_max + 1, xlinks_.p_event_.at("bind_ii").GetVal());
  Vec<double> p_unbind(2 * x_max + 1,
                       xlinks_.p_event_.at("unbind_ii").GetVal());
  double offset(Filaments::x_initial[1] - Filaments::x_initial[0]);
  for (int x{-x_max}; x <= x_max; x++) {
    int x_index{x_max + x};
    double r_x{x * Filaments::site_size + offset};
    double r{sqrt(Square(r_x) + Square(r_y))};
    if (r < xlinks_.r_min_ or r > xlinks_.r_max_) {
      p_bind[x_index] *= 0.0;
      p_unbind[x_index] *= 0.0;
      continue;
    }
    double dr{r - Params::Xlinks::r_0};
    double dE{0.5 * Params::Xlinks::k_spring * Square(dr)};
    p_bind[x_index] *= exp(-(1.0 - _lambda_spring) * dE / Params::kbT);
    p_unbind[x_index] *= exp(_lambda_spring * dE / Params::kbT);
  }
  test_ref_.emplace("bind_ii", p_bind);
  test_ref_.emplace("unbind_ii", p_unbind);
  Vec<Pair<size_t, size_t>> zeros(2 * x_max + 1, {0, 0});
  test_stats_.emplace("bind_ii", zeros);
  test_stats_.emplace("unbind_ii", zeros);
  // Place first xlink head on lower MT; remains static for entire sim
  int i_site{(int)x_max};
  BindingSite *site{&filaments_->protofilaments_[0].sites_[i_site]};
  Protein *xlink{xlinks_.GetFreeEntry()};
  bool executed{xlink->Bind(site, &xlink->head_one_)};
  if (executed) {
    xlinks_.AddToActive(xlink);
    filaments_->FlagForUpdate();
  } else {
    Sys::ErrorExit("ProteinTester::InitializeTestEnvironment()");
  }
  // KMC Event -- Bind_II
  // Add population tracker for potential targets of Bind_II event
  auto is_singly_bound = [](Object *protein) -> Vec<Object *> {
    if (protein->GetNumHeadsActive() == 1) {
      return {protein->GetActiveHead()};
    }
    return {};
  };
  xlinks_.AddPop("bind_ii", is_singly_bound);
  // Helper functions for Bind_II event structure
  auto poisson_bind_ii = [&](double p, int n) {
    for (int x{0}; x < test_stats_.at("bind_ii").size(); x++) {
      int n_entries{(int)xlinks_.sorted_.at("bind_ii").size_};
      test_stats_.at("bind_ii")[x].second += n_entries;
    }
    if (p == 0.0) {
      return 0;
    }
    return SysRNG::SamplePoisson(p);
  };
  auto get_weight_bind_ii = [](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    return head->parent_->GetWeight_Bind_II();
  };
  auto exe_bind_ii = [&](Object *base_head) {
    auto bound_head{dynamic_cast<BindingHead *>(base_head)};
    auto head{bound_head->GetOtherHead()};
    auto site{head->parent_->GetNeighbor_Bind_II()};
    auto executed{head->parent_->Bind(site, head)};
    if (executed) {
      xlinks_.FlagForUpdate();
      filaments_->FlagForUpdate();
      bool still_attached{head->parent_->UpdateExtension()};
      if (!still_attached) {
        return;
      }
      double r_x{head->pos_[0] - bound_head->pos_[0]};
      double offset(Filaments::x_initial[1] - Filaments::x_initial[0]);
      int x{(int)std::round((r_x - offset) / Params::Filaments::site_size)};
      int x_max{int(test_stats_.at("unbind_ii").size() - 1) / 2};
      test_stats_.at("bind_ii")[x + x_max].first++;
    } else {
      Sys::ErrorExit("Bind_II (TEST)");
    }
  };
  // Construct KMC event fr Bind_II
  kmc_.events_.emplace_back("bind_ii", xlinks_.p_event_.at("bind_ii").GetVal(),
                            &xlinks_.sorted_.at("bind_ii").size_,
                            &xlinks_.sorted_.at("bind_ii").entries_,
                            poisson_bind_ii, get_weight_bind_ii, exe_bind_ii);
  // KMC event -- Unbind_II
  auto is_doubly_bound = [&](Object *protein) -> Vec<Object *> {
    // Only ever unbind second head
    if (protein->GetNumHeadsActive() == 2) {
      return {protein->GetHeadTwo()};
    }
    return {};
  };
  xlinks_.AddPop("unbind_ii", is_doubly_bound);
  auto poisson_unbind_ii = [&](double p, int n) {
    if (xlinks_.sorted_.at("unbind_ii").size_ > 0) {
      auto head{xlinks_.sorted_.at("unbind_ii").entries_[0]};
      double r_x{head->pos_[0] - head->GetOtherHead()->pos_[0]};
      double offset(Filaments::x_initial[1] - Filaments::x_initial[0]);
      int x{(int)std::round((r_x - offset) / Params::Filaments::site_size)};
      int x_max{int(test_stats_.at("unbind_ii").size() - 1) / 2};
      test_stats_.at("unbind_ii")[x + x_max].second += 1;
    }
    if (p == 0.0) {
      return 0;
    }
    return SysRNG::SamplePoisson(p);
  };
  auto get_weight_unbind_ii = [](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    return head->GetWeight_Unbind_II();
  };
  auto exe_unbind_ii = [&](Object *base) {
    auto head{dynamic_cast<BindingHead *>(base)};
    bool executed{head->Unbind()};
    if (executed) {
      xlinks_.FlagForUpdate();
      filaments_->FlagForUpdate();
      double r_x{head->pos_[0] - head->GetOtherHead()->pos_[0]};
      double offset(Filaments::x_initial[1] - Filaments::x_initial[0]);
      int x{(int)std::round((r_x - offset) / Params::Filaments::site_size)};
      int x_max{int(test_stats_.at("unbind_ii").size() - 1) / 2};
      test_stats_.at("unbind_ii")[x + x_max].first++;
      // bool head->parent_->UpdateExtension();
    } else {
      Sys::ErrorExit("Unbind_II (TEST)");
    }
  };
  kmc_.events_.emplace_back(
      "unbind_ii", xlinks_.p_event_.at("unbind_ii").GetVal(),
      &xlinks_.sorted_.at("unbind_ii").size_,
      &xlinks_.sorted_.at("unbind_ii").entries_, poisson_unbind_ii,
      get_weight_unbind_ii, exe_unbind_ii);
}