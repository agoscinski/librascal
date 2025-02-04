/**
 * @file   bind_py_structure_manager.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   9 May 2018
 *
 * @brief  File for binding the Structure Managers
 *
 * Copyright  2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "bind_py_structure_manager.hh"
#include "bind_py_representation_calculator.hh"

namespace rascal {
  /* NOTE: In this file SMI stands for StructureManagerImplementation. */

  // Collection of typedefs to simplify the bindings
  template <typename SMI>
  using LayerByOrder = typename SMI::traits::LayerByOrder;

  template <typename SMI, size_t Order>
  struct HelperLayer {
    static constexpr size_t layer{
        compute_cluster_layer<Order>(LayerByOrder<SMI>{})};
  };

  template <typename SMI, size_t Order>
  using ClusterRefkey_t = ClusterRefKey<Order, HelperLayer<SMI, Order>::layer>;

  template <typename SMI, size_t Order>
  using ClusterRef_t =
      typename StructureManager<SMI>::template ClusterRef<Order>;

  template <typename SMI, size_t Order>
  using PyClusterRef =
      py::class_<ClusterRef_t<SMI, Order>, ClusterRefkey_t<SMI, Order>>;

  template <typename StructureManagerImplementation>
  using PyManager = py::class_<StructureManagerImplementation,
                               StructureManager<StructureManagerImplementation>,
                               std::shared_ptr<StructureManagerImplementation>>;

  //! Bind a ClusterRef (to enable the access to properties such as distances)
  template <size_t Order, size_t Layer>
  void add_cluster_ref(py::module & m) {
    std::string cluster_parent_name =
        internal::GetBindingTypeName<ClusterRefKey<Order, Layer>>();

    py::class_<ClusterRefKey<Order, Layer>, ClusterRefBase>(
        m, cluster_parent_name.c_str());
  }

  //! Bind many cluster refs of a given order using template recursion
  template <size_t Order, size_t Layer, size_t LayerEnd>
  struct add_cluster_refs {
    //! starts recursion
    static void static_for(py::module & m) {
      add_cluster_ref<Order, Layer>(m);
      add_cluster_refs<Order, Layer + 1, LayerEnd>::static_for(m);
    }
  };

  //! Stop the recursion
  template <size_t Order, size_t LayerEnd>
  struct add_cluster_refs<Order, LayerEnd, LayerEnd> {
    static void static_for(py::module &) {}
  };

  //! Bind clusters of different orders
  //! (up to quadruplet)
  template <size_t Order, typename SMI>
  decltype(auto) add_cluster(py::module & m) {
    using ClusterRef = ClusterRef_t<SMI, Order>;
    using ClusterRefKey = ClusterRefkey_t<SMI, Order>;

    std::string cluster_name = internal::GetBindingTypeName<SMI>();

    static_assert(Order < 5, R"(Not currently configured to bind clusters
                                higher than quadruplets.)");

    if (Order == 1) {
      cluster_name += std::string(".Center");
    } else if (Order == 2) {
      cluster_name += std::string(".Pair");
    } else if (Order == 3) {
      cluster_name += std::string(".Triplet");
    } else if (Order == 4) {
      cluster_name += std::string(".Quadruplet");
    }

    py::class_<ClusterRef, ClusterRefKey> py_cluster(m, cluster_name.c_str());
    py_cluster
        .def_property_readonly("atom_tag", &ClusterRef::get_atom_tag,
                               py::return_value_policy::reference)
        .def_property_readonly(
            "atom_type",
            [](const ClusterRef & cluster) { return cluster.get_atom_type(); },
            py::return_value_policy::reference)
        .def_property_readonly(
            "index",
            [](const ClusterRef & cluster) { return cluster.get_index(); },
            py::return_value_policy::reference)
        .def_property_readonly("size", &ClusterRef::size,
                               py::return_value_policy::reference)
        .def_property_readonly("position", &ClusterRef::get_position,
                               py::return_value_policy::reference);
    return py_cluster;
  }

  //! Bind iterator and ClusterRef for Order >= 2
  template <typename StructureManagerImplementation, size_t Order>
  decltype(auto) add_iterator(
      py::module & m,
      PyClusterRef<StructureManagerImplementation, Order - 1> & py_cluster) {
    // Order-1 for PyClusterRef because it is the
    // cluster from the previous iteration
    using Child = StructureManagerImplementation;
    using Parent = typename Child::Parent;

    // bind the iteration over clusterRef<Order-1>
    using ClusterRef = typename Parent::template ClusterRef<Order - 1>;
    py_cluster.def(
        "__iter__",
        [](ClusterRef & v) { return py::make_iterator(v.begin(), v.end()); },
        py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */
    auto py_cluster_new = add_cluster<Order, Child>(m);
    return py_cluster_new;
  }

  //! bind iterator and ClusterRef for Order == 1
  template <typename StructureManagerImplementation, size_t Order>
  decltype(auto)
  add_iterator(py::module & m,
               PyManager<StructureManagerImplementation> & manager) {
    using Child = StructureManagerImplementation;
    using Parent = typename Child::Parent;

    // bind the iteration over the manager
    manager.def(
        "__iter__",
        [](Parent & v) { return py::make_iterator(v.begin(), v.end()); },
        py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */
    auto py_cluster = add_cluster<1, Child>(m);
    return py_cluster;
  }

  /**
   * Bind the clusterRef allowing to iterate over the manager, atom, neigh...
   * Use signature overloading to dispatch to the proper function.
   * Use iteration by recursion to iterate from Order to MaxOrder-1 statically
   */
  template <typename StructureManagerImplementation, size_t Order,
            size_t MaxOrder>
  struct add_iterators {
    //! starts recursion
    static void
    static_for(py::module & m,
               PyManager<StructureManagerImplementation> & manager) {
      auto py_cluster =
          add_iterator<StructureManagerImplementation, Order>(m, manager);
      add_iterators<StructureManagerImplementation, Order + 1,
                    MaxOrder>::static_for(m, py_cluster);
    }
    //! following recursion
    static void static_for(
        py::module & m,
        PyClusterRef<StructureManagerImplementation, Order - 1> & py_cluster) {
      auto py_cluster_new =
          add_iterator<StructureManagerImplementation, Order>(m, py_cluster);
      add_iterators<StructureManagerImplementation, Order + 1,
                    MaxOrder>::static_for(m, py_cluster_new);
    }
  };

  //! Stop the recursion
  template <typename StructureManagerImplementation, size_t MaxOrder>
  struct add_iterators<StructureManagerImplementation, MaxOrder, MaxOrder> {
    static void
    static_for(py::module &,
               PyClusterRef<StructureManagerImplementation, MaxOrder - 1> &) {}
  };

  template <typename Manager_t>
  decltype(auto) add_manager(py::module & mod) {
    using Parent = typename Manager_t::Parent;
    std::string manager_name = internal::GetBindingTypeName<Manager_t>();
    py::class_<Manager_t, Parent, std::shared_ptr<Manager_t>> manager(
        mod, manager_name.c_str());
    return manager;
  }

  /**
   * Is used by the `BindAdaptor` to keep track of the names of the managers
   * to ensure that the same manager present in different stacks of managers
   * is not binded twice.
   */
  template <typename Manager_t>
  decltype(auto) add_manager_safe(py::module & mod,
                                  const std::string & manager_name) {
    using Parent = typename Manager_t::Parent;
    py::class_<Manager_t, Parent, std::shared_ptr<Manager_t>> manager(
        mod, manager_name.c_str());
    return manager;
  }

  /**
   * templated function for adding a StructureManager interface
   * to allow using the iteration over the manager in python, the interface
   * of the structure manager need to be bound.
   */
  template <typename StructureManagerImplementation>
  decltype(auto) add_structure_manager_interface(py::module & m) {
    using Child = StructureManagerImplementation;
    using Parent = typename Child::Parent;

    std::string manager_name = "StructureManager_";
    manager_name += internal::GetBindingTypeName<Child>();
    py::class_<Parent, StructureManagerBase, std::shared_ptr<Parent>> manager(
        m, manager_name.c_str());

    return manager;
  }

  //! templated function for adding a StructureManager implementation
  template <typename StructureManagerImplementation>
  decltype(auto) add_structure_manager_implementation(py::module & m,
                                                      py::module & m_internal) {
    using Child = StructureManagerImplementation;
    constexpr static size_t MaxOrder = Child::traits::MaxOrder;

    auto manager = add_manager<Child>(m);
    manager.def(py::init<>());

    // MaxOrder+1 because it stops at Val-1
    add_iterators<Child, 1, MaxOrder + 1>::static_for(m_internal, manager);
    return manager;
  }

  /**
   * Bind the update function when the atomic structure is provided as
   * a set of positions, the corresponding atom_types, the cell vectors and
   * the periodic boundary conditions.
   */
  template <typename StructureManagerImplementation>
  void
  bind_update_unpacked(PyManager<StructureManagerImplementation> & manager) {
    manager.def("update",
                [](StructureManagerImplementation & manager,
                   const py::EigenDRef<const Eigen::MatrixXd> & positions,
                   const py::EigenDRef<const Eigen::VectorXi> & atom_types,
                   const py::EigenDRef<const Eigen::MatrixXd> & cell,
                   const py::EigenDRef<const Eigen::MatrixXi> & pbc) {
                  manager.update(positions, atom_types, cell, pbc);
                },
                py::arg("positions"), py::arg("atom_types"), py::arg("cell"),
                py::arg("pbc"), py::call_guard<py::gil_scoped_release>());
  }

  /**
   * Bind the update function when no atomic structure is provided. It
   * corresponds to the case when several adaptors are stacked on a
   * structure manager centers and then you update keeping the structure as is.
   */
  template <typename StructureManagerImplementation>
  void bind_update_empty(PyManager<StructureManagerImplementation> & manager) {
    manager.def("update",
                [](StructureManagerImplementation & v) { v.update(); },
                py::call_guard<py::gil_scoped_release>());
  }

  /**
   * Binding utility for the construction of adaptors.
   * Uses partial specialization to define how to bind the constructor of each
   * adaptors.
   *
   * Register Adaptors here
   */
  template <template <class> class Adaptor, typename Implementation_t>
  struct BindAdaptor {};

  template <typename Implementation_t>
  struct BindAdaptor<AdaptorStrict, Implementation_t> {
    using Manager_t = AdaptorStrict<Implementation_t>;
    using ImplementationPtr_t = std::shared_ptr<Implementation_t>;
    using PyManager_t = PyManager<Manager_t>;
    static void bind_adaptor_init(PyManager_t & adaptor) {
      adaptor.def(py::init<std::shared_ptr<Implementation_t>, double>(),
                  py::keep_alive<1, 2>());
    }

    static void bind_adapted_manager_maker(const std::string & name,
                                           py::module & m_adaptor) {
      m_adaptor.def(
          name.c_str(),
          [](ImplementationPtr_t & manager, double cutoff) {
            return make_adapted_manager<AdaptorStrict, Implementation_t>(
                manager, cutoff);
          },
          py::arg("manager"), py::arg("cutoff"), py::return_value_policy::copy);
    }
  };

  template <typename Implementation_t>
  struct BindAdaptor<AdaptorCenterContribution, Implementation_t> {
    using Manager_t = AdaptorCenterContribution<Implementation_t>;
    using ImplementationPtr_t = std::shared_ptr<Implementation_t>;
    using PyManager_t = PyManager<Manager_t>;
    static void bind_adaptor_init(PyManager_t & adaptor) {
      adaptor.def(py::init<std::shared_ptr<Implementation_t>>(),
                  py::keep_alive<1, 2>());
    }

    static void bind_adapted_manager_maker(const std::string & name,
                                           py::module & m_adaptor) {
      m_adaptor.def(name.c_str(),
                    [](ImplementationPtr_t & manager) {
                      return make_adapted_manager<AdaptorCenterContribution,
                                                  Implementation_t>(manager);
                    },
                    py::arg("manager"), py::return_value_policy::copy);
    }
  };

  template <typename Implementation_t>
  struct BindAdaptor<AdaptorMaxOrder, Implementation_t> {
    using Manager_t = AdaptorMaxOrder<Implementation_t>;
    using ImplementationPtr_t = std::shared_ptr<Implementation_t>;
    using PyManager_t = PyManager<Manager_t>;
    static void bind_adaptor_init(PyManager_t & adaptor) {
      adaptor.def(py::init<std::shared_ptr<Implementation_t>>(),
                  py::keep_alive<1, 2>());
    }

    static void bind_adapted_manager_maker(const std::string & name,
                                           py::module & m_adaptor) {
      m_adaptor.def(
          name.c_str(),
          [](ImplementationPtr_t & manager) {
            return make_adapted_manager<AdaptorMaxOrder, Implementation_t>(
                manager);
          },
          py::arg("manager"), py::return_value_policy::copy);
    }
  };

  template <typename Implementation_t>
  struct BindAdaptor<AdaptorHalfList, Implementation_t> {
    using Manager_t = AdaptorHalfList<Implementation_t>;
    using ImplementationPtr_t = std::shared_ptr<Implementation_t>;
    using PyManager_t = PyManager<Manager_t>;
    static void bind_adaptor_init(PyManager_t & adaptor) {
      adaptor.def(py::init<std::shared_ptr<Implementation_t>>(),
                  py::keep_alive<1, 2>());
    }

    static void bind_adapted_manager_maker(const std::string & name,
                                           py::module & m_adaptor) {
      m_adaptor.def(
          name.c_str(),
          [](ImplementationPtr_t & manager) {
            return make_adapted_manager<AdaptorHalfList, Implementation_t>(
                manager);
          },
          py::arg("manager"), py::return_value_policy::copy);
    }
  };

  template <typename Implementation_t>
  struct BindAdaptor<AdaptorFullList, Implementation_t> {
    using Manager_t = AdaptorFullList<Implementation_t>;
    using ImplementationPtr_t = std::shared_ptr<Implementation_t>;
    using PyManager_t = PyManager<Manager_t>;
    static void bind_adaptor_init(PyManager_t & adaptor) {
      adaptor.def(py::init<std::shared_ptr<Implementation_t>>(),
                  py::keep_alive<1, 2>());
    }

    static void bind_adapted_manager_maker(const std::string & name,
                                           py::module & m_adaptor) {
      m_adaptor.def(
          name.c_str(),
          [](ImplementationPtr_t & manager) {
            return make_adapted_manager<AdaptorFullList, Implementation_t>(
                manager);
          },
          py::arg("manager"), py::return_value_policy::copy);
    }
  };

  template <typename Implementation_t>
  struct BindAdaptor<AdaptorNeighbourList, Implementation_t> {
    using Manager_t = AdaptorNeighbourList<Implementation_t>;
    using ImplementationPtr_t = std::shared_ptr<Implementation_t>;
    using PyManager_t = PyManager<Manager_t>;
    static void bind_adaptor_init(PyManager_t & adaptor) {
      adaptor.def(py::init<std::shared_ptr<Implementation_t>, double, bool>(),
                  py::arg("manager"), py::arg("cutoff"),
                  py::arg("consider_ghost_neighbours") = false,
                  py::keep_alive<1, 2>());
    }

    static void bind_adapted_manager_maker(const std::string & name,
                                           py::module & m_adaptor) {
      m_adaptor.def(
          name.c_str(),
          [](ImplementationPtr_t manager, double cutoff,
             bool consider_ghost_neighbours) {
            return make_adapted_manager<AdaptorNeighbourList, Implementation_t>(
                manager, cutoff, consider_ghost_neighbours);
          },
          py::arg("manager"), py::arg("cutoff"),
          py::arg("consider_ghost_neighbours") = false,
          py::return_value_policy::copy);
    }
  };

  template <typename Manager_t>
  void bind_make_structure_manager(py::module & m_str_mng) {
    std::string factory_name = "make_structure_manager_";
    factory_name += internal::GetBindingTypeName<Manager_t>();
    m_str_mng.def(factory_name.c_str(), &make_structure_manager<Manager_t>,
                  py::return_value_policy::copy);
  }

  template <template <class> class Adaptor, typename Manager_t>
  void bind_make_adapted_manager(py::module & m_adaptor) {
    std::string factory_name = "make_adapted_manager_";
    factory_name += internal::GetBindingTypeName<Adaptor<Manager_t>>();
    BindAdaptor<Adaptor, Manager_t>::bind_adapted_manager_maker(factory_name,
                                                                m_adaptor);
  }

  template <class Calculator, class ManagerCollection_t,
            class ManagerCollectionBinder>
  void
  bind_feature_matrix_getter(ManagerCollectionBinder & manager_collection) {
    manager_collection.def(
        "get_dense_feature_matrix",
        &ManagerCollection_t::template get_dense_feature_matrix<Calculator>,
        py::call_guard<py::gil_scoped_release>());
  }

  template <typename Manager, template <class> class... Adaptor>
  void bind_structure_manager_collection(py::module & m_str_mng) {
    using ManagerCollection_t = ManagerCollection<Manager, Adaptor...>;
    using Manager_t = typename ManagerCollection_t::Manager_t;
    std::string factory_name = "ManagerCollection_";
    factory_name += internal::GetBindingTypeName<Manager_t>();
    py::class_<ManagerCollection_t> manager_collection(m_str_mng,
                                                       factory_name.c_str());

    // bind constructor
    manager_collection.def(py::init([](std::string & adaptor_inputs_str) {
      // convert to json
      json hypers = json::parse(adaptor_inputs_str);
      return std::make_unique<ManagerCollection_t>(hypers);
    }));

    // bind iteration over the managers
    manager_collection.def("__iter__",
                           [](ManagerCollection_t & v) {
                             return py::make_iterator(v.begin(), v.end());
                           },
                           py::keep_alive<0, 1>());
    // bind [] accessor
    manager_collection.def(
        "__getitem__",
        [](ManagerCollection_t & v, int index) { return v[index]; },
        py::keep_alive<0, 1>());

    /**
     * Binds the `add_structures`. Instead of invoking the targeted function to
     * bind within a lambda function, a pointer-to-member-function is used here.
     *
     * @param ManagerCollection_t::add_structures a pointer-to-member-function
     * of the add_structures function in ManagerCollection which takes an vector
     * of AtomicStructures<3> and returns void.
     */
    manager_collection.def("add_structures",
                           (void (ManagerCollection_t::*)(  // NOLINT
                               const std::vector<AtomicStructure<3>> &)) &
                               ManagerCollection_t::add_structures);
    /**
     * Binds the `add_structures`. Instead of invoking the targeted function to
     * bind within a lambda function, a pointer-to-member-function is used here.
     *
     * @param ManagerCollection_t::add_structures a pointer-to-member-function
     * of the add_structures function in ManagerCollection which takes a string,
     * const int and int, and returns void.
     */
    manager_collection.def(
        "add_structures",
        (void (ManagerCollection_t::*)(const std::string &, int, int)) &
            ManagerCollection_t::add_structures,
        R"(Read a file and extract the structures from start to start + length.)",
        py::arg("filename"), py::arg("start") = 0, py::arg("length") = -1);

    manager_collection.def("__len__",
                           [](ManagerCollection_t & v) { return v.size(); });

    manager_collection.def("get_parameters", [](ManagerCollection_t & v) {
      return std::string(v.get_adaptors_parameters().dump(2));
    });

    bind_feature_matrix_getter<CalculatorSortedCoulomb, ManagerCollection_t>(
        manager_collection);

    // Register Calculator here
    bind_feature_matrix_getter<CalculatorSphericalExpansion,
                               ManagerCollection_t>(manager_collection);
    bind_feature_matrix_getter<CalculatorSphericalInvariants,
                               ManagerCollection_t>(manager_collection);
    bind_feature_matrix_getter<CalculatorSphericalCovariants,
                               ManagerCollection_t>(manager_collection);
  }

  /**
   * Bind a list of adaptors by stacking them using template recursion.
   * @tparam ManagerImplementation a fully typed manager
   * @tparam AdaptorImplementationPack list of adaptor partial type
   */
  template <typename ManagerImplementation,
            template <class> class AdaptorImplementation,
            template <class> class... AdaptorImplementationPack>
  struct BindAdaptorStack {
    using Manager_t = AdaptorImplementation<ManagerImplementation>;
    using Parent = typename Manager_t::Parent;
    constexpr static size_t MaxOrder = Manager_t::traits::MaxOrder;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using ManagerPtr = std::shared_ptr<Manager_t>;
    using type = BindAdaptorStack<Manager_t, AdaptorImplementationPack...>;

    BindAdaptorStack(py::module & m_nl, py::module & m_adaptor,
                     py::module & m_internal, std::set<std::string> & name_list)
        : next_stack{m_nl, m_adaptor, m_internal, name_list} {
      std::string manager_name{internal::GetBindingTypeName<Manager_t>()};
      if (not name_list.count(manager_name)) {
        name_list.insert(manager_name);
        add_structure_manager_interface<Manager_t>(m_internal);
        auto adaptor = add_manager_safe<Manager_t>(m_adaptor, manager_name);
        BindAdaptor<AdaptorImplementation,
                    ManagerImplementation>::bind_adaptor_init(adaptor);
        // bind_update_empty<Manager_t>(adaptor);
        bind_update_unpacked<Manager_t>(adaptor);
        // bind clusterRefs so that one can loop over adaptor
        // MaxOrder+1 because recursion stops at Val-1
        add_iterators<Manager_t, 1, MaxOrder + 1>::static_for(m_internal,
                                                              adaptor);
        // bind the factory function
        bind_make_adapted_manager<AdaptorImplementation, ManagerImplementation>(
            m_nl);
      }
    }

    type next_stack;
  };

  //! End of recursion
  template <typename ManagerImplementation,
            template <class> class AdaptorImplementation>
  struct BindAdaptorStack<ManagerImplementation, AdaptorImplementation> {
    using Manager_t = AdaptorImplementation<ManagerImplementation>;
    using Parent = typename Manager_t::Parent;
    constexpr static size_t MaxOrder = Manager_t::traits::MaxOrder;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using ManagerPtr = std::shared_ptr<Manager_t>;

    BindAdaptorStack(py::module & m_nl, py::module & m_adaptor,
                     py::module & m_internal,
                     std::set<std::string> & name_list) {
      std::string manager_name{internal::GetBindingTypeName<Manager_t>()};
      if (not name_list.count(manager_name)) {
        name_list.insert(manager_name);
        add_structure_manager_interface<Manager_t>(m_internal);
        auto adaptor = add_manager_safe<Manager_t>(m_adaptor, manager_name);
        BindAdaptor<AdaptorImplementation,
                    ManagerImplementation>::bind_adaptor_init(adaptor);
        bind_update_empty<Manager_t>(adaptor);
        bind_update_unpacked<Manager_t>(adaptor);
        // bind clusterRefs so that one can loop over adaptor
        // MaxOrder+1 because recursion stops at Val-1
        add_iterators<Manager_t, 1, MaxOrder + 1>::static_for(m_internal,
                                                              adaptor);
        // bind the factory function
        bind_make_adapted_manager<AdaptorImplementation, ManagerImplementation>(
            m_nl);
      }
    }
  };

  //! Template overloading of the binding of the structure managers
  template <typename StructureManagerImplementation>
  void bind_structure_manager(py::module & mod, py::module & m_internal);

  template <typename StructureManagerCenters>
  void bind_structure_manager(py::module & mod, py::module & m_internal) {
    using Manager_t = StructureManagerCenters;

    // bind parent class
    add_structure_manager_interface<Manager_t>(m_internal);
    // bind implementation class
    auto manager =
        add_structure_manager_implementation<Manager_t>(mod, m_internal);
    //
    bind_update_unpacked<Manager_t>(manager);
  }

  //! Bind the ClusterRef up to order 4 and from Layer 0 to 6
  void bind_cluster_refs(py::module & m_internal) {
    add_cluster_refs<1, 0, 6>::static_for(m_internal);
    add_cluster_refs<2, 0, 6>::static_for(m_internal);
    add_cluster_refs<3, 0, 6>::static_for(m_internal);
    add_cluster_refs<4, 0, 6>::static_for(m_internal);
  }

  template <typename T>
  using ArrayConstRef_t =
      const Eigen::Ref<const Eigen::Array<T, Eigen::Dynamic, 1>>;

  /**
   * bind AtomicStructure class and bind a vector of them so that a vector of
   * AtomicStructure can be passed from python to c++ without copy to the
   * ManagerCollection.
   */
  void bind_atomic_structure(py::module & mod) {
    using AtomicStructure_t = AtomicStructure<3>;
    py::class_<AtomicStructure_t>(mod, "AtomicStructure")
        .def(py::init<>())
        .def("get_positions",
             [](AtomicStructure_t & atomic_structure) {
               return atomic_structure.positions;
             },
             py::return_value_policy::reference_internal)
        .def("get_atom_types",
             [](AtomicStructure_t & atomic_structure) {
               return atomic_structure.atom_types;
             },
             py::return_value_policy::reference_internal)
        .def("get_cell",
             [](AtomicStructure_t & atomic_structure) {
               return atomic_structure.cell;
             },
             py::return_value_policy::reference_internal)
        .def("get_pbc",
             [](AtomicStructure_t & atomic_structure) {
               return atomic_structure.pbc;
             },
             py::return_value_policy::reference_internal);

    using AtomicStructureList_t = std::vector<AtomicStructure<3>>;

    py::class_<AtomicStructureList_t>(mod, "AtomicStructureList")
        .def(py::init<>())
        .def("append",
             [](AtomicStructureList_t & v,
                const py::EigenDRef<const Eigen::MatrixXd> & positions,
                const py::EigenDRef<const Eigen::VectorXi> & atom_types,
                const py::EigenDRef<const Eigen::MatrixXd> & cell,
                const py::EigenDRef<const Eigen::MatrixXi> & pbc) {
               v.emplace_back();
               v.back().set_structure(positions, atom_types, cell, pbc);
             },
             py::arg("positions"), py::arg("atom_types"), py::arg("cell"),
             py::arg("pbc"), py::call_guard<py::gil_scoped_release>())
        .def("append",
             [](AtomicStructureList_t & v,
                const py::EigenDRef<const Eigen::MatrixXd> & positions,
                const py::EigenDRef<const Eigen::VectorXi> & atom_types,
                const py::EigenDRef<const Eigen::MatrixXd> & cell,
                const py::EigenDRef<const Eigen::MatrixXi> & pbc,
                ArrayConstRef_t<bool> center_atoms_mask) {
               v.emplace_back();
               v.back().set_structure(positions, atom_types, cell, pbc);
               v.back().set_atom_property("center_atoms_mask",
                                          center_atoms_mask);
             },
             py::arg("positions"), py::arg("atom_types"), py::arg("cell"),
             py::arg("pbc"), py::arg("center_atoms_mask"),
             py::call_guard<py::gil_scoped_release>())
        .def("__len__",
             [](const AtomicStructureList_t & v) { return v.size(); })
        .def("__iter__",
             [](AtomicStructureList_t & v) {
               return py::make_iterator(v.begin(), v.end());
             },
             py::keep_alive<0, 1>());
  }

  //! Main function to add StructureManagers and their Adaptors
  void add_structure_managers(py::module & m_nl, py::module & m_internal) {
    // Bind StructureManagerBase (needed for virtual inheritance)
    py::class_<StructureManagerBase, std::shared_ptr<StructureManagerBase>>(
        m_internal, "StructureManagerBase");
    // Bind ClusterRefBase (needed for virtual inheritance)
    py::class_<ClusterRefBase>(m_internal, "ClusterRefBase");
    bind_cluster_refs(m_internal);

    py::module m_strc_mng = m_nl.def_submodule("StructureManager");
    m_strc_mng.doc() = "Structure Manager Classes";
    py::module m_adp = m_nl.def_submodule("Adaptor");
    m_adp.doc() = "Adaptor Classes";

    using Manager_t = StructureManagerCenters;
    bind_structure_manager<Manager_t>(m_strc_mng, m_internal);

    // bind the factory function
    bind_make_structure_manager<Manager_t>(m_nl);
    std::set<std::string> name_list{};
    BindAdaptorStack<Manager_t, AdaptorNeighbourList, AdaptorStrict>
        adaptor_stack_1{m_nl, m_adp, m_internal, name_list};

    BindAdaptorStack<Manager_t, AdaptorNeighbourList, AdaptorCenterContribution,
                     AdaptorStrict>
        adaptor_stack_2{m_nl, m_adp, m_internal, name_list};

    // bind the manager collection
    bind_structure_manager_collection<Manager_t, AdaptorNeighbourList,
                                      AdaptorStrict>(m_nl);

    bind_structure_manager_collection<Manager_t, AdaptorNeighbourList,
                                      AdaptorCenterContribution, AdaptorStrict>(
        m_nl);

    bind_atomic_structure(m_nl);
  }
}  // namespace rascal
