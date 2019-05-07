/**
 * file   adaptor_increase_maxorder.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   19 Jun 2018
 *
 * @brief implements an adaptor for structure_managers, which
 * creates a full and half neighbourlist if there is none and
 * triplets/quadruplets, etc. if existent.
 *
 * Copyright  2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_ADAPTOR_INCREASE_MAXORDER_HH_
#define SRC_STRUCTURE_MANAGERS_ADAPTOR_INCREASE_MAXORDER_HH_

#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "rascal_utility.hh"
#include "lattice.hh"
#include "basic_types.hh"

#include <typeinfo>
#include <set>
#include <vector>


namespace rascal {
  namespace internal {
    
    template <typename T>
    constexpr T nb_distances(T order) {
      return (order*(order-1))/2;
    }
  }
  /**
   * Forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorMaxOrder;

  /**
   * Specialisation of traits for increase <code>MaxOrder</code> adaptor
   */
  template <class ManagerImplementation>
  struct StructureManager_traits<AdaptorMaxOrder<ManagerImplementation>> {
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    // inherit strictness
    // constexpr static AdaptorTraits::Strict Strict{
    //    ManagerImplementation::traits::Strict};
    constexpr static bool HasDistances{false};
    // TODO(alex) when implemented put flag on true
    // constexpr static bool HasDistances{true};
    constexpr static bool HasDirectionVectors{
        ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    // New MaxOrder upon construction
    constexpr static size_t MaxOrder{ManagerImplementation::traits::MaxOrder +
                                     1};
    constexpr static AdaptorTraits::NeighbourListType NeighbourListType{
        AdaptorTraits::NeighbourListType::half};
    // TODO(alex) change above to below. Above is done to match the current
    // behaviour of tests, conflicts currently with ANL inits type full list but
    // AMO extends with type half list, resolve this when Behler-Parinello
    // symmetry functions have been implemented
    //
    // constexpr static AdaptorTraits::NeighbourListType NeighbourListType{
    //    ManagerImplementation::traits::NeighbourListType};

    // Extend the layer by one with the new MaxOrder
    using LayerByOrder = typename LayerExtender<
        MaxOrder, typename ManagerImplementation::traits::LayerByOrder>::type;
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Adaptor that increases the MaxOrder of an existing StructureManager. This
   * means, if the manager does not have a neighbourlist, there is nothing this
   * adaptor can do (hint: use adaptor_neighbour_list before and stack this on
   * top), if it exists, triplets, quadruplets, etc. lists are created.
   */
  template <class ManagerImplementation>
  class AdaptorMaxOrder
      : public StructureManager<AdaptorMaxOrder<ManagerImplementation>>,
        public std::enable_shared_from_this<
            AdaptorMaxOrder<ManagerImplementation>> {
   public:
    using Manager_t = AdaptorMaxOrder<ManagerImplementation>;
    using Parent = StructureManager<Manager_t>;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using traits = StructureManager_traits<AdaptorMaxOrder>;
    using AtomRef_t = typename ManagerImplementation::AtomRef_t;
    template <size_t Order>
    using ClusterRef_t =
        typename ManagerImplementation::template ClusterRef<Order>;
    using Vector_ref = typename Parent::Vector_ref;
    using Hypers_t = typename Parent::Hypers_t;

    static_assert(traits::MaxOrder > 2,
                  "ManagerImplementation needs at least a pair list for"
                  " extension.");

    //! Default constructor
    AdaptorMaxOrder() = delete;

    /**
     * Given at least a pair list, this adaptor creates the next Order
     * list. I.e. from pairs to triplets, triplets to quadruplet, etc. Not
     * cutoff is needed, the cutoff is implicitly given by the neighbourlist,
     * which was built
     */
    explicit AdaptorMaxOrder(ImplementationPtr_t manager);

    AdaptorMaxOrder(ImplementationPtr_t manager, std::tuple<>)
        : AdaptorMaxOrder(manager) {}

    AdaptorMaxOrder(ImplementationPtr_t manager,
                    const Hypers_t & /*adaptor_hypers*/)
        : AdaptorMaxOrder(manager) {}

    //! Copy constructor
    AdaptorMaxOrder(const AdaptorMaxOrder & other) = delete;

    //! Move constructor
    AdaptorMaxOrder(AdaptorMaxOrder && other) = default;

    //! Destructor
    virtual ~AdaptorMaxOrder() = default;

    //! Copy assignment operator
    AdaptorMaxOrder & operator=(const AdaptorMaxOrder & other) = delete;

    //! Move assignment operator
    AdaptorMaxOrder & operator=(AdaptorMaxOrder && other) = default;

    /**
     * Updates just the adaptor assuming the underlying manager was
     * updated. this function invokes making triplets, quadruplets,
     * etc. depending on the MaxOrder, pair list has to be present.
     */
    void update_self();

    //! Updates the underlying manager as well as the adaptor
    template <class... Args>
    void update(Args &&... arguments);

    /**
     * Returns the linear indices of the clusters (whose atom indices are stored
     * in counters). For example when counters is just the list of atoms, it
     * returns the index of each atom. If counters is a list of pairs of indices
     * (i.e. specifying pairs), for each pair of indices i,j it returns the
     * number entries in the list of pairs before i,j appears.
     */
    template <size_t Order>
    inline size_t
    get_offset_impl(const std::array<size_t, Order> & counters) const;


    //! returns the distance between atoms in a given pair

    template <size_t Order, size_t Layer>
    inline const std::conditional_t<Order == 2, double,
                           std::array<double, internal::nb_distances(Order)>> &
    get_distance(const ClusterRefKey<Order, Layer> & cluster) {
      if(traits::HasDistances) {
        if (traits::MaxOrder > Order) {
          return this->manager->get_distance(cluster);
        } 
        return get_distances(cluster);
      } 
      throw std::runtime_error("Access on distance"
                               " without underlying AdaptorFilter "
                               " with distance property.");        
    }
    // 1. get_iterator_at_pair_index (the the atom indices the cluster
    // 2. get all pairs from an array<size_t,..>
    template <size_t Order, size_t Layer>
    const inline std::conditional_t<Order == 2, double,
           std::array<double, internal::nb_distances(Order)>> &
      get_distances(const ClusterRefKey<Order, Layer> & cluster) {
      const size_t nb_distances = internal::nb_distances(Order);
      std::conditional_t<Order == 2, double,
           std::array<double, nb_distances>> distances[nb_distances];
      size_t distances_index{0};
      auto & root_manager{cluster.get_manager()};
      for (auto atom_index : cluster.get_atom_indices()) {
        auto && iterator_at_atom_index{root_manager->get_iterator_at(atom_index)};
        auto atom_at_atom_index{*iterator_at_atom_index};
        for (auto pair : atom_at_atom_index) {
          distances.at(distances_index) = this->get_distance(pair);
        }
      }
      return distances;
    }

    //! Returns the number of clusters of size cluster_size
    inline size_t get_nb_clusters(size_t order) const {
      switch (order) {
      case traits::MaxOrder: {
        return this->neighbours.size();
        break;
      }
      default:
        return this->manager->get_nb_clusters(order);
        break;
      }
    }

    //! Returns number of clusters of the original manager
    inline size_t get_size() const { return this->manager->get_size(); }

    //! Returns position of an atom with index atom_index
    inline Vector_ref get_position(const size_t & atom_index) {
      return this->manager->get_position(atom_index);
    }

    //! Returns position of the given atom object (useful for users)
    inline Vector_ref get_position(const AtomRef_t & atom) {
      return this->manager->get_position(atom.get_index());
    }

    //! get atom type from underlying manager
    inline const int & get_atom_type(const int & atom_index) const {
      return this->manager->get_atom_type(atom_index);
    }

    //! get atom type from underlying manager
    inline int & get_atom_type(const int & atom_index) {
      return this->manager->get_atom_type(atom_index);
    }

    //! return atom type
    inline int & get_atom_type(const AtomRef_t & atom) {
      return this->manager->get_atom_type(atom.get_atom_index());
    }

    //! get atom_index of the index-th atom in manager
    inline int get_cluster_neighbour(const Parent &, size_t index) const {
      return this->manager->get_cluster_neighbour(*this->manager, index);
    }

    //! Returns the id of the index-th neighbour atom of a given cluster
    template <size_t Order, size_t Layer>
    inline int
    get_cluster_neighbour(const ClusterRefKey<Order, Layer> & cluster,
                          size_t index) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation only handles up to traits::MaxOrder");

      // necessary helper construct for static branching
      using IncreaseHelper_t =
          internal::IncreaseHelper<Order == (traits::MaxOrder - 1)>;

      if (Order < (traits::MaxOrder - 1)) {
        return IncreaseHelper_t::get_cluster_neighbour(*this->manager, cluster,
                                                       index);
      } else {
        auto && offset = this->offsets[cluster.get_cluster_index(Layer)];
        return this->neighbours[offset + index];
      }
    }

    //! Returns the number of neighbors of a given cluster
    template <size_t Order, size_t Layer>
    inline size_t
    get_cluster_size(const ClusterRefKey<Order, Layer> & cluster) const {
      static_assert(Order < traits::MaxOrder,
                    "this implementation handles only the respective MaxOrder");
      /*
       * Here it is traits::MaxOrder-1, because only the current manager has the
       * right answer to the number of neighbours of the MaxOrder-1 tuple. This
       * is the 'else' case.
       */

      // necessary helper construct for static branching
      using IncreaseHelper_t =
          internal::IncreaseHelper<Order == (traits::MaxOrder - 1)>;

      if (Order < (traits::MaxOrder - 1)) {
        return IncreaseHelper_t::get_cluster_size(*this->manager, cluster);
      } else {
        auto access_index = cluster.get_cluster_index(Layer);
        return this->nb_neigh[access_index];
      }
    }

    //! Get the manager used to build the instance
    ImplementationPtr_t get_previous_manager() {
      return this->manager->get_shared_ptr();
    }

   protected:
    //! Extends the list containing the number of neighbours with a 0
    inline void add_new_entry_to_number_of_neighbours() {
      this->nb_neigh.push_back(0);
    }
    inline void add_new_neighbours_of_cluster(std::vector<size_t> neighbours) {
      this->nb_neigh.push_back(0);
      if (neighbours.size() > 0) {
        for (auto neighbour : neighbours) {
          this->add_neighbour_to_most_recent_added_cluster(neighbour);
        }
      }
    }

    inline void
    add_neighbour_to_most_recent_added_cluster(const int atom_index) {
      // adds `atom_index` to neighbours
      this->neighbours.push_back(atom_index);
      // increases the number of neighbours
      this->nb_neigh.back()++;
    }

    //! Sets the correct offsets for accessing neighbours
    inline void set_offsets() {
      auto n_tuples{nb_neigh.size()};
      if (n_tuples > 0) {
        this->offsets.reserve(n_tuples);
        this->offsets.resize(1);
        for (size_t i{0}; i < n_tuples; ++i) {
          this->offsets.emplace_back(this->offsets[i] + this->nb_neigh[i]);
        }
      }
    }

    //! reference to underlying manager
    ImplementationPtr_t manager;

    //! Construct for reaching the MaxOrder and adding neighbours of at MaxOrder
    template <size_t Order, bool IsOldMaxOrderReached>
    struct AddOrderLoop;

    //! Stores the number of neighbours for every traits::MaxOrder-1-clusters
    std::vector<size_t> nb_neigh{};

    //! Stores all neighbours of traits::MaxOrder-1-clusters
    std::vector<size_t> neighbours{};

    /**
     * Stores the offsets of traits::MaxOrder-1-*clusters for accessing
     * `neighbours`, from where nb_neigh can be counted
     */
    std::vector<size_t> offsets{};
  };

  /* ---------------------------------------------------------------------- */
  //! Constructor of the next level manager
  template <class ManagerImplementation>
  AdaptorMaxOrder<ManagerImplementation>::AdaptorMaxOrder(
      std::shared_ptr<ManagerImplementation> manager)
      : manager{std::move(manager)}, nb_neigh{}, neighbours{}, offsets{} {
    if (traits::MaxOrder < 3) {
      throw std::runtime_error("Increase MaxOrder: No pair list in underlying"
                               " manager.");
    }
  }

  /* ---------------------------------------------------------------------- */
  //! update, involing the update of the underlying manager
  template <class ManagerImplementation>
  template <class... Args>
  void AdaptorMaxOrder<ManagerImplementation>::update(Args &&... arguments) {
    // if sizeof...(arguments) == 0 then the underlying structure
    // is not changed
    if (sizeof...(arguments) > 0) {
      this->set_update_status(false);
    }
    this->manager->update(std::forward<Args>(arguments)...);
  }

  /* ---------------------------------------------------------------------- */
  /* Helper structure to loop through all orders copying the cluster indices
   * list of each order and creating the new MaxOrder list.
   */
  template <class ManagerImplementation>
  template <size_t Order, bool IsOldMaxOrderReached>
  struct AdaptorMaxOrder<ManagerImplementation>::AddOrderLoop {};

  /*
   * Recursion start:
   * The cluster indices lists of the underlying structure manager are
   * recursively copied to this new structure manager until the old max order is
   * reached. Container_t can be StructureManager or ClusterRef.
   * Order is the order of the container, starting with 0 for a structure
   * manager and then corresponds to the ClusterRef's order.
   */
  template <class ManagerImplementation>
  template <size_t Order>
  struct AdaptorMaxOrder<ManagerImplementation>::AddOrderLoop<Order, false> {
    static constexpr int OldMaxOrder{ManagerImplementation::traits::MaxOrder};
    using ClusterRef_t =
        typename ManagerImplementation::template ClusterRef<Order>;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using Container_t =
        std::conditional_t<Order == 0, ImplementationPtr_t, ClusterRef_t>;

    using NextOrderLoop = AddOrderLoop<Order + 1, (Order + 1 == OldMaxOrder)>;

    static void loop(Container_t & container,
                     AdaptorMaxOrder<ManagerImplementation> & manager) {
      for (auto next_cluster : container) {
        auto & next_cluster_indices{
            std::get<(Order + 1) - 1>(manager.cluster_indices_container)};

        auto indices{next_cluster.get_cluster_indices()};
        next_cluster_indices.push_back(indices);

        NextOrderLoop::loop(next_cluster, manager);
      }
    }
  };

  /**
   * Recursion end: The new cluster indices list of order MaxOrder is created.
   * At desired MaxOrder (plus one), here is where the magic happens and the
   * neighbours of the same order are added as the Order+1.  add check for non
   * half neighbour list.
   *
   * TODO: currently, this implementation is not distinguishing between minimal
   * and full lists. E.g. this has to be adjusted to include both, the i- and
   * the j-atoms of each pair as an i-atom in a triplet (center). Should this
   * adaptor be able to build a half neighbour list from a full neighbour list?
   */
  template <class ManagerImplementation>
  template <size_t Order>
  struct AdaptorMaxOrder<ManagerImplementation>::AddOrderLoop<Order, true> {
    static constexpr int OldMaxOrder{ManagerImplementation::traits::MaxOrder};

    using ClusterRef_t =
        typename ManagerImplementation::template ClusterRef<Order>;
    using ClusterRefOrder1_t =
        typename ManagerImplementation::template ClusterRef<1>;

    using traits = typename AdaptorMaxOrder<ManagerImplementation>::traits;

    using strict_yes_t = 
      typename std::integral_constant<AdaptorTraits::Strict, AdaptorTraits::Strict::yes>;
    using strict_no_t = 
      typename std::integral_constant<AdaptorTraits::Strict, AdaptorTraits::Strict::no>;
    using traits_strict_t = 
      typename std::integral_constant<AdaptorTraits::Strict, traits::Strict>;

    using half_t = 
      typename std::integral_constant<AdaptorTraits::NeighbourListType, AdaptorTraits::NeighbourListType::half>;
    using full_t = 
      typename std::integral_constant<AdaptorTraits::NeighbourListType, AdaptorTraits::NeighbourListType::full>;
    using traits_neighbourlisttype_t = 
      typename std::integral_constant<AdaptorTraits::NeighbourListType, traits::NeighbourListType>;

    //! loop through the orders to get to the maximum order, this is agnostic to
    //! the underlying MaxOrder, just goes to the maximum order
    static void loop(ClusterRef_t & cluster,
                     AdaptorMaxOrder<ManagerImplementation> & manager) {
      extend_cluster_indices_container(cluster, manager);
    }

    static inline void extend_cluster_indices_container(
        ClusterRef_t & cluster,
        AdaptorMaxOrder<ManagerImplementation> & manager) {
      auto atom_indices_of_cluster = cluster.get_atom_indices();

      // access to root manager for access to atom pairs
      auto & root_manager{cluster.get_manager()};

      // a set of new neighbours for the cluster, which will be added to extend
      // the cluster for new MaxOrder
      std::set<size_t> atom_indices_in_neighbour_environment_of_cluster{};

      // careful: atom_indices_of_cluster can include ghosts: ghosts have to be
      // ignored, since they to not have a neighbour list themselves, they are
      // only neighbours
      insert_neighbours_of_atom_indices<traits_strict_t>(
            atom_indices_of_cluster,
            atom_indices_in_neighbour_environment_of_cluster,
            root_manager);

      // to remove cluster's atom indices in the cluster's neighbour environment
      std::vector<size_t> neighbours_of_cluster{};
      std::set_difference(
          atom_indices_in_neighbour_environment_of_cluster.begin(),
          atom_indices_in_neighbour_environment_of_cluster.end(),
          atom_indices_of_cluster.begin(), atom_indices_of_cluster.end(),
          std::inserter(neighbours_of_cluster, neighbours_of_cluster.begin()));

      // add an entry for the current clusters' neighbours
      manager.add_new_neighbours_of_cluster(neighbours_of_cluster);
    }
  
    //static inline typename std::enable_if<!traits::Strict>::type insert_neighbours_of_atom_indices(
    template <typename Strictness>
    static inline typename std::enable_if<std::is_same<Strictness,strict_no_t>::value>::type insert_neighbours_of_atom_indices(
        std::array<int, Order> & atom_indices_of_cluster,
        std::set<size_t> & atom_indices_in_neighbour_environment_of_cluster,
        StructureManager<ManagerImplementation> & manager) {
      for (auto atom_index : atom_indices_of_cluster) {
        insert_neighbours_of_atom( 
            atom_index,
            atom_indices_in_neighbour_environment_of_cluster,
            atom_indices_of_cluster.back(),
            manager);
      }
    }

    template <typename Strictness>
    static inline typename std::enable_if<std::is_same<Strictness,strict_yes_t>::value>::type insert_neighbours_of_atom_indices(
        std::array<int, Order> & atom_indices_of_cluster,
        std::set<size_t> & atom_indices_in_neighbour_environment_of_cluster,
        StructureManager<ManagerImplementation> & manager) {
      std::set<size_t> atom_indices_in_neighbour_environment_of_atom{};
      std::set<size_t> current_atom_indices_in_neighbour_environment_of_cluster{};
      insert_neighbours_of_atom( 
          atom_indices_of_cluster.front(),
          atom_indices_in_neighbour_environment_of_cluster,
          atom_indices_of_cluster.back(),
          manager);
      for (auto atom_index : atom_indices_of_cluster) {
        std::set<size_t> current_atom_indices_in_neighbour_environment_of_cluster{};
        insert_neighbours_of_atom( 
            atom_index,
            atom_indices_in_neighbour_environment_of_cluster,
            atom_indices_of_cluster.back(),
            manager); 
        std::set_intersection(
          atom_indices_in_neighbour_environment_of_cluster.begin(),
          atom_indices_in_neighbour_environment_of_cluster.end(),
          atom_indices_in_neighbour_environment_of_atom.begin(),
          atom_indices_in_neighbour_environment_of_atom.end(),
          std::inserter(current_atom_indices_in_neighbour_environment_of_cluster,
            current_atom_indices_in_neighbour_environment_of_cluster.begin()));
        atom_indices_in_neighbour_environment_of_cluster =
            current_atom_indices_in_neighbour_environment_of_cluster;
      }
    }

    static inline void insert_neighbours_of_atom(
        size_t atom_index, std::set<size_t> & set, int index,
        StructureManager<ManagerImplementation> & manager) {
      auto && iterator_at_atom_index{manager.get_iterator_at(atom_index)};
      ClusterRefOrder1_t atom_at_atom_index{*iterator_at_atom_index};
      insert_neighbours_of_atom<traits_neighbourlisttype_t>(
          atom_at_atom_index,
          set,
          index);
    }

    // pushes the neighbours of the atom to the set
    template <typename NeighbourListType>
    static inline typename std::enable_if<std::is_same<NeighbourListType, half_t>::value>::type insert_neighbours_of_atom(
        ClusterRefOrder1_t & atom, std::set<size_t> & set, int index) {
      for (auto pair : atom) {
        auto neighbour_of_atom_at_atom_access_index =
            pair.get_internal_cluster_neighbour_index();
        if (neighbour_of_atom_at_atom_access_index > index) {
          set.insert(neighbour_of_atom_at_atom_access_index);
        }
      }
    }

    template <typename NeighbourListType>
    static inline typename std::enable_if<std::is_same<NeighbourListType, full_t>::value>::type insert_neighbours_of_atom(
        ClusterRefOrder1_t & atom, std::set<size_t> & set, int) {
      for (auto pair : atom) {
        auto neighbour_of_atom_at_atom_access_index =
            pair.get_internal_cluster_neighbour_index();
        set.insert(neighbour_of_atom_at_atom_access_index);
      }
    }
    // TODO(alex) this function does not work because of dereference issues
    //static inline ClusterRefOrder1_t get_atom_at_atom_index_from_manager(
    //    StructureManager<ManagerImplementation> & underlying_manager,
    //    size_t atom_index) {
    //  return *(underlying_manager.get_iterator_at(atom_index));
    //}
  };




  /* ---------------------------------------------------------------------- */
  /**
   * This is the loop, which runs recursively goes to the maximum Order and then
   * increases it by one (i.e. pairs->triplets, triplets->quadruplets, etc.
   */
  template <class ManagerImplementation>
  void AdaptorMaxOrder<ManagerImplementation>::update_self() {
    static_assert(traits::MaxOrder > 2,
                  "No neighbourlist present; extension not possible.");

    internal::for_each(this->cluster_indices_container,
                       internal::ResizePropertyToZero());

    this->nb_neigh.clear();
    this->offsets.clear();
    this->neighbours.clear();

    AddOrderLoop<0, false>::loop(this->manager, *this);

    // correct the offsets for the new cluster order
    this->set_offsets();

    // add correct cluster_indices for the highest order
    auto & max_cluster_indices{
        std::get<traits::MaxOrder - 1>(this->cluster_indices_container)};
    max_cluster_indices.fill_sequence();
  }

  /* ---------------------------------------------------------------------- */
  /**
   * Returns the linear indices of the clusters (whose atom indices are stored
   * in counters). For example when counters is just the list of atoms, it
   * returns the index of each atom. If counters is a list of pairs of indices
   * (i.e. specifying pairs), for each pair of indices i,j it returns the number
   * entries in the list of pairs before i,j appears.
   */
  template <class ManagerImplementation>
  template <size_t Order>
  inline size_t AdaptorMaxOrder<ManagerImplementation>::get_offset_impl(
      const std::array<size_t, Order> & counters) const {
    static_assert(Order < traits::MaxOrder,
                  "this implementation handles only up to the respective"
                  " MaxOrder");
    // Order accessor: 0 - atoms
    //                 1 - pairs
    //                 2 - triplets
    //                 etc.
    // Order is determined by the ClusterRef building iterator, not by the Order
    // of the built iterator

    // necessary construct for static branching
    using IncreaseHelper_t =
        internal::IncreaseHelper<Order == (traits::MaxOrder - 1)>;

    if (Order < (traits::MaxOrder - 1)) {
      // If not accessible at this order, call lower Order offsets from lower
      // order manager or push through to lower levels, if adaptors are stacked.
      return IncreaseHelper_t::get_offset_impl(*this->manager, counters);
    } else {
      // Counters is an array to call parent offset multiplet. This can then be
      // used to access the actual offset for the Order which was built here.
      // It needs to be cast into a smaller one to access the order of this
      // Cluster(Order-1) from the manager below
      std::array<size_t, Order - 1> counters_below{};
      for (size_t c_index{0}; c_index < Order - 1; ++c_index) {
        counters_below[c_index] = counters[c_index];
      }
      // Linear index of the Cluster (Order-1)
      auto i{this->manager->get_offset_impl(counters_below)};
      // Number of cluster in its current iteration
      auto j{counters[Order - 1]};
      auto tuple_index{i + j};
      auto main_offset{this->offsets[tuple_index]};
      return main_offset;
    }
  }
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_ADAPTOR_INCREASE_MAXORDER_HH_
