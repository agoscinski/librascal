/**
 * file   structure_skew_example.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   30 May 2018
 *
 * @brief test neighbourlist for a problematic structures, probably temporary
 *
 * Copyright  2018 markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#include <structure_managers/adaptor_half_neighbour_list.hh>
#include <structure_managers/adaptor_increase_maxorder.hh>
#include <structure_managers/adaptor_neighbour_list.hh>
#include <structure_managers/adaptor_strict.hh>
#include <structure_managers/make_structure_manager.hh>
#include <structure_managers/property.hh>
#include <structure_managers/structure_manager_centers.hh>

// Shorthands for templated types for readability.
using Manager_t = rascal::StructureManagerCenters;
using PairManager_t = rascal::AdaptorNeighbourList<Manager_t>;

using StrictPairManager_t = rascal::AdaptorStrict<PairManager_t>;
using TripletManager_t = rascal::AdaptorMaxOrder<StrictPairManager_t>;

int main() {
  auto manager{rascal::make_structure_manager<Manager_t>()};
  double cutoff{3.5};

  std::string filename{"problematic_structure.json"};

  std::cout << "Reading structure " << filename << std::endl;

  auto pair_manager{rascal::make_adapted_manager<rascal::AdaptorNeighbourList>(
      manager, cutoff, false)};

  auto strict_manager{rascal::make_adapted_manager<rascal::AdaptorStrict>(
      pair_manager, cutoff)};

  strict_manager->update();

  // auto triplet_manager{
  //     rascal::make_adapted_manager<rascal::AdaptorMaxOrder>(strict_manager)};
  // triplet_manager->update(filename);

  // // Iteration over `manager`
  // std::cout << "manager iteration over atoms" << std::endl;
  // for (auto atom : manager) {
  //   std::cout << "atom " << atom.get_atom_index() << " global index "
  //             << atom.get_global_index() << std::endl;
  // }

  // // `pair_manager` provides iteration over atoms and pairs
  // for (auto atom : pair_manager) {
  //   for (auto pair : atom) {
  //     std::cout << "pair (" << atom.get_atom_index() << ", "
  //               << pair.get_atom_index() << " ) global index "
  //               << pair.get_global_index() << std::endl;
  //   }
  // }

  // // `strict_manager` provides iteration over atoms and strict pairs
  // for (auto atom : strict_manager) {
  //   for (auto pair : atom) {
  //     std::cout << "strict pair (" << atom.get_atom_index() << ", "
  //               << pair.get_atom_index() << ") global index "
  //               << pair.get_global_index() << std::endl;
  //   }
  // }

  // // `triplet_manager` provides iteration over atoms, strict pairs and strict
  // // triplets
  // for (auto atom : triplet_manager) {
  //   for (auto pair : atom) {
  //     for (auto triplet : pair) {
  //       std::cout << "triplet (" << atom.get_atom_index() << ", "
  //                 << pair.get_atom_index() << ", " << triplet.get_atom_index()
  //                 << ") global index " << triplet.get_global_index()
  //                 << std::endl;
  //     }
  //   }
  // }
}
