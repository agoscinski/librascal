/**
 * file   test_nl.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   2018
 *
 * @brief Example for Neighbour list
 *
 * Copyright  2018 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "structure_managers/structure_manager_centers.hh"
#include "structure_managers/adaptor_strict.hh"
#include "structure_managers/adaptor_neighbour_list.hh"
#include "structure_managers/adaptor_increase_maxorder.hh"
#include "structure_managers/make_structure_manager.hh"
#include "rascal_utility.hh"
#include "representations/representation_manager_sorted_coulomb.hh"
#include "representations/representation_manager_spherical_expansion.hh"
#include "representations/feature_manager_dense.hh"
#include "basic_types.hh"

#include <Eigen/StdVector>

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>
#include <functional>
#include <initializer_list>

// using namespace std;
using namespace rascal;  // NOLINT

/* Bug happens in the structure manager leve of adaptor max order in its loop
 * function, the Iterator of j_cluster (Order 2) has at some point a index
 * out of the range of the vector
 */
int main() {
  bool verbose{true};
  // bool verbose_rep{false};
  double cutoff{2.};

  // std::string
  // filename{"/local/scratch/goscinsk/lib/librascal/tests/simple_cubic_9.json"};
  std::string filename{"/local/scratch/goscinsk/lib/librascal/tests/"
                       "reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};
  auto stucture_manager{make_structure_manager<StructureManagerCenters>()};
  stucture_manager->update(filename);
  auto pair_manager{
      make_adapted_manager<AdaptorNeighbourList>(stucture_manager, cutoff)};
  pair_manager->update();
  auto manager{make_adapted_manager<AdaptorMaxOrder>(pair_manager)};
  manager->update();

  for (auto && center : manager) {
    if (verbose) {
      std::cout << "################################# 2" << std::endl;
      std::cout << center.get_atom_type() << std::endl;
    }
    for (auto neigh : center) {
      if (verbose) {
        std::cout << neigh.get_atom_type() << ", ";
      }
    }
    if (verbose)
      std::cout << std::endl;
  }
  return (0);
}
