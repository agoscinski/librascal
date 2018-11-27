/**
 * file   test_adaptor_filter.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   26 Nov 2018
 *
 * @brief  Test the filter adaptor
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#include "structure_managers/adaptor_filter.hh"
// TODO: replace the following line after the branch feat/better_tests
// is merged into master
#include "test_structure.hh"
#include "tests.hh"

#include <boost/mpl/list.hpp>

#include <random>

namespace rascal {

  template<class ManagerImplementation, size_t MaxOrder>
  struct FilterFixture
  {
    using Filter_t = AdaptorFilter<ManagerImplementation, MaxOrder>;

    FilterFixture():
      manager{fixture.manager}{}

    ~FilterFixture() = default;

    size_t get_MaxOrder() {return MaxOrder;}

    ManagerFixture<ManagerImplementation> fixture{};
    Filter_t manager;
  };

  using Fixtures = boost::mpl::list<
    FilterFixture<StructureManagerCenters, 1>,
    FilterFixture<StructureManagerLammps, 2> >;

  using FixturesMax1 = boost::mpl::list<
    FilterFixture<StructureManagerCenters, 1>>;
  using FixturesMax2 = boost::mpl::list<
    FilterFixture<StructureManagerLammps, 2>>;
  using FixturesMax3 = boost::mpl::list<>;


  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE(test_adaptor_filter);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(constructor_test, Fix, Fixtures, Fix) {
  }

  /**
   * This test iterates over a Structure manager's atoms, randomly
   * filters them or not using AdaptorFilter, and stores the indices
   * of the retained atoms in a separate backup vector. In a second
   * loop, we check that the content of the AdaptorFilter is identical
   * to the backup vector. The random filtering is meant to catch all
   * corner cases which might have been missed the authors.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(filter_1_test, Fix, FixturesMax1, Fix) {
    std::random_device rd{};
    std::uniform_int_distribution<int> dist(0, 1);
    std::vector<int> atom_indices{};

    for (auto atom: Fix::fixture.manager) {
      const bool include(dist(rd));
      if (include) {
        Fix::manager.add_cluster(atom);
        atom_indices.push_back(atom.get_atom_index());
      }
    }

    size_t counter{0};
    for (auto atom: Fix::manager) {
      BOOST_CHECK_EQUAL(atom.get_atom_index(), atom_indices[counter]);
      counter++;
      const auto & pos_a{atom.get_position()};
      const auto & pos_b{
        this->fixture.manager.get_position(atom.get_atom_index())};
      const auto error{(pos_a-pos_b).norm()};
      BOOST_CHECK_EQUAL(error, 0.);

      const auto & atom_type_a{atom.get_atom_type()};
      const auto & atom_type_b{
        this->fixture.manager.get_atom_type(atom.back())};
      BOOST_CHECK_EQUAL(atom_type_a, atom_type_b);
    }
  }

  /**
   * This test iterates over a Structure manager's pairs, randomly
   * filters them or not using AdaptorFilter, and stores the atom
   * indices of the retained atom pairs in a separate backup
   * vector. In a second loop, we check that the content of the
   * AdaptorFilter is identical to the backup vector. The random
   * filtering is meant to catch all corner cases which might have
   * been missed the authors.
   */
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(filter_2_test, Fix, FixturesMax2, Fix) {
    std::random_device rd{};
    std::uniform_int_distribution<int> dist(0, 1);
    std::vector<std::array<int, 2>> atom_indices{};


    for (auto atom: Fix::fixture.manager) {
      for (auto pair: atom) {
        const bool include{dist(rd)};
        if (include) {
          Fix::manager.add_cluster(pair);
          atom_indices.push_back(pair.get_atom_indices());
        }
      }
    }

    size_t counter{0};
    for (auto atom: Fix::manager) {
      for (auto pair: atom) {
        auto && a{pair.get_atom_indices()};
        auto && b{atom_indices[counter]};
        BOOST_CHECK_EQUAL_COLLECTIONS(a.begin(), a.end(), b.begin(), b.end());

        const auto & pos_a{pair.get_position()};
        const auto & pos_b{
          this->fixture.manager.get_position(pair.get_atom_index())};
        const auto error{(pos_a-pos_b).norm()};
        BOOST_CHECK_EQUAL(error, 0.);

        const auto & atom_type_a{pair.get_atom_type()};
        const auto & atom_type_b{
          this->fixture.manager.get_atom_type(pair.back())};
        BOOST_CHECK_EQUAL(atom_type_a, atom_type_b);
        counter++;
      }
    }
  }


  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // rascal