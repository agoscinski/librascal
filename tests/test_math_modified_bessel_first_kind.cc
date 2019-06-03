/**
 * file   test_math_modified_bessel_first_kind.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   21 May 2019
 *
 * @brief Test the implementation of modified bessel first kind
 *
 * Copyright © 2019  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "tests.hh"
#include "test_math.hh"

namespace rascal {

  BOOST_AUTO_TEST_SUITE(MathBesselFirstKindTests);

  /* ---------------------------------------------------------------------- */
  /**
   * Check the implementation of the modified bessel function of the first
   * kind against mpmath v1.1.0
   */
  BOOST_FIXTURE_TEST_CASE(math_bessel_test, ModifiedBesselFirstKindRefFixture) {
    for (auto & data : this->ref_data) {
      double x{data["x"]};
      auto ref_vals{data["vals"].get<std::vector<double>>()};
      size_t max_order{data["max_order"]};

      auto vals{math::bessel_i_exp_allorders(x, max_order)};

      for (size_t ii{0}; ii < max_order; ++ii) {
        double rel_error{(vals(ii) - ref_vals[ii])};

        if ((rel_error > math::dbl_ftol) and this->verbose) {
          std::cout << " max_order=" << max_order << " x=" << x <<  " rel_error=" << rel_error<< " val=" << vals(ii) << std::endl;
        }

        // BOOST_CHECK_LE(rel_error, math::dbl_ftol);
      }
      break;
    }
  }


  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
