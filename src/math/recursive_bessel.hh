/**
 * file   recursive_bessel.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 * @author  Max Veit <max.veit@epfl.ch>
 *
 * @date   27 May 2019
 *
 * @brief Implementation of the modified spherical bessel of the 1st kind
 *
 * Copyright  2019  Felix Musil, Max Veit COSMO (EPFL), LAMMM (EPFL)
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


#ifndef SRC_MATH_RECURSIVE_BESSEL_HH_
#define SRC_MATH_RECURSIVE_BESSEL_HH_

#include "math/math_interface.hh"
#include "math/math_utils.hh"

namespace rascal {
  namespace math {

      /**
       * Compute the modified spherical Bessel function of the first kind
       *
       * This function, hereafter denoted MBSF, is notated as i_n(x).  Here,
       * n must be an integer and x must be a real number > 0.  All orders n
       * are computed and returned, up to a specified maximum order.
       *
       * @param x     Argument of the MBSF
       * @param l_max Maximum order to compute
       *
       * @return Eigen::Array (1-D) of Bessel function values
       *
       * Since the MSBF blows up exponentially with increasing argument (x),
       * what is returned here is e^{-x}*i_n(x), which is bounded (decays to
       * leading order as 1/x).  If you want i_n(x) itself, multiply the
       * result by e^x (at your own peril).
       */
      Eigen::ArrayXd bessel_i_exp_allorders(double x, size_t l_max) {
        if (l_max < 1) {
          l_max = 1;
        }
        Eigen::ArrayXd function_values{l_max + 1};
        // Note: Since x is not likely to be close to zero, there is no
        // signficant improvement in accuracy from using std::expm1 instead
        function_values(0) = (1. - std::exp(2.*x)) / (2.*x);
        function_values(1) = ((x - 1.) + std::exp(2.*x)*(x + 1.)) / (2.*x*x)
        for (size_t l_order{2}; l_order < l_max; ++l_order) {
          function_values(l_order) =
              function_values(l_order - 1)*(2.*l_order - 1.) / x
              + function_values(l_order - 2);
        }
        return function_values;
      }

      /**
       * Compute the modified spherical Bessel function of the first kind
       *
       * Vectorized version for multiple arguments (xs)
       *
       * @param xs    Eigen::Array of MBSF arguments
       * @param l_max Maximum order to compute
       *
       * @return Eigen::Array (2-D) of Bessel function values.  The different
       *         arguments (xs) go along the rows, while the column indexes
       *         the orders (n-values).
       */
      Eigen::ArrayXXd bessel_i_exp_allorders(
          const Eigen::Ref<const Eigen::ArrayXd>& xs, size_t l_max) {
        if (l_max < 1) {
          l_max = 1;
        }
        Eigen::ArrayXXd function_values{xs.size(), l_max + 1};
        function_values.col(0) = (1. - Eigen::exp(2.*xs)) / (2.*xs);
        function_values.col(1) = ((xs - 1.) + Eigen::exp(2.*xs)*(xs + 1.))
                               / (2.*xs.square());
        for (size_t l_order{2}; l_order < l_max; ++l_order) {
          function_values.col(l_order) =
              function_values.col(l_order - 1)*(2.*l_order - 1.) / xs
              + function_values.col(l_order - 2);
        }
        return function_values;
      }

      /**
       * Cache to avoid recomputation of the MSBFs
       *
       * Just initialize and use get_values() to access the values
       */
      class ModifiedSphericalBesselCache {
       public:
        ModifiedSphericalBesselCache = default;

        precompute(const size_t& l_max, const size_t& n_max) {
          this->bessel_values.resize(n_max, l_max+1);
          this->l_max = l_max;
        }

        /**
         * Compute all the MBSFs for the given x-values up to the given order
         */
        inline calc(const Eigen::Ref<const Eigen::VectorXd>& x_values, size_t l_max) {
          bessel_values = bessel_i_exp_allorders(x_values, this->l_max);
        }

        ~ModifiedSphericalBesselCache = default;

        /**
         * Return a reference to the precomputed Bessel function values
         *
         * @return Eigen::Array (2-D) of Bessel function values.  The different
         *         arguments (xs) go along the rows, while the column indexes
         *         the orders (n-values).
         *         Note that a reference is returned to avoid unnecessary
         *         copies.
         */
        Eigen::Ref<const Eigen::ArrayXXd> get_values() {
          return bessel_values;
        }

        Eigen::ArrayXXd bessel_values{};
        size_t l_max{};

      }

  }  // namespace math
}  // namespace rascal

#endif  // SRC_MATH_RECURSIVE_BESSEL_HH_
