/**
 * file   hyp1f1.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 * @author  Max Veit <max.veit@epfl.ch>
 *
 * @date   27 May 2019
 *
 * @brief Implementation of the gauss legendre quadrature weights and points
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


#ifndef SRC_MATH_GAUSS_LEGENDRE_HH_
#define SRC_MATH_GAUSS_LEGENDRE_HH_

#include "math/math_interface.hh"
#include "math/math_utils.hh"

#include <map>

namespace rascal {
  namespace math {
    /**
     * @params r_st starting point of the integral
     * @params r_nd ending point of the integral
     */
    // 
    inline MatrixX2_t compute_gauss_legendre_points_weigths(const double& r_st, const double& r_nd, const int& order_n) {
      MatrixX2_t point_weight(order_n, 2);
      double z{0.}, pp{0.}, r_midpoint{(r_st+r_nd)*0.5}, r_length{0.5*(r_nd - r_st)};
      for (int ii{1}; ii <= static_cast<int>((order_n+1)/2); ++ii) {
        z = std::cos(math::PI * (ii - 0.25) / (order_n + 0.5));
        double z1{0.};
        while (std::abs(z-z1) > dbl_ftol) {
          double p1{1.}, p2{0.};
          for (int jj{1}; jj <= order_n; ++jj) {
            double p3{p2};
            p2 = p1;
            p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
          }
          pp = order_n*(z*p1-p2)/(z*z-1.0);
          z1 = z;
          z = z1 - p1/pp;
        }

        point_weight(ii-1, 0) = r_midpoint - r_length*z;
        point_weight(order_n+2-ii, 0) = r_midpoint + r_length*z;
        point_weight(ii-1, 1) = (2.0*r_length)/((1.0-z*z)*pp*pp);
        point_weight(order_n+2-ii, 1) = point_weight(ii-1, 1);
      }
      return point_weight;
    }

  }  // namespace math
}  // namespace rascal

#endif  // SRC_MATH_GAUSS_LEGENDRE_HH_