#==============================================================================
# file   CMakeLists.txt
#
# @author Till Junge <till.junge@epfl.ch>
#
# @date   05 Apr 2018
#
# @brief  configuration for compiling structure_managers
#
# @section LICENSE
#
# Copyright  2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
#
# rascal is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3, or (at
# your option) any later version.
#
# rascal is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GNU Emacs; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
#

set (manager_SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/structure_manager_lammps.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/structure_manager_centers.cc
  )

target_sources("${LIBRASCAL_NAME}" PRIVATE ${manager_SRC})
