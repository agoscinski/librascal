#!/usr/bin/env python3

import unittest
import faulthandler

from python_structure_manager_test import (
    TestStructureManagerCenters, TestNL, TestNLStrict
)
from python_representation_calculator_test import (
    TestSortedCoulombRepresentation, TestSphericalExpansionRepresentation,
    TestSphericalInvariantsRepresentation
)

from python_math_test import TestMath
from python_test_sparsify_fps import TestFPS

if __name__ == '__main__':
    faulthandler.enable()

    unittest.main()
