.. _how_to_test:

.. contents::
   :local:


How to write a test
-------------------

Every feature (e.g., functions, class methods and constructors, algorithms) needs its own unit test. Unit tests serve two main purposes, on the one hand, they allow test-driven development (I.e.,  you define a test case and your expected results, then develop your feature. Once you replicate the expected results, your feature is ready) and on the other hand, they help catching regressions, especially in combination with the continuous integration server (It runs all test cases after every commit, and complains if the change causes any of them to fail.

*Rascal* uses the `boost unit testing framework <http://www.boost.org/doc/libs/1_66_0/libs/test/doc/html/index.html>`_ for unit tests of the C++ core library and `unittest module <https://docs.python.org/3/library/unittest.html>`_ of the Python standard library for Python binding tests.

It is instructive to go through the documentations and tutorials for both testing frameworks for details, as the following examples only serve as pointers in the right direction.

Writing a Boost test
^^^^^^^^^^^^^^^^^^^^

Tests can be added to any of the ``test_*.cc`` files in the ``tests`` folder, or you can add a new file containing tests in a new file, as long as it follows the pattern ``test_*.cc`` (after adding a new file, you will have to run ``cmake .`` in the build folder for CMake to pick up the modification).

A test file needs to have the following structure:

.. code-block:: c++

   #include "tests.hh"

   namespace rascal {

     ... // test cases go here

   } // rascal

Any test case in such a file will be added to *Rascal*\' main test suite. It is recommended to group test cases that logically belong together in sub test suites using the ``BOOST_AUTO_TEST_SUITE`` macro. Imagine we write a new sub-suite called ``tutorial_test``

.. code-block:: c++

   #include "tests.hh"

   namespace rascal {

     BOOST_AUTO_TEST_SUITE(tutorial_test);

     ... // test cases go here

     BOOST_AUTO_TEST_SUITE_END()

   } // rascal

The most used types of test cases will very likely be ``BOOST_AUTO_TEST_CASE`` (for straight-forward test cases that do not share common code with other test cases) and ``BOOST_FIXTURE_TEST_CASE_TEMPLATE`` (for testing more involved features which require a setup phase and are parametrised by template parameters, see `Fixtures <http://www.boost.org/doc/libs/1_66_0/libs/test/doc/html/boost_test/tests_organization/fixtures.html>`_ for a detailed discussion)

Writing a ``BOOST_AUTO_TEST_CASE``
..................................
This is as simple as running some function from the library and checking the results with ``BOOST_TEST`` (`here <http://www.boost.org/doc/libs/1_63_0/libs/test/doc/html/boost_test/testing_tools/boost_test_universal_macro.html>`_), e.g.:

.. code-block:: c++

   #include "tests.hh"
   #include "module.hh"

   namespace rascal {

     BOOST_AUTO_TEST_SUITE(tutorial_test);

     BOOST_AUTO_TEST_CASE(f_test) {
       BOOST_TEST(f(12) == 2, "Should have been 2");
     }

     BOOST_AUTO_TEST_SUITE_END()

   } // rascal


Writing a ``BOOST_FIXTURE_TEST_CASE_TEMPLATE``
..............................................

While the previous example was simple, it was also very limited. Frequently, we wish to test multiple properties and methods of an initialised class or data structure, and will do so for multiple template parameters. The following very contrived example creates a so-called templated *test fixture*, defines a list of template instantiations that we wish to test, and runs the test cases on each member of the list.


.. code-block:: c++

   #include "tests.hh"
   #include <boost/mpl/list.hpp>

   namespace rascal {

     BOOST_AUTO_TEST_SUITE(tutorial_test);

     // creation of the test fixture. In practice, this structure would
     // contain data members (here `int val`) that correspond to some data
     // structure of rascal. The constructor (which is required to be a
     // *default constructor*, i.e., without parameters) initialises the
     //structure)
     template <int Dim>
     struct DemoTestFixture {

       static constexpr int dim(){return Dim;}

       DemoTestFixture()
         :val{Dim}
       {}

       int val;
     };

     // create a list of template instantiations to test:
     using fixtures = boost::mpl::list<DemoTestFixture<2>,
                                       DemoTestFixture<3>>;

     // declare a fixture test using the list
     BOOST_FIXTURE_TEST_CASE_TEMPLATE(
       templated_basic_fixture_test, Fix, fixtures, Fix) {
       BOOST_TEST(Fix::val == Fix::dim());
     }

     BOOST_AUTO_TEST_SUITE_END()

   } // rascal



Writing a binding test
^^^^^^^^^^^^^^^^^^^^^^

Similarly to the library tests, binding tests can be added to any of the ``python_*.py`` files in the ``tests`` folder, or you can add a new file containing tests in a new file, as long as it follows the pattern ``python_*.py`` (after adding a new file, you will have to run ``cmake .`` in the build folder for CMake to pick up the modification). Furthermore, if you add a new file (say, ``python_tutorial_test.py``, you will have to import it in the main python test file ``python_binding_tests.py``:

.. code-block:: py

   import python_tutorial_test

The basic unit test tool in Python's ``unittest`` module is the ``unittest.TestCase`` class. New test cases need to inherit from it, define a test case initialisation method ``setUp(self)`` and one or more test methods ``test_*(self)``. Say we create a new test case to test the distance matrix calculation function ``cdist``:


.. code-block:: py

   import unittest
   import numpy as np
   from python_import_rascal import _rascal as pt

   class TestCdist(unittest.TestCase):
       def setUp(self):
           """builds the test case. we'll use it to create the matrices
           coordinate matrices A and B between which we wish to compute
           the distances

           """
           self.A = np.array([[0., 0.],
                              [1., 0.],
                              [2., 0.]])
           nb_A = self.A.shape[0]
           self.B = np.array([[0., 1.],
                              [1., 1.]])
           nb_B = self.B.shape[0]

           # the distance matrix is trivial to compute:

           self.dists_ref = np.empty([nb_A, nb_B])

           for i in range(nb_A):
               for j in range(nb_B):
                   self.dists_ref[i, j] = np.linalg.norm(
                       self.A[i, :] - self.B[j, :])


       def test_cdist(self):
           """feeds the matrices A and B to rascal' cdist function and compares
           the results to the local reference dist_ref
           """
           dists = pt.cdist(self.A, self.B)

           error = np.linalg.norm(dists-self.dists_ref)
           tol = 1e-10
           self.assertLessEqual(error, tol)

How to run the tests
--------------------

All tests are added as targets during  compilation by default. You can run all tests by executing ``ctest`` in the build folder. If the tests fail, you can re-run them verbosely using ``ctest -V`` to get a detailed account of which tests have failed.

