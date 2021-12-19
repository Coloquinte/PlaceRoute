
import unittest

import numpy as np

from net_topology import NetTopology


class TestNetTopology(unittest.TestCase):
    def test_basics(self):
        """
        Test the creation of a simple topology
        """
        simple_nets = [
            [0, 1, 2],
            [0, 2],
            [0, 3],
        ]
        topo = NetTopology.from_nets(simple_nets)
        self.assertEqual(topo.nb_cells, 4)
        self.assertEqual(topo.nb_nets, 3)

    def test_b2b(self):
        """
        Test that the B2B model yields a matrix of the right size
        """
        simple_nets = [
            [0, 1, 2],
            [0, 2],
            [0, 3],
        ]
        topo = NetTopology.from_nets(simple_nets)
        self.assertEqual(topo.nb_cells, 4)
        self.assertEqual(topo.nb_nets, 3)
        mat, y = topo.matrix_b2b(np.array([0.0, 0.01, 2.0, 3.0]), 0.1)
        self.assertEqual(mat.shape, (4, 4))
        self.assertEqual(y.shape, (4,))
 
    def test_b2b_pen(self):
        """
        Test that the B2B model with penalty yields a correct matrix
        """
        simple_nets = [
            [0, 1, 2, 3]
        ]
        topo = NetTopology.from_nets(simple_nets, weights=[0.0])
        self.assertEqual(topo.nb_cells, 4)
        self.assertEqual(topo.nb_nets, 1)
        penalty = 0.01
        pos = np.array([0.0, 1.0, 2.0, 3.0]) 
        mat, y = topo.matrix_b2b(
            pos, 0.1,
            pos, penalty
        )
        np.testing.assert_allclose(mat.todense().A, 2 * penalty * np.eye(4))
        np.testing.assert_allclose(y, 2 * penalty * pos)

    def test_b2b_gradient(self):
        """
        Test that the B2B model has correct gradients on a simple case
        """
        simple_nets = [
            [0, 1, 2],
        ]
        topo = NetTopology.from_nets(simple_nets)
        self.assertEqual(topo.nb_cells, 3)
        self.assertEqual(topo.nb_nets, 1)
        x = np.array([0.0, 1.0, 2.0])
        mat, y = topo.matrix_b2b(x, 0.1)
        self.assertEqual(mat.shape, (3, 3))
        self.assertEqual(y.shape, (3,))
        np.testing.assert_allclose(mat.todense(), np.array([[0.75, -0.5, -0.25], [-0.5, 1, -0.5], [-0.25, -0.5, 0.75]]))
        np.testing.assert_allclose(y, np.zeros(3))
        np.testing.assert_allclose(mat * x, np.array([-1, 0, 1]))

    def test_b2b_gradient_with_offset(self):
        """
        Test that the B2B model has correct gradients on a case with offsets
        """
        nets = [
            [0, 1, 2],
        ]
        offsets= [
            [0.0, 1.0, 2.0],
        ]
        topo = NetTopology.from_nets(nets, offsets)
        self.assertEqual(topo.nb_cells, 3)
        self.assertEqual(topo.nb_nets, 1)
        x = np.array([0.0, 0.0, 0.0])
        mat, y = topo.matrix_b2b(x, 0.1)
        self.assertEqual(mat.shape, (3, 3))
        self.assertEqual(y.shape, (3,))
        np.testing.assert_allclose(mat.todense(), np.array([[0.75, -0.5, -0.25], [-0.5, 1, -0.5], [-0.25, -0.5, 0.75]]))
        np.testing.assert_allclose(mat * x - y, np.array([-1, 0, 1]))



if __name__ == '__main__':
    unittest.main()
