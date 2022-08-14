
import numpy as np
import scipy.sparse
import scipy.sparse.linalg


class NetTopology:
    """
    Local representation of the objective function as the sum of lengths of 1D nets
    """

    def __init__(self):
        self._nb_nets = None
        self._nb_cells = None
        self._net_weight = None
        self._net_limits = None
        self._pin_cells = None
        self._pin_offsets = None

    @staticmethod
    def from_nets(nets, offsets=None, weights=None):
        assert isinstance(nets, list)
        for net in nets:
            assert isinstance(net, list)
        if weights is None:
            weights = [1.0 for net in nets]
        if offsets is None:
            offsets = [[0.0 for pin in net] for net in nets]
        assert len(weights) == len(nets)
        assert len(offsets) == len(nets)
        nb_cells = 0
        for net in nets:
            for pin in net:
                nb_cells = max(pin + 1, nb_cells)
        net_limits = [0]
        pin_cells = []
        pin_offsets = []
        for net, net_offsets in zip(nets, offsets):
            assert len(net) == len(net_offsets)
            for pin, pin_offset in zip(net, net_offsets):
                pin_cells.append(pin)
                pin_offsets.append(pin_offset)
            net_limits.append(net_limits[-1] + len(net))
        ret = NetTopology()
        ret._net_weight = np.array(weights)
        ret._net_limits = np.array(net_limits)
        ret._pin_cells = np.array(pin_cells)
        ret._pin_offsets = np.array(pin_offsets)
        ret._nb_cells = nb_cells
        ret._nb_nets = len(nets)
        ret.check()
        return ret

    @staticmethod
    def from_circuit(circuit):
        x_topo = NetTopology()
        y_topo = NetTopology()
        x_topo._nb_nets = circuit.nb_nets
        y_topo._nb_nets = circuit.nb_nets
        x_topo._nb_cells = circuit.nb_cells
        y_topo._nb_cells = circuit.nb_cells
        x_topo._net_weight = np.ones(circuit.nb_nets, dtype=np.float32)
        y_topo._net_weight = np.ones(circuit.nb_nets, dtype=np.float32)
        x_topo._pin_cells = circuit._pin_cells
        y_topo._pin_cells = circuit._pin_cells
        x_topo._net_limits = circuit._net_limits
        y_topo._net_limits = circuit._net_limits
        # Use the centered offsets (instead of lower-left corner)
        x_topo._pin_offsets = np.array([
            p - 0.5 * circuit._cell_width[c]
            for c, p in zip(circuit._pin_cells, circuit._pin_x)
        ], dtype=np.float32)
        y_topo._pin_offsets = np.array([
            p - 0.5 * circuit._cell_height[c]
            for c, p in zip(circuit._pin_cells, circuit._pin_y)
        ], dtype=np.float32)
        x_topo.check()
        y_topo.check()
        return x_topo, y_topo

    @property
    def nb_nets(self):
        """
        Total number of nets
        """
        return self._nb_nets

    @property
    def nb_cells(self):
        """
        Total number of (movable) cells
        """
        return self._nb_cells

    @property
    def net_weight(self):
        """
        Weight associated with the net
        """
        return self._net_weight

    def pin_cells(self, net):
        """
        Cells for the pins of a given net
        """
        return self._pin_cells[self._net_limits[net]:self._net_limits[net+1]]

    def pin_offsets(self, net):
        """
        Offsets for the pins of a given net
        """
        return self._pin_offsets[self._net_limits[net]:self._net_limits[net+1]]

    def matrix_b2b(self, x, epsilon, target=None, penalty=None):
        """
        Quadratic approximation with the bound-to-bound method
        """
        assert isinstance(x, np.ndarray)
        assert len(x) == self.nb_cells
        data = []
        y = np.zeros(self.nb_cells, dtype=np.float32)
        if target is not None:
            assert penalty is not None
            assert isinstance(target, np.ndarray)
            assert len(target) == self.nb_cells
            for i, t in enumerate(target):
                data.append((i, i, 2.0 * penalty))
                y[i] += 2.0 * t * penalty
        for i in range(self.nb_nets):
            weight = self.net_weight[i]
            pins = self.pin_cells(i)
            offsets = self.pin_offsets(i)
            if len(pins) <= 1:
                continue
            elif len(pins) == 2:
                p1 = pins[0]
                p2 = pins[1]
                o1 = offsets[0]
                o2 = offsets[1]
                w = weight / max(abs(x[p1] + o1 - x[p2] - o2), epsilon)
                data.append((p1, p2, -w))
                data.append((p2, p1, -w))
                data.append((p1, p1, w))
                data.append((p2, p2, w))
                y[p1] += w * (o2 - o1)
                y[p2] += w * (o1 - o2)
            else:
                pin_pos = x[pins] + offsets
                mx = np.argmax(pin_pos)
                mn = np.argmin(pin_pos)
                for p in range(len(pins)):
                    if p != mn and p != mx:
                        w = weight / \
                            max(abs(pin_pos[p] - pin_pos[mn]), epsilon)
                        w /= (len(pins) - 1)
                        data.append((pins[p], pins[mn], -w))
                        data.append((pins[mn], pins[p], -w))
                        data.append((pins[p], pins[p], w))
                        data.append((pins[mn], pins[mn], w))
                        y[p] += w * (offsets[mn] - offsets[p])
                        y[mn] += w * (offsets[p] - offsets[mn])
                    if p != mx:
                        w = weight / \
                            max(abs(pin_pos[p] - pin_pos[mx]), epsilon)
                        w /= (len(pins) - 1)
                        data.append((pins[p], pins[mx], -w))
                        data.append((pins[mx], pins[p], -w))
                        data.append((pins[p], pins[p], w))
                        data.append((pins[mx], pins[mx], w))
                        y[p] += w * (offsets[mx] - offsets[p])
                        y[mx] += w * (offsets[p] - offsets[mx])
        rows, cols, vals = zip(*data)
        mat = scipy.sparse.coo_matrix((vals, (rows, cols)), shape=(
            self.nb_cells, self.nb_cells)).tocsr()
        return mat, y

    def matrix_star(self, x, epsilon):
        """
        Quadratic approximation with the star method
        """
        pass

    def matrix_clique(self, x, epsilon):
        """
        Quadratic approximation with the clique method
        """
        pass

    def matrix_linstar(self, x, epsilon):
        """
        Quadratic approximation with the linstar method
        """
        pass

    def _pin_pos(self, net, x):
        ret = []
        for c, o in zip(self.pin_cells(net), self.pin_offsets(net)):
            ret.append(x[c] + o)
        return np.array(ret)

    def val(self, x):
        """
        Value using the actual objective function
        """
        res = 0.0
        for i, w in enumerate(self.net_weight):
            pos = self._pin_pos(i, x)
            res += w * (np.max(pos) - np.min(pos))
        return res

    def val_lse(self, epsilon):
        """
        Smoothed value with the log-sum-exp approximation
        """
        res = 0.0
        for i, w in enumerate(self.net_weight):
            pos = self._pin_pos(i, x)
            res += w * (_lse(pos) - _lse(-pos))
        return res

    def val_wa(self, epsilon):
        """
        Smoothed value with the weighted average approximation
        """
        res = 0.0
        for i, w in enumerate(self.net_weight):
            pos = self._pin_pos(i, x)
            res += w * (_wa(pos) - _wa(-pos))
        return res

    def grad(self):
        """
        Gradient using the actual objective function
        """
        res = np.zeros(self.nb_cells)
        for i, w in enumerate(self.net_weight):
            pos = self._pin_pos(i, x)
            c = self.pin_cells(i)
            res[c[np.argmin(pos)]] -= w
            res[c[np.argmax(pos)]] += w
        return res

    def grad_lse(self, epsilon):
        """
        Smoothed gradient with the log-sum-exp approximation
        """
        pass

    def grad_wa(self, epsilon):
        """
        Smoothed gradient with the weighted average approximation
        """
        pass

    def optim_b2b(self, epsilon, target, penalty):
        x = target
        print(f"Initial HPWL: {self.val(x):.3E}")
        for i in range(100):
            val_penalty = penalty * np.sum((x - target)**2)
            val_hpwl = self.val(x)
            mat, y = self.matrix_b2b(x, epsilon, target, penalty)
            x, info = scipy.sparse.linalg.cg(mat, y, x0=x)
            print(
                f"Iter #{i+1}: HPWL {val_hpwl:.3E}, penalty {val_penalty:.3E}")
            #import pdb; pdb.set_trace()

    def check(self):
        # No member is None
        assert self._nb_nets is not None
        assert self._nb_cells is not None
        assert self._net_weight is not None
        assert self._net_limits is not None
        assert self._pin_cells is not None
        assert self._pin_offsets is not None

        # Check the dimensions
        assert len(self._net_weight) == self.nb_nets
        assert len(self._net_limits) == self.nb_nets + 1
        assert len(self._pin_cells) == len(self._pin_offsets)

        # Pin cells must be valid
        assert np.all(self._pin_cells >= 0)
        assert np.all(self._pin_cells < self.nb_cells)

        # Valid non-negative weights
        assert np.all(self._net_weight >= 0)

        # Sorted and valid net limits
        assert np.all(np.diff(self._net_limits) >= 0)
        assert self._net_limits[0] == 0
        assert self._net_limits[-1] == len(self._pin_cells)


def _lse(data, epsilon):
    cmax = np.max(data)
    data = data - cmax
    return cmax + epsilon * np.log(np.exp(data / epsilon).sum())


def _wa(data, epsilon):
    cmax = np.max(data)
    data = data - cmax
    e = np.exp(data / epsilon)
    return cmax + (data * e).sum() / e.sum()
