
import numpy as np


class Circuit:
    def __init__(self):
        self._nb_cells = None
        self._nb_nets = None
        self._cell_name = None
        self._cell_height = None
        self._cell_width = None
        self._cell_fixed = None
        self._cell_x = None
        self._cell_y = None
        self._cell_flip_x = None
        self._cell_flip_y = None
        self._net_name = None
        self._net_limits = None
        self._pin_cells = None
        self._pin_x = None
        self._pin_y = None
        self._pl_area = None
        self._row_height = None

    @property
    def nb_cells(self):
        return self._nb_cells

    @property
    def nb_nets(self):
        return self._nb_nets

    @property
    def cell_name(self):
        return self._cell_name

    @property
    def cell_height(self):
        return self._cell_height

    @property
    def cell_width(self):
        return self._cell_width

    @property
    def cell_fixed(self):
        return self._cell_fixed

    def net_pins_cell(self, i):
        return self._pin_cells[self._net_limits[i]:self._net_limits[i+1]]

    def net_pins_x(self, i):
        return self._pin_x[self._net_limits[i]:self._net_limits[i+1]]

    def net_pins_y(self, i):
        return self._pin_y[self._net_limits[i]:self._net_limits[i+1]]

    @property
    def cell_x(self):
        return self._cell_x

    @property
    def cell_y(self):
        return self._cell_y

    @property
    def cell_flip_x(self):
        return self._cell_flip_x

    @property
    def cell_flip_y(self):
        return self._cell_flip_y

    @property
    def pl_area(self):
        return self._pl_area

    @property
    def row_height(self):
        return self._row_height

    def check(self):
        # No member is None
        assert self._nb_cells is not None
        assert self._nb_nets is not None
        assert self._cell_name is not None
        assert self._cell_height is not None
        assert self._cell_width is not None
        assert self._cell_fixed is not None
        assert self._cell_x is not None
        assert self._cell_y is not None
        assert self._cell_flip_x is not None
        assert self._cell_flip_y is not None
        assert self._net_name is not None
        assert self._net_limits is not None
        assert self._pin_cells is not None
        assert self._pin_x is not None
        assert self._pin_y is not None
        assert self._pl_area is not None
        assert self._row_height is not None

        assert isinstance(self._cell_width, np.ndarray)
        assert isinstance(self._cell_height, np.ndarray)
        assert isinstance(self._cell_fixed, np.ndarray)
        assert isinstance(self._cell_x, np.ndarray)
        assert isinstance(self._cell_y, np.ndarray)
        assert isinstance(self._cell_flip_x, np.ndarray)
        assert isinstance(self._cell_flip_y, np.ndarray)
        assert isinstance(self._pin_cells, np.ndarray)
        assert isinstance(self._pin_x, np.ndarray)
        assert isinstance(self._pin_y, np.ndarray)

        assert self._cell_width.dtype == np.int32
        assert self._cell_height.dtype == np.int32
        assert self._cell_fixed.dtype == np.bool
        assert self._cell_x.dtype == np.int32
        assert self._cell_y.dtype == np.int32
        assert self._cell_flip_x.dtype == np.bool
        assert self._cell_flip_y.dtype == np.bool
        assert self._pin_cells.dtype == np.int32
        assert self._pin_x.dtype == np.int32
        assert self._pin_y.dtype == np.int32

        # Check the dimensions
        assert len(self._cell_name) == self.nb_cells
        assert len(self._cell_height) == self.nb_cells
        assert len(self._cell_width) == self.nb_cells
        assert len(self._cell_fixed) == self.nb_cells
        assert len(self._cell_x) == self.nb_cells
        assert len(self._cell_y) == self.nb_cells
        assert len(self._cell_flip_x) == self.nb_cells
        assert len(self._cell_flip_y) == self.nb_cells
        assert len(self._net_limits) == self.nb_nets + 1
        assert len(self._pin_cells) == len(self._pin_x)
        assert len(self._pin_cells) == len(self._pin_y)

        # Pin cells must be valid
        assert np.all(self._pin_cells >= 0)
        assert np.all(self._pin_cells < self.nb_cells)

    def place(self):
        """
        Call the C++ library to place the circuit
        """
        import ctypes
        dll = ctypes.CDLL("./build/libcoloquinte.so")
        c_bool_p = ctypes.POINTER(ctypes.c_char)
        c_int_p = ctypes.POINTER(ctypes.c_int)
        c_float_p = ctypes.POINTER(ctypes.c_float)
        dll.place_ispd(
            self.nb_cells,
            self.nb_nets,
            self._cell_width.ctypes.data_as(c_int_p),
            self._cell_height.ctypes.data_as(c_int_p),
            self._cell_fixed.ctypes.data_as(c_bool_p),
            self._net_limits.ctypes.data_as(c_int_p),
            self._pin_cells.ctypes.data_as(c_int_p),
            self._pin_x.ctypes.data_as(c_int_p),
            self._pin_y.ctypes.data_as(c_int_p),
            self._cell_x.ctypes.data_as(c_int_p),
            self._cell_y.ctypes.data_as(c_int_p),
            self._cell_flip_x.ctypes.data_as(c_bool_p),
            self._cell_flip_y.ctypes.data_as(c_bool_p),
            *self._pl_area
        )

    def benchmark_all_quadratic_models(self):
        for model_type in ["STAR", "BSTAR", "B2B"]:
            for epsilon in [0.1, 1.0, 10.0]:
                self.benchmark_quadratic_models(model_type, 20, epsilon, 0.0)

    def benchmark_quadratic_models(self, model_type, nb_steps, epsilon, relaxation):
        """
        Call the C++ library to o
        """
        model_types = {"STAR": 0, "BSTAR": 1, "B2B": 2}
        assert epsilon > 0.0
        assert 0.0 <= relaxation <= 1.0
        assert nb_steps > 0
        assert model_type in model_types.keys()

        import ctypes
        dll = ctypes.CDLL("./build/libcoloquinte.so")
        c_bool_p = ctypes.POINTER(ctypes.c_char)
        c_int_p = ctypes.POINTER(ctypes.c_int)
        c_float_p = ctypes.POINTER(ctypes.c_float)
        dll.benchmark_quadratic_models(
            self.nb_cells,
            self.nb_nets,
            self._cell_width.ctypes.data_as(c_int_p),
            self._cell_height.ctypes.data_as(c_int_p),
            self._cell_fixed.ctypes.data_as(c_bool_p),
            self._net_limits.ctypes.data_as(c_int_p),
            self._pin_cells.ctypes.data_as(c_int_p),
            self._pin_x.ctypes.data_as(c_int_p),
            self._pin_y.ctypes.data_as(c_int_p),
            self._cell_x.ctypes.data_as(c_int_p),
            self._cell_y.ctypes.data_as(c_int_p),
            self._cell_flip_x.ctypes.data_as(c_bool_p),
            self._cell_flip_y.ctypes.data_as(c_bool_p),
            model_types[model_type],
            nb_steps,
            ctypes.c_float(epsilon),
            ctypes.c_float(relaxation)
        )

    def partition(self):
        """
        Try and partition the circuit using KaHypar's interface
        """
        import kahypar

        hyperedge_indices = [int(i) for i in self._net_limits]
        hyperedges = [int(i) for i in self._pin_cells]

        node_weights = [1 for i in range(self.nb_cells)]
        edge_weights = [1 for i in range(self.nb_nets)]

        k=2
        hypergraph = kahypar.Hypergraph(self.nb_cells, self.nb_nets, hyperedge_indices, hyperedges, k, edge_weights, node_weights)

        context = kahypar.Context()
        context.loadINIconfiguration("/home/alnurn/Documents/kahypar/config/km1_kKaHyPar_eco_sea20.ini")

        context.setK(k)
        context.setEpsilon(0.03)

        kahypar.partition(hypergraph, context)

    def export_hgr(self, filename):
        """
        Export the circuit in hgr format for hypergraph partitioning tools.

        Node areas are not exported
        """
        edges = []
        for i in range(self.nb_nets):
            b = self._net_limits[i]
            e = self._net_limits[i+1]
            pins = [int(self._pin_cells[j]) for j in range(b, e)]
            # Uniquify the pins and ensure that there are at least two
            pins = list(sorted(set(pins)))
            if len(pins) < 2:
                continue
            edges.append(pins)
        self.export_hgr_data(filename, self.nb_cells, edges)

    def export_edge_hgr(self, filename):
        """
        Export the circuit in hgr format to partition the edges, instead of the nodes
        """
        edges = [[] for n in range(self.nb_cells)]
        net_ind = 0
        for i in range(self.nb_nets):
            b = self._net_limits[i]
            e = self._net_limits[i+1]
            pins = [int(self._pin_cells[j]) for j in range(b, e)]
            # Uniquify the pins and ensure that there are at least two
            pins = list(sorted(set(pins)))
            if len(pins) < 2:
                continue
            for p in pins:
                edges[p].append(net_ind)
            net_ind += 1
        self.export_hgr_data(filename, net_ind, edges)


    def export_hgr_data(self, filename, nb_cells, edges):
        with open(filename, "w") as f:
            print(f"{len(edges)} {nb_cells}", file=f)
            for edge in edges:
                pins = " ".join(str(p + 1) for p in edge)
                print(pins, file=f)
