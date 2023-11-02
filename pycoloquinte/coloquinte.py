"""
Coloquinte VLSI placer
"""


import gzip
import lzma
import math
import os
import sys

import coloquinte_pybind
from coloquinte_pybind import (
    CellOrientation,
    CellRowPolarity,
    ColoquinteParameters,
    GlobalPlacerParameters,
    LegalizationParameters,
    DetailedPlacerParameters,
    LegalizationModel,
    NetModel,
    PlacementStep,
    Rectangle,
    Row,
)


def _open_file(name, write=False):
    """
    Open the file with the appropriate decompression method. In read mode, search for compressed versions if the exact name does not exist.
    """
    mode = "wt" if write else "rt"
    if name.endswith(".gz"):
        return gzip.open(name, mode=mode)
    elif name.endswith(".xz") or name.endswith(".lzma"):
        return lzma.open(name, mode=mode)
    elif write:
        return open(name, mode=mode)
    elif os.path.exists(name):
        return open(name, mode=mode)
    elif os.path.exists(name + ".gz"):
        return gzip.open(name + ".gz", mode=mode)
    elif os.path.exists(name + ".xz"):
        return lzma.open(name + ".xz", mode=mode)
    elif os.path.exists(name + ".lzma"):
        return lzma.open(name + ".lzma", mode=mode)
    else:
        raise RuntimeError(f"Could not find file {name}")


def _read_aux(filename):
    if os.path.isdir(filename):
        dir_list = os.listdir(filename)
        all_files = [f for f in dir_list if f.endswith(".aux")]
        default_name = os.path.basename(os.path.normpath(filename)) + ".aux"
        if len(all_files) == 1:
            filename = os.path.join(filename, all_files[0])
        elif len(all_files) > 1 and default_name in all_files:
            filename = os.path.join(filename, default_name)
        else:
            raise RuntimeError(
                f"There should be one file ending with .aux, got {', '.join(all_files)}"
            )
    elif not os.path.exists(filename):
        filename = filename + ".aux"
    dirname = os.path.dirname(filename)
    files = []
    with open(filename) as f:
        for line in f:
            for name in line.split():
                files.append(name)
    node_files = [n for n in files if n.endswith(".nodes")]
    net_files = [n for n in files if n.endswith(".nets")]
    pl_files = [n for n in files if n.endswith(".pl")]
    scl_files = [n for n in files if n.endswith(".scl")]
    if len(node_files) != 1:
        raise RuntimeError("There should be a .nodes file in .aux")
    if len(net_files) != 1:
        raise RuntimeError("There should be a .nets file in .aux")
    if len(pl_files) != 1:
        raise RuntimeError("There should be a .pl file in .aux")
    if len(scl_files) != 1:
        raise RuntimeError("There should be a .scl file in .aux")
    return (
        filename,
        os.path.join(dirname, node_files[0]),
        os.path.join(dirname, net_files[0]),
        os.path.join(dirname, pl_files[0]),
        os.path.join(dirname, scl_files[0]),
    )


def _parse_num_line(line):
    tokens = line.split(":")
    if len(tokens) != 2:
        raise RuntimeError(f"Couldn't interpret <{line}> as <Key : Value>")
    return int(tokens[1].strip())


def _read_nodes(filename):
    nb_nodes = None
    nb_terminals = None
    nodes = []
    with _open_file(filename) as f:
        first_line_found = False
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line.startswith("#"):
                continue
            if line.startswith("UCLA") and not first_line_found:
                first_line_found = True
                continue
            if line.startswith("NumNodes"):
                assert nb_nodes is None
                nb_nodes = _parse_num_line(line)
                continue
            if line.startswith("NumTerminals"):
                assert nb_terminals is None
                nb_terminals = _parse_num_line(line)
                continue
            vals = line.split()
            fixed = False
            obstruction = True
            if "terminal" in vals[1:]:
                fixed = True
            name = vals[0]
            if len(vals) >= 3:
                width, height = vals[1:3]
                width = int(width)
                height = int(height)
            else:
                # Dummy placed cell
                width = 0
                height = 0
                obstruction = False
            nodes.append((name, width, height, fixed, obstruction))
    if nb_nodes is not None:
        assert len(nodes) == nb_nodes
    if nb_terminals is not None:
        assert len([n for n in nodes if n[3]]) == nb_terminals
    names = [n[0] for n in nodes]
    widths = [n[1] for n in nodes]
    heights = [n[2] for n in nodes]
    fixed = [n[3] for n in nodes]
    obstruction = [n[4] for n in nodes]
    return names, widths, heights, fixed, obstruction


def _read_nets(filename, cell_names, cell_widths, cell_heights, sort_entries=False):
    name_dir = dict((name, i) for i, name in enumerate(cell_names))
    nb_nets = None
    nb_pins = None
    net_degree = None
    nets = []
    with _open_file(filename) as f:
        first_line_found = False
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line.startswith("#"):
                continue
            if line.startswith("UCLA") and not first_line_found:
                first_line_found = True
                continue
            if line.startswith("NumNets"):
                assert nb_nets is None
                nb_nets = _parse_num_line(line)
                continue
            if line.startswith("NumPins"):
                assert nb_pins is None
                nb_pins = _parse_num_line(line)
                continue
            line = line.replace(":", " ")
            if line.startswith("NetDegree"):
                vals = line.split()
                assert 2 <= len(vals) <= 3
                net_degree = int(vals[1])
                if len(vals) == 3:
                    name = vals[2]
                else:
                    name = f"n{len(nets)}"
                nets.append((name, net_degree, []))
                continue
            vals = line.split()
            assert len(vals) == 4 or len(vals) == 2
            if len(vals) == 4:
                cell, direction, x, y = vals
                x = float(x)
                y = float(y)
            else:
                cell, direction = vals
                x = 0.0
                y = 0.0
            assert cell in name_dir
            nets[-1][2].append((name_dir[cell], direction, x, y))
    total_pins = 0
    for name, net_degree, pins in nets:
        if net_degree != len(pins):
            raise RuntimeError(
                f"Net degree for {name} is {len(pins)}; expected {net_degree}"
            )
        total_pins += net_degree
    if nb_nets is not None:
        assert len(nets) == nb_nets
    if nb_pins is not None:
        assert total_pins == nb_pins
    if sort_entries:
        # Sort so that same-size nets are contiguous
        nets.sort(key=lambda net: len(net[-1]))
    cell_x_offset = [0.5 * c for c in cell_widths]
    cell_y_offset = [0.5 * c for c in cell_heights]
    ret = []
    for name, _, pins in nets:
        cells = []
        pin_x = []
        pin_y = []
        for cell, _, x, y in pins:
            cells.append(cell)
            pin_x.append(int(round(cell_x_offset[cell] + x)))
            pin_y.append(int(round(cell_y_offset[cell] + y)))
        ret.append((name, cells, pin_x, pin_y))
    return ret


def _read_place(filename, cell_names):
    name_dir = dict((name, i) for i, name in enumerate(cell_names))
    cell_x = [0 for i in cell_names]
    cell_y = [0 for i in cell_names]
    cell_orient = [None for i in cell_names]

    with _open_file(filename) as f:
        first_line_found = False
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line.startswith("#"):
                continue
            if line.startswith("UCLA") and not first_line_found:
                first_line_found = True
                continue
            line = line.replace(":", " ")
            vals = line.split()
            assert len(vals) >= 3
            cell, x, y, orient = vals[:4]
            assert cell in name_dir
            cell_ind = name_dir[cell]
            cell_x[cell_ind] = int(x)
            cell_y[cell_ind] = int(y)
            if orient in CellOrientation.__members__:
                cell_orient[cell_ind] = CellOrientation.__members__[orient]
            else:
                raise RuntimeError(f"Unknown orientation encountered {orient}")
    return cell_x, cell_y, cell_orient


def _read_rows(filename):
    nb_rows = None
    rows = []
    with _open_file(filename) as f:
        lines = [l.strip() for l in f]
        for line in lines:
            if line.startswith("NumRows"):
                assert nb_rows is None
                nb_rows = _parse_num_line(line)
        row_descs = []
        in_row = False
        for line in lines:
            if line.startswith("CoreRow"):
                row_descs.append([])
                in_row = True
            elif line.startswith("End"):
                in_row = False
            elif in_row:
                row_descs[-1].extend(line.replace(":", " ").split())
        for desc in row_descs:
            min_x = None
            min_y = None
            width = None
            height = None
            site_width = 1
            orient = CellOrientation.N
            for i in range(1, len(desc)):
                if desc[i - 1].lower() == "coordinate":
                    min_y = int(desc[i])
                if desc[i - 1].lower() == "subroworigin":
                    min_x = int(desc[i])
                if desc[i - 1].lower() == "numsites":
                    width = int(desc[i])
                if desc[i - 1].lower() == "height":
                    height = int(desc[i])
                if desc[i - 1].lower() == "sitewidth":
                    site_width = int(desc[i])
                if desc[i - 1].lower() == "siteorient":
                    if desc[i] in CellOrientation.__members__:
                        orient = CellOrientation.__members__[desc[i]]

            width *= site_width
            assert min_x is not None
            assert min_y is not None
            assert width is not None
            assert height is not None
            r = Rectangle(min_x, min_x + width, min_y, min_y + height)
            rows.append(Row(r, orient))
        return rows


def _str2bool(v):
    """
    Parse boolean arguments with argparse, from
    https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        import argparse
        raise argparse.ArgumentTypeError('Boolean value expected.')


class Circuit(coloquinte_pybind.Circuit):
    """
    Representation of a circuit for the placement tool
    """

    def __init__(self, nb_cells):
        super(Circuit, self).__init__(nb_cells)
        self._filename = None
        self._cell_name = None
        self._net_name = None

    @staticmethod
    def read_ispd(filename, ignore_obstructions=False):
        """
        Read an ISPD benchmark

        :param filename: Path of the .aux file or its directory
        :param ignore_obstructions: Do not consider preplaced macros as obstacles
        :return: A new Circuit object
        """
        (
            aux_filename,
            node_filename,
            net_filename,
            pl_filename,
            scl_filename,
        ) = _read_aux(filename)
        cell_names, cell_widths, cell_heights, cell_fixed, cell_obstruction = _read_nodes(
            node_filename
        )
        nets = _read_nets(net_filename, cell_names, cell_widths, cell_heights)
        cell_x, cell_y, cell_orient = _read_place(pl_filename, cell_names)
        rows = _read_rows(scl_filename)

        ret = Circuit(len(cell_names))
        ret._filename = os.path.splitext(aux_filename)[0]

        # Setup cell properties
        ret._cell_name = cell_names
        ret.cell_width = cell_widths
        ret.cell_height = cell_heights
        ret.cell_is_fixed = cell_fixed
        ret.cell_is_obstruction = cell_obstruction
        if ignore_obstructions:
            # All fixed cells marked as not obstructions
            ret.cell_is_obstruction = [not f for f in cell_fixed]

        # Setup nets and pins
        ret._net_name = []
        for name, cells, pins_x, pins_y in nets:
            ret._net_name.append(name)
            ret.add_net(cells, pins_x, pins_y)

        # Setup initial cell placement
        ret.cell_x = cell_x
        ret.cell_y = cell_y
        ret.cell_orientation = cell_orient

        # Setup rows
        ret.rows = rows

        # Allow macros to have any orientation
        row_height = ret.row_height
        polarities = ret.cell_row_polarity
        for i in range(ret.nb_cells):
            if cell_heights[i] > 4 * row_height:
                polarities[i] = CellRowPolarity.ANY
            elif cell_heights[i] % row_height != 0:
                polarities[i] = CellRowPolarity.NW
            elif cell_heights[i] % row_height == 0:
                polarities[i] = CellRowPolarity.SAME
        ret.cell_row_polarity = polarities

        ret.check()
        return ret

    def place_global(self, params, callback=None):
        """
        Run the global placement

        :param params: An integer effort (1-9) or the placement parameters
        :param callback: A callback to execute at each placement step
        """
        if not isinstance(params, ColoquinteParameters):
            if not isinstance(params, int):
                raise TypeError("Argument should be an integer effort")
            params = ColoquinteParameters(params)
        super().place_global(params, callback)

    def legalize(self, params, callback=None):
        """
        Run the legalization

        :param params: An integer effort (1-9) or the placement parameters
        :param callback: A callback to execute at each placement step
        """
        if not isinstance(params, ColoquinteParameters):
            if not isinstance(params, int):
                raise TypeError("Argument should be an integer effort")
            params = ColoquinteParameters(params)
        super().legalize(params, callback)

    def place_detailed(self, params, callback=None):
        """
        Run the detailed placement

        :param params: An integer effort (1-9) or the placement parameters
        :param callback: A callback to execute at each placement step
        """
        if not isinstance(params, ColoquinteParameters):
            if not isinstance(params, int):
                raise TypeError("Argument should be an integer effort")
            params = ColoquinteParameters(params)
        super().place_detailed(params, callback)

    def load_placement(self, filename):
        """
        Load a placement solution from an ISPD placement file

        :param filename: Path of the file
        """
        if filename is None:
            return

        cell_x, cell_y, cell_orient = _read_place(filename, self._cell_name)
        self.cell_x = cell_x
        self.cell_y = cell_y
        self.cell_orientation = cell_orient

    def write_placement(self, filename):
        """
        Write the placement result in ISPD file format

        :param filename: Path of the file
        """
        if filename is None:
            if self._filename is None:
                raise RuntimeError("No filename to export placement to")
            filename = os.path.splitext(self._filename)[0] + ".sol.pl"
        with _open_file(filename, True) as f:
            print("UCLA pl 1.0", file=f)
            print("# Created by Coloquinte", file=f)
            print("# https://github.com/Coloquinte/PlaceRoute", file=f)
            print("", file=f)
            cell_x = self.cell_x
            cell_y = self.cell_y
            cell_orientation = self.cell_orientation
            cell_is_fixed = self.cell_is_fixed
            for i in range(self.nb_cells):
                name = self._cell_name[i]
                x = cell_x[i]
                y = cell_y[i]
                orient = cell_orientation[i].name
                if cell_is_fixed[i]:
                    orient += " /FIXED"
                print(f"{name}\t{x}\t{y}\t: {orient}", file=f)

    def write_image(self, filename, macros_only=False, image_width=2048):
        """
        Export an image of the circuit
        """
        img, scale_factor = self._make_image(image_width)
        self._draw_rows(img, scale_factor)
        self._draw_cells(img, True, scale_factor)
        if not macros_only:
            self._draw_cells(img, False, scale_factor)
        self._save_image(img, filename)

    def write_displacement(self, filename, pl1, pl2, image_width=2048):
        """
        Export an image representing the displacement between two solutions
        """
        img, scale_factor = self._make_image(image_width)
        self._draw_rows(img, scale_factor)
        self._draw_cells(img, True, scale_factor)
        self._draw_displacement(img, pl1, pl2, scale_factor)
        self._save_image(img, filename)

    def _save_image(self, img, filename):
        img.save(filename, lossless=True)

    def _draw_area(self):
        placement = self.cell_placement
        rows = self.rows

        min_x_c = min(pl.min_x for pl in placement)
        min_y_c = min(pl.min_y for pl in placement)
        max_x_c = max(pl.max_x for pl in placement)
        max_y_c = max(pl.max_y for pl in placement)

        min_x_r = min(r.min_x for r in rows)
        min_y_r = min(r.min_y for r in rows)
        max_x_r = max(r.max_x for r in rows)
        max_y_r = max(r.max_y for r in rows)

        return (min(min_x_c, min_x_r), min(min_y_c, min_y_r), max(max_x_c, max_x_r), max(max_y_c, max_y_r))

    def _scale_factor(self, image_width):
        min_x, min_y, max_x, max_y = self._draw_area()
        max_allowed_ratio = 8
        fact = max(1.0, (max_x - min_x) / image_width,
                   (max_y - min_y) / (max_allowed_ratio * image_width))
        return fact, math.ceil((max_y - min_y) / fact)

    def _make_image(self, image_width):
        from PIL import Image
        scale_factor, image_height = self._scale_factor(image_width)
        img = Image.new("RGB", (image_width, image_height), "lightgray")
        return img, scale_factor

    def _draw_rows(self, img, scale_factor):
        from PIL import ImageDraw
        min_x, min_y = self._draw_area()[:2]
        draw = ImageDraw.Draw(img)

        rows = self.rows
        for r in rows:
            xmn = round((r.min_x - min_x) / scale_factor)
            xmx = round((r.max_x - min_x) / scale_factor)
            ymn = round((r.min_y - min_y) / scale_factor)
            ymx = round((r.max_y - min_y) / scale_factor)
            rect = [(xmn, ymn), (xmx, ymx)]
            draw.rectangle(rect, fill="white")

    def _draw_cells(self, img, macros, scale_factor):
        from PIL import ImageDraw
        min_x, min_y, max_x, max_y = self._draw_area()
        draw = ImageDraw.Draw(img)

        placement = self.cell_placement
        orientation = self.cell_orientation
        fixed = self.cell_is_fixed
        row_height = self.row_height
        for i, pl in enumerate(placement):
            if fixed[i] != macros:
                continue
            xmn = round((pl.min_x - min_x) / scale_factor)
            xmx = round((pl.max_x - min_x) / scale_factor)
            ymn = round((pl.min_y - min_y) / scale_factor)
            ymx = round((pl.max_y - min_y) / scale_factor)
            rect = [(xmn, ymn), (xmx, ymx)]
            if xmn == xmx or ymn == ymx:
                continue
            if fixed[i]:
                # Fixed cells and macros
                fill = "gray"
                outline = "black"
                width = 1
            elif pl.height > 4 * row_height:
                # Movable macros
                fill = "aqua"
                outline = "mediumblue"
                width = 1
            else:
                # Nice standard cells (or close enough)
                fill = "blue"
                outline = "mediumblue"
                width = 1
            draw.rectangle(rect, fill=fill, outline=outline, width=width)

            # Show the orientation with a small square in the corner
            osize = round(row_height / 5 / scale_factor)
            osize = min(osize, xmx - xmn)
            osize = min(osize, ymx - ymn)
            if osize > 2:  # Only do it for 2px and more
                if orientation[i] in [CellOrientation.N, CellOrientation.FW]:
                    rect = [(xmn, ymn), (xmn+osize, ymn+osize)]
                elif orientation[i] in [CellOrientation.S, CellOrientation.FE]:
                    rect = [(xmx-osize, ymx-osize), (xmx, ymx)]
                elif orientation[i] in [CellOrientation.FS, CellOrientation.E]:
                    rect = [(xmn, ymx-osize), (xmn+osize, ymx)]
                else:
                    rect = [(xmx-osize, ymn), (xmx, ymn+osize)]
                draw.rectangle(rect, fill=outline,
                               outline=outline, width=width)

        return img

    def _draw_displacement(self, img, pl1, pl2, scale_factor):
        from PIL import ImageDraw
        min_x, min_y = self._draw_area()[:2]

        draw = ImageDraw.Draw(img)
        fixed = self.cell_is_fixed
        assert len(pl1) == len(pl2)
        for i in range(self.nb_cells):
            if not fixed[i]:
                x1 = (pl1[i].min_x + pl1[i].max_x) / 2
                x2 = (pl2[i].min_x + pl2[i].max_x) / 2
                y1 = (pl1[i].min_y + pl1[i].max_y) / 2
                y2 = (pl2[i].min_y + pl2[i].max_y) / 2
                x1 = round((x1 - min_x) / scale_factor)
                x2 = round((x2 - min_x) / scale_factor)
                y1 = round((y1 - min_y) / scale_factor)
                y2 = round((y2 - min_y) / scale_factor)
                draw.line([(x1, y1), (x2, y2)], fill="red", width=1)
                draw.arc([x1 - 1, y1 - 1, x1 + 1, y1 + 1],
                         0, 360, fill="black")


def _add_arguments(parser, obj, prefix):
    import argparse
    for name in obj.__dir__():
        if name.startswith("_"):
            continue
        if name in ["check", "seed"]:
            continue
        child = getattr(obj, name)
        arg_type = type(child)
        if arg_type in (int, float):
            parser.add_argument(
                "--" + ".".join(prefix + [name, ]),
                type=arg_type,
                metavar=name.upper(),
                help=argparse.SUPPRESS
            )
        elif arg_type is bool:
            parser.add_argument(
                "--" + ".".join(prefix + [name, ]),
                type=_str2bool,
                metavar=name.upper(),
                help=argparse.SUPPRESS
            )
        elif '__members__' in child.__dir__():
            parser.add_argument(
                "--" + ".".join(prefix + [name, ]),
                choices=list(arg_type.__members__.keys()),
                metavar=name.upper(),
                help=argparse.SUPPRESS
            )
        else:
            _add_arguments(parser, child, prefix + [name,])


def _parse_arguments(args, obj, prefix):
    for name in obj.__dir__():
        if name.startswith("_"):
            continue
        if name in ["check", "seed"]:
            continue
        argname = ".".join(prefix + [name, ])
        child = getattr(obj, name)
        arg_type = type(child)
        if arg_type in (int, float, bool) or "__members__" in child.__dir__():
            val = getattr(args, argname)
            if val is not None:
                new_val = val
                old_val = getattr(obj, name)
                if arg_type not in (int, float, bool):
                    val = arg_type.__members__[val]
                    old_val = old_val.name
                print(
                    f"Overloading {argname} placement parameter "
                    f"({old_val} -> {new_val})"
                )
                setattr(obj, name, val)
        else:
            _parse_arguments(args, child, prefix + [name, ])


def _show_params(obj, tabs):
    for name in obj.__dir__():
        if name.startswith("_"):
            continue
        if name in ["check", "seed"]:
            continue
        child = getattr(obj, name)
        arg_type = type(child)
        if arg_type in (int, float, bool) or "__members__" in child.__dir__():
            if arg_type not in (int, float, bool):
                child = child.name
            p = "\t" * tabs
            print(f"{p}{name}: {child}")
    for name in obj.__dir__():
        if name.startswith("_"):
            continue
        if name in ["check", "seed"]:
            continue
        child = getattr(obj, name)
        arg_type = type(child)
        if arg_type in (int, float, bool) or "__members__" in child.__dir__():
            continue
        p = "\t" * tabs
        print(f"{p}{name}: ")
        _show_params(child, tabs + 1)


class OptimizationCallback:
    """
    Default callback for placement; used to save images, graphs and statistics
    """

    def __init__(self, circuit, prefix, image_width, extension):
        self.circuit = circuit
        self.step = 1
        self.prefix = prefix
        self.image_width = image_width
        self.extension = extension
        self.save_view = True
        self.save_displacement = True
        self.prev_placement = None
        self.history = []

    def __call__(self, step_name):
        if self.save_view:
            filename = f"{self.prefix}_{self.step:04d}_{step_name.name.lower()}.{self.extension}"
            self.circuit.write_image(filename, image_width=self.image_width)
            self._save_graph()
        if self.save_displacement:
            if self.prev_placement is not None:
                filename = f"{self.prefix}_{self.step:04d}_{step_name.name.lower()}_disp.{self.extension}"
                self.circuit.write_displacement(
                    filename, self.prev_placement, self.circuit.cell_placement, image_width=self.image_width)
            self.prev_placement = self.circuit.cell_placement

        self.history.append((self.step, step_name, self.circuit.hpwl()))
        self.step += 1

    def _save_graph(self):
        import matplotlib.pyplot as plt

        filename = f"{self.prefix}_WL.{self.extension}"
        plt.title("Wirelength over time")
        plt.xlabel("Step")
        plt.ylabel("Wirelength")
        for step_name in (
            PlacementStep.LowerBound,
            PlacementStep.UpperBound,
            PlacementStep.Detailed,
        ):
            steps = []
            vals = []
            for step, name, val in self.history:
                if name == step_name:
                    steps.append(step)
                    vals.append(val)
            if len(vals) > 0:
                plt.plot(steps, vals, label=step_name.name)
                plt.legend(loc="upper left")
        plt.savefig(filename)
        plt.clf()


def main():
    """
    Run the whole placement algorithm from the command line
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Place a benchmark circuit from the command line",
        usage="usage: coloquinte [-h] [--effort EFFORT] [--seed SEED] [--load-solution FILE] [--save-solution FILE] instance",
    )
    parser.add_argument("instance", help="Benchmark instance", nargs="?")
    parser.add_argument("--effort", help="Placement effort",
                        type=int, default=3)
    parser.add_argument(
        "--density", help="Expand the cells to reach the target density", type=float)
    parser.add_argument("--seed", help="Random seed", type=int, default=-1)
    parser.add_argument(
        "--load-solution", help="Load initial placement", metavar="FILE"
    )
    parser.add_argument("--save-solution",
                        help="Save final placement", metavar="FILE")
    parser.add_argument(
        "--show-parameters",
        help="Show available tuning parameters and their current value",
        action="store_true",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--report-only", help="Load the circuit and return", action="store_true")
    group.add_argument(
        "--no-global", help="Skip global placement", action="store_true")
    group.add_argument(
        "--only-global",
        help="Run only global placement (no legalization)",
        action="store_true",
    )
    group.add_argument(
        "--no-detailed", help="Skip detailed placement", action="store_true"
    )
    parser.add_argument(
        "--ignore-macros",
        help="Ignore fixed placement obstructions",
        action="store_true",
        dest="ignore_obstructions",
    )

    # Prefix to save images
    parser.add_argument("--save-images", help=argparse.SUPPRESS, type=str)
    # Save intermediate placement images
    parser.add_argument(
        "--save-all-images", help=argparse.SUPPRESS, action="store_true"
    )
    parser.add_argument(
        "--save-displacement", help=argparse.SUPPRESS, action="store_true"
    )
    parser.add_argument(
        "--save-graph", help=argparse.SUPPRESS, action="store_true")
    # Save intermediate placement images
    parser.add_argument(
        "--image-width", help=argparse.SUPPRESS, type=int, default=1080)
    # Save intermediate placement images
    parser.add_argument(
        "--image-extension", help=argparse.SUPPRESS, type=str, default="webp"
    )
    # Export the ISPD benchmark files after reading
    parser.add_argument(
        "--export-ispd", help=argparse.SUPPRESS, type=str
    )

    tuning_options = parser.add_argument_group("tuning options")
    _add_arguments(tuning_options, ColoquinteParameters(), [])
    args = parser.parse_args()

    print(f"Placement effort {args.effort}, seed {args.seed}")
    params = ColoquinteParameters(args.effort, args.seed)
    _parse_arguments(args, params, [])
    params.check()

    if args.show_parameters:
        print("Parameters can be set using command line options. For example, --detailed.nb_passes 2")
        print("Current parameter values:")
        _show_params(params, 1)
        return
    if args.instance is None:
        parser.print_help()
        return

    circuit = Circuit.read_ispd(args.instance, args.ignore_obstructions)
    if args.ignore_obstructions:
        print("Ignoring macros for standard cell placement")
    if args.load_solution is not None:
        print("Loading initial solution")
        circuit.load_placement(args.load_solution)
    if args.density is not None:
        if args.density <= 0.0 or args.density >= 1.0:
            raise RuntimeError(
                "Target density should be strictly between 0 and 1.")
        circuit.expand_cells_to_density(args.density, 0.0)
    print(circuit.report())

    sys.stdout.flush()
    if args.export_ispd is not None:
        print("Exporting as ISPD benchmark")
        circuit.export_ispd(args.export_ispd)
    callback = None
    if args.save_images is not None:
        circuit.write_image(
            f"{args.save_images}_macros.{args.image_extension}", True, args.image_width
        )
        if args.load_solution is not None:
            circuit.write_image(
                f"{args.save_images}_initial.{args.image_extension}", False, args.image_width
            )
        if args.save_all_images or args.save_graph:
            callback = OptimizationCallback(
                circuit, args.save_images, args.image_width, args.image_extension
            )
            callback.save_view = args.save_all_images
            callback.save_displacement = args.save_displacement
            if args.load_solution is not None:
                callback.prev_placement = circuit.cell_placement
    if args.report_only:
        return

    if args.no_global:
        print("Global placement skipped at user's request")
    else:
        circuit.place_global(params, callback)

    sys.stdout.flush()
    if args.only_global:
        print("Legalization and detailed placement skipped at user's request")
    elif args.no_detailed:
        circuit.legalize(params, callback)
        print("Detailed placement skipped at user's request")
    else:
        circuit.place_detailed(params, callback)
    circuit.write_placement(args.save_solution)

    if args.save_images is not None:
        circuit.write_image(
            f"{args.save_images}_placed.{args.image_extension}", False, args.image_width
        )
    if callback is not None:
        callback._save_graph()


__all__ = [
    "Circuit",
    "GlobalPlacerParameters",
    "LegalizationParameters",
    "DetailedPlacerParameters",
    "Rectangle",
    "LegalizationModel",
    "NetModel",
    "CellOrientation",
    "main",
]

if __name__ == "__main__":
    main()
