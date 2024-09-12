import numpy as np


def grp_start_len(a):
    """Given a sorted 1D input array `a`, e.g., [0 0, 1, 2, 3, 4, 4, 4], this
    routine returns the indices where the blocks of equal integers start and
    how long the blocks are.
    """
    # https://stackoverflow.com/a/50394587/353337
    m = np.concatenate([[True], a[:-1] != a[1:], [True]])
    idx = np.flatnonzero(m)
    return idx[:-1], np.diff(idx)


def _column_stack(a, b):
    # https://stackoverflow.com/a/39638773/353337
    return np.stack([a, b], axis=1)


def unique_rows(a):
    # The cleaner alternative `np.unique(a, axis=0)` is slow; cf.
    # <https://github.com/numpy/numpy/issues/11136>.
    b = np.ascontiguousarray(a).view(
        np.dtype((np.void, a.dtype.itemsize * a.shape[1]))
    )
    a_unique, inv, cts = np.unique(b, return_inverse=True, return_counts=True)
    a_unique = a_unique.view(a.dtype).reshape(-1, a.shape[1])
    return a_unique, inv, cts


def cpt_tri_areas_and_ce_ratios(ei_dot_ej):
    """Given triangles (specified by their edges), this routine will return the
    triangle areas and the signed distances of the triangle circumcenters to
    the edge midpoints.
    """
    cell_volumes = 0.5 * np.sqrt(
        ei_dot_ej[2] * ei_dot_ej[0]
        + ei_dot_ej[0] * ei_dot_ej[1]
        + ei_dot_ej[1] * ei_dot_ej[2]
    )

    ce = -ei_dot_ej * 0.25 / cell_volumes[None]

    return cell_volumes, ce


def compute_triangle_circumcenters(X, ei_dot_ei, ei_dot_ej):
    """
    Computes the circumcenters of all given triangles.
    """
    alpha = ei_dot_ei * ei_dot_ej
    alpha_sum = alpha[0] + alpha[1] + alpha[2]
    beta = alpha / alpha_sum[None]
    a = X * beta[..., None]
    cc = a[0] + a[1] + a[2]

    return cc


class VoroBuild(object):
    """
    Class for handling triangular meshes.
    """

    def __init__(self):
        # nothing to do here

        return

    def __del__(self):
        # nothing to do here

        return

    def initVoronoi(self, nodes, cells, sort_cells=False):
        """Initialization.
        """
        if sort_cells:
            cells = np.sort(cells, axis=1)
            cells = cells[cells[:, 0].argsort()]

        self.node_coords = nodes
        self._edge_lengths = None

        assert (
            len(nodes.shape) == 2
        ), "Illegal node coordinates shape {}".format(nodes.shape)
        assert (
            len(cells.shape) == 2 and cells.shape[1] == 3
        ), "Illegal cells shape {}".format(cells.shape)

        # super(VoroBuild, self).__init__(nodes, cells)

        self.node_is_used = np.zeros(len(nodes), dtype=bool)
        self.node_is_used[cells] = True

        self.cells = {"nodes": cells}

        self._interior_ce_ratios = None
        self._control_volumes = None
        self._cell_partitions = None
        self._cv_centroids = None
        self._surface_areas = None
        self.edges = None
        self._cell_circumcenters = None
        self._signed_cell_areas = None
        self.subdomains = {}
        self._is_interior_node = None
        self._is_boundary_node = None
        self.is_boundary_edge = None
        self._is_boundary_facet = None
        self._interior_edge_lengths = None
        self._ce_ratios = None
        self._edges_cells = None
        self._edge_gid_to_edge_list = None
        self._edge_to_edge_gid = None
        self._cell_centroids = None

        self.local_idx = np.array([[1, 2], [2, 0], [0, 1]]).T
        nds = self.cells["nodes"].T
        self.idx_hierarchy = nds[self.local_idx]

        # The inverted local index.
        # This array specifies for each of the three nodes which edge endpoints
        # correspond to it. For the above local_idx, this should give
        #
        #    [[(1, 1), (0, 2)], [(0, 0), (1, 2)], [(1, 0), (0, 1)]]
        #
        self.local_idx_inv = [
            [tuple(i) for i in zip(*np.where(self.local_idx == k))]
            for k in range(3)
        ]

        # Create the corresponding edge coordinates.
        self.half_edge_coords = (
            self.node_coords[self.idx_hierarchy[1]]
            - self.node_coords[self.idx_hierarchy[0]]
        )

        self.ei_dot_ej = np.einsum(
            "ijk, ijk->ij",
            self.half_edge_coords[[1, 2, 0]],
            self.half_edge_coords[[2, 0, 1]],
        )

        e = self.half_edge_coords
        self.ei_dot_ei = np.einsum("ijk, ijk->ij", e, e)

        self.cell_volumes, self.ce_ratios = cpt_tri_areas_and_ce_ratios(
            self.ei_dot_ej
        )

        return

    @property
    def edge_lengths(self):
        if self._edge_lengths is None:
            self._edge_lengths = np.sqrt(self.ei_dot_ei)
        return self._edge_lengths

    @property
    def ce_ratios_per_interior_edge(self):
        if self._interior_ce_ratios is None:
            if "edges" not in self.cells:
                self.create_edges()

            self._ce_ratios = np.zeros(len(self.edges["nodes"]))
            np.add.at(self._ce_ratios, self.cells["edges"].T, self.ce_ratios)
            self._interior_ce_ratios = self._ce_ratios[
                ~self.is_boundary_edge_individual
            ]

        return self._interior_ce_ratios

    @property
    def control_volumes(self):
        if self._control_volumes is None:
            v = self.cell_partitions

            # Summing up the arrays first makes the work for np.add.at
            # lighter.
            ids = self.cells["nodes"].T
            vals = np.array(
                [
                    sum([v[i] for i in np.where(self.local_idx.T == k)[0]])
                    for k in range(3)
                ]
            )
            control_volume_data = [(ids, vals)]

            # sum up from self.control_volume_data
            self._control_volumes = np.zeros(len(self.node_coords))
            for d in control_volume_data:
                # TODO fastfunc
                np.add.at(self._control_volumes, d[0], d[1])

        return self._control_volumes

    @property
    def surface_areas(self):
        if self._surface_areas is None:
            self._surface_areas = self._compute_surface_areas()
        return self._surface_areas

    @property
    def control_volume_centroids(self):
        # This function is necessary, e.g., for Lloyd's
        # smoothing <https://en.wikipedia.org/wiki/Lloyd%27s_algorithm>.
        #
        # The centroid of any volume V is given by
        #
        #   c = \int_V x / \int_V 1.
        #
        # The denominator is the control volume. The numerator can be computed
        # by making use of the fact that the control volume around any vertex
        # v_0 is composed of right triangles, two for each adjacent cell.
        if self._cv_centroids is None:
            _, v = self._compute_integral_x()
            # Again, make use of the fact that edge k is opposite of node k in
            # every cell. Adding the arrays first makes the work for
            # np.add.at lighter.
            ids = self.cells["nodes"].T
            vals = np.array(
                [v[1, 1] + v[0, 2], v[1, 2] + v[0, 0], v[1, 0] + v[0, 1]]
            )
            centroid_data = [(ids, vals)]
            # add it all up
            num_components = centroid_data[0][1].shape[-1]
            self._cv_centroids = np.zeros((len(self.node_coords),
                                           num_components))
            for d in centroid_data:
                # TODO fastfunc
                np.add.at(self._cv_centroids, d[0], d[1])
            # Divide by the control volume
            self._cv_centroids /= self.control_volumes[:, None]

        return self._cv_centroids

    @property
    def signed_cell_areas(self):
        """Signed area of a triangle in 2D.
        """
        # http://mathworld.wolfram.com/TriangleArea.html
        assert (
            self.node_coords.shape[1] == 2
        ), "Signed areas only make sense for triangles in 2D."

        if self._signed_cell_areas is None:
            # One could make p contiguous by adding a copy(), but that's not
            # really worth it here.
            p = self.node_coords[self.cells["nodes"]].T
            # <https://stackoverflow.com/q/50411583/353337>
            self._signed_cell_areas = (
                +p[0][2] * (p[1][0] - p[1][1])
                + p[0][0] * (p[1][1] - p[1][2])
                + p[0][1] * (p[1][2] - p[1][0])
            ) / 2
        return self._signed_cell_areas

    def mark_boundary(self):
        if self.edges is None:
            self.create_edges()

        assert self.is_boundary_edge is not None

        self._is_boundary_node = np.zeros(len(self.node_coords), dtype=bool)
        self._is_boundary_node[
            self.idx_hierarchy[..., self.is_boundary_edge]] = True

        self._is_interior_node = self.node_is_used & ~self.is_boundary_node

        self._is_boundary_facet = self.is_boundary_edge
        return

    @property
    def is_boundary_node(self):
        if self._is_boundary_node is None:
            self.mark_boundary()
        return self._is_boundary_node

    @property
    def is_interior_node(self):
        if self._is_interior_node is None:
            self.mark_boundary()
        return self._is_interior_node

    @property
    def is_boundary_facet(self):
        if self._is_boundary_facet is None:
            self.mark_boundary()
        return self._is_boundary_facet

    def create_edges(self):
        """Set up edge-node and edge-cell relations.
        """
        # Reshape into individual edges.
        # Sort the columns to make it possible for `unique()` to identify
        # individual edges.
        s = self.idx_hierarchy.shape
        a = np.sort(self.idx_hierarchy.reshape(s[0], -1).T)
        a_unique, inv, cts = unique_rows(a)

        assert np.all(
            cts < 3
        ), "No edge has more than 2 cells. Are cells listed twice?"

        self.is_boundary_edge = (cts[inv] == 1).reshape(s[1:])

        self.is_boundary_edge_individual = cts == 1

        self.edges = {"nodes": a_unique}

        # cell->edges relationship
        self.cells["edges"] = inv.reshape(3, -1).T

        self._edges_cells = None
        self._edge_gid_to_edge_list = None

        # Store an index {boundary,interior}_edge -> edge_gid
        self._edge_to_edge_gid = [
            [],
            np.where(self.is_boundary_edge_individual)[0],
            np.where(~self.is_boundary_edge_individual)[0],
        ]
        return

    @property
    def edges_cells(self):
        if self._edges_cells is None:
            self._compute_edges_cells()
        return self._edges_cells

    def _compute_edges_cells(self):
        """This creates interior edge->cells relations. While it's not
        necessary for many applications, it sometimes does come in handy.
        """
        if self.edges is None:
            self.create_edges()

        num_edges = len(self.edges["nodes"])

        counts = np.zeros(num_edges, dtype=int)
        np.add.at(
            counts,
            self.cells["edges"],
            np.ones(self.cells["edges"].shape, dtype=int),
        )

        # <https://stackoverflow.com/a/50395231/353337>
        edges_flat = self.cells["edges"].flat
        idx_sort = np.argsort(edges_flat)
        idx_start, count = grp_start_len(edges_flat[idx_sort])
        res1 = idx_sort[idx_start[count == 1]][:, np.newaxis]
        idx = idx_start[count == 2]
        res2 = np.column_stack([idx_sort[idx], idx_sort[idx + 1]])
        self._edges_cells = [
            [],  # no edges with zero adjacent cells
            res1 // 3,
            res2 // 3,
        ]

        # For each edge, store the number of adjacent cells plus the index into
        # the respective edge array.
        self._edge_gid_to_edge_list = np.empty((num_edges, 2), dtype=int)
        self._edge_gid_to_edge_list[:, 0] = count
        c1 = count == 1
        l1 = np.sum(c1)
        self._edge_gid_to_edge_list[c1, 1] = np.arange(l1)
        c2 = count == 2
        l2 = np.sum(c2)
        self._edge_gid_to_edge_list[c2, 1] = np.arange(l2)
        assert l1 + l2 == len(count)

        return

    @property
    def edge_gid_to_edge_list(self):
        if self._edge_gid_to_edge_list is None:
            self._compute_edges_cells()
        return self._edge_gid_to_edge_list

    @property
    def face_partitions(self):
        # face = edge for triangles.
        # The partition is simply along the center of the edge.
        edge_lengths = self.edge_lengths
        return np.array([0.5 * edge_lengths, 0.5 * edge_lengths])

    @property
    def cell_partitions(self):
        if self._cell_partitions is None:
            # Compute the control volumes. Note that
            #   0.5 * (0.5 * edge_length) * covolume
            # = 0.25 * edge_length**2 * ce_ratio_edge_ratio
            self._cell_partitions = 0.25 * self.ei_dot_ei * self.ce_ratios
        return self._cell_partitions

    @property
    def cell_circumcenters(self):
        if self._cell_circumcenters is None:
            node_cells = self.cells["nodes"].T
            self._cell_circumcenters = compute_triangle_circumcenters(
                self.node_coords[node_cells], self.ei_dot_ei, self.ei_dot_ej
            )
        return self._cell_circumcenters

    @property
    def cell_centroids(self):
        """Computes the centroids (barycenters) of all triangles.
        """
        if self._cell_centroids is None:
            self._cell_centroids = (
                np.sum(self.node_coords[self.cells["nodes"]], axis=1) / 3.0
            )
        return self._cell_centroids

    @property
    def cell_barycenters(self):
        return self.cell_centroids

    @property
    def inradius(self):
        # See <http://mathworld.wolfram.com/Incircle.html>.
        abc = np.sqrt(self.ei_dot_ei)
        return 2 * self.cell_volumes / np.sum(abc, axis=0)

    @property
    def circumradius(self):
        # See <http://mathworld.wolfram.com/Incircle.html>.
        a, b, c = np.sqrt(self.ei_dot_ei)
        return (a * b * c) / np.sqrt(
            (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c)
        )

    @property
    def triangle_quality(self):
        # q = 2 * r_in / r_out
        #   = (-a+b+c) * (+a-b+c) * (+a+b-c) / (a*b*c),
        #
        # where r_in is the incircle radius and r_out the circumcircle radius
        # and a, b, c are the edge lengths.
        a, b, c = np.sqrt(self.ei_dot_ei)
        return (-a + b + c) * (a - b + c) * (a + b - c) / (a * b * c)

    @property
    def angles(self):
        # The cosines of the angles are the negative dot products of the
        # normalized edges adjacent to the angle.
        norms = np.sqrt(self.ei_dot_ei)
        normalized_ei_dot_ej = np.array(
            [
                self.ei_dot_ej[0] / norms[1] / norms[2],
                self.ei_dot_ej[1] / norms[2] / norms[0],
                self.ei_dot_ej[2] / norms[0] / norms[1],
            ]
        )
        return np.arccos(-normalized_ei_dot_ej)

    def _compute_integral_x(self):
        r"""Computes the integral of x,

          \int_V x,

        over all atomic "triangles", i.e., areas cornered by a node, an edge
        midpoint, and a circumcenter.
        """
        # The integral of any linear function over a triangle is the average of
        # the values of the function in each of the three corners, times the
        # area of the triangle.
        right_triangle_vols = self.cell_partitions

        node_edges = self.idx_hierarchy

        corner = self.node_coords[node_edges]
        edge_midpoints = 0.5 * (corner[0] + corner[1])
        cc = self.cell_circumcenters

        average = (corner + edge_midpoints[None] + cc[None, None]) / 3.0

        contribs = right_triangle_vols[None, :, :, None] * average

        return node_edges, contribs

    def _compute_surface_areas(self, cell_ids):
        """For each edge, one half of the the edge goes to each of the end
        points. Used for Neumann boundary conditions if on the boundary of the
        mesh and transition conditions if in the interior.
        """
        # Each of the three edges may contribute to the surface areas of all
        # three vertices. Here, only the two adjacent nodes receive a
        # contribution, but other approaches (e.g., the flat cell corrector),
        # may contribute to all three nodes.
        cn = self.cells["nodes"][cell_ids]
        ids = np.stack([cn, cn, cn], axis=1)

        half_el = 0.5 * self.edge_lengths[..., cell_ids]
        zero = np.zeros([half_el.shape[1]])
        vals = np.stack(
            [
                np.column_stack([zero, half_el[0], half_el[0]]),
                np.column_stack([half_el[1], zero, half_el[1]]),
                np.column_stack([half_el[2], half_el[2], zero]),
            ],
            axis=1,
        )

        return ids, vals

    def compute_curl(self, vector_field):
        r"""Computes the curl of a vector field over the mesh. While the vector
        field is point-based, the curl will be cell-based. The approximation is
        based on

        .. math::
            n\cdot curl(F) = \lim_{A\\to 0} |A|^{-1} <\int_{dGamma}, F> dr;

        see <https://en.wikipedia.org/wiki/Curl_(mathematics)>. Actually, to
        approximate the integral, one would only need the projection of the
        vector field onto the edges at the midpoint of the edges.
        """
        # Compute the projection of A on the edge at each edge midpoint.
        # Take the average of `vector_field` at the endpoints to get the
        # approximate value at the edge midpoint.
        A = 0.5 * np.sum(vector_field[self.idx_hierarchy], axis=0)
        # sum of <edge, A> for all three edges
        sum_edge_dot_A = np.einsum("ijk, ijk->j", self.half_edge_coords, A)

        # Get normalized vector orthogonal to triangle
        z = np.cross(self.half_edge_coords[0], self.half_edge_coords[1])

        # Now compute
        #
        #    curl = z / ||z|| * sum_edge_dot_A / |A|.
        #
        # Since ||z|| = 2*|A|, one can save a sqrt and do
        #
        #    curl = z * sum_edge_dot_A * 0.5 / |A|^2.
        #
        curl = z * (0.5 * sum_edge_dot_A / self.cell_volumes ** 2)[..., None]
        return curl

