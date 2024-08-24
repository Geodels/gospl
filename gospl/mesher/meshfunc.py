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
    """Class for handling triangular meshes.

    .. inheritance-diagram:: MeshTri
    """

    def __init__(self):
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

    def get_vertex_mask(self, subdomain=None):
        if subdomain is None:
            # https://stackoverflow.com/a/42392791/353337
            return np.s_[:]
        if subdomain not in self.subdomains:
            self._mark_vertices(subdomain)
        return self.subdomains[subdomain]["vertices"]

    def get_edge_mask(self, subdomain=None):
        """Get faces which are fully in subdomain.
        """
        if subdomain is None:
            # https://stackoverflow.com/a/42392791/353337
            return np.s_[:]

        if subdomain not in self.subdomains:
            self._mark_vertices(subdomain)

        # A face is inside if all its edges are in.
        # An edge is inside if all its nodes are in.
        is_in = self.subdomains[subdomain]["vertices"][self.idx_hierarchy]
        # Take `all()` over the first index
        is_inside = np.all(is_in, axis=tuple(range(1)))

        if subdomain.is_boundary_only:
            # Filter for boundary
            is_inside = is_inside & self.is_boundary_edge

        return is_inside

    def get_face_mask(self, subdomain):
        """Get faces which are fully in subdomain.
        """
        if subdomain is None:
            # https://stackoverflow.com/a/42392791/353337
            return np.s_[:]

        if subdomain not in self.subdomains:
            self._mark_vertices(subdomain)

        # A face is inside if all its edges are in.
        # An edge is inside if all its nodes are in.
        is_in = self.subdomains[subdomain]["vertices"][self.idx_hierarchy]
        # Take `all()` over all axes except the last two (face_ids, cell_ids).
        n = len(is_in.shape)
        is_inside = np.all(is_in, axis=tuple(range(n - 2)))

        if subdomain.is_boundary_only:
            # Filter for boundary
            is_inside = is_inside & self.is_boundary_facet

        return is_inside

    def get_cell_mask(self, subdomain=None):
        if subdomain is None:
            # https://stackoverflow.com/a/42392791/353337
            return np.s_[:]

        if subdomain.is_boundary_only:
            # There are no boundary cells
            return np.array([])

        if subdomain not in self.subdomains:
            self._mark_vertices(subdomain)

        is_in = self.subdomains[subdomain]["vertices"][self.idx_hierarchy]
        # Take `all()` over all axes except the last one (cell_ids).
        n = len(is_in.shape)
        return np.all(is_in, axis=tuple(range(n - 1)))

    def _mark_vertices(self, subdomain):
        """Mark faces/edges which are fully in subdomain.
        """
        if subdomain is None:
            is_inside = np.ones(len(self.node_coords), dtype=bool)
        else:
            is_inside = subdomain.is_inside(self.node_coords.T).T

            if subdomain.is_boundary_only:
                # Filter boundary
                self.mark_boundary()
                is_inside = is_inside & self.is_boundary_node

        self.subdomains[subdomain] = {"vertices": is_inside}
        return
    
    def update_values(self):
        if self.half_edge_coords is not None:
            # Constructing the temporary arrays
            # self.node_coords[self.idx_hierarchy] can take quite a while here.
            self.half_edge_coords = (
                self.node_coords[self.idx_hierarchy[1]]
                - self.node_coords[self.idx_hierarchy[0]]
            )

        if self.ei_dot_ej is not None:
            self.ei_dot_ej = np.einsum(
                "ijk, ijk->ij",
                self.half_edge_coords[[1, 2, 0]],
                self.half_edge_coords[[2, 0, 1]],
            )

        if self.ei_dot_ei is not None:
            e = self.half_edge_coords
            self.ei_dot_ei = np.einsum("ijk, ijk->ij", e, e)

        if self.cell_volumes is not None or self.ce_ratios is not None:
            self.cell_volumes, self.ce_ratios = cpt_tri_areas_and_ce_ratios(
                self.ei_dot_ej
            )

        self._interior_edge_lengths = None
        self._cell_circumcenters = None
        self._interior_ce_ratios = None
        self._control_volumes = None
        self._cell_partitions = None
        self._cv_centroids = None
        self._surface_areas = None
        self._signed_cell_areas = None
        self._cell_centroids = None
        return

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

    def num_delaunay_violations(self):
        # Delaunay violations are present exactly on the interior edges where
        # the ce_ratio is negative. Count those.
        return np.sum(self.ce_ratios_per_interior_edge < 0.0)

    def flip_until_delaunay(self):
        # If all coedge/edge ratios are positive, all cells are Delaunay.
        if np.all(self.ce_ratios > 0):
            return False

        # If all _interior_ coedge/edge ratios are positive, all cells are
        # Delaunay.
        if self.is_boundary_edge is None:
            self.create_edges()
        self.mark_boundary()
        if np.all(self.ce_ratios[~self.is_boundary_edge] > 0):
            return False

        if self._edges_cells is None:
            self._compute_edges_cells()

        num_flip_steps = 0
        ce_ratios_per_interior_edge = self.ce_ratios_per_interior_edge
        while np.any(ce_ratios_per_interior_edge < 0.0):
            num_flip_steps += 1

            is_flip_interior_edge = ce_ratios_per_interior_edge < 0.0

            interior_edges_cells = self._edges_cells[2]
            adj_cells = interior_edges_cells[is_flip_interior_edge].T

            # Check if there are cells for which more than one edge needs to be
            # flipped. For those, only flip one edge, namely that with the
            # smaller (more negative) ce_ratio.
            cell_gids, num_flips_per_cell = np.unique(adj_cells,
                                                      return_counts=True)
            critical_cell_gids = cell_gids[num_flips_per_cell > 1]
            for cell_gid in critical_cell_gids:
                edge_gids = self.cells["edges"][cell_gid]
                num_adj_cells, edge_id = self._edge_gid_to_edge_list[
                    edge_gids].T
                edge_ids = edge_id[num_adj_cells == 2]
                k = np.argmin(ce_ratios_per_interior_edge[edge_ids])
                is_flip_interior_edge[edge_ids] = False
                is_flip_interior_edge[edge_ids[k]] = True

            self.flip_interior_edges(is_flip_interior_edge)
            ce_ratios_per_interior_edge = self.ce_ratios_per_interior_edge

        return num_flip_steps > 1

    def flip_interior_edges(self, is_flip_interior_edge):
        if self._edges_cells is None:
            self._compute_edges_cells()

        interior_edges_cells = self._edges_cells[2]
        adj_cells = interior_edges_cells[is_flip_interior_edge].T

        edge_gids = self._edge_to_edge_gid[2][is_flip_interior_edge]
        adj_cells = interior_edges_cells[is_flip_interior_edge].T

        # Get the local ids of the edge in the two adjacent cells.
        # Get all edges of the adjacent cells
        ec = self.cells["edges"][adj_cells]
        # Find where the edge sits.
        hits = ec == edge_gids[None, :, None]
        # Make sure that there is exactly one match per cell
        assert np.all(np.sum(hits, axis=2) == 1)
        # translate to lids
        idx = np.empty(hits.shape, dtype=int)
        idx[..., 0] = 0
        idx[..., 1] = 1
        idx[..., 2] = 2
        lids = idx[hits].reshape((2, -1))

        #        3                   3
        #        A                   A
        #       /|\                 / \
        #      / | \               /   \
        #     /  |  \             /  1  \
        #   0/ 0 |   \1   ==>   0/_______\1
        #    \   | 1 /           \       /
        #     \  |  /             \  0  /
        #      \ | /               \   /
        #       \|/                 \ /
        #        V                   V
        #        2                   2
        #
        verts = np.array(
            [
                self.cells["nodes"][adj_cells[0], lids[0]],
                self.cells["nodes"][adj_cells[1], lids[1]],
                self.cells["nodes"][adj_cells[0], (lids[0] + 1) % 3],
                self.cells["nodes"][adj_cells[0], (lids[0] + 2) % 3],
            ]
        )

        self.edges["nodes"][edge_gids] = np.sort(verts[[0, 1]].T, axis=1)
        # No need to touch self.is_boundary_edge,
        # self.is_boundary_edge_individual; we're only flipping interior edges.

        # Do the neighboring cells have equal orientation (both node sets
        # clockwise/counterclockwise?
        equal_orientation = (
            self.cells["nodes"][adj_cells[0], (lids[0] + 1) % 3]
            == self.cells["nodes"][adj_cells[1], (lids[1] + 2) % 3]
        )

        # Set new cells
        self.cells["nodes"][adj_cells[0]] = verts[[0, 1, 2]].T
        self.cells["nodes"][adj_cells[1]] = verts[[0, 1, 3]].T

        # Set up new cells->edges relationships.
        previous_edges = self.cells["edges"][adj_cells].copy()

        i0 = np.ones(equal_orientation.shape[0], dtype=int)
        i0[~equal_orientation] = 2
        i1 = np.ones(equal_orientation.shape[0], dtype=int)
        i1[equal_orientation] = 2

        self.cells["edges"][adj_cells[0]] = np.column_stack(
            [
                np.choose((lids[1] + i0) % 3, previous_edges[1].T),
                np.choose((lids[0] + 2) % 3, previous_edges[0].T),
                edge_gids,
            ]
        )
        self.cells["edges"][adj_cells[1]] = np.column_stack(
            [
                np.choose((lids[1] + i1) % 3, previous_edges[1].T),
                np.choose((lids[0] + 1) % 3, previous_edges[0].T),
                edge_gids,
            ]
        )

        # update is_boundary_edge
        for k in range(3):
            self.is_boundary_edge[
                k, adj_cells] = self.is_boundary_edge_individual[
                self.cells["edges"][adj_cells, k]
            ]

        # Update the edge->cells relationship. It doesn't change for the
        # edge that was flipped, but for two of the other edges.
        confs = [
            (0, 1, np.choose((lids[0] + 1) % 3, previous_edges[0].T)),
            (1, 0, np.choose((lids[1] + i0) % 3, previous_edges[1].T)),
        ]
        for conf in confs:
            c, d, edge_gids = conf
            num_adj_cells, edge_id = self._edge_gid_to_edge_list[edge_gids].T

            k1 = num_adj_cells == 1
            k2 = num_adj_cells == 2
            assert np.all(np.logical_xor(k1, k2))

            # outer boundary edges
            edge_id1 = edge_id[k1]
            assert np.all(self._edges_cells[1][edge_id1]
                          [:, 0] == adj_cells[c, k1])
            self._edges_cells[1][edge_id1, 0] = adj_cells[d, k1]

            # interior edges
            edge_id2 = edge_id[k2]
            is_column0 = self._edges_cells[2][edge_id2][:, 0] == adj_cells[c,
                                                                           k2]
            is_column1 = self._edges_cells[2][edge_id2][:, 1] == adj_cells[c,
                                                                           k2]
            assert np.all(np.logical_xor(is_column0, is_column1))
            #
            self._edges_cells[2][edge_id2[is_column0],
                                 0] = adj_cells[d, k2][is_column0]
            self._edges_cells[2][edge_id2[is_column1],
                                 1] = adj_cells[d, k2][is_column1]

        # Schedule the cell ids for updates.
        update_cell_ids = np.unique(adj_cells.T.flat)
        # Same for edge ids
        k, edge_gids = self._edge_gid_to_edge_list[
            self.cells["edges"][update_cell_ids].flat
        ].T
        update_interior_edge_ids = np.unique(edge_gids[k == 2])

        self._update_cell_values(update_cell_ids, update_interior_edge_ids)
        return

    def _update_cell_values(self, cell_ids, interior_edge_ids):
        """Updates all sorts of cell information for the given cell IDs.
        """
        # update idx_hierarchy
        nds = self.cells["nodes"][cell_ids].T
        self.idx_hierarchy[..., cell_ids] = nds[self.local_idx]

        # update self.half_edge_coords
        self.half_edge_coords[:, cell_ids, :] = np.moveaxis(
            self.node_coords[self.idx_hierarchy[1, ..., cell_ids]]
            - self.node_coords[self.idx_hierarchy[0, ..., cell_ids]],
            0,
            1,
        )

        # update self.ei_dot_ej
        self.ei_dot_ej[:, cell_ids] = np.einsum(
            "ijk, ijk->ij",
            self.half_edge_coords[[1, 2, 0]][:, cell_ids],
            self.half_edge_coords[[2, 0, 1]][:, cell_ids],
        )

        # update self.ei_dot_ei
        e = self.half_edge_coords[:, cell_ids]
        self.ei_dot_ei[:, cell_ids] = np.einsum("ijk, ijk->ij", e, e)

        # update cell_volumes, ce_ratios_per_half_edge
        cv, ce = cpt_tri_areas_and_ce_ratios(self.ei_dot_ej[:, cell_ids])
        self.cell_volumes[cell_ids] = cv
        self.ce_ratios[:, cell_ids] = ce

        if self._interior_ce_ratios is not None:
            self._interior_ce_ratios[interior_edge_ids] = 0.0
            edge_gids = self._edge_to_edge_gid[2][interior_edge_ids]
            adj_cells = self._edges_cells[2][interior_edge_ids]

            is0 = self.cells["edges"][adj_cells[:, 0]][:, 0] == edge_gids
            is1 = self.cells["edges"][adj_cells[:, 0]][:, 1] == edge_gids
            is2 = self.cells["edges"][adj_cells[:, 0]][:, 2] == edge_gids
            assert np.all(
                np.sum(np.column_stack([is0, is1, is2]), axis=1) == 1
            )
            #
            self._interior_ce_ratios[interior_edge_ids[is0]] += self.ce_ratios[
                0, adj_cells[is0, 0]
            ]
            self._interior_ce_ratios[interior_edge_ids[is1]] += self.ce_ratios[
                1, adj_cells[is1, 0]
            ]
            self._interior_ce_ratios[interior_edge_ids[is2]] += self.ce_ratios[
                2, adj_cells[is2, 0]
            ]

            is0 = self.cells["edges"][adj_cells[:, 1]][:, 0] == edge_gids
            is1 = self.cells["edges"][adj_cells[:, 1]][:, 1] == edge_gids
            is2 = self.cells["edges"][adj_cells[:, 1]][:, 2] == edge_gids
            assert np.all(
                np.sum(np.column_stack([is0, is1, is2]), axis=1) == 1
            )
            #
            self._interior_ce_ratios[interior_edge_ids[is0]] += self.ce_ratios[
                0, adj_cells[is0, 1]
            ]
            self._interior_ce_ratios[interior_edge_ids[is1]] += self.ce_ratios[
                1, adj_cells[is1, 1]
            ]
            self._interior_ce_ratios[interior_edge_ids[is2]] += self.ce_ratios[
                2, adj_cells[is2, 1]
            ]

        if self._signed_cell_areas is not None:
            # One could make p contiguous by adding a copy(), but that's not
            # really worth it here.
            p = self.node_coords[self.cells["nodes"][cell_ids]].T
            # <https://stackoverflow.com/q/50411583/353337>
            self._signed_cell_areas[cell_ids] = (
                +p[0][2] * (p[1][0] - p[1][1])
                + p[0][0] * (p[1][1] - p[1][2])
                + p[0][1] * (p[1][2] - p[1][0])
            ) / 2

        # TODO update those values
        self._cell_centroids = None
        self._edge_lengths = None
        self._cell_circumcenters = None
        self._control_volumes = None
        self._cell_partitions = None
        self._cv_centroids = None
        self._surface_areas = None
        self.subdomains = {}
        return
