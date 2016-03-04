"""
Class Optimize implements the algorithm from the paper.

Almost everything can be exported to pdf/txt files.
Here we create the mesh and call FEniCS/scipy solver.
"""

import subprocess
from dolfin import Mesh, MeshEditor, refine, MeshFunction, plot
from solver2 import Solver
import numpy as np


def Triangle(left, bottom):
    """ Triangular mesh with just one triangular cell. """
    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, 2, 2)
    editor.init_vertices(3)
    editor.init_cells(1)
    editor.add_vertex(0, 0, 0)
    editor.add_vertex(1, bottom, 0)
    editor.add_vertex(2, 0, left)
    editor.add_cell(0, 0, 1, 2)
    editor.close()
    return mesh


def trunc_triangle(n, side):
    """ Build the complement of an excluded set (in grid coordinates). """
    if len(side) < n + 1:
        raise Exception("More points please")
    shape = [(0, side[0])]
    current_point = 0
    last_y = side[0]
    while current_point + side[current_point] < n:
        if side[current_point] > last_y:
            shape.append((current_point, last_y))
            shape.append((current_point, side[current_point]))
            last_y = side[current_point]
        current_point += 1
    if side[current_point] > last_y:
        shape.append((current_point, last_y))
    shape.append((current_point, n - current_point))
    shape.append((0, n))
    return shape


def augment(side, n, x, y):
    """ Enlarge excluded domain. """
    increment = [0] * x
    increment.extend([y] * (n - x + 1))
    result = [max(increment[u], side[u]) for u in range(len(side))]
    return result


def capture_command(command):
    """ Wrapper for command execution with output capturing. """
    arg_list = ["bash", "-lc", command]

    proc = subprocess.Popen(arg_list, stdout=subprocess.PIPE)
    out = proc.communicate()[0]
    return out


class Optimize:

    """ Stores all eigenproblem parameters and implements the algorithm. """

    def __init__(self, n=16, divide=1, threshold=12.3,
                 left_side_num=3, left_side_denom=4):
        """ Store parameters, initialize data exports (latex, txt). """
        # left_side_num/denom are either height, assuming base is 1,
        #   or left and bottom side lengths
        self.n = n
        self.threshold = threshold / left_side_denom ** 2
        self.left_side = [0] * (n + 1)
        self.bottom_side = [0] * (n + 1)
        self.left_side_len = left_side_num * 1.0 / left_side_denom
        self.left_side_num = left_side_num
        self.left_side_denom = left_side_denom
        self.folder = "results"
        self.name = "domains_{}_{}_{}".format(n, divide, threshold)
        self.filename = self.folder + "/" + self.name + ".tex"
        self.matrixZname = self.folder + "/matrices_Z/" + self.name
        self.matrixQname = self.folder + "/matrices_Q/" + self.name
        self.genericname = self.folder + "/" + self.name
        self.textname = self.folder + "/" + self.name + ".txt"
        f = open(self.textname, 'w')
        f.close()
        self.latex = "cd " + self.folder + "; pdflatex --interaction=batchmode" \
            + " {0}.tex; rm {0}.aux {0}.log".format(self.name)
        self.shift = 0
        # common mesh, ready to cut/apply Dirichlet BC
        self.mesh = Triangle(self.left_side_num, left_side_denom)
        while self.mesh.size(2) < self.n * self.n:
            self.mesh = refine(self.mesh)
        for i in range(divide):
            self.mesh = refine(self.mesh)
        boundary = MeshFunction("size_t", self.mesh, 1)
        boundary.set_all(0)
        self.plotted = plot(
            boundary,
            prefix='animation/animation',
            scalarbar=False,
            window_width=1024,
            window_height=1024)
        print 'Grid size: ', self.n ** 2
        print 'Mesh size: ', self.mesh.size(2)
        # setup solver
        self.solver = Solver(self.mesh, left_side_num, left_side_denom)
        self.mass = self.solver.B
        self.stiff = self.solver.A
        print np.unique(self.stiff.data)
        self.save_matrices(0)
        # save dofs
        f = open(self.genericname+'_dofmap.txt', 'w')
        for vec in self.solver.dofs:
            print >>f, vec
        f.close()

    def save_matrices(self, ind):
        """
        Export scipy matrices to files.

        Save both rational and integer forms. Rational matrices are based on
        a x b mesh instead of (a/b) x 1 mesh.
        """
        M = self.stiff.tocoo()
        f = open(self.matrixZname + '_stiff_{}.txt'.format(ind), 'w')
        print >>f, zip(M.row, M.col, M.data)
        f.close()
        f = open(self.matrixQname + '_stiff_{}.txt'.format(ind), 'w')
        print >>f, zip(M.row, M.col, M.data/self.solver.scaleA)
        f.close()
        M = self.mass.tocoo()
        f = open(self.matrixZname + '_mass_{}.txt'.format(ind), 'w')
        print >>f, zip(M.row, M.col, M.data)
        f.close()
        f = open(self.matrixQname + '_mass_{}.txt'.format(ind), 'w')
        print >>f, zip(M.row, M.col, M.data/self.solver.scaleB)
        f.close()

    def init_left(self, left_side=None):
        """ Initialize left exclusion. """
        if left_side is not None:
            self.left_side = [int(s) for s in left_side]
            l = len(self.left_side)
            if l < self.n + 1:
                self.left_side.extend([self.n - l] * (self.n - l + 1))
        else:
            self.left_side = [0] * (self.n + 1)
        self.find_corners(self.left_side)

    def init_bottom(self, bottom_side=None):
        """ Initialize bottom exclusion. """
        if bottom_side is not None:
            self.bottom_side = bottom_side
            l = len(self.bottom_side)
            if l < self.n + 1:
                self.bottom_side.extend([self.n - l] * (self.n - l + 1))
        else:
            self.bottom_side = [0] * (self.n + 1)
        self.find_corners(self.bottom_side)

    def upper_domain(self, bottom_side=None, raw=False):
        """
        Return upper polygon vertices.

        If raw, then in grid coordinates, otherwise scaled with side legths.
        """
        # bottom_side might be augmented
        if bottom_side is None:
            # otherwise use from object
            bottom_side = self.bottom_side
        raw_points = trunc_triangle(self.n, bottom_side)
        if raw:
            return raw_points
        output = []
        for (u, v) in raw_points:
            x = float(u - self.shift * 0.5) / self.n
            y = float(v - self.shift * (self.n - 1.0)) / \
                self.n * self.left_side_len
            output.extend([x, y])
        self.find_corners(bottom_side)
        return output

    def lower_domain(self, left_side=None, raw=False):
        """
        Return lower polygon vertices.

        If raw, then in grid coordinates, otherwise scaled with side legths.
        """
        # left_side might be augmented
        if left_side is None:
            # otherwise use from object
            left_side = self.left_side
        raw_points = trunc_triangle(self.n, left_side)
        if raw:
            return raw_points
        output = []
        for (u, v) in raw_points:
            x = float(v - self.shift * (self.n - 1)) / self.n
            y = float(u - self.shift * 0.5) / self.n * self.left_side_len
            output.extend([x, y])
        self.find_corners(left_side)
        return output

    def update_left_side(self, old):
        """ Search for the optimal rectangle to be added to left side. """
        (x, y) = (0, 1)
        best_score = 0  # best score so far
        for step in range(self.n + 1):
            if self.left_side[y] < x:
                test_side = augment(self.bottom_side, self.n, x, y)
                if self.left_side_test(test_side):
                    # save a few things from the solver in left_side_test
                    tempmass = self.solver.BB
                    tempstiff = self.solver.AA
                    tempraweigs = self.solver.raweigs
                    tempresidual = self.solver.residual
                    test_complement = augment(self.left_side, self.n, y, x)
                    self.find_corners(test_complement)
                    test_eig, _ = self.solver.solve(
                        lambda x: self.outside_left(x), plotted=self.plotted)
                    test_score = test_eig * self.left_side_denom ** 2
                    if test_score > best_score:
                        best_score = test_score
                        best_rect = (y, x)
                        # save solver results
                        self.mass = tempmass
                        self.stiff = tempstiff
                        self.raweigs = tempraweigs
                        self.lowerbound = self.templower*self.left_side_denom**2
                        self.residual = tempresidual
                        self.testraweigs = self.solver.raweigs
                        self.testresidual = self.solver.residual
                        self.teststiff = self.solver.AA
                        self.testmass = self.solver.BB
                    x += 1
                else:
                    y += 1
            else:  # success: alredy in boundary.
                x += 1
        return best_score, best_rect

    def left_side_test(self, test_side):
        """ Use solver to check if above threshold. """
        self.find_corners(test_side)
        self.templower, _ = self.solver.solve(lambda x: self.outside_bottom(x),
                                              plotted=self.plotted)
        return self.templower > self.threshold

    def update_bottom_side(self, old):
        """ Search for the optimal rectangle to be added to bottom side. """
        (x, y) = (1, 0)
        best_score = 0  # best score so far
        for step in range(self.n + 1):
            if self.bottom_side[x] < y:
                test_side = augment(self.left_side, self.n, y, x)
                if self.bottom_side_test(test_side):
                    # found good rectangle
                    # save a few things from the solver in left_side_test
                    tempmass = self.solver.BB
                    tempstiff = self.solver.AA
                    tempraweigs = self.solver.raweigs
                    tempresidual = self.solver.residual
                    test_complement = augment(self.bottom_side, self.n, x, y)
                    self.find_corners(test_complement)
                    test_eig, _ = self.solver.solve(
                        lambda x: self.outside_bottom(x), plotted=self.plotted)
                    test_score = test_eig * self.left_side_denom ** 2
                    if test_score > best_score:
                        best_score = test_score
                        best_rect = (x, y)
                        # update FEM matrices for the best rectangle
                        self.mass = tempmass
                        self.stiff = tempstiff
                        # update raw eigenvalues and residuals
                        self.raweigs = tempraweigs
                        self.lowerbound = self.templower*self.left_side_denom**2
                        self.residual = tempresidual
                        self.testraweigs = self.solver.raweigs
                        self.testresidual = self.solver.residual
                        self.teststiff = self.solver.AA
                        self.testmass = self.solver.BB
                    y += 1
                else:
                    x += 1
            else:  # success: alredy in boundary.
                y += 1
        return best_score, best_rect

    def bottom_side_test(self, test_side):
        """ Use solver to check if above threshold. """
        self.find_corners(test_side)
        self.templower, _ = self.solver.solve(lambda x: self.outside_left(x),
                                              plotted=self.plotted)
        return self.templower > self.threshold

    def find_corners(self, side):
        """ Find all corners on the boundary of the exclusion. """
        self.corners = [
            (i, side[i])
            for i in range(1, len(side)) if side[i] > side[i - 1]]

    def outside_bottom(self, p):
        """ Check if p is outside the domain (to mark Dirichlet condition). """
        for corner in self.corners:
            if self.n * p[0] >= self.left_side_denom * corner[0] and \
                    self.n * p[1] <= self.left_side_num * corner[1]:
                return True
        return p[1] <= 0

    def outside_left(self, p):
        """ Check if p is outside the domain (to mark Dirichlet condition). """
        for corner in self.corners:
            if self.n * p[0] <= self.left_side_denom * corner[1] and \
                    self.n * p[1] >= self.left_side_num * corner[0]:
                return True
        return False
