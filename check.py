#!/usr/bin/env python
"""
Check the calculations using nonuniform mesh and polygonal domains.

This time we use PETSc matrices and SLEPc eigensolver.
"""
from dolfin import Mesh, has_cgal, SubDomain, Point, FunctionSpace, Function, \
    TrialFunction, TestFunction, PETScMatrix, inner, dx, assemble, grad, \
    FacetFunction, DirichletBC, plot, Constant, SLEPcEigenSolver, Vector
import numpy as np
import sys

# CGAL or mshr
# use generate_mesh instead of Mesh for CSG
if not has_cgal():
    try:
        from mshr import Polygon, generate_mesh
        print "Using mshr, instead of CGAL"
    except:
        print "No CGAL and no mshr. Will be hard to generate meshes."
        exit(0)
else:
    from dolfin import Polygon
    # mshr uses this function to get mesh, while CGAL uses Mesh constructor
    generate_mesh = Mesh


try:
    gridsize = int(sys.argv[1])
except:
    gridsize = 64
try:
    ref = int(sys.argv[2])
except:
    ref = 0

FILE = open('results/domains_{}_{}_12.25.txt'.format(gridsize, ref), 'r')
save = open(
    'results/domains_{}_{}_12.25_checked.txt'.format(gridsize, ref), 'w')

if ref < 0:
    ref = 0

neumann = None
lengths = None


class NeumannBC (SubDomain):

    """ Neumann boundary conditons for the boundary of the triangle. """

    def inside(self, x, on):
        """ Sum of distances to endpoints should equal interval length. """
        if not on:
            return on
        vectors = neumann - np.array(x)
        errors = np.sum(np.sqrt(np.sum(vectors * vectors, axis=2)),
                        axis=1) - lengths
        return np.any(errors < 1E-14)


class OnBoundary (SubDomain):

    """ Mark all boundary DOFs. Later unmark few Neumann sides. """

    def inside(self, x, on):
        """ x unused, check if on boundary. """
        return on

ind = 0
for poly in FILE:
    try:
        poly = eval(poly)
        ind += 1
    except:
        continue
    print >>save, "\n\n ## Domain {} ## ".format(ind)
    print >>save, 'Polygon: ', poly
    vertices = poly[:-1] + poly[-1]
    # neumann edges
    neumann = poly[-2:-1] + poly[-1] + poly[0:1]
    neumann = [(p[0] / float(gridsize), 93 * p[1] / float(gridsize) / 128)
               for p in neumann]
    neumann = np.array([[p, q] for p, q in zip(neumann[:-1], neumann[1:])])
    lengths = [np.sqrt((p - q).dot(p - q)) for p, q in neumann]
    points = [Point(p[0] / float(gridsize), 93 * p[1] / float(gridsize) / 128)
              for p in vertices]
    try:
        mesh = generate_mesh(Polygon(points), 50*gridsize / 64*2**ref)
    except:
        mesh = generate_mesh(Polygon(points[::-1]), 50*gridsize / 64*2**ref)
    print >>save, "Mesh size: ", mesh.size(2)
    # diameter is smaller than this
    H = mesh.hmax()
    print >>save, 'Largest cell diameter: ', H

    #
    # solver
    #
    V = FunctionSpace(mesh, "CR", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    A = PETScMatrix()
    B = PETScMatrix()
    a = inner(grad(u), grad(v)) * dx
    b = u * v * dx
    assemble(a, tensor=A)
    assemble(b, tensor=B)

    # Neumann boundary marked with True
    nBC = NeumannBC()
    on = OnBoundary()
    boundary = FacetFunction("size_t", mesh, 1)
    on.mark(boundary, 2)
    nBC.mark(boundary, 0)
    bc = DirichletBC(V, 0, boundary, 2)
    plt = plot(boundary)
    # apply BC
    vec = assemble(Constant(0) * v * dx)
    bc.zero_columns(A, vec)
    bc.apply(A)
    bc.zero_columns(B, vec)
    bc.zero(B)
    # this is the same as B applied to a vector of ones, hence gives the
    # diagonal of B for diagonal B
    diag = assemble(Constant(1) * v * dx).array()
    diag = 1 / diag
    # eigensolver
    eigensolver = SLEPcEigenSolver(A, B)
    eigensolver.parameters["spectrum"] = "smallest real"
    eigensolver.parameters["spectral_shift"] = 1E-10
    eigensolver.parameters["spectral_transform"] = "shift-and-invert"
    eigensolver.parameters["problem_type"] = "gen_hermitian"
    eigensolver.solve(2)
    # first eigenvalue
    print >>save, " ## First eigenvalue ## "
    eigv, _, eigf, _ = eigensolver.get_eigenpair(0)
    v = Vector()
    B.mult(eigf, v)
    print >>save, 'B norm of eigenvector: ', v.array().dot(eigf)
    u = Function(V)
    u.vector()[:] = eigf
    print >>save, 'L2 norm of eigenfunction: ', assemble(u * u * dx)
    v2 = Vector()
    A.mult(eigf, v2)
    res = v2.array() - eigv * v.array()
    resnorm = np.sqrt(res.dot(diag * res))
    print >>save, 'Algebraic residual: ', resnorm
    eig = (eigv - resnorm) / (1 + 0.1932 * (eigv - resnorm) * H ** 2)
    print >>save, 'Numerical eigenvalue: ', eigv
    print >>save, 'Lower bound for eigenvalue: ', eig
    # second eigenvalue
    print >>save, " ## Second eigenvalue ## "
    eigv, _, eigf, _ = eigensolver.get_eigenpair(1)
    v = Vector()
    B.mult(eigf, v)
    print >>save, 'B norm of eigenvector: ', v.array().dot(eigf)
    v2 = Vector()
    A.mult(eigf, v2)
    res = v2.array() - eigv * v.array()
    resnorm = np.sqrt(res.dot(diag * res))
    print >>save, 'Algebraic residual: ', resnorm
    eig = (eigv - resnorm) / (1 + 0.1932 * (eigv - resnorm) * H ** 2)
    print >>save, 'Numerical eigenvalue: ', eigv
    if eig < 12.25:
        print >>save, "Error: eigenvalue below threshold"
        save.close()
        print "Error: eigenvalue below threshold"
        exit(1)


save.close()
FILE.close()
