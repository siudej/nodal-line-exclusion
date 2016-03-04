"""
Mixed Laplace eigenvalue solver for subdomains of a right triangle.

Based on:
 - FEniCS mesh generation and matrix assembly
 - scipy.sparse eigenvalue solver and boundary condition handling

Mesh and FEM matrices are created only ones. Then Dirichlet boundary
conditions are applied to interior DOFs to create a subdomain. We extract
slices from mass and stiffness matrices that correspond to non-Dirichlet DOFs.
"""
from dolfin import parameters, SubDomain, FunctionSpace, TrialFunction, \
    TestFunction, MeshFunction, inner, dx, uBLASSparseMatrix, assemble, grad
parameters.linear_algebra_backend = "uBLAS"

import scipy.sparse as sps
import scipy.sparse.linalg as ssl
import numpy as np


class Dirichlet(SubDomain):

    """ Subdomain for Dirichlet boundary from function. """

    def init(self, dirichlet):
        """ argument is function returning True if x on Dirichlet BC. """
        self.dirichlet = dirichlet

    def inside(self, x, on_boundary):
        """ on_boundary is not used, since we apply BC inside. """
        return self.dirichlet(x)


class Solver(object):

    """
    First order FEM with CR elements.

    Nonconforming CR elements give lower bounds for eigenvalues
    after apropriate postprocessing is applied.
    """

    def __init__(self, mesh, num, denom, method="CR"):
        """ Assemble matrices and cache matrices. """
        self.V = FunctionSpace(mesh, method, 1)
        # save coordinates of DOFs
        self.dofs = np.reshape(self.V.dofmap().tabulate_all_coordinates(mesh),
                               (-1, 2))
        u = TrialFunction(self.V)
        v = TestFunction(self.V)
        self.boundary = MeshFunction("size_t", mesh, 1)
        # assemble matrices
        a = inner(grad(u), grad(v)) * dx
        b = u * v * dx
        self.A = uBLASSparseMatrix()
        self.A = assemble(a, tensor=self.A)
        self.A.compress()
        self.B = uBLASSparseMatrix()
        self.B = assemble(b, tensor=self.B)
        self.B.compress()
        size = mesh.size(2)
        # cell diameter calculation
        self.H2 = (num ** 2 + denom ** 2) * 1.0 / size
        # print "Theoretical cell diameter: ", sqrt(self.H2)
        # print "FEniCS calculated: ", mesh.hmax()

        # Matrices have rational entries. We can rescale to get integers
        # and save as scipy matrices
        scaleB = 6.0*size/num/denom
        self.scaleB = scaleB
        self.B *= scaleB
        r, c, val = self.B.data()
        val = np.round(val)
        self.B = sps.csr_matrix((val, c, r))
        # check if B is diagonal
        assert len(self.B.diagonal()) == self.B.nnz
        # find B inverse
        self.Binv = sps.csr_matrix((1.0/val, c, r))
        scaleA = 1.0*num*denom/2
        self.scaleA = scaleA
        self.A *= scaleA
        r, c, val = self.A.data()
        self.A = sps.csr_matrix((np.round(val), c, r))

        self.scale = scaleA/scaleB
        print 'scaling A by: ', scaleA
        print 'scaling B by: ', scaleB
        print 'eigenvalues scale by:', self.scale

    def solve(self, dirichlet=lambda x: False, plotted=None):
        """
        Find eigenvalues given a boundary condition function.

        dirichlet is a function returning True if x is on Dirichlet BC.
        """
        self.boundary.set_all(0)
        d = Dirichlet()
        d.init(dirichlet)
        d.mark(self.boundary, 1)
        if plotted is not None:
            plotted.plot(self.boundary)
        #    plotted.write_png()

        # indices for non-Dirichlet rows/columns
        indices = np.nonzero(np.apply_along_axis(
            lambda x: not dirichlet(x), 1, self.dofs))[0]
        # remove Dirichlet rows and columns from A
        self.AA = (self.A[indices, :]).tocsc()[:, indices]
        # remove Dirichlet rows and columns from B
        self.BB = (self.B[indices, :]).tocsc()[:, indices]
        self.BBinv = (self.Binv[indices, :]).tocsc()[:, indices]
        # solve using scipy
        eigs, eigfs = ssl.eigsh(self.AA, k=2, M=self.BB, sigma=0, which='LM')
        self.raweigs = [eigs[0], eigs[1]]
        eig = eigs[0]
        # turn eigf into a function
        eigf = np.array(eigfs[:, 0]).flatten()
        # we were solving only for some dofs, rest is 0
        # u = Function(self.V)
        # u.vector()[:] = 0
        # u.vector()[indices] = eigf
        # find algebraic residual
        # both L2 norm of u and B norm of eigf should equal 1
        # print assemble(u*u*dx), self.BB.dot(eigf).dot(eigf)
        res = self.AA.dot(eigf) - eig * self.BB.dot(eigf)
        resnorm = np.sqrt(self.BBinv.dot(res).dot(res))
        self.residual = [resnorm, 0]
        # apply Carstensen-Gedicke transformations to eig
        #
        # kappa^2 less than 0.1932
        eig = (eig - resnorm) / (1 + 0.1932 * (eig - resnorm) *
                                 self.H2 / self.scale)
        # scale back the eigenvalue
        eig = eig/self.scale
        # find residual for the second eigenvalue (for gap calculations)
        eigf = np.array(eigfs[:, 1]).flatten()
        res = self.AA.dot(eigf) - eigs[1] * self.BB.dot(eigf)
        resnorm = np.sqrt(self.BBinv.dot(res).dot(res))
        self.residual[1] = resnorm
        # return (eig, u)  # pair (eigenvalue,eigenfunctions)
        return (eig, None)  # pair (eigenvalue,eigenfunctions)
