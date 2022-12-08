import numpy as np
from scipy.linalg import solve

class Solver():
    """
    A class to represent a numerical solver.
    
    ...

    Attributes
    ----------
    None

    Methods
    -------
    solveApproximateMatricesSadiku(A, ui, b, target):
        solves given matrices
    solveApproximateMatricesScipy(A, b):
        solves given matrices
    solveApproximateMatricesNumpy(A, b):
        solves given matrices
    solveNumericalMatricesSadiku(A, ui, b, target):
        solves given matrices
    solveNumericalMatricesScipy(A, b):
        solves given matrices
    solveNumericalMatricesNumpy(A, b):
        solves given matrices
    solveApproximateBoundaryMatricesSadiku(A, ui, b, target):
        solves given matrices
    solveApproximateBoundaryMatricesScipy(A, b):
        solves given matrices
    solveApproximateBoundaryMatricesNumpy(A, b):
        solves given matrices
    solveNumericalMatricesApproxBoundSadiku(A, ui, b, target):
        solves given matrices
    solveNumericalMatricesApproxBoundScipy(A, b):
        solves given matrices
    solveNumericalMatricesApproxBoundNumpy(A, b):
        solves given matrices
    solveNumericalMatricesApproxIterBoundSadiku(A, ui, b, target):
        solves given matrices
    solveNumericalMatricesApproxIterBoundScipy(A, b):
        solves given matrices
    solveNumericalMatricesApproxIterBoundNumpy(A, b):
        solves given matrices
    solveIterBoundaryMatricesSadiku(A, ui, b, target):
        solves given matrices
    solveIterBoundaryMatricesScipy(A, b):
        solves given matrices
    solveIterBoundaryMatricesNumpy(A, b):
        solves given matrices
    solveminiScipy(A, b):
        solves given matrices
    solveminiNumpy(A, b):
        solves given matrices
    solveminiSadiku(A, ui, b, target):
        solves given matrices
    """
    
    def __init__(self):
        """
        Constructs all the necessary attributes for the Solver object.

        Parameters
        ----------
        None
        """

        pass
        

    def solveApproximateMatricesSadiku(self, A, ui, b, target):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        ui : list
            initial guess
        b : list
            forcing matrix
        target : float
            target tolerence
        Returns
        -------
        None
        """

        u = ui.copy()
        error = 1

        while (error > target):
            up = u.copy()
            for i in range(len(u)):
                Az = A[i,:].copy()
                Az[i] = 0
                u[i] = (np.matmul(-Az,up)+b[i])/A[i,i]
                
            error = np.amax(np.absolute((u-up)/up*100))

        self.ua = u
        self.u = u

    def solveApproximateMatricesScipy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.ua = solve(A, b)
        self.u = self.ua.copy()

    def solveApproximateMatricesNumpy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.ua = np.linalg.lstsq(A, b)[0]
        self.u = self.ua.copy()

    def solveNumericalMatricesSadiku(self, A, ui, b, target):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        ui : list
            initial guess
        b : list
            forcing matrix
        target : float
            target tolerence
        Returns
        -------
        None
        """

        u = ui.copy()
        error = 1

        while (error > target):
            up = u.copy()
            for i in range(len(u)):
                Az = A[i,:].copy()
                Az[i] = 0
                u[i] = (np.matmul(-Az,up)+b[i])/A[i,i]
                
            error = np.amax(np.absolute((u-up)/up*100))

        self.un = u
        self.u = u

    def solveNumericalMatricesScipy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.un = solve(A, b)
        self.u = self.un.copy()

    def solveNumericalMatricesNumpy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.un = np.linalg.lstsq(A, b)[0]
        self.u = self.un.copy()

    def solveApproximateBoundaryMatricesSadiku(self, A, ui, b, target):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        ui : list
            initial guess
        b : list
            forcing matrix
        target : float
            target tolerence
        Returns
        -------
        None
        """

        u = ui.copy()
        error = 1

        while (error > target):
            up = u.copy()
            for i in range(len(u)):
                Az = A[i,:].copy()
                Az[i] = 0
                u[i] = (np.matmul(-Az,up)+b[i])/A[i,i]
                
            error = np.amax(np.absolute((u-up)/up*100))

        self.uab = u
        self.u = u

    def solveApproximateBoundaryMatricesScipy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.uab = solve(A, b)
        self.u = self.uab.copy()

    def solveApproximateBoundaryMatricesNumpy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.uab = np.linalg.lstsq(A, b)[0]
        self.u = self.uab.copy()

    def solveNumericalMatricesApproxBoundSadiku(self, A, ui, b, target):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        ui : list
            initial guess
        b : list
            forcing matrix
        target : float
            target tolerence
        Returns
        -------
        None
        """

        u = ui.copy()
        error = 1

        while (error > target):
            up = u.copy()
            for i in range(len(u)):
                Az = A[i,:].copy()
                Az[i] = 0
                u[i] = (np.matmul(-Az,up)+b[i])/A[i,i]
                
            error = np.amax(np.absolute((u-up)/up*100))

        self.unab = u
        self.u = u

    def solveNumericalMatricesApproxBoundScipy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.unab = solve(A, b)
        self.u = self.unab.copy()

    def solveNumericalMatricesApproxBoundNumpy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.unab = np.linalg.lstsq(A, b)[0]
        self.u = self.unab.copy()

    def solveNumericalMatricesApproxIterBoundSadiku(self, A, ui, b, target):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        ui : list
            initial guess
        b : list
            forcing matrix
        target : float
            target tolerence
        Returns
        -------
        None
        """

        u = ui.copy()
        error = 1

        while (error > target):
            up = u.copy()
            for i in range(len(u)):
                Az = A[i,:].copy()
                Az[i] = 0
                u[i] = (np.matmul(-Az,up)+b[i])/A[i,i]
                
            error = np.amax(np.absolute((u-up)/up*100))

        self.unib = u
        self.u = u

    def solveNumericalMatricesApproxIterBoundScipy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.unib = solve(A, b)
        self.u = self.unib.copy()

    def solveNumericalMatricesApproxIterBoundNumpy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.unib = np.linalg.lstsq(A, b)[0]
        self.u = self.unib.copy()

    def solveIterBoundaryMatricesSadiku(self, A, ui, b, target):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        ui : list
            initial guess
        b : list
            forcing matrix
        target : float
            target tolerence
        Returns
        -------
        None
        """

        u = ui.copy()
        error = 1

        while (error > target):
            up = u.copy()
            for i in range(len(u)):
                Az = A[i,:].copy()
                Az[i] = 0
                u[i] = (np.matmul(-Az,up)+b[i])/A[i,i]
                
            error = np.amax(np.absolute((u-up)/up*100))

        self.uib = u
        self.u = u

    def solveIterBoundaryMatricesScipy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.uib = solve(A, b)
        self.u = self.uib.copy()

    def solveIterBoundaryMatricesNumpy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.uib = np.linalg.lstsq(A, b)[0]
        self.u = self.uib.copy()
    
    def solveminiScipy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """

        self.miniu = solve(A, b)

    def solveminiNumpy(self, A, b):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        b : list
            forcing matrix
        Returns
        -------
        None
        """
        
        self.miniu = np.linalg.lstsq(A, b)[0]
        
    def solveminiSadiku(self, A, ui, b, target):
        """
        calculates solution of system of equations represented in matrix form.

        Parameters
        ----------
        A : list
            coeffecient matrix
        ui : list
            initial guess
        b : list
            forcing matrix
        target : float
            target tolerence
        Returns
        -------
        None
        """

        u = ui.copy()
        error = 1

        while (error > target):
            up = u.copy()
            for i in range(len(u)):
                Az = A[i,:].copy()
                Az[i] = 0
                u[i] = (np.matmul(-Az,up)+b[i])/A[i,i]
                
            error = np.amax(np.absolute((u-up)/up*100))

        self.miniu = u