import numpy as np

from . import Domain

class Line(Domain.Domain):
    def __init__(self, width, height, position, relativeMeshDensity, angle, isChild, isParent, childDomains, parentDomains, non):
        """
        Constructs all the necessary attributes for the geometry object.

        Parameters
        ----------
        shape: str
            the name of the shape of geometry
        """
        
        Domain.Domain.__init__(self, width, height, position, relativeMeshDensity, angle, isChild, isParent, childDomains, parentDomains, non)
        self.shape = "line"

    def show(self, plt, lW):
        """
        adds the current geometry to plot

        Parameters
        ----------
        plt : object
            matplotlib object
        lw : float
            width of line stroke

        Returns
        -------
        None
        """

        # Generating boundary nodes for plotting
        w = self.width
        v0 = np.array([[0,0],[w,0]])
        v0 = np.transpose(v0)
        # Rotate the boundary points
        thr = self.angle
        R = np.array([[np.cos(np.radians(thr)), -np.sin(np.radians(thr))],
        [np.sin(np.radians(thr)), np.cos(np.radians(thr))]])
        vr = np.dot(R,v0)
        # Reposition the boundary points
        position = np.array(self.position)
        v = position +  np.transpose(vr)
        # Plot the boundary points
        plt.plot(v[[0,1],0], v[[0,1],1], 'k', lw = lW)

    def calculateArea(self):
        """
        calculates area.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.area = 0

    def calculateMaxFreeNodes(self, domains):
        """
        calculates maximum number of free nodes allowed

        Parameters
        ----------
        domains : list
            list of all geometry objects

        Returns
        -------
        None
        """

        self.maxFreeNodes = 0

    def genVertices(self, plt, lW, markerShape):
        """
        gets vertices of the geometry.

        Parameters
        ----------
        plt : object
            matplotlib object
        lw : float
            width of line stroke
        markerShape : char
            shape of marker

        Returns
        -------
        None
        """

        # Generating vertices
        w = self.width
        v0 = np.array([[0,0], [w,0]])
        v0 = np.transpose(v0)
        # Rotate the boundary points
        thr = self.angle
        R = np.array([[np.cos(np.radians(thr)), -np.sin(np.radians(thr))],
        [np.sin(np.radians(thr)), np.cos(np.radians(thr))]])
        vr = np.dot(R,v0)
        # Reposition the boundary points
        position = np.array(self.position)
        v = np.transpose(position +  np.transpose(vr))
        self.verticesX = v[0]
        self.verticesY = v[1]
        # Plot the vertices
        plt.scatter(self.verticesX, self.verticesY, c=None, marker=markerShape, s = lW*1000)

    def genBoundaryNodes(self, plt, lW, markerShape):
        """
        gets boundary nodes.

        Parameters
        ----------
        plt : object
            matplotlib object
        lw : float
            width of line stroke
        markerShape : char
            shape of marker

        Returns
        -------
        None
        """

        self.boundaryNodesX = np.linspace(self.verticesX[0],self.verticesX[1],self.widthNodes)
        self.boundaryNodesY = np.linspace(self.verticesY[0],self.verticesY[1],self.widthNodes)
        plt.scatter(self.boundaryNodesX, self.boundaryNodesY, c=None, marker=markerShape, s = lW*1000)

    def genRandFreeNodes(self, initialDomains, plt, lW, markerShape):
        """
        gets random free nodes inside a domain.

        Parameters
        ----------
        initialDomains : list
            list of all geometry objects
        plt : object
            matplotlib object
        lw : float
            width of line stroke
        markerShape : char
            shape of marker

        Returns
        -------
        None
        """

        self.freeNodesX = np.array([])
        self.freeNodesY = np.array([])
        self.allNodesX = self.boundaryNodesX
        self.allNodesY = self.boundaryNodesY

    def genTriFreeNodes(self, initialDomains, plt, lW, markerShape):
        """
        gets uniformly distributed triangular free nodes inside a domain.

        Parameters
        ----------
        initialDomains : list
            list of all geometry objects
        plt : object
            matplotlib object
        lw : float
            width of line stroke
        markerShape : char
            shape of marker

        Returns
        -------
        None
        """

        self.freeNodesX = np.array([])
        self.freeNodesY = np.array([])
        self.allNodesX = self.boundaryNodesX
        self.allNodesY = self.boundaryNodesY