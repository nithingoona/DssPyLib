import numpy as np

from . import Domain

class Point(Domain.Domain):
    def __init__(self, width, height, position, relativeMeshDensity, angle, isChild, isParent, childDomains, parentDomains, non):
        """
        Constructs all the necessary attributes for the geometry object.

        Parameters
        ----------
        shape: str
            the name of the shape of geometry
        """
        
        Domain.Domain.__init__(self, width, height, position, relativeMeshDensity, angle, isChild, isParent, childDomains, parentDomains, non)
        self.shape = "point"

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

        # Position of the point
        pos = np.array(self.position)
        # Plot the point
        plt.plot(pos[0], pos[1], 'k.', ms = lW*10)

    def calculateRelativeWidth(self, globalDomain):
        self.relativeWidth = 0

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
        self.verticesX = np.array(self.position[0])
        self.verticesY = np.array(self.position[1])
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

        self.boundaryNodesX = self.verticesX
        self.boundaryNodesY = self.verticesY
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
