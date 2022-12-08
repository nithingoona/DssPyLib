import numpy as np
import random

from . import Domain


class Rectangle(Domain.Domain):
    def __init__(self, width, height, position, relativeMeshDensity, angle, isChild, isParent, childDomains, parentDomains, non):
        """
        Constructs all the necessary attributes for the geometry object.

        Parameters
        ----------
        shape: str
            the name of the shape of geometry
        """

        Domain.Domain.__init__(self, width, height, position, relativeMeshDensity, angle, isChild, isParent, childDomains, parentDomains, non)
        self.shape = "rectangle"

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
        h = self.height
        v0 = np.array([[-w/2, -h/2], [w/2, -h/2], [w/2, h/2], [-w/2, h/2]])
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
        for i in range(3):
            plt.plot(v[[i,i+1],0], v[[i,i+1],1], 'k', lw = lW)
        plt.plot(v[[0,3],0], v[[0,3],1], 'k', lw = lW)

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

        self.area = self.width*self.height

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

        ab = domains[0].non
        a = domains[0].width
        b = domains[0].height
        mtn = round(1.2*(ab-2)**2*b/a)

        area = self.area
        if self.isParent:
            for childDomain in self.childDomains:
                area = area - domains[childDomain].area

        ga = domains[0].area
        rd = self.relativeMeshDensity[1]
        self.maxFreeNodes = round(mtn*area/ga*rd**2)

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
        h = self.height
        v0 = np.array([[-w/2, -h/2], [w/2, -h/2], [w/2, h/2], [-w/2, h/2]])
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

        nw = self.widthNodes
        nh = round(nw*self.height/self.width)
        xb = np.linspace(self.verticesX[0],self.verticesX[1],nw)
        yb = np.linspace(self.verticesY[0],self.verticesY[1],nw)
        xb = np.append(xb[0:-1], np.linspace(self.verticesX[1],self.verticesX[2],nh)[0:-1])
        yb = np.append(yb[0:-1], np.linspace(self.verticesY[1],self.verticesY[2],nh)[0:-1])
        xb = np.append(xb, np.linspace(self.verticesX[2],self.verticesX[3],nw)[0:-1])
        yb = np.append(yb, np.linspace(self.verticesY[2],self.verticesY[3],nw)[0:-1])
        xb = np.append(xb, np.linspace(self.verticesX[3],self.verticesX[0],nh)[0:-1])
        yb = np.append(yb, np.linspace(self.verticesY[3],self.verticesY[0],nh)[0:-1])

        self.boundaryNodesX = xb
        self.boundaryNodesY = yb
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

        hbyw = self.height/self.width
        mnf = self.maxFreeNodes
        thr = self.angle
        R = np.array([[np.cos(np.radians(thr)), -np.sin(np.radians(thr))],
        [np.sin(np.radians(thr)), np.cos(np.radians(thr))]])

        fix = self.boundaryNodesX
        fiy = self.boundaryNodesY
        h = self.width/self.widthNodes

        if self.isParent: # Add child nodes as fixed nodes
            for cd in self.childDomains:
                fix = np.append(fix, initialDomains[cd].allNodesX)
                fiy = np.append(fiy, initialDomains[cd].allNodesY)

        nop = 0 # Initializing no. of points inserted
        count = 0 # Count for while loop
        fx = np.array([0]) # Initializing free node positions
        fy = np.array([0])

        if hbyw >= 1:
            maxit = 1000*mnf*hbyw

            while nop < mnf and count < maxit:
                # Get random points
                xr = random.uniform(0,1) - 0.5
                yr = random.uniform(0,1) - 0.5

                if xr > -1/(2*hbyw) and xr < 1/(2*hbyw):
                    # Scale the random point
                    xr = xr*self.height
                    yr = yr*self.height
                    # Rotate the random point
                    p = np.array([xr,yr])
                    p = np.dot(R,p)
                    # Reposition the point
                    p = self.position + p
                    minDis = min(np.sqrt((np.append(fix, fx)-p[0])**2+(np.append(fiy, fy)-p[1])**2))

                    if minDis > h:
                        fx = np.append(fx, p[0])
                        fy = np.append(fy, p[1])

                    count = count +1

        else:
            maxit = 1000*mnf/hbyw

            while nop < mnf and count < maxit:
                # Get random points
                xr = random.uniform(0,1) - 0.5
                yr = random.uniform(0,1) - 0.5

                if yr > -1/(2/hbyw) and yr < 1/(2/hbyw):
                    # Scale the random point
                    xr = xr*self.width
                    yr = yr*self.width
                    # Rotate the random point
                    p = np.array([xr,yr])
                    p = np.dot(R,p)
                    # Reposition the point
                    p = self.position + p
                    minDis = min(np.sqrt((np.append(fix, fx)-p[0])**2+(np.append(fiy, fy)-p[1])**2))

                    if minDis > h:
                        fx = np.append(fx, p[0])
                        fy = np.append(fy, p[1])

                    count = count +1

        self.freeNodesX = fx[1:None]
        self.freeNodesY = fy[1:None]
        self.allNodesX = np.concatenate((self.boundaryNodesX, self.freeNodesX))
        self.allNodesY = np.concatenate((self.boundaryNodesY, self.freeNodesY))

        for k in self.childDomains:
            if initialDomains[k].shape == "point":
                self.allNodesX = np.concatenate((self.allNodesX, [initialDomains[k].allNodesX]))
                self.allNodesY = np.concatenate((self.allNodesY, [initialDomains[k].allNodesY]))
            else:
                self.allNodesX = np.concatenate((self.allNodesX, initialDomains[k].allNodesX))
                self.allNodesY = np.concatenate((self.allNodesY, initialDomains[k].allNodesY))

        plt.scatter(self.freeNodesX, self.freeNodesY, c=None, marker=markerShape, s = lW*1000)

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

        hbyw = self.height/self.width
        thr = self.angle
        R = np.array([[np.cos(np.radians(thr)), -np.sin(np.radians(thr))],
        [np.sin(np.radians(thr)), np.cos(np.radians(thr))]])

        fix = self.boundaryNodesX
        fiy = self.boundaryNodesY
        h = self.width/self.widthNodes/1.25*self.relativeMeshDensity[0]/self.relativeMeshDensity[1]

        if self.isParent: # Add child nodes as fixed nodes
            for cd in self.childDomains:
                fix = np.append(fix, initialDomains[cd].allNodesX)
                fiy = np.append(fiy, initialDomains[cd].allNodesY)

        fx = np.array([]) # Initializing free node positions
        fy = np.array([])

        if hbyw >= 1:

            wn = round(self.relativeWidth*self.relativeMeshDensity[1]*initialDomains[0].non)
            wn = round(hbyw*wn*0.8)
            xh = np.linspace(-0.788, 0.788, 2*wn-1)
            yh = 0*xh

            for i in range(wn-1):
                xd = np.linspace(-0.788/(2*wn-2)*(2*wn-i-3), 0.788/(2*wn-2)*(2*wn-i-3), 2*wn-i-2)
                xh = np.concatenate([xh, xd, xd])
                yp = 0*xd + 1.576/(2*wn-2)*np.sqrt(3)/2*(i+1)
                yn = -yp
                yh = np.concatenate([yh, yp, yn])

            for i in range(len(xh)):

                if xh[i] > -1/(2*hbyw) and xh[i] < 1/(2*hbyw) and yh[i] > -0.5 and yh[i] < 0.5:

                    xr = xh[i]
                    yr = yh[i]
                    # Scale the random point
                    xr = xr*self.height
                    yr = yr*self.height
                    # Rotate the random point
                    p = np.array([xr,yr])
                    p = np.dot(R,p)
                    # Reposition the point
                    p = self.position + p
                    minDis = min(np.sqrt((np.append(fix, fx)-p[0])**2+(np.append(fiy, fy)-p[1])**2))

                    if minDis > h:
                        fx = np.append(fx, p[0])
                        fy = np.append(fy, p[1])

        else:
            wn = round(self.relativeWidth*self.relativeMeshDensity[1]*initialDomains[0].non*0.8)
            xh = np.linspace(-0.788, 0.788, 2*wn-1)
            yh = 0*xh

            for i in range(wn-1):
                xd = np.linspace(-0.788/(2*wn-2)*(2*wn-i-3), 0.788/(2*wn-2)*(2*wn-i-3), 2*wn-i-2)
                xh = np.concatenate([xh, xd, xd])
                yp = 0*xd + 1.576/(2*wn-2)*np.sqrt(3)/2*(i+1)
                yn = -yp
                yh = np.concatenate([yh, yp, yn])

            for i in range(len(xh)):

                if yh[i] > -1/(2/hbyw) and yh[i] < 1/(2/hbyw) and xh[i] > -0.5 and xh[i] < 0.5:

                    xr = xh[i]
                    yr = yh[i]
                    # Scale the random point
                    xr = xr*self.width
                    yr = yr*self.width
                    # Rotate the random point
                    p = np.array([xr,yr])
                    p = np.dot(R,p)
                    # Reposition the point
                    p = self.position + p
                    minDis = min(np.sqrt((np.append(fix, fx)-p[0])**2+(np.append(fiy, fy)-p[1])**2))

                    if minDis > h:
                        fx = np.append(fx, p[0])
                        fy = np.append(fy, p[1])

        self.freeNodesX = fx[1:None]
        self.freeNodesY = fy[1:None]
        self.allNodesX = np.concatenate((self.boundaryNodesX, self.freeNodesX))
        self.allNodesY = np.concatenate((self.boundaryNodesY, self.freeNodesY))

        for k in self.childDomains:
            if initialDomains[k].shape == "point":
                self.allNodesX = np.concatenate((self.allNodesX, [initialDomains[k].allNodesX]))
                self.allNodesY = np.concatenate((self.allNodesY, [initialDomains[k].allNodesY]))
            else:
                self.allNodesX = np.concatenate((self.allNodesX, initialDomains[k].allNodesX))
                self.allNodesY = np.concatenate((self.allNodesY, initialDomains[k].allNodesY))

        plt.scatter(self.freeNodesX, self.freeNodesY, c=None, marker=markerShape, s = lW*1000)

