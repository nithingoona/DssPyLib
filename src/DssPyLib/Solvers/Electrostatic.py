from . import Solver
import numpy as np

class Electrostatic(Solver.Solver):
    """
    A subclass to represent a specific Solver.
    
    ...

    Attributes
    ----------
    name : str
        name of the solver
    physicalConstant : float
        physical constant associated with the solver
    proportionalityConstant : float
        proportionality constant associated with the solver
    
    Methods
    -------
    getApproximateMatrices(mesh, initialDomains):
        gets approximate matrices
    getLinNumericalMatrix(mesh):
        gets linear numerical matrices
    getNumericalMatrix(mesh):
        gets non linear numerical matrices
    getApproximateBoundaryMatrices(mesh, initialDomains):
        gets approximate matrices on boundary
    getIterBoundaryMatrices(mesh, initialDomains, bt):
        gets approximate matrices on boundary for iterative update
    getMiniLinNumericalMatrix(minix, miniy, miniTri):
        gets linear numerical matrices
    findVectorFields(mesh, solver):
        finds gradient vector fields
    """

    def __init__(self):
        """
        A subclass to represent a specific Solver.
    
        Parameters
        ----------
        name : str
            name of the solver
        physicalConstant : float
            physical constant associated with the solver
        proportionalityConstant : float
            proportionality constant associated with the solver
        """

        Solver.Solver.__init__(self)
        self.name = "electrostatic"
        self.physicalConstant = 8.8541878128*10**-12
        self.proportionalityConstant = 1/(2*np.pi*self.physicalConstant)

    def getApproximateMatrices(self, mesh, initialDomains):
        """
        gets approximate matrices.

        Parameters
        ----------
        mesh : object
            the mesh object
        initialDomains : list
            list of geometries

        Returns
        -------
        None
        """

        x = mesh.rx
        y = mesh.ry
        tri = mesh.tri
        nn = len(x)
        ne = len(tri[:,0])
        pc = self.proportionalityConstant
        A = np.zeros((nn,nn))
        b = np.zeros((nn,1))

        for i in range(ne):

            print("Percentage of elements covered", i/ne*100, end='\r')

            # 0-1 edge
            k = 0
            l = 1
            Ex = 0
            Ey = 0

            for domain in initialDomains:
                if domain.source != 0:
                    if domain.shape == "point" or domain.shape == "line":
                        nodes = domain.boundaryList
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (x[tri[i,k]] + x[tri[i,l]])/2 - px
                            Ry = (y[tri[i,k]] + y[tri[i,l]])/2 - py
                            Ex = Ex + mesh.nodeSource[node]*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + mesh.nodeSource[node]*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

                    else:
                        nodes = np.concatenate((domain.boundaryList, domain.freeList)).astype(int)
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (x[tri[i,k]] + x[tri[i,l]])/2 - px
                            Ry = (y[tri[i,k]] + y[tri[i,l]])/2 - py
                            Ex = Ex + mesh.nodeSource[node]*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + mesh.nodeSource[node]*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

            rx = x[tri[i,l]] - x[tri[i,k]]
            ry = y[tri[i,l]] - y[tri[i,k]]
            Er = np.dot([Ex[0], Ey[0]], [rx, ry]/np.sqrt(rx**2+ry**2))
            A[tri[i,k], tri[i,k]] = A[tri[i,k], tri[i,k]] +1
            A[tri[i,k], tri[i,l]] = -1
            b[tri[i,k]]           = b[tri[i,k]] + Er*np.sqrt(rx**2+ry**2)

            # 1-2 edge
            k = 1
            l = 2
            Ex = 0
            Ey = 0

            for domain in initialDomains:
                if domain.source != 0:
                    if domain.shape == "point" or domain.shape == "line":
                        nodes = domain.boundaryList
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (x[tri[i,k]] + x[tri[i,l]])/2 - px
                            Ry = (y[tri[i,k]] + y[tri[i,l]])/2 - py
                            Ex = Ex + mesh.nodeSource[node]*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + mesh.nodeSource[node]*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

                    else:
                        nodes = np.concatenate((domain.boundaryList, domain.freeList)).astype(int)
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (x[tri[i,k]] + x[tri[i,l]])/2 - px
                            Ry = (y[tri[i,k]] + y[tri[i,l]])/2 - py
                            Ex = Ex + mesh.nodeSource[node]*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + mesh.nodeSource[node]*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

            rx = x[tri[i,l]] - x[tri[i,k]]
            ry = y[tri[i,l]] - y[tri[i,k]]
            Er = np.dot([Ex[0], Ey[0]], [rx, ry]/np.sqrt(rx**2+ry**2))
            A[tri[i,k], tri[i,k]] = A[tri[i,k], tri[i,k]] +1
            A[tri[i,k], tri[i,l]] = -1
            b[tri[i,k]]           = b[tri[i,k]] + Er*np.sqrt(rx**2+ry**2)

            # 2-0 edge
            k = 2
            l = 0
            Ex = 0
            Ey = 0

            for domain in initialDomains:
                if domain.source != 0:
                    if domain.shape == "point" or domain.shape == "line":
                        nodes = domain.boundaryList
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (x[tri[i,k]] + x[tri[i,l]])/2 - px
                            Ry = (y[tri[i,k]] + y[tri[i,l]])/2 - py
                            Ex = Ex + mesh.nodeSource[node]*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + mesh.nodeSource[node]*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

                    else:
                        nodes = np.concatenate((domain.boundaryList, domain.freeList)).astype(int)
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (x[tri[i,k]] + x[tri[i,l]])/2 - px
                            Ry = (y[tri[i,k]] + y[tri[i,l]])/2 - py
                            Ex = Ex + mesh.nodeSource[node]*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + mesh.nodeSource[node]*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

            rx = x[tri[i,l]] - x[tri[i,k]]
            ry = y[tri[i,l]] - y[tri[i,k]]
            Er = np.dot([Ex[0], Ey[0]], [rx, ry]/np.sqrt(rx**2+ry**2))
            A[tri[i,k], tri[i,k]] = A[tri[i,k], tri[i,k]] +1
            A[tri[i,k], tri[i,l]] = -1
            b[tri[i,k]]           = b[tri[i,k]] + Er*np.sqrt(rx**2+ry**2)

        print()
        self.Aa = A
        self.ba = b

    def getLinNumericalMatrix(self, mesh):
        """
        gets linear numerical matrices

        Parameters
        ----------
        mesh : object
            the mesh object

        Returns
        -------
        None
        """

        x = mesh.rx
        y = mesh.ry
        tri = mesh.tri
        nn = len(x)
        A = np.zeros((nn,nn))

        count = 0
        ne = len(tri[:,0])
        for i in range(ne):
            print("percentage of linear elements covered", count/ne*100, end='\r')
            count = count +1
            ex = np.array([x[tri[i,0]], x[tri[i,1]], x[tri[i,2]]])
            ey = np.array([y[tri[i,0]], y[tri[i,1]], y[tri[i,2]]])
            p = np.array([ey[1]-ey[2], ey[2]-ey[0], ey[0]-ey[1]])
            q = np.array([ex[2]-ex[1], ex[0]-ex[2], ex[1]-ex[0]])
            ea = 0.5*abs(p[1]*q[2]-q[1]*p[2])
            P = np.zeros([len(p),len(p)])
            Q = P.copy()
            for j in range(len(p)):
                P[j] = p*p[j]
                Q[j] = q*q[j]

            a = (P+Q)/(4*ea)
            for k in range(3):
                for j in range(3):
                    A[tri[i,k],tri[i,j]] = A[tri[i,k],tri[i,j]].copy() + a[k,j]

        print()
        self.AnLin = A.copy()

    def getNumericalMatrix(self, mesh):
        """
        gets non linear numerical matrices
        Parameters
        ----------
        mesh : object
            the mesh object

        Returns
        -------
        None
        """

        x = mesh.rx
        y = mesh.ry
        tri = mesh.tri
        nn = len(x)
        A = np.zeros((nn,nn))
        er = mesh.elementMaterial

        count = 0
        ne = len(tri[:,0])
        for i in range(ne):
            print("percentage of non-linear elements covered", count/ne*100, end='\r')
            count = count +1
            ex = np.array([x[tri[i,0]], x[tri[i,1]], x[tri[i,2]]])
            ey = np.array([y[tri[i,0]], y[tri[i,1]], y[tri[i,2]]])
            p = np.array([ey[1]-ey[2], ey[2]-ey[0], ey[0]-ey[1]])
            q = np.array([ex[2]-ex[1], ex[0]-ex[2], ex[1]-ex[0]])
            ea = 0.5*abs(p[1]*q[2]-q[1]*p[2])
            P = np.zeros([len(p),len(p)])
            Q = P.copy()
            for j in range(len(p)):
                P[j] = p*p[j]
                Q[j] = q*q[j]

            a = er[i]*(P+Q)/(4*ea)
            for k in range(3):
                for j in range(3):
                    A[tri[i,k],tri[i,j]] = A[tri[i,k],tri[i,j]].copy() + a[k,j]

        print()
        self.An = A.copy()

    def getApproximateBoundaryMatrices(self, mesh, initialDomains):
        """
        gets approximate matrices on boundary

        Parameters
        ----------
        mesh : object
            the mesh object
        initialDomains : list
            list of geometries

        Returns
        -------
        None
        """

        x = mesh.rx
        y = mesh.ry
        tri = mesh.tri
        pc = self.proportionalityConstant

        # Boundary nodes
        gc = initialDomains[0].childDomains
        lbn = initialDomains[0].boundaryList
        lbn = np.append([lbn], [initialDomains[0].freeList])
        lbn = np.append([lbn], [initialDomains[gc[0]].boundaryList])
        lbn = lbn.astype(int)
        bx = x[lbn].copy()
        by = y[lbn].copy()

        # Boundary elements
        belements = initialDomains[0].elementList
        btri = tri[belements].copy()
        bne = len(btri[:,0])
        bnn = len(bx)

        # Mapping btri and tri

        for j in range(bne):
            for k in range(3):
                gn = btri[j,k]

                for l in range(bnn):
                    if gn == lbn[l]:
                        btri[j,k] = l
                        break
        
        # Initializing matrices
        A = np.zeros([bnn,bnn])
        b = np.zeros([bnn,1])

        # For every boundary element
        for i in range(bne):

            print("Percentage of boundary elements covered", i/bne*100, end='\r')

            # 0-1 edge
            k = 0
            l = 1
            Ex = 0
            Ey = 0

            for domain in initialDomains:
                if domain.source != 0:
                    if domain.shape == "point" or domain.shape == "line":
                        nodes = domain.boundaryList
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (bx[btri[i,k]] + bx[btri[i,l]])/2 - px
                            Ry = (by[btri[i,k]] + by[btri[i,l]])/2 - py
                            Ex = Ex + mesh.nodeSource[node]*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + mesh.nodeSource[node]*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

                    else:
                        nodes = np.concatenate((domain.boundaryList, domain.freeList)).astype(int)
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (bx[btri[i,k]] + bx[btri[i,l]])/2 - px
                            Ry = (by[btri[i,k]] + by[btri[i,l]])/2 - py
                            Ex = Ex + mesh.nodeSource[node]*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + mesh.nodeSource[node]*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

            rx = bx[btri[i,l]] - bx[btri[i,k]]
            ry = by[btri[i,l]] - by[btri[i,k]]
            Er = np.dot([Ex[0], Ey[0]], [rx, ry]/np.sqrt(rx**2+ry**2))
            A[btri[i,k], btri[i,k]] = A[btri[i,k], btri[i,k]] +1
            A[btri[i,k], btri[i,l]] = -1
            b[btri[i,k]]           = b[btri[i,k]] + Er*np.sqrt(rx**2+ry**2)

            # 1-2 edge
            k = 1
            l = 2
            Ex = 0
            Ey = 0

            for domain in initialDomains:
                if domain.source != 0:
                    if domain.shape == "point" or domain.shape == "line":
                        nodes = domain.boundaryList
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (bx[btri[i,k]] + bx[btri[i,l]])/2 - px
                            Ry = (by[btri[i,k]] + by[btri[i,l]])/2 - py
                            Ex = Ex + mesh.nodeSource[node]*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + mesh.nodeSource[node]*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

                    else:
                        nodes = np.concatenate((domain.boundaryList, domain.freeList)).astype(int)
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (bx[btri[i,k]] + bx[btri[i,l]])/2 - px
                            Ry = (by[btri[i,k]] + by[btri[i,l]])/2 - py
                            Ex = Ex + mesh.nodeSource[node]*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + mesh.nodeSource[node]*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

            rx = bx[btri[i,l]] - bx[btri[i,k]]
            ry = by[btri[i,l]] - by[btri[i,k]]
            Er = np.dot([Ex[0], Ey[0]], [rx, ry]/np.sqrt(rx**2+ry**2))
            A[btri[i,k], btri[i,k]] = A[btri[i,k], btri[i,k]] +1
            A[btri[i,k], btri[i,l]] = -1
            b[btri[i,k]]           = b[btri[i,k]] + Er*np.sqrt(rx**2+ry**2)

            # 2-0 edge
            k = 2
            l = 0
            Ex = 0
            Ey = 0

            for domain in initialDomains:
                if domain.source != 0:
                    if domain.shape == "point" or domain.shape == "line":
                        nodes = domain.boundaryList
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (bx[btri[i,k]] + bx[btri[i,l]])/2 - px
                            Ry = (by[btri[i,k]] + by[btri[i,l]])/2 - py
                            Ex = Ex + mesh.nodeSource[node]*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + mesh.nodeSource[node]*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

                    else:
                        nodes = np.concatenate((domain.boundaryList, domain.freeList)).astype(int)
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (bx[btri[i,k]] + bx[btri[i,l]])/2 - px
                            Ry = (by[btri[i,k]] + by[btri[i,l]])/2 - py
                            Ex = Ex + mesh.nodeSource[node]*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + mesh.nodeSource[node]*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

            rx = bx[btri[i,l]] - bx[btri[i,k]]
            ry = by[btri[i,l]] - by[btri[i,k]]
            Er = np.dot([Ex[0], Ey[0]], [rx, ry]/np.sqrt(rx**2+ry**2))
            A[btri[i,k], btri[i,k]] = A[btri[i,k], btri[i,k]] +1
            A[btri[i,k], btri[i,l]] = -1
            b[btri[i,k]]           = b[btri[i,k]] + Er*np.sqrt(rx**2+ry**2)

        print()
        self.Aab = A
        self.bab = b
        mesh.bx  = bx
        mesh.by  = by
        mesh.btri= btri

    def getIterBoundaryMatrices(self, mesh, initialDomains, bt):
        """
        gets approximate matrices on boundary for iterative update

        Parameters
        ----------
        mesh : object
            the mesh object
        initialDomains : list
            list of geometries
        bt : list

        Returns
        -------
        None
        """

        x = mesh.rx
        y = mesh.ry
        tri = mesh.tri
        ep0 = self.physicalConstant
        pc = self.proportionalityConstant

        # Boundary nodes
        gc = initialDomains[0].childDomains
        lbn = initialDomains[0].boundaryList
        lbn = np.append([lbn], [initialDomains[0].freeList])
        lbn = np.append([lbn], [initialDomains[gc[0]].boundaryList])
        lbn = lbn.astype(int)
        bx = x[lbn].copy()
        by = y[lbn].copy()

        # Boundary elements
        belements = initialDomains[0].elementList
        btri = tri[belements].copy()
        bne = len(btri[:,0])
        bnn = len(bx)

        # Mapping btri and tri

        for j in range(bne):
            for k in range(3):
                gn = btri[j,k]

                for l in range(bnn):
                    if gn == lbn[l]:
                        btri[j,k] = l
                        break
        
        # Initializing matrices
        A = np.zeros([bnn,bnn])
        b = np.zeros([bnn,1])

        # For every boundary element
        for i in range(bne):

            print("Percentage of boundary elements covered", i/bne*100, end='\r')

            # 0-1 edge
            k = 0
            l = 1
            Ex = 0
            Ey = 0

            for domain in initialDomains:
                if domain.source != 0 or domain.material != 1:
                    if domain.shape == "point" or domain.shape == "line":
                        nodes = domain.boundaryList
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (bx[btri[i,k]] + bx[btri[i,l]])/2 - px
                            Ry = (by[btri[i,k]] + by[btri[i,l]])/2 - py
                            Ex = Ex + bt[node]*ep0*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + bt[node]*ep0*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

                    else:
                        nodes = np.concatenate((domain.boundaryList, domain.freeList)).astype(int)
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (bx[btri[i,k]] + bx[btri[i,l]])/2 - px
                            Ry = (by[btri[i,k]] + by[btri[i,l]])/2 - py
                            Ex = Ex + bt[node]*ep0*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + bt[node]*ep0*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

            rx = bx[btri[i,l]] - bx[btri[i,k]]
            ry = by[btri[i,l]] - by[btri[i,k]]
            Er = np.dot([Ex[0], Ey[0]], [rx, ry]/np.sqrt(rx**2+ry**2))
            A[btri[i,k], btri[i,k]] = A[btri[i,k], btri[i,k]] +1
            A[btri[i,k], btri[i,l]] = -1
            b[btri[i,k]]           = b[btri[i,k]] + Er*np.sqrt(rx**2+ry**2)

            # 1-2 edge
            k = 1
            l = 2
            Ex = 0
            Ey = 0

            for domain in initialDomains:
                if domain.source != 0 or domain.material != 1:
                    if domain.shape == "point" or domain.shape == "line":
                        nodes = domain.boundaryList
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (bx[btri[i,k]] + bx[btri[i,l]])/2 - px
                            Ry = (by[btri[i,k]] + by[btri[i,l]])/2 - py
                            Ex = Ex + bt[node]*ep0*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + bt[node]*ep0*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

                    else:
                        nodes = np.concatenate((domain.boundaryList, domain.freeList)).astype(int)
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (bx[btri[i,k]] + bx[btri[i,l]])/2 - px
                            Ry = (by[btri[i,k]] + by[btri[i,l]])/2 - py
                            Ex = Ex + bt[node]*ep0*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + bt[node]*ep0*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

            rx = bx[btri[i,l]] - bx[btri[i,k]]
            ry = by[btri[i,l]] - by[btri[i,k]]
            Er = np.dot([Ex[0], Ey[0]], [rx, ry]/np.sqrt(rx**2+ry**2))
            A[btri[i,k], btri[i,k]] = A[btri[i,k], btri[i,k]] +1
            A[btri[i,k], btri[i,l]] = -1
            b[btri[i,k]]           = b[btri[i,k]] + Er*np.sqrt(rx**2+ry**2)

            # 2-0 edge
            k = 2
            l = 0
            Ex = 0
            Ey = 0

            for domain in initialDomains:
                if domain.source != 0 or domain.material != 1:
                    if domain.shape == "point" or domain.shape == "line":
                        nodes = domain.boundaryList
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (bx[btri[i,k]] + bx[btri[i,l]])/2 - px
                            Ry = (by[btri[i,k]] + by[btri[i,l]])/2 - py
                            Ex = Ex + bt[node]*ep0*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + bt[node]*ep0*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

                    else:
                        nodes = np.concatenate((domain.boundaryList, domain.freeList)).astype(int)
                        for node in nodes:
                            px = x[node]
                            py = y[node]
                            Rx = (bx[btri[i,k]] + bx[btri[i,l]])/2 - px
                            Ry = (by[btri[i,k]] + by[btri[i,l]])/2 - py
                            Ex = Ex + bt[node]*ep0*pc*Rx/(np.sqrt(Rx**2 + Ry**2))**2
                            Ey = Ey + bt[node]*ep0*pc*Ry/(np.sqrt(Rx**2 + Ry**2))**2

            rx = bx[btri[i,l]] - bx[btri[i,k]]
            ry = by[btri[i,l]] - by[btri[i,k]]
            Er = np.dot([Ex[0], Ey[0]], [rx, ry]/np.sqrt(rx**2+ry**2))
            A[btri[i,k], btri[i,k]] = A[btri[i,k], btri[i,k]] +1
            A[btri[i,k], btri[i,l]] = -1
            b[btri[i,k]]           = b[btri[i,k]] + Er*np.sqrt(rx**2+ry**2)

        print()
        self.Aib = A
        self.bib = b
        mesh.bx  = bx
        mesh.by  = by
        mesh.btri= btri

    def getMiniLinNumericalMatrix(self, minix, miniy, miniTri):
        """
        gets linear numerical matrices

        Parameters
        ----------
        minix : list
            x positions of minimesh
        miniy : list
            y positions of minimesh
        miniTri : list
            connectivity matrix of minimesh
        Returns
        -------
        None
        """

        x = minix
        y = miniy
        tri = miniTri
        nn = len(x)
        A = np.zeros((nn,nn))

        count = 0
        ne = len(tri[:,0])
        for i in range(ne):
            # print("percentage of linear elements covered", count/ne*100, end='\r')
            count = count +1
            ex = np.array([x[tri[i,0]], x[tri[i,1]], x[tri[i,2]]])
            ey = np.array([y[tri[i,0]], y[tri[i,1]], y[tri[i,2]]])
            p = np.array([ey[1]-ey[2], ey[2]-ey[0], ey[0]-ey[1]])
            q = np.array([ex[2]-ex[1], ex[0]-ex[2], ex[1]-ex[0]])
            ea = 0.5*abs(p[1]*q[2]-q[1]*p[2])
            P = np.zeros([len(p),len(p)])
            Q = P.copy()
            for j in range(len(p)):
                P[j] = p*p[j]
                Q[j] = q*q[j]

            a = (P+Q)/(4*ea)
            for k in range(3):
                for j in range(3):
                    A[tri[i,k],tri[i,j]] = A[tri[i,k],tri[i,j]].copy() + a[k,j]

        # print()
        self.miniAn = A.copy()

    def findVectorFields(self, mesh, solver):
        """
        finds gradient vector fields

        Parameters
        ----------
        mesh : object
            the mesh object
        solver : object
            the solver object

        Returns
        -------
        None
        """

        x = mesh.rx
        y = mesh.ry
        u = solver.u.copy()

        Ex = np.zeros([len(x),1])
        Ey = np.zeros([len(x),1])
        farnodeAttNames = mesh.farnodeAttNames.copy()
        neiNodesTn = getattr(mesh, farnodeAttNames[0])

        for i in range(len(x)):
            tNodes = neiNodesTn[i]
            for node in tNodes:
                rx = x[node]-x[i]
                ry = y[node]-y[i]
                # print(rx, ry, (np.sqrt(rx**2 + ry**2))**2)
                Ex[i] = Ex[i] + (u[i] - u[node]) * rx / (np.sqrt(rx**2 + ry**2))**2
                Ey[i] = Ey[i] + (u[i] - u[node]) * ry / (np.sqrt(rx**2 + ry**2))**2

        solver.Ex = Ex
        solver.Ey = Ey