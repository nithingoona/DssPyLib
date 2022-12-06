import matplotlib.pyplot as plt
import numpy as np

from scipy.spatial import Delaunay as delaunay

class Mesh:
    """
    A class to represent mesh.
    
    ...

    Attributes
    ----------
    ix : list
        x positions of nodes
    iy : list
        y positions of nodes
    x : list
        x positions of nodes
    y : list
        y positions of nodes
    tri : list
        connectivity matrix
    

    Methods
    -------
    getNeiElements():
        generates list of neighbouring elements for each node
    getNeiNodes():
        generates list of neighbouring nodes for each node
    """

    def __init__(self, ix, iy, tri):
        """
        Constructs all the necessary attributes for the geometry object.

        Parameters
        ----------
        ix : list
            x positions of nodes
        iy : list
            y positions of nodes
        x : list
            x positions of nodes
        y : list
            y positions of nodes
        tri : list
            connectivity matrix
        """

        self.ix = ix
        self.iy = iy
        self.x = ix
        self.y = iy
        self.tri = tri

    def getNeiElements(self):
        """
        generates list of neighbouring elements for each node.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        nn = len(self.ix)
        ne = len(self.tri[:,0])
        neiElements = []

        for i in range(nn):
            ele = []
    
            for j in range(ne):

                if self.tri[j,0] == i or self.tri[j,1] == i or self.tri[j,2] == i:
                    ele.append(j)

            neiElements.append(ele)

        self.neiElements = neiElements

    def getNeiNodes(self):
        """
        generates list of neighbouring nodes for each node.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        nn = len(self.ix)
        neiNodes = []

        for i in range(nn):

            nod = []
            for j in self.neiElements[i]:

                for k in self.tri[j]:

                    if k != i: 

                        nod.append(k)

            nod = np.unique(np.array(nod))
            nod = nod.tolist()

            neiNodes.append(nod)

        self.neiNodes = neiNodes

def generateNodes(initialDomains, meshingMethod, figSize, showPlot):
    """
    gets nodes.

    Parameters
    ----------
    initialDomains : list
        list of geometries
    meshingMethod : str
        method of meshing
    figSize : list
        size of figure
    showPlot : boolean
        show plot if True

    Returns
    -------
    None
    """

    markers = ["+", "_", "*", "^", "x", "o"]
    plt.close("all")
    # Plotting the geometry again
    if figSize == "Default":
        figSize = initialDomains[0].figSize

    fig, ax = plt.subplots(figsize= figSize)
    w = initialDomains[0].width
    h = initialDomains[0].height
    [xp,yp] = initialDomains[0].position
    dx = w/100
    sf = 1.5
    plt.xlim(xp-w/2-dx, xp+w/2+dx)
    plt.ylim(yp-h/2-dx, yp+h/2+dx)

    lW = (figSize[0]/initialDomains[0].defaultFigSize[0])*0.5
    for domain in initialDomains:
        domain.show(plt, lW)

    # Generating additional data required for meshing

    # Relative width
    for domain in initialDomains:
        domain.calculateRelativeWidth(initialDomains[0].width)

    # Width Nodes
    for domain in initialDomains:
        domain.calculateWidthNodes(initialDomains[0].non)

    # Area
    for domain in initialDomains:
        domain.calculateArea()

    # Maximum number of free nodes
    for domain in initialDomains:
        domain.calculateMaxFreeNodes(initialDomains)

    # Actual meshing starts from here

    # Generate vertices
    count = 0
    for domain in initialDomains:
        domain.genVertices(plt, lW, markers[count % len(markers)])
        count = count + 1

    # Generate boundary nodes
    count = 0
    for domain in initialDomains:
        domain.genBoundaryNodes(plt, lW, markers[count % len(markers)])
        count = count + 1

    # Generate internal free nodes

    for domain in initialDomains:
        domain.isMeshed = False


    for j in range(len(initialDomains)): # Atleast one domain will be meshed in the next loop
                                        # This loop is to take care of the worst case scenario
        count = 0
        for domain in initialDomains: 
            
            if domain.isParent: # Dealing with parents with childs

                if (not domain.isMeshed): # Skip if already meshed

                    # Assume the childs are not meshed
                    domain.childMeshed = False
                    mc = 0 # Initialise no. of meshed childs
                    cd = domain.childDomains

                    # Check if all child domains are meshed
                    for k in cd:
                        if initialDomains[k].isMeshed:
                            mc = mc+1

                    if mc == len(cd):
                        domain.childMeshed = True

                    # Start meshing if all childs are meshed
                    if domain.childMeshed:
                        print("Meshing domain number", count)
                        if meshingMethod == "tri":
                            domain.genTriFreeNodes(initialDomains, plt, lW, markers[count % len(markers)])
                        elif meshingMethod == "rand":
                            domain.genRandFreeNodes(initialDomains, plt, lW, markers[count % len(markers)])
                        domain.isMeshed = True
                        
            else: # Dealing with pure child domains (Pure childs will be meshed first)

                if (not domain.isMeshed): # Skip if already meshed

                    print("Meshing domain number", count)
                    if meshingMethod == "tri":
                        domain.genTriFreeNodes(initialDomains, plt, lW, markers[count % len(markers)])
                    elif meshingMethod == "rand":
                        domain.genRandFreeNodes(initialDomains, plt, lW, markers[count % len(markers)])
                    domain.isMeshed = True

            count = count +1

    plt.gca().set_aspect("equal")
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(figSize[0]/initialDomains[0].defaultFigSize[0])
    
    plt.xlabel("x - axis [meter]", fontsize=figSize[0]/initialDomains[0].defaultFigSize[0]*1.5*20*sf)
    plt.ylabel("y - axis [meter]", fontsize=figSize[0]/initialDomains[0].defaultFigSize[0]*1.5*20*sf)
    plt.xticks(fontsize= figSize[0]/initialDomains[0].defaultFigSize[0]*20*sf)
    plt.yticks(fontsize= figSize[0]/initialDomains[0].defaultFigSize[0]*20*sf)
    plt.grid()
    plt.tick_params(axis='both', width=figSize[0]/initialDomains[0].defaultFigSize[0])
    plt.tick_params(axis='both', length=figSize[0]/initialDomains[0].defaultFigSize[0]*5)
    plt.savefig("Figures/1_Nodes.png")
    if showPlot == True:
        plt.show()

def createMesh(initialDomains):
    """
    gets mesh.

    Parameters
    ----------
    initialDomains : list
        list of geometries
    
    Returns
    -------
    None
    """

    # initialDomains = joblib.load("Data/initialDomains.sav")
    print("Creating mesh")
    # Generate boundary, free node list, x, y positions of all nodes

    x = np.array([])
    y = np.array([])
    nn = -1
    for domain in initialDomains:

        if domain.shape == "point":
            domain.boundaryList = np.array([nn+1])
            nn = nn+1
            domain.freeList = np.array([])
            x =  np.concatenate((x, [domain.boundaryNodesX]))
            y =  np.concatenate((y, [domain.boundaryNodesY]))

        elif domain.shape == "line":
            domain.boundaryList = np.arange(nn+1, nn+len(domain.boundaryNodesX)+1, 1)
            nn = domain.boundaryList[-1]
            domain.freeList = np.array([])
            x =  np.concatenate((x, domain.boundaryNodesX))
            y =  np.concatenate((y, domain.boundaryNodesY))

        else:
            domain.boundaryList = np.arange(nn+1, nn+len(domain.boundaryNodesX)+1, 1)
            nn = domain.boundaryList[-1]
            if len(domain.freeNodesX) == 0:
                domain.freeList = np.array([])
            else:
                domain.freeList = np.arange(nn+1, nn+len(domain.freeNodesX)+1, 1)
                nn = domain.freeList[-1]

            x =  np.concatenate((x, domain.boundaryNodesX, domain.freeNodesX))
            y =  np.concatenate((y, domain.boundaryNodesY, domain.freeNodesY))

    xy = np.transpose([x, y])
    tri = delaunay(xy)
    mesh = Mesh(x, y, tri.simplices)

    # Generate element list in each domain

    ne = len(tri.simplices[:,0])
    nn = len(mesh.ix)

    for domain in initialDomains:
        elementList = [] # Initializing element list

        if domain.shape != "point" and domain.shape != "line":

            if domain.parentDomains == []: # Checking for global domain to also include boundary nodes for element extraction
                if len(domain.freeNodesX) == 0:
                    domainNodes = domain.boundaryList
                else:
                    domainNodes = np.concatenate(domain.boundaryList, domain.freeList)
            else:
                domainNodes = domain.freeList
                boundaryNodes = domain.boundaryList
                
            domainNodes = domainNodes.astype(int)

            for cd in domain.childDomains:
                if initialDomains[cd].shape == "line":
                    domainNodes = np.concatenate((domainNodes, initialDomains[cd].boundaryList))

            for j in range(ne):

                for k in domainNodes:

                    if tri.simplices[j,0] == k or tri.simplices[j,1] == k or tri.simplices[j,2] == k:
                        elementList.append(j)
                        break

        domain.elementList = elementList

    # To assign interior elements formed by only the boundary nodes

    actualElementList = np.arange(0,ne,1).tolist()
    currentElementList = np.array([]).tolist()
    for domain in initialDomains:
        currentElementList = currentElementList + domain.elementList

    remainingElements = actualElementList.copy()

    for element in actualElementList:
        if element in currentElementList:
            remainingElements.remove(element)

    for domain in initialDomains:
        boundaryNodes = domain.boundaryList
        for element in remainingElements:
            for node in boundaryNodes:
                if node in mesh.tri[element]:
                    domain.elementList = domain.elementList + [element]
                    break

    # Extracting neighbour information
    print("Extracting neighbour information")
    # Elements
    mesh.getNeiElements()
    # Nodes
    mesh.getNeiNodes()
    return mesh

def relaxMesh(initialDomains, mesh):
    """
    relaxes the mesh.

    Parameters
    ----------
    initialDomains : list
        list of geometries
    mesh : object
        the mesh object

    Returns
    -------
    None
    """

    print("Relaxing mesh")
    
    rx = mesh.ix
    ry = mesh.iy
    nn = len(rx)

    Fx = np.linspace(0, 0, len(rx))
    Fy = Fx

    rMin = initialDomains[0].width/initialDomains[0].non

    for rm in range(1000):

        print("Percentage of relaxation", rm/1000*100, end='\r')

        for domain in initialDomains:

            if domain.shape != "point" and domain.shape != "line":

                frn = domain.freeList

                for j in frn:
                    nein = mesh.neiNodes[j]
                    Fx[j] = 0
                    Fy[j] = 0

                    for k in nein:
                        stiff = 1
                        rmag = np.sqrt((rx[k]-rx[j])**2+(ry[k]-ry[j])**2)
                        Fx[j] = Fx[j] + stiff*(rx[k]-rx[j])/rmag**0.1
                        Fy[j] = Fy[j] + stiff*(ry[k]-ry[j])/rmag**0.1

                for j in frn:
                    if Fx[j] != 0 and Fy[j] != 0:
                        rx[j] = rx[j] + 0.002*rMin*Fx[j]/np.sqrt(Fx[j]**2+Fy[j]**2)
                        ry[j] = ry[j] + 0.002*rMin*Fy[j]/np.sqrt(Fx[j]**2+Fy[j]**2)

    print()

    mesh.rx = rx
    mesh.ry = ry
    mesh.x = rx
    mesh.y = ry

def assignSources(initialDomains, mesh, solver):
    """
    assigns sources.

    Parameters
    ----------
    initialDomains : list
        list of geometries
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

    nn = len(x)
    # nodeSource is for both plotting and matrix assignment
    # elementSource may serve some purpose in different physics
    if solver.name == 'structural':

        nodeSource = np.zeros((nn,2))
        ne = len(mesh.tri[:,0])
        elementSource = np.zeros((ne,2))
        print("Assigning sources")

        for domain in initialDomains:

            node = domain.sourcePosition
            nodeSource[node] = domain.source

    else:
        nodeSource = np.zeros((nn,1))
        ne = len(mesh.tri[:,0])
        elementSource = np.zeros((ne,1))
        print("Assigning sources")

        for domain in initialDomains:

            if domain.shape == "point" or domain.shape == "line":
                nodes = domain.boundaryList
                for node in nodes:
                    nodeSource[node] = domain.source/len(nodes)
                    elements = mesh.neiElements[node]
                    for element in elements:
                        elementSource[element] = elementSource[element] + domain.source/len(elements)

            else:
                elements = domain.elementList
                for element in elements:
                    ex = []
                    ey = []
                    for node in mesh.tri[element]:
                        ex.append(x[node])
                        ey.append(y[node])

                    p = [ey[1]-ey[2], ey[2]-ey[0], ey[0]-ey[1]]
                    q = [ex[2]-ex[1], ex[0]-ex[2], ex[1]-ex[0]]
                    ea = 0.5*abs(p[1]*q[2]-q[1]*p[2])

                    if ea == 0:
                        print("element area is zero at: ", element)

                    elementSource[element] = domain.source*ea

                for element in elements:

                    for node in mesh.tri[element]:
                        nodeSource[node] = nodeSource[node] + elementSource[element]/3

        mesh.elementSource = elementSource
                    
    mesh.nodeSource = nodeSource
    
def assignMaterials(initialDomains, mesh):
    """
    assigns materials.

    Parameters
    ----------
    initialDomains : list
        list of geometries
    mesh : object
        the mesh object

    Returns
    -------
    None
    """

    # nodeMaterial is for plotting
    # elementMaterial is for matrix computation

    x = mesh.rx
    nn = len(x)
    nodeMaterial = np.ones((nn,1))
    ne = len(mesh.tri[:,0])
    elementMaterial = np.ones((ne,1))
    print("Assigning material properties")

    for domain in initialDomains:

        if domain.shape == "point" or domain.shape == "line":
            nodes = domain.boundaryList
            for node in nodes:
                nodeMaterial[node] = domain.material

        else:
            nodes = np.concatenate((domain.boundaryList, domain.freeList)).astype(int)
            for node in nodes:
                nodeMaterial[node] = domain.material

        elements = domain.elementList
        for element in elements:
            elementMaterial[element] = domain.material

    mesh.nodeMaterial = nodeMaterial
    mesh.elementMaterial = elementMaterial

def getNeighbourInformation(mesh, maxTn):
    """
    gets neighbour information for DSS.

    Parameters
    ----------
    mesh : object
        the mesh object
    maxTn : int
        truncation number, all the lower numbered information will be automatically generated

    Returns
    -------
    None
    """
        
    tri = mesh.tri
    ne = len(tri[:,0])
    nn = len(mesh.rx)

    print("Generating neighbour element and node list upto Truncation Number", maxTn)
    eleAttNames = []
    nodeAttNames = []
    farnodeAttNames = []

    for i in range(maxTn+1): # Initializing all truncated element list
        globals()['neiElementsTn%s' % i] = np.empty([nn,0]).tolist()
        setattr(mesh, 'neiElementsTn%s' % i, globals()['neiElementsTn%s' % i])
        eleAttNames.append('neiElementsTn%s' % i)

    for i in range(nn):
        for j in range(ne):
            if tri[j,0] == i or tri[j,1] == i or tri[j,2] == i:
                mesh.neiElementsTn0[i].append(j)

    for i in range(maxTn+1): # Initializing all truncated node list
        globals()['neiNodesTn%s' % i] = np.empty([nn,0]).tolist()
        setattr(mesh, 'neiNodesTn%s' % i, globals()['neiNodesTn%s' % i])
        nodeAttNames.append('neiNodesTn%s' % i)

    for i in range(nn):
        mesh.neiNodesTn0[i].append(i) # Redundant. However, completes the set.

    for i in range(1, maxTn+1): # Initializing all far neighbour list
        globals()['farneiNodesTn%s' % i] = np.empty([nn,0]).tolist()
        setattr(mesh, 'farneiNodesTn%s' % i, globals()['farneiNodesTn%s' % i])
        farnodeAttNames.append('farneiNodesTn%s' % i)

    for tn in range(1,maxTn+1):
        neiNodesTn = np.empty([nn,0]).tolist()
        neiElementsTn = np.empty([nn,0]).tolist()
        farneiNodesTn = np.empty([nn,0]).tolist()

        for i in range(nn):
            existingNodes = getattr(mesh, nodeAttNames[tn-1])[i]
            newNodes = []

            for node in existingNodes:
                potElements = getattr(mesh, eleAttNames[0])[node]
                for k in potElements:
                    for l in range(3):
                        newNodes.append(tri[k,l])

            newNodes = np.unique(newNodes)
            newNodesList = newNodes.tolist()

            for enode in existingNodes:
                if enode in newNodes:
                    newNodesList.remove(enode)

            newNodes = newNodesList

            totalNodes = existingNodes.copy()
            for node in newNodes:
                totalNodes.append(node)

            neiNodesTn[i] = totalNodes.copy()
            farneiNodesTn[i] = newNodes.copy()

            truncatedElements = []

            for node in existingNodes:
                potElements = getattr(mesh, eleAttNames[0])[node]
                for k in potElements:
                    truncatedElements.append(k)

            truncatedElements = np.unique(truncatedElements)

            neiElementsTn[i] = truncatedElements.tolist()

        print("Current truncation number", tn, end='\r')
        setattr(mesh, nodeAttNames[tn], neiNodesTn)
        setattr(mesh, eleAttNames[tn], neiElementsTn)
        setattr(mesh, farnodeAttNames[tn-1], farneiNodesTn)

    print(end= '\r')
    mesh.eleAttNames = eleAttNames
    mesh.nodeAttNames = nodeAttNames
    mesh.farnodeAttNames = farnodeAttNames

            
            




