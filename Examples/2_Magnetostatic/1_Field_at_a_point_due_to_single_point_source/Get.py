# Import packages
import joblib
# import sys
import numpy as np

def vectorField(location, isNode = False, nodeNumber = -1):
    """
    calculates vector fields

    Parameters
    ----------
    location : int or list
        location of point where vector field is to be calculated
    isNode : boolean
        to mention if selected point is node or not
    nodeNumber : int
        number of the node if node is selected

    Returns
    -------
    field : list
        list of x and y components of vector fields
    """

    if isinstance(location, int):
        domains = joblib.load("Data/initialDomains.sav")
        if domains[location].shape == 'point':
            nodeNumber = domains[location].boundaryList[0]
            return(vectorField(np.array(domains[location].position), True, nodeNumber))
        elif domains[location-1].shape == 'line' and domains[location].shape == 'line' and domains[location+1].shape == 'line':
            print("Field along a line is in progress")
            mesh = joblib.load("Data/mesh.sav")
            solver = joblib.load("Data/solver.sav")
            x = mesh.x
            y = mesh.y
            tri = mesh.tri
            u = solver.u

            l1 = domains[location-1].boundaryList[0]
            l2 = domains[location].boundaryList[0]
            l3 = domains[location+1].boundaryList[0]

            field = []
            for i in range(1,len(l2)-1):
                Exp = (l2[i]-l2[i+1])/(np.sqrt((x[l2[i]] - x[l2[i+1]])**2+(y[l2[i]] - y[l2[i+1]])**2))
                Exn = (l2[i-1]-l2[i])/(np.sqrt((x[l2[i-1]] - x[l2[i]])**2+(y[l2[i-1]] - y[l2[i]])**2))
                Ex = (Exp+Exn)/2
                Eyp = (l2[i]-l1[i])/(np.sqrt((x[l2[i]] - x[l1[i]])**2+(y[l2[i]] - y[l1[i]])**2))
                Eyn = (l3[i]-l2[i])/(np.sqrt((x[l3[i]] - x[l2[i]])**2+(y[l3[i]] - y[l2[i]])**2))
                Ey = (Eyp+Eyn)/2
                field.append([Ex, Ey])

            return field

        else:
            print("Field along the line cannot be obtained for the selected domain", location)
    elif isinstance(location, np.ndarray):
        domains = joblib.load("Data/initialDomains.sav")
        mesh = joblib.load("Data/mesh.sav")
        solver = joblib.load("Data/solver.sav")
        x = mesh.x
        y = mesh.y
        tri = mesh.tri
        u = solver.u

        if isNode:
            field = [0, 0]
            # node = domains[domainNumber].boundaryList[0]
            eleAttNames = mesh.eleAttNames.copy()
            neiElements = getattr(mesh, eleAttNames[0])[nodeNumber]
            elementfields = []
            for element in neiElements:
                A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
                for k in range(3):
                    A[k] = [1, x[tri[element][k]], y[tri[element][k]]]
                b = [[u[tri[element][0]][0]], [u[tri[element][1]][0]], [u[tri[element][2]][0]]]
                X = np.linalg.lstsq(A, b)[0]
                elementfield = [-X[1][0], -X[2][0]]
                elementfields.append(elementfield)
            elementfields = np.array(elementfields)
            field = sum(elementfields)/len(elementfields)
        else:
            field = []
            noOfRows = location.shape[0]
            for i in range(noOfRows):
                field.append([0, 0])

            for j in range(noOfRows):
                isInsideGlobalBoundary = False
                for i in range(len(tri)):
                    if isInside(x[tri[i][0]], y[tri[i][0]], x[tri[i][1]], y[tri[i][1]], x[tri[i][2]], y[tri[i][2]], location[j][0], location[j][1]):
                        theElement = i
                        isInsideGlobalBoundary = True
                        break

                if isInsideGlobalBoundary:
                    A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
                    for k in range(3):
                        A[k] = [1, x[tri[theElement][k]], y[tri[theElement][k]]]
                    b = [[u[tri[theElement][0]][0]], [u[tri[theElement][1]][0]], [u[tri[theElement][2]][0]]]
                    X = np.linalg.lstsq(A, b)[0]
                    field[j] = [-X[1][0], -X[2][0]]
                else:
                    field[j] = [[0, 0]]
                    print("Position", location[j], "is out of Global Boundary")

        if solver.name == "magnetostatic":
            print("Turn field 90 degrees")
        return field

def area(x1, y1, x2, y2, x3, y3):
    """
    calculates area

    Parameters
    ----------
    x1 : float
        x position of vertix 1
    y1 : float
        y position of vertix 1
    x2 : float
        x position of vertix 2
    y2 : float
        y position of vertix 2
    x3 : float
        x position of vertix 3
    y3 : float
        y position of vertix 3
    
    Returns
    -------
    area : float
        area of triangle
    """

    return abs((x1 * (y2 - y3) + x2 * (y3 - y1)+ x3 * (y1 - y2))/ 2.0)

def isInside(x1, y1, x2, y2, x3, y3, x, y):
    """
    calculates if a point is inside a triangle or not

    Parameters
    ----------
    x1 : float
        x position of vertix 1
    y1 : float
        y position of vertix 1
    x2 : float
        x position of vertix 2
    y2 : float
        y position of vertix 2
    x3 : float
        x position of vertix 3
    y3 : float
        y position of vertix 3
    x : float
        x position of the point
    y : float
        y position of the point
    
    Returns
    -------
    isInside : boolean
        returns true if point is inside
    """

    A = area(x1, y1, x2, y2, x3, y3)
    A1 = area(x, y, x2, y2, x3, y3)
    A2 = area(x1, y1, x, y, x3, y3)
    A3 = area(x1, y1, x2, y2, x, y)
    AT = A1 + A2 + A3
    if(AT < A*1.00001): # To avoid floating point error
        return True
    else:
        return False

def force(domainNumber):
    """
    calculates force on a domain

    Parameters
    ----------
    domainNumber: int
        number of domain for which force is to be calculated

    Returns
    -------
    Force : list
        list of x and y components of vector Force
    """

    print("Get.force is under progress")
    domains = joblib.load("Data/initialDomains.sav")
    mesh = joblib.load("Data/mesh.sav")
    solver = joblib.load("Data/solver.sav")
    x = mesh.x
    y = mesh.y
    nodeSource = mesh.nodeSource
    nodes = np.concatenate((domains[domainNumber].boundaryList, domains[domainNumber].freeList)).astype(int)
    nodeForce = []
    for node in nodes:
        nodeVectorField = vectorField(np.array([x[node], y[node]]), True, node)
        if solver.name == "electrostatic":
            nodeForce.append(nodeVectorField*nodeSource[node])
        elif solver.name == "magnetostatic":
            print("turn magnetic force 90 degrees")
            nodeForce.append(nodeVectorField*nodeSource[node])
    Force = sum(nodeForce)
    return Force

