# Import packages
import joblib
import os
from DssPyLib.Domains.Rectangle import Rectangle
from DssPyLib.Domains.Ellipse import Ellipse
from DssPyLib.Domains.Line import Line
from DssPyLib.Domains.Point import Point

import DssPyLib.PlotUtils as PlotUtils
import Output

def Geometry(figSize = "Default", showPlot = False):
    """
    user input geometry.

    Parameters
    ----------
    figSize : list
        size of figure
    showPlot : boolean
        show plot if true

    Returns
    -------
    None
    """

    initialDomains = []
    # The first domain with index 0 is always a Global domain and should be initiated with non

    # Domain 0
    initialDomains.append(Rectangle(2, # Width
                                    1, # Height
                                    [0, 0], # Position
                                    [1, 1], # Relative Mesh Density
                                    0, # Angle
                                    False, # is_child
                                    True, # is_parent
                                    [2], # Child domains
                                    [], # Parent domains
                                    50)) # non, Predefined no. of nodes

    # Domain 1
    initialDomains.append(Point(0, # Width
                                0, # Height
                                [0, 0], # Position
                                [1, 1], # Relative Mesh Density
                                0, # Angle
                                True, # is_child
                                False, # is_parent
                                [], # Child domains
                                [2], # Parent domains
                                0)) # Predefined no. of nodes

    # Domain 2
    h = initialDomains[0].width/initialDomains[0].non
    initialDomains.append(Rectangle(initialDomains[0].width - 2*h, # Width
                                initialDomains[0].height - 2*h, # Height
                                initialDomains[0].position, # Position
                                [1.2, 1.2], # Relative Mesh Density
                                0, # Angle
                                True, # is_child
                                True, # is_parent
                                [1], # Child domains
                                [0], # Parent domains
                                0)) # Predefined no. of nodes

    if figSize == "Default":
        figSize = (20, 20 * initialDomains[0].height/ initialDomains[0].width)
    
    initialDomains[0].defaultFigSize = (20, 20 * initialDomains[0].height/ initialDomains[0].width)
    initialDomains[0].figSize = figSize
    print("The shape of global domain is", initialDomains[0].shape, "(From input.py)")

    if not os.path.exists('Data'):
        os.mkdir('Data')

    if not os.path.exists('Figures'):
        os.mkdir('Figures')

    joblib.dump(initialDomains, "Data/initialDomains.sav")
    PlotUtils.plotGeometry(initialDomains, "Figures/1_Geometry.png", figSize, showPlot)

def Sources(colorbarPosition = "right", figSize = "Default", showPlot = False):
    """
    user input sources.

    Parameters
    ----------
    colorbarPosition : str
        position of color bar
    figSize : list
        size of figure
    showPlot : boolean
        show plot if true

    Returns
    -------
    None
    """

    initialDomains = joblib.load("Data/initialDomains.sav")

    initialDomains[0].source = 0
    initialDomains[1].source = 1
    initialDomains[2].source = 0

    joblib.dump(initialDomains, "Data/initialDomains.sav")
    Output.assigedSources(colorbarPosition, figSize, showPlot)

def Materials(colorbarPosition = "right", figSize = "Default", showPlot = False):
    """
    user input materials.

    Parameters
    ----------
    colorbarPosition : str
        position of color bar
    figSize : list
        size of figure
    showPlot : boolean
        show plot if true

    Returns
    -------
    None
    """

    # Material property of points and lines will not be considered as material property is bulk property 
    initialDomains = joblib.load("Data/initialDomains.sav")

    initialDomains[0].material = 1
    initialDomains[1].material = 1
    initialDomains[2].material = 1

    joblib.dump(initialDomains, "Data/initialDomains.sav")
    Output.assigedMaterials(colorbarPosition, figSize, showPlot)

def Solver(solverName):
    """
    user input solver.

    Parameters
    ----------
    solverName : str
        name of the solver

    Returns
    -------
    None
    """

    if (solverName == "electrostatic"):
        from DssPyLib.Solvers.Electrostatic import Electrostatic
        solver = Electrostatic()
    elif (solverName == "magnetostatic"):
        from DssPyLib.Solvers.Magnetostatic import Magnetostatic
        solver = Magnetostatic()
    joblib.dump(solver, "Data/solver.sav")

if __name__ == "__main__":
    Geometry()