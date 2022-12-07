# DssPyLib

DssPyLib is a Python library for calculating 2-D integral and finite element numerical solutions for Poisson equation with simple non-overlapping shapes. DssPyLib provides support to apply integral Dirichelet boundary conditions for improved accuracy in open boundary problems. It also provides addition support to apply Distributed Source Scheme to reduce error due to presence of field sources. DssPyLib supports calculation of vector fields at any point inside the problem region and forces on the sources.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install DssPyLib.

```bash
pip install DssPyLib
```

## Usage

The working directory of any project using DssPyLib should contain main.py, Input.py, Output.py, and Get.py files which can be simply copied and edited from the Example projects.

### main.py

```python
# main.py
import Input
import Output
import Get
import numpy as np

def main():
    Input.Geometry(figSize = "Default", showPlot = False)

    Output.nodes(meshingMethod = "rand", figSize = "Default", showPlot = False)
    Output.mesh(meshNumbering = "numbered", figSize = "Default", showPlot = False)
    Output.relaxedMesh(meshNumbering = "numbered", figSize = "Default", showPlot = False)

    Input.Solver("electrostatic")
    Input.Sources(colorbarPosition = "right", figSize = "Default", showPlot = False)
    Input.Materials(colorbarPosition = "right", figSize = "Default", showPlot = False)

    Output.neighbourInformation(truncationNumber = 6)
    Output.solution(solutionName = "approximate", solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = False, withDifference = True)

    position = np.array([[0, 0]])
    field = Get.vectorField(position)
    print(field)
    domainNumber = 1
    force = Get.force(domainNumber)
    print(force)

if __name__ == "__main__":
    main()
```

### Input.py

```python
import joblib
import os
from DssPyLib.Domains.Rectangle import Rectangle
from DssPyLib.Domains.Ellipse import Ellipse
from DssPyLib.Domains.Line import Line
from DssPyLib.Domains.Point import Point

import DssPyLib.PlotUtils as PlotUtils
import Output

def Geometry(figSize = "Default", showPlot = False):
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
                                    10)) # non, Predefined no. of nodes

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

    if not os.path.exists('Data'):
        os.mkdir('Data')

    if not os.path.exists('Figures'):
        os.mkdir('Figures')

    joblib.dump(initialDomains, "Data/initialDomains.sav")
    PlotUtils.plotGeometry(initialDomains, "Figures/1_Geometry.png", figSize, showPlot)

def Sources(colorbarPosition = "right", figSize = "Default", showPlot = False):
    initialDomains = joblib.load("Data/initialDomains.sav")

    initialDomains[0].source = 0
    initialDomains[1].source = 1
    initialDomains[2].source = 0

    joblib.dump(initialDomains, "Data/initialDomains.sav")
    Output.assigedSources(colorbarPosition, figSize, showPlot)

def Materials(colorbarPosition = "right", figSize = "Default", showPlot = False):
    # Material property of points and lines will not be considered as material property is bulk property
    initialDomains = joblib.load("Data/initialDomains.sav")

    initialDomains[0].material = 1
    initialDomains[1].material = 1
    initialDomains[2].material = 1

    joblib.dump(initialDomains, "Data/initialDomains.sav")
    Output.assigedMaterials(colorbarPosition, figSize, showPlot)

def Solver(solverName):
    if (solverName == "electrostatic"):
        from DssPyLib.Solvers.Electrostatic import Electrostatic
        solver = Electrostatic()
    elif (solverName == "magnetostatic"):
        from DssPyLib.Solvers.Magnetostatic import Magnetostatic
        solver = Magnetostatic()
    joblib.dump(solver, "Data/solver.sav")

if __name__ == "__main__":
    Geometry()

```

Users are encouraged to modify Get.py to extract additional information from the saved data in the Data directory.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)