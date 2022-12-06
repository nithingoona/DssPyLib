# DssPyLib

DssPyLib is a Python library for calculating 2-D integral and finite element numerical solutions for Poisson equation with simple non-overlapping shapes. DssPyLib provides support to apply integral Dirichelet boundary conditions for improved accuracy in open boundary problems. It also provides addition support to apply Distributed Source Scheme to reduce error due to presence of field sources. DssPyLib supports calculation of vector fields at any point inside the problem region and forces on the sources.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install DssPyLib.

```bash
pip install DssPyLib
```

## Usage

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

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)