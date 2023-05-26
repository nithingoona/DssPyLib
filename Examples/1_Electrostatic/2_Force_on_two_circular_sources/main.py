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

    actualForce = 2*9*10**9
    print(actualForce)
    f1 = Get.force(1)
    f2 = Get.force(2)
    averageForce = (f1[0]-f2[0])/2
    error = (averageForce-actualForce)/actualForce*100
    print(f1)
    print(f2)
    print(error)

    Output.solution(solutionName = "numerical", solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = False, withDifference = True)

    print(actualForce)
    f1 = Get.force(1)
    f2 = Get.force(2)
    averageForce = (f1[0]-f2[0])/2
    error = (averageForce-actualForce)/actualForce*100
    print(f1)
    print(f2)
    print(error)

    Output.solution(solutionName = "numericalApproxBound", solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = False, withDifference = True)

    print(actualForce)
    f1 = Get.force(1)
    f2 = Get.force(2)
    averageForce = (f1[0]-f2[0])/2
    error = (averageForce-actualForce)/actualForce*100
    print(f1)
    print(f2)
    print(error)

    Output.dssSolution(solutionName = "numerical", truncationNumber = 5, solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = True, withDifference = True)

    print(actualForce)
    f1 = Get.force(1)
    f2 = Get.force(2)
    averageForce = (f1[0]-f2[0])/2
    error = (averageForce-actualForce)/actualForce*100
    print(f1)
    print(f2)
    print(error)

    Output.dssSolution(solutionName = "numericalApproxBound", truncationNumber = 5, solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = False, withDifference = True)

    print(actualForce)
    f1 = Get.force(1)
    f2 = Get.force(2)
    averageForce = (f1[0]-f2[0])/2
    error = (averageForce-actualForce)/actualForce*100
    print(f1)
    print(f2)
    print(error)


if __name__ == "__main__":
    main()