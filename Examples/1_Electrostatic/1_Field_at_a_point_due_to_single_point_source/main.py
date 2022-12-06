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

    tht = range(0, 360, 20)
    rad = 0.3
    ep0 = 8.8541878128*10**(-12)

    E_magnitude_analytical = 1/(2*np.pi*ep0*rad)
    print("Magnitude of Analytical Electric FIeld is", E_magnitude_analytical)
    for th in tht:
        x = rad*np.cos(np.deg2rad(th))
        y = rad*np.sin(np.deg2rad(th))
        position = np.array([[x, y]])
        field = Get.vectorField(position)
        print(field, (np.linalg.norm(field)-E_magnitude_analytical)/E_magnitude_analytical*100)

    Output.solution(solutionName = "numerical", solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = False, withDifference = True)

    for th in tht:
        x = rad*np.cos(np.deg2rad(th))
        y = rad*np.sin(np.deg2rad(th))
        position = np.array([[x, y]])
        field = Get.vectorField(position)
        print(field, (np.linalg.norm(field)-E_magnitude_analytical)/E_magnitude_analytical*100)

    Output.solution(solutionName = "numericalApproxBound", solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = False, withDifference = True)

    for th in tht:
        x = rad*np.cos(np.deg2rad(th))
        y = rad*np.sin(np.deg2rad(th))
        position = np.array([[x, y]])
        field = Get.vectorField(position)
        print(field, (np.linalg.norm(field)-E_magnitude_analytical)/E_magnitude_analytical*100)

    Output.solution(solutionName = "numericalIterBound", solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = False, withDifference = True)

    for th in tht:
        x = rad*np.cos(np.deg2rad(th))
        y = rad*np.sin(np.deg2rad(th))
        position = np.array([[x, y]])
        field = Get.vectorField(position)
        print(field, (np.linalg.norm(field)-E_magnitude_analytical)/E_magnitude_analytical*100)

    Output.dssSolution(solutionName = "numerical", truncationNumber = 5, solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = True, withDifference = True)

    for th in tht:
        x = rad*np.cos(np.deg2rad(th))
        y = rad*np.sin(np.deg2rad(th))
        position = np.array([[x, y]])
        field = Get.vectorField(position)
        print(field, (np.linalg.norm(field)-E_magnitude_analytical)/E_magnitude_analytical*100)

    Output.dssSolution(solutionName = "numericalApproxBound", truncationNumber = 5, solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = False, withDifference = True)

    for th in tht:
        x = rad*np.cos(np.deg2rad(th))
        y = rad*np.sin(np.deg2rad(th))
        position = np.array([[x, y]])
        field = Get.vectorField(position)
        print(field, (np.linalg.norm(field)-E_magnitude_analytical)/E_magnitude_analytical*100)

    Output.dssSolution(solutionName = "numericalIterBound", truncationNumber = 5, solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = False, withDifference = True)

    for th in tht:
        x = rad*np.cos(np.deg2rad(th))
        y = rad*np.sin(np.deg2rad(th))
        position = np.array([[x, y]])
        field = Get.vectorField(position)
        print(field, (np.linalg.norm(field)-E_magnitude_analytical)/E_magnitude_analytical*100)


if __name__ == "__main__":
    main()