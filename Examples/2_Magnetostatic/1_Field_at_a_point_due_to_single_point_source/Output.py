# Import packages
import joblib

import DssPyLib.Mesh as Mesh
import DssPyLib.PlotUtils as PlotUtils
import DssPyLib.Solvers.Solve as Solve


def nodes(meshingMethod = "tri", figSize = "Default", showPlot = False):
    """
    generates nodes.

    Parameters
    ----------
    meshingMethod : str
        method of meshing
    figSize : list
        size of figure
    showPlot : boolean
        show plot if true

    Returns
    -------
    None
    """

    initialDomains = joblib.load("Data/initialDomains.sav")
    Mesh.generateNodes(initialDomains, meshingMethod, figSize, showPlot)
    joblib.dump(initialDomains, "Data/initialDomains.sav")
    
def mesh(meshNumbering = "unNumbered", figSize = "Default", showPlot = False):
    """
    generates mesh.

    Parameters
    ----------
    meshNumbering : str
        string mentioning numbered or unNumbered mesh
    figSize : list
        size of figure
    showPlot : boolean
        show plot if true

    Returns
    -------
    None
    """

    initialDomains = joblib.load("Data/initialDomains.sav")
    mesh = Mesh.createMesh(initialDomains)
    joblib.dump(initialDomains, "Data/initialDomains.sav")
    joblib.dump(mesh, "Data/mesh.sav")
    PlotUtils.plotMesh(initialDomains, mesh, "Figures/2_Initial_mesh.png", figSize, showPlot, meshNumbering)

def relaxedMesh(meshNumbering = "unNumbered", figSize = "Default", showPlot = False):
    """
    relaxes mesh.

    Parameters
    ----------
    meshNumbering : str
        string mentioning numbered or unNumbered mesh
    figSize : list
        size of figure
    showPlot : boolean
        show plot if true

    Returns
    -------
    None
    """

    initialDomains = joblib.load("Data/initialDomains.sav")
    mesh = joblib.load("Data/mesh.sav")
    Mesh.relaxMesh(initialDomains, mesh)
    joblib.dump(initialDomains, "Data/initialDomains.sav")
    joblib.dump(mesh, "Data/mesh.sav")
    PlotUtils.plotMesh(initialDomains, mesh, "Figures/2_Relaxed_mesh.png", figSize, showPlot, meshNumbering)

def assigedSources(colorbarPosition, figSize, showPlot):
    """
    assigns sources.

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
    mesh = joblib.load("Data/mesh.sav")
    solver = joblib.load("Data/solver.sav")
    Mesh.assignSources(initialDomains, mesh, solver)
    joblib.dump(mesh, "Data/mesh.sav")
    PlotUtils.plotContour(initialDomains, mesh, solver, "Source", "Figures/3_Sources.png", colorbarPosition, figSize, showPlot)

def assigedMaterials(colorbarPosition, figSize, showPlot):
    """
    assigns materials.

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
    mesh = joblib.load("Data/mesh.sav")
    solver = joblib.load("Data/solver.sav")
    Mesh.assignMaterials(initialDomains, mesh)
    joblib.dump(mesh, "Data/mesh.sav")
    PlotUtils.plotContour(initialDomains, mesh, solver, "Material", "Figures/3_Materials.png", colorbarPosition, figSize, showPlot)

def neighbourInformation(truncationNumber = 5):
    """
    generates neighbour information.

    Parameters
    ----------
    truncationNumber : int
        truncation number    

    Returns
    -------
    None
    """

    mesh = joblib.load("Data/mesh.sav")
    Mesh.getNeighbourInformation(mesh, truncationNumber)
    joblib.dump(mesh, "Data/mesh.sav")

def solution(solutionName = "numerical", solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = False, withDifference = False):
    """
    generates solution.

    Parameters
    ----------
    solutionName : str
        name of solution
    solverMethod : str
        name of solver
    figSize : list
        size of figure
    colorbarPosition : str
        position of colorbar
    showPlot : boolean
        show plot if true
    withField : boolean
        show vector field if true
    withDifference : boolean
        show difference in solution if true

    Returns
    -------
    None
    """

    initialDomains = joblib.load("Data/initialDomains.sav")
    mesh = joblib.load("Data/mesh.sav")
    solver = joblib.load("Data/solver.sav")
    Solve.getMatrices(initialDomains, mesh, solver, solutionName)
    Solve.solve(initialDomains, mesh, solver, solutionName, solverMethod)
    showSolution(mesh, solver, solutionName, figSize, colorbarPosition, showPlot, withField)
    if withDifference:
        showDifference(mesh, solver, solutionName, figSize, showPlot)
    joblib.dump(mesh, "Data/mesh.sav")
    joblib.dump(solver, "Data/solver.sav")

def dssSolution(solutionName = "numerical", truncationNumber = 5, solverMethod = "Scipy", figSize = "Default", colorbarPosition = "right", showPlot = False, withField = False, withDifference = False):
    """
    generates dss solution.

    Parameters
    ----------
    solutionName : str
        name of solution
    truncationNumber : int
        truncation number 
    solverMethod : str
        name of solver
    figSize : list
        size of figure
    colorbarPosition : str
        position of colorbar
    showPlot : boolean
        show plot if true
    withField : boolean
        show vector field if true
    withDifference : boolean
        show difference in solution if true

    Returns
    -------
    None
    """

    initialDomains = joblib.load("Data/initialDomains.sav")
    mesh = joblib.load("Data/mesh.sav")
    solver = joblib.load("Data/solver.sav")
    Solve.getDistributedSources(initialDomains, mesh, solver, truncationNumber)
    Solve.getMatrices(initialDomains, mesh, solver, solutionName)
    Solve.solveDSS(initialDomains, mesh, solver, solutionName, solverMethod)
    showDSSSolution(mesh, solver, solutionName, figSize, colorbarPosition, showPlot, withField)
    showDistributedSources(initialDomains, mesh, solver, solutionName, colorbarPosition, figSize, showPlot)
    if withDifference:
        showDSSDifference(mesh, solver, solutionName, figSize, showPlot)
    joblib.dump(mesh, "Data/mesh.sav")
    joblib.dump(solver, "Data/solver.sav")

def showSolution(mesh, solver, solutionName, figSize, colorbarPosition, showPlot, withField):
    """
    shows solution.

    Parameters
    ----------
    mesh : object
        mesh object
    solver : object
        solver object
    solutionName : str
        name of solution
    figSize : list
        size of figure
    colorbarPosition : str
        position of colorbar
    showPlot : boolean
        show plot if true
    withField : boolean
        show vector field if true

    Returns
    -------
    None
    """

    domains = joblib.load("Data/initialDomains.sav")
    if (solutionName == "approximate"):
        if (withField):
            PlotUtils.plotContourField(domains, "Figures/4_Approximate_solution", mesh, solver, solutionName, figSize, colorbarPosition, showPlot)
        else:
            solution = solver.ua
            PlotUtils.plotSurf(domains, "Figures/4_Approximate_solution", mesh, solver, solution, solutionName, figSize, showPlot)
    elif (solutionName == "numerical"):
        if (withField):
            PlotUtils.plotContourField(domains, "Figures/4_Numerical_solution", mesh, solver, solutionName, figSize, colorbarPosition, showPlot)
        else:
            solution = solver.un
            PlotUtils.plotSurf(domains, "Figures/4_Numerical_solution", mesh, solver, solution, solutionName, figSize, showPlot)
    elif (solutionName == "numericalApproxBound"):
        if (withField):
            PlotUtils.plotContourField(domains, "Figures/4_Numerical_solution_approx_boundary", mesh, solver, solutionName, figSize, colorbarPosition, showPlot)
        else:
            solution = solver.unab
            PlotUtils.plotSurf(domains, "Figures/4_Numerical_solution_approx_boundary", mesh, solver, solution, solutionName, figSize, showPlot)
    elif (solutionName == "numericalIterBound"):
        if (withField):
            PlotUtils.plotContourField(domains, "Figures/4_Numerical_solution_iter_boundary", mesh, solver, solutionName, figSize, colorbarPosition, showPlot)
        else:
            solution = solver.unib
            PlotUtils.plotSurf(domains, "Figures/4_Numerical_solution_iter_boundary", mesh, solver, solution, solutionName, figSize, showPlot)

def showDSSSolution(mesh, solver, solutionName, figSize, colorbarPosition, showPlot, withField):
    """
    shows dss solution.

    Parameters
    ----------
    mesh : object
        mesh object
    solver : object
        solver object
    solutionName : str
        name of solution
    figSize : list
        size of figure
    colorbarPosition : str
        position of colorbar
    showPlot : boolean
        show plot if true
    withField : boolean
        show vector field if true

    Returns
    -------
    None
    """

    domains = joblib.load("Data/initialDomains.sav")
    if (solutionName == "approximate"):
        print("DSS does not support approximate solution")
    elif (solutionName == "numerical"):
        if (withField):
            PlotUtils.plotContourField(domains, "Figures/4_Numerical_DSS_solution", mesh, solver, solutionName, figSize, colorbarPosition, showPlot)
        else:
            solution = solver.un
            PlotUtils.plotSurf(domains, "Figures/4_Numerical_DSS_solution", mesh, solver, solution, solutionName, figSize, showPlot)
    elif (solutionName == "numericalApproxBound"):
        if (withField):
            PlotUtils.plotContourField(domains, "Figures/4_Numerical_DSS_solution_approx_boundary", mesh, solver, solutionName, figSize, colorbarPosition, showPlot)
        else:
            solution = solver.unab
            PlotUtils.plotSurf(domains, "Figures/4_Numerical_DSS_solution_approx_boundary", mesh, solver, solution, solutionName, figSize, showPlot)
    elif (solutionName == "numericalIterBound"):
        if (withField):
            PlotUtils.plotContourField(domains, "Figures/4_Numerical_DSS_solution_iter_boundary", mesh, solver, solutionName, figSize, colorbarPosition, showPlot)
        else:
            solution = solver.unib
            PlotUtils.plotSurf(domains, "Figures/4_Numerical_DSS_solution_iter_boundary", mesh, solver, solution, solutionName, figSize, showPlot)

def showDifference(mesh, solver, solutionName, figSize, showPlot):
    """
    shows difference in solution.

    Parameters
    ----------
    mesh : object
        mesh object
    solver : object
        solver object
    solutionName : str
        name of solution
    figSize : list
        size of figure
    showPlot : boolean
        show plot if true

    Returns
    -------
    None
    """

    domains = joblib.load("Data/initialDomains.sav")
    if (solutionName == "approximate"):
        print("Difference is not supported for approximate solution")
    elif (solutionName == "numerical"):
        appxsolution = solver.ua
        solution = solver.un
        PlotUtils.plotDifference(domains, "Figures/5_Difference_in_Numerical_solution", mesh, solver, appxsolution, solution, figSize, showPlot)
    elif (solutionName == "numericalApproxBound"):
        appxsolution = solver.ua
        solution = solver.unab
        PlotUtils.plotDifference(domains, "Figures/5_Difference_in_Numerical_solution_approx_boundary", mesh, solver, appxsolution, solution, figSize, showPlot)
    elif (solutionName == "numericalIterBound"):
        appxsolution = solver.ua
        solution = solver.unib
        PlotUtils.plotDifference(domains, "Figures/5_Difference_in_Numerical_solution_iter_boundary", mesh, solver, appxsolution, solution, figSize, showPlot)

def showDSSDifference(mesh, solver, solutionName, figSize, showPlot):
    """
    shows difference in dss solution.

    Parameters
    ----------
    mesh : object
        mesh object
    solver : object
        solver object
    solutionName : str
        name of solution
    figSize : list
        size of figure
    showPlot : boolean
        show plot if true

    Returns
    -------
    None
    """

    domains = joblib.load("Data/initialDomains.sav")
    if (solutionName == "approximate"):
        print("Difference is not supported for approximate solution")
    elif (solutionName == "numerical"):
        appxsolution = solver.ua
        solution = solver.un
        PlotUtils.plotDifference(domains, "Figures/5_Difference_in_Numerical_DSS_solution", mesh, solver, appxsolution, solution, figSize, showPlot)
    elif (solutionName == "numericalApproxBound"):
        appxsolution = solver.ua
        solution = solver.unab
        PlotUtils.plotDifference(domains, "Figures/5_Difference_in_Numerical_DSS_solution_approx_boundary", mesh, solver, appxsolution, solution, figSize, showPlot)
    elif (solutionName == "numericalIterBound"):
        appxsolution = solver.ua
        solution = solver.unib
        PlotUtils.plotDifference(domains, "Figures/5_Difference_in_Numerical_DSS_solution_iter_boundary", mesh, solver, appxsolution, solution, figSize, showPlot)

def showDistributedSources(initialDomains, mesh, solver, solutionName, colorbarPosition, figSize, showPlot):
    """
    shows distributed sources.

    Parameters
    ----------
    initialDomains : list
        list of geometry objects
    mesh : object
        mesh object
    solver : object
        solver object
    solutionName : str
        name of solution
    colorbarPosition : str
        position of colorbar
    figSize : list
        size of figure
    showPlot : boolean
        show plot if true

    Returns
    -------
    None
    """
    
    if (solutionName == "approximate"):
        print("Difference not support approximate solution")
    else:
        PlotUtils.plotContour(initialDomains, mesh, solver, "SourceDSS", "Figures/6_Difference_in_Distributed_sources.png", colorbarPosition, figSize, showPlot)

