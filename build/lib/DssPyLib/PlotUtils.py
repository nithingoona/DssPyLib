from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'
import matplotlib.pyplot as plt

import matplotlib.tri as mtri
import numpy as np
from matplotlib.ticker import FormatStrFormatter

def plotGeometry(domains, filename, figSize, showPlot):
    """
    plots geometry.

    Parameters
    ----------
    domains : list
        list of geometries
    filename : str
        name of file
    figSize : list
        size of figure
    showPlot : boolean
        show plot if True

    Returns
    -------
    None
    """

    plt.close("all")
    fig, ax = plt.subplots(figsize=domains[0].figSize)
    w = domains[0].width
    h = domains[0].height
    [xp,yp] = domains[0].position
    dx = w/100
    sf = 1.5
    plt.xlim(xp-w/2-dx, xp+w/2+dx)
    plt.ylim(yp-h/2-dx, yp+h/2+dx)

    lW = figSize[0]/domains[0].defaultFigSize[0]
    for domain in domains:
        domain.show(plt, lW)

    plt.gca().set_aspect("equal")
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(figSize[0]/domains[0].defaultFigSize[0])

    plt.xlabel("x - axis [meter]", fontsize=domains[0].figSize[0]*1.5*sf)
    plt.ylabel("y - axis [meter]", fontsize=domains[0].figSize[0]*1.5*sf)
    plt.xticks(fontsize= domains[0].figSize[0]*sf)
    plt.yticks(fontsize= domains[0].figSize[0]*sf)
    plt.grid()
    plt.tick_params(axis='both', width=figSize[0]/domains[0].defaultFigSize[0])
    plt.tick_params(axis='both', length=figSize[0]/domains[0].defaultFigSize[0]*5)
    plt.savefig(filename)
    if showPlot == True:
        plt.show()

def plotMesh(domains, mesh, filename, figSize, showPlot, meshNumbering):
    """
    plots mesh.

    Parameters
    ----------
    domains : list
        list of geometries
    mesh : object
        mesh object
    filename : str
        name of file
    figSize : list
        size of figure
    showPlot : boolean
        show plot if True
    meshNumbering : boolean
        show node and element numbers if true

    Returns
    -------
    None
    """

    plt.close("all")
    if figSize == "Default":
        figSize = domains[0].figSize

    fig, ax = plt.subplots(figsize=figSize)
    w = domains[0].width
    h = domains[0].height
    [xp,yp] = domains[0].position
    dx = w/100
    sf = 1.5
    plt.xlim(xp-w/2-dx, xp+w/2+dx)
    plt.ylim(yp-h/2-dx, yp+h/2+dx)

    plt.triplot(mesh.x, mesh.y, mesh.tri, c = 'k', lw = figSize[0]/domains[0].defaultFigSize[0])

    if (meshNumbering == "unNumbered"):
        print("Skipping Element and Node numbering")
    elif (meshNumbering == "numbered"):
        # Display element and node numbering
        fs = figSize[0]/domains[0].non * 15
        ne = len(mesh.tri[:,0])
        nn = len(mesh.x)
        for i in range(ne):
            ex = (mesh.x[mesh.tri[i,0]] + mesh.x[mesh.tri[i,1]] + mesh.x[mesh.tri[i,2]])/3
            ey = (mesh.y[mesh.tri[i,0]] + mesh.y[mesh.tri[i,1]] + mesh.y[mesh.tri[i,2]])/3
            ex = ex/(domains[0].width+2*dx)
            ey = ey/(domains[0].height+2*dx)
            ex = ex + 0.5
            ey = ey + 0.5
            ax.text(ex, ey, f"%d" % i, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=fs, color='r')

        for i in range(nn):
            nx = mesh.x[i]
            ny = mesh.y[i]
            nx = nx/(domains[0].width+2*dx)
            ny = ny/(domains[0].height+2*dx)
            nx = nx + 0.5
            ny = ny + 0.5
            ax.text(nx, ny, f"%d" % i, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=fs, color='b')

    plt.gca().set_aspect("equal")
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(figSize[0]/domains[0].defaultFigSize[0])
    
    plt.xlabel("x - axis [meter]", fontsize=figSize[0]/domains[0].defaultFigSize[0]*1.5*20*sf)
    plt.ylabel("y - axis [meter]", fontsize=figSize[0]/domains[0].defaultFigSize[0]*1.5*20*sf)
    plt.xticks(fontsize= figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    plt.yticks(fontsize= figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    plt.grid()
    plt.tick_params(axis='both', width=figSize[0]/domains[0].defaultFigSize[0])
    plt.tick_params(axis='both', length=figSize[0]/domains[0].defaultFigSize[0]*5)
    plt.savefig(filename)
    if showPlot == True:
        plt.show()

def plotContour(domains, mesh, solver, plotType, filename, colorbarPosition, figSize, showPlot):
    """
    plots solution.

    Parameters
    ----------
    domains : list
        list of geometries
    mesh : object
        mesh object
    solver : object
        the solver object
    plotType : str
        type of plot
    filename : str
        name of file
    colorbarPosition : str
        position of color bar
    figSize : list
        size of figure
    showPlot : boolean
        show plot if True

    Returns
    -------
    None
    """

    plt.close("all")
    if figSize == "Default":
        figSize = domains[0].figSize

    fig, ax = plt.subplots(figsize=figSize)
    w = domains[0].width
    h = domains[0].height
    [xp,yp] = domains[0].position
    dx = w/100
    sf = 1.5
    plt.xlim(xp-w/2-dx, xp+w/2+dx)
    plt.ylim(yp-h/2-dx, yp+h/2+dx)

    figRatio = domains[0].defaultFigSize[1]/domains[0].defaultFigSize[0]

    if (plotType == "Source"):
        triang = mtri.Triangulation(mesh.x, mesh.y, mesh.tri)
        plt.tricontourf(triang, mesh.nodeSource[:,0], 990,  cmap = "viridis")
        if solver.name == "electrostatic":
            if colorbarPosition == "top" or colorbarPosition == "bottom":
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.044/figRatio, format='%.0e')
            else:
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.046*figRatio, format='%.0e')
            cbar.set_label(label="Charge Density [C/m]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf)
        elif solver.name == "magnetostatic":
            if colorbarPosition == "top" or colorbarPosition == "bottom":
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.044/figRatio, format='%.0e')
            else:
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.046*figRatio, format='%.0e')
            cbar.set_label(label="Electric Current [Amp]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf)
        else:
            print("plotContour does not support Source type of the specified solver ", solver.name)
        
    elif (plotType == "Material"):
        triang = mtri.Triangulation(mesh.x, mesh.y, mesh.tri)
        plt.tricontourf(triang, mesh.nodeMaterial[:,0], 990,  cmap = "viridis")
        if solver.name == "electrostatic":
            if colorbarPosition == "top" or colorbarPosition == "bottom":
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.044/figRatio, format='%.0e')
            else:
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.046*figRatio, format='%.0e')
            cbar.set_label(label="Relative Permittivity",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf)
        elif solver.name == "magnetostatic":
            if colorbarPosition == "top" or colorbarPosition == "bottom":
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.044/figRatio, format='%.0e')
            else:
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.046*figRatio, format='%.0e')
            cbar.set_label(label="Relative Permeability",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf)
        else:
            print("plotContour does not support Material type of the specified solver ", solver.name)

    elif (plotType == "InducedSources"):
        triang = mtri.Triangulation(mesh.x, mesh.y, mesh.tri)
        ep0 = solver.physicalConstant
        if solver.name == "electrostatic":
            plt.tricontourf(triang, solver.bi[:,0]*ep0, 990,  cmap = "viridis")
            if colorbarPosition == "top" or colorbarPosition == "bottom":
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.044/figRatio, format='%.0e')
            else:
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.046*figRatio, format='%.0e')
            cbar.set_label(label="Charge Density [C/m]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf)
        elif solver.name == "magnetostatic":
            plt.tricontourf(triang, solver.bi[:,0]/ep0, 990,  cmap = "viridis")
            if colorbarPosition == "top" or colorbarPosition == "bottom":
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.044/figRatio, format='%.0e')
            else:
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.046*figRatio, format='%.0e')
            cbar.set_label(label="Electric Current [Amp]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf)
        else:
            print("plotContour does not support Source type of the specified solver ", solver.name)

    elif (plotType == "SourceDSS"):
        triang = mtri.Triangulation(mesh.x, mesh.y, mesh.tri)
        plt.tricontourf(triang, mesh.nodeSourceDss[:,0]-mesh.nodeSource[:,0], 990,  cmap = "viridis")
        if solver.name == "electrostatic":
            if colorbarPosition == "top" or colorbarPosition == "bottom":
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.044/figRatio, format='%.0e')
            else:
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.046*figRatio, format='%.0e')
            cbar.set_label(label="Charge Density [C/m]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf)
        elif solver.name == "magnetostatic":
            if colorbarPosition == "top" or colorbarPosition == "bottom":
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.044/figRatio, format='%.0e')
            else:
                cbar = plt.colorbar(location = colorbarPosition, fraction=0.046*figRatio, format='%.0e')
            cbar.set_label(label="Electric Current [Amp]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf)
        else:
            print("plotContour does not support Source type of the specified solver ", solver.name)

    plt.gca().set_aspect("equal")
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(figSize[0]/domains[0].defaultFigSize[0])
    
    plt.xlabel("x - axis [meter]", fontsize=figSize[0]/domains[0].defaultFigSize[0]*1.5*20*sf)
    plt.ylabel("y - axis [meter]", fontsize=figSize[0]/domains[0].defaultFigSize[0]*1.5*20*sf)
    plt.xticks(fontsize= figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    plt.yticks(fontsize= figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    plt.grid()
    plt.tick_params(axis='both', width=figSize[0]/domains[0].defaultFigSize[0])
    plt.tick_params(axis='both', length=figSize[0]/domains[0].defaultFigSize[0]*5)
    cbar.ax.tick_params(labelsize=figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    plt.savefig(filename)
    if showPlot == True:
        plt.show()

def plotContourField(domains, filename, mesh, solver, solutionName, figSize, colorbarPosition, showPlot):
    """
    plots solution.

    Parameters
    ----------
    domains : list
        list of geometries
    filename : str
        name of file
    mesh : object
        mesh object
    solver : object
        the solver object
    solutionName : str
        name of solution
    figSize : list
        size of figure
    colorbarPosition : str
        position of color bar
    showPlot : boolean
        show plot if True

    Returns
    -------
    None
    """

    solution = solver.u
    x = mesh.rx
    y = mesh.ry
    plt.close("all")
    if figSize == "Default":
        figSize = domains[0].figSize

    fig, ax = plt.subplots(figsize=figSize)
    w = domains[0].width
    h = domains[0].height
    [xp,yp] = domains[0].position
    dx = w/100
    sf = 1.5
    plt.xlim(xp-w/2-dx, xp+w/2+dx)
    plt.ylim(yp-h/2-dx, yp+h/2+dx)
    figRatio = domains[0].defaultFigSize[1]/domains[0].defaultFigSize[0]

    triang = mtri.Triangulation(mesh.x, mesh.y, mesh.tri)
    plt.tricontourf(triang, solution[:,0]-solution[:,0][0], 990,  cmap =  "viridis")

    if colorbarPosition == "top" or colorbarPosition == "bottom":
        cbar = plt.colorbar(location = colorbarPosition, fraction=0.044/figRatio, format='%.0e')
    else:
        cbar = plt.colorbar(location = colorbarPosition, fraction=0.046*figRatio, format='%.0e')

    if solver.name == "electrostatic":
        cbar.set_label(label="Potential [Volt]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf)

    elif solver.name == "magnetostatic":
        cbar.set_label(label="Potential [Tesla meter]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf)

    plt.quiver(x, y, solver.Ex, solver.Ey)
    plt.gca().set_aspect("equal")
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(figSize[0]/domains[0].defaultFigSize[0])
    
    plt.xlabel("x - axis [meter]", fontsize=figSize[0]/domains[0].defaultFigSize[0]*1.5*20*sf)
    plt.ylabel("y - axis [meter]", fontsize=figSize[0]/domains[0].defaultFigSize[0]*1.5*20*sf)
    plt.xticks(fontsize= figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    plt.yticks(fontsize= figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    plt.tick_params(axis='both', width=figSize[0]/domains[0].defaultFigSize[0])
    plt.tick_params(axis='both', length=figSize[0]/domains[0].defaultFigSize[0]*5)
    cbar.ax.tick_params(labelsize=figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    plt.savefig(filename)
    if showPlot == True:
        plt.show()

def plotSurf(domains, filename, mesh, solver, solution, solutionName, figSize, showPlot):
    """
    plots solution.

    Parameters
    ----------
    domains : list
        list of geometries
    filename : str
        name of file
    mesh : object
        mesh object
    solver : object
        the solver object
    solution : list
        the solution
    solutionName : str
        name of solution
    figSize : list
        size of figure
    showPlot : boolean
        show plot if True

    Returns
    -------
    None
    """

    x = mesh.rx
    y = mesh.ry
    plt.close("all")
    if figSize == "Default":
        figSize = domains[0].figSize

    sf = 1.5
    fig = plt.figure(figsize=figSize)
    ax = plt.axes(projection='3d')
    ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(y)))
    ts = ax.plot_trisurf(x, y, np.transpose(solution)[0]-np.transpose(solution)[0][0], cmap = "viridis", edgecolor = 'grey', linewidth=0, antialiased=False)

    ax.set_xlabel("x - axis [meter]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf, labelpad=30*figSize[0]/domains[0].defaultFigSize[0])
    ax.set_ylabel("y - axis [meter]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf, labelpad=30*figSize[0]/domains[0].defaultFigSize[0])
    ax.zaxis.set_label_position("top")
    if solver.name == "electrostatic":
        ax.set_zlabel("Potential [Volt]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf, labelpad=30*figSize[0]/domains[0].defaultFigSize[0]) #, format='%.0e'

    elif solver.name == "magnetostatic":
        ax.set_zlabel("Potential [Tesla meter]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf, labelpad=30*figSize[0]/domains[0].defaultFigSize[0])

    ax.view_init(30, -110)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    plt.grid()
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.0e'))
    ax.xaxis.set_tick_params(labelsize=figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    ax.yaxis.set_tick_params(labelsize=figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    ax.zaxis.set_tick_params(labelsize=figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    plt.savefig(filename)
    if showPlot == True:
        plt.show()

def plotDifference(domains, filename, mesh, solver, appxsolution, solution, figSize, showPlot):
    """
    plots difference solution.

    Parameters
    ----------
    domains : list
        list of geometries
    filename : str
        name of file
    mesh : object
        mesh object
    solver : object
        the solver object
    appxsolution : list
        the approximate solution
    solution : list
        the solution
    solutionName : str
        name of solution
    figSize : list
        size of figure
    showPlot : boolean
        show plot if True

    Returns
    -------
    None
    """

    x = mesh.rx
    y = mesh.ry
    plt.close("all")
    if figSize == "Default":
        figSize = domains[0].figSize
    sf = 1.5
    fig = plt.figure(figsize=figSize)
    ax = plt.axes(projection='3d')
    ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(y)))
    difference = (np.transpose(appxsolution)[0]-np.transpose(appxsolution)[0][0])-(np.transpose(solution)[0]-np.transpose(solution)[0][0])
    ts = ax.plot_trisurf(x, y, difference, cmap = "viridis", edgecolor = 'grey', linewidth=0, antialiased=False)

    ax.set_xlabel("x - axis [meter]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf, labelpad=30*figSize[0]/domains[0].defaultFigSize[0])
    ax.set_ylabel("y - axis [meter]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf, labelpad=30*figSize[0]/domains[0].defaultFigSize[0])

    if solver.name == "electrostatic":
        ax.set_zlabel("Potential [Volt]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf, labelpad=30*figSize[0]/domains[0].defaultFigSize[0])

    elif solver.name == "magnetostatic":
        ax.set_zlabel("Potential [Tesla meter]",size=figSize[0]/domains[0].defaultFigSize[0]*30*sf, labelpad=30*figSize[0]/domains[0].defaultFigSize[0])

    ax.view_init(30, -110)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    plt.grid()
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.0e'))
    ax.xaxis.set_tick_params(labelsize=figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    ax.yaxis.set_tick_params(labelsize=figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    ax.zaxis.set_tick_params(labelsize=figSize[0]/domains[0].defaultFigSize[0]*20*sf)
    plt.savefig(filename)
    if showPlot == True:
        plt.show()
