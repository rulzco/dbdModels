import numpy as np
# from scipy import optimize as op
from scipy.optimize import fsolve

x_min = -0.02    # Domain length in x-axis min
x_max = 0.7     # Domain length in x-axis max
y_min = -0.03   # Domain length in y-axis min
y_max = 0.14    # Domain length in y-axis min

L_expE = 0.01            # Length of exposed electrode
t_expE = 6.0e-5         # Thickness of exposed electrode
L_grdE = 0.01            # Length of grounded electrode
t_grdE = 6.0e-5         # Thickness of grounded electrode
t_Diel = 1.1e-4        # Dielectric between electrodes thickness
L_Gap = 5.0e-4           # Gap between electrodes


cellToCellRatio = 1.05

cellsYDielectric = 18   # Cells y-axis thickness of deielectric
cellsYGrounded = 9     # Cells y-axis thickness of grounded electrode
cellsYExposed = 9       # Cells y-axis thickness of exposed electrode

cellsYTop = 100          # Cells y-axis exposed up to "Top" boundary
cellsYBottom = 50       # Cells y-axis grounded down to "Bottom" boundary

cellsXGap = 15          # Cells x-axis gap (min dx=L/cellsXGap)
cellsXExposed = 50      # Cells x-axis exposed electrode
CellsXGrounded = 50     # Cells x-axis grounded electrode
cellsXRHS = 100          # Cells x-axis RHS of grounded
cellsXLHS = 9          # Cells x-axis LHS of exposed

# Compute the total expansion ratio


def RatioElectrodes():

    dx0 = L_Gap / cellsXGap  # initial ratio

    def f(x):
        return x**cellsXExposed - L_expE / dx0 * x + \
            L_expE / dx0 - 1

    def g(x):
        return x**CellsXGrounded - L_grdE / dx0 * x + \
            L_grdE / dx0 - 1
# Solve for the total exp ratio
    ratiosExp = fsolve(f, [1, 100], maxfev=1000)
    ratiosGrd = fsolve(g, [1, 100], maxfev=1000)

# Check roots convergence for correct total exp ratio
    for i in range(len(ratiosExp)):
        if ratiosExp[i] != 1.0:
            tempLen = (dx0 * (1 - ratiosExp[i]**cellsXExposed)
                       / (1 - ratiosExp[i]))
            cellToCellRatioExp = ratiosExp[i]
            if str(np.isclose(L_expE, tempLen)) == str('True'):
                print('\nCell to cell ratio x-axis exposed: ',
                      '{:6.8f}'.format(cellToCellRatioExp))

    for i in range(len(ratiosGrd)):
        if ratiosGrd[i] != 1.0:
            tempLen = (dx0 * (1 - ratiosGrd[i]**CellsXGrounded)
                       / (1 - ratiosGrd[i]))
            cellToCellRatioGrd = ratiosGrd[i]
            if str(np.isclose(L_grdE, tempLen)) == str('True'):
                print('Cell to cell ratio x-axis grounded: ',
                      '{:6.8f}'.format(cellToCellRatioGrd))

    dxNExposed = dx0 * cellToCellRatioExp ** (cellsXExposed-1)
    dxNGrounded = dx0 * cellToCellRatioGrd ** (CellsXGrounded-1)
    totRatioE = dxNExposed / dx0
    totRatioG = dxNGrounded / dx0

    print('Delta x0', '{:6.8f}'.format(dx0))
    print('Delta-xN exposed: ', '{:6.8f}'.format(dxNExposed))
    print('Delta-xN grounded: ', '{:6.8f}'.format(dxNGrounded))
    print('Total ratio x-axis grounded: ', '{:6.8f}'.format(totRatioG))
    print('Total ratio x-axis exposed: ', '{:6.8f}'.format(1 / totRatioE))
    print('\n')
    return dxNExposed, dxNGrounded


def RatioX(cellToCellRatio):

    dx0 = L_Gap / cellsXGap  # initial ratio
    # dx0RHS = L_expE / cellsXExposed
    # dx0LHS = L_expE / CellsXGrounded
    dxNExposed, dxNGrounded = RatioElectrodes()
    dx0RHS = dxNGrounded
    dx0LHS = dxNExposed

    xLenghtRHS = (x_max - L_grdE - L_Gap)
    xLenghtLHS = abs(x_min + L_expE)

    def f(x):
        return x**cellsXRHS - xLenghtRHS / dx0RHS * x + \
            xLenghtRHS / dx0RHS - 1

    def g(x):
        return x**cellsXLHS - xLenghtLHS / dx0LHS * x + \
            xLenghtLHS / dx0LHS - 1

# Solve for the total exp ratio
    xRHS = fsolve(f, [1, 100], maxfev=1000)
    xLHS = fsolve(g, [1, 100], maxfev=1000)

# Check roots convergence for correct total exp ratio
    for i in range(len(xRHS)):
        if xRHS[i] != 1.0:
            tempLen = (dx0RHS * (1 - xRHS[i]**cellsXRHS)
                       / (1 - xRHS[i]))
            cellToCellRatioRHS = xRHS[i]
            if str(np.isclose(xLenghtRHS, tempLen)) == str('True'):
                print('\nCell to cell ratio x-axis RHS: ',
                      '{:6.8f}'.format(xRHS[i]))

    for i in range(len(xLHS)):
        if xLHS[i] != 1.0:
            tempLen = (dx0LHS * (1 - xLHS[i]**cellsXLHS)
                       / (1 - xLHS[i]))
            cellToCellRatioLHS = xLHS[i]
            if str(np.isclose(xLenghtLHS, tempLen)) == str('True'):
                print('Cell to cell ratio x-axis LHS: ',
                      '{:6.8f}'.format(xLHS[i]))

    dxNrhs = dx0RHS * cellToCellRatioRHS ** (cellsXRHS-1)
    dxNlhs = dx0LHS * cellToCellRatioLHS ** (cellsXLHS-1)
    totRatioRHS = dxNrhs / dx0RHS
    totRatioLHS = dxNlhs / dx0LHS

    print('Delta-x Gap: ', '{:6.8f}'.format(dx0))
    print('Delta-x exposed: ', '{:6.8f}'.format(L_expE / cellsXExposed))
    print('Delta-x grounded: ', '{:6.8f}'.format(L_grdE / CellsXGrounded))
    print('Delta-xN RHS: ', '{:6.8f}'.format(dxNrhs))
    print('Delta-xN LHS: ', '{:6.8f}'.format(dxNlhs))
    print('Total ratio x-axis RHS: ', '{:6.8f}'.format(totRatioRHS))
    print('Total ratio x-axis LHS: ', '{:6.8f}'.format(1 / totRatioLHS))

    print('RHS real length: ', '{:6.8f}'.format(xLenghtRHS))
    print('LHS real length: ', '{:6.8f}'.format(xLenghtLHS))
    print('\n')


def RatioY(cellToCellRatio):

    dy0top = t_expE / cellsYExposed
    dy0bot = t_grdE / cellsYGrounded

    yLenghtTop = (y_max - t_expE)
    yLenghtBot = abs(y_min + t_Diel + t_grdE)

    def f(x):
        return x**cellsYTop - yLenghtTop / dy0top * x + \
            yLenghtTop / dy0top - 1

    def g(x):
        return x**cellsYBottom - yLenghtBot / dy0bot * x + \
            yLenghtBot / dy0bot - 1

    yTop = fsolve(f, [1, 100], maxfev=1000)
    yBot = fsolve(g, [1, 100], maxfev=1000)

    for i in range(len(yTop)):
        if yTop[i] != 1.0:
            tempLen = (dy0top * (1 - yTop[i]**cellsYTop)
                       / (1 - yTop[i]))
            cellToCellRatioTop = yTop[i]
            if str(np.isclose(yLenghtTop, tempLen)) == str('True'):
                print('\nCell to cell ratio y-axis Top: ',
                      '{:6.8f}'.format(yTop[i]))

    for i in range(len(yBot)):
        if yBot[i] != 1.0:
            tempLen = (dy0bot * (1 - yBot[i]**cellsYBottom)
                       / (1 - yBot[i]))
            cellToCellRatioBot = yBot[i]
            if str(np.isclose(yLenghtBot, tempLen)) == str('True'):
                print('Cell to cell ratio y-axis Bot: ',
                      '{:6.8f}'.format(yBot[i]))

    dyNTop = dy0top * cellToCellRatioTop ** (cellsYTop-1)
    dyNBot = dy0bot * cellToCellRatioBot ** (cellsYBottom-1)

    totRatioTop = dyNTop / dy0top
    totRatioBot = dyNBot / dy0bot

    print('Delta-y dielectric: ',
          '{:6.8f}'.format(t_Diel / cellsYDielectric))
    print('Delta-y exposed: ', '{:6.8f}'.format(dy0top))
    print('Delta-y grounded: ', '{:6.8f}'.format(dy0bot))
    print('Delta-yN Top: ', '{:6.8f}'.format(dyNTop))
    print('Delta-yN Bot: ', '{:6.8f}'.format(dyNBot))
    print('Total ratio y-axis exposed up to Top: ',
          '{:6.8f}'.format(totRatioTop))
    print('Total ratio y-axis grounded to Bottom:',
          '{:6.8f}'.format(1 / totRatioBot))

    print('Top real length: ', '{:6.8f}'.format(yLenghtTop))
    print('Bottom real length: ', '{:6.8f}'.format(yLenghtBot))
# compute the number of cells bases on cell-to-cell ratio


def cells(cellToCellRatio):

    print('\n')
    dx0_gap = float(input('Delta-x gap: '))
    dx0_exposed = float(input('Delta-x exposed: '))
    dx0_grounded = float(input('Delta-x grounded: '))

    print('\n')
    dy0_gap = float(input('Delta-y gap: '))
    dy0_exposed = float(input('Delta-y exposed: '))
    dy0_grounded = float(input('Delta-y grounded: '))

    cellsXGap = dx0_gap * L_Gap
    CellsXGrounded = dx0_grounded * L_grdE
    cellsXExposed = dx0_exposed * L_expE
    cellsXRHS = np.log10((x_max - L_Gap - L_grdE) / dx0_grounded *
                         (cellToCellRatio - 1) + 1) / np.log10(cellToCellRatio)
    cellsXLHS = np.log10(np.abs(x_min + L_expE) / dx0_exposed *
                         (cellToCellRatio - 1) + 1) / np.log10(cellToCellRatio)

    cellsYDielectric = dy0_gap * t_Diel
    CellsYGrounded = dy0_grounded * t_grdE
    cellsYExposed = dy0_exposed * t_expE
    cellsYTop = np.log10((y_max - t_expE) / dy0_exposed *
                         (cellToCellRatio - 1) + 1) / np.log10(cellToCellRatio)
    cellsYBottom = (np.log10(np.abs(y_min + t_Diel + t_grdE) / dy0_grounded *
                             (cellToCellRatio - 1) + 1) /
                    np.log10(cellToCellRatio))

    print('\n cells x-axis exposed', '{:6.8f}'.format(cellsXExposed))
    print('cells x-axis gap', '{:6.8f}'.format(cellsXGap))
    print('cells x-axis grounded''{:6.8f}'.format(CellsXGrounded))
    print('cells x-axis RHS', '{:6.8f}'.format(cellsXRHS))
    print('cells x-axis LHS''{:6.8f}'.format(cellsXLHS))

    print('\n cells y-axis exposed', '{:6.8f}'.format(cellsYExposed))
    print('cells y-axis dielectric', '{:6.8f}'.format(cellsYDielectric))
    print('cells y-axis grounded''{:6.8f}'.format(CellsYGrounded))
    print('cells y-axis Top', '{:6.8f}'.format(cellsYTop))
    print('cells y-axis Bottom''{:6.8f}'.format(cellsYBottom))


# Ratios = RatioElectrodes()

RatioX = RatioX(cellToCellRatio)

RatioY = RatioY(cellToCellRatio)

# cells = cells(cellToCellRatio)
