import numpy as np
import math as m

def coordinateImport():
    '''
    This function does not take inputs, but reads a text file with coordinates of the track, with the first column the
    x coordinates and the second, the y coordinates. The function output are two arrays of coordinates (x and y) which
    are reversed.
    '''
    trackcoordinates = open('skidpan.txt', 'r')  # open file to read
    #set empty array
    xArr = np.array([])
    yArr = np.array([])

    for line in trackcoordinates:
        xcoord, ycoord = [float(a) for a in line.strip().split()]
        xArr = np.append(xArr, xcoord)
        yArr = np.append(yArr, ycoord)

    # reverse arrays so that the track calculations are reversed
    xCoords = xArr[::-1]
    yCoords = yArr[::-1]
    return xCoords, yCoords


def effectiveCurve(xCoords, yCoords):
    '''
    The curvecalculator function takes the x and y coordinate arrays as inputs and creates initial and final values for
    the effective radius between three points. The for loop creates the effective radius for all other values in the
    track coordinates. The array of effective radii for each point on the track coordinates is output.
    '''
    # initial length values using the first two points and the last point.
    a0 = ((xCoords[1]-xCoords[-1])**2 + (yCoords[1] - yCoords[-1])**2)**0.5
    b0 = ((xCoords[0]-xCoords[-1])**2 + (yCoords[0] - yCoords[-1])**2)**0.5
    c0 = ((xCoords[1]-xCoords[0])**2 + (yCoords[1] - yCoords[0])**2)**0.5
    effectiveRadiusArr = np.array([])

    if c0 + b0 == a0:  # if the value of c0 + b0 is equal to a0 all three points are in a straight line, therefore
        # recorded as string 'straight'
        effectiveRadius0 = 'straight'
        effectiveRadiusArr = np.append(effectiveRadiusArr,effectiveRadius0)
    else:  # if they are not in a straight line, then the effective radius of the three points is estimated using trig,
        # cosine rule and circle theorem mathematics
        effectiveRadius0 = a0/(2*(np.sin(m.pi - abs(np.arccos(abs(((b0**2)+(c0**2)-(a0**2))/(2*b0*c0))))*180/m.pi)))
        effectiveRadiusArr = np.append(effectiveRadiusArr, effectiveRadius0)

    for i in range(1, len(xCoords)-1):
        # loop for points i-1,i and i+1 (outside the first and final positions).
        a = ((xCoords[i+1]-xCoords[i-1])**2 + (yCoords[i+1] - yCoords[i-1])**2)**0.5
        b = ((xCoords[i]-xCoords[i-1])**2 + (yCoords[i] - yCoords[i-1])**2)**0.5
        c = ((xCoords[i+1]-xCoords[i])**2 + (yCoords[i+1] - yCoords[i])**2)**0.5
        if b + c == a:  # if the value of c0 + b0 is equal to a0 all three points are in a straight line, therefore
        # recorded as string 'straight'
            effectiveRadius = 'straight'
            effectiveRadiusArr = np.append(effectiveRadiusArr,effectiveRadius)
        else:  # if they are not in a straight line, then the effective radius of the three points is estimated using trig,
        # cosine rule and circle theorem mathematics
            intstepresult = abs(((b**2)+(c**2)-(a**2))/(2*b*c))
            thirdstep = a/(2*np.sin(m.pi - abs(np.arccos(intstepresult)))*180/m.pi)
            effectiveRadius = a/(2*(np.sin(180 - abs(np.arccos(abs(intstepresult)*(m.pi/180))))*(m.pi/180)))
            effectiveRadiusArr = np.append(effectiveRadiusArr, effectiveRadius)
    # final length values using the first point and the two last points.
    aend = ((xCoords[0]-xCoords[-2])**2 + (yCoords[0] - yCoords[-2])**2)**0.5
    bend = ((xCoords[-1]-xCoords[-2])**2 + (yCoords[-1] - yCoords[-2])**2)**0.5
    cend = ((xCoords[0]-xCoords[-1])**2 + (yCoords[0] - yCoords[-1])**2)**0.5
    if cend + bend == aend:  # if the value of c0 + b0 is equal to a0 all three points are in a straight line, therefore
        # recorded as string 'straight'
        effectiveRadiusEnd = 'straight'
        effectiveRadiusArr = np.append(effectiveRadiusArr, effectiveRadiusEnd)
    else:  # if they are not in a straight line, then the effective radius of the three points is estimated using trig,
        # cosine rule and circle theorem mathematics
        effectiveRadiusEnd = aend/(2*(np.sin(m.pi - np.arccos(((bend**2)+(cend**2)-(aend**2))/(2*bend*cend))))*(180/m.pi))
        effectiveRadiusArr = np.append(effectiveRadiusArr, effectiveRadiusEnd)
    return effectiveRadiusArr
