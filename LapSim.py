import numpy as np
import matplotlib.pyplot as plt
import Engine
import Track

def vehicleValues():
    '''
    This function reads a text file with various vehicle parameters. It returns the following values: tyre friction
    coefficient, mass, air density, frontal area, drag coefficient, efficiency, acceleration due to gravity, rolling
    resistance coefficient, area of lift, final drive ratio, wheel diameter.
    '''
    values = open('vehiclevalues.txt', 'r')  # open file to read
    for line in values:  # reading and assigning values
        tyreFrictionCoeff, mass, density, frontalArea, dragCoefficient, efficiency, g, rollingResistanceCoefficient, \
        liftArea, liftCoefficient, finalDriveRatio, wheelDiameter = [float(a) for a in line.strip().split()]
    return tyreFrictionCoeff, mass, density, frontalArea, dragCoefficient, efficiency, g, rollingResistanceCoefficient,\
           liftArea, liftCoefficient, finalDriveRatio, wheelDiameter


def maxCornerVel(tyreFrictionCoeff, normalForce, mass, effectiveRadius, airDen, Area, cd):
    '''
    The maxCornerVel function takes inputs of the tyre friction coefficient, the normal Force acting on the vehicle,
    the mass of the vehicle, the effective radius, the air density, the frontal area and the drag coefficient. The
    function returns a value of velocity which represents the maximum the vehicle can navigate the corner.
    '''
    vel = (((tyreFrictionCoeff*normalForce)**2)/(mass/(float(effectiveRadius)**2)+((0.5*airDen*Area*cd))**2))**0.25
    return vel


def acceleration(Torque, gearRatio, finalDriveRatio, efficiency, wheelDiameter, density, dragCoefficient, frontalArea,
                 vel, mass, g, coeffRollingResistance):
    '''
    this function takes a torque values, a gear ratio, final drive ratio, the efficiency of the vehicle powertrain,
    wheel diameter, air density, drag coefficient, frontal area, velocity, mass, acceleration due to gravity and the
    coefficient of rolling resistance. The function returns an a value of acceleration.
    '''
    accelForce = (Torque*gearRatio*finalDriveRatio*efficiency)/(wheelDiameter/2)
    aeroForce = (density*dragCoefficient*frontalArea*(vel**2))/2
    rollingResistance = mass*g*coeffRollingResistance

    a = (accelForce-aeroForce-rollingResistance)/mass
    return a


def normalForce(mass, g, density, cl, liftArea, vel):
    '''
    normalForce function takes inputs of mass, acceleration due to gravity, the lift coefficient(should be negative due
    to downforce, the area of lift, and the velocity. The function outputs a force value.
    '''
    F = float(mass*g - (density*cl*liftArea*(vel**2))/2)
    return F


def calculator(xCoords, yCoords, RPMGearVelArr, initialVel):
    '''
    This function takes the array of x coordinates, array of y coordinates and an array with five columns of:
    Power, RPM, Velocity, Torque, Gear Ratio. The initial normal force is calculated using the maximum velocity of the
    engine performance. The for loop iterates over the trackcoordinates, calculating the acceleration from the previous
    point, with a limit on the speed due to the maximum corner velocity. Assuming the acceleration is constant over the
    distance between the points the velocity can be calculated. The output is the time of the lap and a velocity
    distance graph.
    '''

    effectiveRadiusArr = Track.effectiveCurve(xCoords, yCoords)  # obtain track coordinates from Track script
    u = initialVel
    vArr = np.array([u])  # create velocity array with first value
    torqueInterpArr = np.array([0])  # create Torque array with 0 as the first value
    entranceVelArr = np.array([])  # create empty array for entrance velocity
    distArr = np.array([0])
    count = 0  # set count equal to 0
    for i in range(0, len(xCoords)-1):

        count = count+1  # increase count
        hyp = ((xCoords[i+1]-xCoords[i])**2 + (yCoords[i+1] - yCoords[i])**2)**0.5  # calculate hypotenuse of points
        dist = distArr[-1]+hyp  # cumulative addition of hyp to the distance
        distArr = np.append(distArr, dist)

        u = vArr[[-1]]  # u is equal to the final value in vArr
        torqueInterp = np.interp(u, RPMGearVelArr[:, 2], RPMGearVelArr[:, 3])  # interpolate the torque for the
        # given velocity
        torqueInterpArr = np.append(torqueInterpArr, torqueInterp)
        gearRatioInterp = np.interp(u, RPMGearVelArr[:, 2], RPMGearVelArr[:, 4])  # approximation of interpolating the
        # gear values
        # using function acceleration to calculate the acceleration
        a = acceleration(torqueInterp, gearRatioInterp, finalDriveRatio, efficiency, wheelDiameter, airDensity,
                         dragCoefficient, frontalArea, u, mass, g, rollingResistanceCoefficient)
        v = (u**2 + 2*hyp*a)**0.5  # assuming constant acceleration over the distance between the points
        if effectiveRadiusArr[i] != 'straight':  # if the effective radius is not 'straight', then the max corner
            # velocity may limit the velocity of the vehicle. If the velocity is larger than the maximum value it is
            # assumed the vehicle will travel at the maximum possible velocity
            float(effectiveRadiusArr[i])
            corner = maxCornerVel(tyreFrictionCoefficient, normalForce(mass, g, airDensity, liftCoefficient,
                                liftArea, u), mass, effectiveRadiusArr[i], airDensity, frontalArea, dragCoefficient)
            if v > corner:
                v = corner
        vArr = np.append(vArr, v)
        entranceVelArr = np.append(entranceVelArr, v)
    return distArr, vArr


def accelerationRestCalculator(xCoordsReversed, yCoordsReversed, RPMGearVelArr):
    '''
    accelerationRestCalculator takes the array of x coordinates, array of y coordinates and an array with five columns
    of: Power, RPM, Velocity, Torque, Gear Ratio. The initial normal force is calculated using the current velocity of
    zero. The for loop iterates over the trackcoordinates (which have been reversed back to 'forward', calculating
    the acceleration from the previous point, with a limit on the speed due to the maximum corner velocity. Assuming
    the acceleration is constant over the distance between the points the velocity can be calculated. The output is the
    time of the lap and a velocity distance graph.
    '''
    xCoords = xCoordsReversed[::-1]  # reverse the input arrays so that they are 'forward'
    yCoords = yCoordsReversed[::-1]
    effectiveRadiusArr = Track.effectiveCurve(xCoords, yCoords)
    u = 0  # initial velocity is set to zero
    vArr = np.array([u])  # create velocity array with an initial value (of 0)
    # create empty arrays
    torqueInterpArr = np.array([])
    entranceVelArr = np.array([])
    distArr = np.array([0])
    timeArr = np.array([0])

    for i in range(0, len(xCoords)-1):
        hyp = ((xCoords[i+1]-xCoords[i])**2 + (yCoords[i+1] - yCoords[i])**2)**0.5  # calculate hypotenuse
        dist = distArr[-1]+hyp  # cumulative addition of hyp to the distance
        distArr = np.append(distArr, dist)
        u = vArr[-1]  # u is equal to the final value in vArr
        torqueInterp = np.interp(u, RPMGearVelArr[:, 2], RPMGearVelArr[:, 3])  # interpolate the torque for the
        # given velocity
        gearRatioInterp = np.interp(u, RPMGearVelArr[:, 2], RPMGearVelArr[:, 4])  # approximation of interpolating the
        # gear values
        # using function acceleration to calculate the acceleration
        a = acceleration(torqueInterp, gearRatioInterp, finalDriveRatio, efficiency, wheelDiameter, airDensity,
                         dragCoefficient, frontalArea, u, mass, g, rollingResistanceCoefficient)
        v = (u**2 + 2*hyp*a)**0.5  # assuming constant acceleration over the distance between the points
        if effectiveRadiusArr[i] != 'straight':  # if the effective radius is not 'straight', then the max corner
            # velocity may limit the velocity of the vehicle. If the velocity is larger than the maximum value it is
            # assumed the vehicle will travel at the maximum possible velocity
            float(effectiveRadiusArr[i])
            corner = maxCornerVel(tyreFrictionCoefficient, normalForce(mass, g, airDensity, liftCoefficient,
                                liftArea, u), mass, effectiveRadiusArr[i], airDensity, frontalArea, dragCoefficient)
            if v > corner:
                v = corner

        vArr = np.append(vArr, v)
        torqueInterpArr = np.append(torqueInterpArr, torqueInterp)
        entranceVelArr = np.append(entranceVelArr, v)


    averageSpd = np.mean(vArr)
    return distArr, vArr


def join(reverseVel, accelDist, accelVel):
    '''
    The aim of this function is to produce a single array when the acceleration curve meets the reverse calculated
    array. The inputs of the function are: a velocity array of the reversed calculator, a distance array from the
    accelerating curve and a velocity array from the accelerating curve. The function then merges the two curves when
    the difference in the velocity is less than 1e-4. The function outputs a single distance array and a single
    velocity array.
    '''
    reverseVel = reverseVel[::-1]  # make the reverse velocity array forwards for comparison
    velocityArr = np.array([])
    for i in range(0, len(reverseVel)):
        bool = 0  # create boolean and set to 0
        if abs(reverseVel[i] - accelVel[i]) < 1e-4:  # determine whether the difference in arrays is within the
            # tolerance to join
            velocityArr = np.append(velocityArr, reverseVel[i])
            bool = 1
        elif bool == 1:
            velocityArr = np.append(velocityArr, reverseVel[i])
        else:
            velocityArr = np.append(velocityArr, accelVel[i])
    meanSpd = np.mean(abs(velocityArr))
    time = accelDist[-1]/meanSpd
    AA = len(accelDist)
    BB = len(velocityArr)
    # plt.plot(accelDist, velocityArr, color='blue')
    # plt.title('Distance Velocity graph')
    # plt.xlabel('distance/m')
    # plt.ylabel('velocity/ms^-1')
    # plt.show()
    return accelDist, velocityArr


def lapfunction(dist, vel):

    boolean = True  # boolean for error catch

    while boolean:  # while loop for error catch
        try:
            lapNumber = int(input('Number of laps'))  # input number
            if lapNumber <= 0:  # used to check if lapNumber is positive
                print('Number of laps must be a positive integer.')
            else:  # exits while loop if satisfied
                boolean = False
        except ValueError:  # ValueError raised
            print('Number of laps must be a positive integer.')

    distArr = np.array(dist)
    velArr = np.array(vel)
    if lapNumber == 1:  # if only a single lap
        plt.plot(distArr, velArr)
        plt.show()
        return distArr, velArr
    else:  # if lap is greater than one append more laps to the existing array - using the latest value of
        # velocity for the initial velocity of each lap
        lapOneTime = distArr[-1]/np.mean(velArr)
        print('lap time for lap 1          : ', '{:.4f}'.format(lapOneTime), 'seconds')
        for n in range(0, lapNumber-1):
            calculatedArray = calculator(Track.coordinateImport()[0], Track.coordinateImport()[1],
                                         Engine.enginePerformance(Engine.PowerCurve()[0], Engine.PowerCurve()[1],
                                            Engine.PowerCurve()[2], finalDriveRatio, wheelDiameter), velArr[-1])
            calculatedArray = calculatedArray[::-1]
            lapDist = calculatedArray[1][-1]
            lapTime = lapDist/np.mean(calculatedArray[0])
            dist = distArr[-1] + calculatedArray[1]

            distArr = np.append(distArr, dist)
            velArr = np.append(velArr, calculatedArray[0])
            if 9 < n+2 < 100:  # if loop for formatting depending on number of figures in the lap number
                print('lap time for lap', n+2, '        : ', '{:.6f}'.format(lapTime), 'seconds')
            elif n+2 < 10:
                print('lap time for lap', n+2, '         : ', '{:.6f}'.format(lapTime), 'seconds')
            else:
                print('lap time for lap', n+2, '       : ', '{:.6f}'.format(lapTime), 'seconds')
    aveSpd = np.mean(velArr)  # find average speed of vehicle
    totaldist = distArr[-1]  # obtain the total distance covered
    time = totaldist/aveSpd  # calculate time
    print('Average Velocity of all laps: ', '{:.6f}'.format(aveSpd), 'metres per second',
          '({:.6f}'.format(aveSpd*2.23694), 'mph)')
    print('Total Distance of all laps  : ', '{:.6f}'.format(totaldist), 'metres')
    return distArr, velArr


if __name__ == '__main__':

    tyreFrictionCoefficient = float(vehicleValues()[0])
    mass = float(vehicleValues()[1])
    airDensity = float(vehicleValues()[2])
    frontalArea = float(vehicleValues()[3])
    dragCoefficient = float(vehicleValues()[4])
    efficiency = float(vehicleValues()[5])
    g = float(vehicleValues()[6])  # ms**-2
    rollingResistanceCoefficient = float(vehicleValues()[7])
    liftArea = float(vehicleValues()[8])
    liftCoefficient = float(vehicleValues()[9])
    finalDriveRatio = float(vehicleValues()[10])
    wheelDiameter = float(vehicleValues()[11])

    engineResults = Engine.enginePerformance(Engine.PowerCurve()[0], Engine.PowerCurve()[1], Engine.PowerCurve()[2],
                                             finalDriveRatio, wheelDiameter)
    fig, ax1 = plt.subplots(figsize=(10, 5))  # plot engine results (Power, RPM and Torque against velocity)

    ax1.set_xlabel('Velocity [m/s]')
    ax1.set_ylabel('Power [bhp]')
    plt.plot(engineResults[:, 2], engineResults[:, 0], color='tab:red',label='Power')
    plt.legend(loc='lower right')
    ax2 = ax1.twinx()
    ax2.set_ylabel('RPM [RPM]', color='tab:blue')
    plt.plot(engineResults[:, 2], engineResults[:, 1], color='tab:blue', label='RPM')
    ax2.tick_params(axis='y', labelcolor='tab:blue')
    plt.title('Engine Power and Engine RPM against Vehicle Velocity')
    plt.legend(loc='upper left')
    plt.show()

    plt.plot(engineResults[:, 2], engineResults[:, 3], color='tab:green')
    plt.title('Engine Torque against Vehicle Velocity')
    plt.ylabel('Torque [Nm]')
    plt.xlabel('Vehicle Velocity [m/s]')
    plt.show()

    coordinates = Track.coordinateImport()  # plot track coordinates
    plt.plot(coordinates[0], coordinates[1])
    plt.xlabel('x')
    plt.axis('equal')
    plt.ylabel('y')
    plt.title('Track')
    plt.show()

    a = (calculator(Track.coordinateImport()[0], Track.coordinateImport()[1], Engine.enginePerformance
    (Engine.PowerCurve()[0], Engine.PowerCurve()[1], Engine.PowerCurve()[2], finalDriveRatio, wheelDiameter), 25))

    b = (accelerationRestCalculator(Track.coordinateImport()[0], Track.coordinateImport()[1], Engine.enginePerformance
    (Engine.PowerCurve()[0], Engine.PowerCurve()[1], Engine.PowerCurve()[2], finalDriveRatio, wheelDiameter)))
    # join(a[1], b[0], b[1])
    totaldistandvel = lapfunction(join(a[1], b[0], b[1])[0], join(a[1], b[0], b[1])[1])
    aveSpd = np.mean(totaldistandvel[1])
    time = totaldistandvel[0][-1]/aveSpd
    print('Total Race Time             : ', '{:.6f}'.format(time), 'seconds')

    plt.plot(totaldistandvel[0], totaldistandvel[1])
    plt.title('Lap Simulator Total Distance against Velocity')
    plt.xlabel('Distance [m]')
    plt.ylabel('Velocity [m/s]')
    plt.show()
