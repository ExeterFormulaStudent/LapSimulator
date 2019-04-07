import numpy as np
import math as m

def powerCurve():
    '''
    This function takes no inputs, but reads a text file of coordinates of Engine Horsepower in one column, and Engine
    RPM in the second column. The function also additionally calculates the Engine Torque from the horsepower. The
    function returns three arrays: power array, RPM array and the torque array.
    '''
    hpRPM = open('2017engine.txt', 'r')  # text file to read
    # set up empty arrays
    powerArr = np.array([])
    RPMArr = np.array([])
    hpArr = np.array([])
    torqueArray = np.array([])

    for line in hpRPM:
        try:
            RPM, hp = [float(a) for a in line.strip().split()]
            RPMArr = np.append(RPMArr, RPM)
            hpArr = np.append(hpArr, hp)
            powerArr = np.append(powerArr, hp)
        except ValueError:
            print('engine values are not numerical values.')
    #  for loop to calculate Torque at each RPM
    for i in range(0, len(RPMArr)):
        torqueNewtonMeters = (hpArr[i]*5252)/(RPMArr[i])
        torqueArray = np.append(torqueArray, torqueNewtonMeters)

    return powerArr, RPMArr, torqueArray


def enginePerformance(PowerArray, RPMArray, TorqueArray, finalDriveRatio, wheelDiameter):
    '''
    the function gearGraph takes a power array, RPM array and torque array as inputs. the function also opens and reads
    a text file consisting of the gear ratios, the RPM of the gear change from the given gear. A function of
    gearCalculation is called for each gear. The results are then plotted, with the output and array with 5 columns:
    Power, RPM, Velocity, Torque, Gear Ratio.
    '''
    gearRatios = open('gearratios.txt', 'r')  # file to read
    # set up empty arrays
    gearRatioArr = np.array([])
    gearChangeArr = np.array([])
    RPMGearVelArr = np.array([])
    # startLen is set to 0
    startLen = 0

    for line in gearRatios:  # strip and splitting the values in the file to read
        try:  # error catch of text file
            gearRatio, gearChange = [float(a) for a in line.strip().split()]
            gearRatioArr = np.append(gearRatioArr, gearRatio)
            gearChangeArr = np.append(gearChangeArr, gearChange)
        except ValueError:
            print('Gear Values given are not numerical values.')
    for n in range(0, len(gearRatioArr)):  # for loop for each gear, to output the 5 column array
        a = gearCalculation(PowerArray, RPMArray, TorqueArray, finalDriveRatio, wheelDiameter, gearRatioArr[n], gearChangeArr[n], startLen)
        RPMGearVelArr = np.append(RPMGearVelArr, a[0])
        startLen = a[1]  # set new startLen
    RPMGearVelArr = np.reshape(RPMGearVelArr, (int(len(RPMGearVelArr)/5), 5))  # reshape array to have 5 columns

    return RPMGearVelArr


def gearCalculation(PowerArray, RPMArray, TorqueArray, finalDriveRatio, wheelDiameter, gearRatio, gearChange, startLen):
    '''
    gearCalculation takes inputs of a power array, rpm array, torque array, a final drive ratio, the wheel diameter,
    gear ratio array, gear change array and a start length (therefore when a gear change occurs the rpm is not reset to
    the lowest value). The function outputs an array with 5 columns: Power, RPM, Velocity, Torque, Gear Ratio and a
    startLen value for the next iteration to start on the same rpm.
    '''

    gearArr = np.array([])
    for a in range(startLen, len(RPMArray)):  # for loop from the starLen (governed by previous gear change) and the end
        #  of the RPM array
        gearVel = ((((RPMArray[a]/60)/gearRatio)/finalDriveRatio)/(m.pi*wheelDiameter))  # calculation for the velocity
        # at the particular gear ratio
        ArrOne = np.array([PowerArray[a], RPMArray[a], gearVel, TorqueArray[a], gearRatio])  # set up array with initial values
        gearArr = np.append(gearArr, ArrOne)
        if RPMArray[a] >= gearChange:  # if the current RPM is greater than or equal to the gear change value, return
            # formatted array
            fifthLengthLast = startLen  # previous startLen is saved
            fifthLength = int(len(gearArr)/5)
            gearArr = np.reshape(gearArr, (fifthLength, 5))  # reshape the array into 5 columns
            startLen = fifthLength+fifthLengthLast  # new startLen is equal to the length of the array looped plus
            # the previous startLen value
            return gearArr, startLen  # output the array and the value of startLen
    fifthLength = int(len(gearArr)/5)
    gearArr = np.reshape(gearArr, (fifthLength, 5))  # reshape the array into 5 columns
    startLen = fifthLength+startLen

    return gearArr, startLen  # output the array and the value of startLen

if __name__ == '__main__':
    ab = enginePerformance(PowerCurve()[0],PowerCurve()[1],PowerCurve()[2],3.45,0.513)

