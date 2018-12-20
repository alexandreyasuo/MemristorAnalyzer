import exitroutines
# exit and cleanup

import circlefit
# circle fittings

import csv
# needed to read raw data in CSV

import os
import glob
import sys
# needed for directory manipulation

import configparser
# for config file

import numpy as np
# numerical functions

import io
# needed for separating bitstreams

from openpyxl import Workbook
# needed to save to XLSX

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# plotting

# Initialization of parameters
config = configparser.ConfigParser()

#counter variables
warning = 0
verbose = 10

if os.path.exists('config.ini'):
    config.read('config.ini')
else:
    print('WARNING: Configuration file not found. Initiating config file with default values')
    warning += 1
    configfile = open('config.ini', 'w')
    config['Directories'] = {'WorkingDir' : os.path.dirname(os.path.realpath(__file__))}
    config['Settings'] = {'Verbose' : 10}
    config.write(configfile)
    configfile.flush()
    configfile.close()
    print('Done')

if 'WorkingDir' in config['Directories']:
    workdirectory = config['Directories']['WorkingDir']
else:
    print('ERROR: Working directory configuration not in configuration. Is the ini file correct and not corrupted?')
    exitroutines.getout(warning, -1)

try:
    verbose = (int)(config['Settings']['Verbose'])
except:
    print('WARNING: Verbose not in configuration file. Assigning default value: ', verbose)
    warning += 1

print('Changing to working directory...')
try:
    os.chdir(os.path.join(workdirectory, 'data'))
except:
    print('ERROR: Could not change to working directory. Does directory exist?')
    exitroutines.getout(warning, -1)

numberofcsvs = [i for i in glob.glob('*.{}'.format('csv'))]

if numberofcsvs == []:
    print('Nothing to do, no CSV files found')
    exitroutines.getout(warning, 1)

frequencies = []
forwardareas = []
backwardareas = []
voltages = []
currents = []
memvolts = []
memamps = []
propermem = []

figurenumber = 1

for i in numberofcsvs:
    try:
        currentfile = open(i,'r')
        text = currentfile.read()
        print('Currently working with file: ', currentfile.name)
        measurements = text.split('Step')[1:]
        for mezzies in measurements:
            tempfile = open('Temp.txt', 'w')
            tempfile.write(mezzies)
            tempfile.close()
            # trouble using buffers. Just write each simulation to a temp file and scan that
            # TODO: horribly not elegant. Can I go back to using buffers and make it work?

            data = np.genfromtxt('Temp.txt', skip_header=1, skip_footer=0, names=['Time', 'Voltage', 'Current'])
            # separate V and I. Skip the "Step" header. No footers, skip time column.

            dataNum = len(data['Voltage'])
            # read amount of simulation points.

            forwardLoop = 0.0
            backwardLoop = 0.0
            # total hysteresis for quadrant 1 and 3 respectively.
                
            partialInt1 = 0.0
            partialInt2 = 0.0
            # partial integral for "going" swing and "returning" swing.
            
            negative = 0
            # has the data turned negative?
                
            loopway = 0
            # is it "going" or "returning"?
            # initial conditions.

            currentflux = 0
            currentcharge = 0
            # hopefully at beginning of cycle

            currentchargeintegral = 0

            charge = []
            flux = []

            chargeintegral = []

            resistances = []
            ivplane = []

            angles = []

            voltagephasor_angle = []

            currentphasor_angle = []

            charge.append(currentcharge)
            flux.append(currentflux)

            chargeintegral.append(currentchargeintegral)

            initialtime = data['Time'][0]
            finaltime = data['Time'][len(data['Time'])-1]

            if verbose > 6:
                print('Frequency detected: ', "{:.2e}".format(1/(finaltime-initialtime)), '[Hz]')

            maxV = data['Voltage'][0]
            # to detect turning from "going" to "returning"

            resistances.append(data['Voltage'][0]/data['Current'][0])

            voltagephasor_angle.append(np.arccos(data['Voltage'][0]/np.max(data['Voltage'])))
            currentphasor_angle.append(np.arccos(data['Current'][0]/np.max(data['Current'])))

   
            if(data['Voltage'][0] > 0):
                ivplane.append(np.sqrt(data['Voltage'][0]**2 + data['Current'][0]**2))
            else:
                ivplane.append(-np.sqrt(data['Voltage'][0]**2 + data['Current'][0]**2))

            if(data['Voltage'][0] == 0):
                angles.append(np.pi/2)
            else:
                angles.append(np.arctan(data['Voltage'][0]/data['Current'][0]))
                
            for j in range(dataNum - 1):

                if (verbose > 9):

                    plt.figure()

                    plt.ylim(-1,1)
                    plt.xlim(-1,1)

                    origin = [0], [0]

                    voltagephasor = np.sin(voltagephasor_angle[j]), np.sin(currentphasor_angle[j])
                    currentphasor = np.cos(voltagephasor_angle[j]), np.cos(currentphasor_angle[j])

                    plt.quiver(*origin, voltagephasor, currentphasor, color=['r','b'], scale=5)
                    
                    plt.savefig('Phasor' + str(j) +'.png')
                    plt.close()

                currentchargeintegral = currentchargeintegral + currentcharge*(data['Time'][j+1] - data['Time'][j])
                chargeintegral.append(currentchargeintegral)

                currentcharge = currentcharge + (data['Time'][j+1] - data['Time'][j])*(data['Current'][j+1] - data['Current'][j])/2 + (data['Time'][j+1] - data['Time'][j])*data['Current'][j]
                currentflux = currentflux + (data['Time'][j+1] - data['Time'][j])*(data['Voltage'][j+1] - data['Voltage'][j])/2 + (data['Time'][j+1] - data['Time'][j])*data['Voltage'][j]

                resistances.append(data['Voltage'][j]/data['Current'][j])

                voltagephasor_angle.append(np.arccos(np.abs(data['Voltage'][j]/np.max(data['Voltage']))))
                currentphasor_angle.append(np.arccos(np.abs(data['Current'][j]/np.max(data['Current']))))

                if(data['Voltage'][j] > 0):
                    ivplane.append(np.sqrt(data['Voltage'][j]**2 + data['Current'][j]**2))
                else:
                    ivplane.append(-np.sqrt(data['Voltage'][j]**2 + data['Current'][j]**2))

                if(data['Voltage'][j] == 0):
                    angles.append(np.pi/2)
                else:
                    angles.append(np.arctan(data['Voltage'][j]/data['Current'][j]))
                
                charge.append(currentcharge)
                flux.append(currentflux)
                
                if negative == 0:
                    if data['Voltage'][j] >= maxV:
                        maxV = data['Voltage'][j]
                                        # update max V while "going"
                    else:
                        loopway = 1
                                        # not "going", so obviously "returning"
                    if loopway == 0:
                        partialInt1 = partialInt1 + (data['Voltage'][j+1] - data['Voltage'][j])*(data['Current'][j+1] - data['Current'][j])/2 + (data['Voltage'][j+1] - data['Voltage'][j])*data['Current'][j]
                                        # integrate "going"
                    else:
                        partialInt2 = partialInt2 + (data['Voltage'][j+1] - data['Voltage'][j])*(data['Current'][j+1] - data['Current'][j])/2 + (data['Voltage'][j+1] - data['Voltage'][j])*data['Current'][j]
                                        # integrate "returning"
                    if data['Voltage'][j] < 0 and j != 0:
                        negative = 1
                        partialInt2 = -partialInt2
                        forwardLoop = partialInt1 - partialInt2
                        partialInt1 = 0.0
                        partialInt2 = 0.0
                        maxV = 0
                        loopway = 0
                                        # first quadrant done, calculate and go to third quadrant now
                else:
                    if np.abs(data['Voltage'][j]) >= maxV:
                        maxV = np.abs(data['Voltage'][j])
                                        # same, while "going" update maxV
                    else:
                        loopway = 1
                                        # else "returning"
                    if loopway == 0:
                        partialInt1 = partialInt1 + (data['Voltage'][j+1] - data['Voltage'][j])*data['Current'][j]
                                        # integrate "going"
                    else:
                        partialInt2 = partialInt2 + (data['Voltage'][j+1] - data['Voltage'][j])*data['Current'][j]
                                        # integrate "returning"

            partialInt2 = -partialInt2
            backwardLoop = partialInt1 - partialInt2
            if forwardLoop < 0:
                forwardDirection = 'CCW'
            else:
                forwardDirection = 'CW'
            if backwardLoop < 0:
                backwardDirection = 'CCW'
            else:
                backwardDirection = 'CW'
        

            frequencies.append(1/(finaltime-initialtime))
            forwardareas.append(forwardLoop)
            backwardareas.append(backwardLoop)

            voltagecircle = circlefit.leastsquares(data['Voltage'],resistances)
            currentcircle = circlefit.leastsquares(data['Current'],resistances)

            anglescircle = circlefit.leastsquares(ivplane, angles)

            ivplanecircle = circlefit.leastsquares(ivplane,resistances)

            fittedradiusvoltage = voltagecircle[0]*voltagecircle[6]
            fittedradiuscurrent = currentcircle[0]*currentcircle[6]
            fittedradiusiv = ivplanecircle[0]*ivplanecircle[6]
            fittedradiusangles = anglescircle[0]*anglescircle[6]

            R1 = voltagecircle[7]
            R2 = currentcircle[7]
            R3 = ivplanecircle[7]
            R4 = anglescircle[7]

            excitationvoltage = voltagecircle[0]*voltagecircle[5]
            excitationcurrent = currentcircle[0]*currentcircle[5]
            excitationiv = ivplanecircle[0]*ivplanecircle[5]

            cirfitvoltage_x = []
            cirfitvoltage_y = []
            cirfitcurrent_x = []
            cirfitcurrent_y = []

            cirfitiv_x = []
            cirfitiv_y = []

            anglesiv_x = []
            anglesiv_y = []

            voltages.append(excitationvoltage)
            currents.append(excitationcurrent)

            memvolts.append(fittedradiusvoltage/(excitationvoltage*(finaltime-initialtime)))
            memamps.append(fittedradiuscurrent/(excitationcurrent*(finaltime-initialtime)))

            propermem.append(fittedradiusiv/(excitationiv*(finaltime-initialtime)))

            for k in range(101):
                theta = k*2*np.pi/100
                cirfitvoltage_x.append((voltagecircle[1]+np.cos(theta)*voltagecircle[0])*voltagecircle[5]+voltagecircle[3])
                cirfitvoltage_y.append((voltagecircle[2]+np.sin(theta)*voltagecircle[0])*voltagecircle[6]+voltagecircle[4])
                cirfitcurrent_x.append((currentcircle[1]+np.cos(theta)*currentcircle[0])*currentcircle[5]+currentcircle[3])
                cirfitcurrent_y.append((currentcircle[2]+np.sin(theta)*currentcircle[0])*currentcircle[6]+currentcircle[4])

                cirfitiv_x.append((ivplanecircle[1]+np.cos(theta)*ivplanecircle[0])*ivplanecircle[5]+ivplanecircle[3])
                cirfitiv_y.append((ivplanecircle[2]+np.sin(theta)*ivplanecircle[0])*ivplanecircle[6]+ivplanecircle[4])

                anglesiv_x.append((anglescircle[1]+np.cos(theta)*anglescircle[0])*anglescircle[5]+anglescircle[3])
                anglesiv_y.append((anglescircle[2]+np.sin(theta)*anglescircle[0])*anglescircle[6]+anglescircle[4])

            if (verbose > 7):
                print('Foward loop: ', "{:.2e}".format(np.abs(forwardLoop)), '[W] ', forwardDirection)
                print('Backward loop: ', "{:.2e}".format(np.abs(backwardLoop)), '[W] ', backwardDirection)
                #print('Memristance (by voltage): ', "{:.2e}".format(fittedradiusvoltage/(excitationvoltage*(finaltime-initialtime))), '[Ch] --- R = ', R1)
                #print('Memristance (by current): ', "{:.2e}".format(fittedradiuscurrent/(excitationcurrent*(finaltime-initialtime))), '[Ch] --- R = ', R2)
                print('Memristance (IV method): ', "{:.2e}".format(fittedradiusiv/(excitationiv*(finaltime-initialtime))), '[Ch] --- R = ', R3)
            
            if (verbose > 8):
                
                plt.figure(figsize=(8,8))
                plt.subplot(221)
                plt.plot(data['Voltage'], data['Current'])
                title = str('I x V\nForward loop area: '+str("{:.2e}".format(np.abs(forwardLoop)))+'[W] '+forwardDirection+'\nBackward loop area: '+str("{:.2e}".format(np.abs(backwardLoop)))+ '[W] '+ backwardDirection)
                plt.title(title)
                plt.ylabel('Current [A]')
                plt.xlabel('Voltage [V]')
                plt.axhline(color = 'k', linestyle = '--', linewidth = 0.5)
                plt.axvline(color = 'k', linestyle = '--', linewidth = 0.5)
                plt.tight_layout()

                plt.subplot(222)
                plt.plot(flux, charge)
                title = str('Frequency detected: ' + str("{:.2e}".format(1/(finaltime-initialtime))) + '[Hz]')
                plt.title(title)
                plt.ylabel('Charge [C]')
                plt.xlabel('Flux [Wb]')
                plt.tight_layout()

                plt.subplot(223)
                plt.plot(ivplane, resistances)
                plt.plot(cirfitiv_x, cirfitiv_y, color = 'k', linestyle = '--', linewidth = 0.5)
                title = str('Memristance: '+ str("{:.2e}".format(fittedradiusiv/(excitationiv*(finaltime-initialtime)))) + '[Ch]\nR = ' + str(R3))
                plt.title(title)
                plt.ylabel('Resistance [Ohm]')
                plt.xlabel('IV Plane')
                plt.tight_layout()

                plt.subplot(224)
                plt.plot(ivplane, angles)
                plt.plot(anglesiv_x, anglesiv_y, color = 'k', linestyle = '--', linewidth = 0.5)
                title = str('Memristance: ' + str("{:.2e}".format(fittedradiusiv/(excitationiv*(finaltime-initialtime)))) + '[Ch]\nR = ' + str(R4))
                plt.title(title)
                plt.ylabel('Phase angle [radians]')
                plt.xlabel('IV Plane')
                plt.tight_layout()
                
                plt.savefig('Cycle' + str(figurenumber) +'PartA.png')
                plt.close()

            if (verbose > 9):
                
                plt.figure(figsize=(8,8))
                plt.subplot(221)
                plt.plot(resistances, charge)
                title = str('ohm x q')
                plt.title(title)
                plt.ylabel('Charge [C]')
                plt.xlabel('Voltage [V]')
                #plt.axhline(color = 'k', linestyle = '--', linewidth = 0.5)
                #plt.axvline(color = 'k', linestyle = '--', linewidth = 0.5)
                plt.tight_layout()

                plt.subplot(222)
                plt.plot(data['Voltage'], chargeintegral)
                title = str('intq x V')
                plt.title(title)
                plt.ylabel('Charge Integral [Cs]')
                plt.xlabel('Voltage [V]')
                plt.tight_layout()

                plt.subplot(223)
                plt.plot(flux, chargeintegral)
                title = str('intq x phi')
                plt.title(title)
                plt.ylabel('Charge Integral [Cs]')
                plt.xlabel('Flux [Wb]')
                plt.tight_layout()

                plt.subplot(224)
                plt.plot(ivplane, angles)
                plt.plot(anglesiv_x, anglesiv_y, color = 'k', linestyle = '--', linewidth = 0.5)
                title = str('Memristance: ' + str("{:.2e}".format(fittedradiusiv/(excitationiv*(finaltime-initialtime)))) + '[Ch]\nR = ' + str(R4))
                plt.title(title)
                plt.ylabel('Phase angle [radians]')
                plt.xlabel('IV Plane')
                plt.tight_layout()
                
                plt.savefig('Cycle' + str(figurenumber) +'PartB.png')
                plt.close()

                '''

                np.savetxt('Flux' + str(figurenumber) +'.csv', flux, delimiter=",")
                np.savetxt('Charge' + str(figurenumber) +'.csv', charge, delimiter=",")

                

                fig = plt.figure()
                ax = fig.gca(projection='3d')
                ax.plot(data['Voltage'], -data['Current'], angles)
                plt.show()
                #for l in range(0,90):
                    #ax.view_init(90-l,l/2)
                    #plt.savefig("movie"+str(l)+".png")                
                
                input("Press Enter to continue...")

                '''

            figurenumber += 1
                
        currentfile.close()
    except TypeError:
        print('WARNING: Could not open or process CSV file ', i)
        warning += 1

deltaarea = np.diff(np.abs(forwardareas)) / np.abs(forwardareas[:-1])
deltaf = np.diff(frequencies) / np.abs(frequencies[:-1])

memperarea = (deltaf) / deltaarea
#print(memperarea)

plt.figure(figsize=(8,8))
plt.subplot(221)
plt.semilogx(frequencies, forwardareas)
title = str('Forward loop area x f')
plt.title(title)
plt.ylabel('Forward loop area [W]')
plt.xlabel('Frequency [Hz]')
plt.tight_layout()

plt.subplot(222)
plt.semilogx(frequencies, backwardareas)
title = str('Backward loop area x f')
plt.title(title)
plt.ylabel('Backward loop area [W]')
plt.xlabel('Frequency [Hz]')
plt.tight_layout()

plt.subplot(223)
plt.plot(propermem)
title = str('Memristance (IV method)\nAverage: ' + str("{:.2e}".format(np.average(propermem))))
plt.title(title)
plt.ylabel('Memristance [Ch]')
plt.xlabel('Count')
plt.tight_layout()

plt.subplot(224)
plt.plot(memperarea)
title = str('Memristance (Area method)\nAverage: ' + str("{:.2e}".format(np.average(propermem))))
plt.title(title)
plt.ylabel('Memristance [Ch]')
plt.xlabel('Count')
plt.tight_layout()
                
plt.savefig('Final.png')
plt.close()

'''
plt.plot(memperarea)
title = str('Memristance (Area method)\nAverage: ' + str("{:.2e}".format(np.average(propermem))))
plt.title(title)
plt.ylabel('Memristance [Ch]')
plt.xlabel('Count')
plt.tight_layout()
plt.show()
input("Press Enter to continue...")
plt.close()
'''

exitroutines.getout(warning, 1)
