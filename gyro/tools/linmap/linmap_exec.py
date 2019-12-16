#!/usr/bin/env python

# USAGE: 
#
#    python linmap_exec.py <directory> <nProc> <start_over>

import sys
import os
import shutil
import subprocess
import numpy as np

#----------------------------------------
# Capture command-line arguments
#

directory = sys.argv[1]
nProc = int(sys.argv[2])
start_over = int(sys.argv[3])
start_over = bool(start_over)


#To solve python rounding issue, i.e. floats are not represented by exact values (important for comparisons).
#MinFloat = sys.float_info.epsilon
MinFloat = pow(10, -sys.float_info.dig)

#Go to directory
os.chdir(directory)

print("Starting gyro_linmap in directory : %s" % os.getcwd())
print("**********************************************************")


#-------------------------------------------------------
# Read input parameters from file

inputParams = np.genfromtxt('input.linmap.gen', dtype='str')
VariableValue = inputParams[:,0].astype(np.float)
VariableName = inputParams[:,1]
NumVariables = len(VariableValue)

i=0
while i < NumVariables:
	if VariableName[i] == 'DET_TOLERANCE':
		DET_TOLERANCE = VariableValue[i]
	elif VariableName[i] == 'ERROR_TOLERANCE':
		ERROR_TOLERANCE = VariableValue[i]
	elif VariableName[i] == 'WI_STABLE_LIMIT':
		WI_STABLE_LIMIT = VariableValue[i]
	elif VariableName[i] == 'startWR':
		startWR = VariableValue[i]
	elif VariableName[i] == 'startWI':
		startWI = VariableValue[i]
	elif VariableName[i] == 'RADIUSstart':
		RADIUSstart = VariableValue[i]
	elif VariableName[i] == 'RADIUSmin':
		RADIUSmin = VariableValue[i]
	elif VariableName[i] == 'RADIUSmax':
		RADIUSmax = VariableValue[i]
	elif VariableName[i] == 'RADIUSstep':
		RADIUSstep = VariableValue[i]
	elif VariableName[i] == 'L_Ystart':
		L_Ystart = VariableValue[i]
	elif VariableName[i] == 'L_Ymin':
		L_Ymin = VariableValue[i]
	elif VariableName[i] == 'L_Ymax':
		L_Ymax = VariableValue[i]
	elif VariableName[i] == 'L_Ystep':
		L_Ystep = VariableValue[i]
	else:
		print("Got false variable") 
	i += 1

print("Parameters:")
print("DET_TOLERANCE: %E" % DET_TOLERANCE)
print("ERROR_TOLERANCE: %E" % ERROR_TOLERANCE)
print("WI_STABLE_LIMIT: %E" % WI_STABLE_LIMIT)
print("startWR: %E" % startWR)
print("startWI: %E" % startWI)
print("RADIUSstart: %E" % RADIUSstart)
print("RADIUSmin: %E" % RADIUSmin)
print("RADIUSmax: %E" % RADIUSmax)
print("RADIUSstep: %E" % RADIUSstep)
print("L_Ystart: %E" % L_Ystart)
print("L_Ymin: %E" % L_Ymin)
print("L_Ymax: %E" % L_Ymax)
print("L_Ystep: %E" % L_Ystep)
print("**********************************************************")
#-------------------------------------------------------

#----------------------------------------
NstepsR = int((RADIUSmax - RADIUSmin) / RADIUSstep)
NstepsL_Y= int((L_Ymax - L_Ymin) / L_Ystep)

print("Steps in RADIUS = %i" % NstepsR)
print("Steps in L_Y = %i" % NstepsL_Y)
#----------------------------------------
#Set up first simulation
os.chdir(directory)
if not os.path.exists(directory + '/input.gyro.orig'):
	print("No input.gyro.orig found in %s" % directory)
	sys.exit(0)

shutil.copyfile('input.gyro.orig', 'input.gyro')

input_gyro = open('input.gyro', 'a')
input_gyro.write("\n")
input_gyro.write("FIELDEIGEN_WR=%E\n" % startWR)
input_gyro.write("FIELDEIGEN_WI=%E\n" % startWI)
input_gyro.close()
#----------------------------------------

if start_over :
	print("Starting from scratch!")
	if os.path.exists(directory + '/radius_krho_fieldeigen_omega.out'):
		os.remove("radius_krho_fieldeigen_omega.out")
	if os.path.exists(directory + '/fieldeigen_gbflux.out'):
		os.remove("fieldeigen_gbflux.out")
	if os.path.exists(directory + '/rho.out'):
		os.remove("rho.out")
else :
	print("Restarting!")

#----------------------------------------
print("Your tolerance is |det| < %E	error < %E" % (DET_TOLERANCE,  ERROR_TOLERANCE))
print("An eigenmode is considered stable if gamma < %F" % WI_STABLE_LIMIT)
#----------------------------------------

restartRadiusWR=startWR
restartRadiusWI=startWI

resultsWR=1
resultsWI=1

currentRADIUS=RADIUSstart
#Start by scanning downwards in RADIUS
downwardsRadius=True
radiusbreak=False

while not (currentRADIUS < RADIUSmin - MinFloat or currentRADIUS > RADIUSmax + MinFloat):
	radiusbreak=False
	currentL_Y=L_Ystart
	#start by scanning downwards in L_Y
	downwardsL_Y=True
	
	while not (currentL_Y < (L_Ymin - MinFloat) or currentL_Y > (L_Ymax + MinFloat)):
		converged=False
		trials=0
		#Loop to find a solution
		while not converged :
			print("**********************************************************")
			print("RADIUS = %F" % currentRADIUS)
			print("L_Y = %F" % currentL_Y)
			input_gyro = open('input.gyro', 'a')
			input_gyro.write("\n")
			input_gyro.write("RADIUS=%F\n" % currentRADIUS)
			input_gyro.write("L_Y=%F\n" % currentL_Y)
			input_gyro.close()
			#Run GYRO
			subprocess.call(["gyro", "-e", ".", "-n", "%i" %nProc])
			fieldeigen_out = open('fieldeigen.out', 'r')
			lines = fieldeigen_out.readlines()
			results = list(map(float, lines[-1].split())) #Last line
			fieldeigen_out.close()
			resultsWR=results[0]
			resultsWI=results[1]
			det=results[2]
			error=results[3]

			if det <= (DET_TOLERANCE + MinFloat) and error <= (ERROR_TOLERANCE + MinFloat):
				print("Converged!")
				converged=True
				break

			trials += 1
			if trials >= 3 : #Give up searching
				if currentL_Y == L_Ystart :
					radiusbreak=True
				break
			
			if currentL_Y == L_Ystart :

				if trials == 1 : #Try reducing step in RADIUS to half
					print("Try reducing step length in RADIUS to half")
					if downwardsRadius :
						currentRADIUS = currentRADIUS + 0.5000*RADIUSstep
					else :
						currentRADIUS = currentRADIUS - 0.5000*RADIUSstep
				if trials == 2 : #Try reducing step in RADIUS to 1/4th
					print("Try reducing step length in RADIUS to 1/4th")
					if downwardsRadius :
						currentRADIUS = currentRADIUS + 0.2500*RADIUSstep
					else :
						currentRADIUS = currentRADIUS - 0.2500*RADIUSstep

			else :

				if trials == 1 : #Try reducing step in L_Y to half
					print("Try reducing step length in L_Y to half")
					if downwardsL_Y :
						currentL_Y = currentL_Y + 0.5000*L_Ystep
					else :
						currentL_Y = currentL_Y - 0.5000*L_Ystep
				if trials == 2 : #Try reducing step in L_Y to 1/4th
					print("Try reducing step length in L_Y to 1/4th")
					if downwardsL_Y :
						currentL_Y = currentL_Y + 0.2500*L_Ystep
					else :
						currentL_Y = currentL_Y - 0.2500*L_Ystep

		if converged : #Write eigenvalue to output files if converged
			if resultsWI >= (WI_STABLE_LIMIT - MinFloat) :
				print("Writing to output")
				radius_krho_fieldeigen_omega = open('radius_krho_fieldeigen_omega.out', 'a')
				radius_krho_fieldeigen_omega.write("%F\t%F\t%E\t%E\t%E\t%E\n" % (currentRADIUS, currentL_Y, results[0], results[1], results[2], results[3]))
				radius_krho_fieldeigen_omega.close()

				out_gyro_gbflux = open('out.gyro.gbflux', 'r')
				lines = out_gyro_gbflux.readlines()
				fieldeigen_gbflux = open('fieldeigen_gbflux.out', 'a')
				fieldeigen_gbflux.write("%s\n" % lines[-1]) #Last line
				fieldeigen_gbflux.close()
				out_gyro_gbflux.close()
				
				#####NOT YET IMPLEMENTED, CHECK IF NEEDED?
				#profiles_gen -i input.profiles -loc_rad "$currentRADIUS" | grep "rho " >> rho.out
				#rho
		else :
			print("Couldn't find solution! Stop following branch")
			if radiusbreak :
				break

			if downwardsL_Y : #Start scanning upwards
				downwardsL_Y=False
				currentL_Y = L_Ystart + L_Ystep
				shutil.copyfile('input.gyro.orig', 'input.gyro')
				input_gyro = open('input.gyro', 'a')
				input_gyro.write("\n")
				input_gyro.write("FIELDEIGEN_WR=%E\n" % restartRadiusWR)
				input_gyro.write("FIELDEIGEN_WI=%E\n" % restartRadiusWI)
				input_gyro.close()
				continue
			else :
				break #Can't continue on branch


		if resultsWI < (WI_STABLE_LIMIT - MinFloat) : #Mode has become stable
			print("Branch is stable! Stop following")
			if downwardsL_Y : #Start scanning upwards
				downwardsL_Y=False
				currentL_Y = L_Ystart + L_Ystep
				shutil.copyfile('input.gyro.orig', 'input.gyro')
				input_gyro = open('input.gyro', 'a')
				input_gyro.write("\n")
				input_gyro.write("FIELDEIGEN_WR=%E\n" % restartRadiusWR)
				input_gyro.write("FIELDEIGEN_WI=%E\n" % restartRadiusWI)
				input_gyro.close()
				continue
			else :
				break #Can't continue on branch

		if currentL_Y == L_Ystart : #save solution for later
			restartRadiusWR=resultsWR
			restartRadiusWI=resultsWI
		
		shutil.copyfile('input.gyro.orig', 'input.gyro')
		input_gyro = open('input.gyro', 'a')
		input_gyro.write("\n")
		input_gyro.write("FIELDEIGEN_WR=%E\n" % resultsWR)
		input_gyro.write("FIELDEIGEN_WI=%E\n" % resultsWI)
		input_gyro.close()

		if (currentL_Y - L_Ystep) < (L_Ymin - MinFloat) : #Change to scanning upwards
			downwardsL_Y = False
			currentL_Y = L_Ystart
			input_gyro = open('input.gyro', 'a')
			input_gyro.write("FIELDEIGEN_WR=%E\n" % restartRadiusWR)
			input_gyro.write("FIELDEIGEN_WI=%E\n" % restartRadiusWI)
			input_gyro.close()

		if downwardsL_Y :
			currentL_Y = currentL_Y - L_Ystep
		else :
			currentL_Y = currentL_Y + L_Ystep

	if radiusbreak :
		if downwardsRadius :
			downwardsRadius = False
			currentRADIUS = RADIUSstart + RADIUSstep
			shutil.copyfile('input.gyro.orig', 'input.gyro')
			input_gyro = open('input.gyro', 'a')
			input_gyro.write("\n")
			input_gyro.write("FIELDEIGEN_WR=%E\n" % startWR)
			input_gyro.write("FIELDEIGEN_WI=%E\n" % startWI)
			input_gyro.close()
			continue
		else :
			break #Can't continue on branch, stop execution

	shutil.copyfile('input.gyro.orig', 'input.gyro')

	input_gyro = open('input.gyro', 'a')
	input_gyro.write("\n")
	input_gyro.write("FIELDEIGEN_WR=%E\n" % restartRadiusWR)
	input_gyro.write("FIELDEIGEN_WI=%E\n" % restartRadiusWI)
	input_gyro.close()

	if (currentRADIUS - RADIUSstep) < (RADIUSmin - MinFloat) : #Change to scanning upwards
		downwardsRadius = False
		currentRADIUS = RADIUSstart
		input_gyro = open('input.gyro', 'a')
		input_gyro.write("FIELDEIGEN_WR=%E\n" % startWR)
		input_gyro.write("FIELDEIGEN_WI=%E\n" % startWI)
		input_gyro.close()

	if downwardsRadius :
		currentRADIUS = currentRADIUS - RADIUSstep 
	else :
		currentRADIUS = currentRADIUS + RADIUSstep

sys.exit(0)

