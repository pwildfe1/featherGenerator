
import eGeometry.interpolateGeo as iG
import eGeometry.exportPolyline as eP
import eGeometry.importGeo as iReader
import random as r
import math as m
import numpy as np
import os
import pandas as pd
import json

"""

instruction class packages all of the parameters needed for each new generation of barbs

INPUT:
angle(deg) =  the angle you want the barb to come off the stem
sprout_ratio(0-1) = the length of the new barb as a ratio of the stem
density(0-1) = the density of further barb locations on this barb
entropy(0-1) = the variation affecting every other train 0-1
curl(0-1)  = the amount the barb curls off of the stem)

"""

class Instruction:


	def __init__(self, angle, sprout_ratio, density, entropy, curl, thinning, stub_ratio):

		self.ang = angle
		self.length = sprout_ratio
		self.density = density
		self.entropy = entropy
		self.curl = curl
		self.thinning = thinning
		self.ratio = stub_ratio


"""

Stem class stores the points, curve and tapering radius data of the previous barb or original path

INPUT:
pts = numpy array of 3D points
radius = the start radius of this completed barb
gradient = the rate of the radius's descent towards at least 1/4 the initial thickness 

"""


class Stem:

	def __init__(self,pts,radius=3,gradient=.5):

		self.pts = pts
		self.crv = iG.interpCrv(self.pts)
		self.t = []
		self.radius = radius
		self.radii = []
		self.length = 0
		self.updateTan()
		self.genRadii(gradient)

	def updateTan(self):

		self.t = []
		self.length = 0

		for i in range(self.pts.shape[0]):

			if i<self.pts.shape[0]-1:
				t = self.pts[i+1] - self.pts[i]
				self.length = self.length + np.linalg.norm(t)
				self.t.append(t/np.linalg.norm(t))

		self.t.append(self.t[-1])
		self.t = np.array(self.t)

	def genRadii(self,gradient=.5):

		start = int(self.pts.shape[0]*.1)

		for i in range(self.pts.shape[0]):

			if i<start:
				section_radius = self.radius
			else:
				f = 1-m.pow((i-start)/(self.pts.shape[0]-start-1), gradient)+.25
				if f > 1: f = 1
				section_radius = self.radius*f

			self.radii.append(section_radius)


"""

Barb class stores stem it is branching off of, the index of the point on the stem and instructions for how to shape
the barbs form

INPUT:
stem(Stem class) = Stem class that the barb is jutting out from
index(int) = The index of the point in the stem that the barb is jutting out from 
reso(default = 200) = The number of points initially determining the form of the barb

"""

class Barb:

	def __init__(self, stem, index, manual, reso=200, threeD=False):

		self.stem = stem
		self.index = index
		self.ang = manual.ang
		self.entropy = manual.entropy
		self.reso = reso
		self.thinning = manual.thinning

		if self.ang<0: 
			reverse = -1
		else:
			reverse = 1

		# barb location on stem
		self.stem_location = self.index/self.stem.pts.shape[0]

		# the radius is taken from the stem at the sprout point and reduced by a factor of .75
		self.radius = self.stem.radii[self.index]*.75

		# the actual angle is randomized greatly by the entropy
		self.ang = self.ang + 2*(r.random()-.5)*self.ang*m.pow(self.entropy,2)*.5

		# the actual barb length is a factor of the stem length, radius and entropy
		self.length = self.stem.length*manual.length*(1.1-m.pow(self.stem_location,3))*(1 + .1*(r.random()-.5)*self.entropy)

		# the actual curl is a fraction of the barb length and affected by entropy aswell
		self.curl = manual.curl*self.length*(1 - 2*(r.random()-.5)*self.entropy)

		otherTan = np.subtract(self.stem.pts[self.index],self.stem.pts[0]) # vector based on stem plane
		self.plane = np.cross(otherTan,self.stem.t[self.index]) #cross two vectors on the stem plane to get the stem plane 
		
		# the stem plane helps remove the need for a fixed axis, by basing branch movement on the relative plane of the stem
		self.plane = self.plane/np.linalg.norm(self.plane)

		# x-axis for this barb is along the barbs direction (rotating self.ang from the stem)
		self.x = np.array(vecRotate(self.stem.t[self.index],self.ang,self.plane))
		# y-axis for this barb is rotated -90 degrees from x-axis based on the stem curves plane
		self.y = np.array(vecRotate(self.x,-90,self.plane))

		if np.dot(self.y,self.stem.t[self.index]) < 0:
			self.y = -self.y

		##### ROTATES THE BARB IN THE 3RD DIMENSION
		
		torque = r.random()

		self.torsion = self.ang*r.random()

		if torque<.5:
			self.torsion = -self.torsion

		if threeD:
			self.x = np.array(vecRotate(self.x,self.torsion,self.stem.t[self.index]))
			self.y = np.array(vecRotate(self.y,self.torsion,self.stem.t[self.index])) 
		
		#####

		self.pts = [] # will store the final points of the barb
		self.t = [] # will store the tangents of the barb
		self.setPts()

	####
	# setPts: finalizes the position of the points and the barb by applying a tappering off the stem and a wave eq
	####

	def setPts(self,reso=200):

		self.pts = []
		st = self.stem.pts[self.index]
		glide = self.stem.t[self.index]*self.stem.radii[self.index]

		for i in range(self.reso):

			pt = st + self.x*self.length/(self.reso-1)*i
			f = m.pow(i/(self.reso/10),.5)
			if f<1:
				pt = pt + glide*f
			else:
				pt = pt + glide

			self.pts.append(pt)

		self.pts = np.array(self.pts)
		self.wavePoints()
		self.updateTan()

	####
	# wavePoints: moves the points with a wave equation to add curl and flow to hard straight lines.
	# The wavelength is kept below 1.25 to ensure no wonky waves and for the most part a simple arc
	####

	def wavePoints(self,dir=1):

		waveLength = self.entropy*1.5 #the entropy results in a wavelength that is greater and therefore a more complete arc
		newPts = []

		for i in range(len(self.pts)):
			# the greater the curl the more exaggerated the arc
			f = self.curl*m.sin(2*m.pi*waveLength*i/len(self.pts))
			# y = self.stem.t[self.index]*f
			y = self.y*f
			self.pts[i] = self.pts[i]+y


	####
	# updateTan: updates the tangents of the barb to match the points
	####

	def updateTan(self):

		self.t = []

		for i in range(self.pts.shape[0]):

			if i<self.pts.shape[0]-1:
				t = self.pts[i+1] - self.pts[i]
				self.t.append(t/np.linalg.norm(t))

		self.t.append(-self.t[-1])
		self.t = np.array(self.t)


	####
	# genPath: creates the final Stem based on the barbs fully placed points
	####

	def genPath(self):

		self.path = Stem(self.pts,self.radius,self.thinning)

		return self.path




####
# INITIALIZE: places new barbs on an existing path (Stem)
# path (Stem class) = Stem that the barbs will branch off of 
# manual (Instruction class) = stores all the data for the new barbs
#
# returns = [the new barbs contributing to next generations, the new barbs not contributing to next generations]
####


def initialize(path, manual, reso=100):

	newBarbs = []
	stubs = []

	b_interval = int(1/(manual.density+.0001)) # manual.density translates to b_interval to place barbs
	density = manual.density # stores original barb density

	for i in range(path.pts.shape[0]):

		index = i + int(path.pts.shape[0]/5*manual.entropy)

		if i%b_interval==0 and i/path.pts.shape[0]>.1 and i/path.pts.shape[0]<.9:

			# decreases barb density as we move further up the stem
			manual.density = density*(1-.0125*m.pow(i/path.pts.shape[0],3))

			# makes first side barb
			newBarbs.append(Barb(path,i,manual,reso))

			# switches up angle for mirrored barb
			manual.ang = -manual.ang
			
			# makes second side barb
			newBarbs.append(Barb(path,i,manual,reso))

			# returns angle to original side
			manual.ang = -manual.ang

		# if there are stubs form another set of barbs with 0 density to keep them from contributing to growth

		if manual.ratio!=0 and i%b_interval!=0 and i/path.pts.shape[0]>.1 and i%4==0:

				stub_manual = Instruction(50, manual.ratio, 0, 0, 0, 1, 0) # 0 density prevents further growth
				stubs.append(Barb(path,i,stub_manual,reso))

				stub_manual = Instruction(-50, manual.ratio, 0, 0, 0, 1, 0) # 0 density prevents further growth
				stubs.append(Barb(path,i,stub_manual,reso))

	return [newBarbs,stubs]


def vecRotate(vec,ang,axis):

	cos = m.cos(m.pi/180*ang)
	sin = m.sin(m.pi/180*ang)
	v = vec
	u = [axis[0]/np.linalg.norm(axis),axis[1]/np.linalg.norm(axis),axis[2]/np.linalg.norm(axis)]
	R1,R2,R3 = [] , [] , []
	c = 1-cos

	R1.append(cos+m.pow(u[0],2)*c)
	R1.append(u[0]*u[1]*c-u[2]*sin)
	R1.append(u[0]*u[2]*c+u[1]*sin)

	R2.append(u[1]*u[0]*c+u[2]*sin)
	R2.append(cos+m.pow(u[1],2)*c)
	R2.append(u[1]*u[2]*c-u[0]*sin)

	R3.append(u[2]*u[0]*c-u[1]*sin)
	R3.append(u[2]*u[1]*c+u[0]*sin)
	R3.append(cos+m.pow(u[2],2)*c)

	x = v[0]*R1[0] + v[1]*R1[1] + v[2]*R1[2]
	y = v[0]*R2[0] + v[1]*R2[1] + v[2]*R2[2]
	z = v[0]*R3[0] + v[1]*R3[1] + v[2]*R3[2]

	return [x,y,z]


####
# WRITEPATHS: writes every path's information to a .csv file including chaning radius information
# paths (list) = array of Stem class objects 
# dest (str) = the location of the resulting .csv file
####

def writePaths(paths,dest='RHINO/all_paths.csv'):

	f = open(dest,'w')

	for i in range(len(paths)):

		for j in range(paths[i].pts.shape[0]):

			f.write(str(paths[i].pts[j][0]))
			f.write(',')
			f.write(str(paths[i].pts[j][1]))
			f.write(',')
			f.write(str(paths[i].pts[j][2]))
			f.write(',')
			f.write(str(paths[i].radii[j]))

			if j<paths[i].pts.shape[0]-1:

				f.write(' ')

		if i < len(paths)-1:

			f.write('\n')

	f.close()



def Main():

	path_pts = iReader.importCrvsCSV('RHINO/paths.csv')[0]
	rachis = Stem(path_pts)

	# Correct first path radius to normalize for any length

	initialRadius = rachis.length/90 # decrease the denominator to increase the radius
	rachis = Stem(path_pts,radius=initialRadius)
	
	#

	manuals = []
	generations = 4

	#### HYPER-PARAMETERS ####

	# each list represents a different trait and each entry in the list represents that trait for that generation

	angles =       [60, 90, 70, 60]
	long_ratio =   [.2, .5, .5, .5]
	barb_density = [.03, .03, .04, .05]
	entropy =      [.7, .5, .6, .6]
	curl =         [.075, .1, .1, .1]
	thinning =     [1, .75, .5, .5]
	short_ratio =  [.05, 0, 0, 0]

	####################


	 
	# loop sets up the instructions for every generation based on user input above

	for n in range(generations):

		step = Instruction(angles[n], long_ratio[n], barb_density[n], entropy[n], curl[n], thinning[n], short_ratio[n])
		manuals.append(step)

	all_barbs = [] # stores all barbs
	all_paths = [rachis] # stores all stems

	curr_barbs = [] # stores current generation barbs
	curr_paths = [rachis] # stores current generation stems

	for n in range(generations):

		# loop initalizes new barbs based on every current path

		for i in range(len(curr_paths)):

			barbs = initialize(curr_paths[i], manuals[n])
			curr_barbs.extend(barbs[0])

			# save short barbs path information separately because the barb  will not be used for next generations

			for j in range(len(barbs[1])):
				all_paths.append(barbs[1][j].genPath())

		curr_paths = []

		# loop records new paths based on every new barb in this generation

		for i in range(len(curr_barbs)):

			curr_paths.append(curr_barbs[i].genPath())
			all_paths.append(curr_paths[-1])
			all_barbs.append(curr_barbs[i])

		curr_barbs = []

		# saves a checkpoint at every generation so work can be resumed
		writePaths(all_paths, 'RHINO/generation_0' + str(n) + '.csv') 

	writePaths(all_paths) # writes all path information, including radius information, to 'RHINO' folder




Main()