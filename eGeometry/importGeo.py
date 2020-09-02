import numpy as np
import os
import sys
from scipy.interpolate import CubicSpline


"""
interpCrv takes a series of points and forms separate equations for the x,y,z coordinates based on index numbers
parameters:
- pts: list of points for function generation
- reso: the number of points evaluated using the equations generated
"""

"""
interpSrf takes a series of contours (interpolated curves) and uses their points to form a grid that defines the UV coordinates of an interpolated srf
parameters:
- contours: interpolated curves
- resolution U: number of points between contours that define the surface grid
- resolution V: number of points along the contours that define the surface grid
- distDivide (optional, default = False): whether or not the surface grid is normalized based on distance 

"""

def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()


class importOBJ:

	def __init__ (self, loc):

		self.v = []
		self.f = []
		self.ft = []
		self.fn = []
		self.vn = []
		self.file = loc

		self.parse()

	def parse(self):

		f = open(self.file,'r')
		lines = f.read().split('\n')

		for i in range(len(lines)):

			line = lines[i].split(' ')

			if 'v' in line:
				point = []
				for j in range(len(line)):
					if j>0: 
						point.append(float(line[j]))
				if len(point)==3:
					self.v.append(point)

			if 'vn' in line:
				norm = []
				for j in range(len(line)):
					if j>0: 
						norm.append(float(line[j]))
				if len(norm)==3:
					self.vn.append(norm)

			if 'f' in line:
				face = []
				ftxt = []
				fnorm = []
				for j in range(len(line)):
					entry = line[j].split('/')
					if len(entry)==3:
						if is_integer(entry[0]): face.append(int(entry[0]))
						if is_integer(entry[1]): fxt.append(int(entry[1]))
						if is_integer(entry[2]): fnorm.append(int(entry[2]))
				if len(face)>2:
					self.f.append(face)
				if len(ftxt)>2:
					self.ft.append(ftxt)
				if len(fnorm)>2:
					self.fn.append(fnorm)

		self.v = np.array(self.v)
		self.vn = np.array(self.vn)
		self.f = np.array(self.f)
		self.ft = np.array(self.ft)
		self.fn = np.array(self.fn)


def importCrvOBJ(loc):

	f = open(loc,'r')
	allLines = f.read().split('\n')

	crvs = []
	verts = []

	broken = False

	for i in range(len(allLines)):

		line = allLines[i].split(' ')

		if line[0] == 'v':

			if '\\' in line:

				pt = [float(line[1]),float(line[2])]
				pt.append(float(allLines[i+1]))

				broken = True

			else:

				pt = [float(line[1]),float(line[2]),float(line[3])]
				verts.append(pt)
				broken = False

		elif len(verts)>0 and broken == False:

			crvs.append(verts)
			verts = []

	return [np.array(crvs),np.array(verts)]


def importCrvsCSV(dir):

	f = open(dir,'r')
	lines = f.read().split('\n')
	profiles = []

	for i in range(len(lines)):
		
		crvPts = []
		pts = lines[i].split(' ')
		
		for j in range(len(pts)):

			point = []
			entry = pts[j].split(',')

			for k in range(len(entry)):
				point.append(float(entry[k]))
			crvPts.append(point)

		profiles.append(np.array(crvPts))

	return profiles
