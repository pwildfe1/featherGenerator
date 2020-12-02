import numpy as np
import os
import sys
from scipy.interpolate import CubicSpline


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


def genCircleSection(cnt,vector,radius,axis=[0,0,1],reso=8):
    
	pts = []
	vec = vector/np.linalg.norm(vector)
	v = vecRotate(vec,90,axis)

	for i in range(reso):
		mv = np.array(vecRotate(v,(360/reso)*i,vector))
		mv = mv * radius
		pts.append(cnt + mv)
	pts.append(pts[0])

	return pts



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


class interpCrv:

	def __init__ (self, pts, reso = 200):

		self.data = pts
		self.pts = pts
		self.tans = []
		self.reso = reso
		self.crvLen = self.updateLength()
		self.updateCrv()


	def updateCrv(self):

		indexes = np.arange(self.pts.shape[0])

		if np.linalg.norm(self.pts[0]-self.pts[-1]) == 0:

			self.fx = CubicSpline(indexes,self.pts[:,0],bc_type='periodic')
			self.fy = CubicSpline(indexes,self.pts[:,1],bc_type='periodic')
			self.fz = CubicSpline(indexes,self.pts[:,2],bc_type='periodic')

		else:

			self.fx = CubicSpline(indexes,self.pts[:,0])
			self.fy = CubicSpline(indexes,self.pts[:,1])
			self.fz = CubicSpline(indexes,self.pts[:,2])

		self.cntPoint()


	def updateLength(self):

		crvLen = 0

		for i in range(len(self.pts)-1):

			vec = np.subtract(np.array(self.pts[i]),np.array(self.pts[i+1]))
			crvLen = crvLen + np.linalg.norm(vec)
		
		self.crvLen = crvLen

		return self.crvLen

	# updateTans() recalculates the tangents derivatives at the spline control points

	def updateTans(self):

		#self.updateCrv()

		indexes = np.arange(self.pts.shape[0])

		tans = []

		for i in range(len(indexes)):

			tan = [self.fx(indexes[i],1),self.fy(indexes[i],1),self.fz(indexes[i],1)]
			mag = np.linalg.norm(np.array(tan))
			tans.append([tan[0]/mag, tan[1]/mag, tan[2]/mag])

		self.tans = np.array(tans)


	# genPts(param) function generates the curve based a number of points (reso). Each dimension gets its own spline equation [f(x),f(y),f(z)] 


	def genPts(self,reso):

		step = (len(self.data)-1)/reso
		indexes = np.arange(0,step*reso+step,step)

		pts = []
		tans = []

		for i in range(len(indexes)):

			pts.append([self.fx(indexes[i]),self.fy(indexes[i]),self.fz(indexes[i])])
			tan = [self.fx(indexes[i],1),self.fy(indexes[i],1),self.fz(indexes[i],1)]
			mag = np.linalg.norm(np.array(tan))
			tans.append([tan[0]/mag, tan[1]/mag, tan[2]/mag])

		self.pts = np.array(pts)
		self.tans = np.array(tans)

		return pts


	# genUniform(param) function generates the curve based a number of points (reso) that are equidistant from one another on the overall curve (approximately)
	# before this function can be called the self.pts array must be created.
	# parameters:
	# param = value from 0-1 determining where along the curve you want to evaluate a point

	def genUniform(self,segments):

		pts = []

		print(segments)

		for i in range(segments+1):

			param = self.evalDistParam(i/segments)
			pts.append(self.evalParam(param))

		self.pts = np.array(pts)

		self.updateTans()

		return pts


	# evalDistParam(normParam) function gets point closest to normalized parameter on curve
	# before this function can be called the self.pts array must be created.
	# parameters:
	# param = value from 0-1 determining where along the curve you want to evaluate a point


	def evalDistParam(self,normParam):

		self.updateLength()

		dist = self.crvLen*normParam

		param = 0

		total = 0

		for i in range(len(self.pts)-1):

			vec = np.subtract(self.pts[i],self.pts[i+1])

			total = np.linalg.norm(vec) + total
			
			if total > dist:

				diff = total - dist
				mag = np.linalg.norm(vec)
				factor = diff/mag
				if factor>1:
					print('shit')
				param = (i + factor)/len(self.pts)
				break

			elif total == dist:

				param = i/len(self.pts)
				break


		return param

	#evalParam points retrieves a point on the spline based on the curves existing x y and z functions. If you want to evaluate based on a unitized
	#see genDistPts

	def evalParam(self,param):

		index = param*(len(self.data)-1)

		result = [self.fx(index),self.fy(index),self.fz(index)]

		return result

	#cntPoint finds the approximate center point of the curve based on the cntrl points

	def cntPoint(self):

		points = self.genPts(self.reso)
		sum = [0,0,0]

		for i in range(len(points)):

			for j in range(len(points[i])):

				sum[j] = sum[j]+points[i][j]

		self.cnt = np.array(sum)*1/len(points)

		return self.cnt

	#closestPt(pt) finds the approx location and index of the closest point on the curve to a point given by the user
	#pt = test point to find closest point to.

	def closestPt(self,pt):

		points = self.genPts(self.reso)
		index = 0
		minimum = np.linalg.norm(np.array(points[0]) - pt)

		#print(np.full((self.pts.shape),pt))

		for i in range(len(points)):

			dist = np.linalg.norm(np.array(points[i]) - pt)
			if dist<minimum:
				minimum = dist
				index = i

		return [index,np.array(points[index])]

	#recordCrv(dir) writes the point information that defines the curve to a .csv file at a specified directory
	#dir = output directory

	def recordCrv(self,dir):

		#self.genPts(200)
		f = open(dir,'w')

		for i in range(len(self.pts)):
			for j in range(len(self.pts[i])):
				f.write(str(self.pts[i][j]))
				if j<len(self.pts[i])-1:
					f.write(',')
			if i<len(self.pts)-1:
				f.write('\n')

		f.close()

	#mvCrv(vec) moves the entire curve by a user specified vector
	#vec = vector you want to move the entire curve by

	def mvCrv(self,vec):

		for i in range(len(self.pts)):
			self.pts[i] = np.array([self.pts[i][0]+vec[0],self.pts[i][1]+vec[1],self.pts[i][2]+vec[2]])


	#rotCrv(ang,origin,axi) rotates an entire curve around an origin point and axis by a specified angle
	#ang = angle to rotate curve by
	#origin = center of rotation
	#axi = axis you are rotating the curve around


	def rotCrv(self,ang,origin,axi):
		
		for i in range(len(self.pts)):
			vec = np.subtract(self.pts[i],origin)
			vec = np.array(vecRotate(vec,ang,axi))
			self.pts[i] = np.add(vec,origin)


	#offsetCrv(plane,value) offsets a curve in a specified plane by a certain value
	#plane = a numpy vector defining the plane normal
	#value = the magnitude the curve will be offset by (can be - or + to offset in and out)


	def offsetCrv(self,plane,value):

		oldVec = np.cross(self.tans[0],plane)
		oldVec = oldVec/np.linalg.norm(oldVec)

		for i in range(len(self.pts)):

			vec = np.cross(self.tans[i],plane)
			vec = vec/np.linalg.norm(vec)

			if np.dot(vec,oldVec)<0:
				factor = -value
			else:
				factor = value
			
			self.pts[i] = np.add(self.pts[i],vec*factor)
			oldVec = vec

		return self.pts


	#scaleCrv(origin,factor) scales the entire curve from an origin point by a specified factor in the [x,y,z] directions
	#origin = origin of scale
	#factor = list of three values defining how much the curve will be scaled in the x,y,z directions


	def scaleCrv(self,origin,factor):

		vecs = []

		for i in range(len(self.pts)):

			vecs.append([self.pts[i][0]-origin[0],self.pts[i][1]-origin[1],self.pts[i][2]-origin[2]])
			vecs[-1] = [vecs[-1][0]*factor[0], vecs[-1][1]*factor[1], vecs[-1][2]*factor[2]]

			self.pts[i] = np.add(origin,vecs[i])

		return self.pts


	#reverseCrv() reverses the direction of the curve (start point and end point)


	def reverseCrv(self):

		newPts = []

		for i in range(len(self.pts)):

			newPts.append(self.pts[len(self.pts)-1-i])

		self.pts = newPts


"""
interpSrf takes a series of contours (interpolated curves) and uses their points to form a grid that defines the UV coordinates of an interpolated srf
parameters:
- contours: interpolated curves
- resolution U: number of points between contours that define the surface grid
- resolution V: number of points along the contours that define the surface grid
- distDivide (optional, default = False): whether or not the surface grid is normalized based on distance 

"""

class interpSrf:


	def __init__(self,contours,resolutionU,resolutionV,distDivide = False):

		self.vCrv = contours
		self.uCrv = []
		self.resoU = resolutionU
		self.resoV = resolutionV
		self.divByDist = distDivide

		self.U = []
		self.V = []
		self.faces = []
		self.v = []
		self.vn = []
		self.vt = []
		self.utan = []
		self.vtan = []

		crosses = []
		columns = []

		# START: divides the original contours into a a set of resolution V points

		# step = self.vCrv[0].crvLen/self.resoV

		for i in range(len(self.vCrv)):

			if self.divByDist == True:

				columns.append(self.vCrv[i].genUniform(self.resoV))

			else:

				columns.append(self.vCrv[i].genPts(self.resoV-1))

		rowPts = []

		

		for i in range(len(columns[0])):

			cross = []

			for j in range(len(columns)):

				cross.append(columns[j][i])

			crosses.append(interpCrv(np.array(cross)))
			rowPts.append(cross)

		# END


		# START: divides the original contours into a a set of resolution U points


		for i in range(self.resoV):

			#step = crosses[i].crvLen/self.resoU

			row = crosses[i].genPts(self.resoU)

			#row = crosses[i].genUniform(step)

			self.U.append(row)

		# END

		self.genVCoord()
		self.boundingBox()


	# updates the V coordinates 


	def genVCoord(self):

		self.V = []

		for i in range(len(self.U[0])):

			row = []

			for j in range(len(self.U)):

				row.append(self.U[j][i])

			self.V.append(row)

	# calculates the bounding box around a interpolated surface. [lower corners 1-4 counterclockwise then upper corners 1-4 counterclockwise]

	def boundingBox(self):


		minX,minY,minZ = self.U[0][0][0],self.U[0][0][1],self.U[0][0][2]
		maxX,maxY,maxZ = self.U[0][0][0],self.U[0][0][1],self.U[0][0][2]

		for i in range(len(self.U)):
			for j in range(len(self.U[i])):

				pt = self.U[i][j]

				if pt[0]<minX: minX = pt[0]
				if pt[0]>maxX: maxX = pt[0]
				if pt[1]<minY: minY = pt[1]
				if pt[1]>maxY: maxY = pt[1]
				if pt[2]<minZ: minZ = pt[2]
				if pt[2]>maxZ: maxZ = pt[2]

		self.bounds = [[minX,minY,minZ],[maxX,minY,minZ],[maxX,maxY,minZ],[minX,maxY,minZ]]
		self.bounds.extend([[minX,minY,maxZ],[maxX,minY,maxZ],[maxX,maxY,maxZ],[minX,maxY,maxZ]])


	# correctAlign(thres,level,align): Aligns all the surface verticies within a specified tolerance to a defined plane in the X,Y or Z direction
	# thres = tolerance between vertex location and the level that you want to align
	# level = the X,Y,Z value at which vertices will be corrected
	# align =  X Y or Z (what plane you want to align to)

	def correctAlign(self,thres,level,align):

		correctedLocs = []

		for i in range(len(self.U)):
			for j in range(len(self.U[i])):

				if abs(self.U[i][j][align] - level)<thres:
					
					self.U[i][j][align] = self.bounds[0][align]

		self.genVCoord()
		self.boundingBox()

	def swapUV(self):

		self.U = self.V

		resoU = self.resoU
		self.resoU = self.resoV
		self.resoV = resoU

		self.genVCoord()
		self.boundingBox()


	def recordPts(self):

		f = open('gridPts.csv','w')

		for i in range(len(self.U)):
			for j in range(len(self.U[i])):

				f.write(str(self.U[i][j][0]))
				f.write(',')
				f.write(str(self.U[i][j][1]))
				f.write(',')
				f.write(str(self.U[i][j][2]))
				f.write('\n')

		f.close()

	# createFaces() generates the faces, vertex normals and vertex texture values that would be used to export a mesh 

	def createFaces(self):

		for i in range(len(self.U)):
			for j in range(len(self.U[i])):
				
				self.v.append(self.U[i][j])

				if i<len(self.U)-1:
					vecX  = np.subtract(self.U[i+1][j],self.U[i][j])
				else:
					vecX = np.subtract(self.U[i][j],self.U[i-1][j])

				if j<len(self.U[i])-1:
					vecY = np.subtract(self.U[i][j+1],self.U[i][j])
				else:
					vecY = np.subtract(self.U[i][j],self.U[i][j-1])

				
				vNorm = np.cross(vecX,vecY)

				self.vn.append(np.multiply(vNorm,1/np.linalg.norm(vNorm)))

				self.vt.append(np.array([i/(len(self.U)-1),j/(len(self.U[i])-1)]))

				self.utan.append(vecX/np.linalg.norm(vecX))

				self.vtan.append(vecY/np.linalg.norm(vecY))


		for i in range(len(self.U)-1):
			for j in range(len(self.U[i])):

				if j>0:

					face = []
					face.append(i*len(self.U[i])+j)
					face.append(i*len(self.U[i])+j+1)
					face.append((i+1)*len(self.U[i])+j+1)
					face.append((i+1)*len(self.U[i])+j)

					self.faces.append(face)

	# exportMesh(output) uses a standard .obj format to export the interpsrf as a quad mesh
	# output = the output .obj file

	def exportMesh(self,output):

		self.createFaces()

		out = open(output,'w')

		out.write('# Rhino')
		out.write('\n')
		out.write('\n')
		out.write('o object_1')
		out.write('\n')

		for i in range(len(self.v)):

			out.write('v ')
			out.write(str(self.v[i][0]))
			out.write(' ')
			out.write(str(self.v[i][1]))
			out.write(' ')
			out.write(str(self.v[i][2]))
			out.write('\n')

		for i in range(len(self.v)):

			out.write('vt ')
			out.write(str(self.vt[i][0]))
			out.write(' ')
			out.write(str(self.vt[i][1]))
			out.write('\n')

		for i in range(len(self.v)):

			out.write('vn ')
			out.write(str(self.vn[i][0]))
			out.write(' ')
			out.write(str(self.vn[i][1]))
			out.write(' ')
			out.write(str(self.vn[i][2]))
			out.write('\n')

		for i in range(len(self.faces)):

			out.write('f ')
			out.write(str(self.faces[i][0]))
			out.write(' ')
			out.write(str(self.faces[i][1]))
			out.write(' ')
			out.write(str(self.faces[i][2]))
			out.write(' ')
			out.write(str(self.faces[i][3]))
			out.write('\n')

		out.close()


	# offsetSrf(thickness) offsets the U,V grid based on the vertex normals and a specified thickness that can be (- or +)


	def offsetSrf(self,thickness,solid=True):

		self.createFaces()

		U,V = [],[]

		srfCrvs = []
		faces = []
		norms = []
		tans = []
		vertNum = len(self.v)

		for i in range(len(self.U)):
			
			row = []

			for j in range(len(self.U[i])):

				index = i*len(self.U[i])+j

				pt = np.add(self.U[i][j],list(self.vn[index]*thickness))

				if solid:

					self.v.append(pt)
					self.vn.append(list(np.multiply(self.vn[index],-1)))
					self.vt.append(self.vt[index])

					if i==0 and j<len(self.U[i])-1:

						faces.append([index+1,index+2,index+2+vertNum,index+vertNum+1])

					if i==len(self.U)-1 and j<len(self.U[i])-1:

						faces.append([index+1,index+2,index+2+vertNum,index+vertNum+1])

					if j==0 and i<len(self.U)-1:

						faces.append([index+1,index+len(self.U[i])+1,index+len(self.U[i])+vertNum+1,index+vertNum+1])

					if j==len(self.U[i])-1 and i<len(self.U)-1:

						faces.append([index+1,index+len(self.U[i])+1,index+len(self.U[i])+vertNum+1,index+vertNum+1])

				else:

					self.v[index] = pt
					self.vn[index] = list(np.multiply(self.vn[index],-1))
					self.vt[index] = self.vt[index]


				row.append(pt)

			U.append(row)
			srfCrvs.append(interpCrv(np.array(row)))

		if solid:

			faceNum = len(self.faces)

			for i in range(faceNum):

				newFace = []

				for j in range(len(self.faces[i])):

					newFace.append(self.faces[i][j]+vertNum)

				if newFace not in self.faces:

					self.faces.append(newFace)

			for i in range(len(faces)):

				if faces[i] not in self.faces:

					self.faces.append(faces[i])

		temp = []

		for i in range(len(self.faces)):

			if self.faces[i] not in temp:

				temp.append(self.faces[i])
		
		self.faces = temp


		return srfCrvs


	def offsetVerts(self,thickness, dim3 = True):

		midU = int(len(self.U)/2)

		baseVerts = len(self.U)*len(self.U[0])

		for i in range(len(self.U)):

			midV = int(len(self.U[i])/2)

			for j in range(len(self.U[i])):

				index = i*len(self.U[i])+j
				index02 = index + baseVerts

				if dim3:

					self.v[index] = list(np.add(self.v[index],np.multiply(self.vn[index],thickness)))
					self.v[index02] = list(np.add(self.v[index02],np.multiply(self.vn[index02],thickness)))

				if j<midV:

					vec = np.subtract(self.U[i][j+1],self.U[i][j])
					vec = vec/np.linalg.norm(vec)
					vec = vec*(j-midV)/midV*thickness
					self.v[index] = list(np.add(self.v[index],vec))

					if dim3:
						self.v[index02] = list(np.add(self.v[index02],vec))

				if j>midV:

					vec = np.subtract(self.U[i][j-1],self.U[i][j])
					vec = vec/np.linalg.norm(vec)
					vec = vec*(midV-j)/midV*thickness
					self.v[index] = list(np.add(self.v[index],vec))

					if dim3:
						self.v[index02] = list(np.add(self.v[index02],vec))

				if i<midU:

					vec = np.subtract(self.U[i+1][j],self.U[i][j])
					vec = vec/np.linalg.norm(vec)
					vec = vec*(j-midU)/midU*thickness
					self.v[index] = list(np.add(self.v[index],vec))

					if dim3:
						self.v[index02] = list(np.add(self.v[index02],vec))

				if i>midU:

					vec = np.subtract(self.U[i-1][j],self.U[i][j])
					vec = vec/np.linalg.norm(vec)
					vec = vec*(j-midU)/midU*thickness
					self.v[index] = list(np.add(self.v[index],vec))

					if dim3:
						self.v[index02] = list(np.add(self.v[index02],vec))

				#self.U[i][j] = self.v[index]

		return self.v

	def createBoxes(self,atts,hAtts,dimU,dimV,dimW,height,limit=10,limitH=15,closed=True):

		u = []
		v = []
		w = []

		columns = []
		cols = []
		rows = []
		crosses = []

		layers = []

		for i in range(len(self.U)):

			columns.append(interpCrv(np.array(self.U[i])))
			cols.append(columns[-1].genPts(dimV))

		for i in range(len(cols[0])):
			cross = []
			for j in range(len(cols)):
				cross.append(cols[j][i])
			crosses.append(interpCrv(np.array(cross)))

		for i in range(len(crosses)):
			row = crosses[i].genPts(dimV)
			rows.append(row)


		for n in range(dimW):
			
			nxtRows = []

			for i in range(len(cols)):
				
				nxtRow = []
				prevN = np.array([0,0,0])

				for j in range(len(cols[i])):

					if j==len(cols[i])-1:
						vec01 = np.subtract(cols[i][j],cols[i][j-1])
					else:
						vec01 = np.subtract(cols[i][(j+1)%len(cols[i])],cols[i][j])

					if i==len(cols)-1:
						vec02 = np.subtract(cols[i][j],cols[i-1][j])
					else:
						vec02 = np.subtract(cols[(i+1)%len(cols)][j],cols[i][j])
					
					norm = np.cross(vec01,vec02)

					if i==len(cols)-1 or i==0:
						norm[2] = 0

					if j!=0 and np.dot(prevN,norm)<-.7:
						norm = -norm

					prevN = norm/np.linalg.norm(norm)
					norm = norm/np.linalg.norm(norm)

					if n>1:
						rise = 2.5*height
					else:
						rise = height

					nxtRow.append(np.add(cols[i][j],norm*n*rise/(dimW-1)))

				nxtRows.append(nxtRow)

			layers.append(nxtRows)


		cells = []
		upperCells = []
		finalCells = []
		cellAreas = []
		crvs = []
		num = 0

		if self.closed:
			closeAdjustment = 0
		else:
			closeAdjustment = -1
		

		for i in range(len(layers)):
			for j in range(len(layers[i])+closeAdjustment):
				for k in range(len(layers[i][j])-1):

					pt01 = layers[i][j][k]
					pt02 = layers[i][j][k+1]
					pt03 = layers[i][(j+1)%len(layers[i])][k+1]
					pt04 = layers[i][(j+1)%len(layers[i])][k]
					botPts = [pt01,pt02,pt03,pt04]

					pt05 = layers[(i+1)%len(layers)][j][k]
					pt06 = layers[(i+1)%len(layers)][j][k+1]
					pt07 = layers[(i+1)%len(layers)][(j+1)%len(layers[i])][k+1]
					pt08 = layers[(i+1)%len(layers)][(j+1)%len(layers[i])][k]
					topPts = [pt05,pt06,pt07,pt08]

					cellAreas.append(np.linalg.norm(np.subtract(pt01,pt02))*np.linalg.norm(np.subtract(pt04,pt01)))

					newBotPts = []
					newTopPts = []
					rises = []

					for p in range(len(botPts)):

						pt = botPts[p]
						factor = 0
						count = 0

						for q in range(len(hAtts)):

							attPt = hAtts[q].closestPt(pt)[1]
							vec = np.subtract(attPt,pt)

							if np.linalg.norm(vec)<limitH:
								factor = factor + 1-np.linalg.norm(vec)/limitH
								count = count + 1

						if count>0:
							factor = factor/count

						riseNorm = np.subtract(topPts[p],pt)

						rises.append(riseNorm*factor*2)
						newBotPts.append(np.add(rises[-1],pt))
						newTopPts.append(np.add(rises[-1],topPts[p]))


					botPts = newBotPts
					topPts = newTopPts


					if i>0:

						newTopPts = []
						valid = False
						edges = []
						rises = []

						for p in range(len(botPts)):

							factor = 0
							count = 0

							pt = botPts[p]

							for q in range(len(atts)):
								attPt = atts[q].closestPt(pt)[1]
								vec = np.subtract(attPt,pt)

								if np.linalg.norm(vec)<limit:
									factor = factor + np.linalg.norm(vec)/limit
									count = count + 1
									valid = True

							if count>0:
								factor = 1 - factor/count

							
							edges.append(np.subtract(botPts[(p+1)%len(botPts)],pt)*.5)
							rises.append(np.subtract(topPts[p],pt)*factor)
							newTopPts.append(np.add(rises[-1],pt))


						if valid:

							cell = []
							cell.extend(botPts)
							cell.extend(newTopPts)
							crvs.append(cell)
							#finalCells.append(cell)
							upperCells.append(cell)

					else:

						cell = []
						cell.extend(botPts)
						cell.extend(topPts)
						crvs.append(cell)
						finalCells.append(cell)

					print(str(num/((len(layers))*len(layers[i])*(len(layers[i][j])))*100) + ' %')
					num = num + 1

		self.finalCells = finalCells
		self.cellMembers = crvs
		self.upperCells = upperCells
		self.cellU = dimU
		self.cellV = dimV


	def divideReso(self,domU,domV,thresX,thresY,loopLimit=4):

		for n in range(loopLimit):

			allSubCells = []
			maxCellX = 0
			maxCellY = 0

			for i in range(len(self.finalCells)):

				midB , midT , cntB , cntT= [] , [] , np.array([0,0,0]) , np.array([0,0,0])

				for j in range(int(len(self.finalCells[i])/2)):

					midB.append(np.add(self.finalCells[i][j],self.finalCells[i][(j+1)%4])/2)
					cntB = cntB + midB[-1]
					index = j+5
					if j == 3:
						index = 4
					midT.append(np.add(self.finalCells[i][j+4],self.finalCells[i][index])/2)
					cntT = cntT + midT[-1]

				distX = np.linalg.norm(np.subtract(self.finalCells[i][4],self.finalCells[i][7]))
				distY = np.linalg.norm(np.subtract(self.finalCells[i][4],self.finalCells[i][5]))

				cntT = cntT/4
				cntB = cntB/4

				divideX = False
				divideY = False


				if i%self.cellV<self.cellV*domV[1] and i%self.cellV>=self.cellV*domV[0] and n==0:

					if int(i/self.cellV)<self.cellU*domU[1] and int(i/self.cellV)>self.cellU*domU[0]:

						divideX = True
						divideY = True

				if distX > thresX: 
					divideX = True
				if distY > thresY: 
					divideY = True

				if distX > maxCellX: 
					maxCellX = distX
				if distY > maxCellY: 
					maxCellY = distY


				if divideX and divideY:

					subCells = []

					bot = [self.finalCells[i][0],midB[0],cntB,midB[3]]
					top = [self.finalCells[i][4],midT[0],cntT,midT[3]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[0],self.finalCells[i][1],midB[1],cntB]
					top = [midT[0],self.finalCells[i][5],midT[1],cntT]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [cntB,midB[1],self.finalCells[i][2],midB[2]]
					top = [cntT,midT[1],self.finalCells[i][6],midT[2]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[3],cntB,midB[2],self.finalCells[i][3]]
					top = [midT[3],cntT,midT[2],self.finalCells[i][7]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					allSubCells.extend(subCells)

				if divideX and divideY==False:

					subCells = []

					bot = [self.finalCells[i][0],self.finalCells[i][1],midB[1],midB[3]]
					top = [self.finalCells[i][4],self.finalCells[i][5],midT[1],midT[3]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[1],self.finalCells[i][2],self.finalCells[i][3],midB[3]]
					top = [midT[1],self.finalCells[i][6],self.finalCells[i][7],midT[3]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					allSubCells.extend(subCells)


				if divideY and divideX==False:

					subCells = []

					bot = [self.finalCells[i][0],midB[0],midB[2],self.finalCells[i][3]]
					top = [self.finalCells[i][4],midT[0],midT[2],self.finalCells[i][7]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[0],self.finalCells[i][1],self.finalCells[i][2],midB[2]]
					top = [midT[0],self.finalCells[i][5],self.finalCells[i][6],midT[2]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					allSubCells.extend(subCells)

				if divideX==False and divideY==False:

					allSubCells.append(self.finalCells[i])

			self.finalCells = allSubCells

			allSubCells = []

			print(maxCellX)
			print(maxCellY)

			if maxCellX<thresX and maxCellY<thresY: 
				break

		return self.finalCells


	def divideResoUpper(self,thresX,thresY,loopLimit=4):

		for n in range(loopLimit):

			allSubCells = []
			maxCellX = 0
			maxCellY = 0

			for i in range(len(self.upperCells)):

				midB , midT , cntB , cntT= [] , [] , np.array([0,0,0]) , np.array([0,0,0])

				for j in range(int(len(self.upperCells[i])/2)):

					midB.append(np.add(self.upperCells[i][j],self.upperCells[i][(j+1)%4])/2)
					cntB = cntB + midB[-1]
					index = j+5
					if j == 3:
						index = 4
					midT.append(np.add(self.upperCells[i][j+4],self.upperCells[i][index])/2)
					cntT = cntT + midT[-1]

				distX = np.linalg.norm(np.subtract(self.upperCells[i][4],self.upperCells[i][7]))
				distY = np.linalg.norm(np.subtract(self.upperCells[i][4],self.upperCells[i][5]))

				cntT = cntT/4
				cntB = cntB/4

				divideX = False
				divideY = False

				if distX > thresX: 
					divideX = True
				if distY > thresY: 
					divideY = True

				if distX > maxCellX: 
					maxCellX = distX
				if distY > maxCellY: 
					maxCellY = distY


				if divideX and divideY:

					subCells = []

					bot = [self.upperCells[i][0],midB[0],cntB,midB[3]]
					top = [self.upperCells[i][4],midT[0],cntT,midT[3]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[0],self.upperCells[i][1],midB[1],cntB]
					top = [midT[0],self.upperCells[i][5],midT[1],cntT]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [cntB,midB[1],self.upperCells[i][2],midB[2]]
					top = [cntT,midT[1],self.upperCells[i][6],midT[2]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[3],cntB,midB[2],self.upperCells[i][3]]
					top = [midT[3],cntT,midT[2],self.upperCells[i][7]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					allSubCells.extend(subCells)

				if divideX and divideY==False:

					subCells = []

					bot = [self.upperCells[i][0],self.upperCells[i][1],midB[1],midB[3]]
					top = [self.upperCells[i][4],self.upperCells[i][5],midT[1],midT[3]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[1],self.upperCells[i][2],self.upperCells[i][3],midB[3]]
					top = [midT[1],self.upperCells[i][6],self.upperCells[i][7],midT[3]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					allSubCells.extend(subCells)


				if divideY and divideX==False:

					subCells = []

					bot = [self.upperCells[i][0],midB[0],midB[2],self.upperCells[i][3]]
					top = [self.upperCells[i][4],midT[0],midT[2],self.upperCells[i][7]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[0],self.upperCells[i][1],self.upperCells[i][2],midB[2]]
					top = [midT[0],self.upperCells[i][5],self.upperCells[i][6],midT[2]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					allSubCells.extend(subCells)

				if divideX==False and divideY==False:

					allSubCells.append(self.upperCells[i])

			self.upperCells = allSubCells

			allSubCells = []

			print(maxCellX)
			print(maxCellY)

			if maxCellX<thresX and maxCellY<thresY: 
				break

		return self.upperCells


	def exportJSON(self,dir):

		crvs = self.cellMembers
		finalCells = self.finalCells
		upperCells = self.upperCells

		f = open('testCrvs.csv','w')

		for i in range(len(crvs)):

			for j in range(len(crvs[i])):

				for k in range(len(crvs[i][j])):

					f.write(str(crvs[i][j][k]))
					if k<len(crvs[i][j])-1:
						f.write(',')

				if j<len(crvs[i])-1:
					f.write(' ')

			if i<len(crvs)-1:
				f.write('\n')

		f.close()

		f = open(dir,'w')

		f.write('[')

		for i in range(len(finalCells)):

			f.write('[')

			if finalCells[i] != None:

				for j in range(len(finalCells[i])):

					f.write('[')

					for k in range(len(finalCells[i][j])):

						f.write(str(finalCells[i][j][k]))
						if k<len(finalCells[i][j])-1:
							f.write(',')
						else:
							f.write(']')

					if j<len(finalCells[i])-1:
						f.write(',')
					else:
						f.write(']')

			if i<len(finalCells)-1:
				f.write(',')
			else:
				f.write(']')

		f.close()

		########### CREATE UPPER BOXES #############

		f = open(dir.split('.')[0] + '_h.json','w')
		f.write('[')

		for i in range(len(upperCells)):
			f.write('[')
			for j in range(len(upperCells[i])):
				f.write('[')
				for k in range(len(upperCells[i][j])):
					f.write(str(upperCells[i][j][k]))
					if k<len(upperCells[i][j])-1:
						f.write(',')
					else:
						f.write(']')

				if j<len(upperCells[i])-1:
					f.write(',')
				else:
					f.write(']')

			if i<len(upperCells)-1:
				f.write(',')
			else:
				f.write(']')
		f.close()

		return [finalCells,upperCells]




class interpPipe:

	def __init__ (self, path ,radius):

		self.pts = path.pts
		self.t = path.tans
		self.sections =  []
		self.radius = radius

	def genSections (self, taper = True, gradient = .5, start = .1):

		self.radii = []

		for i in range(len(self.pts)):

			if taper and i/len(self.pts)>start:

				f = 1 - .5*m.pow(i/len(self.pts)-start,gradient)
			
			else:
			
				f = 1

			if f>1: f=1
			
			self.radii.append(self.radius*f)

			self.sections.append(genCircleSection(self.pts[i],self.t[i],radius))

	def genPipe(self, reso = 8, cap = True):

		if len(self.sections) == 0:

			self.genSections()

		self.srf = interpSrf(self.sections,len(self.pts),reso)
		self.srf.createFaces()

		if cap:
			self.srf.v[self.U[-1]]




def exportCrvsOBJ(dir, curves):

	f = open(dir,'w')

	f.write('# Rhino\n\n')

	for i in range(len(curves)):
		f.write('o object_' + str(i+1) + '\n')

		for j in range(len(curves[i])):
			f.write('v ')

			for k in range(len(curves[i][j])):

				f.write(str(curves[i][j][k]))

				if k<len(curves[i][j])-1:
					f.write(' ')

			if j < len(curves[i])-1:
				f.write('\n')

		f.write('\n')
		f.write('cstype bspline\ndeg 1\n')
		f.write('curv')

		length = np.linalg.norm(np.subtract(curves[i][j][0],curves[i][j][-1]))

		f.write(' ' + str(0))
		f.write(' ' + str(length))

		f.write(' ' + str(2*i+1) + ' ' + str(2*i+2))

		f.write('\n')
		f.write('param u')

		f.write(' ' + str(0))
		f.write(' ' + str(0))
		f.write(' ' + str(length))
		f.write(' ' + str(length))

		f.write('\n')
		f.write('end')

		if i<len(curves)-1:
			f.write('\n')

	f.close()