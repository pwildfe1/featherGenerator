import rhinoscriptsyntax as rs
import math as m
import sys

def vecRotate(vec,ang,axis):
    cos = m.cos(m.pi/180*ang)
    sin = m.sin(m.pi/180*ang)
    v = vec
    u = rs.VectorUnitize(axis)
    #u = [axis[0]/np.linalg.norm(axis),axis[1]/np.linalg.norm(axis),axis[2]/np.linalg.norm(axis)]
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
# genCircleSection: creates a circle of points around a single point on a line
# 
# cnt(3D vector) = point on the line
# vector(3D vector) = tangent at point on line
# radius(float) = the radius of the circle formed
# axis(default=[0,0,1]) = used to find the plane at that point on the line
# reso(default=8) = the number of points defining the circle
####

def genCircleSection(cnt,vector,radius,axis=[0,0,1],reso=8):
    pts = []
    vec = rs.VectorUnitize(vector)
    v = vecRotate(vec,90,axis)
    for i in range(reso):
        mv = vecRotate(v,(360/reso)*i,vector)
        mv = rs.VectorScale(mv,radius)
        pts.append(rs.PointAdd(cnt,mv))
    pts.append(pts[0])
    return pts


####
# meshLoftSections: creates a quad based mesh by lofting between equal groups of points (sections)
# 
# sections (array) = an array of points to be lofted between
####

def meshLoftSections(sections,closed=True):
    faces = []
    v = []
    for i in range(len(sections)):
        v.extend(sections[i])
    for i in range(len(sections)-1):
        for j in range(len(sections[0])-1):
            index0 = i*len(sections[0]) + j
            index1 = i*len(sections[0]) + (j + 1)%len(sections)
            index2 = (i + 1)*len(sections[0]) + (j + 1)%len(sections)
            index3 = (i + 1)*len(sections[0]) + j
            faces.append([index0,index1,index2,index3])
    if closed:
        st_cnt, en_cnt = [0,0,0], [0,0,0]
        for i in range(len(sections[0])):
            st_cnt = rs.PointAdd(st_cnt,sections[0][i])
            en_cnt = rs.PointAdd(en_cnt,sections[-1][i])
        st_cnt = rs.VectorScale(st_cnt,1/len(sections[0]))
        en_cnt = rs.VectorScale(en_cnt,1/len(sections[0]))
        init = 0
        v.append(st_cnt)
        for j in range(len(sections[0])-1):
            index0 = init + j
            index1 = len(v)-1
            index2 = init + (j+1)%(len(sections[0])-1)
            index3 = init + j
            faces.append([index0,index1,index2])
        #init = len(sections)*(len(sections[0])-1) + 1
        init = len(v)-1-len(sections[0])
        v.append(en_cnt)
        for j in range(len(sections[0])-1):
            index0 = init + j
            index1 = len(v)-1
            index2 = init + j + 1
            index3 = init + j
            faces.append([index0,index1,index2])
    mesh = rs.AddMesh(v,faces)


####
# importPaths: imports path information from the feather_v2.py file 
# 
# loc(str) = file location
####
1
def importPaths(loc):
    f = open(loc,'r')
    paths = f.read().split('\n')
    crvs = []
    profiles = []
    for i in range(len(paths)):
        entries = paths[i].split(' ')
        pts = []
        radii = []
        for j in range(len(entries)):
            # EACH ENTRY HAS 4 PARTS: X,Y,Z and Radius
            entry = entries[j].split(',')
            point = []
            point.append(float(entry[0]))
            point.append(float(entry[1]))
            point.append(float(entry[2]))
            pts.append(point)
            if len(entry)>3:
                radii.append(float(entry[3]))
        crvs.append(pts)
        profiles.append(radii)
    return [crvs,profiles]


def Main(mesh=True):
    info = importPaths('all_paths.csv')
    crvs = info[0]
    radii = info[1]
    rs.EnableRedraw(False)
    for i in range(len(crvs)):
        sections = []
        for j in range(len(crvs[i])):
            if j<len(crvs[i])-1:
                t = rs.VectorCreate(crvs[i][j+1],crvs[i][j])
            else:
                t = rs.VectorCreate(crvs[i][-1],crvs[i][-2])
            if mesh:
                sections.append(genCircleSection(crvs[i][j],t,radii[i][j]))
        if mesh:
            meshLoftSections(sections)
        else:
            rs.AddCurve(crvs[i])

Main(True)