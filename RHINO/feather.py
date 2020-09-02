import rhinoscriptsyntax as rs
import math as m
import random as r

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

def vecToPath(st,vec,reso=10,curl=[0,0,0]):
    pts = []
    for i in range(reso):
        pt = rs.PointAdd(st,rs.VectorScale(vec,i/(reso-1)))
        pt = rs.PointAdd(pt,rs.VectorScale(curl,m.pow(i/10,2)))
        pts.append(pt)
    return pts

class Zipper:
    def __init__(self,top_barb,bot_barb,gap=.2,angle=30):
        self.top = top_barb
        self.bot = bot_barb
        self.gap = gap
        self.top_teeth = []
        self.bot_teeth = []
        self.bot.createBarbules(-angle,.01)
        self.top.createBarbules(angle,.01)
        self.reso = self.bot.barbules_r[0].reso
        self.genTeeth()
    #
    #
    def genTeeth(self):
        for i in range(len(self.top.pts)):
            vec = rs.VectorCreate(self.bot.pts[i],self.top.pts[i])
            vec = rs.VectorUnitize(vec)
            length = rs.Distance(self.bot.pts[i],self.top.pts[i])*(1-self.gap)
            top_tooth = []
            bot_tooth = []
            for j in range(self.reso):
                y = rs.VectorScale(vec,j/(self.reso-1)*length)
                top_tooth.append(rs.PointAdd(self.top.pts[i],y))
                bot_tooth.append(rs.PointAdd(self.bot.pts[i],rs.VectorScale(y,-1)))
            self.top_teeth.append(top_tooth)
            self.bot_teeth.append(bot_tooth)
    #
    #
    def modBarbules(self,factor=1):
        #
        for i in range(len(self.bot.barbules_r)):
            self.top.barbules_r[i].blendPts(self.top_teeth[i],factor)
            self.bot.barbules_l[i].blendPts(self.bot_teeth[i],factor)
        #
    #
    #


class Barb:
    def __init__(self,stem,index,length,angle,entropy,interval=20,reso=100,taper=20):
        self.stem = stem
        self.index = index
        self.ang = angle
        self.entropy = entropy
        self.interval = interval
        self.reso = reso
        self.radius = self.stem.radii[self.index]
        self.length = length
        self.bLength = length/5
        self.taper = self.length/taper
        self.x = vecRotate(self.stem.t[index],self.ang,[0,0,1])
        self.y = vecRotate(self.x,-90,[0,0,1])
        self.lock=-1
        self.ang = angle
        stPt = self.stem.pts[self.index]
        self.pts = []
        glide = rs.VectorScale(self.stem.t[self.index],self.taper)
        for i in range(self.reso):
            pt = rs.PointAdd(stPt,rs.VectorScale(self.x,self.length/(self.reso-1)*i))
            f = 1-m.pow(i/(self.reso/10),.5)
            if f>0:
                pt = rs.PointAdd(pt,rs.VectorScale(glide,-f))
            self.pts.append(pt)
        if self.entropy>0:
            self.wavePoints()
        self.t = []
        self.updateTan()
        self.barbules_r = []
        self.barbules_l = []
        self.sections = []
        self.genSections()
    #
    #
    def resetPts(self):
        stPt = self.stem.pts[self.index]
        self.pts = []
        glide = rs.VectorScale(self.stem.t[self.index],self.taper)
        for i in range(self.reso):
            pt = rs.PointAdd(stPt,rs.VectorScale(self.x,self.length/(self.reso-1)*i))
            f = 1-m.pow(i/(self.reso/10),3)
            if f>0:
                pt = rs.PointAdd(pt,rs.VectorScale(glide,-f))
            self.pts.append(pt)
        if self.entropy>0:
            self.wavePoints()
        self.t = []
        self.updateTan()
    #
    #
    def blendPts(self,goals,factor=1):
        #
        self.length = rs.Distance(goals[0],goals[-1])
        self.resetPts()
        #
        for i in range(len(goals)):
            vec = rs.VectorCreate(goals[i],self.pts[i])
            vec = rs.VectorScale(vec,factor)
            self.pts[i] = rs.PointAdd(self.pts[i],vec)
        #
        self.t = []
        self.updateTan()
        return self.pts
    #
    #
    def wavePoints(self,dir=1):
        #
        curl = self.entropy*self.stem.radii[self.index]*dir
        waveLength = self.entropy*1.5
        newPts = []
        for i in range(len(self.pts)):
            f = curl*m.sin(2*m.pi*waveLength*i/len(self.pts))
            y = rs.VectorScale(self.y,f)
            self.pts[i] = rs.PointAdd(self.pts[i],y)
    #
    #
    def updateTan(self):
        #
        self.t = []
        for i in range(len(self.pts)):
            #
            if i<len(self.pts)-1:
                t = rs.VectorCreate(self.pts[i+1],self.pts[i])
                t = rs.VectorUnitize(t)
                self.t.append(t)
            else:
                t = rs.VectorCreate(self.pts[i],self.pts[i-1])
                t = rs.VectorUnitize(t)
                self.t.append(t)
    #
    #
    def createBarbules(self,sub_entropy,angle=70,reso=20):
        #
        self.barbules_r = []
        self.barbules_l = []
        #
        # Barb(stem,index,length,angle,entropy,interval=20,reso=100)
        #
        f = 1
        for i in range(len(self.pts)):
            #
            if i>0 and i<len(self.pts)-1:
                #
                ang = angle + 2*(r.random()-.5)*30*self.entropy
                #
                if self.lock!=0:
                    self.barbules_r.append(Barb(self,i,self.bLength*f,ang,sub_entropy,reso))
                if self.lock!=1:
                    self.barbules_l.append(Barb(self,i,self.bLength*f,-ang,sub_entropy,reso))
    #
    #
    def genBarbules(self,offset,interval=4,gradient=.5):
        #
        for i in range(len(self.barbules_r)):
            inter = i+offset
            self.genSections(.1)
            if inter<len(self.barbules_r) and inter%interval==0 and i>0:
                self.barbules_r[i].genMesh(.1,.75)
                self.barbules_l[i].genMesh(.1,.75)
    #
    #
    def genSections(self,start=.1,gradient=.5,reduction=.75):
        #
        self.sections=[]
        self.radii = []
        radius = self.radius*reduction
        #
        for i in range(len(self.pts)):
            #
            if i/len(self.pts)<start:
                section_radius = radius
            else:
                f = 1-m.pow((i-start)/(len(self.pts)-start-1), gradient)+.25
                if f > 1: f = 1
                section_radius = radius*f
            #
            self.sections.append(genCircleSection(self.pts[i],self.t[i],section_radius))
            self.radii.append(section_radius)
    #
    #
    def genMesh(self,start,gradient=.5,reduction=.75):
        self.genSections(start,gradient,reduction)
        self.mesh = meshLoftSections(self.sections)


class Rachis:
    def __init__(self,path,radius,interval,reso=200):
        #
        self.radius = radius
        self.interval = interval
        self.stub_interval = 5
        self.radii = []
        #
        self.sections = []
        self.branches_r = []
        self.branches_l = []
        self.stubs_r = []
        self.stubs_l = []
        #
        self.pts = rs.DivideCurve(path,reso)
        self.t = []
        self.updateTan()
        self.zip_r = []
        self.zip_l = []
        self.genSections()
    #
    #
    def updateTan(self):
        #
        self.t = []
        for i in range(len(self.pts)):
            #
            if i<len(self.pts)-1:
                t = rs.VectorCreate(self.pts[i+1],self.pts[i])
                t = rs.VectorUnitize(t)
                self.t.append(t)
            else:
                t = rs.VectorCreate(self.pts[i],self.pts[i-1])
                t = rs.VectorUnitize(t)
                self.t.append(t)
    #
    #
    def createBarbs(self,start,bLength,sLength,entropy,gradient=1):
        #
        self.stubs_r = []
        self.stubs_l = []
        self.branches_r = []
        self.branches_l = []
        #
        # Barb(stem,index,length,angle,entropy,reso=100)
        #
        f = 1
        gap = .2
        start = int(len(self.pts)*start)
        #
        for i in range(len(self.pts)):
            #
            if i>start and i%self.stub_interval == 0 and i%self.interval != 0:
                self.stubs_r.append(Barb(self,i,sLength,45,0,reso=10))
                self.stubs_l.append(Barb(self,i,sLength,-45,0,reso=10))
                #
            if i>start and i%self.interval == 0 and i!=len(self.pts)-1:
                f = 1-m.pow((i-start)/(len(self.pts)-start-1), gradient)+.4
                if f > 1: f = 1
                if f*bLength>sLength:
                    self.branches_r.append(Barb(self,i,bLength*f,50,entropy))
                    self.branches_l.append(Barb(self,i,bLength*f,-50,entropy))
        for i in range(len(self.branches_r)-1):
            self.zip_r.append(Zipper(self.branches_r[i],self.branches_r[i+1],gap,-30))
            self.zip_l.append(Zipper(self.branches_l[i],self.branches_l[i+1],gap,30))
    #
    #
    def genBarbs(self,gradient=.5):
        #
        for i in range(len(self.branches_r)):
            self.branches_r[i].genMesh(0,.5)
            self.branches_l[i].genMesh(0,.5)
        #
        for i in range(len(self.stubs_l)):
            self.stubs_r[i].genMesh(.1,1)
            self.stubs_l[i].genMesh(.1,1)
    #
    #
    def genSections(self,start=.1,gradient=.5):
        #
        self.sections=[]
        self.radii = []
        radius = self.radius
        start = int(start*len(self.pts))
        #
        for i in range(len(self.pts)):
            #
            if i<start:
                section_radius = radius
            else:
                f = 1-m.pow((i-start)/(len(self.pts)-start-1), gradient)+.25
                if f > 1: f = 1
                section_radius = radius*f
            #
            self.sections.append(genCircleSection(self.pts[i],self.t[i],section_radius))
            self.radii.append(section_radius)
    #
    #
    def genMesh(self,start,gradient=.5):
        #
        self.genSections(start,gradient)
        self.mesh = meshLoftSections(self.sections)
    #
    #


def Main():
    path = rs.GetObject("please select feather path",rs.filter.curve)
    rs.EnableRedraw(False)
    radius = 3
    barb_interval = 10
    start = .1
    #
    barb_len = 60
    stub_len = 10
    grade_barb = .75
    entropy_barb = .5
    grade_thick = .15
    #
    barbule_entropy = .5
    #
    rachis = Rachis(path,radius,barb_interval)
    rachis.createBarbs(start,barb_len,stub_len,entropy_barb,grade_barb)
    rachis.genMesh(.1,grade_thick)
    rachis.genBarbs()
    #
    interlock = len(rachis.branches_r)/10
    if interlock<.5:
        interlock = 0
    if interlock>.5:
        interlock = 1
    #
    for i in range(len(rachis.branches_r)):
        #
        barb_r = rachis.branches_r[i]
        barb_l = rachis.branches_l[i]
        #
        if interlock>.75:
            if i<len(rachis.branches_r)-1:
                if i==0:
                    rachis.zip_r[i].bot.lock=0
                    rachis.zip_l[i].bot.lock=0
                rachis.zip_r[i].modBarbules(interlock)
                rachis.zip_l[i].modBarbules(interlock)
            interval = 4
        else:
            rachis.branches_r[i].createBarbules(barbule_entropy)
            rachis.branches_l[i].createBarbules(barbule_entropy)
            interval = 40
        #
        rachis.branches_r[i].genBarbules(2*(i%2),interval)
        rachis.branches_l[i].genBarbules(2*(i%2),interval)

Main()