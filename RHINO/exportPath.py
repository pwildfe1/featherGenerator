import rhinoscriptsyntax as rs
import math as m
import sys


def extractCrvs(crvs,reso,dir):
    f = open(dir,'w')
    for i in range(len(crvs)):
        divPts = rs.DivideCurve(crvs[i],reso)
        for j in range(len(divPts)):
            pt = divPts[j]
            f.write(str(pt[0]))
            f.write(',')
            f.write(str(pt[1]))
            f.write(',')
            f.write(str(pt[2]))
            if j < len(divPts)-1:
                f.write(' ')
        if i<len(crvs)-1:
            f.write('\n')
    f.close()


def Main():
    paths = rs.GetObjects("please select path curves",rs.filter.curve)
    if paths != None:
        extractCrvs(paths,200,'paths.csv')

Main()