import math as m
import numpy as np
import sys

def exportPolyline(members,dir):

    if dir.split('.')[-1] == 'obj':

        f = open(dir,'w')

        f.write('# Rhino\n\n')

        for i in range(len(members)):
            f.write('o object_' + str(i+1) + '\n')

            for j in range(len(members[i])):
                f.write('v ')

                for k in range(len(members[i][j])):

                    f.write(str(members[i][j][k]))

                    if k<len(members[i][j])-1:
                        f.write(' ')

                if j < len(members[i])-1:
                    f.write('\n')

            f.write('\n')
            f.write('cstype bspline\ndeg 1\n')
            f.write('curv')

            length = np.linalg.norm(np.subtract(members[i][j][0],members[i][j][-1]))

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

            if i<len(members)-1:
                f.write('\n')

        f.close()

        print(dir)

    if dir.split('.')[-1] == 'csv':

        f = open(dir,'w')

        for i in range(len(members)):

            for j in range(len(members[i])):

                for k in range(len(members[i][j])):

                    f.write(str(members[i][j][k]))

                    if k<len(members[i][j])-1:
                        f.write(',')

                if j<len(members[i])-1:

                    f.write(' ')

            if i < len(members)-1:

                f.write('\n')

        f.close()

        print(dir)