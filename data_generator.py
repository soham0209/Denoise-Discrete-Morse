import numpy as np
import os
import sys
import vtk
from tqdm import tqdm


class bbox:
    def __init__(self):
        self.xmax = float("-inf")
        self.xmin = float("inf")
        self.ymax = float("-inf")
        self.ymin = float("inf")
        self.zmax = float("-inf")
        self.zmin = float("inf")

    def getWidth(self):
        width = [self.xmax-self.xmin, self.ymax-self.ymin, self.zmax-self.zmin]
        return width

    def getbox(self, vertices):
        for vert in vertices:
            if vert[0] > self.xmax:
                self.xmax = vert[0]
            if vert[0] < self.xmin:
                self.xmin = vert[0]
            if vert[1] > self.ymax:
                self.ymax = vert[1]
            if vert[1] < self.ymin:
                self.ymin = vert[1]
            if vert[2] > self.zmax:
                self.zmax = vert[2]
            if vert[2] < self.zmin:
                self.zmin = vert[2]

    def __str__(self):
        return 'Xmin %s\nXmax %s\n Ymin %s\n Ymax %s\n Zmin %s\n Zmax %s' % (self.xmin, self.xmax, self.ymin, self.ymax,
                                                                             self.zmin, self.zmax)


def gaussian(X, var=1.0):
    v = np.sum(np.multiply(X, X) / var, axis=1)
    v = np.exp(-v)
    val = np.sum(v)
    return val

def writedata3Dv2(pts, mn_pts,  vertexfile, std_dv=0.1):
    fieldvalues = list()
    print('Total vertices read ', len(pts))
    std_dev = std_dv
    var = std_dev ** 2
    for p in tqdm(pts):
        v = mn_pts - p
        val = gaussian(v, var)
        fieldvalues.append(val)
        vertexfile.write(str(p[0]) + ' ' + str(p[1]) + ' ' + str(p[2]) + '\n')
    vertexfile.close()
    return fieldvalues


def createfromFile(fname, bounds: bbox):
    filepts = np.loadtxt(fname)
    return filepts


cwd = os.getcwd()
dataName = ''
sigma_dict = {'bimba':0.0151, 'botijo':0.0227, 'MotherChild':12.2478}
if len(sys.argv) > 1:
    dataName = sys.argv[1]
if dataName not in sigma_dict:
    raise Exception('Dataset does not exist. Please try between bimba/botijo/MotherChild')
std = 2.0 *  sigma_dict[dataName]
if len(sys.argv) > 2:
    std = float(sys.argv[2])


dirname = cwd + '/' + dataName

if not os.path.exists(dirname):
    os.mkdir(dirname)


datafilename = dirname + '/' + dataName + '.txt'
vertfilename = dirname + '/' + dataName + '_vert.txt'
bbox1 = bbox()

field_to_be_calculated = dirname + '/' + dataName + '_vert_refined.txt'
mean_pts = createfromFile(dataName + '/' + dataName + '_vert_Q.txt', bbox1)


pts = np.loadtxt(field_to_be_calculated)
vertfile = open(vertfilename, 'w')


values = writedata3Dv2(np.array(pts), mean_pts, vertfile, std_dv=std)
datafile = open(datafilename, 'w')
datafile.write('3\n' + str(len(values)) + '\n1\n1\n')
values = np.array(values) / np.max(values)
np.savetxt(datafile, values, fmt='%.6f')
datafile.close()
