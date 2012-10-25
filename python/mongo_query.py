#! /usr/bin/python

# to import: mongoimport -d test -c miniFOF --type csv --file miniFOF.csv --headerline

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pymongo import Connection

connection = Connection('localhost', 27017)
db = connection['test']
collection = db['bigFOF']

result = collection.find({ "snapnum" : 85 }, { 'x' : 1, 'y' : 1, 'z' : 1 })
rows = result.count()

xs = numpy.zeros(rows)
ys = numpy.zeros(rows)
zs = numpy.zeros(rows)

dump = open('minidump.csv', 'w')

for i, row in enumerate(result):
	xs[i] = row['x']
	ys[i] = row['y']
	zs[i] = row['z']
	dump.write(str(xs[i]) + ',' + str(ys[i]) + ',' + str(zs[i]) + '\n')

dump.close()


#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(xs, ys, zs, c='r', marker='o')

#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Z')

#plt.show()
connection.close()
