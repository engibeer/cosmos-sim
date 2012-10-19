#! /usr/bin/python

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pymongo import Connection

connection = Connection('localhost', 27017)
db = connection['test']
collection = db['miniFOF']

result = collection.find({ "snapnum" : 85 }, { 'x' : 1, 'y' : 1, 'z' : 1 })
rows = result.count()

xs = numpy.zeros(rows)
ys = numpy.zeros(rows)
zs = numpy.zeros(rows)

for i, row in enumerate(result):
	xs[i] = row['x']
	ys[i] = row['y']
	zs[i] = row['z']

print max(xs)
print min(xs)
print
print max(ys)
print min(ys)
print
print max(zs)
print min(zs)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs, ys, zs, c='r', marker='o')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()
connection.close()
