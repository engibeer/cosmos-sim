#! /usr/bin/python

import psycopg2
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

conn = psycopg2.connect("dbname=dark user=postgres password=lesbians")
cur = conn.cursor()

query = """
select x,y,z from "miniMDR1"."FOF" where snapnum = 85
"""

"""
with massive_halo as (
   select x,y,z from "miniMDR1"."FOF" where snapnum=85 order by np desc limit 1
  )
   select f.x, f.y, f.z, f.snapnum from massive_halo mh, "miniMDR1"."FOF" f
   where f.snapnum = 85
   and f.x between (mh.x - 500) and (mh.x + 500)
   and f.y between (mh.y - 500) and (mh.y + 500)
   and f.z between (mh.z - 500) and (mh.z + 500)
"""

cur.execute(query)
arr = numpy.array(cur.fetchall())

xs = arr[:,0]
ys = arr[:,1]
zs = arr[:,2]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs, ys, zs, c='r', marker='o')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()

cur.close()
conn.close()
