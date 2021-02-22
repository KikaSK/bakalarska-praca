# -*- coding: utf-8 -*-
"""Vykreslovanie.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1odaj3o6dFoWspupxTegM8RmvqxL4Svp1
"""

! pip3 install pygnuplot
! pip3 install numpy
! pip3 install PyWavefront

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math

fig = plt.figure()
ax = fig.gca(projection='3d')

# Dáta pre funkciu
precision = 0.001
X = np.arange(-2, 2, precision)
Y = np.arange(-2, 2, precision)
X, Y = np.meshgrid(X, Y)

# Funkcia
Z = X**2+Y**2

class Triangle:
  def __init__(self, point_a, point_b, point_c):
    self.a = point_a
    self.b = point_b
    self.c = point_c

# Vykreslenie funkcie
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

def is_triangle(T):
  A = T.a
  B = T.b
  C = T.c

  distAB = math.sqrt( (A[0]-B[0])**2 + (A[1]-B[1])**2 + (A[2]-B[2])**2 )
  distCB = math.sqrt( (C[0]-B[0])**2 + (C[1]-B[1])**2 + (C[2]-B[2])**2 )
  distAC = math.sqrt( (A[0]-C[0])**2 + (A[1]-C[1])**2 + (A[2]-C[2])**2 )
  
  if( (distAB + distCB > distAC) and (distAB + distAC > distCB) and (distCB + distAC > distAB)):
    return True
  return False

# Rozsah osi z
ax.set_zlim(-10, 10)

# Počet deliacich bodov na osi z
ax.zaxis.set_major_locator(LinearLocator(10)) 

# Počet desatinných miest deliacich bodov na osi z
ax.zaxis.set_major_formatter(FormatStrFormatter('%.2f'))


# Veľkosť grafu
fig.set_size_inches(18.5, 10.5)

# Na boku grafu, asociácia farby grafu s funkčnou hodnotou
fig.colorbar(surf, shrink=0.5, aspect=5)

# Pole trojuholníkov
# pr. triangles = [Triangle((1, 1, 1), (2, 2, 2), (1, 3, 4)), Triangle((2, 3, 4), (9, 9, 9), (3, 4, 5))]
my_triangles = []

# Načítanie súradníc trojuholníkov zo súboru troj.txt
with open("troj.txt", mode='r') as f:
  for l in f.readlines():
    val = [int(i) for i in l.split()]
    if len(val) != 9:
      continue
    
    T = Triangle((val[0], val[1], val[2]), (val[3], val[4], val[5]), (val[6], val[7], val[8]));
    
    # skontrolujeme, či zadané trojuholníky spĺňajú trojuholníkovú nerovnosť
    assert is_triangle(T), "Input points don't correspond to triangle."
    
    print(val)
    my_triangles.append(T)

# Vytvorenie polí x,y,z
# Tieto polia berie ako parameter funkcia plot_trisurf()
# V poli x sa nachádzajú všetky x-ové súradnice, v poli y všetky y-ové, atď.
# Napr. pre dva trojuholníky dané vrcholmi 
# T1 - (0, 1, 4), (5, 4, 7), (2, 3, 1)
# T2 - (3, 4, 1), (6, 5, 7), (2, 1, 9)
# budú polia vyzerať nasledovne: x = (0, 5, 2, 3, 6, 2), y = (1, 4, 3, 4, 5, 1), z = (4, 7, 1, 1, 7, 9)
def create_input(triangles):
  x, y, z = [], [], []
  for t in triangles:
    x += [t.a[0], t.b[0], t.c[0]]
    y += [t.a[1], t.b[1], t.c[1]]
    z += [t.a[2], t.b[2], t.c[2]]
  return x, y, z

x, y, z = create_input(my_triangles)

tri_indexes = [(3 * i, 3 * i + 1, 3 * i + 2) for i in range(len(x)//3)]
ax.plot_trisurf(x, y, z, triangles = tri_indexes)

plt.show()


