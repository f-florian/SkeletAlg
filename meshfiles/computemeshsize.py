#!python3
import math

triangleArea = 0.25
totalArea = 24
ratio = totalArea/triangleArea

def size2Triangles(meshsize):
    return ratio / meshsize**2

def triangles2Size(triangles):
    return math.sqrt(ratio / triangles)

# sizes =[0.09  , 0.0450, 0.0337, 0.0225, 0.0168, 0.0112, 0.0084, 0.0056, 0.0042]
# print(list(map(size2Triangles, sizes)))

triangles = [
    100,
    300,
    700,
    1e3,
    2e3,
    3e3,
    4e3,
    5e3,
    7e3,
    1e4,
    2e4,
    3e4,
    5e4,
    7e4,
    1e5,
    3e5,
    5e5,
    7e5,
    1e6,
    2e6,
    3e6] 
sizes = (list(map(triangles2Size, triangles)))
for j in sizes:
    print(j)

