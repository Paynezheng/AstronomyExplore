"""
Auther: ScutRookie
modified from Trigradient Demo in matplotlib.com

"""

import matplotlib.pyplot as plt
import numpy as np


# 势能，当然，我不知道计算公式,就假设是x^2+y^2吧(圆)
def potential(x,y):
    return x*x+y*y;

# 生成向量
x_axis = np.arange(-2.0,2.0,0.01)
# 生成矩阵
x,y=np.meshgrid(x_axis,x_axis)
#print(x_axis)
#print(y_axis)
# 这是个二维全0矩阵
ling = np.zeros(np.shape(x)) 
# for i in range(0,len[x])
#     for j in range(0,len[x[0]):
#         ling[i][j] = potential(i,j)

# 建图
plt.figure(figsize=(10,10))
# 这是填充颜色的
# plt.contourf(x,y,potential(x,y))
lines = np.arange(0.5,8.0,0.25)
plt.contour(x,y,potential(x,y),lines)

# 这个是生成矢量场的
Sparse_vector = np.arange(-2.0,2.0,0.5)
x_1,y_1=np.meshgrid(Sparse_vector,Sparse_vector)
U = 2*x_1
V = 2*y_1 
plt.quiver(x_1, y_1, U, V, potential(x_1,y_1))
plt.show()

