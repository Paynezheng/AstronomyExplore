"""
Auther: ScutRookie

"""

import matplotlib.pyplot as plt
import numpy as np
import math as mt



#-------------------------------------------------------------------
# 输入 参数
# 将star1 star2放在(0,0)两侧
#--------------------------------------------------------------------
star1 = float(input("Enter Star1(>0):"))
star2 = float(input("Enter Star2(>0):"))  
d = float(input("Enter d(>0):"))

m = star1+star2
x1 = -star2*d/m
x2 = star1*d/m
#----------------------------------------------------------------------
# 生成坐标系
#---------------------------------------------------------------------

# 生成向量(横坐标范围)
x_axis = np.arange(-5.0,5.0,0.01)
y_axis = np.arange(-5.0,5.0,0.01)
# 生成矩阵(两个都是二维矩阵，表示一个点的横纵坐标)
x_col,y_row=np.meshgrid(x_axis,y_axis)


#---------------------------------------------------------------------
# 计算势能
#----------------------------------------------------------------------
# omage=1
def potential(x,y):
    if y == 0 and x == x1:
        return 0
    if y == 0 and x == x2:
        return 0
    return -star1/(((x-x1)**2+y**2)**0.5)-star2/(((x-x2)**2+y**2)**0.5)-m*(x**2+y**2)/(d**3)/2

## FIXME: 这里的势能计算可能有问题(公式不对？)

# 这是个二维全0矩阵
ling = np.zeros(np.shape(x_col)) 
for i in range(0,len(x_col)):
    for j in range(0,len(x_col[0])):
        x_temp = x_col[i][j]
        y_temp = y_row[i][j]
        # 给二维的全0矩阵赋值
        ling[i][j]=potential(x_temp,y_temp)

#print(ling)

#---------------------------------------------------------------------
# 建图，绘制等值线
#----------------------------------------------------------------------

plt.figure(figsize=(10,10))
# lines用来指定划线的位置，可以不指定
lines = np.arange(-10,0,0.1)
plt.contour(x_col,y_row,ling,lines)


#---------------------------------------------------------------------
# 绘制矢量场(箭头)，需要计算加速度
#----------------------------------------------------------------------


Sparse_vector = np.arange(-5.0,5.0,0.25) #这个矩阵关系到矢量的密度
x_1,y_1=np.meshgrid(Sparse_vector,Sparse_vector)
def accelerationU(x,y):
    if y == 0 and x == -d*star2/(star2+star1):
        return 0
    if y == 0 and x == d*star1/(star2+star1):
        return 0
    return -star1/(((x-x1)**2+y**2)**1.5)*(x-x1)-star2/(((x-x2)**2+y**2)**1.5)*(x-x2)+m*x/(d**3)

def accelerationV(x,y):
    if y == 0 and x == -d*star2/(star2+star1):
        return 0
    if y == 0 and x == d*star1/(star2+star1):
        return 0
    return -star1/(((x-x1)**2+y**2)**1.5)*y-star2/(((x-x2)**2+y**2)**1.5)*y+m*x/(d**3)


U = np.zeros(np.shape(x_1)) 
V = np.zeros(np.shape(x_1))
for i in range(0,len(x_1)):
    for j in range(0,len(x_1[0])):
        x_temp = x_1[i][j]
        y_temp = y_1[i][j]
        # 给二维的全0矩阵赋值
        U[i][j]=accelerationU(x_temp,y_temp)
        V[i][j]=accelerationV(x_temp,y_temp)


# #这个是生成矢量场的函数调用

################################################
# x_1:
# y_1:
# U:
# V:
################################################
plt.quiver(x_1, y_1, U, V)#, potential(x_1,y_1))
plt.show()

