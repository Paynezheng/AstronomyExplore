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

#----------------------------------------------------------------------
# 生成坐标系
#---------------------------------------------------------------------

# 生成向量(横坐标范围)
x_axis = np.arange(-2.0,2.0,0.1)
y_axis = np.arange(-2.0,2.0,0.1)
# 生成矩阵(两个都是二维矩阵，表示一个点的横纵坐标)
x_col,y_row=np.meshgrid(x_axis,y_axis)


#---------------------------------------------------------------------
# 计算势能
#----------------------------------------------------------------------
# omage=1
def potential(x,y):

    if y == 0 and x == -d*star2/(star2+star1):
        return 0
    if y == 0 and x == d*star1/(star2+star1):
        return 0
    
    # |r-r1|
    R_r1 = mt.sqrt((x+d*star2/(star2+star1))*(x+d*star2/(star2+star1))+y*y)
    # |r-r2|
    R_r2 = mt.sqrt((x-d*star1/(star2+star1))*(x-d*star1/(star2+star1))+y*y)
    # omega的二次方
    omega_2 = (star1+star2)/(d*d*d)
    return -star1/R_r1-star2/R_r2-omega_2*(x*x+y*y)/2

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
# lines = np.arange(-10294370.15547379,-11000000.15547379,5000)
plt.contour(x_col,y_row,ling)#,lines)


#---------------------------------------------------------------------
# 绘制矢量场(箭头)，需要计算加速度
#----------------------------------------------------------------------

# #依然按照这样计算
# ling = np.zeros(np.shape(x_col)) 
# for i in range(0,len(x_col)):
#     for j in range(0,len(x_col[0])):
#         x_temp = x_col[i][j]
#         y_temp = y_row[i][j]
#         # 给二维的全0矩阵赋值
#         ling[i][j]=potential(x_temp,y_temp)


# #这个是生成矢量场的函数调用
# Sparse_vector = np.arange(-2.0,2.0,0.5) #这个矩阵关系到矢量的密度
# x_1,y_1=np.meshgrid(Sparse_vector,Sparse_vector)
# #矢量的x,y
# U = 2*x_1
# V = 2*y_1 
# plt.quiver(x_1, y_1, U, V, potential(x_1,y_1))
plt.show()

