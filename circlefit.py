from scipy import optimize
from numpy import *

def leastsquares(x,y):

    def calc_R(c):
        """ calculate the distance of each 2D points from the center c=(xc, yc) """
        return sqrt((x-c[0])**2 + (y-c[1])**2)

    def f_2(c):
        """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
        Ri = calc_R(c)
        return Ri - Ri.mean()



    x_mean = mean(x)
    y_mean = mean(y)

    x = x - x_mean
    y = y - y_mean

    max_x = amax(x)
    max_y = amax(y)

    x_m = 0
    y_m = 0

    x = x/max_x
    y = y/max_y
    
    center_estimate = x_m, y_m
    center_2, ier = optimize.leastsq(f_2, center_estimate)

    xc_2, yc_2 = center_2
    Ri_2       = calc_R(center_2)
    R_2        = Ri_2.mean()
    residu_2   = sum((Ri_2 - R_2)**2)
    residu2_2  = sum((Ri_2**2-R_2**2)**2)

    return R_2, xc_2, yc_2, x_mean, y_mean, max_x, max_y, residu_2
