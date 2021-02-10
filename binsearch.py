def f(x):
    return x-0.0003
def binsearch(s,l,tolerance):
    mid = (s+l)/2
    i = 1

    #while i<100000:
    while abs(mid**i-mid**(i-1))>tolerance:
        if f(mid)==0:
            return mid
        elif f(s)*f(mid) >0:
            s=mid
        else:
            l=mid
        i+=1
        mid=(s+l)/2
    return mid

print(binsearch(100,-23,0.0000000000000000000000001)) 