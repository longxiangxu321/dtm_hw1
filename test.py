import math
def get_area(*args):
    t_area = 0
    for i in range(len(args)-1):
        t_area = t_area + (args[i][0] * args[i + 1][1]) - (args[i+1][0] * args[i][1] )

    t_area = t_area + args[-1][0] * args[0][1] - args[0][0] * args[-1][1]
    t_area = abs(t_area) / 2
    return t_area

def tri(p0,p1,p2):
    a = math.dist(p0,p1)
    b = math.dist(p0,p2)
    c = math.dist(p1,p2)
    s = (a + b + c) / 2
    area = (s*(s-a)*(s-b)*(s-c)) ** 0.5
    return area

print(get_area([0, 4], [2, 0], [0, -3]))
