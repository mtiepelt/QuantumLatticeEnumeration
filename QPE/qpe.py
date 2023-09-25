from sage.all import exp
import numpy as np

p1 = .25
p2 = .5

def p1_bound(a, t):
    assert t > 0
    return (1 + p1*(exp(t) - 1))/exp(a*t)

def p2_bound(a, t):
    assert t < 0
    # return (1 + p2*(exp(t) - 1))/exp(a*t)
    return (exp(t) + p2*(1 - exp(t)))/exp(a*t)


from sage.all import oo
def find_min(f, ran):
    min = oo
    arg = None
    for x in ran:
        val = f(x)
        if val < min:
            min = val
            arg = x
    return min, arg


for a in np.arange(0.368, 0.370, 0.000001): # after looking at larger but coarser ranges, and zoomin in
    min_p1, arg_p1 = find_min(lambda t: p1_bound(a, t), np.arange(0.01, 1., 0.01))
    min_p2, arg_p2 = find_min(lambda t: p2_bound(a, -t), np.arange(0.01, 1., 0.01))
    # min_p1, arg_p1 = find_min(lambda t: p1_bound(a, t), [.5])
    # min_p2, arg_p2 = find_min(lambda t: p2_bound(a, -t), [.5])

    for eps in range(18, 21):
        if min_p1**eps < 1/2 and min_p2**eps < 1/2:
            print(f"alpha = {round(a,8)}, epsilon = {eps}, minimising p1 ({min_p1**eps}) at {arg_p1} and p2 ({min_p2**eps}) at {arg_p2}")
            break

