from math import sqrt

def mag(a):
    return sqrt(a[0]**2+a[1]**2 +a[2]**2)

def idx(t, nc):
    return t[2]*nc*nc + t[1]*nc + t[0]

def no_neg(alpha, rcut):
    results = []
    nc = 10/(rcut/alpha)
    for x in range(1, alpha+2):
        for y in range(1, alpha+2):
            for z in range(1, alpha+2):
                
                if (x-1)**2+(y-1)**2+(z-1)**2 < alpha**2:
                    results.append((x,y,z))

    for i in range(1, alpha+1):
        results.append((0,0,i))
        results.append((0,i,0))
        results.append((i,0,0))

    print(results)
    results = [idx(_,nc) for _ in results]
    results += [0.0]
    results.sort()

    return results

def second_way(alpha, rcut):
    results = []
    nc = 10/(rcut/alpha)
    for x in range(1, alpha+2):
        for y in range(1, alpha+2):
            for z in range(1, alpha+2):
                
                if (x-1)**2+(y-1)**2+(z-1)**2 < alpha**2:
                    results.append((x,y,z))
                    results.append((x,y,-z))
                    results.append((x,-y,z))
                    results.append((x,-y,-z))
                    results.append((-x,y,z))
                    results.append((-x,y,-z))
                    results.append((-x,-y,z))
                    results.append((-x,-y,-z))


    for i in range(1, alpha+1):
        results.append((0,0,i))
        results.append((0,i,0))
        results.append((i,0,0))
        results.append((0,0,-i))
        results.append((0,-i,0))
        results.append((-i,0,0))

    print(results)
    results = [idx(_,nc) for _ in results]
    results += [0.0]
    results.sort()

    return results

second_way(1,1)
no_neg(1,2)

for r in range(1,10):
    for a in range(1,10):
        nm = (2*a+1)**3
        #re = second_way(a,r)
        #if len(re) > nm:
        #    print(f'rcut: {r}, a: {a}, NM: {nm}, RE:{len(re)}')
