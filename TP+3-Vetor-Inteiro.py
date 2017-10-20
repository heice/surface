
# coding: utf-8

# # Résolution de problèmes multiobjectif par agrégation pondérée
#
# ## Utilisation de Scipy
#
# La fonction *minimize* fournie dans Scipy a un fonctionnement similaire à fmincon dans Matlab:

# In[1]:

from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import math
from pymopt import *
# Pour avoir les figures intégrées dans la page html:
get_ipython().magic(u'matplotlib inline')

# In[2]:

#Constants
Xmax = 20
Ymax = 12
PR = [0, 6.5]
PE = [0, 1.5]
lp = [0, 4]
l = [4, 4, 3, 2, 2, 2, 2, 2, 3, 2, 3, 3]

# In[3]:

#Fonctions
def dist(x1, y1, x2, y2):
    return sqrt((x1-x2)**2+(y1-y2)**2)

def croise2(x1, y1, lx1, ly1, x2, y2, lx2, ly2):
    #exception
    if (x1+lx1)<(x2+lx2) and x1>x2 and y1>y2 and (y1+ly1)<(y2+ly2) :
        return lx1*ly1
    else :
        if (x1+lx1)>(x2+lx2) and x1<x2 and y1<y2 and (y1+ly1)>(y2+ly2) :
            return lx2*ly2
        else :
                if (x1 <= x2) and (x1+ lx1) > x2 :
                    if ( (y1 >= y2) and y1 < (y2 + ly2)) :
                        return ((y2 + ly2) - y1)*((x1+ lx1) - x2)
                    else :
                        if (y1 < y2) and (y1 + ly1) > y2 :
                            return ((y1 + ly1) - y2)*((x1+ lx1) - x2)
                        else :
                            return 0

                else :
                    if (x1 > x2) and x1 < (x2 + lx2) :
                        if ( (y1 > y2) and y1 < (y2 + ly2)) :
                            return ((y2 + ly2) - y1)*((x2 + lx2) - x1)
                        else :
                            if (y1 < y2) and (y1 + ly1) > y2 :
                                return ((y1 + ly1) - y2)*((x2 + lx2) - x1)
                            else :
                                return 0
                    else :
                        return 0

# In[4]:

# Définition de la fonction objectif:
def f1(x):
    somme = 0
    for i in range(2, 5):
        somme = somme + dist(x[0],x[1],x[i*2],x[i*2+1]) + dist(x[2],x[3],x[i*2],x[i*2+1])
    somme = somme + dist(x[4],x[5],x[6],x[7]) + dist(x[8],x[9],x[6],x[7])
    return somme

def f2(x):
    somme = 0
    for i in range(6):
        for j in range(i+1,6):
            somme = somme + dist(x[i*2],x[i*2+1],x[j*2],x[j*2+1])
    return -somme

# In[5]:

# Bornes inférieures
lb=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
# Bornes supérieures
ub=[16, 8, 17, 10, 18, 10, 18, 10, 18, 10, 17, 9]
# Contraintes g(x) <= 0
C=[lambda x: -x[0]+l[0]+3,   lambda x: -x[2]+l[2]+3,   lambda x: -x[4]+l[4]+3,   lambda x: -x[6]+l[6]+3,   lambda x: -x[8]+l[8]+3,  lambda x: -x[10]+l[10]+3]

C += [lambda x: croise2(x[0],x[1],l[0],l[1],x[2],x[3],l[2],l[3])]
C += [lambda x: croise2(x[0],x[1],l[0],l[1],x[4],x[5],l[4],l[5])]
C += [lambda x: croise2(x[0],x[1],l[0],l[1],x[6],x[7],l[6],l[7])]
C += [lambda x: croise2(x[0],x[1],l[0],l[1],x[8],x[9],l[8],l[9])]
C += [lambda x: croise2(x[0],x[1],l[0],l[1],x[10],x[11],l[10],l[11])]

C += [lambda x: croise2(x[2],x[3],l[2],l[3],x[4],x[5],l[4],l[5])]
C += [lambda x: croise2(x[2],x[3],l[2],l[3],x[6],x[7],l[6],l[7])]
C += [lambda x: croise2(x[2],x[3],l[2],l[3],x[8],x[9],l[8],l[9])]
C += [lambda x: croise2(x[2],x[3],l[2],l[3],x[10],x[11],l[10],l[11])]

C += [lambda x: croise2(x[4],x[5],l[4],l[5],x[6],x[7],l[6],l[7])]
C += [lambda x: croise2(x[4],x[5],l[4],l[5],x[8],x[9],l[8],l[9])]
C += [lambda x: croise2(x[4],x[5],l[4],l[5],x[10],x[11],l[10],l[11])]

C += [lambda x: croise2(x[6],x[7],l[6],l[7],x[8],x[9],l[8],l[9])]
C += [lambda x: croise2(x[6],x[7],l[6],l[7],x[10],x[11],l[10],l[11])]

C += [lambda x: croise2(x[8],x[9],l[8],l[9],x[10],x[11],l[10],l[11])]

# Instanciation du problème
pb = Problem([f1,f2],C,lb,ub)

# In[6]:

# Définition du nombre d'individus dans la population
pb_size=100
# Définition du nombre de générations à calculer
ngen=100
# Définition des probabilités de croissance et mutation
pcross=0.99
pmut=0.01

# In[7]:

# Résolution du problème et affichage du minimum trouvé
sols=nsga2(pb,pb_size, ngen, pcross, pmut)

# Affichage de la dernière population
fig_y2(sols)

# In[16]:

s = sols[99]
print(s)
rects=[]
facecolor='r'
edgecolor='None'
alpha=0.5

for i in range(int(len(s)/2)):
    x_i=s[2*i]
    y_i=s[2*i+1]
    l_i = l[2*i]
    L_i = l[2*i+1]
    rects.append(Rectangle((x_i,y_i),l_i,L_i))
    print(x_i,y_i,l_i,L_i)

pc = PatchCollection(rects, facecolor=facecolor, alpha=alpha, edgecolor=edgecolor)
fig, ax = plt.subplots(1)
ax.add_collection(pc)
ax.set_xlim([0,20])
ax.set_ylim([0,12])
plt.show()

patterns = ['-', '+', 'x', 'o', 'O', '.']  # more patterns
fig4 = plt.figure()
ax4 = fig4.add_subplot(111, aspect='equal')
for p in [
    patches.Rectangle(
        ((s)[i*2],(s)[i*2+1]),
        l[i*2],
        l[i*2+1],
        hatch=patterns[i],
        fill=False
    ) for i in range(len(patterns))
]:
    ax4.add_patch(p)
ax4.set_xlim([0,20])
ax4.set_ylim([0,12])
