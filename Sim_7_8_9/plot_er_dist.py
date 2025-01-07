import matplotlib.pyplot as plt

#### Code to plot the data of Fig.7
#### Author: Vladlen Galetsky, updated: 7/1/2025
#################################################
"""
Renormalization is needed to start counting the storage time at t=0 for all codes to be comparable, despite until storage different codes take different initialization, merging and splitting times.



Parameters:
tc0: Time steps for unencoded Bell state
tc1: Time for BS d=3
tc2: Time for unrotated S d=3
tc3: Time for rotated S d=3

tc3b: Extention of storage time for rotated S d=3 to have more iterations (to cover the BS d=3 time)
tc2b: Extention of storage time for unrotated S d=3 to have more iterations (to cover the BS d=3 time)
Vec_Y2b: Logical error rate for unrotated S d=3 (extension)
Vec_Y3b: Logical error rate for rotated S d=3 (extension)

Std_p2b: Standard deviation of logical error rate for unrotated S d=3 (extended)
Std_p3b: Standard deviation of logical error rate for rotated S d=3 (extended)

Std_m2b: Mean logical error rate for unrotated S d=3 (extension)
Std_m3b: Mean logical error rate for rotated S d=3 (extension)


S_Bell: Error of unencoded bell state

Std_m: Standard deviation (lower bound) of logical error rate for BS d=3 
Std_m2: Standard deviation (lower bound) of logical error rate for unrotated S d=3 
Std_m3: Standard deviation (lower bound) of logical error rate for rotated S d=3 

Std_p: Standard deviation (upper bound) of logical error rate for BS d=3 
Std_p2: Standard deviation (upper bound) of logical error rate for unrotated S d=3 
Std_p3: Standard deviation (upper bound) of logical error rate for rotated S d=3 

Vec_Y: Logical error rate BS d=3
Vec_Y2: Logical error rate unrotated S d=3
Vec_Y3: Logical error rate rotated S d=3

"""

#Plot distance parameters between scheme 1 and scheme 2
########## 5km 1 0.5 -ready

tc0 = [9.014322970853389, 11.416322970853388, 13.818322970853389, 16.220322970853392, 18.622322970853393, 21.02432297085339, 23.42632297085339, 25.828322970853392, 28.230322970853393, 30.632322970853394, 33.03432297085339]
tc1 = [62.22732297085339, 71.59532297085339, 80.96332297085338, 90.33132297085338, 99.69932297085339, 109.06732297085338, 118.43532297085339, 127.80332297085339, 137.1713229708534, 146.5393229708534, 155.9073229708534]
tc2 = [26.82732297085339, 30.29532297085339, 33.76332297085339, 37.23132297085339, 40.699322970853395, 44.16732297085339, 47.63532297085339, 51.10332297085339, 54.571322970853394, 58.03932297085339, 61.50732297085339]
tc3 = [26.82732297085339, 30.29532297085339, 33.76332297085339, 37.23132297085339, 40.699322970853395, 44.16732297085339, 47.63532297085339, 51.10332297085339, 54.571322970853394, 58.03932297085339, 61.50732297085339]

S_Bell = [0.037786, 0.044804, 0.05132, 0.05783, 0.064702, 0.071116, 0.077003, 0.083029, 0.089796, 0.095361, 0.10141]

Std_p = [0.36459425, 0.38722156, 0.40831963, 0.42777770000000004, 0.44542032, 0.46314355, 0.47883664000000004, 0.49323675, 0.5082686599999999, 0.5210774, 0.53402083]
Std_p2 = [0.16398743, 0.17408858, 0.18230273000000002, 0.19058058, 0.19760259, 0.20540072, 0.21102579999999999, 0.21871778, 0.22550125, 0.23315578, 0.24028931]
Std_p3 = [0.14217279000000002, 0.15131782000000002, 0.15784481, 0.16426775, 0.17086305, 0.17792859, 0.18379295, 0.19001604, 0.19654625, 0.20266163, 0.20853531]

Std_m = [0.36101283, 0.38360637, 0.40466490000000005, 0.42410982, 0.44173621999999996, 0.45944057, 0.47512528000000004, 0.48951544, 0.50454575, 0.5173609, 0.53030232]
Std_m2 = [0.16124656, 0.17127773999999998, 0.17943539, 0.18766592999999998, 0.19464685, 0.20239889000000003, 0.20799389999999998, 0.21565101, 0.22240403, 0.2300205, 0.23712321]
Std_m3 = [0.13958361, 0.14865948, 0.15513891, 0.16152222, 0.16806859, 0.17508897, 0.18092170000000002, 0.18710376, 0.19360054000000002, 0.19968060999999998, 0.2055201]

Vec_Y = [0.362802, 0.385414, 0.406491, 0.425943, 0.443578, 0.461292, 0.476981, 0.491376, 0.506406, 0.51922, 0.532161]
Vec_Y2 = [0.162614, 0.17268, 0.180866, 0.18912, 0.196123, 0.203897, 0.209506, 0.217182, 0.22395, 0.231586, 0.238704]
Vec_Y3 = [0.140875, 0.149986, 0.156489, 0.162892, 0.169463, 0.176506, 0.182354, 0.188556, 0.195071, 0.201168, 0.207025]




tc0 = [x - tc0[0] for x in tc0]
tc1 = [x - tc1[0] for x in tc1]
tc2 = [x - tc2[0] for x in tc2]
tc3 = [x - tc3[0] for x in tc3]



fig, ax = plt.subplots(1, 1)


ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"v",linestyle='-',color="black",label="Protocol 1: BS [[18,2,3]] ")

ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"v",linestyle='-',color="blue",label="Protocol 1: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"v",linestyle='-',color="green",label="Protocol 1: rotated S [[18,2,3]]")

#ax.plot(tc0,S_Bell,linestyle='--',color="red", label="Unencoded Bell-State")
#scheme 2, 1,0.5 D=5 km


tc0     = [28.22761485426695, 30.62961485426695, 33.03161485426695, 35.433614854266956, 37.83561485426696, 40.23761485426695, 42.63961485426695, 45.04161485426695, 47.443614854266954, 49.845614854266955, 52.24761485426695]
tc1     = [129.47384456280085, 138.84184456280084, 148.20984456280087, 157.57784456280086, 166.94584456280086, 176.31384456280085, 185.68184456280085, 195.04984456280084, 204.41784456280084, 213.78584456280086, 223.15384456280086]
tc2     = [94.07384456280084, 97.54184456280085, 101.00984456280084, 104.47784456280084, 107.94584456280084, 111.41384456280083, 114.88184456280084, 118.34984456280084, 121.81784456280084, 125.28584456280085, 128.75384456280085]
tc3     = [94.07384456280084, 97.54184456280085, 101.00984456280084, 104.47784456280084, 107.94584456280084, 111.41384456280083, 114.88184456280084, 118.34984456280084, 121.81784456280084, 125.28584456280085, 128.75384456280085]

S_Bell  = [0.02411, 0.030737, 0.038085, 0.044376, 0.051483, 0.057833, 0.064266, 0.070834, 0.076782, 0.082875, 0.089249]


Std_p   = [0.410659, 0.43240657, 0.45273990000000003, 0.47202175, 0.48890946, 0.50468483, 0.52049611, 0.53414001, 0.54811738, 0.55923221, 0.57143003]
Std_p2  = [0.22152629999999998, 0.22952812, 0.23669926000000002, 0.24357884, 0.24930975, 0.2575807, 0.26317124, 0.27014848999999996, 0.2765809, 0.28365687, 0.28963648999999997]
Std_p3  = [0.22498068, 0.23082487, 0.23737172, 0.24311601000000002, 0.24824060999999997, 0.25470737, 0.26008654999999997, 0.26578512, 0.2705627, 0.27611391999999996, 0.2822171]

Std_m   = [0.40699367999999997, 0.4287345, 0.44904336, 0.46829909000000003, 0.48519204, 0.50097695, 0.51676848, 0.53042204, 0.54442996, 0.5555448900000001, 0.5677631999999999]
Std_m2  = [0.21844547, 0.22641023000000002, 0.23355438, 0.24039104, 0.24609942, 0.25433516, 0.25989504, 0.26685535, 0.27326745, 0.28031352000000004, 0.28626126]
Std_m3  = [0.22188235, 0.22769604000000002, 0.23421304999999998, 0.23993606, 0.24503816, 0.25147677, 0.25682932, 0.26250453999999995, 0.26726475, 0.27279384, 0.27887625]

Vec_Y   = [0.408827, 0.43057, 0.450891, 0.47016, 0.48705, 0.50283, 0.518633, 0.532282, 0.546275, 0.557389, 0.569597]
Vec_Y2  = [0.219982, 0.227966, 0.235125, 0.241982, 0.247702, 0.255956, 0.261531, 0.268499, 0.274922, 0.281983, 0.287947]
Vec_Y3  = [0.223429, 0.229258, 0.23579, 0.241524, 0.246637, 0.253089, 0.258456, 0.264143, 0.268911, 0.274452, 0.280544]


tc1=tc1[:4]
Std_m=Std_m[:4]
Std_p=Std_p[:4]
Vec_Y=Vec_Y[:4]



tc0 = [x - tc0[0] for x in tc0]
tc1 = [x - tc1[0] for x in tc1]
tc2 = [x - tc2[0] for x in tc2]
tc3 = [x - tc3[0] for x in tc3]

################################################
ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"x",linestyle='-',color="black",label="Protocol 2: BS [[18,2,3]]")

ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"x",linestyle='-',color="blue",label="Protocol 2: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"x",linestyle='-',color="green",label="Protocol 2: rotated S [[18,2,3]]")




plt.xlim(0,tc1[-1])
plt.ylim(0,Vec_Y[-1])

ax.grid(linestyle='--')
plt.legend(fontsize = 9)

props = dict(facecolor='white', alpha=0.5)

ax.text(0.410, 0.7, "$T_{1}, T_{2}, 0.5p_{err_{H}}, 0.5p_{err_{CX}}, 0.5p_{err_{M}}, D=5 km$", transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

# Manually add the first legend back to the plot

ax.set_xlabel("Storage time [$\mu$s]",fontsize=14)
ax.set_ylabel("Logical Error Rate",fontsize=14)
plt.show()





#############################################################################################
# d=10 km 1, 0.5, scheme 1 -ready
tc0     = [9.014322970853389, 11.416322970853388, 13.818322970853389, 16.220322970853392, 18.622322970853393, 21.02432297085339, 23.42632297085339, 25.828322970853392, 28.230322970853393, 30.632322970853394, 33.03432297085339]
tc1     = [62.22732297085339, 71.59532297085339, 80.96332297085338, 90.33132297085338, 99.69932297085339, 109.06732297085338, 118.43532297085339, 127.80332297085339, 137.1713229708534, 146.5393229708534, 155.9073229708534]
tc2     = [26.82732297085339, 30.29532297085339, 33.76332297085339, 37.23132297085339, 40.699322970853395, 44.16732297085339, 47.63532297085339, 51.10332297085339, 54.571322970853394, 58.03932297085339, 61.50732297085339]
tc3     = [26.82732297085339, 30.29532297085339, 33.76332297085339, 37.23132297085339, 40.699322970853395, 44.16732297085339, 47.63532297085339, 51.10332297085339, 54.571322970853394, 58.03932297085339, 61.50732297085339]


Std_p  = [0.52944239, 0.54175689, 0.5542398000000001, 0.5648295600000001, 0.5737806400000001, 0.58341461, 0.59208402, 0.60118003, 0.60826982, 0.61608617, 0.62291587]
Std_p2   = [0.32780458, 0.33467821000000003, 0.34150863, 0.34568532, 0.35101396999999995, 0.35541954, 0.36219643, 0.36597637, 0.37117699, 0.37631408, 0.3804169]
Std_p3  = [0.2848449, 0.29081614, 0.29631618, 0.30059381, 0.30675256, 0.31153728, 0.31523890000000004, 0.32083213, 0.32565164, 0.33115448999999997, 0.33490171]

Std_m   = [0.52573088, 0.53805982, 0.55054223, 0.5611473, 0.57010047, 0.5797515999999999, 0.58843326, 0.59753577, 0.60463593, 0.61247228, 0.6193143000000001]
Std_m2  = [0.3243227, 0.33118010999999997, 0.33799465, 0.34215490000000004, 0.34747468, 0.35185934, 0.35863568, 0.36239398, 0.36759585, 0.37270654999999997, 0.37680786]
Std_m3  = [0.28149257, 0.28743588, 0.29292882000000003, 0.29719398999999996, 0.30332874, 0.30810943, 0.3117908, 0.3173655, 0.32216928, 0.32767057, 0.3313947]

Vec_Y   = [0.527586, 0.539907, 0.552392, 0.56299, 0.571942, 0.581584, 0.590259, 0.599358, 0.606454, 0.614281, 0.621115]
Vec_Y2  = [0.326063, 0.332927, 0.33975, 0.343918, 0.349243, 0.353638, 0.360415, 0.364183, 0.369384, 0.374509, 0.378611]
Vec_Y3  = [0.283167, 0.289124, 0.29462, 0.298892, 0.305038, 0.309822, 0.313513, 0.319098, 0.323909, 0.329411, 0.333146]


tc1=tc1[:4]
Std_m=Std_m[:4]
Std_p=Std_p[:4]
Vec_Y=Vec_Y[:4]

tc0 = [x - tc0[0] for x in tc0]
tc1 = [x - tc1[0] for x in tc1]
tc2 = [x - tc2[0] for x in tc2]
tc3 = [x - tc3[0] for x in tc3]

fig, ax = plt.subplots(1, 1)


ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"v",linestyle='-',color="black",label="Protocol 1: BS [[18,2,3]] ")

ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"v",linestyle='-',color="blue",label="Protocol 1: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"v",linestyle='-',color="green",label="Protocol 1: rotated S [[18,2,3]]")


### 1 0.5 scheme 2 -ready

tc0     = [52.2442297085339, 54.6462297085339, 57.0482297085339, 59.450229708533904, 61.852229708533905, 64.2542297085339, 66.65622970853389, 69.0582297085339, 71.46022970853389, 73.8622297085339, 76.26422970853389]
tc1     = [201.5236891256017, 210.8916891256017, 220.2596891256017, 229.6276891256017, 238.9956891256017, 248.3636891256017, 257.7316891256017, 267.0996891256017, 276.4676891256017, 285.8356891256017, 295.2036891256017]
tc2     = [166.1236891256017, 169.5916891256017, 173.0596891256017, 176.52768912560168, 179.9956891256017, 183.4636891256017, 186.9316891256017, 190.3996891256017, 193.8676891256017, 197.3356891256017, 200.8036891256017]
tc3     = [166.1236891256017, 169.5916891256017, 173.0596891256017, 176.52768912560168, 179.9956891256017, 183.4636891256017, 186.9316891256017, 190.3996891256017, 193.8676891256017, 197.3356891256017, 200.8036891256017]

S_Bell  = [0.024065, 0.031197, 0.037571, 0.044351, 0.051474, 0.058184, 0.06463, 0.070642, 0.077361, 0.083379, 0.089472]

Vec_Y2  = [0.474326, 0.478859, 0.482533, 0.486108, 0.490185, 0.49376, 0.496931, 0.501298, 0.505396, 0.508543, 0.512014]

Std_p   = [0.5658574300000001, 0.58021373, 0.5932604399999999, 0.605746, 0.61657127, 0.6276443100000001, 0.63678742, 0.6457836800000001, 0.65337453, 0.65998563, 0.66689726]
Std_p2  = [0.47618413, 0.48071131, 0.48439125, 0.48796987, 0.49204343, 0.49561863, 0.49878581, 0.50315714, 0.50725255, 0.51039758, 0.51387552]
Std_p3  = [0.47590831, 0.48047142, 0.48252843, 0.48535177, 0.48813443, 0.49082546000000005, 0.4944998, 0.49833710999999997, 0.50095857, 0.50356692, 0.50743697]

Std_m   = [0.56217786, 0.5765340600000001, 0.58961046, 0.6021164, 0.61295352, 0.62404989, 0.6332149699999999, 0.6422162299999999, 0.64982696, 0.6564598, 0.663392]
Std_m2  = [0.47246663, 0.47700696000000004, 0.48067418, 0.48424667, 0.48832817, 0.49190096, 0.49507773, 0.49943989, 0.50353942, 0.50668841, 0.5101540800000001]
Std_m3  = [0.47220352000000004, 0.47676298, 0.47880892999999997, 0.48164007000000003, 0.48440908, 0.48712071999999995, 0.49078379, 0.49461809999999995, 0.49723738, 0.49984635, 0.50371172]

Vec_Y   = [0.564018, 0.578375, 0.591437, 0.603933, 0.614764, 0.625849, 0.635003, 0.644002, 0.651603, 0.658225, 0.665146]
Vec_Y2  = [0.474326, 0.478859, 0.482533, 0.486108, 0.490185, 0.49376, 0.496931, 0.501298, 0.505396, 0.508543, 0.512014]
Vec_Y3  = [0.474056, 0.478618, 0.480668, 0.483496, 0.486273, 0.488973, 0.492641, 0.496477, 0.499097, 0.501706, 0.505574]




tc1=tc1[:4]
Std_m=Std_m[:4]
Std_p=Std_p[:4]
Vec_Y=Vec_Y[:4]

tc0 = [x - tc0[0] for x in tc0]
tc1 = [x - tc1[0] for x in tc1]
tc2 = [x - tc2[0] for x in tc2]
tc3 = [x - tc3[0] for x in tc3]

################################################
ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"x",linestyle='-',color="black",label="Protocol 2: BS [[18,2,3]]")

ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"x",linestyle='-',color="blue",label="Protocol 2: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"x",linestyle='-',color="green",label="Protocol 2: rotated S [[18,2,3]]")




plt.xlim(0,tc1[-1])
plt.ylim(0,Vec_Y[-1])

ax.grid(linestyle='--')
plt.legend(fontsize = 9)

props = dict(facecolor='white', alpha=0.5)

ax.text(0.025, 0.425, "$T_{1}, T_{2}, 0.5p_{err_{H}}, 0.5p_{err_{CX}}, 0.5p_{err_{M}}, D=10 km$", transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

# Manually add the first legend back to the plot

ax.set_xlabel("Storage time [$\mu$s]",fontsize=14)
ax.set_ylabel("Logical Error Rate",fontsize=14)
plt.show()

##########################################################################################################################################
