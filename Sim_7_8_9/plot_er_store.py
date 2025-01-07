import matplotlib.pyplot as plt

#### Code to plot the data of Fig.9
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


#Protocol 1:
#DATA: 1, 1

tc0=[9.014322970853389, 11.416322970853388, 13.818322970853389, 16.220322970853392, 18.622322970853393, 21.02432297085339, 23.42632297085339, 
        25.828322970853392, 28.230322970853393, 30.632322970853394, 33.03432297085339, 
        35.436322970853396, 37.8383229708534, 40.24032297085339, 42.64232297085339, 
        45.04432297085339, 47.446322970853394, 49.848322970853395, 52.2503229708534, 
        54.6523229708534, 57.05432297085339, 59.45632297085339, 61.85832297085339, 
        64.2603229708534, 66.6623229708534, 69.0643229708534, 71.46632297085338, 
        73.86832297085338, 76.27032297085339, 78.67232297085339, 81.07432297085339, 
        83.47632297085339, 85.87832297085339, 88.28032297085339]
tc1 = [62.22732297085339, 71.59532297085339, 80.96332297085338, 90.33132297085338]
tc2 = [26.82732297085339, 30.29532297085339, 33.76332297085339, 37.23132297085339]
tc3 = [26.82732297085339, 30.29532297085339, 33.76332297085339, 37.23132297085339]


# Third line of output
tc3b = [37.23132297085339, 40.699322970853395, 44.16732297085339, 47.63532297085339, 
        51.10332297085339, 54.571322970853394, 58.03932297085339, 61.50732297085339, 
        64.97532297085338, 68.44332297085339, 71.91132297085339, 75.37932297085338, 
        78.84732297085338, 82.31532297085339, 85.78332297085338, 89.25132297085338]
tc2b = [37.23132297085339, 40.699322970853395, 44.16732297085339, 47.63532297085339, 
        51.10332297085339, 54.571322970853394, 58.03932297085339, 61.50732297085339, 
        64.97532297085338, 68.44332297085339, 71.91132297085339, 75.37932297085338, 
        78.84732297085338, 82.31532297085339, 85.78332297085338, 89.25132297085338]





Vec_Y2b = [0.143997, 0.162055, 0.17882, 0.196208, 0.21268, 0.227711, 0.243013, 0.257899]
Vec_Y3b = [0.132769, 0.147905, 0.162752, 0.176846, 0.190881, 0.20475, 0.218256, 0.230756]

Std_m2b = [0.14269341, 0.1606891, 0.17739669, 0.1947349, 0.21116284, 0.22615493, 0.24142251, 0.25627138]
Std_p2b = [0.14530606, 0.16342854, 0.18024789000000002, 0.19768739999999999, 0.21420349, 0.22927088, 0.24460695999999998, 0.25952955]

Std_m3b = [0.13151137, 0.1465854, 0.16138084, 0.17543082000000002, 0.18942196, 0.20325179000000002, 0.21672528, 0.22919593]
Std_p3b = [0.13403395999999998, 0.14923198999999998, 0.16412849, 0.17826775, 0.19234618, 0.20625311999999998, 0.21979223, 0.23232054000000002]






S_Bell = [0.053712, 0.060229, 0.06685, 0.07319, 0.07925, 0.085423, 0.091353, 0.097113, 0.103165, 0.109168, 
           0.114742, 0.120809, 0.126236, 0.131305, 0.137377, 0.141835, 0.147891, 
           0.152257, 0.157975, 0.162959, 0.168246, 0.172594, 0.176868, 0.18223, 
           0.186996, 0.191399, 0.196247, 0.200043, 0.204976, 0.209454, 0.212958, 
           0.217491, 0.220637, 0.226561]

Vec_Y2 = [0.085994, 0.107501, 0.125772, 0.144056]

Std_p = [0.40294846, 0.449072, 0.48871673, 0.52361931]
Std_p2 = [0.08704084, 0.10865712, 0.12700909, 0.14536735999999997]
Std_p3 = [0.08445239, 0.10142486, 0.11856541000000001, 0.133362]

Std_m = [0.3993027, 0.44538342999999997, 0.48500462, 0.51990249]
Std_m2 = [0.08495475999999999, 0.10635178999999999, 0.12454187, 0.14275260999999997]
Std_m3 = [0.08239539, 0.09919205, 0.11617175, 0.13084284]

Vec_Y = [0.401125, 0.447228, 0.48686, 0.521761]
Vec_Y2 = [0.085994, 0.107501, 0.125772, 0.144056]
Vec_Y3 = [0.08342, 0.100305, 0.117365, 0.132099]



tc2b = [x - tc2[0] for x in tc2b]
tc3b = [x - tc3[0] for x in tc3b]
tc3b=tc3b[:len(Std_m3b)]
tc2b=tc2b[:len(Std_m2b)]

tc0 = [x - tc0[0] for x in tc0]
tc1 = [x - tc1[0] for x in tc1]
tc2 = [x - tc2[0] for x in tc2]
tc3 = [x - tc3[0] for x in tc3]





fig, ax = plt.subplots(1, 1)


ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"v",linestyle='-',color="black",label="Protocol 1: BS [[18,2,3]] ")


ax.fill_between(tc2b, Std_m2b, Std_p2b,alpha=0.2,color="blue")
ax.plot(tc2b,Vec_Y2b,"v",linestyle='-',color="blue")

ax.fill_between(tc3b, Std_m3b, Std_p3b,alpha=0.2,color="green")
ax.plot(tc3b,Vec_Y3b,"v",linestyle='-',color="green")


ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"v",linestyle='-',color="blue",label="Protocol 1: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"v",linestyle='-',color="green",label="Protocol 1: rotated S [[18,2,3]]")

ax.plot(tc0,S_Bell,linestyle='--',color="red", label="Unencoded Bell-State")







####################################################



Vec_Y2b = [0.150552, 0.167952, 0.184631, 0.201689, 0.217316, 0.232907, 0.247927, 0.262195]
Vec_Y3b = [0.1391, 0.154594, 0.168235, 0.18348, 0.196989, 0.210191, 0.224211, 0.235409]

Std_m2b = [0.14922767, 0.16656394, 0.18318931, 0.20019995000000002, 0.21578276000000002, 0.23133838, 0.24632309, 0.26056216]
Std_p2b = [0.15188332999999998, 0.16934562, 0.18607727, 0.20318341, 0.21885404, 0.23448160999999998, 0.24953554, 0.26383017]

Std_m3b = [0.13781452, 0.15325294, 0.16684749, 0.18204188000000002, 0.19551389000000002, 0.20867727, 0.22266342, 0.23383632000000001]
Std_p3b = [0.14039274, 0.15593996, 0.16962786, 0.184924, 0.19847073, 0.21171098000000002, 0.22576362, 0.23698708]

Vec_Y2 = [0.093352, 0.112898, 0.131316, 0.149848]

Std_p = [0.41949309999999995, 0.46368463, 0.5017878, 0.53464613]
Std_p2 = [0.09443773, 0.11408027, 0.13257354999999998, 0.15117881]
Std_p3 = [0.09324151, 0.11003136999999999, 0.12486813000000001, 0.14081117]

Std_m = [0.41582403, 0.45997978, 0.49806763, 0.5309283100000001]
Std_m2 = [0.09227353999999999, 0.11172275, 0.13006291, 0.14852328]
Std_m3 = [0.09109151, 0.10771343, 0.12242067, 0.13822901999999998]

Vec_Y = [0.417658, 0.461832, 0.499928, 0.532788]
Vec_Y2 = [0.093352, 0.112898, 0.131316, 0.149848]
Vec_Y3 = [0.092163, 0.108869, 0.123641, 0.139517]


###################################

ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"x",linestyle='-',color="black",label="Protocol 2: BS [[18,2,3]]")


ax.fill_between(tc2b, Std_m2b, Std_p2b,alpha=0.2,color="blue")
ax.plot(tc2b,Vec_Y2b,"x",linestyle='-',color="blue")

ax.fill_between(tc3b, Std_m3b, Std_p3b,alpha=0.2,color="green")
ax.plot(tc3b,Vec_Y3b,"x",linestyle='-',color="green")


ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"x",linestyle='-',color="blue",label="Protocol 2: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"x",linestyle='-',color="green",label="Protocol 2: rotated S [[18,2,3]]")


plt.xlim(0,tc1[-1])
plt.ylim(0,Vec_Y[-1])

ax.grid(linestyle='--')
plt.legend(fontsize = 9)

props = dict(facecolor='white', alpha=0.5)

ax.text(0.020, 0.975, "$T_{1}, T_{2}, p_{err_{H}}, p_{err_{CX}}, p_{err_{M}}$", transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

# Manually add the first legend back to the plot

ax.set_xlabel("Storage time [$\mu$s]",fontsize=14)
ax.set_ylabel("Logical Error Rate",fontsize=14)
plt.show()



#DATA 1,0.5




Vec_Y2b = [0.076113, 0.084562, 0.093467, 0.102241, 0.110913, 0.119517, 0.128091, 0.136003]
Vec_Y3b = [0.068507, 0.075809, 0.083639, 0.091311, 0.09929, 0.106632, 0.114013, 0.121121]

Std_m2b = [0.07512949000000001, 0.08352933, 0.0923866, 0.10111849, 0.10974766, 0.1183142, 0.12685189, 0.13473229]
Std_p2b = [0.07710441999999999, 0.08560174000000001, 0.09455461999999999, 0.10337049000000001, 0.11208483, 0.12072661, 0.12933722, 0.13728003]

Std_m3b = [0.06757205000000001, 0.07482838, 0.08261418, 0.09024227, 0.09818183999999999, 0.10548836, 0.11283347, 0.11991128]
Std_p3b = [0.06945051, 0.076798, 0.084672, 0.0923873, 0.100405, 0.10778305, 0.11519908, 0.12233769]
S_Bell = [0.037825, 0.044567, 0.051324, 0.057622, 0.06422, 0.071134, 0.077116, 0.082924, 0.089733, 0.095752, 0.100968, 
           0.107004, 0.113168, 0.117825, 0.124207, 0.129847, 0.135966, 0.140697, 0.145362, 
           0.151104, 0.155753, 0.160867, 0.166104, 0.170848, 0.175569, 0.180656, 0.185341, 
           0.190703, 0.194464, 0.198901, 0.202864, 0.206952, 0.211828, 0.215924]

Std_p = [0.20318251999999998, 0.23586277, 0.26647415999999996, 0.29524262]
Std_p2 = [0.04737994, 0.057942839999999995, 0.06735548, 0.0766323]
Std_p3 = [0.044638800000000006, 0.05359247, 0.06184155, 0.06937524]

Std_m = [0.20019778, 0.23270891, 0.26318821000000003, 0.29186192]
Std_m2 = [0.04581242, 0.05621684, 0.06550187, 0.07466564]
Std_m3 = [0.04311525, 0.05192956, 0.060060220000000004, 0.06749816]

Vec_Y = [0.201687, 0.234283, 0.264829, 0.29355]
Vec_Y2 = [0.046592, 0.057076, 0.066425, 0.075645]
Vec_Y3 = [0.043873, 0.052757, 0.060947, 0.068433]


fig, ax = plt.subplots(1, 1)


ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"v",linestyle='-',color="black",label="Protocol 1: BS [[18,2,3]] ")


ax.fill_between(tc2b, Std_m2b, Std_p2b,alpha=0.2,color="blue")
ax.plot(tc2b,Vec_Y2b,"v",linestyle='-',color="blue")

ax.fill_between(tc3b, Std_m3b, Std_p3b,alpha=0.2,color="green")
ax.plot(tc3b,Vec_Y3b,"v",linestyle='-',color="green")


ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"v",linestyle='-',color="blue",label="Protocol 1: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"v",linestyle='-',color="green",label="Protocol 1: rotated S [[18,2,3]]")

ax.plot(tc0,S_Bell,linestyle='--',color="red", label="Unencoded Bell-State")


#1, 0.5


Vec_Y2b = [0.074432, 0.083809, 0.09213, 0.100637, 0.109272, 0.118127, 0.127436, 0.135884]
Vec_Y3b = [0.070727, 0.078932, 0.086557, 0.094255, 0.101862, 0.108783, 0.11594, 0.123519]

Std_m2b = [0.0734582, 0.08278209, 0.09105907, 0.09952234, 0.10811469, 0.11693060000000001, 0.12619965, 0.13461128]
Std_p2b = [0.07541319, 0.08484262, 0.09320878, 0.10175933, 0.11043597, 0.11932963, 0.12867833, 0.13716255]

Std_m3b = [0.06977788, 0.07793338000000001, 0.08551518, 0.09317117, 0.10074213, 0.10762874, 0.11475236999999999, 0.1222972]
Std_p3b = [0.07168411999999999, 0.07993807000000001, 0.08760669, 0.09534595, 0.10298958999999999, 0.10994484, 0.11713454, 0.12474778]

Std_p = [0.2141014, 0.2458879, 0.27578982, 0.30348364]
Std_p2 = [0.047106470000000004, 0.056606779999999995, 0.06614585, 0.07544895]
Std_p3 = [0.04799287, 0.05579432, 0.06430759, 0.07225402]

Std_m = [0.21105919, 0.24268876, 0.27248036999999997, 0.30006325]
Std_m2 = [0.04554412, 0.05490129, 0.06431002, 0.07349735]
Std_m3 = [0.04641696, 0.05410002, 0.06249451, 0.0703415]

Vec_Y = [0.212577, 0.244286, 0.274134, 0.301771]
Vec_Y2 = [0.046321, 0.05575, 0.065224, 0.074469]
Vec_Y3 = [0.047201, 0.054943, 0.063397, 0.071294]



################################################
ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"x",linestyle='-',color="black",label="Protocol 2: BS [[18,2,3]]")


ax.fill_between(tc2b, Std_m2b, Std_p2b,alpha=0.2,color="blue")
ax.plot(tc2b,Vec_Y2b,"x",linestyle='-',color="blue")

ax.fill_between(tc3b, Std_m3b, Std_p3b,alpha=0.2,color="green")
ax.plot(tc3b,Vec_Y3b,"x",linestyle='-',color="green")


ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"x",linestyle='-',color="blue",label="Protocol 2: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"x",linestyle='-',color="green",label="Protocol 2: rotated S [[18,2,3]]")


second_legend = plt.legend("here", loc='lower right')


plt.xlim(0,tc1[-1])
plt.ylim(0,Vec_Y[-1])

ax.grid(linestyle='--')
plt.legend(fontsize = 9)

props = dict(facecolor='white', alpha=0.5)

ax.text(0.020, 0.975, "$T_{1}, T_{2}, 0.5p_{err_{H}}, 0.5p_{err_{CX}}, 0.5p_{err_{M}}$", transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

# Manually add the first legend back to the plot

ax.set_xlabel("Storage time [$\mu$s]",fontsize=14)
ax.set_ylabel("Logical Error Rate",fontsize=14)
plt.show()




#DATA 1,0.1



Vec_Y2b = [0.031656, 0.034679, 0.038038, 0.041142, 0.044481, 0.047905, 0.050789, 0.054016]
Vec_Y3b = [0.027321, 0.029639, 0.032191, 0.034688, 0.037569, 0.039945, 0.042017, 0.044861]

Std_m2b = [0.03100941, 0.03400286, 0.03733125, 0.0404078, 0.04371856, 0.04711538, 0.04997717, 0.05317947]
Std_p2b = [0.03231106, 0.03536367, 0.03875331, 0.041884870000000005, 0.04525199, 0.04870322, 0.051609050000000004, 0.05486083]

Std_m3b = [0.02671937, 0.029012740000000002, 0.03153872, 0.03401157, 0.0368659, 0.03922114, 0.0412752, 0.044095800000000004]
Std_p3b = [0.02793128, 0.030273900000000003, 0.03285218, 0.035372940000000005, 0.038280669999999996, 0.04067736, 0.04276751, 0.04563424]

S_Bell = [0.024958, 0.032138, 0.038714,0.045674, 0.052135, 0.059071, 0.06555, 0.07129, 0.077928, 0.083876, 0.090375, 0.096337, 0.101981, 0.10809, 0.113989, 0.119664, 0.125795, 0.130506, 0.135507, 0.140855, 0.146439, 0.151882, 0.156029, 0.161869, 0.16622, 0.171413, 0.176969, 0.181014, 0.186102, 0.190854, 0.195361, 0.200635, 0.203821, 0.208177]


Vec_Y2 = [0.021664, 0.024709, 0.028295, 0.031484]
Std_p = [0.0568579, 0.06574039999999999, 0.07375313, 0.0816531]
Std_p2 = [0.022209720000000002, 0.02529037, 0.02891604, 0.03213784]
Std_p3 = [0.0197115, 0.02243427, 0.025212650000000003, 0.02784139]

Std_m = [0.05514968, 0.06390926, 0.07182252, 0.07962644]
Std_m2 = [0.02112716, 0.02413618, 0.02768277, 0.03083879]
Std_m3 = [0.0186914, 0.021346470000000003, 0.024060130000000002, 0.02663144]

Vec_Y = [0.056, 0.064821, 0.072784, 0.080636]
Vec_Y2 = [0.021664, 0.024709, 0.028295, 0.031484]
Vec_Y3 = [0.019197, 0.021886, 0.024632, 0.027232]


fig, ax = plt.subplots(1, 1)


ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"v",linestyle='-',color="black",label="Protocol 1: BS [[18,2,3]] ")


ax.fill_between(tc2b, Std_m2b, Std_p2b,alpha=0.2,color="blue")
ax.plot(tc2b,Vec_Y2b,"v",linestyle='-',color="blue")

ax.fill_between(tc3b, Std_m3b, Std_p3b,alpha=0.2,color="green")
ax.plot(tc3b,Vec_Y3b,"v",linestyle='-',color="green")


ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"v",linestyle='-',color="blue",label="Protocol 1: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"v",linestyle='-',color="green",label="Protocol 1: rotated S [[18,2,3]]")

ax.plot(tc0,S_Bell,linestyle='--',color="red", label="Unencoded Bell-State")



#Data 1,0.1


Vec_Y2b = [0.02748, 0.030593, 0.034203, 0.037144, 0.040554, 0.043595, 0.046869, 0.049825]
Vec_Y3b = [0.027386, 0.02993, 0.032644, 0.035014, 0.037504, 0.040454, 0.042612, 0.045007]

Std_m2b = [0.02687692, 0.02995663, 0.03353148, 0.036445080000000005, 0.0398248, 0.04283998000000001, 0.04608733, 0.049020339999999996]
Std_p2b = [0.02809183, 0.031238119999999998, 0.034883029999999995, 0.03785153, 0.04129193, 0.044358339999999996, 0.04765888, 0.050638]

Std_m3b = [0.02678372, 0.02930093, 0.03198714, 0.03433506, 0.03680171, 0.03972528, 0.041864910000000005, 0.044240510000000004]
Std_p3b = [0.02799695, 0.030567650000000002, 0.03330951, 0.03570162, 0.038214910000000005, 0.04119085, 0.04336768, 0.04578175]



S_Bell2 = [0.031813, 0.038846, 0.045734, 0.052156, 0.058952, 0.065148, 0.07143, 0.077791, 0.084476, 0.090749, 0.096596, 0.101925, 0.10822, 0.113454, 0.11984, 0.125499, 0.130024, 0.136176, 0.140948, 0.146399, 0.151847, 0.156255, 0.162413, 0.167005, 0.171424, 0.176397, 0.181577, 0.185982, 0.190661, 0.194472, 0.198995]
S_Bell = [0.010903, 0.017992, 0.024973, 0.031565]

Vec_Y2 = [0.017426, 0.020826, 0.024104, 0.027401]
Std_p = [0.05326627, 0.061382019999999995, 0.06966488, 0.07785825]
Std_p2 = [0.01791693, 0.02136109, 0.02467844, 0.02801225]
Std_p3 = [0.02040031, 0.02283401, 0.02545816, 0.02781525]

Std_m = [0.051608089999999995, 0.05960803, 0.0677835, 0.07587607]
Std_m2 = [0.01694398, 0.02029958, 0.02353831, 0.02679848]
Std_m3 = [0.01936251, 0.02173668, 0.024300509999999997, 0.02660535]

Vec_Y = [0.052433, 0.060491, 0.06872, 0.076863]
Vec_Y2 = [0.017426, 0.020826, 0.024104, 0.027401]
Vec_Y3 = [0.019877, 0.022281, 0.024875, 0.027206]
###############################################################



################################################
ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"x",linestyle='-',color="black",label="Protocol 2: BS [[18,2,3]]")


ax.fill_between(tc2b, Std_m2b, Std_p2b,alpha=0.2,color="blue")
ax.plot(tc2b,Vec_Y2b,"x",linestyle='-',color="blue")

ax.fill_between(tc3b, Std_m3b, Std_p3b,alpha=0.2,color="green")
ax.plot(tc3b,Vec_Y3b,"x",linestyle='-',color="green")


ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"x",linestyle='-',color="blue",label="Protocol 2: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"x",linestyle='-',color="green",label="Protocol 2: rotated S [[18,2,3]]")


second_legend = plt.legend("here", loc='lower right')


plt.xlim(0,tc1[-1])
plt.ylim(0,Vec_Y[-1])

ax.grid(linestyle='--')
plt.legend(fontsize = 9)

props = dict(facecolor='white', alpha=0.5)

ax.text(0.020, 0.975, "$T_{1}, T_{2}, 0.1p_{err_{H}}, 0.1p_{err_{CX}}, 0.1p_{err_{M}}$", transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

# Manually add the first legend back to the plot

ax.set_xlabel("Storage time [$\mu$s]",fontsize=14)
ax.set_ylabel("Logical Error Rate",fontsize=14)
plt.show()


#DATA 2,1


Vec_Y2b = [0.099828, 0.114202, 0.127619, 0.141452, 0.153564, 0.166316, 0.178592, 0.191357]
Vec_Y3b = [0.093031, 0.105018, 0.116908, 0.128011, 0.13965, 0.150511, 0.160817, 0.171767]

Std_m2b = [0.09871676, 0.11302232000000001, 0.126382, 0.14015695, 0.1522285, 0.16493417000000002, 0.17717182999999997, 0.18989537]
Std_p2b = [0.10094615, 0.11538868, 0.12886181, 0.14275131, 0.1549064, 0.16770253, 0.18001798000000002, 0.19282367]

Std_m3b = [0.09195519, 0.10388084, 0.1157164, 0.12677175, 0.13836293, 0.14918508, 0.15945493, 0.17036403]
Std_p3b = [0.09411435, 0.10616176, 0.11810645, 0.12925695, 0.14094273000000002, 0.1518432, 0.16218560999999998, 0.1731762]

S_Bell = [0.043875, 0.047039, 0.050384,0.053827, 0.056255, 0.060017, 0.063825, 0.066712, 0.069556, 0.073166, 0.075445, 0.079075, 0.082505, 0.085583, 0.088509, 0.091837, 0.094531, 0.097325, 0.100168, 0.103468, 0.105796, 0.108499, 0.11162, 0.115069, 0.117564, 0.120376, 0.123247, 0.125949, 0.129165, 0.131558, 0.133926, 0.137072, 0.140273, 0.142308]


Vec_Y2 = [0.05447, 0.070967, 0.084964, 0.099724]
Std_p = [0.3358182, 0.38102905, 0.42036265, 0.45411776000000004]
Std_p2 = [0.05531732, 0.07192628999999999, 0.08600505, 0.10084314]
Std_p3 = [0.05562355, 0.06923225999999999, 0.08201227999999999, 0.0940809]

Std_m = [0.33230815999999996, 0.37741664, 0.41669522999999997, 0.45040846999999995]
Std_m2 = [0.053630989999999996, 0.07001594, 0.08393081, 0.09861164]
Std_m3 = [0.05393278, 0.06735611, 0.07998286, 0.09191913]

Vec_Y = [0.334062, 0.379221, 0.418528, 0.452263]
Vec_Y2 = [0.05447, 0.070967, 0.084964, 0.099724]
Vec_Y3 = [0.054774, 0.06829, 0.080994, 0.092996]


fig, ax = plt.subplots(1, 1)


ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"v",linestyle='-',color="black",label="Protocol 1: BS [[18,2,3]] ")


ax.fill_between(tc2b, Std_m2b, Std_p2b,alpha=0.2,color="blue")
ax.plot(tc2b,Vec_Y2b,"v",linestyle='-',color="blue")

ax.fill_between(tc3b, Std_m3b, Std_p3b,alpha=0.2,color="green")
ax.plot(tc3b,Vec_Y3b,"v",linestyle='-',color="green")


ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"v",linestyle='-',color="blue",label="Protocol 1: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"v",linestyle='-',color="green",label="Protocol 1: rotated S [[18,2,3]]")

ax.plot(tc0,S_Bell,linestyle='--',color="red", label="Unencoded Bell-State")




#Data 2,1


Vec_Y2b = [0.103085, 0.116342, 0.130877, 0.143504, 0.156994, 0.1692, 0.181814, 0.19359]
Vec_Y3b = [0.095521, 0.107046, 0.118649, 0.129927, 0.141585, 0.152311, 0.162591, 0.173089]

Std_m2b = [0.10195842, 0.11515348, 0.12962539, 0.14220293, 0.15564253, 0.16780826999999998, 0.18038052, 0.19212089999999998]
Std_p2b = [0.10421883, 0.1175378, 0.13213608, 0.14481176999999998, 0.15835264999999998, 0.17059748000000002, 0.18325369, 0.19506483]

Std_m3b = [0.09443222, 0.10589842999999999, 0.11744841, 0.12868017, 0.14029214, 0.15097584, 0.16122341, 0.17168388]
Std_p3b = [0.09661732, 0.10820078, 0.11985677, 0.13118201000000002, 0.1428837, 0.15365296, 0.16396566, 0.17449999]


Vec_Y2 = [0.058848, 0.072836, 0.088353, 0.102767]
Std_p = [0.3460475, 0.38958695, 0.42843408, 0.46196174]
Std_p2 = [0.05972767, 0.07380605999999999, 0.08941335, 0.10390031]
Std_p3 = [0.05895995, 0.07156856, 0.08398419, 0.09635473]

Std_m = [0.34251437, 0.38596934000000005, 0.42476216, 0.45826554999999997]
Std_m2 = [0.05797629, 0.07187309, 0.08730016, 0.10164129]
Std_m3 = [0.05722017, 0.06966180000000001, 0.08193173, 0.09417248]

Vec_Y = [0.344279, 0.387777, 0.426597, 0.460113]
Vec_Y2 = [0.058848, 0.072836, 0.088353, 0.102767]
Vec_Y3 = [0.058086, 0.070611, 0.082954, 0.09526]




################################################
ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"x",linestyle='-',color="black",label="Protocol 2: BS [[18,2,3]]")


ax.fill_between(tc2b, Std_m2b, Std_p2b,alpha=0.2,color="blue")
ax.plot(tc2b,Vec_Y2b,"x",linestyle='-',color="blue")

ax.fill_between(tc3b, Std_m3b, Std_p3b,alpha=0.2,color="green")
ax.plot(tc3b,Vec_Y3b,"x",linestyle='-',color="green")


ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"x",linestyle='-',color="blue",label="Protocol 2: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"x",linestyle='-',color="green",label="Protocol 2: rotated S [[18,2,3]]")


second_legend = plt.legend("here", loc='lower right')


plt.xlim(0,tc1[-1])
plt.ylim(0,Vec_Y[-1])

ax.grid(linestyle='--')
plt.legend(fontsize = 9)

props = dict(facecolor='white', alpha=0.5)

ax.text(0.020, 0.975, "$2T_{1}, 2T_{2}, p_{err_{H}}, p_{err_{CX}}, p_{err_{M}}$", transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

# Manually add the first legend back to the plot

ax.set_xlabel("Storage time [$\mu$s]",fontsize=14)
ax.set_ylabel("Logical Error Rate",fontsize=14)
plt.show()
#########################################################




#DATA 2,0.5



Vec_Y2b = [0.0409, 0.046929, 0.052367, 0.057973, 0.063806, 0.068874, 0.074813, 0.07986]
Vec_Y3b = [0.038689, 0.043429, 0.047978, 0.052951, 0.057823, 0.062588, 0.067095, 0.071624]

Std_m2b = [0.04016796, 0.04614679, 0.05154288, 0.05710779, 0.0629002, 0.06793666000000001, 0.07383856, 0.07885466000000001]
Std_p2b = [0.04164051, 0.04771959, 0.053198959999999997, 0.05884609, 0.06471931, 0.06981928999999999, 0.07579522, 0.08087219]

Std_m3b = [0.03797586, 0.04267603, 0.047187730000000004, 0.052122720000000004, 0.05695933, 0.061691300000000004, 0.06616881, 0.07066978]
Std_p3b = [0.03941074, 0.04419066000000001, 0.048776660000000006, 0.05378686, 0.05869463, 0.06349315, 0.06802916, 0.0725866]

S_Bell = [0.027411, 0.031158, 0.034749, 0.03766, 0.04134, 0.04434, 0.047976, 0.051126, 0.054283, 0.057904, 0.060724, 0.064371, 0.067743, 0.070855, 0.073944, 0.07709, 0.080057, 0.083246, 0.08612, 0.08998, 0.092843, 0.095584, 0.098486, 0.10139, 0.104022, 0.107227, 0.110285, 0.112955, 0.115912, 0.117942, 0.120863, 0.12416, 0.1268, 0.12979]

Vec_Y2 = [0.023145, 0.029848, 0.035418, 0.041239]
Std_p = [0.13990276999999998, 0.16539222, 0.18929018, 0.2111772]
Std_p2 = [0.02370812, 0.03048524, 0.03610928, 0.0419822]
Std_p3 = [0.023808669999999997, 0.02916258, 0.03433825, 0.039067769999999995]

Std_m = [0.13732751999999998, 0.16263032, 0.18638225, 0.20815124]
Std_m2 = [0.0225906, 0.02921947, 0.03473516, 0.040504019999999995]
Std_m3 = [0.02268811, 0.027923959999999998, 0.03299631, 0.037640690000000004]

Vec_Y = [0.138611, 0.164008, 0.187833, 0.209661]
Vec_Y2 = [0.023145, 0.029848, 0.035418, 0.041239]
Vec_Y3 = [0.023244, 0.028539, 0.033663, 0.03835]


fig, ax = plt.subplots(1, 1)


ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"v",linestyle='-',color="black",label="Protocol 1: BS [[18,2,3]] ")


ax.fill_between(tc2b, Std_m2b, Std_p2b,alpha=0.2,color="blue")
ax.plot(tc2b,Vec_Y2b,"v",linestyle='-',color="blue")

ax.fill_between(tc3b, Std_m3b, Std_p3b,alpha=0.2,color="green")
ax.plot(tc3b,Vec_Y3b,"v",linestyle='-',color="green")


ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"v",linestyle='-',color="blue",label="Protocol 1: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"v",linestyle='-',color="green",label="Protocol 1: rotated S [[18,2,3]]")

ax.plot(tc0,S_Bell,linestyle='--',color="red", label="Unencoded Bell-State")



#Data 2,0.5


Vec_Y2b = [0.040194, 0.045753, 0.051874, 0.057108, 0.063011, 0.068477, 0.073475, 0.079524]
Vec_Y3b = [0.03874, 0.04324, 0.047897, 0.053522, 0.05795, 0.062747, 0.067691, 0.072158]

Std_m2b = [0.03946775, 0.04498015, 0.05105382, 0.056248690000000004, 0.06211142, 0.06754273, 0.07250912, 0.07852063000000001]
Std_p2b = [0.04092876, 0.04653425, 0.05270242, 0.057975160000000005, 0.06391836000000001, 0.06941932, 0.07444814, 0.08053531]

Std_m3b = [0.03802646, 0.04248824, 0.04710697, 0.05268898, 0.05708609, 0.06184875, 0.06675983, 0.07120021]
Std_p3b = [0.039462050000000005, 0.0440004, 0.0486951, 0.05436366, 0.05882228, 0.0636533, 0.06863002, 0.07312315]

S_Bell = [0.020246, 0.02401, 0.02752, 0.031128, 0.034473, 0.037886, 0.041768, 0.044256, 0.04799, 0.051445, 0.054105, 0.057557, 0.061082, 0.06436, 0.067638, 0.070205, 0.073154, 0.077131, 0.080449, 0.082937, 0.086368, 0.08924, 0.092351, 0.095107, 0.09832, 0.101788, 0.104559, 0.106838, 0.110193, 0.112973, 0.115944, 0.118507, 0.121474, 0.124221]

Vec_Y2 = [0.023094, 0.028983, 0.034661, 0.040384]
Std_p = [0.143414, 0.16845558, 0.19180938, 0.21434739000000003]
Std_p2 = [0.02365688, 0.02961163, 0.03534577, 0.041120330000000004]
Std_p3 = [0.02440557, 0.0291026, 0.034503160000000005, 0.03932075]

Std_m = [0.14081622, 0.16567783, 0.18889113000000002, 0.21130301999999998]
Std_m2 = [0.02254, 0.02836295, 0.03398475, 0.03965593]
Std_m3 = [0.02327109, 0.02786601, 0.03315933, 0.037887870000000004]

Vec_Y = [0.142112, 0.167063, 0.190347, 0.212822]
Vec_Y2 = [0.023094, 0.028983, 0.034661, 0.040384]
Vec_Y3 = [0.023834, 0.02848, 0.033827, 0.0386]



ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"x",linestyle='-',color="black",label="Protocol 2: BS [[18,2,3]]")


ax.fill_between(tc2b, Std_m2b, Std_p2b,alpha=0.2,color="blue")
ax.plot(tc2b,Vec_Y2b,"x",linestyle='-',color="blue")

ax.fill_between(tc3b, Std_m3b, Std_p3b,alpha=0.2,color="green")
ax.plot(tc3b,Vec_Y3b,"x",linestyle='-',color="green")


ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"x",linestyle='-',color="blue",label="Protocol 2: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"x",linestyle='-',color="green",label="Protocol 2: rotated S [[18,2,3]]")


second_legend = plt.legend("here", loc='lower right')


plt.xlim(0,tc1[-1])
plt.ylim(0,Vec_Y[-1])

ax.grid(linestyle='--')
plt.legend(fontsize = 9,loc="center left")

props = dict(facecolor='white', alpha=0.5)

ax.text(0.020, 0.975, "$2T_{1}, 2T_{2}, 0.5p_{err_{H}}, 0.5p_{err_{CX}}, 0.5p_{err_{M}}$", transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

# Manually add the first legend back to the plot

ax.set_xlabel("Storage time [$\mu$s]",fontsize=14)
ax.set_ylabel("Logical Error Rate",fontsize=14)
plt.show()

#######################################################





#DATA 2,0.1

Vec_Y2b = [0.01094, 0.012345, 0.013346, 0.014814, 0.016117, 0.017213, 0.018461, 0.019621]
Vec_Y3b = [0.009773, 0.010742, 0.011819, 0.012587, 0.013713, 0.014804, 0.016016, 0.016762]

Std_m2b = [0.010557879999999999, 0.01193895, 0.01292376, 0.014369360000000001, 0.01565332, 0.016734009999999997, 0.01796501, 0.01910999]
Std_p2b = [0.01133109, 0.01276006, 0.01377717, 0.015267530000000001, 0.01658966, 0.01770093, 0.01896586, 0.02014094]

Std_m3b = [0.00941184, 0.010363299999999999, 0.01142173, 0.01217695, 0.013285209999999999, 0.014359549999999999, 0.01555372, 0.01628924]
Std_p3b = [0.01014318, 0.01112975, 0.01222527, 0.013006, 0.01414974, 0.01525742, 0.01648722, 0.01724373]

S_Bell = [0.014182, 0.017927, 0.021433, 0.024872, 0.028463, 0.031982, 0.035652, 0.038933, 0.041991, 0.045678, 0.048689, 0.052718, 0.054976, 0.058641, 0.062191, 0.065007, 0.068438, 0.071435, 0.07445, 0.078084, 0.081462, 0.084259, 0.087034, 0.090151, 0.09295, 0.09652, 0.099532, 0.101889, 0.105397, 0.108027, 0.111001, 0.113845, 0.1159, 0.119147]

Vec_Y2 = [0.007304, 0.008651, 0.009813, 0.011052]
Std_p = [0.02244741, 0.02580503, 0.02907645, 0.03290525]
Std_p2 = [0.00762508, 0.00899978, 0.01018401, 0.011445120000000001]
Std_p3 = [0.00695563, 0.00811523, 0.00904875, 0.01011768]

Std_m = [0.021359490000000002, 0.024639619999999997, 0.02784015, 0.03159137]
Std_m2 = [0.0069920500000000005, 0.008311260000000001, 0.00945108, 0.010667870000000001]
Std_m3 = [0.0063514700000000006, 0.007461850000000001, 0.008358280000000001, 0.009387340000000001]

Vec_Y = [0.021899, 0.025218, 0.028454, 0.032244]
Vec_Y2 = [0.007304, 0.008651, 0.009813, 0.011052]
Vec_Y3 = [0.006649, 0.007784, 0.008699, 0.009748]


#Data 2,0.1
fig, ax = plt.subplots(1, 1)


ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"v",linestyle='-',color="black",label="Protocol 1: BS [[18,2,3]] ")


ax.fill_between(tc2b, Std_m2b, Std_p2b,alpha=0.2,color="blue")
ax.plot(tc2b,Vec_Y2b,"v",linestyle='-',color="blue")

ax.fill_between(tc3b, Std_m3b, Std_p3b,alpha=0.2,color="green")
ax.plot(tc3b,Vec_Y3b,"v",linestyle='-',color="green")


ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"v",linestyle='-',color="blue",label="Protocol 1: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"v",linestyle='-',color="green",label="Protocol 1: rotated S [[18,2,3]]")

ax.plot(tc0,S_Bell,linestyle='--',color="red", label="Unencoded Bell-State")




Vec_Y2b = [0.009381, 0.010536, 0.011808, 0.013085, 0.014344, 0.015486, 0.01688, 0.017928]
Vec_Y3b = [0.009163, 0.010321, 0.011486, 0.012394, 0.013253, 0.014217, 0.015147, 0.016255]

Std_m2b = [0.009027209999999999, 0.010161, 0.01141099, 0.01266711, 0.01390644, 0.01503141, 0.016405490000000002, 0.017439]
Std_p2b = [0.00974382, 0.01092, 0.012213950000000001, 0.01351186, 0.01479051, 0.01594953, 0.01736341, 0.01842592]

Std_m3b = [0.008813379999999999, 0.00994983, 0.01109438, 0.01198711, 0.01283241, 0.01378145, 0.01469745, 0.01578952]
Std_p3b = [0.00952171, 0.010701229999999999, 0.01188661, 0.01280992, 0.01368257, 0.01466151, 0.0156055, 0.01672941]


Vec_Y2 = [0.005777, 0.006906, 0.008121, 0.009364]
Std_p = [0.021345939999999997, 0.02443775, 0.028192389999999998, 0.03170459]
Std_p2 = [0.00606325, 0.00721831, 0.00845909, 0.00972648]
Std_p3 = [0.00658534, 0.00763008, 0.00866114, 0.009461639999999999]

Std_m = [0.020284919999999998, 0.02330288, 0.0269743, 0.03041407]
Std_m2 = [0.005499850000000001, 0.00660272, 0.0077919799999999996, 0.009010549999999999]
Std_m3 = [0.00599775, 0.00699699, 0.00798591, 0.00875541]

Vec_Y = [0.020811, 0.023866, 0.027579, 0.031055]
Vec_Y2 = [0.005777, 0.006906, 0.008121, 0.009364]
Vec_Y3 = [0.006287, 0.007309, 0.008319, 0.009104]


ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
ax.plot(tc1,Vec_Y,"x",linestyle='-',color="black",label="Protocol 2: BS [[18,2,3]]")


ax.fill_between(tc2b, Std_m2b, Std_p2b,alpha=0.2,color="blue")
ax.plot(tc2b,Vec_Y2b,"x",linestyle='-',color="blue")

ax.fill_between(tc3b, Std_m3b, Std_p3b,alpha=0.2,color="green")
ax.plot(tc3b,Vec_Y3b,"x",linestyle='-',color="green")


ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
ax.plot(tc2,Vec_Y2,"x",linestyle='-',color="blue",label="Protocol 2: S [[18,2,3]]")

ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
ax.plot(tc3,Vec_Y3,"x",linestyle='-',color="green",label="Protocol 2: rotated S [[18,2,3]]")


second_legend = plt.legend("here", loc='lower right')


plt.xlim(0,tc1[-1])
plt.ylim(0,Vec_Y[-1])

ax.grid(linestyle='--')
plt.legend(fontsize = 9)

props = dict(facecolor='white', alpha=0.5)

ax.text(0.020, 0.1, "$2T_{1}, 2T_{2}, 0.1p_{err_{H}}, 0.1p_{err_{CX}}, 0.1p_{err_{M}}$", transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

# Manually add the first legend back to the plot

ax.set_xlabel("Storage time [$\mu$s]",fontsize=14)
ax.set_ylabel("Logical Error Rate",fontsize=14)
plt.show()
##########################################################################################################################################################################






