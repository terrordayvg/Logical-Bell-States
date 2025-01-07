import itertools
import stim
import pytest
import sys 
import sinter 
from typing import List
import matplotlib.pyplot as plt

import sinter
import stim
from typing import List
import matplotlib.pyplot as plt
import numpy as np

from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.quantum_info import state_fidelity
import math
import sys
from typing import Callable, TypeVar, List, Any, Iterable, Optional, TYPE_CHECKING, Dict, Union, Literal, Tuple
from typing import cast

import numpy as np

from sinter._probability_util import fit_binomial, shot_error_rate_to_piece_error_rate, Fit
from qiskit import QuantumCircuit, transpile
#from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram
from qiskit_aer.noise import NoiseModel
from qiskit import transpile, assemble
#from qiskit_ibm_provider import IBMProvider


############### PROTOCOL 2 ###########################################
####Author: Vladlen G. updated: 7/1/2025
#### This code is responsible for figures 7,8,9 
#### - To obtain 7,9 export the data by varying the parameters below
#### - To obtain 8 erase the identation on lines of this code: 1030-1100 lines
"""
Parameters: 
dist: distance in meters, eg: 1000 = 1000 meters
l_p: loss ratio of the overall meaurements and gate errors
t_p: T1 and T2 ratio during the protocol procedure
iterations: Amount of QEC iterations whith is represented with increased storage time
max_shots: To increase the accuracy of the results increase the amount of montecarlo shots, currently at 10^6

Requirements:

-Update path in sys.path.insert(1, '...'), to include calibration2.csv which containts the IBM data for the gate noise.
-Have all imports up to date

"""
sys.path.insert(1, 'Include here the path ')



TCurveId = TypeVar('TCurveId')

MARKERS: str = "ov*sp^<>8PhH+xXDd|" * 100
T = TypeVar('T')
TVal = TypeVar('TVal')
TKey = TypeVar('TKey')
import gen


def s_after(s, delim):
    return s.partition(delim)[2]

def test(S):
    try:
        fs=float(S)
        return fs 
    except:
        return

def Noise_model2(idle_t,los_p,t_l):
  

    #backend = provider.get_backend("ibm_brisbane")
    #props = backend.properties()

    #Import from a file of ibm_brisbane 127 qubits
    import csv
    Errd=[]
    Amp_pha=[]

    CX_err=[]
    H_err=[]
    M_err=[]

    Loc_H=[]
    Loc_CNOT=[]

    Read_err=[]
    p_l=[]
    #Gate time
    #t_CX= 0.660 #micros  #660ns + buffer
    #t_H=  0.120  #micros  #120ns
    #t_M=  4     #micros   4000ns
    t_CX=0.533
    t_H=0.06
    t_M=1.22


    A=0.00654
    B=0.00486
    m=1
    p_l=((A/(A+B))*(1-np.exp(-(A+B)*m)))

    #backend = provider.get_backend('ibmq_16_melbourne')
    #print(backend.properties())
    cross_talk=1.3             #Based on mean from fig4. Software Mitigation of Crosstalk on Noisy Intermediate-Scale Quantum Computers - random circuits on processor need to be created to study the coupling error relation between cx

    #Import calibration data: Amplitude attenuation + dephasing, logical local error, Readout qubit error
    with open("calibration2.csv") as fp:
        reader = csv.reader(fp, delimiter=",", quotechar='"')
    # next(reader, None)  # skip the headers
        
        i=0
        for row in reader:
            data_read = row
            #print(data_read)
            if (i==1):
                #print(data_read)
                #calibration.csv
                #               Q,                      T1,                T2 (micros),    Readout error,      ID,                 Sx,                 X,                      ECR,                gate time (2 qubit ECR) data_read[13]
                #Errd.append([float(data_read[0]),float(data_read[1]),float(data_read[2]),float(data_read[5]),float(data_read[9]),float(data_read[10]),float(data_read[11]),float(data_read[12])])

                #calibration2.csv
                #               Q,                      T1,                T2 (micros),    Readout error,      ID,                 Z angle,              Sx,                    X,                      ECR,                gate time (2 qubit ECR) data_read[13]
                #Errd.append([float(data_read[0]),float(data_read[1]),float(data_read[2]),float(data_read[5]),float(data_read[9]),float(data_read[10]), float(data_read[11]),float(data_read[12]),float(data_read[13])])



                #print(Errd)

                #Damping + dephasing, pauli channel (px,py,pz)
                #CX_err.append([(1-np.exp(-t_CX/float(data_read[1])*t_l))/4,(1-np.exp(-t_CX/float(data_read[1])*t_l))/4,(1-np.exp(-t_CX/float(data_read[2])*t_l))/2 - (1-np.exp(-t_CX/float(data_read[1])*t_l))/4])
               # ([1-((1-np.exp(-t_CX/float(data_read[1])*t_l))/4 + (1-np.exp(-t_CX/float(data_read[1])*t_l))/4 + (1-np.exp(-t_CX/float(data_read[2])*t_l))/2 - (1-np.exp(-t_CX/float(data_read[1])*t_l))/4),(1-np.exp(-t_CX/float(data_read[1])*t_l))/4,(1-np.exp(-t_CX/float(data_read[1])*t_l))/4,(1-np.exp(-t_CX/float(data_read[2])*t_l))/2 - (1-np.exp(-t_CX/float(data_read[1])*t_l))/4])
                #H_err.append([(1-np.exp(-t_H/float(data_read[1])*t_l))/4,(1-np.exp(-t_H/float(data_read[1])*t_l))/4,(1-np.exp(-t_H/float(data_read[2])*t_l))/2 - (1-np.exp(-t_H/float(data_read[1])*t_l))/4])
                #M_err.append([(1-np.exp(-t_M/float(data_read[1])*t_l))/4,(1-np.exp(-t_M/float(data_read[1])*t_l))/4,(1-np.exp(-t_M/float(data_read[2])*t_l))/2 - (1-np.exp(-t_M/float(data_read[1])*t_l))/4])
                
                CX_err.append([
                1 - ((1 - np.exp(-t_CX / (float(data_read[1]) * t_l))) / 4 +
            (1 - np.exp(-t_CX / (float(data_read[1]) * t_l))) / 4 +
            (1 - np.exp(-t_CX / (float(data_read[2]) * t_l))) / 2 -
            (1 - np.exp(-t_CX / (float(data_read[1]) * t_l))) / 4),
        (1 - np.exp(-t_CX / (float(data_read[1]) * t_l))) / 4,
        (1 - np.exp(-t_CX / (float(data_read[1]) * t_l))) / 4,
        (1 - np.exp(-t_CX / (float(data_read[2]) * t_l))) / 2 - (1 - np.exp(-t_CX / (float(data_read[1]) * t_l))) / 4
                ])



               
               # CX_err.append([(1 - np.exp(-t_CX / (float(data_read[1]) * t_l))) / 4,
               # (1 - np.exp(-t_CX / (float(data_read[1]) * t_l))) / 4,
               # (1 - np.exp(-t_CX / (float(data_read[2]) * t_l))) / 2 - (1 - np.exp(-t_CX / (float(data_read[1]) * t_l))) / 4])



                H_err.append([(1 - np.exp(-t_H / (float(data_read[1]) * t_l))) / 4,
               (1 - np.exp(-t_H / (float(data_read[1]) * t_l))) / 4,
               (1 - np.exp(-t_H / (float(data_read[2]) * t_l))) / 2 - (1 - np.exp(-t_H / (float(data_read[1]) * t_l))) / 4])

                M_err.append([(1 - np.exp(-t_M / (float(data_read[1]) * t_l))) / 4,
               (1 - np.exp(-t_M / (float(data_read[1]) * t_l))) / 4,
               (1 - np.exp(-t_M / (float(data_read[2]) * t_l))) / 2 - (1 - np.exp(-t_M / (float(data_read[1]) * t_l))) / 4])


                #Local gate error 

                Loc_H.append([float(data_read[11])*los_p]) #10

                

                S=s_after(data_read[13], ":") #12
                if test(S)!=None:
                    Loc_CNOT.append(test(S)*los_p)
                     
                #Readout error
                Read_err.append(float(data_read[5])*los_p)

                #print(CX_err)
                #print(H_err)
                #Amp_pha.append([])
            i=1

    ##########################################################################
    #Return vectors for errors:
    return CX_err,H_err,Loc_H,Loc_CNOT,Read_err,M_err,p_l



def Noise_model(idle_t,los_p,t_l):
  



    #Import from a file of ibm_brisbane 127 qubits
    import csv
    Errd=[]
    Amp_pha=[]

    ID_err=[]
    H_err=[]
    M_err=[]

    Loc_H=[]
    Loc_CNOT=[]

    Read_err=[]
    p_l=[]
    #Gate time
    t_CX= 0.660 #micros  #660ns + buffer
    t_H=  0.120 #micros  #120ns
    t_M=  4     #micros   4000ns

    A=0.00654
    B=0.00486
    m=1
    p_l=((A/(A+B))*(1-np.exp(-(A+B)*m)))

    #backend = provider.get_backend('ibmq_16_melbourne')
    cross_talk=1.3             #Based on mean from fig4. Software Mitigation of Crosstalk on Noisy Intermediate-Scale Quantum Computers - random circuits on processor need to be created to study the coupling error relation between cx

    #Import calibration data: Amplitude attenuation + dephasing, logical local error, Readout qubit error
    with open("calibration2.csv") as fp:
        reader = csv.reader(fp, delimiter=",", quotechar='"')
    # next(reader, None)  # skip the headers
        
        i=0
        for row in reader:
            data_read = row
            #print(data_read)
            if (i==1):
                #               Q,                      T1,                T2 (micros),    Readout error,      ID,                 Sx,                 X,                      ECR,                gate time (2 qubit ECR) data_read[13]
                #Errd.append([float(data_read[0]),float(data_read[1]),float(data_read[2]),float(data_read[5]),float(data_read[9]),float(data_read[10]),float(data_read[11]),float(data_read[12])])

                #print(Errd)

                #dephasing decoherence on waiting time - pauli channel (px,py,pz)
                ID_err.append([
    (1 - np.exp(-idle_t / (float(data_read[1]) * t_l))) / 4,
    (1 - np.exp(-idle_t / (float(data_read[1]) * t_l))) / 4,
    (1 - np.exp(-idle_t / (float(data_read[2]) * t_l))) / 2 - (1 - np.exp(-idle_t / (float(data_read[1]) * t_l))) / 4
])                #Local#append([1-((1-np.exp(-idle_t/float(data_read[1]*t_l)))/4 + (1-np.exp(-idle_t/float(data_read[1]*t_l)))/4 + (1-np.exp(-idle_t/float(data_read[2]*t_l)))/2 - (1-np.exp(-idle_t/float(data_read[1]*t_l)))/4),(1-np.exp(-t_CX/float(data_read[1]*t_l)))/4,(1-np.exp(-t_CX/float(data_read[1]*t_l)))/4,(1-np.exp(-t_CX/float(data_read[2]*t_l)))/2 - (1-np.exp(-t_CX/float(data_read[1]*t_l)))/4])
               
             
            i=1


    #Amplitude + dephasing ----------------------------------------------------------------------------------------------------------------------------------------------------------------------



    #IBMProvider.save_account("b0b806108cf6b0d6561dd1cf5626b521536d0f054c294143cf6778620042106be73424c02a31898269f75dcfe5db0dd72f11f23d3ce4a590ec5ab246741e9349",overwrite=True) #POL
    #IBMQ.enable_account("")

    #provider = IBMProvider("31ce3dfece14c5c94b26f5fd59b804d2ade3ba778a1a1926f4848b65311a52ce79d8f6be87c243b18c06930231fd6fe399f4ab3939fa1e0b158aa7074f465e5a")
    #print(provider)
    #backend=provider.get_backend('ibm_brisbane')#

# Build noise model from backend properties
   # backend = FakeVigo()
   
    #print(noise_model.__dict__.keys())
    #print("basis_gates")
    #print(noise_model._basis_gates)
    #print("noise_intructions")
    #print(noise_model._noise_instructions)
    #print("noise_qubits")

    #print(noise_model._noise_qubits)
    #print("default_quantum_errors")

    #print(noise_model._default_quantum_errors)
    #print("local_quantum_errros")

    #print(noise_model._local_quantum_errors)
    #print("nonlocal_quantum_errors")

    #print(noise_model._nonlocal_quantum_errors)
    #print("default_redout_error")

    #print(noise_model._default_readout_error)
    #print("local_redout_errors")

    #print(noise_model._local_readout_errors[(1,)])
    #print("custom_noise_passes")

    #print(noise_model._custom_noise_passes)
    #print(" ")

    ##########################################################################
    #Return vectors for errors:
    return ID_err

#30nanoseconds both qubits
def Scheduler_source(j,dist,T_p):
    #Scheduler source has frequency of generation time of the circuit.
    #Sends photon-photon with that frequency frequency 720ns + poisson distribution for the frequency.
    #Time of travel c/
    #Sends driving pulse storage -

    tH=30*pow(10,-9)
    #Initialization of the photon
    T_p[0]=tH
    T_p[1]=tH
    #T_p[0]=target T_p[1]=control, [0] delay travel, [1] storage normal

    #delay time
    return  T_p


from scipy.constants import c as clight
def Scheduler_memory(j,dist,T_p,shots):
    #Fiber optic delay time of photon depending on the distance.
    c=clight
    n=1.44 #refractive index in glass
    vel=c/n 
    ttrav=dist/vel #photon delay time

    ##### Driving pulse for storage
    #Lets consider 60ns for storage pulse
    #Lets consider 60ns for readout pulse

    #Position in memory in terms of storage times
    A=30*pow(10,-9)          #interval betweeen sending (frequency of the source 30nanoseconds)
    tpul_store=30*pow(10,-9) #storage waiting time to input states

    T_p[0]=(j)*(A+tpul_store)+ttrav #Time of storage during the control
    T_p[1]=(j)*(A+tpul_store) #Time of storage during the target
    #print("Signal_in_storage times: Control - Target: "+str(shots-j))
    #print(T_p[j][0],T_p[j][1])

    return T_p

#Additional photon time to wait to be sent back in storage
def Scheduler_returnS(j,dist,T_p,shots):
    #Here it gets tricky, you need to return the signal back to from Bob->Alice, the mem qubits in Alice will decohere more, but, we optimize by sending first the "longest" stored qubits to avoid decoherence.
    c=clight #speed of light in vacum
    n=1.44 #refractive index in glass
    vel=c/n 
    ttrav=dist/vel

    tpul_read=30*pow(10,-9) #recover time 60ns
    T_p[0]=(shots-j-1)*(tpul_read)+ttrav+T_p[0] #Control
    T_p[1]=(shots-j-1)*(tpul_read)+T_p[1] #Target
    #print("Signal_out_storage times: Control - Target: "+str(shots-j))
    #print(T_p[j][0],T_p[j][1])

    return T_p

def Scheduler_Ver(j,dist,T_p,shots):

    tH=30*pow(10,-9)
    tCX=30*pow(10,-9)

    #Initialization of the photon
    T_p[0]=T_p[0]+tH+tCX
    T_p[1]=T_p[0]+tH+tCX
    #print("Signal_verification times: Control - Target: "+str(shots-j))
    #print(T_p[j][0],T_p[j][1])

    return T_p



def Fiber_loss2(dist,lamb,j):
    #We implement this by losing the qubits
    #Measuring the state and disapearing from a memory, how do you inform Alice about that, if this information is secret
    #dist = [km]
    #lamb = [db/km]
    #p_resonant_p=1-0.06           #probability of resonance if happens
    p_transm=pow(10,-dist*lamb/10) # in fiber optic case 
    #p_fconversion=0.43 #Lets assume a convertion of 43% of photons from optical-microwave and vice-versa from article 20 (fit) - optimal case with conversion 1588nm-637nm (assume 1538nm)
    
    ############### detectoion
    p_detection=0.90  #Assume fiber optic detector probability to find your photon #best parameters possible 0.9 here
    p_dark_counts=1-(1-np.exp(-0.00000025))

    p_channel=p_transm*p_detection*p_dark_counts  # A ------- B -------- A  2-3-4-5
    random.seed(66122+j)
    s = np.random.uniform(0,1,1)
    #Measure if it was lost, and create a tag "NM=Not measured" during the time window
    if (s>=p_channel):
        return 1
    return 0




import random


def Rate(store_t,lamb):
    data=20000 #1 KHZ =blocks of [16]+[16] qubits
    Count=0
    Count_T=[]
    Count_tot=0
    for dist in range(1000,10000,5000):
        for j in range(data):
            for k in range(9):
                Count=Fiber_loss(dist*pow(10,-3),lamb,k,j) + Read_loss(store_t*pow(10,-9),store_t*pow(10,-9),k,j)+Count

            if (Count == 0):
                Count_tot=Count_tot+1
            Count=0


        Count_T.append(10*np.log10(Count_tot/data))
        Count_tot=0


    Count2=0
    Count_T2=[]
    Count_tot2=0
    for dist in range(1000,10000,5000):
        print(dist)
        for j in range(data):
            for k in range(1):
                Count2=Fiber_loss(dist*pow(10,-3),lamb,k,j+24124) + Read_loss(store_t*pow(10,-9),store_t*pow(10,-9),k,j+33341)+Count2

            if (Count2 == 0):
                Count_tot2=Count_tot2+1
            Count2=0


        Count_T2.append(10*np.log10(Count_tot2/data))
        Count_tot2=0

    print(Count_T,Count_T2)
    return Count_T,Count_T2

def Fiber_loss(dist,lamb,k,j):
    #We implement this by losing the qubits
    #Measuring the state and disapearing from a memory, how do you inform Alice about that, if this information is secret
    #dist = [km]
    #lamb = [db/km]
    #p_resonant_p=1-0.06           #probability of resonance if happens
    p_transm=pow(10,-dist*lamb/10) # in fiber optic case 
    #p_fconversion=0.43 #Lets assume a convertion of 43% of photons from optical-microwave and vice-versa from article 20 (fit) - optimal case with conversion 1588nm-637nm (assume 1538nm)
    
    ############### detectoion
    p_detection=0.90  #Assume fiber optic detector probability to find your photon #best parameters possible 0.9 here
    p_dark_counts=1-(1-np.exp(-0.00000025))
    p_channel=p_transm*p_detection*p_dark_counts  # A ------- B -------- A  2-3-4-5
    

    random.seed(1014+j*k)
    s = np.random.uniform(0,1,1)

    
    #print(p_channel)
    #Measure if it was lost, and create a tag "NM=Not measured" during the time window
    if (s>=p_channel):
        return 1

    random.seed(152155+j*k)
    s2 = np.random.uniform(0,1,1)
    #Measure if it was lost, and create a tag "NM=Not measured" during the time window
    if (s2>=p_channel):
        return 1

    return 0




def Read_loss(T1a,T1b,k,j):

    R1=0.96
    R2=0.999
    alpha=1                                         #Alpha*L=Alpha
    lamb=5000                                       #comb FWHM
    F=40                                            #comb finesse
    lambd=(2*np.pi*lamb)/(np.sqrt(8*np.log(2)))
    alpha_d=((alpha)/F)*np.sqrt(np.pi/(4*np.log(2)))

    #t=np.arange(0,0.001,0.0001)
    t=T1a #100mus
    nu=(4*(alpha_d)*(alpha_d)*np.exp(-2*(alpha_d))*(1-R1)*(1-R1)*R2*np.exp(-t*t*lambd*lambd))/pow((1-np.sqrt(R1*R2)*np.exp(-alpha_d)),4)


    t2=T1b #100mus
    nu2=(4*(alpha_d)*(alpha_d)*np.exp(-2*(alpha_d))*(1-R1)*(1-R1)*R2*np.exp(-t2*t2*lambd*lambd))/pow((1-np.sqrt(R1*R2)*np.exp(-alpha_d)),4)
    #We implement this by losing the qubits
    #Measuring the state and disapearing from a memory, how do you inform Alice about that, if this information is secret
    #dist = [km]
    #lamb = [db/km]
    #p_transm=pow(10,-dist*lamb/10) # in fiber optic case with conversion 1588nm-637nm (assume 1538nm)
    #p_fconversion=0.43 #Lets assume a convertion of 43% of photons from optical-microwave and vice-versa from article 20 (fit) - optimal case
    
    ############### detectoion
    #p_detection=0.85  #Assume fiber optic detector probability to find your photon #best parameters possible 0.9 here


    p_channel=nu # A ------- B -------- A  2-3-4-5 ------ if p_resonant double excitation happens, we abort the photon
    p_channel2=nu2 # A ------- B -------- A  2-3-4-5 ------ if p_resonant double excitation happens, we abort the photon

    #print(p_channel,p_channel2)
    random.seed(10116+j*k)
    s = np.random.uniform(0,1,1)    #Target qubit measurement
    if (s>=p_channel):
        return 1

    random.seed(1052266+j*k)
    s2 = np.random.uniform(0,1,1)
    if (s2>=p_channel2):
        return 1 

    return 0



def better_sorted_str_terms(val: Any) -> Any:
    """A function that orders "a10000" after "a9", instead of before.

    Normally, sorting strings sorts them lexicographically, treating numbers so
    that "1999999" ends up being less than "2". This method splits the string
    into a tuple of text pairs and parsed number parts, so that sorting by this
    key puts "2" before "1999999".

    Because this method is intended for use in plotting, where it's more
    important to see a bad result than to see nothing, it returns a type that
    tries to be comparable to everything.

    Args:
        val: The value to convert into a value with a better sorting order.

    Returns:
        A custom type of object with a better sorting order.

    Examples:
        >>> import sinter
        >>> items = [
        ...    "distance=199999, rounds=3",
        ...    "distance=2, rounds=3",
        ...    "distance=199999, rounds=199999",
        ...    "distance=2, rounds=199999",
        ... ]
        >>> for e in sorted(items, key=sinter.better_sorted_str_terms):
        ...    print(e)
        distance=2, rounds=3
        distance=2, rounds=199999
        distance=199999, rounds=3
        distance=199999, rounds=199999
    """

    if isinstance(val, tuple):
        return tuple(better_sorted_str_terms(e) for e in val)
    if not isinstance(val, str):
        return val
    terms = split_by(val, lambda c: c in '.0123456789')
    result = []
    for term in terms:
        term = ''.join(term)
        if '.' in term:
            try:
                term = float(term)
            except ValueError:
                try:
                    term = tuple(int(e) for e in term.split('.'))
                except ValueError:
                    pass
        else:
            try:
                term = int(term)
            except ValueError:
                pass
        result.append(term)
    return tuple(LooseCompare(e) for e in result)

def group_by(items: Iterable[TVal],
             *,
             key: Callable[[TVal], TKey],
             ) -> Dict[TKey, List[TVal]]:
    """Groups items based on whether they produce the same key from a function.

    Args:
        items: The items to group.
        key: Items that produce the same value from this function get grouped together.

    Returns:
        A dictionary mapping outputs that were produced by the grouping function to
        the list of items that produced that output.

    Examples:
        >>> import sinter
        >>> sinter.group_by([1, 2, 3], key=lambda i: i == 2)
        {False: [1, 3], True: [2]}

        >>> sinter.group_by(range(10), key=lambda i: i % 3)
        {0: [0, 3, 6, 9], 1: [1, 4, 7], 2: [2, 5, 8]}
    """

    result: Dict[TKey, List[TVal]] = {}

    for item in items:
        curve_id = key(item)
        result.setdefault(curve_id, []).append(item)

    return result

def plot_custom(
        *,
        stats: 'Iterable[sinter.TaskStats]',
        x_func: Callable[['sinter.TaskStats'], Any],
        y_func: Callable[['sinter.TaskStats'], Union['sinter.Fit', float, int]],
        group_func: Callable[['sinter.TaskStats'], TCurveId] = lambda _: None,
        filter_func: Callable[['sinter.TaskStats'], Any] = lambda _: True,
        plot_args_func: Callable[[int, TCurveId, List['sinter.TaskStats']], Dict[str, Any]] = lambda index, group_key, group_stats: dict(),
        line_fits: Optional[Tuple[Literal['linear', 'log', 'sqrt'], Literal['linear', 'log', 'sqrt']]] = None,
) -> None:
    """Plots error rates in curves with uncertainty highlights.

    Args:
        ax: The plt.Axes to plot onto. For example, the `ax` value from `fig, ax = plt.subplots(1, 1)`.
        stats: The collected statistics to plot.
        x_func: The X coordinate to use for each stat's data point. For example, this could be
            `x_func=lambda stat: stat.json_metadata['physical_error_rate']`.
        y_func: The Y value to use for each stat's data point. This can be a float or it can be a
            sinter.Fit value, in which case the curve will follow the fit.best value and a
            highlighted area will be shown from fit.low to fit.high.
        group_func: Optional. When specified, multiple curves will be plotted instead of one curve.
            The statistics are grouped into curves based on whether or not they get the same result
            out of this function. For example, this could be `group_func=lambda stat: stat.decoder`.
        filter_func: Optional. When specified, some curves will not be plotted.
            The statistics are filtered and only plotted if filter_func(stat) returns True.
            For example, `filter_func=lambda s: s.json_metadata['basis'] == 'x'` would plot only stats
            where the saved metadata indicates the basis was 'x'.
        plot_args_func: Optional. Specifies additional arguments to give the the underlying calls to
            `plot` and `fill_between` used to do the actual plotting. For example, this can be used
            to specify markers and colors. Takes the index of the curve in sorted order and also a
            curve_id (these will be 0 and None respectively if group_func is not specified). For example,
            this could be:

                plot_args_func=lambda index, group_key, group_stats: {
                    'color': (
                        'red'
                        if group_key == 'decoder=pymatching p=0.001'
                        else 'blue'
                    ),
                }
        line_fits: Defaults to None. Set this to a tuple (x_scale, y_scale) to include a dashed line
            fit to every curve. The scales determine how to transform the coordinates before
            performing the fit, and can be set to 'linear', 'sqrt', or 'log'.
    """

    # Backwards compatibility to when the group stats argument wasn't present.
    import inspect
    if len(inspect.signature(plot_args_func).parameters) == 2:
        old_plot_args_func = cast(Callable[[int, TCurveId], Any], plot_args_func)
        plot_args_func = lambda a, b, _: old_plot_args_func(a, b)

    filtered_stats: List['sinter.TaskStats'] = [
        stat
        for stat in stats
        if filter_func(stat)
    ]

    curve_groups = group_by(filtered_stats, key=group_func)
    for k, curve_id in enumerate(sorted(curve_groups.keys(), key=better_sorted_str_terms)):
        this_group_stats = sorted(curve_groups[curve_id], key=x_func)

        xs = []
        ys = []
        xs_range = []
        ys_low = []
        ys_high = []
        saw_fit = False
        for stat in this_group_stats:
            num_kept = stat.shots - stat.discards
            if num_kept == 0:
                continue
            x = float(x_func(stat))
            y = y_func(stat)
            if isinstance(y, Fit):
                xs_range.append(x)
                ys_low.append(y.low)
                ys_high.append(y.high)
                saw_fit = True
                y = y.best
            if not math.isnan(y):
                xs.append(x)
                ys.append(y)
        curve_id="BC code d=4"
        kwargs: Dict[str, Any] = dict(plot_args_func(k, curve_id, this_group_stats))
        kwargs.setdefault('marker', MARKERS[k])
        if curve_id is not None:
            kwargs.setdefault('label', "BC code d=4")
            kwargs.setdefault('color', f'C{k}')
        kwargs.setdefault('color', 'black')
        #ax.plot(xs, ys, **kwargs)

        print("-------------------------------------")
        print("Real output:")
        print(xs,ys)

        if line_fits is not None and len(set(xs)) >= 2:
            x_scale, y_scale = line_fits
            fit_xs = _rescale(xs, x_scale, False)
            fit_ys = _rescale(ys, y_scale, False)

            from scipy.stats import linregress
            line_fit = linregress(fit_xs, fit_ys)

            x0 = fit_xs[0]
            x1 = fit_xs[-1]
            dx = x1 - x0
            x0 -= dx*10
            x1 += dx*10
            if x0 < 0 <= fit_xs[0] > x0 and x_scale == 'sqrt':
                x0 = 0

            out_xs = np.linspace(x0, x1, 1000)
            out_ys = out_xs * line_fit.slope + line_fit.intercept
            out_xs = _rescale(out_xs, x_scale, True)
            out_ys = _rescale(out_ys, y_scale, True)

            line_kwargs = kwargs.copy()
            line_kwargs.pop('marker', None)
            line_kwargs.pop('label', "BC code d=4")
            line_kwargs['linestyle'] = '--'
            line_kwargs.setdefault('linewidth', 1)
            line_kwargs['linewidth'] /= 2
            ax.plot(out_xs, out_ys, **line_kwargs)

        if saw_fit:
            fit_kwargs = kwargs.copy()
            fit_kwargs.setdefault('zorder', 0)
            fit_kwargs.setdefault('alpha', 1)
            fit_kwargs['zorder'] -= 100
            fit_kwargs['alpha'] *= 0.25
            fit_kwargs.pop('marker', None)
            fit_kwargs.pop('linestyle', None)
            fit_kwargs.pop('label', "BC code d=4")
            #ax.fill_between(xs_range, ys_low, ys_high, **fit_kwargs)

            print("-------------------------------------")
            print("Area ranges")
            print(xs_range,ys_low,ys_high)
            return xs, ys, ys_low, ys_high

def plot_error_rate2(
        *,
        stats: 'Iterable[sinter.TaskStats]',
        x_func: Callable[['sinter.TaskStats'], Any],
        failure_units_per_shot_func: Callable[['sinter.TaskStats'], Any] = lambda _: 1,
        failure_values_func: Callable[['sinter.TaskStats'], Any] = lambda _: 1,
        group_func: Callable[['sinter.TaskStats'], TCurveId] = lambda _: None,
        filter_func: Callable[['sinter.TaskStats'], Any] = lambda _: True,
        plot_args_func: Callable[[int, TCurveId, List['sinter.TaskStats']], Dict[str, Any]] = lambda index, group_key, group_stats: dict(),
        highlight_max_likelihood_factor: Optional[float] = 1e3,
        line_fits: Optional[Tuple[Literal['linear', 'log', 'sqrt'], Literal['linear', 'log', 'sqrt']]] = None,
) -> None:
    """Plots error rates in curves with uncertainty highlights.

    Args:
        ax: The plt.Axes to plot onto. For example, the `ax` value from `fig, ax = plt.subplots(1, 1)`.
        stats: The collected statistics to plot.
        x_func: The X coordinate to use for each stat's data point. For example, this could be
            `x_func=lambda stat: stat.json_metadata['physical_error_rate']`.
        failure_units_per_shot_func: How many error chances there are per shot. This rescales what the
            logical error rate means. By default, it is the logical error rate per shot, but this allows
            you to instead make it the logical error rate per round. For example, if the metadata
            associated with a shot has a field 'r' which is the number of rounds, then this can be
            achieved with `failure_units_per_shot_func=lambda stats: stats.metadata['r']`.
        failure_values_func: How many independent ways there are for a shot to fail, such as
            the number of independent observables in a memory experiment. This affects how the failure
            units rescaling plays out (e.g. with 1 independent failure the "center" of the conversion
            is at 50% whereas for 2 independent failures the "center" is at 75%).
        group_func: Optional. When specified, multiple curves will be plotted instead of one curve.
            The statistics are grouped into curves based on whether or not they get the same result
            out of this function. For example, this could be `group_func=lambda stat: stat.decoder`.
        filter_func: Optional. When specified, some curves will not be plotted.
            The statistics are filtered and only plotted if filter_func(stat) returns True.
            For example, `filter_func=lambda s: s.json_metadata['basis'] == 'x'` would plot only stats
            where the saved metadata indicates the basis was 'x'.
        plot_args_func: Optional. Specifies additional arguments to give the the underlying calls to
            `plot` and `fill_between` used to do the actual plotting. For example, this can be used
            to specify markers and colors. Takes the index of the curve in sorted order and also a
            curve_id (these will be 0 and None respectively if group_func is not specified). For example,
            this could be:

                plot_args_func=lambda index, curve_id: {'color': 'red'
                                                        if curve_id == 'pymatching'
                                                        else 'blue'}

        highlight_max_likelihood_factor: Controls how wide the uncertainty highlight region around curves is.
            Must be 1 or larger. Hypothesis probabilities at most that many times as unlikely as the max likelihood
            hypothesis will be highlighted.
        line_fits: Defaults to None. Set this to a tuple (x_scale, y_scale) to include a dashed line
            fit to every curve. The scales determine how to transform the coordinates before
            performing the fit, and can be set to 'linear', 'sqrt', or 'log'.
    """
    if highlight_max_likelihood_factor is None:
        highlight_max_likelihood_factor = 1
    if not (highlight_max_likelihood_factor >= 1):
        raise ValueError(f"not (highlight_max_likelihood_factor={highlight_max_likelihood_factor} >= 1)")

    def y_func(stat: 'sinter.TaskStats') -> Union[float, 'sinter.Fit']:
        result = fit_binomial(
            num_shots=stat.shots - stat.discards,
            num_hits=stat.errors,
            max_likelihood_factor=highlight_max_likelihood_factor,
        )

        pieces = failure_units_per_shot_func(stat)
        values = failure_values_func(stat)
        result = Fit(
            low=shot_error_rate_to_piece_error_rate(result.low, pieces=pieces, values=values),
            best=shot_error_rate_to_piece_error_rate(result.best, pieces=pieces, values=values),
            high=shot_error_rate_to_piece_error_rate(result.high, pieces=pieces, values=values),
        )

        if stat.errors == 0:
            result = Fit(low=result.low, high=result.high, best=float('nan'))

        if highlight_max_likelihood_factor == 1:
            return result.best
        return result

    print("Dataset:============================================")
    print(x_func)
    print(y_func)

    X,Y,Y_l,Y_h=plot_custom(
        stats=stats,
        x_func=x_func,
        y_func=y_func,
        group_func=group_func,
        filter_func=filter_func,
        plot_args_func=plot_args_func,
        line_fits=line_fits,
    )

    #print(X,Y,Y_l,Y_h)
    return X,Y,Y_l,Y_h


@pytest.mark.parametrize('width,height,basis,rounds', itertools.product(
    [4, 6],
    [6, 12],
    ['X', 'Z'],
    [1, 5, 6],
))
def test_make_bacon_shor_xx_lattice_surgery_circuit(width: int, height: int, basis: str, rounds: int):
    circuit = make_bacon_shor_xx_lattice_surgery_circuit(
        width=width,
        height=height,
        basis=basis,
        rounds=rounds,
    )
    circuit = gen.NoiseModel.uniform_depolarizing(1e-3).noisy_circuit(circuit)
    circuit.detector_error_model()

    expected_determined = circuit.num_detectors + circuit.num_observables
    assert gen.count_determined_measurements_in_circuit(circuit) == expected_determined

    assert circuit.num_ticks == (rounds + 2) * 4 + 1

    expected_distance = min(width // 2, rounds) if basis == 'X' else height // 2
    actual_distance = len(circuit.shortest_graphlike_error())
    assert actual_distance == expected_distance

    assert gen.gates_used_by_circuit(circuit) <= {
        'R',
        'M',
        'RX',
        'MX',
        'MXX',
        'MZZ',

        'TICK',
        'DEPOLARIZE1',
        'DEPOLARIZE2',
        'DETECTOR',
        'X_ERROR',
        'Z_ERROR',
        'OBSERVABLE_INCLUDE',
        'QUBIT_COORDS',
        'SHIFT_COORDS',
    }


def test_gate_counts():
    circuit = make_bacon_shor_xx_lattice_surgery_circuit(
        width=5,
        height=5,
        basis='X',
        rounds=3,
    )
    no_noise = gen.NoiseRule(after={})
    circuit = gen.NoiseModel(
        idle_depolarization=1e-3,
        gate_rules={
            'RX': 0.01,
            'MX': 0.01,
            'MPP': 0.01,
        }
    ).noisy_circuit(circuit)

    assert gen.gate_counts_for_circuit(circuit) == {
        'MXX': 215,
        'MZZ': 200,
        'DEPOLARIZE1': 170,  # Counts idles.

        'RX': 50,
        'MX': 50,

        'DETECTOR': 102,
        'QUBIT_COORDS': 50,
        'SHIFT_COORDS': 5,
        'OBSERVABLE_INCLUDE': 3,
        'TICK': 21,
    }



#print(circuit)
from scipy.constants import c as clight


def count_determined_measurements_in_circuit(circuit: stim.Circuit) -> int:
    """Simulates the circuit, counting how many measurements were determined.

    In most cases, for a quantum error correcting code, the result should be
    related to the number of detectors plus the number of observables declared
    in the circuit.
    """

    num_determined_measurements = 0
    sim = stim.TableauSimulator()
    n = circuit.num_qubits

    def run_block(block: stim.Circuit, reps: int):
        nonlocal num_determined_measurements
        for _ in range(reps):
            for inst in block:
                if isinstance(inst, stim.CircuitRepeatBlock):
                    run_block(inst.body_copy(), inst.repeat_count)
                elif inst.name == 'M' or inst.name == 'MR':
                    args = inst.gate_args_copy()
                    for t in inst.targets_copy():
                        assert t.is_qubit_target
                        known = sim.peek_z(t.value) != 0
                        num_determined_measurements += known
                        sim.do(stim.CircuitInstruction(inst.name, [t.value], args))
                elif inst.name == 'MX' or inst.name == 'MRX':
                    args = inst.gate_args_copy()
                    for t in inst.targets_copy():
                        assert t.is_qubit_target
                        known = sim.peek_x(t.value) != 0
                        num_determined_measurements += known
                        sim.do(stim.CircuitInstruction(inst.name, [t.value], args))
                elif inst.name == 'MY' or inst.name == 'MRY':
                    args = inst.gate_args_copy()
                    for t in inst.targets_copy():
                        assert t.is_qubit_target
                        known = sim.peek_y(t.value) != 0
                        num_determined_measurements += known
                        sim.do(stim.CircuitInstruction(inst.name, [t.value], args))
                elif inst.name == 'MPP':
                    args = inst.gate_args_copy()
                    targets = inst.targets_copy()
                    start = 0
                    while start < len(targets):
                        end = start + 1
                        while end < len(targets) and targets[end].is_combiner:
                            end += 2

                        p = stim.PauliString(n)
                        for t in targets[start:end:2]:
                            if t.is_x_target:
                                p[t.value] = 'X'
                            elif t.is_y_target:
                                p[t.value] = 'Y'
                            elif t.is_z_target:
                                p[t.value] = 'Z'
                            else:
                                raise NotImplementedError(f'{t=} {inst=}')

                        known = sim.peek_observable_expectation(p) != 0
                        num_determined_measurements += known
                        sim.do(stim.CircuitInstruction(inst.name, targets[start:end], args))

                        start = end
                else:
                    sim.do(inst)

    run_block(circuit, 1)
    return num_determined_measurements


if __name__ == '__main__':
    dist=1000
    c=clight
    n=1.44 #refractive index in glass
    vel=c/n 
    print(dist/vel)
    ttrav=(dist/vel)*1000000 #photon delay time
    lambd=0.17



    print("Traveling time:")
    print(ttrav)
    #Storage time decoherence and dephasing error
    l_p=0.5 #cut down the loss by 50%
    t_p=1 #T1 and T2 parameter 
    ID_err=Noise_model(ttrav,l_p,t_p)
    CX_err,H_err,Loc_H,Loc_CNOT,Read_err,M_err,p_l=Noise_model2(ttrav,l_p,t_p)

    tc0=[]
    tc1=[]
    tc2=[]
    tc3=[]

    tc0b=[]
    tc1b=[]
    tc2b=[]
    tc3b=[]

    tc0c=[]
    tc1c=[]
    tc2c=[]
    tc3c=[]

    tH=0.060
    tCX=0.533
    tM=1.216

    iterations=5
    Xl=np.arange(1,iterations,1)
    #Fiber optic errorf
    #fig = plt.figure(figsize=(10,10))
    #ax = fig.add_subplot(111)

    #d=5

    #noise=0.001
   # circuit=stim.Circuit.generated(
   #         "surface_code:unrotated_memory_x",
   #         rounds=3,
   #         distance=3,
   #         after_clifford_depolarization=noise,
   #         after_reset_flip_probability=noise,
   #         before_measure_flip_probability=noise,
   #         before_round_data_depolarization=noise,
   #     )
   # print(" ")
   # print(" ----------------------------------------------- ")
   # print(circuit.diagram)
    #rect = [0.35,0.1,0.6,0.6]
    ##for store_t in range(6760,6700,6760):
     #   Count,Count2,Count3=Rate(store_t,lambd)
     #   print(Count)
        
     #   ax.plot(Xl,Count,label='W_corr|t='+str(store_t/1000)+r'$\mu$s')
     #   ax.plot(Xl,Count2,"--",label='No_corr|t='+str(store_t/1000)+r'$\mu$s')
     #   ax.plot(Xl,Count2,"--",label='No_corr|t='+str(store_t/1000)+r'$\mu$s')

     #   ax.set_xlabel("Distance [km]")
     #   ax.set_ylabel("Loss of raw rate [dB]")
     #   ax.grid(which='major')
        #axa.plot(X2,Y2,"--", color="black",label='Google_D3')
     #   ax.legend(loc ="upper left")
    
   # plt.show()

    #for dist in range(1000,10000,1000):
    
    #    Count=Rate(dist,lambd)
    #    print(Count)
        
    #rect = [0.35,0.1,0.6,0.6]
    #    ax.plot(Xl,Count,label='dist='+str(dist/1000)+' km')
    #    ax.set_xlabel(r"Storage time [$\mu$s]")
    #    ax.set_ylabel("Loss of raw rate [dB]")
    #    ax.grid(which='major')
    #axa.plot(X2,Y2,"--", color="black",label='Google_D3')
    #    ax.legend(loc ="upper left")
    
    #plt.show()




    init=1
      #Variables for True X and Y values vectors and their standard deviations in upper and lower bounds
    Vec_Y=[]
    Vec_X=np.arange(init,iterations)        #Iteration vector
    Std_p=[]
    Std_m=[]

    Vec_Y2=[]
    Vec_X2=np.arange(init,iterations)        #Iteration vector
    Std_p2=[]
    Std_m2=[]

    Vec_Y3=[]
    Vec_X3=np.arange(init,iterations)        #Iteration vector
    Std_p3=[]
    Std_m3=[]

    S_Bell=[]
    S_Bell2=[]

    Ad=[]
    Ad2=[]

    for j in range(init,iterations):
        wait=(0.533*2+0.06*2+1.216)*j
        t0c=wait+ttrav+tH+tCX+tM

        Ad.append(t0c)

    for j in range(iterations-1,35):
        wait=(0.533*2+0.06*2+1.216)*j
        
        t0c=wait+ttrav+tH+tCX+tM

        Ad2.append(t0c)

    print("Processing time calculation for every code in protocol 2:----------------")
    print(Ad)
    print(Ad2)
    for j in range(init,iterations):
        #Times for equivalent gate duration 0.66 cx, 0.06 H, 4 measuring
        #wait=(0.660*2+0.06*2+4)*j

        wait=(0.533*2+0.06*2+1.216)*j
        t0=wait+ttrav+tH+tCX+tM

        t1cyc=(tM+2*tCX+2*tH)*2+(tM+2*tCX)*2
        t1=t1cyc+t1cyc*3+ttrav*3+t1cyc+t1cyc*j+tM
        

        t2cyc=2*tH+4*tCX+tM
        t2=t2cyc+t2cyc*3+ttrav*3+t2cyc+t2cyc*j+tM
        t3=t2cyc+t2cyc*3+ttrav*3+t2cyc+t2cyc*j+tM

        print(" ")
        print("tcycle, tmerge, t travv,t total")
        print(t1cyc, t1cyc+t1cyc*3,ttrav*3,t1cyc+t1cyc*3+ttrav*3+t1cyc+t1cyc*1+tM)
        print(" ")
        
        
        ID_err2=Noise_model(wait,l_p,t_p)

        print("CX idle----------------------------------------------------------------")
        print(CX_err[1][0],CX_err[1][1],CX_err[1][2])
        print("H idle----------------------------------------------------------------")
        print(H_err[1][0],H_err[1][1],H_err[1][2])
        print("M idle---------------------------------------------------------------")
        print(M_err[1][0],M_err[1][1],M_err[1][2])
        print("CX_loc| H_loc------depolarization")
        print(Loc_CNOT[2],Loc_H[2][0])
        print(" ")
        print("Read_err")
        print(Read_err[2])




        
        circsurf_notr=stim.Circuit('''
    QUBIT_COORDS(0, 0) 0
    QUBIT_COORDS(1, 0) 1
    QUBIT_COORDS(2, 0) 2
    QUBIT_COORDS(3, 0) 3
    QUBIT_COORDS(4, 0) 4
    QUBIT_COORDS(0, 1) 5
    QUBIT_COORDS(1, 1) 6
    QUBIT_COORDS(2, 1) 7
    QUBIT_COORDS(3, 1) 8
    QUBIT_COORDS(4, 1) 9
    QUBIT_COORDS(0, 2) 10
    QUBIT_COORDS(1, 2) 11
    QUBIT_COORDS(2, 2) 12
    QUBIT_COORDS(3, 2) 13
    QUBIT_COORDS(4, 2) 14
    QUBIT_COORDS(0, 3) 15
    QUBIT_COORDS(1, 3) 16
    QUBIT_COORDS(2, 3) 17
    QUBIT_COORDS(3, 3) 18
    QUBIT_COORDS(4, 3) 19
    QUBIT_COORDS(0, 4) 20
    QUBIT_COORDS(1, 4) 21
    QUBIT_COORDS(2, 4) 22
    QUBIT_COORDS(3, 4) 23
    QUBIT_COORDS(4, 4) 24
    QUBIT_COORDS(5, 0) 25
    QUBIT_COORDS(5, 1) 26
    QUBIT_COORDS(5, 2) 27
    QUBIT_COORDS(5, 3) 28
    QUBIT_COORDS(5, 4) 29
    QUBIT_COORDS(6, 0) 30
    QUBIT_COORDS(7, 0) 31
    QUBIT_COORDS(8, 0) 32
    QUBIT_COORDS(9, 0) 33
    QUBIT_COORDS(10, 0) 34
    QUBIT_COORDS(6, 1) 35
    QUBIT_COORDS(7, 1) 36
    QUBIT_COORDS(8, 1) 37
    QUBIT_COORDS(9, 1) 38
    QUBIT_COORDS(10, 1) 39
    QUBIT_COORDS(6, 2) 40
    QUBIT_COORDS(7, 2) 41
    QUBIT_COORDS(8, 2) 42
    QUBIT_COORDS(9, 2) 43
    QUBIT_COORDS(10, 2) 44
    QUBIT_COORDS(6, 3) 45
    QUBIT_COORDS(7, 3) 46
    QUBIT_COORDS(8, 3) 47
    QUBIT_COORDS(9, 3) 48
    QUBIT_COORDS(10, 3) 49
    QUBIT_COORDS(6, 4) 50
    QUBIT_COORDS(7, 4) 51
    QUBIT_COORDS(8, 4) 52
    QUBIT_COORDS(9, 4) 53
    QUBIT_COORDS(10, 4) 54

    #data qubits into |+> state
    RX 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54
    TICK
    
    #x,z stabilizer qubits
    R 1 3 5 7 9 11 13 15 17 19 21 23
    R 31 33 35 37 39 41 43 45 47 49 51 53
    R 25 27 29
    TICK
    #Merging------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    H 1 3 11 13 21 23 25 27 29 31 33 41 43 51 53
    TICK
    CX 1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47 29 50 27 40 25 30 28 19 26 9
    TICK
    CX 1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49 27 28 25 26
    TICK
    CX 11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49 29 28 27 26
    TICK
    CX 1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49 29 24 27 14 25 4 26 35 28 45
    TICK
    H 1 3 11 13 21 23 25 27 29 31 33 41 43 51 53
    TICK
    MR 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53

    DETECTOR(1, 0, 0) rec[-27]
    DETECTOR(1, 2, 0) rec[-22]
    DETECTOR(1, 4, 0) rec[-17]

    DETECTOR(3, 0, 0) rec[-26]
    DETECTOR(3, 2, 0) rec[-21]
    DETECTOR(3, 4, 0) rec[-16]

    DETECTOR(5, 0, 0) rec[-15]
    DETECTOR(5, 2, 0) rec[-14]
    DETECTOR(5, 4, 0) rec[-13]

    DETECTOR(7, 0, 0) rec[-12]
    DETECTOR(7, 2, 0) rec[-7]
    DETECTOR(7, 4, 0) rec[-2]

    DETECTOR(9, 0, 0) rec[-11]
    DETECTOR(9, 2, 0) rec[-6]
    DETECTOR(9, 4, 0) rec[-1]


    REPEAT 3 {
        #Delay for information to send----------------------------------------------------------------------------------------------------------------------------------------------------------------------
        PAULI_CHANNEL_1('''+str(ID_err[1][0])+''', '''+str(ID_err[1][1])+''', '''+str(ID_err[1][2])+''') 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54

        
        H 1 3 11 13 21 23 25 27 29 31 33 41 43 51 53
        PAULI_CHANNEL_1('''+str(H_err[1][0])+''', '''+str(H_err[1][1])+''', '''+str(H_err[1][2])+''') 0 2 4 5 6 7 8 9 10 12 14 15 16 17 18 19 20 22 24 26 28 30 32 34 35 36 37 38 39 40 42 44 45 46 47 48 49 50 52 54
        DEPOLARIZE1('''+str(Loc_H[2][0])+''') 1 3 11 13 21 23 25 27 29 31 33 41 43 51 53
        TICK

        CX 1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47 29 50 27 40 25 30 28 19 26 9

        PAULI_CHANNEL_1('''+str(CX_err[1][0])+''', '''+str(CX_err[1][1])+''', '''+str(CX_err[1][2])+''') 0 10 20 39
        DEPOLARIZE2('''+str(Loc_CNOT[2])+''') 1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47 29 50 27 40 25 30 28 19 26 9
        TICK

        CX 1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49 27 28 25 26
        PAULI_CHANNEL_1('''+str(CX_err[1][0])+''', '''+str(CX_err[1][1])+''', '''+str(CX_err[1][2])+''') 0 2 4 21 23 29 30 32 34 51 53
        DEPOLARIZE2('''+str(Loc_CNOT[2])+''') 1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49 27 28 25 26
        TICK

        CX 11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49 29 28 27 26
        PAULI_CHANNEL_1('''+str(CX_err[1][0])+''', '''+str(CX_err[1][1])+''', '''+str(CX_err[1][2])+''') 1 3 20 22 24 25 31 33 50 52 54
        DEPOLARIZE2('''+str(Loc_CNOT[2])+''') 11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49 29 28 27 26
        TICK

        CX 1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49 29 24 27 14 25 4 26 35 28 45
        PAULI_CHANNEL_1('''+str(CX_err[1][0])+''', '''+str(CX_err[1][1])+''', '''+str(CX_err[1][2])+''') 5 15 34 44 54
        DEPOLARIZE2('''+str(Loc_CNOT[2])+''') 1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49 29 24 27 14 25 4 26 35 28 45
        TICK



        H 1 3 11 13 21 23 25 27 29 31 33 41 43 51 53
        PAULI_CHANNEL_1('''+str(H_err[1][0])+''', '''+str(H_err[1][1])+''', '''+str(H_err[1][2])+''') 1 3 11 13 21 23 25 27 29 31 33 41 43 51 53
        DEPOLARIZE1('''+str(Loc_H[2][0])+''') 0 2 4 5 6 7 8 9 10 12 14 15 16 17 18 19 20 22 24 26 28 30 32 34 35 36 37 38 39 40 42 44 45 46 47 48 49 50 52 54

        TICK

        X_ERROR('''+str(Read_err[2])+''') 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53
        MR 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53
        PAULI_CHANNEL_1('''+str(M_err[1][0])+''', '''+str(M_err[1][1])+''', '''+str(M_err[1][2])+''') 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54




                 


        SHIFT_COORDS(0, 0, 1)
        DETECTOR(1, 0, 0) rec[-27] rec[-54]
        DETECTOR(3, 0, 0) rec[-26] rec[-53]
        DETECTOR(0, 1, 0) rec[-25] rec[-52]
        DETECTOR(2, 1, 0) rec[-24] rec[-51]
        DETECTOR(4, 1, 0) rec[-23] rec[-50]
        DETECTOR(1, 2, 0) rec[-22] rec[-49]
        DETECTOR(3, 2, 0) rec[-21] rec[-48]
        DETECTOR(0, 3, 0) rec[-20] rec[-47]
        DETECTOR(2, 3, 0) rec[-19] rec[-46]
        DETECTOR(4, 3, 0) rec[-18] rec[-45] 
        DETECTOR(1, 4, 0) rec[-17] rec[-44]
        DETECTOR(3, 4, 0) rec[-16] rec[-43]

        DETECTOR(5, 0, 0) rec[-15] rec[-42] 
        DETECTOR(5, 2, 0) rec[-14] rec[-41]
        DETECTOR(5, 4, 0) rec[-13] rec[-40]

        DETECTOR(7, 0, 0) rec[-12] rec[-39]
        DETECTOR(9, 0, 0) rec[-11] rec[-38]
        DETECTOR(6, 1, 0) rec[-10] rec[-37]


        DETECTOR(8, 1, 0) rec[-9] rec[-36]
        DETECTOR(10, 1, 0) rec[-8] rec[-35]
        DETECTOR(7, 2, 0) rec[-7] rec[-34]
        DETECTOR(9, 2, 0) rec[-6] rec[-33]
        DETECTOR(6, 3, 0) rec[-5] rec[-32]
        DETECTOR(8, 3, 0) rec[-4] rec[-31]
        DETECTOR(10, 3, 0) rec[-3] rec[-30] 
        DETECTOR(7, 4, 0) rec[-2] rec[-29]
        DETECTOR(9, 4, 0) rec[-1] rec[-28]


    }
    ##################################################################################
    #Measurement of the data qubits in the intermediate section
    X_ERROR('''+str(Read_err[2])+''') 26 28
    MZ 26 28    
    PAULI_CHANNEL_1('''+str(M_err[1][0])+''', '''+str(M_err[1][1])+''', '''+str(M_err[1][2])+''') 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54


    ##################################################################################

    R 1 3 5 7 9 11 13 15 17 19 21 23
    R 31 33 35 37 39 41 43 45 47 49 51 53

    TICK
    H 1 3 11 13 21 23 31 33 41 43 51 53
    PAULI_CHANNEL_1('''+str(H_err[1][0])+''', '''+str(H_err[1][1])+''', '''+str(H_err[1][2])+''') 0 2 4 5 6 7 8 9 10 12 14 15 16 17 18 19 20 22 24 25 27 29 30 32 34 35 36 37 38 39 40 42 44 45 46 47 48 49 50 52 54
    DEPOLARIZE1('''+str(Loc_H[2][0])+''') 1 3 11 13 21 23 31 33 41 43 51 53

    TICK
    CX 1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47 

    PAULI_CHANNEL_1('''+str(CX_err[1][0])+''', '''+str(CX_err[1][1])+''', '''+str(CX_err[1][2])+''') 0 9 10 19 20 25 27 29 30 39 40 49 50
    DEPOLARIZE2('''+str(Loc_CNOT[2])+''') 1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47  
    TICK
    CX 1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49 

    PAULI_CHANNEL_1('''+str(CX_err[1][0])+''', '''+str(CX_err[1][1])+''', '''+str(CX_err[1][2])+''') 0 2 4 21 23 25 27 29 30 32 34 51 53
    DEPOLARIZE2('''+str(Loc_CNOT[2])+''') 1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49  
    TICK

    CX 11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49
    PAULI_CHANNEL_1('''+str(CX_err[1][0])+''', '''+str(CX_err[1][1])+''', '''+str(CX_err[1][2])+''') 1 3 20 22 24 25 27 29 31 33 50 52 54
    DEPOLARIZE2('''+str(Loc_CNOT[2])+''') 11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49 
    TICK
    CX 1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49
    PAULI_CHANNEL_1('''+str(CX_err[1][0])+''', '''+str(CX_err[1][1])+''', '''+str(CX_err[1][2])+''') 4 5 14 15 24 25 27 29 34 35 44 45 54
    DEPOLARIZE2('''+str(Loc_CNOT[2])+''') 1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49  
    TICK
    H 1 3 11 13 21 23 31 33 41 43 51 53
    #PAULI_CHANNEL_1('''+str(H_err[1][0])+''', '''+str(H_err[1][1])+''', '''+str(H_err[1][2])+''') 0 2 4 5 6 7 8 9 10 12 14 15 16 17 18 19 20 22 24 25 27  29 30 32 34 35 36 37 38 39 40 42 44 45 46 47 48 49 50 52 54
    DEPOLARIZE1('''+str(Loc_H[2][0])+''') 1 3 11 13 21 23 31 33 41 43 51 53

    TICK
    X_ERROR('''+str(Read_err[2])+''') 1 3 5 7 9 11 13 15 17 19 21 23 31 33 35 37 39 41 43 45 47 49 51 53
    MR 1 3 5 7 9 11 13 15 17 19 21 23 31 33 35 37 39 41 43 45 47 49 51 53
    PAULI_CHANNEL_1('''+str(M_err[1][0])+''', '''+str(M_err[1][1])+''', '''+str(M_err[1][2])+''') 0 2 4 6 8 10 12 14 16 18 20 22 24 25 27 29 30 32 34 36 38 40 42 44 46 48 50 52 54

     #25 27 29 


    DETECTOR(1, 0, 0) rec[-24] rec[-53]
    DETECTOR(3, 0, 0) rec[-23] rec[-52]
    DETECTOR(0, 1, 0) rec[-22] rec[-51]
    DETECTOR(2, 1, 0) rec[-21] rec[-50]
    DETECTOR(4, 1, 0) rec[-20] rec[-49] rec[-26]
    DETECTOR(1, 2, 0) rec[-19] rec[-48]
    DETECTOR(3, 2, 0) rec[-18] rec[-47]
    DETECTOR(0, 3, 0) rec[-17] rec[-46]
    DETECTOR(2, 3, 0) rec[-16] rec[-45]
    DETECTOR(4, 3, 0) rec[-15] rec[-44] rec[-25] 
    DETECTOR(1, 4, 0) rec[-14] rec[-43]
    DETECTOR(3, 4, 0) rec[-13] rec[-42]

    DETECTOR(7, 0, 0) rec[-12] rec[-38]
    DETECTOR(9, 0, 0) rec[-11] rec[-37]
    DETECTOR(6, 1, 0) rec[-10] rec[-36] rec[-26]
    DETECTOR(8, 1, 0) rec[-9] rec[-35]
    DETECTOR(10, 1, 0) rec[-8] rec[-34]
    DETECTOR(7, 2, 0) rec[-7] rec[-33]
    DETECTOR(9, 2, 0) rec[-6] rec[-32]
    DETECTOR(6, 3, 0) rec[-5] rec[-31] rec[-25]
    DETECTOR(8, 3, 0) rec[-4] rec[-30]
    DETECTOR(10, 3, 0) rec[-3] rec[-29] 
    DETECTOR(7, 4, 0) rec[-2] rec[-28]
    DETECTOR(9, 4, 0) rec[-1] rec[-27]



    R 1 3 5 7 9 11 13 15 17 19 21 23
    R 31 33 35 37 39 41 43 45 47 49 51 53

    REPEAT '''+str(j)+''' {
        
        TICK
        H 1 3 11 13 21 23 31 33 41 43 51 53
        PAULI_CHANNEL_1('''+str(H_err[1][0])+''', '''+str(H_err[1][1])+''', '''+str(H_err[1][2])+''') 0 2 4 5 6 7 8 9 10 12 14 15 16 17 18 19 20 22 24 25 27 29 30 32 34 35 36 37 38 39 40 42 44 45 46 47 48 49 50 52 54
        DEPOLARIZE1('''+str(Loc_H[2][0])+''') 1 3 11 13 21 23 31 33 41 43 51 53
        
        TICK
        CX 1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47 
        PAULI_CHANNEL_1('''+str(CX_err[1][0])+''', '''+str(CX_err[1][1])+''', '''+str(CX_err[1][2])+''') 1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47 
        DEPOLARIZE2('''+str(Loc_CNOT[2])+''') 1 2 11 12 21 22 3 4 13 14 23 24 6 5 16 15 8 7 18 17 31 32 41 42 51 52 33 34 43 44 53 54 36 35 46 45 38 37 48 47  
        
        TICK
        CX 1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49 
        PAULI_CHANNEL_1('''+str(CX_err[1][0])+''', '''+str(CX_err[1][1])+''', '''+str(CX_err[1][2])+''') 0 2 4 7 9 14 15 25 27 29 30 32 34 36 37 40 42 44 46 48 50 51 53
        DEPOLARIZE2('''+str(Loc_CNOT[2])+''') 1 6 11 16 3 8 13 18 10 5 20 15 12 7 22 17 14 9 24 19 31 36 41 46 33 38 43 48 40 35 50 45 42 37 52 47 44 39 54 49  
        
        TICK
        CX 11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49 
        PAULI_CHANNEL_1('''+str(CX_err[1][0])+''', '''+str(CX_err[1][1])+''', '''+str(CX_err[1][2])+''') 1 3 4 9 20 22 24 25 27 29 31 33 34 35 36 39 40 42 44 45 49 50 52 54
        DEPOLARIZE2('''+str(Loc_CNOT[2])+''') 11 6 21 16 13 8 23 18 0 5 10 15 2 7 12 17 4 9 14 19 41 36 51 46 43 38 53 48 30 35 40 45 32 37 42 47 34 39 44 49 
        
        TICK
        CX 1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49 
        PAULI_CHANNEL_1('''+str(CX_err[1][0])+''', '''+str(CX_err[1][1])+''', '''+str(CX_err[1][2])+''') 4 5 14 15 24 25 27 29 34 35 44 45 54
        DEPOLARIZE2('''+str(Loc_CNOT[2])+''') 1 0 11 10 21 20 3 2 13 12 23 22 6 7 16 17 8 9 18 19 31 30 41 40 51 50 33 32 43 42 53 52 36 37 46 47 38 39 48 49  
        
        TICK
        H 1 3 11 13 21 23 31 33 41 43 51 53
        #PAULI_CHANNEL_1('''+str(H_err[1][0])+''', '''+str(H_err[1][1])+''', '''+str(H_err[1][2])+''') 0 2 4 5 6 7 8 9 10 12 14 15 16 17 18 19 20 22 24 25 27 29 30 32 34 35 36 37 38 39 40 42 44 45 46 47 48 49 50 52 54
        DEPOLARIZE1('''+str(Loc_H[2][0])+''') 1 3 11 13 21 23 31 33 41 43 51 53
        
        TICK
        X_ERROR('''+str(Read_err[2])+''') 1 3  5  7  9 11 13 15 17 19 21 23 31 33 35 37 39 41 43 45 47 49 51 53
        MR 1 3  5  7  9 11 13 15 17 19 21 23 31 33 35 37 39 41 43 45 47 49 51 53
        PAULI_CHANNEL_1('''+str(M_err[1][0])+''', '''+str(M_err[1][1])+''', '''+str(M_err[1][2])+''') 0 2 4 6 8 10 12 14 16 18 20 22 24 25 29 30 32 34 36 38 40 42 44 46 48 50 52 54
         

        SHIFT_COORDS(0, 0, 1)
        DETECTOR(1, 0, 0) rec[-24] rec[-48]
        DETECTOR(3, 0, 0) rec[-23] rec[-47]
        DETECTOR(0, 1, 0) rec[-22] rec[-46]
        DETECTOR(2, 1, 0) rec[-21] rec[-45]
        DETECTOR(4, 1, 0) rec[-20] rec[-44]
        DETECTOR(1, 2, 0) rec[-19] rec[-43]
        DETECTOR(3, 2, 0) rec[-18] rec[-42]
        DETECTOR(0, 3, 0) rec[-17] rec[-41]
        DETECTOR(2, 3, 0) rec[-16] rec[-40]
        DETECTOR(4, 3, 0) rec[-15] rec[-39] 
        DETECTOR(1, 4, 0) rec[-14] rec[-38]
        DETECTOR(3, 4, 0) rec[-13] rec[-37]

      

        DETECTOR(7, 0, 0) rec[-12] rec[-36]
        DETECTOR(9, 0, 0) rec[-11] rec[-35]
        DETECTOR(6, 1, 0) rec[-10] rec[-34]
        DETECTOR(8, 1, 0) rec[-9] rec[-33]
        DETECTOR(10, 1, 0) rec[-8] rec[-32]
        DETECTOR(7, 2, 0) rec[-7] rec[-31]
        DETECTOR(9, 2, 0) rec[-6] rec[-30]
        DETECTOR(6, 3, 0) rec[-5] rec[-29]
        DETECTOR(8, 3, 0) rec[-4] rec[-28]
        DETECTOR(10, 3, 0) rec[-3] rec[-27] 
        DETECTOR(7, 4, 0) rec[-2] rec[-26]
        DETECTOR(9, 4, 0) rec[-1] rec[-25]


    }



    X_ERROR('''+str(Read_err[2])+''') 0 2 4 6 8  10 12 14 16 18 20 22 24 30 32 34 36 38 40  42 44 46 48 50 52 54
    MX 0 2 4 6 8  10 12 14 16 18 20 22 24 30 32 34 36 38 40 42 44 46 48 50 52 54


    
        
    DETECTOR(1, 0, 1) rec[-23] rec[-25] rec[-26] rec[-50]                   
    DETECTOR(1, 2, 1) rec[-18] rec[-20] rec[-21] rec[-23] rec[-45]
    DETECTOR(1, 4, 1) rec[-15] rec[-16] rec[-18] rec[-40]
    
    DETECTOR(3, 0, 1) rec[-22] rec[-24] rec[-25] rec[-49]
    DETECTOR(3, 2, 1) rec[-17] rec[-19] rec[-20] rec[-22] rec[-44]
    DETECTOR(3, 4, 1) rec[-14] rec[-15] rec[-17] rec[-39]

    DETECTOR(7, 0, 1) rec[-10] rec[-12] rec[-13] rec[-38]
    DETECTOR(7, 2, 1) rec[-5] rec[-7] rec[-8] rec[-10] rec[-33]
    DETECTOR(7, 4, 1) rec[-2] rec[-3] rec[-5] rec[-28]

    DETECTOR(9, 0, 1) rec[-9] rec[-11] rec[-12] rec[-37]
    DETECTOR(9, 2, 1) rec[-4] rec[-6] rec[-7] rec[-9] rec[-32]
    DETECTOR(9, 4, 1) rec[-1] rec[-2] rec[-4] rec[-27]


    OBSERVABLE_INCLUDE(0) rec[-16] rec[-21] rec[-26] #rec[-3] rec[-8] rec[-13] 
    OBSERVABLE_INCLUDE(1) rec[-3] rec[-8] rec[-13] 


''')    
####################################################################################################################################################################################################################        
        circsurf_notr2=stim.Circuit("""
QUBIT_COORDS(1, 1) 1
QUBIT_COORDS(2, 0) 2
QUBIT_COORDS(3, 1) 3
QUBIT_COORDS(5, 1) 5
QUBIT_COORDS(1, 3) 8
QUBIT_COORDS(2, 2) 9
QUBIT_COORDS(3, 3) 10
QUBIT_COORDS(4, 2) 11
QUBIT_COORDS(5, 3) 12
QUBIT_COORDS(6, 2) 13
QUBIT_COORDS(0, 4) 14
QUBIT_COORDS(1, 5) 15
QUBIT_COORDS(2, 4) 16
QUBIT_COORDS(3, 5) 17
QUBIT_COORDS(4, 4) 18
QUBIT_COORDS(5, 5) 19
QUBIT_COORDS(4, 6) 25

#data
QUBIT_COORDS(7, 1) 0
QUBIT_COORDS(7, 3) 4
QUBIT_COORDS(7, 5) 7

#stabilizers
QUBIT_COORDS(6, 0) 26
QUBIT_COORDS(8, 2) 30
QUBIT_COORDS(6, 4) 32
QUBIT_COORDS(8, 6) 33


QUBIT_COORDS(9, 1) 27
QUBIT_COORDS(10, 0) 28
QUBIT_COORDS(11, 1) 29
QUBIT_COORDS(13, 1) 31
QUBIT_COORDS(9, 3) 34
QUBIT_COORDS(10, 2) 35
QUBIT_COORDS(11, 3) 36
QUBIT_COORDS(12, 2) 37
QUBIT_COORDS(13, 3) 38
QUBIT_COORDS(14, 2) 39
QUBIT_COORDS(8, 4) 40
QUBIT_COORDS(9, 5) 41
QUBIT_COORDS(10, 4) 42
QUBIT_COORDS(11, 5) 43
QUBIT_COORDS(12, 4) 44
QUBIT_COORDS(13, 5) 45
QUBIT_COORDS(12, 6) 51



RX 1 3 5 8 10 12 15 17 19 27 29 31 34 36 38 41 43 45 0 4 7
R 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51 26 30 32 33 
TICK
H 2 11 16 25 28 37 42 51 26 30 32 33
TICK
CX 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44 26 0 32 7 30 34 4 13


TICK
CX 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44 26 5 32 19 30 4 0 13


TICK
CX 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39 30 27 32 4 33 41 7 40


TICK
CX 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39 32 12 30 0 33 7 4 40


TICK
H 2 11 16 25 28 37 42 51 26 30 32 33

TICK
MR 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51 26 30 32 33

DETECTOR(2, 0, 0) rec[-20]
DETECTOR(2, 4, 0) rec[-15]
DETECTOR(4, 2, 0) rec[-18]
DETECTOR(4, 6, 0) rec[-13]

DETECTOR(6, 0, 0) rec[-4]
DETECTOR(8, 2, 0) rec[-3]
DETECTOR(6, 4, 0) rec[-2]
DETECTOR(8, 6, 0) rec[-1]

DETECTOR(10, 0, 0) rec[-12]  
DETECTOR(10, 4, 0) rec[-7]  
DETECTOR(12, 2, 0) rec[-10] 
DETECTOR(12, 6, 0) rec[-5] 




REPEAT 3 {

    #Delay for information to send----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    PAULI_CHANNEL_1("""+str(ID_err[1][0])+""", """+str(ID_err[1][1])+""", """+str(ID_err[1][2])+""") 0 1 2 3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 51


    TICK
    H 2 11 16 25 28 37 42 51 26 30 32 33
    DEPOLARIZE2("""+str(Loc_H[2][0])+""") 2 11 16 25 28 37 42 51 26 30 32 33
    PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 0 1 3 4 5 6 7 8 9 10 12 13 14 15 17 18 19 27 29 31 34 35 36 38 39 40 41 43 44 45
    TICK
    CX 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44 26 0 32 7 30 34 4 13
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44 26 0 32 7 30 34 4 13
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 1 5 6 8 20 21 22 23 24 25 27 31 33 39 46 47 48 49 50 51


    TICK
    CX 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44 26 5 32 19 30 4 0 13
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44 26 5 32 19 30 4 0 13

    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 5 13 17 19 25 31 39 43 45 51

    TICK
    CX 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39 30 27 32 4 33 41 7 40
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39 30 27 32 4 33 41 7 40

    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 1 2 3 6 14 15 20 21 22 23 24 26 28 29 46 47 48 49 50

    TICK
    CX 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39 32 12 30 0 33 7 4 40
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39 32 12 30 0 33 7 4 40

    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 2 6 14 15 19 20 21 22 23 24 26 28 38 41 45 46 47 48 49 50

    TICK
    H 2 11 16 25 28 37 42 51 26 30 32 33
    DEPOLARIZE1("""+str(Loc_H[2][0])+""") 2 11 16 25 28 37 42 51
    PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 0 1 3 4 5 6 7 8 9 10 12 13 14 15 17 18 19 20 21 22 23 24 26 27 29 30 31 32 33 34 35 36 38 39 40 41 43 44 45 46 47 48 49 50

    TICK
    X_ERROR("""+str(Read_err[2])+""") 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51 26 30 32 33
    MR 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51 26 30 32 33
    PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 0 1 3 4 5 6 7 8 10 12 15 17 19 20 21 22 23 24 27 29 31 34 36 38 41 43 45 46 47 48 49 50

    SHIFT_COORDS(0, 0, 1)

    DETECTOR(2, 0, 0) rec[-20] rec[-40]
    DETECTOR(2, 2, 0) rec[-19] rec[-39]
    DETECTOR(4, 2, 0) rec[-18] rec[-38]
    DETECTOR(6, 2, 0) rec[-17] rec[-37]
    DETECTOR(0, 4, 0) rec[-16] rec[-36]
    DETECTOR(2, 4, 0) rec[-15] rec[-35]
    DETECTOR(4, 4, 0) rec[-14] rec[-34]
    DETECTOR(4, 6, 0) rec[-13] rec[-33]

    DETECTOR(10, 0, 0) rec[-12] rec[-32]   
    DETECTOR(10, 2, 0) rec[-11] rec[-31]   
    DETECTOR(12, 2, 0) rec[-10] rec[-30]  
    DETECTOR(14, 2, 0) rec[-9] rec[-29]  
    DETECTOR(8, 4, 0) rec[-8] rec[-28]   
    DETECTOR(10, 4, 0) rec[-7] rec[-27]   
    DETECTOR(12, 4, 0) rec[-6] rec[-26]  
    DETECTOR(12, 6, 0) rec[-5] rec[-25]

    DETECTOR(6, 0, 0) rec[-4] rec[-24]   
    DETECTOR(8, 2, 0) rec[-3] rec[-23]   
    DETECTOR(6, 4, 0) rec[-2] rec[-22]  
    DETECTOR(8, 6, 0) rec[-1] rec[-21]

}

#################### splitting ################################################################
MZ 0 4 7
PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54

#Reset auxiliary qubits
R 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51 
#

TICK
H 2 11 16 25 28 37 42 51
DEPOLARIZE2("""+str(Loc_H[2][0])+""") 2 11 16 25 28 37 42 51
PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 1 3 5 8 9 10 12 13 14 15 17 18 19 27 29 31 34 35 36 38 39 40 41 43 44 45 
TICK
CX 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44
DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44
PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 5 13 1 19 25 31 39 27 45 51

TICK
CX 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44
DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44

PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 5 13 17 19 25 31 39 43 45 51

TICK
CX 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39
DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39

PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 1 2 3 14 15 27 28 29 40 41

TICK
CX 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39
DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39

PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 1 3 5 8 9 10 12 13 14 15 17 18 19 27 29 31 34 35 36 38 39 40 41 43 44 45 

TICK
H 2 11 16 25 28 37 42 51
DEPOLARIZE1("""+str(Loc_H[2][0])+""") 2 11 16 25 28 37 42 51
PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 1 3 5 8 9 10 12 13 14 15 17 18 19 27 29 31 34 35 36 38 39 40 41 43 44 45 
TICK
X_ERROR("""+str(Read_err[2])+""") 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51
MR 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51
PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 1 3 5 8 10 12 15 17 19 27 29 31 34 36 38 41 43 45
SHIFT_COORDS(0, 0, 1)

DETECTOR(2, 0, 0) rec[-16] rec[-39]
DETECTOR(2, 2, 0) rec[-15] rec[-38]
DETECTOR(4, 2, 0) rec[-14] rec[-37]
DETECTOR(6, 2, 0) rec[-13] rec[-36] rec[-18] rec[-19] 
DETECTOR(0, 4, 0) rec[-12] rec[-35]
DETECTOR(2, 4, 0) rec[-11] rec[-34]
DETECTOR(4, 4, 0) rec[-10] rec[-33]
DETECTOR(4, 6, 0) rec[-9] rec[-32]

DETECTOR(10, 0, 0) rec[-8] rec[-31]   
DETECTOR(10, 2, 0) rec[-7] rec[-30]   
DETECTOR(12, 2, 0) rec[-6] rec[-29]  
DETECTOR(14, 2, 0) rec[-5] rec[-28]  
DETECTOR(8, 4, 0) rec[-4] rec[-27] rec[-17] rec[-18] 
DETECTOR(10, 4, 0) rec[-3] rec[-26]   
DETECTOR(12, 4, 0) rec[-2] rec[-25]  
DETECTOR(12, 6, 0) rec[-1] rec[-24]

R 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51 

REPEAT """+str(j)+""" {
    TICK
    H 2 11 16 25 28 37 42 51
    DEPOLARIZE2("""+str(Loc_H[2][0])+""") 2 11 16 25 28 37 42 51
    PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 1 3 5 8 9 10 12 13 14 15 17 18 19 27 29 31 34 35 36 38 39 40 41 43 44 45 
    TICK
    CX 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 2 3 16 17 11 12 15 14 10 9 19 18 28 29 42 43 37 38 41 40 36 35 45 44
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 5 13 1 19 25 31 39 27 45 51

    TICK
    CX 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 2 1 16 15 11 10 8 14 3 9 12 18 28 27 42 41 37 36 34 40 29 35 38 44

    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 5 13 17 19 25 31 39 43 45 51

    TICK
    CX 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 16 10 11 5 25 19 8 9 17 18 12 13 42 36 37 31 51 45 34 35 43 44 38 39

    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 1 2 3 14 15 27 28 29 40 41

    TICK
    CX 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 16 8 11 3 25 17 1 9 10 18 5 13 42 34 37 29 51 43 27 35 36 44 31 39

    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 1 3 5 8 9 10 12 13 14 15 17 18 19 27 29 31 34 35 36 38 39 40 41 43 44 45 

    TICK
    H 2 11 16 25 28 37 42 51
    DEPOLARIZE1("""+str(Loc_H[2][0])+""") 2 11 16 25 28 37 42 51
    PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 1 3 5 8 9 10 12 13 14 15 17 18 19 27 29 31 34 35 36 38 39 40 41 43 44 45 
    TICK
    X_ERROR("""+str(Read_err[2])+""") 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51
    MR 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51
    PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 1 3 5 8 10 12 15 17 19 27 29 31 34 36 38 41 43 45
    SHIFT_COORDS(0, 0, 1)
    DETECTOR(2, 0, 0) rec[-16] rec[-32]
    DETECTOR(2, 2, 0) rec[-15] rec[-31]
    DETECTOR(4, 2, 0) rec[-14] rec[-30]
    DETECTOR(6, 2, 0) rec[-13] rec[-29]
    DETECTOR(0, 4, 0) rec[-12] rec[-28]
    DETECTOR(2, 4, 0) rec[-11] rec[-27]
    DETECTOR(4, 4, 0) rec[-10] rec[-26]
    DETECTOR(4, 6, 0) rec[-9] rec[-25]

    DETECTOR(10, 0, 0) rec[-8] rec[-24]   
    DETECTOR(10, 2, 0) rec[-7] rec[-23]   
    DETECTOR(12, 2, 0) rec[-6] rec[-22]  
    DETECTOR(14, 2, 0) rec[-5] rec[-21]  
    DETECTOR(8, 4, 0) rec[-4] rec[-20]   
    DETECTOR(10, 4, 0) rec[-3] rec[-19]   
    DETECTOR(12, 4, 0) rec[-2] rec[-18]  
    DETECTOR(12, 6, 0) rec[-1] rec[-17]
}

X_ERROR("""+str(Read_err[2])+""") 1 3 5 8 10 12 15 17 19 27 29 31 34 36 38 41 43 45
MX 1 3 5 8 10 12 15 17 19 27 29 31 34 36 38 41 43 45
PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 2 9 11 13 14 16 18 25 28 35 37 39 40 42 44 51

DETECTOR(2, 0, 1) rec[-17] rec[-18] rec[-34]#
DETECTOR(2, 4, 1) rec[-11] rec[-12] rec[-14] rec[-15] rec[-29]#

DETECTOR(4, 2, 1) rec[-13] rec[-14] rec[-16] rec[-17] rec[-32]
DETECTOR(4, 6, 1) rec[-10] rec[-11] rec[-27] #


DETECTOR(10, 0, 1) rec[-8] rec[-9] rec[-26] #
DETECTOR(10, 4, 1) rec[-2] rec[-3] rec[-5] rec[-6] rec[-21] #

DETECTOR(12, 2, 1) rec[-4] rec[-5] rec[-7] rec[-8] rec[-24]#
DETECTOR(12, 6, 1) rec[-1] rec[-2] rec[-19] #

OBSERVABLE_INCLUDE(0) rec[-12] rec[-15] rec[-18] #rec[-3] rec[-6] rec[-9] 

OBSERVABLE_INCLUDE(1) rec[-3] rec[-6] rec[-9] 
               
""")
#######################################################################################################################################################################################################################################
        circuit2=stim.Circuit("""
            QUBIT_COORDS(0,1) 0
            QUBIT_COORDS(0,2) 1

            H 0
            DEPOLARIZE1("""+str(Loc_H[2][0])+""") 0
            PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 1

            CX 0 1
            DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 1


            #Travel distance

            #TICK
            #PAULI_CHANNEL_1("""+str(ID_err[1][0])+""", """+str(ID_err[1][1])+""", """+str(ID_err[1][2])+""") 0 1 

            #Waiting time in a quantum memory + travelling

            TICK 
            PAULI_CHANNEL_1("""+str(ID_err2[1][0])+""", """+str(ID_err2[1][1])+""", """+str(ID_err2[1][2])+""") 0 1 
            

            TICK
            X_ERROR("""+str(Read_err[2])+""") 0 1
            M 0 1
            DETECTOR rec[-1] rec[-2]

            """)

        circuit1=stim.Circuit("""
QUBIT_COORDS(0, 0) 0
QUBIT_COORDS(0, 1) 1
QUBIT_COORDS(0, 2) 2
QUBIT_COORDS(0, 3) 3
QUBIT_COORDS(0, 4) 4
QUBIT_COORDS(1, 0) 5
QUBIT_COORDS(1, 1) 6
QUBIT_COORDS(1, 2) 7
QUBIT_COORDS(1, 3) 8
QUBIT_COORDS(1, 4) 9
QUBIT_COORDS(2, 0) 10
QUBIT_COORDS(2, 1) 11
QUBIT_COORDS(2, 2) 12
QUBIT_COORDS(2, 3) 13
QUBIT_COORDS(2, 4) 14
QUBIT_COORDS(3, 0) 15
QUBIT_COORDS(3, 1) 16
QUBIT_COORDS(3, 2) 17
QUBIT_COORDS(3, 3) 18
QUBIT_COORDS(3, 4) 19
QUBIT_COORDS(4, 0) 20
QUBIT_COORDS(4, 1) 21
QUBIT_COORDS(4, 2) 22
QUBIT_COORDS(4, 3) 23
QUBIT_COORDS(4, 4) 24
QUBIT_COORDS(5, 0) 25
QUBIT_COORDS(5, 1) 26
QUBIT_COORDS(5, 2) 27
QUBIT_COORDS(5, 3) 28
QUBIT_COORDS(5, 4) 29
QUBIT_COORDS(6, 0) 30
QUBIT_COORDS(6, 1) 31
QUBIT_COORDS(6, 2) 32
QUBIT_COORDS(6, 3) 33
QUBIT_COORDS(6, 4) 34
QUBIT_COORDS(7, 0) 35
QUBIT_COORDS(7, 1) 36
QUBIT_COORDS(7, 2) 37
QUBIT_COORDS(7, 3) 38
QUBIT_COORDS(7, 4) 39
QUBIT_COORDS(8, 0) 40
QUBIT_COORDS(8, 1) 41
QUBIT_COORDS(8, 2) 42
QUBIT_COORDS(8, 3) 43
QUBIT_COORDS(8, 4) 44
QUBIT_COORDS(9, 0) 45
QUBIT_COORDS(9, 1) 46
QUBIT_COORDS(9, 2) 47
QUBIT_COORDS(9, 3) 48
QUBIT_COORDS(9, 4) 49
#ignore creation in first phase of merging
RX 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
TICK
MPP X0*X5 X1*X6 X2*X7 X3*X8 X4*X9 X10*X15 X11*X16 X12*X17 X13*X18 X14*X19 X30*X35 X31*X36 X32*X37 X33*X38 X34*X39 X40*X45 X41*X46 X42*X47 X43*X48 X44*X49
TICK
MPP X5*X10 X6*X11 X7*X12 X8*X13 X9*X14 X15*X20 X16*X21 X17*X22 X18*X23 X19*X24 X25*X30 X26*X31 X27*X32 X28*X33 X29*X34 X35*X40 X36*X41 X37*X42 X38*X43 X39*X44
TICK
MPP Z0*Z1 Z2*Z3 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z15*Z16 Z17*Z18 Z20*Z21 Z22*Z23 Z25*Z26 Z27*Z28 Z30*Z31 Z32*Z33 Z35*Z36 Z37*Z38 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48
TICK
MPP Z1*Z2 Z3*Z4 Z6*Z7 Z8*Z9 Z11*Z12 Z13*Z14 Z16*Z17 Z18*Z19 Z21*Z22 Z23*Z24 Z26*Z27 Z28*Z29 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49
DETECTOR(0.5, 0, 0) rec[-80]
DETECTOR(0.5, 1, 0) rec[-79]
DETECTOR(0.5, 2, 0) rec[-78]
DETECTOR(0.5, 3, 0) rec[-77]
DETECTOR(0.5, 4, 0) rec[-76]
DETECTOR(1.5, 0, 0) rec[-60]
DETECTOR(1.5, 1, 0) rec[-59]
DETECTOR(1.5, 2, 0) rec[-58]
DETECTOR(1.5, 3, 0) rec[-57]
DETECTOR(1.5, 4, 0) rec[-56]
DETECTOR(2.5, 0, 0) rec[-75]
DETECTOR(2.5, 1, 0) rec[-74]
DETECTOR(2.5, 2, 0) rec[-73]
DETECTOR(2.5, 3, 0) rec[-72]
DETECTOR(2.5, 4, 0) rec[-71]
DETECTOR(3.5, 0, 0) rec[-55]
DETECTOR(3.5, 1, 0) rec[-54]
DETECTOR(3.5, 2, 0) rec[-53]
DETECTOR(3.5, 3, 0) rec[-52]
DETECTOR(3.5, 4, 0) rec[-51]
DETECTOR(5.5, 0, 0) rec[-50]
DETECTOR(5.5, 1, 0) rec[-49]
DETECTOR(5.5, 2, 0) rec[-48]
DETECTOR(5.5, 3, 0) rec[-47]
DETECTOR(5.5, 4, 0) rec[-46]
DETECTOR(6.5, 0, 0) rec[-70]
DETECTOR(6.5, 1, 0) rec[-69]
DETECTOR(6.5, 2, 0) rec[-68]
DETECTOR(6.5, 3, 0) rec[-67]
DETECTOR(6.5, 4, 0) rec[-66]
DETECTOR(7.5, 0, 0) rec[-45]
DETECTOR(7.5, 1, 0) rec[-44]
DETECTOR(7.5, 2, 0) rec[-43]
DETECTOR(7.5, 3, 0) rec[-42]
DETECTOR(7.5, 4, 0) rec[-41]
DETECTOR(8.5, 0, 0) rec[-65]
DETECTOR(8.5, 1, 0) rec[-64]
DETECTOR(8.5, 2, 0) rec[-63]
DETECTOR(8.5, 3, 0) rec[-62]
DETECTOR(8.5, 4, 0) rec[-61]
SHIFT_COORDS(0, 0, 1)

TICK #complete
#DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 5 1 6 2 7 3 8 4 9 10 15 11 16 12 17 13 18 14 19 20 25 21 26 22 27 23 28 24 29 30 35 31 36 32 37 33 38 34 39 40 45 41 46 42 47 43 48 44 49 #q0+q1
#DEPOLARIZE1("""+str(Loc_H[2][0])+""") 0 1 2 3 4 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44 #q0
#PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 45 46 47 48 49
MPP X0*X5 X1*X6 X2*X7 X3*X8 X4*X9 X10*X15 X11*X16 X12*X17 X13*X18 X14*X19 X20*X25 X21*X26 X22*X27 X23*X28 X24*X29 X30*X35 X31*X36 X32*X37 X33*X38 X34*X39 X40*X45 X41*X46 X42*X47 X43*X48 X44*X49
#PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 45 46 47 48 49 #q1
#PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 45 46 47 48 49
#DEPOLARIZE1("""+str(Loc_H[2][0])+""") 0 1 2 3 4 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44
#DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 5 1 6 2 7 3 8 4 9 10 15 11 16 12 17 13 18 14 19 20 25 21 26 22 27 23 28 24 29 30 35 31 36 32 37 33 38 34 39 40 45 41 46 42 47 43 48 44 49

TICK #complete
#DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 5 10 6 11 7 12 8 13 9 14 15 20 16 21 17 22 18 23 19 24 25 30 26 31 27 32 28 33 29 34 35 40 36 41 37 42 38 43 39 44
#DEPOLARIZE1("""+str(Loc_H[2][0])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39
#PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49
#PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44
MPP X5*X10 X6*X11 X7*X12 X8*X13 X9*X14 X15*X20 X16*X21 X17*X22 X18*X23 X19*X24 X25*X30 X26*X31 X27*X32 X28*X33 X29*X34 X35*X40 X36*X41 X37*X42 X38*X43 X39*X44
#DEPOLARIZE1("""+str(Loc_H[2][0])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 #q0
#DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 5 10 6 11 7 12 8 13 9 14 15 20 16 21 17 22 18 23 19 24 25 30 26 31 27 32 28 33 29 34 35 40 36 41 37 42 38 43 39 44 #q0+q1
#PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44 #q1+missing
#PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44 #
#PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 #missing #


TICK #complete
#DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48
#PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 4 9 14 19 24 29 34 39 44 49 #missing 
MPP Z0*Z1 Z2*Z3 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z15*Z16 Z17*Z18 Z20*Z21 Z22*Z23 Z25*Z26 Z27*Z28 Z30*Z31 Z32*Z33 Z35*Z36 Z37*Z38 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48
#PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 0 2 5 7 10 12 15 17 20 22 25 27 30 32 35 37 40 42 45 47 4 9 14 19 24 29 34 39 44 49   #q0+missing
#PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 4 9 14 19 24 29 34 39 44 49 #missing 
#DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48

TICK #complete


#DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49
#PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 5 10 15 20 25 30 35 40 45 #missing 
MPP Z1*Z2 Z3*Z4 Z6*Z7 Z8*Z9 Z11*Z12 Z13*Z14 Z16*Z17 Z18*Z19 Z21*Z22 Z23*Z24 Z26*Z27 Z28*Z29 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49
#PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 1 3 6 8 11 13 16 18 21 23 26 28 31 33 36 38 41 43 46 48 0 5 10 15 20 25 30 35 40 45   #q0+missing
#PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 5 10 15 20 25 30 35 40 45 #missing 
#DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49


DETECTOR(0.5, 0, 0) rec[-165] rec[-164] rec[-163] rec[-162] rec[-161] rec[-85] rec[-84] rec[-83] rec[-82] rec[-81]
DETECTOR(1.5, 0, 0) rec[-145] rec[-144] rec[-143] rec[-142] rec[-141] rec[-60] rec[-59] rec[-58] rec[-57] rec[-56]
DETECTOR(2.5, 0, 0) rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-80] rec[-79] rec[-78] rec[-77] rec[-76]
DETECTOR(3.5, 0, 0) rec[-140] rec[-139] rec[-138] rec[-137] rec[-136] rec[-55] rec[-54] rec[-53] rec[-52] rec[-51]
DETECTOR(5.5, 0, 0) rec[-135] rec[-134] rec[-133] rec[-132] rec[-131] rec[-50] rec[-49] rec[-48] rec[-47] rec[-46]
DETECTOR(6.5, 0, 0) rec[-155] rec[-154] rec[-153] rec[-152] rec[-151] rec[-70] rec[-69] rec[-68] rec[-67] rec[-66]
DETECTOR(7.5, 0, 0) rec[-130] rec[-129] rec[-128] rec[-127] rec[-126] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41]
DETECTOR(8.5, 0, 0) rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-65] rec[-64] rec[-63] rec[-62] rec[-61]
DETECTOR(0, 0.5, 0) rec[-125] rec[-123] rec[-121] rec[-119] rec[-117] rec[-115] rec[-113] rec[-111] rec[-109] rec[-107] rec[-40] rec[-38] rec[-36] rec[-34] rec[-32] rec[-30] rec[-28] rec[-26] rec[-24] rec[-22]
DETECTOR(0, 1.5, 0) rec[-105] rec[-103] rec[-101] rec[-99] rec[-97] rec[-95] rec[-93] rec[-91] rec[-89] rec[-87] rec[-20] rec[-18] rec[-16] rec[-14] rec[-12] rec[-10] rec[-8] rec[-6] rec[-4] rec[-2]
DETECTOR(0, 2.5, 0) rec[-124] rec[-122] rec[-120] rec[-118] rec[-116] rec[-114] rec[-112] rec[-110] rec[-108] rec[-106] rec[-39] rec[-37] rec[-35] rec[-33] rec[-31] rec[-29] rec[-27] rec[-25] rec[-23] rec[-21]
DETECTOR(0, 3.5, 0) rec[-104] rec[-102] rec[-100] rec[-98] rec[-96] rec[-94] rec[-92] rec[-90] rec[-88] rec[-86] rec[-19] rec[-17] rec[-15] rec[-13] rec[-11] rec[-9] rec[-7] rec[-5] rec[-3] rec[-1]
OBSERVABLE_INCLUDE(2) rec[-75] rec[-74] rec[-73] rec[-72] rec[-71]
SHIFT_COORDS(0, 0, 1)

TICK
REPEAT 3 {

    PAULI_CHANNEL_1("""+str(ID_err[1][0])+""", """+str(ID_err[1][1])+""", """+str(ID_err[1][2])+""") 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51
   

    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 5 1 6 2 7 3 8 4 9 10 15 11 16 12 17 13 18 14 19 20 25 21 26 22 27 23 28 24 29 30 35 31 36 32 37 33 38 34 39 40 45 41 46 42 47 43 48 44 49 #q0+q1
    DEPOLARIZE1("""+str(Loc_H[2][0])+""") 0 1 2 3 4 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44 #q0
    PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 45 46 47 48 49
    MPP("""+str(Read_err[2])+""") X0*X5 X1*X6 X2*X7 X3*X8 X4*X9 X10*X15 X11*X16 X12*X17 X13*X18 X14*X19 X20*X25 X21*X26 X22*X27 X23*X28 X24*X29 X30*X35 X31*X36 X32*X37 X33*X38 X34*X39 X40*X45 X41*X46 X42*X47 X43*X48 X44*X49
    PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 45 46 47 48 49 #q1
    PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 45 46 47 48 49
    DEPOLARIZE1("""+str(Loc_H[2][0])+""") 0 1 2 3 4 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 5 1 6 2 7 3 8 4 9 10 15 11 16 12 17 13 18 14 19 20 25 21 26 22 27 23 28 24 29 30 35 31 36 32 37 33 38 34 39 40 45 41 46 42 47 43 48 44 49

    TICK
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 5 10 6 11 7 12 8 13 9 14 15 20 16 21 17 22 18 23 19 24 25 30 26 31 27 32 28 33 29 34 35 40 36 41 37 42 38 43 39 44
    DEPOLARIZE1("""+str(Loc_H[2][0])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49
    PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44
    MPP("""+str(Read_err[2])+""") X5*X10 X6*X11 X7*X12 X8*X13 X9*X14 X15*X20 X16*X21 X17*X22 X18*X23 X19*X24 X25*X30 X26*X31 X27*X32 X28*X33 X29*X34 X35*X40 X36*X41 X37*X42 X38*X43 X39*X44
    DEPOLARIZE1("""+str(Loc_H[2][0])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 #q0
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 5 10 6 11 7 12 8 13 9 14 15 20 16 21 17 22 18 23 19 24 25 30 26 31 27 32 28 33 29 34 35 40 36 41 37 42 38 43 39 44 #q0+q1 
    PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44 #q1+missing
    PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44 #
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 #missing #

    TICK
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 4 9 14 19 24 29 34 39 44 49 #missing 
    MPP("""+str(Read_err[2])+""") Z0*Z1 Z2*Z3 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z15*Z16 Z17*Z18 Z20*Z21 Z22*Z23 Z25*Z26 Z27*Z28 Z30*Z31 Z32*Z33 Z35*Z36 Z37*Z38 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48
    PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 0 2 5 7 10 12 15 17 20 22 25 27 30 32 35 37 40 42 45 47 4 9 14 19 24 29 34 39 44 49   #q0+missing
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 4 9 14 19 24 29 34 39 44 49 #missing 
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48

    TICK
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 5 10 15 20 25 30 35 40 45 #missing 
    MPP("""+str(Read_err[2])+""") Z1*Z2 Z3*Z4 Z6*Z7 Z8*Z9 Z11*Z12 Z13*Z14 Z16*Z17 Z18*Z19 Z21*Z22 Z23*Z24 Z26*Z27 Z28*Z29 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49
    PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 1 3 6 8 11 13 16 18 21 23 26 28 31 33 36 38 41 43 46 48 0 5 10 15 20 25 30 35 40 45   #q0+missing
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 5 10 15 20 25 30 35 40 45 #missing 
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49

    DETECTOR(0.5, 0, 0) rec[-170] rec[-169] rec[-168] rec[-167] rec[-166] rec[-85] rec[-84] rec[-83] rec[-82] rec[-81]
    DETECTOR(1.5, 0, 0) rec[-145] rec[-144] rec[-143] rec[-142] rec[-141] rec[-60] rec[-59] rec[-58] rec[-57] rec[-56]
    DETECTOR(2.5, 0, 0) rec[-165] rec[-164] rec[-163] rec[-162] rec[-161] rec[-80] rec[-79] rec[-78] rec[-77] rec[-76]
    DETECTOR(3.5, 0, 0) rec[-140] rec[-139] rec[-138] rec[-137] rec[-136] rec[-55] rec[-54] rec[-53] rec[-52] rec[-51]
    DETECTOR(4.5, 0, 0) rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-75] rec[-74] rec[-73] rec[-72] rec[-71]
    DETECTOR(5.5, 0, 0) rec[-135] rec[-134] rec[-133] rec[-132] rec[-131] rec[-50] rec[-49] rec[-48] rec[-47] rec[-46]
    DETECTOR(6.5, 0, 0) rec[-155] rec[-154] rec[-153] rec[-152] rec[-151] rec[-70] rec[-69] rec[-68] rec[-67] rec[-66]
    DETECTOR(7.5, 0, 0) rec[-130] rec[-129] rec[-128] rec[-127] rec[-126] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41]
    DETECTOR(8.5, 0, 0) rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-65] rec[-64] rec[-63] rec[-62] rec[-61]
    DETECTOR(0, 0.5, 0) rec[-125] rec[-123] rec[-121] rec[-119] rec[-117] rec[-115] rec[-113] rec[-111] rec[-109] rec[-107] rec[-40] rec[-38] rec[-36] rec[-34] rec[-32] rec[-30] rec[-28] rec[-26] rec[-24] rec[-22]
    DETECTOR(0, 1.5, 0) rec[-105] rec[-103] rec[-101] rec[-99] rec[-97] rec[-95] rec[-93] rec[-91] rec[-89] rec[-87] rec[-20] rec[-18] rec[-16] rec[-14] rec[-12] rec[-10] rec[-8] rec[-6] rec[-4] rec[-2]
    DETECTOR(0, 2.5, 0) rec[-124] rec[-122] rec[-120] rec[-118] rec[-116] rec[-114] rec[-112] rec[-110] rec[-108] rec[-106] rec[-39] rec[-37] rec[-35] rec[-33] rec[-31] rec[-29] rec[-27] rec[-25] rec[-23] rec[-21]
    DETECTOR(0, 3.5, 0) rec[-104] rec[-102] rec[-100] rec[-98] rec[-96] rec[-94] rec[-92] rec[-90] rec[-88] rec[-86] rec[-19] rec[-17] rec[-15] rec[-13] rec[-11] rec[-9] rec[-7] rec[-5] rec[-3] rec[-1]
    SHIFT_COORDS(0, 0, 1)
    DEPOLARIZE2(0.01) 1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49
    DEPOLARIZE1(0.01) 0 5 10 15 20 25 30 35 40 45
    TICK
}



DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 5 1 6 2 7 3 8 4 9 10 15 11 16 12 17 13 18 14 19 30 35 31 36 32 37 33 38 34 39 40 45 41 46 42 47 43 48 44 49 #q0+q1
DEPOLARIZE1("""+str(Loc_H[2][0])+""") 0 1 2 3 4 10 11 12 13 14 30 31 32 33 34 40 41 42 43 44 #q0
PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 20 21 22 23 24 25 26 27 28 29
PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 45 46 47 48 49 20 21 22 23 24 25 26 27 28 29
MPP("""+str(Read_err[2])+""") X0*X5 X1*X6 X2*X7 X3*X8 X4*X9 X10*X15 X11*X16 X12*X17 X13*X18 X14*X19 X30*X35 X31*X36 X32*X37 X33*X38 X34*X39 X40*X45 X41*X46 X42*X47 X43*X48 X44*X49
PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 45 46 47 48 49 20 21 22 23 24 25 26 27 28 29#q1
PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 45 46 47 48 49 20 21 22 23 24 25 26 27 28 29
PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 20 21 22 23 24 25 26 27 28 29
DEPOLARIZE1("""+str(Loc_H[2][0])+""") 0 1 2 3 4 10 11 12 13 14 30 31 32 33 34 40 41 42 43 44 #q0
DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 5 1 6 2 7 3 8 4 9 10 15 11 16 12 17 13 18 14 19 30 35 31 36 32 37 33 38 34 39 40 45 41 46 42 47 43 48 44 49 #q0+q1


TICK
DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 5 10 6 11 7 12 8 13 9 14 15 20 16 21 17 22 18 23 19 24 25 30 26 31 27 32 28 33 29 34 35 40 36 41 37 42 38 43 39 44
DEPOLARIZE1("""+str(Loc_H[2][0])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39
PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49
PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44
MPP("""+str(Read_err[2])+""") X5*X10 X6*X11 X7*X12 X8*X13 X9*X14 X15*X20 X16*X21 X17*X22 X18*X23 X19*X24 X25*X30 X26*X31 X27*X32 X28*X33 X29*X34 X35*X40 X36*X41 X37*X42 X38*X43 X39*X44
DEPOLARIZE1("""+str(Loc_H[2][0])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 #q0
DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 5 10 6 11 7 12 8 13 9 14 15 20 16 21 17 22 18 23 19 24 25 30 26 31 27 32 28 33 29 34 35 40 36 41 37 42 38 43 39 44 #q0+q1 
PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44 #q1+missing
PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44 #
PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 #missing #


TICK
DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48
PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 4 9 14 19 24 29 34 39 44 49 #missing 
MPP("""+str(Read_err[2])+""") Z0*Z1 Z2*Z3 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z15*Z16 Z17*Z18 Z20*Z21 Z22*Z23 Z25*Z26 Z27*Z28 Z30*Z31 Z32*Z33 Z35*Z36 Z37*Z38 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48
PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 0 2 5 7 10 12 15 17 20 22 25 27 30 32 35 37 40 42 45 47 4 9 14 19 24 29 34 39 44 49   #q0+missing
PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 4 9 14 19 24 29 34 39 44 49 #missing 
DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48


TICK
DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49
PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 5 10 15 20 25 30 35 40 45 #missing 
MPP("""+str(Read_err[2])+""") Z1*Z2 Z3*Z4 Z6*Z7 Z8*Z9 Z11*Z12 Z13*Z14 Z16*Z17 Z18*Z19 Z21*Z22 Z23*Z24 Z26*Z27 Z28*Z29 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49
PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 1 3 6 8 11 13 16 18 21 23 26 28 31 33 36 38 41 43 46 48 0 5 10 15 20 25 30 35 40 45   #q0+missing
PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 5 10 15 20 25 30 35 40 45 #missing 
DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49

DETECTOR(0.5, 0, 0) rec[-165] rec[-164] rec[-163] rec[-162] rec[-161] rec[-80] rec[-79] rec[-78] rec[-77] rec[-76]
DETECTOR(1.5, 0, 0) rec[-140] rec[-139] rec[-138] rec[-137] rec[-136] rec[-60] rec[-59] rec[-58] rec[-57] rec[-56]
DETECTOR(2.5, 0, 0) rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-75] rec[-74] rec[-73] rec[-72] rec[-71]
DETECTOR(3.5, 0, 0) rec[-135] rec[-134] rec[-133] rec[-132] rec[-131] rec[-55] rec[-54] rec[-53] rec[-52] rec[-51]
DETECTOR(5.5, 0, 0) rec[-130] rec[-129] rec[-128] rec[-127] rec[-126] rec[-50] rec[-49] rec[-48] rec[-47] rec[-46]
DETECTOR(6.5, 0, 0) rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-70] rec[-69] rec[-68] rec[-67] rec[-66]
DETECTOR(7.5, 0, 0) rec[-125] rec[-124] rec[-123] rec[-122] rec[-121] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41]
DETECTOR(8.5, 0, 0) rec[-145] rec[-144] rec[-143] rec[-142] rec[-141] rec[-65] rec[-64] rec[-63] rec[-62] rec[-61]
DETECTOR(0, 0.5, 0) rec[-120] rec[-118] rec[-116] rec[-114] rec[-112] rec[-40] rec[-38] rec[-36] rec[-34] rec[-32]
DETECTOR(0, 1.5, 0) rec[-100] rec[-98] rec[-96] rec[-94] rec[-92] rec[-20] rec[-18] rec[-16] rec[-14] rec[-12]
DETECTOR(0, 2.5, 0) rec[-119] rec[-117] rec[-115] rec[-113] rec[-111] rec[-39] rec[-37] rec[-35] rec[-33] rec[-31]
DETECTOR(0, 3.5, 0) rec[-99] rec[-97] rec[-95] rec[-93] rec[-91] rec[-19] rec[-17] rec[-15] rec[-13] rec[-11]
DETECTOR(5, 0.5, 0) rec[-110] rec[-108] rec[-106] rec[-104] rec[-102] rec[-30] rec[-28] rec[-26] rec[-24] rec[-22]
DETECTOR(5, 1.5, 0) rec[-90] rec[-88] rec[-86] rec[-84] rec[-82] rec[-10] rec[-8] rec[-6] rec[-4] rec[-2]
DETECTOR(5, 2.5, 0) rec[-109] rec[-107] rec[-105] rec[-103] rec[-101] rec[-29] rec[-27] rec[-25] rec[-23] rec[-21]
DETECTOR(5, 3.5, 0) rec[-89] rec[-87] rec[-85] rec[-83] rec[-81] rec[-9] rec[-7] rec[-5] rec[-3] rec[-1]
SHIFT_COORDS(0, 0, 1)

#Transport and storage of information
TICK

REPEAT """+str(j)+"""{
    

    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 5 1 6 2 7 3 8 4 9 10 15 11 16 12 17 13 18 14 19 30 35 31 36 32 37 33 38 34 39 40 45 41 46 42 47 43 48 44 49 #q0+q1
    DEPOLARIZE1("""+str(Loc_H[2][0])+""") 0 1 2 3 4 10 11 12 13 14 30 31 32 33 34 40 41 42 43 44 #q0
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 20 21 22 23 24 25 26 27 28 29
    PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 45 46 47 48 49 20 21 22 23 24 25 26 27 28 29
    MPP("""+str(Read_err[2])+""") X0*X5 X1*X6 X2*X7 X3*X8 X4*X9 X10*X15 X11*X16 X12*X17 X13*X18 X14*X19 X30*X35 X31*X36 X32*X37 X33*X38 X34*X39 X40*X45 X41*X46 X42*X47 X43*X48 X44*X49
    PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 45 46 47 48 49 20 21 22 23 24 25 26 27 28 29#q1
    PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 45 46 47 48 49 20 21 22 23 24 25 26 27 28 29
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 20 21 22 23 24 25 26 27 28 29
    DEPOLARIZE1("""+str(Loc_H[2][0])+""") 0 1 2 3 4 10 11 12 13 14 30 31 32 33 34 40 41 42 43 44 #q0
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 5 1 6 2 7 3 8 4 9 10 15 11 16 12 17 13 18 14 19 30 35 31 36 32 37 33 38 34 39 40 45 41 46 42 47 43 48 44 49 #q0+q1
    
    TICK
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 5 10 6 11 7 12 8 13 9 14 15 20 16 21 17 22 18 23 19 24 25 30 26 31 27 32 28 33 29 34 35 40 36 41 37 42 38 43 39 44
    DEPOLARIZE1("""+str(Loc_H[2][0])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49
    PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44
    MPP("""+str(Read_err[2])+""") X5*X10 X6*X11 X7*X12 X8*X13 X9*X14 X15*X20 X16*X21 X17*X22 X18*X23 X19*X24 X25*X30 X26*X31 X27*X32 X28*X33 X29*X34 X35*X40 X36*X41 X37*X42 X38*X43 X39*X44
    DEPOLARIZE1("""+str(Loc_H[2][0])+""") 5 6 7 8 9 15 16 17 18 19 25 26 27 28 29 35 36 37 38 39 #q0
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 5 10 6 11 7 12 8 13 9 14 15 20 16 21 17 22 18 23 19 24 25 30 26 31 27 32 28 33 29 34 35 40 36 41 37 42 38 43 39 44 #q0+q1 
    PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44 #q1+missing
    PAULI_CHANNEL_1("""+str(H_err[1][0])+""", """+str(H_err[1][1])+""", """+str(H_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 10 11 12 13 14 20 21 22 23 24 30 31 32 33 34 40 41 42 43 44 #
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 1 2 3 4 45 46 47 48 49 #missing #

    
    TICK
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 4 9 14 19 24 29 34 39 44 49 #missing 
    MPP("""+str(Read_err[2])+""") Z0*Z1 Z2*Z3 Z5*Z6 Z7*Z8 Z10*Z11 Z12*Z13 Z15*Z16 Z17*Z18 Z20*Z21 Z22*Z23 Z25*Z26 Z27*Z28 Z30*Z31 Z32*Z33 Z35*Z36 Z37*Z38 Z40*Z41 Z42*Z43 Z45*Z46 Z47*Z48
    PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 0 2 5 7 10 12 15 17 20 22 25 27 30 32 35 37 40 42 45 47 4 9 14 19 24 29 34 39 44 49   #q0+missing
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 4 9 14 19 24 29 34 39 44 49 #missing 
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 0 1 2 3 5 6 7 8 10 11 12 13 15 16 17 18 20 21 22 23 25 26 27 28 30 31 32 33 35 36 37 38 40 41 42 43 45 46 47 48

    TICK
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 5 10 15 20 25 30 35 40 45 #missing 
    MPP("""+str(Read_err[2])+""") Z1*Z2 Z3*Z4 Z6*Z7 Z8*Z9 Z11*Z12 Z13*Z14 Z16*Z17 Z18*Z19 Z21*Z22 Z23*Z24 Z26*Z27 Z28*Z29 Z31*Z32 Z33*Z34 Z36*Z37 Z38*Z39 Z41*Z42 Z43*Z44 Z46*Z47 Z48*Z49
    PAULI_CHANNEL_1("""+str(M_err[1][0])+""", """+str(M_err[1][1])+""", """+str(M_err[1][2])+""") 1 3 6 8 11 13 16 18 21 23 26 28 31 33 36 38 41 43 46 48 0 5 10 15 20 25 30 35 40 45   #q0+missing
    PAULI_CHANNEL_1("""+str(CX_err[1][0])+""", """+str(CX_err[1][1])+""", """+str(CX_err[1][2])+""") 0 5 10 15 20 25 30 35 40 45 #missing 
    DEPOLARIZE2("""+str(Loc_CNOT[2])+""") 1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19 21 22 23 24 26 27 28 29 31 32 33 34 36 37 38 39 41 42 43 44 46 47 48 49

    
    DETECTOR(0.5, 0, 0) rec[-160] rec[-159] rec[-158] rec[-157] rec[-156] rec[-80] rec[-79] rec[-78] rec[-77] rec[-76]#
    DETECTOR(1.5, 0, 0) rec[-140] rec[-139] rec[-138] rec[-137] rec[-136] rec[-60] rec[-59] rec[-58] rec[-57] rec[-56]#
    DETECTOR(2.5, 0, 0) rec[-155] rec[-154] rec[-153] rec[-152] rec[-151] rec[-75] rec[-74] rec[-73] rec[-72] rec[-71]#
    DETECTOR(3.5, 0, 0) rec[-135] rec[-134] rec[-133] rec[-132] rec[-131] rec[-55] rec[-54] rec[-53] rec[-52] rec[-51]#
    DETECTOR(5.5, 0, 0) rec[-130] rec[-129] rec[-128] rec[-127] rec[-126] rec[-50] rec[-49] rec[-48] rec[-47] rec[-46]#
    DETECTOR(6.5, 0, 0) rec[-150] rec[-149] rec[-148] rec[-147] rec[-146] rec[-70] rec[-69] rec[-68] rec[-67] rec[-66]#
    DETECTOR(7.5, 0, 0) rec[-125] rec[-124] rec[-123] rec[-122] rec[-121] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41]#
    DETECTOR(8.5, 0, 0) rec[-145] rec[-144] rec[-143] rec[-142] rec[-141] rec[-65] rec[-64] rec[-63] rec[-62] rec[-61]#
    DETECTOR(0, 0.5, 0) rec[-120] rec[-118] rec[-116] rec[-114] rec[-112] rec[-40] rec[-38] rec[-36] rec[-34] rec[-32]#
    DETECTOR(0, 1.5, 0) rec[-100] rec[-98] rec[-96] rec[-94] rec[-92] rec[-20] rec[-18] rec[-16] rec[-14] rec[-12]#
    DETECTOR(0, 2.5, 0) rec[-119] rec[-117] rec[-115] rec[-113] rec[-111] rec[-39] rec[-37] rec[-35] rec[-33] rec[-31]#
    DETECTOR(0, 3.5, 0) rec[-99] rec[-97] rec[-95] rec[-93] rec[-91] rec[-19] rec[-17] rec[-15] rec[-13] rec[-11]#
    DETECTOR(5, 0.5, 0) rec[-110] rec[-108] rec[-106] rec[-104] rec[-102] rec[-30] rec[-28] rec[-26] rec[-24] rec[-22]#
    DETECTOR(5, 1.5, 0) rec[-90] rec[-88] rec[-86] rec[-84] rec[-82] rec[-10] rec[-8] rec[-6] rec[-4] rec[-2]#
    DETECTOR(5, 2.5, 0) rec[-109] rec[-107] rec[-105] rec[-103] rec[-101] rec[-29] rec[-27] rec[-25] rec[-23] rec[-21]#
    DETECTOR(5, 3.5, 0) rec[-89] rec[-87] rec[-85] rec[-83] rec[-81] rec[-9] rec[-7] rec[-5] rec[-3] rec[-1]#

}

TICK
MX("""+str(Read_err[2])+""") 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
DETECTOR(0.5, 0, 0) rec[-130] rec[-129] rec[-128] rec[-127] rec[-126] rec[-50] rec[-49] rec[-48] rec[-47] rec[-46] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41]
DETECTOR(1.5, 0, 0) rec[-110] rec[-109] rec[-108] rec[-107] rec[-106] rec[-45] rec[-44] rec[-43] rec[-42] rec[-41] rec[-40] rec[-39] rec[-38] rec[-37] rec[-36]
DETECTOR(2.5, 0, 0) rec[-125] rec[-124] rec[-123] rec[-122] rec[-121] rec[-40] rec[-39] rec[-38] rec[-37] rec[-36] rec[-35] rec[-34] rec[-33] rec[-32] rec[-31]
DETECTOR(3.5, 0, 0) rec[-105] rec[-104] rec[-103] rec[-102] rec[-101] rec[-35] rec[-34] rec[-33] rec[-32] rec[-31] rec[-30] rec[-29] rec[-28] rec[-27] rec[-26]
DETECTOR(5.5, 0, 0) rec[-100] rec[-99] rec[-98] rec[-97] rec[-96] rec[-25] rec[-24] rec[-23] rec[-22] rec[-21] rec[-20] rec[-19] rec[-18] rec[-17] rec[-16]
DETECTOR(6.5, 0, 0) rec[-120] rec[-119] rec[-118] rec[-117] rec[-116] rec[-20] rec[-19] rec[-18] rec[-17] rec[-16] rec[-15] rec[-14] rec[-13] rec[-12] rec[-11]
DETECTOR(7.5, 0, 0) rec[-95] rec[-94] rec[-93] rec[-92] rec[-91] rec[-15] rec[-14] rec[-13] rec[-12] rec[-11] rec[-10] rec[-9] rec[-8] rec[-7] rec[-6]
DETECTOR(8.5, 0, 0) rec[-115] rec[-114] rec[-113] rec[-112] rec[-111] rec[-10] rec[-9] rec[-8] rec[-7] rec[-6] rec[-5] rec[-4] rec[-3] rec[-2] rec[-1]
OBSERVABLE_INCLUDE(0) rec[-30] rec[-29] rec[-28] rec[-27] rec[-26]
OBSERVABLE_INCLUDE(1) rec[-25] rec[-24] rec[-23] rec[-22] rec[-21]
               
""")
    
        

        print("----------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        print(j)
        print(t0,t1,t2,t3)

        tc0.append(t0)
        tc1.append(t1)
        tc2.append(t2)
        tc3.append(t3)



        
        # Verification if there are detectors missing or something wrong going on:--------------------------
        # Check conditions for each code:

        print("Count detector sizes, expected and true :=====================================================")


        a=count_determined_measurements_in_circuit(circuit1)
        expected_determined = circuit1.num_detectors + circuit1.num_observables
        print(expected_determined,a) #expected, true

        a=count_determined_measurements_in_circuit(circsurf_notr)
        expected_determined = circsurf_notr.num_detectors + circsurf_notr.num_observables
        print(expected_determined,a) #expected, true

        a=count_determined_measurements_in_circuit(circsurf_notr2)
        expected_determined = circsurf_notr2.num_detectors + circsurf_notr2.num_observables
        print(expected_determined,a) #expected, true


        sampler2 = circuit2.compile_detector_sampler()
        S_Bell.append(np.sum(sampler2.sample(shots=10**6)) / 10**6)


        print("Distances of error :=====================================================")
        actual_distance1 = len(circsurf_notr.shortest_graphlike_error())
        actual_distance2 = len(circsurf_notr2.shortest_graphlike_error())
        actual_distance3 = len(circuit1.shortest_graphlike_error())


        print("distance 3?")
        print(actual_distance1,actual_distance2,actual_distance3)
        if (actual_distance1 and actual_distance2 and actual_distance3 !=3):
            raise ValueError('Distance !=3, increase/decrease lattice dimensions')


        print(" ------------------------------------------------------------------------ ")
        print("End of verification...")


        surface_code_tasks = [
            sinter.Task(
            circuit = circuit1,
            json_metadata={'d': 4, 'p': 0.01},
            )#,


      #  sinter.Task(
       # circuit = circuit1,
        #json_metadata={'d': 3, 'r': 3 * 3, 'p': 0.1},
    #)
    
        ]

        surface_code_tasks2 = [
            sinter.Task(
            circuit = circsurf_notr,
            json_metadata={'d': 4, 'p': 0.01},
            )#,


      #  sinter.Task(
       # circuit = circuit1,
        #json_metadata={'d': 3, 'r': 3 * 3, 'p': 0.1},
    #)
    
        ]

        surface_code_tasks3 = [
            sinter.Task(
            circuit = circsurf_notr2,
            json_metadata={'d': 4, 'p': 0.01},
            )#,
            

      #  sinter.Task(
       # circuit = circuit1,
        #json_metadata={'d': 3, 'r': 3 * 3, 'p': 0.1},
    #)
    
        ]
        #print("TASK--------------------------------------------------------------------")
        #print(surface_code_tasks)
#Surface code statistics associated with monte carlo sampling using the decoder to the errors per iteration
        collected_stats: List[sinter.TaskStats] = sinter.collect(
        num_workers=4,
        tasks=surface_code_tasks,
        decoders=['pymatching'],
        max_shots=1_000_000,
    #max_errors=7_000,
        print_progress=True,
        count_detection_events=True



            )

        collected_stats2: List[sinter.TaskStats] = sinter.collect(
        num_workers=4,
        tasks=surface_code_tasks2,
        decoders=['pymatching'],
        max_shots=1_000_000,
    #max_errors=7_000,
        print_progress=True,
        count_detection_events=True


        
            )

        collected_stats3: List[sinter.TaskStats] = sinter.collect(
        num_workers=4,
        tasks=surface_code_tasks3,
        decoders=['pymatching'],
        max_shots=1_000_000,
    #max_errors=7_000,
        print_progress=True,
        count_detection_events=True
            )

        print("Iteration number "+str(j))
        print("Metadata----------------------------------------------")
        print(sinter.CSV_HEADER)
        for sample in collected_stats:
            print(sample.to_csv_line())
        
        X,Y,Y_l,Y_h=plot_error_rate2(
    stats=collected_stats,
    x_func=lambda stat: stat.json_metadata['p'],
    group_func=lambda stat: stat.json_metadata['d'],
    #failure_units_per_shot_func=lambda stat: stat.json_metadata['r'],
        )   

        Vec_Y.append(Y[0])
        Std_p.append(Y_h[0])
        Std_m.append(Y_l[0])


        X2,Y2,Y_l2,Y_h2=plot_error_rate2(
    stats=collected_stats2,
    x_func=lambda stat: stat.json_metadata['p'],
    group_func=lambda stat: stat.json_metadata['d'],
    #failure_units_per_shot_func=lambda stat: stat.json_metadata['r'],
        )   

        Vec_Y2.append(Y2[0])
        Std_p2.append(Y_h2[0])
        Std_m2.append(Y_l2[0])



        X3,Y3,Y_l3,Y_h3=plot_error_rate2(
    stats=collected_stats3,
    x_func=lambda stat: stat.json_metadata['p'],
    group_func=lambda stat: stat.json_metadata['d'],
    #failure_units_per_shot_func=lambda stat: stat.json_metadata['r'],
        )   

        Vec_Y3.append(Y3[0])
        Std_p3.append(Y_h3[0])
        Std_m3.append(Y_l3[0])



     
    print("Data output parameters")
    print(" ")
    print(tc0,tc1,tc2,tc3)
    print(" ")
    print(S_Bell)
    print(" ")
    print(Vec_Y2)
    print(" ")
    print(Std_p,Std_p2,Std_p3)
    print(" ")
    print(Std_m,Std_m2,Std_m3)
    print(" ")
    print(Vec_Y,Vec_Y2,Vec_Y3)



    ####################################################################
    #plot without offset, all codes should start QEC cycles of the storage at t=0 (this is included in the plot.py file)

    fig, ax = plt.subplots(1, 1)
    ax.fill_between(tc1, Std_m, Std_p,alpha=0.2,color="black")
    ax.plot(tc1,Vec_Y,"v",color="black")
    ax.plot(tc1,Vec_Y,color="black",label="BS [[18,2,3]]")

    ax.fill_between(tc2, Std_m2, Std_p2,alpha=0.2,color="blue")
    ax.plot(tc2,Vec_Y2,"v",color="blue")
    ax.plot(tc2,Vec_Y2,color="blue",label="S [[18,2,3]]")

    ax.fill_between(tc3, Std_m3, Std_p3,alpha=0.2,color="green")
    ax.plot(tc3,Vec_Y3,"v",color="green")
    ax.plot(tc3,Vec_Y3,color="green",label="rotated S [[18,2,3]]")

    ax.plot(tc0,S_Bell,"v",color="red")
    ax.plot(tc0,S_Bell,color="red",label="Unencoded Bell-State")

    ax.plot(tc0c,S_Bell2,"v",alpha=0.3,color="red")
    ax.plot(tc0c,S_Bell2,linestyle='--',alpha=0.3,color="red")

    ax.grid(linestyle='--')
    plt.legend(loc="lower right")
    ax.set_xlabel("Processing time [$\mu$s]")
    ax.set_ylabel("Logical Error Rate")
    plt.show()

