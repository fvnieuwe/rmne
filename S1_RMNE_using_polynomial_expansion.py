# RMNE calculation using polynomial expansion

from sympy import symbols as symbols
from sympy import expand as expand
from functools import reduce
from operator import mul

def RMNE(p1,p2,p3):
    x = symbols('x')
    expr = reduce(mul, [(i+j*x+k*x**2) for i,j,k in zip(p1,p2,p3)], 1)
    expr = expand(expr)

    sympyPE=[]
    for r in range(len(p1) * 2 + 1):
        sympyPE.append(expr.coeff(x, r))
    return sympyPE

def PE0Lon(observed_list,nonobserved_list):
    """
    Function calculates the PEL0 for each locus and puts the results in list PE0L_list
    PEL0 is the RMNE value for the locus assuming that exactly 0 drop-outs have occurred.
    The function takes as input 2 list of lists: 
     -list of the frequencies in the population of the observed alleles for each locus 
       ([[freq. of observed alleles in locus 1],[freq. of observed in locus 2],...])
     -list of the frequencies in the population of the nonobserved alleles for each locus (is not actually used in the calculation)
       ([[freq. of nonobserved alleles in locus 1],[freq. of nonobserved in locus 2],...])
     the loci need to be in the same order on both lists of lists
    """
    PE0L_list=[]
    for listperlocus_o,listperlocus_n in zip(observed_list,nonobserved_list):
        PE0L=(sum(listperlocus_o))**2
        if PE0L==0: # when there are no observed alleles, the locus is ignored and PE0L is 1
            PE0L=1
        PE0L_list.append(PE0L)
    return PE0L_list


def PE1Lon(observed,nonobserved):
    """
    Function calculates the PEL1 for each locus and puts the results in list PE1L_list
    PEL1 is the RMNE value for the locus assuming that exactly 1 drop-out has occurred.
    The function takes as input 2 list of lists: 
     -list of the frequencies in the population of the observed alleles for each locus 
       ([[freq. of observed alleles in locus 1],[freq. of observed in locus 2],...])
     -list of the frequencies in the population of the nonobserved alleles for each locus
       ([[freq. of nonobserved alleles in locus 1],[freq. of nonobserved in locus 2],...])
     the loci need to be in the same order on both lists
    """
    PE1L_list=[]
    for listperlocus_o,listperlocus_n in zip(observed,nonobserved):
        PE1L=(sum(listperlocus_o)*sum(listperlocus_n)*2)
        PE1L_list.append(PE1L)
    return PE1L_list


def PE2Lon(observed,nonobserved):
    """
    Function calculates the PEL2 for each locus and puts the results in list PE2L_list
    PEL2 is the RMNE value for the locus assuming that exactly 2 drop-outs have occurred.
    The function takes as input 2 list of lists: 
     -list of the frequencies in the population of the observed alleles for each locus 
       ([[freq. of observed alleles in locus 1],[freq. of observed in locus 2],...])
     -list of the frequencies in the population of the nonobserved alleles for each locus
       ([[freq. of nonobserved alleles in locus 1],[freq. of nonobserved in locus 2],...])
     the loci need to be in the same order on both lists o
    """
    PE2L_list=[]
    for listperlocus_o,listperlocus_n in zip(observed,nonobserved):
        PE2L=(sum(listperlocus_n))**2
        so=sum(listperlocus_o)
        if so==0: # when there are no observed alleles, the locus is ignored and PE2L is 0
            PE2L=0
        PE2L_list.append(PE2L)
    return PE2L_list

# Profile data
observed = [
    [0.24762, 0.23452],
    [0.21548, 0.18929],
    [0.15476, 0.09286],
    [0.12143000000000001, 0.05595],
    [0.19405000000000003, 0.006],
    [0.07024, 0.23095000000000002],
    [0.19524000000000002, 0.16667],
    [0.0117, 0.006],
    [0.12738000000000002, 0.31667],
    [0.28300000000000003, 0.17],
    [0.059],
    [0.2, 0.12],
    [0.28600000000000003],
    [0.371, 0.006]
]
nonobserved = [
    [0.006, 0.006, 0.14643, 0.20476000000000003, 0.14048000000000002, 0.01429, 0.006, 0.006],
    [0.006, 0.006, 0.12262, 0.14524, 0.3119, 0.010710000000000003, 0.006, 0.006, 0.006],
    [0.006, 0.006, 0.006, 0.006, 0.006, 0.02857, 0.006, 0.2119, 0.006, 0.2631, 0.027380000000000005, 0.08333, 0.01667, 0.07262, 0.006, 0.03214, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006],
    [0.006, 0.006, 0.00833, 0.006, 0.010710000000000003, 0.12857000000000002, 0.15476, 0.006, 0.15595, 0.15, 0.13214, 0.04048, 0.01667, 0.00952, 0.00833, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006],
    [0.006, 0.006, 0.006, 0.006, 0.1, 0.11071, 0.27857000000000004, 0.2131, 0.08452, 0.0131, 0.006, 0.006],
    [0.006, 0.02619, 0.00714, 0.07738, 0.13690000000000002, 0.32738, 0.1, 0.02381, 0.006, 0.006, 0.006],
    [0.006, 0.006, 0.01667, 0.006, 0.04762, 0.006, 0.13214, 0.006, 0.17143000000000003, 0.006, 0.00952, 0.006, 0.13452, 0.006, 0.07619, 0.006, 0.0369, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006],
    [0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.037, 0.006, 0.049, 0.006, 0.039, 0.006, 0.063, 0.006, 0.073, 0.006, 0.059, 0.006, 0.054000000000000006, 0.022, 0.049, 0.015, 0.017, 0.024, 0.006, 0.037, 0.034, 0.034, 0.046, 0.078, 0.083, 0.073, 0.0386, 0.027000000000000003, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006],
    [0.006, 0.006, 0.006, 0.006, 0.01905, 0.044050000000000006, 0.29762000000000005, 0.16667, 0.027380000000000005, 0.006, 0.006, 0.006],
    [0.006, 0.006, 0.006, 0.006, 0.027999999999999997, 0.34400000000000003, 0.14300000000000002, 0.027999999999999997, 0.006, 0.006, 0.006],
    [0.006, 0.006, 0.08900000000000001, 0.13, 0.053, 0.091, 0.006, 0.135, 0.066, 0.125, 0.051, 0.055999999999999994, 0.12, 0.006, 0.006, 0.012999999999999998, 0.006, 0.006],
    [0.006, 0.046, 0.025999999999999995, 0.085, 0.015, 0.021, 0.126, 0.015, 0.006, 0.12, 0.006, 0.088, 0.078, 0.04, 0.016, 0.006, 0.006, 0.006],
    [0.006, 0.006, 0.17, 0.36100000000000004, 0.044, 0.051, 0.006, 0.029000000000000005, 0.05, 0.006, 0.006, 0.006],
    [0.006, 0.006, 0.006, 0.006, 0.139, 0.006, 0.006, 0.05, 0.34800000000000003, 0.073, 0.009000000000000001, 0.006, 0.006]
]

# Main program (example)
p1=PE0Lon(observed,nonobserved)
p2=PE1Lon(observed,nonobserved)
p3=PE2Lon(observed,nonobserved)

import time
t = time.process_time()
results=RMNE(p1,p2,p3)
elapsed_time = time.process_time() - t

for result in enumerate(results):
    print(result[0],result[1])
    
print("\n")
print("Total time: "+str(elapsed_time))