# coding: utf-8

# In[492]:

from operator import mul
import itertools
import collections
from functools import reduce


# In[493]:

def combinations_withmax2repeats(items, n):
    if n==0: yield []
    else:
        for i in range(len(items)):
            for cc in combinations_withmax2repeats(items[i:],n-1): # repetition of elements is allowed
                if len(cc)>1 and cc[0] == cc[1] and cc[0] == items[i]: # when an element is repeated more than 2 times, this combination is skipped
                    continue
                yield [items[i]]+cc                    


# In[494]:

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


# In[495]:

def PE1Lon(observed,nonobserved):
    """
    Function calculates the PEL0 for each locus and puts the results in list PE1L_list
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


# In[496]:

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


# In[497]:

def calculate_combinations(combination,PE0Ls,PE1Ls,PE2Ls,PE0Lsreduced):
    """
    Calculate the RMNE value based on the observed and nonobserved alleles and a certain combination of drop-outs
    
    This function takes as input:
     - A combination consisting of a list of loci where drop-outs have occurred.
       When multiple drop-outs have occurred in the same locus, the number of the locus is multiple times in this list
     - A list with the PE0Ls from each locus
     - A list with the PE1Ls from each locus
     - A list with the PE2Ls from each locus
     - PE0Lsreduced: the RMNE calculation as if there are no drop-outs 
    """
    #PE=0
    loci_with_2dropouts=[] # list of loci where exact 2 drop-outs have occurred considering the combination
    loci_with_1dropout=[] # list of loci where exact 2 drop-outs have occurred considering the combination
    
    for x, y in collections.Counter(combination).items(): # for each locus with at least a drop-out, determine how many drop-outs are in the locus
        if y==2: # when there are exact 2 drop-outs in a locus, put this locus in the loci_with_2dropouts list
            loci_with_2dropouts.append(x)
        elif y==1: # when there are exact 1 drop-outs in a locus, put this locus in the loci_with_1dropout list
            loci_with_1dropout.append(x)
    
    # in stead of calculating the RMNE value for the profile and the drop-outs in the combination, 
    # the RMNE value with no drop-outs is used as a basis value which is corrected for each locus with one or 2 drop-outs
    PEc=PE0Lsreduced
    for c in loci_with_1dropout:
            PEc=(PEc/PE0Ls[int(c)])*PE1Ls[int(c)]
    for c in loci_with_2dropouts: 
            PEc=(PEc/PE0Ls[int(c)])*PE2Ls[int(c)]
    #PE+=PEc
    return(PEc)
    


# In[498]:

def PE(PE0Ls,PE1Ls,PE2Ls,DO):
    """
    Generate all possible combinations of exactly DO drop-outs in the number of loci. 
    The number of loci is determined based on the length of the PE0Ls list, which contains the PE0Ls from each locus.
        
    This function takes as input:
     - A list with the PE0Ls from each locus
     - A list with the PE1Ls from each locus
     - A list with the PE2Ls from each locus
     The loci need to be in the same order in the lists 
     - The number of drop-outs
    """
    
    n_loci=list(range(0,len(PE0Ls))) # make a list with a number for each locus in the analysis based on the length of the PE0Ls list
    
    #calculate PE as if there are no drop-outs
    PE0Lsreduced=reduce(mul, PE0Ls, 1) # Caculate the RMNE value with 0 drop-outs
    
    if DO==0:
        return (DO,PE0Lsreduced) # when there are no drop-outs, return the RMNE value for 0 drop-outs
    else: 
        # When there are drop-outs, make all combinations of DO drop-outs in the number of loci. 
        # All combinations are generated of DO elements from n_loci elements. 
        # Each element in the list can be used two times in a combination.
        combinations=combinations_withmax2repeats(n_loci,DO)
    
    #For each combination: calculate the PE value of that combination. The final PE value is the sum of PE values from all these combinations
    PE=0
    for combination in combinations:
        PE+=calculate_combinations(combination,PE0Ls,PE1Ls,PE2Ls,PE0Lsreduced)
            
    return (DO,PE)


# In[499]:

def PElist(PE0Ls,PE1Ls,PE2Ls,allowedDO):
    """
    This function calculates the PE (see PE function) for a range of drop-outs starting from 0 to the number of allowedDO
    Returns a list of tuples with the number of drop-outs and the corresponding PE value
    """
    values=[]
    for i in range(0,allowedDO+1):
        values.append(PE(PE0Ls,PE1Ls,PE2Ls,i))
    
    return (values)


# In[503]:

def PE_star(x):
    return PE(*x)


# In[504]:

def PElist_multi_proc(PE0Ls,PE1Ls,PE2Ls,allowedDO):
    from multiprocessing import Pool
    """
    This function calculates the PE (see PE function) for a range of drop-outs starting from 0 to the number of allowedDO
    Returns a list of tuples with the number of drop-outs and the corresponding PE value.
    Each PE is caclulated in a separate process.
    """
    PE_list=[]
    pool = Pool()
    changing_arg=range(0,allowedDO+1)
    values=pool.map(PE_star, zip(itertools.repeat(PE0Ls),itertools.repeat(PE1Ls),itertools.repeat(PE2Ls),changing_arg))
    pool.close()
    pool.join()
    return (values)


# In[490]:

#Main program (example)
#observed=[[0.24762, 0.23452], [0.21548, 0.18929], [0.15476, 0.09286], [0.12143000000000001, 0.05595], [0.19405000000000003, 0.006], [0.07024, 0.23095000000000002], [0.19524000000000002, 0.16667], [0.0117, 0.006], [0.12738000000000002, 0.31667], [0.28300000000000003, 0.17], [0.059], [0.2, 0.12], [0.28600000000000003], [0.371, 0.006]]
#nonobserved=[[0.006, 0.006, 0.14643, 0.20476000000000003, 0.14048000000000002, 0.01429, 0.006, 0.006], [0.006, 0.006, 0.12262, 0.14524, 0.3119, 0.010710000000000003, 0.006, 0.006, 0.006], [0.006, 0.006, 0.006, 0.006, 0.006, 0.02857, 0.006, 0.2119, 0.006, 0.2631, 0.027380000000000005, 0.08333, 0.01667, 0.07262, 0.006, 0.03214, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006], [0.006, 0.006, 0.00833, 0.006, 0.010710000000000003, 0.12857000000000002, 0.15476, 0.006, 0.15595, 0.15, 0.13214, 0.04048, 0.01667, 0.00952, 0.00833, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006], [0.006, 0.006, 0.006, 0.006, 0.1, 0.11071, 0.27857000000000004, 0.2131, 0.08452, 0.0131, 0.006, 0.006], [0.006, 0.02619, 0.00714, 0.07738, 0.13690000000000002, 0.32738, 0.1, 0.02381, 0.006, 0.006, 0.006], [0.006, 0.006, 0.01667, 0.006, 0.04762, 0.006, 0.13214, 0.006, 0.17143000000000003, 0.006, 0.00952, 0.006, 0.13452, 0.006, 0.07619, 0.006, 0.0369, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006], [0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.037, 0.006, 0.049, 0.006, 0.039, 0.006, 0.063, 0.006, 0.073, 0.006, 0.059, 0.006, 0.054000000000000006, 0.022, 0.049, 0.015, 0.017, 0.024, 0.006, 0.037, 0.034, 0.034, 0.046, 0.078, 0.083, 0.073, 0.0386, 0.027000000000000003, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006], [0.006, 0.006, 0.006, 0.006, 0.01905, 0.044050000000000006, 0.29762000000000005, 0.16667, 0.027380000000000005, 0.006, 0.006, 0.006], [0.006, 0.006, 0.006, 0.006, 0.027999999999999997, 0.34400000000000003, 0.14300000000000002, 0.027999999999999997, 0.006, 0.006, 0.006], [0.006, 0.006, 0.08900000000000001, 0.13, 0.053, 0.091, 0.006, 0.135, 0.066, 0.125, 0.051, 0.055999999999999994, 0.12, 0.006, 0.006, 0.012999999999999998, 0.006, 0.006], [0.006, 0.046, 0.025999999999999995, 0.085, 0.015, 0.021, 0.126, 0.015, 0.006, 0.12, 0.006, 0.088, 0.078, 0.04, 0.016, 0.006, 0.006, 0.006], [0.006, 0.006, 0.17, 0.36100000000000004, 0.044, 0.051, 0.006, 0.029000000000000005, 0.05, 0.006, 0.006, 0.006], [0.006, 0.006, 0.006, 0.006, 0.139, 0.006, 0.006, 0.05, 0.34800000000000003, 0.073, 0.009000000000000001, 0.006, 0.006]]
#p1=PE0Lon(observed,nonobserved)
#p2=PE1Lon(observed,nonobserved)
#p3=PE2Lon(observed,nonobserved)
#import time
#for i in range(0,29):
#    t = time.process_time()
#    v=PE(p1,p2,p3,i)
#    RMNE+=v[-1]
#    elapsed_time = time.process_time() - t
#    print(i,v,RMNE,elapsed_time)

