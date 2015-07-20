#----------------------------------------------------------#
# MarkovEqnGenSimple.py
# Authors: Stuart Campbell and Fred Lionetti
# 
# The following python code generates MATLAB functions 
# dynamically which contain vectors of coefficients and other
# data structures which ultimately desribe the system of ODEs
# representing the Markov model of myofilament activation and 
# force production.  
#
# This version is 'simple' because it eliminates explicit depiction
# of XB binding, and instead assumes that this does not directly
# alter activation of neighboring RUs.  The result is a drastic reduction
# in the total number of Markov states
#----------------------------------------------------------#



from string import *
import copy
import time
import numpy


def tobase(base,number):
    def tb(b,n,result=''):
        if n == 0: return result
        else: return tb(b,n/b,str(n%b)+result)
    if number == 0:
        return '0'
    if number > 0:
        return tb(base, number)
    if number < 0:
        return '-' + tb(base, -1*number)

def findall(mylist, key):
    found = []
    lastfound = -1
    try:
        while True:
            found.append(mylist.index(key, lastfound+1))
            lastfound = found[-1]
    except ValueError:
        return found

def unique_name(name):
    rname = copy.copy(name)
    rname.reverse()

    num_0 = name.count(0)
    num_1 = name.count(1)
    num_2 = name.count(2)
    potential_pivots = [x for x in (num_0, num_1, num_2) if x != 0]
    min_pivots = min(potential_pivots)
    if num_0 == min_pivots:
        pivots = findall(name, 0)
        rpivots = findall(rname, 0)
    elif num_1 == min_pivots:
        pivots = findall(name, 1)
        rpivots = findall(rname, 1)
    elif num_2 == min_pivots:
        pivots = findall(name, 2)
        rpivots = findall(rname, 2)
    
    variations = []
    for pivot in pivots:
        new_variation = name[pivot:] + name[:pivot]
        if new_variation in variations:
            continue
        variations.append(new_variation)
    for pivot in rpivots:
        new_variation = rname[pivot:] + rname[:pivot]
        if new_variation in variations:
            continue
        variations.append(new_variation)
    return min(variations)

def pivotMask(mask):
    #rotate mask so it starts with a 1
    pivot = mask.index(1)
    mask = mask[pivot:] + mask[:pivot]
    return mask

def matches(mask, target):
    target_sum = sum(target)
    if target_sum < sum(mask):
        return False
    if target_sum == len(target):
        return True
    #print 'hi'

    name = target
    rname = copy.copy(name)
    rname.reverse()

    pivots = findall(name, 1)
    rpivots = findall(rname, 1)

    
    variations = []
    for pivot in pivots:
        new_variation = name[pivot:] + name[:pivot]
        if new_variation in variations:
            continue
        variations.append(new_variation)
    for pivot in rpivots:
        new_variation = rname[pivot:] + rname[:pivot]
        if new_variation in variations:
            continue
        variations.append(new_variation)

    for variation in variations:
        a = numpy.array(mask, numpy.bool)
        b = numpy.array(variation, numpy.bool)
        if ((a&b)==a).all():
            return True
    return False

def collapseStateID(a):         # Turn 2's in a state.id into 1's (for k_tr experiments)
    b = ''
    for i, val in enumerate(a):
        if val == '2': 
            b = b + '1'
        else:
            b = b + a[i]
    return b
    
def collapseRU(a):         # Turn 1's in a state.id into 0's, and 3's into 1's (for NEM-S1)
    b = ''
    for i, val in enumerate(a):
        if val == '2': 
            b = b + '1'
        elif val == '1':
            b = b + '0'
        else:
            b = b + a[i]
    return b
    
def diffMasks(mask, id):            # Checks to see if id is contained in the set of states where at least all the 'on' RUs in mask are also 'on'
    id = collapseStateID(id)
    sum = 0
    for i in mask: sum = sum + i
    if not sum: return True         # When a mask is all zeros, this function always returns true
    for x in range(len(mask)):       
        if mask[x] and not id[x]:
            return True
    return False
    
def listOfListsIndex2MATLAB(a):         # Adds 1 to each entry of a list of lists to reconcile to MATLAB indexing
    for i, val in enumerate(a):
        for j, val2 in enumerate(val):
            a[i][j] = val2 + 1
    return a
    
def index2MATLAB(a):                      # Adds 1 to each entry of a list to reconcile to MATLAB indexing
    b = [0]*len(a)
    for i, val in enumerate(a):
        b[i] = val + 1
    return b
    
def getTransType(curid, newVal, i):          # Create transtype based on neighbors and newVal
    if i+1 == len(curid): b = atoi(curid[0])
    else: b = atoi(curid[i+1])

    if i == 0: a = atoi(curid[-1])
    else: a = atoi(curid[i-1])
    
    if a == 2: a = a+1
    if b == 2: b = b+1

    transType = a + b
    
    if newVal == '2':                 #If this is an XB binding transition,
        transType = transType + 10  #Offset by 10 to indicate an XB transType (not RU)
    
    return transType
    
def getNewVal(id1, id2):                 # Given two state id, determines what i and newVal for the transition are
    RU = [i for i, ch in enumerate(id1) if ch != id2[i]]
    newVal = id2[RU[0]]     # RU should only ever have one entry
    return RU[0], newVal
        
def getfileptr(varname, numRUs):
    filename = "EqnData/%s%d.txt"%(varname,numRUs)
    f = open(filename, 'w')
    return f
    
def getfileptrmdex(varname, numRUs, mdex):
    filename = "EqnData/%s%d_mdex%d.txt"%(varname,numRUs,mdex)
    f = open(filename, 'w')
    return f


class Transition:
    def __init__(self):
        self.m = 1
    def __str__(self):
        return "<transition: m = %d, n = %d, iPre = %d, iPost = %d, transType = %d>"%(self.m, self.n, self.iPre, self.iPost, self.transType)
    def __repr__(self):
        return self.__str__()

class State:
    def __init__(self, id):
        # precondition: id is the unique key, and is a string
        self.id = id
        self.alpha = len([x for x in id if x != '0']) / float(len(id))
        self.beta = len([x for x in id if x == '2']) / float(len(id))
        
    def __str__(self):
        return "<state: %s, alpha=%.2f, beta=%.2f>"%(self.id, self.alpha, self.beta)

    def __repr__(self):
        return self.__str__()

class MarkovChainGenerator:
    def __init__(self, numRUs):
        self.numRUs = numRUs
        self.transitions = {}
        self.allstates = ()
        
        #initialize root
        self.initStuff(numRUs)
        #root = [0,]*numRUs # State((0,)*numRUs)
        root = '0'*numRUs
        self.all_states = [root]
        self.j = 0  # Initialize transition indexer
        self.addStates(root, 0, '0', '1')
        #self.addStates(root, 0, '1', '2')

        forward_transitions = copy.deepcopy(self.transitions)
        forward_allstates = copy.deepcopy(self.all_states)

        print "finished forward state generation..."
        #print "Number of transitions: ", len(self.transitions)

        #reverse states
        self.transitions = {}    # reset for reverse pass
        root = '1'*numRUs
        self.all_states = [root]            # reset for reverse pass
        self.j = 0  # Initialize transition indexer
        #self.addStates(root, 0, '2', '1', True)
        self.addStates(root, 0, '1', '0', True)
        reverse_transitions = copy.deepcopy(self.transitions)
        reverse_allstates = copy.deepcopy(self.all_states)
        
        print "finished reverse state generation..."
        
        # print "len(states):", len(self.all_states)

        # Place reverse redundancies in the final transitions structure.
        for t in forward_transitions:
            trans = forward_transitions[t]
            rtrans = reverse_transitions[(t[1], t[0])]
            trans.n = rtrans.m
            # print t, "-->", trans

        print "num states:", len(reverse_allstates)
        print "num transitions:", len(forward_transitions.keys())
        
        # Store final versions of transitions and all_states
        self.transitions = forward_transitions
        self.allstates = [State(id) for id in forward_allstates]
        
        #for i, state in enumerate(self.allstates):
            #print i, state
        
        self.generateOutput()
        
        
    def initStuff(self, N):
        datamap = {}
        for i in range(2**N):
            name = tobase(2, i).zfill(N)
            iname = [int(x) for x in name]
            uname = unique_name(iname)
            #print "name:", name, type(name), uname
            datamap[name] = ''.join(['%d'%d for d in uname])
        #self.all_states = []
        self.datamap = datamap
        
    def addStates(self, state_id, state_index, target, newVal, reverse = False):
        transitions = self.transitions
        all_states = self.all_states

        [potentials, rawids] = self.generator(state_id, target, newVal)
        all_states = self.all_states
        
        for i_p, p_id in enumerate(potentials):
            newState = False
            try:
                p_index = all_states.index(p_id)
            except ValueError:
                p_index = -1
            #print "p_index?", p_index
            
            if p_index == -1:
               all_states.append(p_id)
               p_index = len(all_states) - 1
               newState = True
            elif (state_id,p_id) in transitions:
               transitions[(state_id, p_id)].m += 1
               continue  # NOTE: Only continue state generation with new, unique states!
            #create new transition (s,p)
            #print "state_id:", state_id
            #print "p_id:", p_id
            newTrans = Transition()
            transitions[(state_id, p_id)] = newTrans
            newTrans.iPre = state_index
            newTrans.iPost = p_index
            j = self.j
            newTrans.j = j
            self.j = j + 1 # Increment transition counter
            [RU, newVal] = getNewVal(state_id, rawids[i_p])
            newTrans.transType = getTransType(state_id, newVal, RU)
            


            if newState:
                if not reverse:
                    #recursively follow branch
                    self.addStates(p_id, p_index, '0', '1', reverse)
                    #self.addStates(p_id, p_index, '1', '2', reverse)
                else:
                    #self.addStates(p_id, p_index, '2', '1', reverse)
                    self.addStates(p_id, p_index, '1', '0', reverse)

    def generator(self, id, target, newVal):
        results = []
        #transTypes = []
        rawids = []
        for i, val in enumerate(id):
            if val == target:
                curid = id[:i] + newVal + id[i+1:]
                results.append(self.datamap[curid])
                rawids.append(curid)
            
        return results, rawids


    def generateOutput(self):
        
        #create sorted list by J's
        sortedKeys = self.transitions.keys()
        sortedKeys.sort(lambda a,b:cmp(self.transitions[a].j, self.transitions[b].j))
        
        alphas = [s.alpha for s in self.allstates]
        f = getfileptr('alphas', self.numRUs)
        print >>f, joinfields(map(str,alphas),'\n')
        f.close()
        
        m = [self.transitions[t].m for t in sortedKeys]
        f = getfileptr('m', self.numRUs)
        print >>f, joinfields(map(str,m),'\n')
        f.close()
        
        n = [self.transitions[t].n for t in sortedKeys]
        f = getfileptr('n', self.numRUs)
        print >>f, joinfields(map(str,n),'\n')
        f.close()
        
        iPre = [self.transitions[t].iPre+1 for t in sortedKeys]  # MATLAB indexing shift
        f = getfileptr('iPre', self.numRUs)
        print >>f, joinfields(map(str,iPre),'\n')         
        f.close()
        
        iPost = [self.transitions[t].iPost+1 for t in sortedKeys]  # MATLAB indexing shift
        f = getfileptr('iPost', self.numRUs)
        print >>f, joinfields(map(str,iPost),'\n')   
        f.close()
    
        transTypes = [self.transitions[t].transType for t in sortedKeys] 
        f = getfileptr('transTypes', self.numRUs)
        print >>f, joinfields(map(str, index2MATLAB(transTypes)),' ')  # MATLAB indexing shift
        f.close()
        
        mask_IDs = [val.id for val in self.allstates]
        f = getfileptr('mask_IDs', self.numRUs)
        print >>f, joinfields(mask_IDs, '\n')
        f.close()
        
        



def main():
    for i in range(3,16):
        startt = time.clock()
        MarkovChainGenerator(i)
        print "Elapsed time for numRUs = %d is: %.3f seconds"%(i, time.clock() - startt)


    

if __name__ == "__main__":
    main()
