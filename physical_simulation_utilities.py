import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
def plot_counts(counts):
    heights=[]
    labels=[]
    for c in counts:
        if c[0]>0:
            heights.append(c[0])
            labels.append(c[1])
    plt.bar(labels,heights)
    plt.show()
def get_log2(n):
    i=0
    while 2**i<n:
        i += 1
    return i
def element_to_qubits(dim,element):
    bin = format(element,'b')
    bin = str.rjust(bin,dim,'0')
    bin = bin[::-1]
    lbin = list(bin)
    qubits = np.zeros(dim, dtype=int)
    for q, b in enumerate(lbin):
        qubits[q]=b
    return qubits
def state_mask(dim,state):
    str1 = list('                                ')
    for j in range(2**dim):
        if np.around(np.abs(state[j]),3)>.001:
            str1[j]='1'
        else:
            str1[j]='0'
    return "".join(str1)+'\n01234567890123456789012345678901'
def same_state(state1,state2):
    for index, amp in enumerate(state1):
        if np.abs(amp)>0:
            break
    if state2[index]==0:
        return False
    phase_factor = state2[index]/state1[index]
    for index, amp in enumerate(state1):
        if state2[index] != phase_factor*state1[index]:
            return False
    return True
def inner_product(psi1,psi2):
    return np.dot(np.conj(psi1),psi2)
def visualize(matrix, title=None):
    dim = get_log2(matrix.shape[0])
    print('01234567890123456789012345678901')
    if title is not None:
        print(title)
    for row in range(matrix.shape[0]):
        print(sign_mask(dim,matrix[row,:]))
    print()
def state_sign_mask(dim,state, title=None):
    str1 = list('                                ')
    for j in range(2**dim):
        if np.real(state[j])==0:
            str1[j]='0'
        elif np.real(state[j])>0:
            str1[j]='+'
        else:
            str1[j]='-'
    str2 = list('                                ')
    for j in range(2**dim):
        if np.imag(state[j])==0:
            str2[j]='0'
        elif np.imag(np.real(state[j]))>0:
            str2[j]='+'
        else:
            str2[j]='-'
    if title is None:
        retval = "".join(str1)+'\n' + "".join(str2)+'\n01234567890123456789012345678901'+'\n'
    else:
        retval = title + '\n' + "".join(str1)+'\n' + "".join(str2)+'\n01234567890123456789012345678901'+'\n'
    return retval
def sign_mask(dim,state, title=None):
    str1 = list('                                ')
    for j in range(2**dim):
        if state[j]==0:
            str1[j]='0'
        elif state[j]>0:
            str1[j]='+'
        else:
            str1[j]='-'
    return "".join(str1)
def display_state(state):
    dim =get_log2(len(state))
    print(state_sign_mask(dim,state))
def initial_state(dim,qubit):
    state = np.zeros(2**dim)
    if qubit==0:
        state[0]=1
    else:
        state[1]=1
    return state
def state_norm(state):
    n=np.dot(np.conj(state),state)
    return np.sqrt(n)
def measure_qubit(qubit, state):
    dim = get_log2(len(state))
    m = measurement_operator0(dim,qubit)
    state0 = np.matmul(m,state)
    if state_norm(state0)==0:
        state0=None
        prob0=0
    else:
        prob0=state_norm(state0)**2
        state0 /= state_norm(state0)
    
    m = measurement_operator1(dim,qubit)
    state1 = np.matmul(m,state)
    if state_norm(state1)==0:
        state1=None
        prob1=0
    else:
        prob1=state_norm(state1)**2
        state1 /= state_norm(state1)
    return [[prob0,state0],[prob1,state1]]
def get_counts(state):
    dim = get_log2(len(state))
    counts=[]
    for i in range(len(state)):
        state1 = np.zeros(len(state))
        state1[i]=1
        proj = np.outer(state1,state1)
        state2 = np.matmul(proj,state)
        n = state_norm(state2)
        p = n**2
        if p==0:
            state2=None
        else:
            state2 /n
        string = bin(i)[2:]
        string = string.rjust(dim,'0')
        # string = string[::-1]
        counts.append([p,string])
    return counts
def count_to_bits(count):
    dim = len(count)
    bits = np.zeros(dim, dtype=int)
    for i in range(dim):
        bits[i]=0
        if count[dim-i-1]=='1':
            bits[i]=1
    return bits
def n_particle_state(gate, num_particles):
    return linalg.kron(np.eye(2**(num_particles-1)),gate)[:,0]
def one_particle_state(gate):
    return gate[:,0]


