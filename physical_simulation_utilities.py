""" Utilities for symbolic entanglement code """
import numpy as np
from scipy import linalg

def get_log2(length):
    """ Helper for getting dimenstion of a state"""
    dim=0
    while 2**dim<length:
        dim += 1
    return dim
def same_state(state1,state2):
    """ Utility to check if two states are the ame"""
    cross = state1[1]*state2[0] - state1[0]*state2[1]
    return np.abs(cross)<1e-6
def inner_product(psi1,psi2):
    """ Inner product of two states"""
    return np.dot(np.conj(psi1),psi2)
def visualize(matrix, title=None):
    """ Utility to visualize a matrix """
    dim = get_log2(matrix.shape[0])
    print('01234567890123456789012345678901')
    if title is not None:
        print(title)
    for row in range(matrix.shape[0]):
        print(sign_mask(dim,matrix[row,:]))
    print()
def sign_mask(dim,state):
    """ Utility for visualizing a state """
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
        elif np.image(state[j])>0:
            str2[j]='+'
        else:
            str2[j]='-'
    return "".join(str1)
def state_sign_mask(dim,state, title=None):
    """ Utility to visualize a state """
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
        retval = title + '\n' + "".join(str1)+'\n' + \
            "".join(str2)+'\n01234567890123456789012345678901'+'\n'
    return retval
def display_state(state):
    """ Display the state """
    dim =get_log2(len(state))
    print(state_sign_mask(dim,state))
# def state_norm(state):
#     n=np.dot(np.conj(state),state)
#     return np.sqrt(n)
def n_particle_state(gate, num_particles):
    """ Create an n-particle state by applying a gate """
    return linalg.kron(np.eye(2**(num_particles-1)),gate)[:,0]
def one_particle_state(gate):
    """ Create a one-particle state by applying a gate """
    return gate[:,0]
def dagger(density):
    """ The 'dagger' of a state """
    return np.transpose(np.conj(density))
