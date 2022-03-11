""" Teleporation and quantum repeater - without dm """
import numpy as np
from scipy import linalg
from physical_simulation_utilities import inner_product, get_log2
from gates import x_gate, y_gate, z_gate, make_bell,  entangle_gate
def get_state(state,bits01,bitsab,entangle_gate1,entangler):
    """ Reconstruct teleported/repeated state """
    dim = get_log2(len(state))
    reference_state = np.zeros(2**dim)
    reference_state[0]=1
    gates = [\
        [np.eye(2), z_gate(), x_gate(), y_gate()],\
         [z_gate(),np.eye(2), y_gate(), x_gate()],\
         [x_gate(), y_gate(), np.eye(2), z_gate()],\
          [y_gate(), x_gate(), z_gate(), np.eye(2)]]
    xgate = linalg.kron(x_gate(),np.eye(2**(dim-1)))
    psi = np.matmul(entangle_gate1,np.matmul(entangler,reference_state))
    psi1 = np.matmul(xgate,psi)
    amp0=2*inner_product(psi, state)
    amp1=2*inner_product(psi1,state)
    gate = gates[bits01][bitsab]
    return np.matmul(gate, np.array([amp0,amp1]))
def get_states(dim,state):
    """ Get states for a variety of measurement outcomes """
    states=[]
    for bit0 in [0,1]:
        for bit1 in [0,1]:
            entangler=entangle_gate(dim,2,3,bit0,bit1)
            for bita in [0,1]:
                for bitb in [0,1]:
                    states.append(\
                        get_state(state,2*bit0+bit1,\
                        2*bita+bitb,\
                        entangle_gate(dim,0,1,bita,bitb),\
                        entangler))
    return states
def quantum_repeater(state):
    """ Quantum repeater circuit """
    dim=get_log2(len(state))

    gate = entangle_gate(2,0,1,1,1)
    entangler1=linalg.kron(np.eye(4),linalg.kron(gate,np.eye(2)))
    entangler2=linalg.kron(gate,np.eye(8))
    state = np.matmul(entangler1,state)
    state = np.matmul(entangler2,state)
    return get_states(dim,state)
def teleport(state):
    """ Teleportation circuit """
    dim=get_log2(len(state))
    entangler = entangle_gate(2,0,1,1,1)
    entangler=linalg.kron(entangler,np.eye(2)) # entangle qubits 1 and 2, leave 0 unchanged
    xgate = linalg.kron(x_gate(),np.eye(2**(dim-1)))
    state = np.matmul(entangler,state)
    states=[]
    for bit0 in [0,1]:
        for bit1 in [0,1]:
            psi=make_bell(dim,0,1,bit0,bit1)
            psi1 = np.matmul(xgate, psi)
            amp0=2*inner_product(psi, state)
            amp1=2*inner_product(psi1,state)
            state2 = np.array([amp0,amp1])

            if [bit0,bit1]==[0,0]: #Phi+
                gate = y_gate()
            elif [bit0,bit1]==[0,1]: #Phi-
                gate = x_gate()
            elif [bit0,bit1]==[1,0]: #Psi+
                gate = z_gate()
            else:
                gate = np.eye(2)
            state2 = np.matmul(gate,state2)
            states.append(state2)
    return states
