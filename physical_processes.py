import numpy as np
from scipy import linalg
from physical_simulation_utilities import inner_product, get_log2, state_sign_mask, visualize
from gates import x_gate, y_gate, z_gate, make_bell,  entangle_gate
def quantum_repeater(state):
    dim=get_log2(len(state))

    e = entangle_gate(2,0,1,1,1)
    g=linalg.kron(np.eye(4),linalg.kron(e,np.eye(2)))
    visualize(g,'first entangle')
    g2=linalg.kron(e,np.eye(8))
    print(state_sign_mask(dim,state, 'initial'))
    state = np.matmul(g,state)
    print(state_sign_mask(dim, state, 'after first entangle'))
    state = np.matmul(g2,state)
    xgate = linalg.kron(x_gate(),np.eye(2**(dim-1)))
    # print(state_sign_mask(dim,psi))
    reference_state = np.zeros(2**dim)
    reference_state[0]=1
    states=[]
    for x in [0,1]:
        for y in [0,1]:
            g1=entangle_gate(dim,2,3,x,y)
            for xx in [0,1]:
                for yy in [0,1]:
                    g2=entangle_gate(dim,0,1,xx,yy)
                    state2 = np.matmul(g1,reference_state)
                    psi = np.matmul(g2,state2)
                    # print(state_sign_mask(dim,psi))
                    psi1 = np.matmul(xgate,psi)
                    a=2*inner_product(psi, state)
                    b=2*inner_product(psi1,state)
                    state2 = np.array([a,b])         
                    if [x,y]==[0,0]: #Phi+apbp
                        if [xx,yy]==[0,0]: # Phi+
                            gate = np.eye(2)
                        elif [xx,yy]==[0,1]: #Phi-
                            gate = z_gate()
                        elif [xx,yy]==[1,0]: #Psi+
                            gate = x_gate()
                        else: #Psi-
                            gate = y_gate()
                    elif [x,y]==[0,1]: #Phi-apbp
                        if [xx,yy]==[0,0]: #Phi+
                            gate = z_gate()
                        elif [xx,yy]==[0,1]: #Phi-
                            gate = np.eye(2)
                        elif [xx,yy]==[1,0]: #Psi+
                            gate = y_gate()
                        else: #Psi-
                            gate = x_gate()
                    elif [x,y]==[1,0]: #Psi+apbp
                        if [xx,yy]==[0,0]: #Phi+
                            gate = x_gate()
                        elif [xx,yy]==[0,1]: #Phi-
                            gate = y_gate()
                        elif [xx,yy]==[1,0]: #Psi+
                            gate = np.eye(2)
                        else: #Psi-
                            gate=z_gate()
                    else: #Psi - apbp
                        if [xx,yy]==[0,0]: #Phi+
                            gate=y_gate()
                        elif [xx,yy]==[0,1]: #Phi-
                            gate=x_gate()
                        elif [xx,yy]==[1,0]: #Psi+
                            gate=z_gate()
                        else: #Psi-
                            gate=np.eye(2)
                    state2 = np.matmul(gate, state2)
                    states.append(state2)
    return states

def teleport(state):
    dim=get_log2(len(state))
    e = entangle_gate(2,0,1,1,1)
    g1=linalg.kron(e,np.eye(2)) # entangle qubits 1 and 2, leave 0 unchanged
    xgate = linalg.kron(x_gate(),np.eye(2**(dim-1)))
    state = np.matmul(g1,state)
    states=[]
    for x in [0,1]:
        for y in [0,1]:
            psi=make_bell(dim,0,1,x,y)
            psi1 = np.matmul(xgate, psi)
            a=2*inner_product(psi, state)
            b=2*inner_product(psi1,state)
            state2 = np.array([a,b])         

            if [x,y]==[0,0]: #Phi+
                gate = y_gate()
            elif [x,y]==[0,1]: #Phi-
                gate = x_gate()
            elif [x,y]==[1,0]: #Psi+
                gate = z_gate()
            else:
                gate = np.eye(2)
            state2 = np.matmul(gate,state2)
            states.append(state2)
    return states