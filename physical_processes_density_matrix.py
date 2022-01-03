import numpy as np
from scipy import linalg
from physical_simulation_utilities import inner_product, get_log2, dagger, visualize, state_sign_mask
from density import update, measurement_outcome, get_state, partial_trace_second_subsystem
from gates import x_gate, y_gate, z_gate, make_bell,  entangle_gate, make_x, identity_gate, entangle_gate
def teleport(rho):
    dim = get_log2(rho.shape[0])

    # Get projection operators for all possible measurment outcomes
    measurements = np.empty([2,2, 2**dim,2**dim])
    for x in [0,1]:
        for y in [0,1]:
            psi = identity_gate(dim)[:,0] #Initial state - all zeros
            gate = entangle_gate(dim,0,1,x,y)
            psi0 = np.matmul(gate,psi)
            # Toggle the last qubit
            psi1 = np.matmul(make_x(dim,dim-1), psi0)
            # Measurement projection operator
            m= np.outer(np.conj(psi0),psi0) + np.outer(np.conj(psi1), psi1)
            measurements[x,y] =m

    e = entangle_gate(2,0,1,1,1)
    g1=linalg.kron(e,np.eye(2)) # entangle qubits 1 and 2, leave 0 unchanged
    rho=update(rho, g1)
    rhos=[]
    
    # Cycle through all possible outcomes
    for x in [0,1]:
        for y in [0,1]:
            rho_meas, p = measurement_outcome(rho,measurements[x,y])
            rho_meas = partial_trace_second_subsystem(rho_meas, 1, 2)
            # Probability of getting this outcome
            if p>0:
                
                # Pick a gate (unitary operation) to reconstruct the teleported qubit
                if [x,y]==[0,0]: #Phi+
                    gate = y_gate()
                elif [x,y]==[0,1]: #Phi-
                    gate = x_gate()
                elif [x,y]==[1,0]: #Psi+
                    gate = z_gate()
                else:
                    gate = np.eye(2)
                rho_meas1 = update(rho_meas, gate)
                rhos.append(rho_meas1)
    return rhos

def quantum_repeater(rho):
    
    dim=get_log2(rho.shape[0])
    
    # Get projection operators for measurements of all qubits
    measurements = np.empty([2,2,2,2,2**dim,2**dim])
    for x in [0,1]:
        for y in [0,1]:
            for xx in [0,1]:
                for yy in [0,1]:
                    
                    # Bell basis qubits 1 and 2
                    psi_i = identity_gate(dim)[:,0]
                    gate = entangle_gate(dim,1,2,x,y)
                    psi0 = np.matmul(gate, psi_i)
                    
                    # Bel basis qubits 3 and 4
                    gate = entangle_gate(dim, 3, 4, xx, yy)
                    psi1 = np.matmul(gate, psi_i)
                    
                    # Target qubit 4
                    psi1 = np.matmul(make_x(dim, dim-1), psi1)
                    m = np.outer(np.conj(psi0), psi0) + np.outer(np.conj(psi1), psi1)
                    measurements[x,y,xx,yy]=m
                        
    # Bell Psi-
    e = entangle_gate(2,0,1,1,1)
    g=linalg.kron(np.eye(4),linalg.kron(e,np.eye(2))) #Entangle 1 and 2
    g2=linalg.kron(e,np.eye(8)) #Entangle 3 and 4
    
    rho = update(rho,g)
    rho = update(rho,g2)
    
    # Cycle through all possible measurement outcomes
    rhos=[]
    for x in [0,1]:
        for y in [0,1]:
            for xx in [0,1]:
                for yy in [0,1]:
                    rho_meas, p = measurement_outcome(rho, measurements[x,y,xx,yy])

                    # probability of getting this measurement outcome
                    if p> 0:
                        rho_meas1 = partial_trace_second_subsystem(rho_meas, 1,dim-1)

                        # Based on measurement outcome, pick a gate (unitary operation) to transform the
                        # final qubit into a replica of the transmitted qubit
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
                        rho_meas1 = update(rho_meas1,gate)
                        rhos.append(rho_meas1)
    return rhos

