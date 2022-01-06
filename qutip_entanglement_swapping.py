import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from physical_simulation_utilities import same_state
import sys
import qutip
from qutip import *
from qutip.measurement import *
from qutip.qip.operations import hadamard_transform, cnot
for initial_state in [basis(2,0), basis(2,1)]:
#   Make five-particle state
    q=tensor(initial_state, basis(2,0),basis(2,0), basis(2,0), basis(2,0))

    # Entangle data qubit (one) with link qubit (two) in the Psi minus state
    x1 = tensor(identity(2), sigmax(), identity(2), identity(2), identity(2))
    x2 = tensor(identity(2), identity(2), sigmax(), identity(2), identity(2))
    h = tensor(identity(2), hadamard_transform(), identity(2), identity(2), identity(2))
    not1 = tensor(identity(2), cnot(), identity(2), identity(2))
    psi_link = not1*h*x2*x1*q
    
    # Entangle qubits 3 and 4 in Psi minus Bell state
    x1 = tensor(identity(2), identity(2), identity(2), sigmax(), identity(2))
    x2 = tensor(identity(2), identity(2), identity(2), identity(2), sigmax())
    h = tensor(identity(2), identity(2), identity(2), hadamard_transform(), identity (2))
    not1 = tensor(identity(2), identity(2), identity(2), cnot())
    psi_link = not1 * h * x2 * x1 * psi_link
        
    # Bell basis for measurements
    phi_plus = bell_state('00')
    phi_minus = bell_state('01')
    psi_plus = bell_state('10')
    psi_minus = bell_state('11')
    for i1, psi1 in enumerate([phi_plus, phi_minus, psi_plus, psi_minus]):
        for i2, psi2 in enumerate([phi_plus, phi_minus, psi_plus, psi_minus]):
            
            # Measure qubits 2 and 3 in the Bell basis
            proj = psi1.proj()
            measure1 = tensor(identity(2), identity(2),proj, identity(2))
            psi_link1 = measure1*psi_link

            # Measure qubits 0 and 1 in the Bell basis
            proj = psi2.proj()
            measure1 = tensor(proj, identity(2), identity(2), identity(2))
            psi_link1 = measure1 * psi_link1
        
            # Receiving station extracts the state of the link bit
            
            # Get density matrix so that we can extract the last qubit
            dm = ket2dm(psi_link1)
            
            # Trace over the first four qubits
            dm = dm.ptrace([4])
            
            # Probabilities of the outcomes of measurement i1, i2. Only one outcome has nonzero probability.
            probs = dm.eigenstates()[0]
            states = dm.eigenstates()[1]
            for [p, state] in zip(probs, states):
                if p>0:
                    
                    # Depending on the measurement outcome, Bob applies a gate to his particle
                    if i1==0:
                        if i2 == 0:
                            gate = None
                        elif i2 == 1:
                            gate = sigmaz()
                        elif i2 == 2:
                            gate = sigmax()
                        else:
                            gate = sigmay()
                    elif i1==1:
                        if i2 == 0:
                            gate = sigmaz()
                        elif i2 == 1:
                            gate = None
                        elif i2 == 2:
                            gate = sigmay()
                        else:
                            gate = sigmax()

                    elif i1 == 2:
                        if i2 == 0:
                            gate = sigmax()
                        elif i2 == 1:
                            gate = sigmay()
                        elif i2 == 2:
                            gate = None
                        else:
                            gate = sigmaz()
                    elif i1 == 3:
                        if i2 == 0:
                            gate = sigmay()
                        elif i2 == 1:
                            gate = sigmax()
                        elif i2 == 2:
                            gate = sigmaz()
                        else:
                            gate = None
                    if gate is not None:
                        state = gate * state
                        
                    # Verify that receiving station successfully reconstructed the state of particle one
                    print(f'{i1} {i2}')
                    assert same_state(initial_state.full()[:,0], state.full()[:,0])


