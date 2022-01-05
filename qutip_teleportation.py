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
    q=tensor(initial_state,basis(2,0),basis(2,0))
    x1 = tensor(identity(2),sigmax(), identity(2))
    x2 = tensor(identity(2), identity(2),sigmax())
    h=tensor(identity(2),hadamard_transform(), identity (2))
    not1=tensor(identity (2),cnot())
    psi_23 = h*x2*x1*q
    psi_23 = not1*psi_23
    phi_plus = bell_state('00')
    phi_minus = bell_state('01')
    psi_plus = bell_state('10')
    psi_minus = bell_state('11')
    for i, psi in enumerate([phi_plus, phi_minus, psi_plus, psi_minus]):
        proj = psi.proj()
        measure = tensor(proj, identity(2))
        psi_123 = measure*psi_23
        dm = psi_123.proj()
        dm = dm.ptrace([2])
        probs = dm.eigenstates()[0]
        states = dm.eigenstates()[1]
        for [p, state] in zip(probs, states):
            if p>0:
                if i==0:
                    state = sigmay()*state
                elif i==1:
                    state = sigmax()*state
                elif i==2:
                    state = sigmaz()*state
                assert same_state(initial_state.full()[:,0], state.full()[:,0])


