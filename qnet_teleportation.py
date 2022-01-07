import qalgebra as qa
from sympy import sqrt
import sympy as sp
import numpy as np
hs = [qa.LocalSpace(i, dimension=2) for i in range(5)]
def hadamard_transform(hs):
    return (
        qa.LocalSigma(0, 0, hs=hs)  # LocalSigma(i,j) = |i⟩⟨j|
        + qa.LocalSigma(0, 1, hs=hs)
        + qa.LocalSigma(1, 0, hs=hs)
        - qa.LocalSigma(1, 1, hs=hs)
    ) / sqrt(2)
def ket(*indices, hs):
    local_spaces = hs.local_factors # for multiple indices, hs should be product space
    local_kets = [ls.basis_state(i) for (ls, i) in zip(local_spaces, indices)]
    return qa.TensorKet.create(*local_kets)
def ketbra(indices_ket, indices_bra, hs):
    return qa.KetBra(ket(*indices_ket, hs=hs), ket(*indices_bra, hs=hs))
def cnot(hs1, hs2):
    hs = hs1*hs2
    return (
        ketbra((0,0), (0,0), hs)
        + ketbra((0,1), (0,1), hs)
        + ketbra((1,0), (1,1), hs)
        + ketbra((1, 1), (1,0), hs)
    )
def sigmax(hs):
    return qa.LocalSigma(0, 1, hs=hs) + qa.LocalSigma(1, 0, hs=hs)
def sigmay(hs):
    return -sp.I*qa.LocalSigma(0,1,hs=hs) + sp.I*qa.LocalSigma(1,0,hs=hs)
def sigmaz(hs):
    return qa.LocalSigma(0,0, hs=hs)-qa.LocalSigma(1,1, hs=hs)
def bell_state(state, hs):
    bs = {
        '00': ket(0, 0, hs=hs) + ket(1, 1, hs=hs),
        '01': ket(0, 0, hs=hs) - ket(1, 1, hs=hs),
        '10': ket(0, 1, hs=hs) + ket(1, 0, hs=hs),
        '11': ket(0, 1, hs=hs) - ket(1, 0, hs=hs),
    }
    return bs[state] / sqrt(2)
def proj(state):
    return (state * state.dag()).expand()

initial_state = (hs[0].basis_state(0)+hs[0].basis_state(1))/sp.sqrt(2) # fist one in the loop
expected_final = (hs[2].basis_state(0)+hs[2].basis_state(1))/sp.sqrt(2)
dm_expected_final = proj(expected_final)

# initial_state = hs[0].basis_state(1)  # fist one in the loop
q = initial_state * hs[1].basis_state(0) * hs[2].basis_state(0)
x1 = sigmax(hs[1])

x2 = sigmax(hs[2])

h = hadamard_transform(hs[1])

not1 = cnot(hs[1], hs[2])

psi_23 = (h * x2 * x1 * q).expand()
psi_23 = (not1 * psi_23).expand()
phi_plus = bell_state('00', hs[0]*hs[1])
phi_minus = bell_state('01', hs[0]*hs[1])
psi_plus= bell_state('10', hs[0]*hs[1])
psi_minus = bell_state('11', hs[0]*hs[1])
def do_loop():
    # I'm just putting it in a function because I can't split a loop over
    # multiple cells
    result = []
    for i, psi in enumerate([phi_plus, phi_minus, psi_plus, psi_minus]):
        measure1 = proj(psi)
        psi_123 = measure1 * psi_23
        dm = proj(psi_123)
        dm = qa.OperatorTrace(dm, over_space=hs[1])
        dm = qa.OperatorTrace(dm, over_space=hs[0]).doit()

        if i==0:
            dm = sigmay(hs[2]).dag()*dm*sigmay(hs[2])
        elif i==1:
            dm = sigmax(hs[2]).dag()*dm*sigmax(hs[2])
        elif i==2:
            dm = sigmaz(hs[2]).dag()*dm*sigmaz(hs[2])
        result.append(dm)
        assert (dm_expected_final/4).expand_in_basis()==dm.expand_in_basis()
    return result

result=do_loop()
# print(result)