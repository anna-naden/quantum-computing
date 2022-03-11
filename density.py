"""
Stuff to deal with density matrices
"""
import numpy as np
from physical_simulation_utilities import dagger, state_sign_mask, get_log2
def partial_trace(rho,first,second):
    """
    Partial trace over first subsystem
    n=#qubits in first subsyste
    m=# qubits in second subsystem
    n+m=dim
    """
    first=2**first
    second=2**second
    return np.trace(rho.reshape(first,second,first,second), axis1=0, axis2=2)
def partial_trace_second_subsystem(rho,first,second):
    """
    Partial trace over second subsystem
    n=#qubits in first subsyste
    m=# qubits in second subsystem
    n+m=dim
    """
    first=2**first
    second=2**second
    return np.trace(rho.reshape(first,second,first,second), axis1=1, axis2=3)
def update(rho, unitary):
    """
    Apply time evolution operator to density matrix
    """
    return np.matmul(unitary,np.matmul(rho,dagger(unitary)))
def ispure(rho):
    """
    Check if density matrix represents a pure state
    """
    return np.trace(rho)==1
def valid_density(rho):
    """
    Check if density matrix is valid
    """
    assert np.imag(np.trace(rho))==0
    assert np.real(np.trace(rho)) <= 1
    assert np.real(np.trace(rho)) > 0.0
    assert (rho==dagger(rho)).all()
    assert rho.shape[0]==rho.shape[1]
    assert len(rho.shape)==2
    return \
        np.imag(np.trace(rho))==0 and \
        np.real(np.trace(rho)) <= 1 and \
        np.real(np.trace(rho)) > 0.0 and \
        (rho==dagger(rho)).all() and \
        rho.shape[0]==rho.shape[1] and \
        len(rho.shape)==2
def measurement_outcome(rho, measurement):
    """
    Compute probability of a measurement outcome and the corresponding projected state
    as a density matrix
    """
    valid_density(rho)
    prod = np.matmul(dagger(measurement),np.matmul(measurement, rho))
    prob = np.trace(prod)
    prod2 = np.matmul(measurement, np.matmul(rho, dagger(measurement)))
    if prob>0:
        rho = prod2/prob
        valid_density(rho)
    else:
        rho=None
    return rho, prob
def dump_rho(rho, title=None):
    """
    Dumpy density matrix
    """
    dim = get_log2(rho.shape[0])
    eigen = np.linalg.eig(rho)
    if title is not None:
        print(title)
    for i, prob in enumerate(eigen[0]):
        if prob>0:
            state = eigen[1][i]
            print(f'{prob}\n{state_sign_mask(dim,state)}')
def get_state(rho):
    """
    Extract state vector from density matrix
    """
    eigen = np.linalg.eig(rho)
    state = None
    num_states=0
    for i, prob in enumerate(eigen[0]):
        if np.round(prob,9)>0:
            state = eigen[1][:,i]
            ret_prob=prob
            num_states += 1
    assert num_states==1
    assert np.round(ret_prob,5)==1
    return state
