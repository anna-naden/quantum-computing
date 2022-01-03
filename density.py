import numpy as np
from physical_simulation_utilities import dagger, state_sign_mask, get_log2
def partial_trace(rho,n,m):
    """
    n=#qubits to trace over
    m=#untraced qubits
    n+m=dim
    
    """
    n=2**n
    m=2**m
    return np.trace(rho.reshape(n,m,n,m), axis1=0, axis2=2)
def update(rho, u):
    return np.matmul(u,np.matmul(rho,dagger(u)))
def ispure(rho):
    return np.trace(rho)==1
def valid_density(rho):
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
    dim = get_log2(rho.shape[0])
    eigen = np.linalg.eig(rho)
    if title is not None:
        print(title)
    for i, p in enumerate(eigen[0]):
        if p>0:
            state = eigen[1][i]
            print(f'{p}\n{state_sign_mask(dim,state)}')
def get_state(rho):
    dim = get_log2(rho.shape[0])
    eigen = np.linalg.eig(rho)
    state = None
    n=0
    for i, p in enumerate(eigen[0]):
        if np.round(p,9)>0:
            state = eigen[1][:,i]
            p2=p
            n += 1
    assert n==1
    assert np.round(p2,5)==1
    return state

