""" quantum circuit gates """
import numpy as np
from scipy import linalg
def x_gate():
    """ x gate """
    return np.array([[0,1],[1,0]])
def y_gate():
    """ y gate """
    gate = [[0,complex(0,-1)],[complex(0,1),0]]
    return np.array(gate)
def z_gate():
    """ z gate """
    gate = [[1,0],[0,-1]]
    return np.array(gate)
def cnot_gate():
    """ c not gate """
    gate = [[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]]
    return np.array(gate)
def cnot_gate_reverse():
    """ c not gate reverse """
    gate = [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]
    return np.array(gate)
def h_gate():
    """ Hadamard """
    gate = [[1,1],[1,-1]]/np.sqrt(2)
    return np.array(gate)

def make_x(dim,subspace1):
    """ make x gate for n1-st subspace """
    xgate = x_gate()
    subspace2=dim-subspace1-1
    identity1=np.eye(2**subspace1)
    identity2=np.eye(2**subspace2)
    return linalg.kron(linalg.kron(identity2,xgate),identity1)
def make_h(dim,subspace1):
    """ Make hadamard gate in subspace """
    hgate = h_gate()
    subspace2=dim-subspace1-1
    identity1=np.eye(2**subspace1)
    identity2=np.eye(2**subspace2)
    return linalg.kron(linalg.kron(identity2,hgate),identity1)
def make_z(dim,subspace1):
    """ Make z gate in subspace """
    hgate = z_gate()
    subspace2=dim-subspace1-1
    identity1=np.eye(2**subspace1)
    identity2=np.eye(2**subspace2)
    return linalg.kron(linalg.kron(identity2,hgate),identity1)
def make_cnot(dim,subspace1,subspace2):
    """ Make cnot gate in subspace """
    if dim==5 and subspace1==1 and subspace2==4:
        matrix=np.zeros((2**dim)**2).reshape(2**dim,2**dim)
        rows=np.array([0,1,18,19,4,5,22,23,8,9,26,27,12,13,30,31,\
            16,17,2,3,20,21,6,7,24,25,10,11,28,29,14,15])
        assert len(rows)==2**dim
        for col in range(2**dim):
            matrix[rows[col],col]=1
        assert (np.matmul(matrix,matrix)==np.eye(2**dim)).any()
        return matrix
    if subspace2 not in (subspace1+1, subspace1-1):
        raise Exception('cnot not implemented for subspace1={subspace1} subspace2={subspace2}')
    if subspace1<subspace2:
        cnot = cnot_gate()
    else:
        cnot=cnot_gate_reverse()
    if subspace1<subspace2:
        dim1=subspace1
        dim2=dim-dim1-2
    else:
        dim1=subspace2
        dim2=dim-dim1-2
    assert dim1>=0
    assert dim2>=0
    identity1 = np.eye(2**dim1)
    identity2 = np.eye(2**dim2)
    prod = cnot
    if dim2!=0:
        prod = linalg.kron(identity2,cnot)
    if dim1 != 0:
        prod = linalg.kron(prod,identity1)
    matrix = prod
    assert (np.matmul(matrix,matrix)==np.eye(2**dim)).any()
    return matrix
def make_cnot2(dim,subspace1,subspace2):
    """ Make CNOT in subspace"""
    cnot = cnot_gate()
    assert subspace2==subspace1+1
    subspacea=subspace1
    subspaceb=dim-subspacea-2
    assert subspacea>=0
    assert subspaceb>=0
    identity1 = np.eye(2**subspacea)
    identity2 = np.eye(2**subspaceb)
    prod = linalg.kron(identity2,cnot)
    prod = linalg.kron(prod,identity1)
    matrix = prod
    assert (np.matmul(matrix,matrix)==np.eye(2**dim)).any()
    return matrix
def make_bell(dim,subspace1,subspace2, bit0, bit1):
    """ Make Bell state in subspace """
    xgate1 = make_x(dim,subspace1)
    xgate2 = make_x(dim,subspace2)
    gate = np.eye(2**dim)
    if bit1==1:
        gate = np.matmul(xgate1, gate)
    if bit0==1:
        gate = np.matmul(xgate2,gate)
    hgate = make_h(dim,subspace1)
    cnot = make_cnot(dim,subspace1,subspace2)
    gate = np.matmul(hgate,gate)
    gate = np.matmul(cnot,gate)
    return gate[:,0]

def bell_measure_matrix():
    """
    Phi+ (1,0,0,0) (00)
    Phi- (0,1,0,0) (10)
    Psi+ (0,0,1,0) (01)
    Psi- (0,0,0,1) (11)
    """
    return np.array(\
        [[1,0,0,1],\
        [1,0,0,-1],\
        [0,1,1,0],\
        [0,-1,1,0]])\
        /np.sqrt(2)
def bell_state(bit0,bit1):
    """
    Create Bell state
    """
    assert bit0 in [0,1]
    assert bit1 in [0,1]
    row = 2*bit0+bit1
    return bell_measure_matrix()[row,:]
def entangle_gate(dim,subspace1,subspace2, bit0, bit1):
    """
    Create gate for entangling states
    """
    xgate1 = make_x(dim,subspace1)
    xgate2 = make_x(dim,subspace2)
    gate = np.eye(2**dim)
    if bit1==1:
        gate = np.matmul(xgate1, gate)
    if bit0==1:
        gate = np.matmul(xgate2,gate)
    hgate = make_h(dim,subspace1)
    cnot = make_cnot2(dim,subspace1,subspace2)
    gate = np.matmul(hgate,gate)
    gate = np.matmul(cnot,gate)
    return gate
def phase_gate():
    """ Phase gate """
    return np.array([\
    [1,0],\
    [0,np.exp(complex(0,1)*np.pi/4)]])
def identity_gate(dim=1):
    """ Identity gate """
    return np.eye(2**dim)
