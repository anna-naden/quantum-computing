import numpy as np
from scipy import linalg
def x_gate():
    return np.array([[0,1],[1,0]])
def y_gate():
    gate = [[0,complex(0,-1)],[complex(0,1),0]]
    return np.array(gate)
def z_gate():
    gate = [[1,0],[0,-1]]
    return np.array(gate)
def cnot_gate():
    gate = [[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]]
    return np.array(gate)
def cnot_gate_reverse():
    gate = [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]
    return np.array(gate)
def h_gate():
    # Hadamard
    gate = [[1,1],[1,-1]]/np.sqrt(2)
    return np.array(gate)

def make_x(dim,q):
    x = x_gate()
    n1=q
    n2=dim-n1-1
    i1=np.eye(2**n1)
    i2=np.eye(2**n2)
    return linalg.kron(linalg.kron(i2,x),i1)
def make_h(dim,q):
    h = h_gate()
    n1=q
    n2=dim-n1-1
    i1=np.eye(2**n1)
    i2=np.eye(2**n2)
    return linalg.kron(linalg.kron(i2,h),i1)
def make_z(dim,q):
    h = z_gate()
    n1=q
    n2=dim-n1-1
    i1=np.eye(2**n1)
    i2=np.eye(2**n2)
    return linalg.kron(linalg.kron(i2,h),i1)
def make_cnot(dim,q1,q2):
    if dim==5 and q1==1 and q2==4:
        matrix=np.zeros((2**dim)**2).reshape(2**dim,2**dim)
        rows=np.array([0,1,18,19,4,5,22,23,8,9,26,27,12,13,30,31,\
            16,17,2,3,20,21,6,7,24,25,10,11,28,29,14,15])
        assert(len(rows)==2**dim)
        for col in range(2**dim):
            matrix[rows[col],col]=1
        assert (np.matmul(matrix,matrix)==np.eye(2**dim)).any()
        return matrix
    elif q2 != q1+1 and q2 != q1-1:
        raise Exception('cnot not implemented for q1={q1} q2={q2}')
    if q1<q2:
        cnot = cnot_gate()
    else:
        cnot=cnot_gate_reverse()
    if q1<q2:
        n1=q1
        n2=dim-n1-2
    else:
        n1=q2
        n2=dim-n1-2
    assert n1>=0
    assert n2>=0
    i1 = np.eye(2**n1)
    i2 = np.eye(2**n2)
    prod = cnot
    if n2!=0:
        prod = linalg.kron(i2,cnot)
    if n1 != 0:
        prod = linalg.kron(prod,i1)
    matrix = prod
    assert (np.matmul(matrix,matrix)==np.eye(2**dim)).any()
    return matrix
def make_cnot2(dim,q1,q2):
    cnot = cnot_gate()
    assert q2==q1+1
    n1=q1
    n2=dim-n1-2
    assert n1>=0
    assert n2>=0
    i1 = np.eye(2**n1)
    i2 = np.eye(2**n2)
    prod = linalg.kron(i2,cnot)
    prod = linalg.kron(prod,i1)
    matrix = prod
    assert (np.matmul(matrix,matrix)==np.eye(2**dim)).any()
    return matrix
def make_bell(dim,q1,q2, x, y):
    x1 = make_x(dim,q1)
    x2 = make_x(dim,q2)
    gate = np.eye(2**dim)
    if y==1:
        gate = np.matmul(x1, gate)
    if x==1:
        gate = np.matmul(x2,gate)
    h = make_h(dim,q1)
    cnot = make_cnot(dim,q1,q2)
    gate = np.matmul(h,gate)
    gate = np.matmul(cnot,gate)
    return gate[:,0]

def bell_measure_matrix():
    # Phi+ (1,0,0,0) (00)
    # Phi- (0,1,0,0) (10)
    # Psi+ (0,0,1,0) (01)
    # Psi- (0,0,0,1) (11)
    return np.array(\
        [[1,0,0,1],\
        [1,0,0,-1],\
        [0,1,1,0],\
        [0,-1,1,0]])\
        /np.sqrt(2)
def bell_state(x,y):
    assert x in [0,1]
    assert y in [0,1]
    row = 2*x+y
    return bell_measure_matrix()[row,:]
def entangle_gate(dim,q1,q2, x, y):
    x1 = make_x(dim,q1)
    x2 = make_x(dim,q2)
    gate = np.eye(2**dim)
    if y==1:
        gate = np.matmul(x1, gate)
    if x==1:
        gate = np.matmul(x2,gate)
    h = make_h(dim,q1)
    cnot = make_cnot2(dim,q1,q2)
    gate = np.matmul(h,gate)
    gate = np.matmul(cnot,gate)
    return gate
def phase_gate():
    return np.array([\
    [1,0],\
    [0,np.exp(complex(0,1)*np.pi/4)]])
def identity_gate(dim=1):
    return np.eye(2**dim)


