# https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.70.1895

from qiskit import QuantumCircuit, execute
from qiskit.providers.aer import QasmSimulator
from qiskit.tools.visualization import plot_histogram
import qiskit
import matplotlib.pyplot as plt

counts=[]
results=[]
for input_state in [0,1,2]:
    data =''
    measure=''
    if input_state == 1:
        data = 'x input; //Input qubit |1>\n'
        measure = ''        
    elif input_state == 2:
        data = 'h input; //Hadamard gate: input qubit = |0>+|1>\n'
        measure = 'h target; //Output qubit'
    transporter = \
    """
    OPENQASM 2.0;
    include "qelib1.inc";
    
    creg c[2]; //Classical bits
    creg c_output[1]; //Measurement of the teleported qubit
    qreg input[1]; //Qubit/particle state to be teleported
    qreg link[1];
    qreg target[1]; // Reconstructed state

    gate entangle() q1,q2 //Bell basis
    {
        x q1;
        x q2;
        h q1;
        cx q1,q2;
    }

    gate bell_measure() q1,q2
    {
        x q2;
        cx q1,q2;
        x q2;
        x q1;
        h q1;
        z q1;
    }
    

    entangle() link[0], target[0];

    
    """ + data + \
    """
    bell_measure() input[0], link[0];
    measure input[0]->c[0];
    measure link[0]->c[1];
    
    if (c==0) y target[0];
    if (c==1) x target[0];
    if (c==2) z target[0];
    """ + measure + \
    """
    measure target[0]->c_output[0];
    """
    
    qc= QuantumCircuit.from_qasm_str(transporter)
    backend=qiskit.BasicAer.get_backend('qasm_simulator')
    result = execute(qc, backend=backend,shots=1000).result()
    counts.append(result.get_counts())
for c in counts:
    plot_histogram(c)
plt.show()