from qiskit import QuantumCircuit, execute
from qiskit.providers.aer import QasmSimulator
from qiskit.tools.visualization import plot_histogram
import qiskit
import matplotlib.pyplot as plt

counts=[]
results=[]
for input_state in [0,1,2]:
    input_data =''
    measure=''
    if input_state == 1:
        input_data = 'x input; //Input qubit |1>\n'
        measure = ''        
    elif input_state == 2:
        input_data = 'h input; //Hadamard gate: input qubit = |0>+|1>\n'
        measure = 'h output; //Output qubit'
    repeater = \
    """
    OPENQASM 2.0;
    include "qelib1.inc";
    
    creg c[4]; // Results of two Bell measurements
    creg c_output[1]; 
    qreg input[1]; // The qubit to be transmitted
    qreg link[1]; // Link to repeater
    qreg repeater1[1];
    qreg repeater2[1];
    qreg output[1];// The reconstructed qubit

    gate entangle() q1,q2
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
    

    entangle() link,repeater1;
    entangle() repeater2,output;

    
    """ + input_data + \
    """
    bell_measure() repeater1, repeater2;
    measure repeater1->c[2];
    measure repeater2->c[3];
    
    bell_measure() input,link;
    measure input->c[0];
    measure link->c[1];
    
    // Reconstruct transmitted qubit based on measurement outcome
    if (c==1) z output;
    if (c==2) x output;
    if (c==3) y output;
    if (c==4) z output;
    if (c==6) y output;
    if (c==7) x output;
    if (c==8) x output;
    if (c==9) y output;
    if (c==11) z output;
    if (c==12) y output;
    if (c==13) x output;
    if (c==14) z output;
    """ + \
    measure + \
    """
    measure output->c_output;
    """
    
    qc= QuantumCircuit.from_qasm_str(repeater)
    backend=qiskit.BasicAer.get_backend('qasm_simulator')
    result = execute(qc, backend=backend,shots=6553).result()
    counts.append(result.get_counts())
    results.append(result)
for c in counts:
    plot_histogram(c)
plt.show()