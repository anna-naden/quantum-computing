import qiskit
from qiskit import QuantumCircuit, execute
from qiskit.providers.aer import QasmSimulator
from qiskit.tools.visualization import plot_histogram
from qiskit import IBMQ
from qiskit.tools.monitor import job_monitor

circuit = qiskit.QuantumCircuit(2,2)
circuit.cx(0,1)
circuit.cx(1,0)
circuit.cx(0,1)
circuit.measure(0,0)
circuit.draw()

counts=[]
for input_state in [0,1]:
    initialization=''
    if input_state == 1:
        initialization = 'x q[0]; //qubit to be cloned |1>\n'
    cloner1 = """
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[5];
    creg c[2];
    x q[2];  //Set reference qubit 2 to |1>
    x q[4];  //Set reference qubit 4 to |1>

    """ + initialization + \
    """
    measure q[0] -> c[0];
    if (c==1) CX q[1],q[2]; //swap reference qubits 1 and 2 if qubit to be cloned was |1>
    if (c==1) CX q[2],q[1];
    if (c==1) CX q[1],q[2];
    
    if (c==1) CX q[3],q[4]; //swap reference qubits 3 and 4 if qubit to be cloned was |1>
    if (c==1) CX q[4],q[3];
    if (c==1) CX q[3],q[4];
    measure q[1]->c[0];
    measure q[3]->c[1];
    """

    qc= QuantumCircuit.from_qasm_str(cloner1)
    result = execute(qc, backend=qiskit.BasicAer.get_backend('qasm_simulator'),shots=1024).result()
    counts.append(result.get_counts())
