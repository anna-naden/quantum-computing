import qiskit
from qiskit import QuantumCircuit, execute
from qiskit.tools.visualization import plot_histogram
import matplotlib.pyplot as plt

counts=[]
for input_state in [0,1]:
    particle_1=''
    if input_state == 1:
        particle_1 = 'x q[0]; //Particle 1 spin down\n'
    cloner1 = """
    OPENQASM 2.0;
    include "qelib1.inc";

    qreg q[3]; //Three particles
    creg c[3]; //Three classical bits
    
    h q[1]; // Alice entangles particles two and three to get Psi_{23}^-
    cx q[1], q[2];
    x q[1];

    z q[1];
    
    """ + particle_1 + \
    """
    
    x q[0]; // Measurement: Bell operator basis
    x q[1];
    cx q[0],q[1];
    x q[1];
    x q[0];
    h q[0];
    
    
    measure q[0]->c[0]; //Alice measurement
    measure q[1]->c[1];

    if (c == 0) z q[2]; //Eqn 6 - Bob's transformation to get the teleported qubit
    if (c == 3) x q[2]; // |11>
    if (c == 2) y q[2]; // |10>
    
    measure q[2]->c[2]; // Measure the teleported qubit
    """

    qc= QuantumCircuit.from_qasm_str(cloner1)
    backend=qiskit.BasicAer.get_backend('qasm_simulator')
    result = execute(qc, backend=backend,shots=1000).result()
    counts.append(result.get_counts())
plot_histogram(counts[0])
plot_histogram(counts[1])
plt.show()
exit()