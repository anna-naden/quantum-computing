import numpy as np
# from scipy import linalg
# import scipy
from new_physical_processes_density_matrix import teleport, quantum_repeater
from density import get_state, partial_trace
from physical_simulation_utilities import same_state, n_particle_state, \
    one_particle_state
from gates import x_gate, y_gate, z_gate, h_gate, phase_gate, identity_gate
# The gates used to initialize the particles
initializers = [identity_gate(), x_gate(), y_gate(), z_gate(), h_gate(), phase_gate()]

processes = [[teleport,3] , [quantum_repeater, 5]] # The physical processes
for [process, num_particles] in processes:
    for gate in initializers: # Simulate the physical process on all initial states
        initial_state = n_particle_state(gate, num_particles)
        initial_state = np.outer(np.conj(initial_state), initial_state)

        #Get and check all its possible outcomes
        for final_state in process(initial_state):
            n=num_particles-1 #number of qubits to trace
            m=1 # number of qubits to keep
            assert same_state(get_state(final_state), one_particle_state(gate))