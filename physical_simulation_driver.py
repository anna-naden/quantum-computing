from physical_processes import teleport, quantum_repeater
from physical_simulation_utilities import same_state, n_particle_state, one_particle_state
from gates import x_gate, y_gate, z_gate, h_gate, phase_gate, identity_gate

# The gates used to initialize the particles
initializers = [identity_gate(), x_gate(), y_gate(), z_gate(), h_gate(), phase_gate()]

processes = [[teleport,3] , [quantum_repeater, 5]] # The physical processes
for [process, num_particles] in processes:
    for gate in initializers: # Simulate the physical process on all initial states

        #Get all its possible outcomes
        for final_state in process(n_particle_state(gate, num_particles)): 
                assert same_state(one_particle_state(gate), final_state)