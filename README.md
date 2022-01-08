# Teleportation and a Quantum Repeater Using Entanglement Swapping
## qnet_teleportation.py
Using the qalgebra system to demonstration teleportation using computer algebra rather than numerical calculation.
## ent_swap_algebra.py
Use symbolic computation (computer algebra) to determine which gates to apply in qalgebra_entanglement_swapping.py
## qalgebra_entanglement_swapping.py
Like in a quantum repeater, using computer algebra rather than multiplying matrices of numbers
## qnet_requirements.lst
The pip requirements for usign qabebra and qnet
## qutip_teleportation.py
Simulates and demonstrates teleportation
## qutip_teleportation_nondeterministic
The qutip library function "measure_observable" is used to simulate a Bell measurement, and returns a nonzero probability, with the associted projected state, if the particle was measured to be in the pre-specified state. The outcome of this function call is nondeterministic.
## qutip_entanglement_swapping.py
Simulates and demonstrates entanglement swapping
## physical_simulation_utilities.py
Function same_state() returns True if both arguments represent the same state (up to a phase factor)
## physical-simulation-requirements.lst
I use numpy version 1.18 because the latest version is incompatible with qutip
