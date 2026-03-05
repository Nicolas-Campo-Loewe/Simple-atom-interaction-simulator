# ATOMIC DYNAMICS SIMULATOR

**General explanation**

This project simulates atomic interactions based on simplified physical and chemical principles. The simulation models the motion of atoms, their interactions, and the formation of chemical bonds under different energetic and environmental conditions.

Atoms in the system possess dipole moments and kinetic energy, both of which influence their behavior and interactions. The kinetic energy of the system can be adjusted dynamically during runtime using the arrow keys (↑ / ↓), effectively simulating changes in temperature. Increasing the temperature raises particle motion and collision frequency, while lowering it stabilizes the system and facilitates bond formation.

When atoms approach each other, the simulation evaluates the possibility of forming chemical bonds in order to complete their valence electron shell. Depending on the electronegativity difference between atoms, different types of bonds may form:

- Covalent bonds
- Ionic bonds
- Metallic bonds

In addition, covalent interactions may result in single, double, or triple bonds, depending on the valence requirements of the atoms involved. The strength and stability of each bond depend on multiple factors such as electronegativity differences, bonding type, and local energetic conditions within the simulation.

The physical behavior of atoms and the chemical rules governing bonding are implemented separately in the project, to keep the system modular and easy to modify. Key parameters and rules can also be customized in their corresponding files: `physics.py`, `chemics.py`. This structure allows users to experiment with different modeled scenarios.


**How to Run**

1. Fork or clone this repository.
2. Open a terminal in the project directory.
3. Run the simulation with `main.py`


**Code demonstration**

https://www.youtube.com/watch?v=AUlwlbtnwDM
