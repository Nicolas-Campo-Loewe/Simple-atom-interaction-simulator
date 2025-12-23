import random
import math

from physics import (
    compute_nonbonded_force,
    compute_bond_force,
    compute_dipole_torque,
    integrate,
    wrap,
    WORLD_WIDTH,
    WORLD_HEIGHT,
)
from chemistry import (
    prebond_damp,
    try_form_bond,
    try_break_bond,
    reaction_addition,
    reaction_substitution,
    reaction_dissociation,
    remove_bond_between,
)

__all__ = ["Simulation", "Particle"]


class Particle:
    def __init__(self, x, y, vx, vy, p_type):
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.type = p_type
        self.bonds = []


class SpatialGrid:
    def __init__(self, cell_size):
        self.cell_size = cell_size
        self.grid = {}

    def clear(self):
        self.grid.clear()

    def insert(self, idx, x, y):
        cx = int(x // self.cell_size)
        cy = int(y // self.cell_size)
        self.grid.setdefault((cx, cy), []).append(idx)

    def neighbors(self, x, y):
        cx = int(x // self.cell_size)
        cy = int(y // self.cell_size)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                cell = (cx + dx, cy + dy)
                if cell in self.grid:
                    for idx in self.grid[cell]:
                        yield idx


class Simulation:
    def __init__(self, n_types=5, particles_per_type=35, temperature=0.8):
        self.temperature = temperature
        self.n_types = n_types
        self.particles = []
        self.bonds = []
        self.grid = SpatialGrid(cell_size=42.0)
        self._create_particles(particles_per_type)

    def _create_particles(self, n_per_type):
        for t in range(self.n_types):
            for _ in range(n_per_type):
                x = random.uniform(0, WORLD_WIDTH)
                y = random.uniform(0, WORLD_HEIGHT)
                vx = random.uniform(-2.2, 2.2)
                vy = random.uniform(-2.2, 2.2)
                self.particles.append(Particle(x, y, vx, vy, t))

    def step(self):
        parts = self.particles
        bonds = self.bonds
        n = len(parts)
        T = self.temperature

        self.grid.clear()
        for i, p in enumerate(parts):
            self.grid.insert(i, p.x, p.y)

        fx = [0.0] * n
        fy = [0.0] * n

        for i, p_i in enumerate(parts):
            for j in self.grid.neighbors(p_i.x, p_i.y):
                if j <= i:
                    continue
                if j in p_i.bonds:
                    continue

                p_j = parts[j]
                f = compute_nonbonded_force(p_i, p_j)
                fx[i] += f[0]
                fy[i] += f[1]
                fx[j] -= f[0]
                fy[j] -= f[1]

                self._process_reactions(i, j, T)

        bonds_to_break = []
        for b_idx, b in enumerate(bonds):
            b.age += 1  # rampa física

            i, j = b.i, b.j
            p_i, p_j = parts[i], parts[j]

            fbx, fby, broken_phys = compute_bond_force(p_i, p_j, b)
            fx[i] += fbx
            fy[i] += fby
            fx[j] -= fbx
            fy[j] -= fby

            tfx_i, tfy_i, tfx_j, tfy_j = compute_dipole_torque(p_i, p_j, b.dipole)
            fx[i] += tfx_i
            fy[i] += tfy_i
            fx[j] += tfx_j
            fy[j] += tfy_j

            if broken_phys or try_break_bond(b, parts, T):
                bonds_to_break.append(b_idx)

        for b_idx in reversed(bonds_to_break):
            b = bonds[b_idx]
            remove_bond_between(b.i, b.j, bonds, parts)

        # geometría covalente a baja T (torque angular ligero)
        if T < 1.0:
            self._apply_lowT_geometry_forces(fx, fy, T)

        forces = list(zip(fx, fy))
        integrate(parts, forces, T)

    def _process_reactions(self, i, j, T):
        parts = self.particles
        bonds = self.bonds

        dx = parts[j].x - parts[i].x
        dy = parts[j].y - parts[i].y
        dx, dy = wrap(dx, dy)
        r2 = dx * dx + dy * dy

        if r2 > 34.0 * 34.0:
            return

        prebond_damp(i, j, parts, bonds, T)

        if try_form_bond(i, j, parts, bonds, T):
            return

        reaction_addition(i, j, bonds, parts, T)
        reaction_substitution(i, j, bonds, parts, T)
        reaction_dissociation(i, bonds, parts, T)

    def _apply_lowT_geometry_forces(self, fx, fy, T):

        parts = self.particles

        target_angles = {
            2: math.pi,
            3: 2.0 * math.pi / 3.0,
            4: 109.5 * math.pi / 180.0,
        }

        strength = 0.015 * (1.0 - T)  # más fuerte cuanto más frío

        for i, p in enumerate(parts):
            neigh = p.bonds
            k = len(neigh)
            if k < 2 or k > 4:
                continue

            target = target_angles.get(k)
            if target is None:
                continue

            vecs = []
            for j in neigh:
                pj = parts[j]
                dx = pj.x - p.x
                dy = pj.y - p.y
                dx, dy = wrap(dx, dy)
                r2 = dx*dx + dy*dy
                if r2 < 1e-6:
                    continue
                r = math.sqrt(r2)
                vecs.append((dx/r, dy/r, j))

            if len(vecs) < 2:
                continue

            for a in range(len(vecs)):
                for b in range(a+1, len(vecs)):
                    ax, ay, ja = vecs[a]
                    bx, by, jb = vecs[b]
                    dot = max(-1.0, min(1.0, ax*bx + ay*by))
                    ang = math.acos(dot)
                    diff = ang - target

                    txa, tya = -ay, ax
                    txb, tyb = -by, bx

                    fx[i] += -strength * diff * (txa + txb)
                    fy[i] += -strength * diff * (tya + tyb)

    def get_particles(self):
        return self.particles

    def get_bonds(self):
        return self.bonds

    def increase_temperature(self, amount=0.1):
        self.temperature += amount

    def decrease_temperature(self, amount=0.1):
        self.temperature = max(0.05, self.temperature - amount)
