import math
import random

DEFAULT_TEMPERATURE = 1.0
TIME_STEP = 0.2
MAX_SPEED = 9.5

THERMOSTAT_STRENGTH = 0.02
LANGEVIN_GAMMA = 0.04       

NONBONDED_REPULSION = 140.0
NONBONDED_SOFT_REPULSION = 6.0
NONBONDED_RANGE = 42.0

BOND_SPRING_CONSTANT = 0.32
BOND_DAMPING = 0.03
BOND_BREAK_THRESHOLD = 1.9     # factor sobre longitud de equilibrio
BOND_FORCE_CAP = 6.5           # evita "patadas" numéricas
BOND_SOFTSTART_STEPS = 12      # rampa de entrada de fuerza

DIPOLE_STRENGTH = 0.22
DIPOLE_VISUAL_SCALE = 18.0

WORLD_WIDTH = 800
WORLD_HEIGHT = 800


def wrap(dx, dy):
    """Condiciones periódicas (toro)."""
    if dx > WORLD_WIDTH * 0.5:
        dx -= WORLD_WIDTH
    elif dx < -WORLD_WIDTH * 0.5:
        dx += WORLD_WIDTH

    if dy > WORLD_HEIGHT * 0.5:
        dy -= WORLD_HEIGHT
    elif dy < -WORLD_HEIGHT * 0.5:
        dy += WORLD_HEIGHT

    return dx, dy


def limit_speed(vx, vy, max_speed=MAX_SPEED):
    v2 = vx * vx + vy * vy
    if v2 > max_speed * max_speed:
        s = max_speed / math.sqrt(v2)
        return vx * s, vy * s
    return vx, vy

def compute_nonbonded_force(p_i, p_j):
    dx = p_j.x - p_i.x
    dy = p_j.y - p_i.y
    dx, dy = wrap(dx, dy)

    r2 = dx * dx + dy * dy
    if r2 < 1e-9:
        return (0.0, 0.0)

    r = math.sqrt(r2)
    if r > NONBONDED_RANGE:
        return (0.0, 0.0)

    nx = dx / r
    ny = dy / r

    if r < 7.5:
        F = NONBONDED_REPULSION / (r2 + 0.08)
        return (-nx * F, -ny * F)

    rep = NONBONDED_SOFT_REPULSION / r
    F = -rep
    return (nx * F, ny * F)

def compute_bond_force(p_i, p_j, bond):
    dx = p_j.x - p_i.x
    dy = p_j.y - p_i.y
    dx, dy = wrap(dx, dy)

    r2 = dx * dx + dy * dy
    if r2 < 1e-9:
        return (0.0, 0.0, False)

    r = math.sqrt(r2)
    nx = dx / r
    ny = dy / r

    soften = 1.0
    if getattr(bond, "age", 0) < BOND_SOFTSTART_STEPS:
        soften = (bond.age + 1) / BOND_SOFTSTART_STEPS

    k_eff = BOND_SPRING_CONSTANT * bond.order * soften

    cooper = 1.0 + 0.18 * min(len(p_i.bonds), len(p_j.bonds))
    k_eff *= cooper

    stretch = r - bond.length
    F_spring = k_eff * stretch

    dvx = p_j.vx - p_i.vx
    dvy = p_j.vy - p_i.vy
    F_damp = BOND_DAMPING * (dvx * nx + dvy * ny) * bond.order

    F_total = F_spring + F_damp

    if F_total > BOND_FORCE_CAP:
        F_total = BOND_FORCE_CAP
    elif F_total < -BOND_FORCE_CAP:
        F_total = -BOND_FORCE_CAP

    fx = nx * F_total
    fy = ny * F_total

    broken = (r > bond.length * BOND_BREAK_THRESHOLD)
    return (fx, fy, broken)

def compute_dipole_torque(p_i, p_j, dipole_vector):
    dx = p_j.x - p_i.x
    dy = p_j.y - p_i.y
    dx, dy = wrap(dx, dy)

    r2 = dx * dx + dy * dy
    if r2 < 1e-9:
        return (0.0, 0.0, 0.0, 0.0)

    r = math.sqrt(r2)
    nx = dx / r
    ny = dy / r

    tx = -ny
    ty = nx

    desired = dipole_vector[0] * tx + dipole_vector[1] * ty
    torque = desired * DIPOLE_STRENGTH

    fx_i = tx * torque
    fy_i = ty * torque
    fx_j = -fx_i
    fy_j = -fy_i

    return (fx_i, fy_i, fx_j, fy_j)

def apply_langevin(p, temperature, bond_count, avg_bond_order):

    gamma = LANGEVIN_GAMMA

    sigma = math.sqrt(2.0 * gamma * temperature) * bonding_factor

    p.vx += -gamma * p.vx + random.gauss(0.0, sigma)
    p.vy += -gamma * p.vy + random.gauss(0.0, sigma)


def rescale_velocities(particles, temperature):
    if not particles:
        return

    Ek = 0.0
    for p in particles:
        Ek += p.vx * p.vx + p.vy * p.vy

    T_current = Ek / len(particles)
    if T_current < 1e-9:
        return

    scale_target = math.sqrt(temperature / T_current)
    scale = 1.0 + (scale_target - 1.0) * THERMOSTAT_STRENGTH

    for p in particles:
        p.vx *= scale
        p.vy *= scale


def integrate(particles, forces, temperature):
    dt = TIME_STEP

    for i, p in enumerate(particles):
        fx, fy = forces[i]
        p.vx += fx * dt
        p.vy += fy * dt

        p.vx, p.vy = limit_speed(p.vx, p.vy)

    for p in particles:
        bond_count = len(p.bonds)
        if bond_count == 0:
            avg_order = 0.0
        else:
            avg_order = 1.5
        apply_langevin(p, temperature, bond_count, avg_order)

    rescale_velocities(particles, temperature)

    for p in particles:
        p.x += p.vx * dt
        p.y += p.vy * dt

        if p.x < 0:
            p.x += WORLD_WIDTH
        elif p.x >= WORLD_WIDTH:
            p.x -= WORLD_WIDTH

        if p.y < 0:
            p.y += WORLD_HEIGHT
        elif p.y >= WORLD_HEIGHT:
            p.y -= WORLD_HEIGHT
