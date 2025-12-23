import math
import random
from physics import BOND_BREAK_THRESHOLD, wrap

TARGET_SHELL = 5
VALENCE_ELECTRONS = [1, 2, 3, 4, 5] 

ELECTRONEG = [1.0, 1.6, 2.1, 2.6, 3.1]

TYPE_PREF = {
    0: {"max_bonds": 1, "mode": "terminal",   "prefer_multiple": False, "prefer_branch": False},
    1: {"max_bonds": 2, "mode": "linear",     "prefer_multiple": False, "prefer_branch": False},
    2: {"max_bonds": 3, "mode": "branched",   "prefer_multiple": False, "prefer_branch": True},
    3: {"max_bonds": 4, "mode": "multi",      "prefer_multiple": True,  "prefer_branch": True},
    4: {"max_bonds": 1, "mode": "inertish",   "prefer_multiple": False, "prefer_branch": False},
}

MAX_BONDS_BASE = 4

T_GEOM_MAX = 0.9
T_IONIC_MIN = 1.6

CAPTURE_RADIUS = 34.0
BOND_FORM_RADIUS = 26.0

REL_V2_BASE = 6.0
REL_V2_STABIL_BONUS = 8.0

P_BOND_BASE = 0.28             
P_BOND_DEFICIT_GAIN = 0.16
P_BOND_LOW_T_BOOST = 0.20
P_BOND_HIGH_T_BOOST = 0.10

P_DISSOC_BASE = 0.005
P_DISSOC_HIGH_T = 0.05


class Bond:
    def __init__(self, i, j, bond_type="covalent", order=1, length=17.0):
        self.i = i
        self.j = j
        self.type = bond_type
        self.order = order
        self.length = length
        self.dipole = (0.0, 0.0)
        self.donor = None
        self.acceptor = None
        self.age = 0  # soft-start en física

def shell_fill(p_idx, particles, bonds):
    p = particles[p_idx]
    fill = VALENCE_ELECTRONS[p.type]

    for b in bonds:
        if b.i != p_idx and b.j != p_idx:
            continue

        if b.type == "covalent":
            fill += b.order
        elif b.type == "ionic":
            if b.acceptor == p_idx:
                fill += b.order
            elif b.donor == p_idx:
                fill -= b.order
        elif b.type == "metallic":
            fill += 0.5 * b.order

    return fill


def shell_deficit(p_idx, particles, bonds):
    return max(0.0, TARGET_SHELL - shell_fill(p_idx, particles, bonds))


def dynamic_max_bonds(p_idx, particles, bonds):
    p = particles[p_idx]
    base = TYPE_PREF.get(p.type, {}).get("max_bonds", MAX_BONDS_BASE)
    d = shell_deficit(p_idx, particles, bonds)
    # si tiene déficit alto, puede temporalmente superar el base (hasta 5)
    extra = int(math.ceil(d * 0.7))
    return max(0, min(5, base + extra))

def classify_bond_type_temp(p_i, p_j, temperature):
    diff = abs(ELECTRONEG[p_i.type] - ELECTRONEG[p_j.type])

    if temperature <= T_GEOM_MAX:
        return "covalent" if diff < 1.6 else "ionic"

    if temperature >= T_IONIC_MIN:
        return "ionic" if diff > 0.9 else "metallic"

    if diff < 0.35:
        return "metallic"
    elif diff < 1.2:
        return "covalent"
    else:
        return "ionic"

def dipole_from_EN(p_i, p_j):
    e1, e2 = ELECTRONEG[p_i.type], ELECTRONEG[p_j.type]
    dx = p_j.x - p_i.x
    dy = p_j.y - p_i.y
    dx, dy = wrap(dx, dy)
    r2 = dx*dx + dy*dy
    if r2 < 1e-9:
        return (0.0, 0.0)
    r = math.sqrt(r2)
    nx, ny = dx/r, dy/r
    return (nx, ny) if e2 > e1 else (-nx, -ny)

def prebond_damp(i, j, particles, bonds, temperature):
    p_i = particles[i]
    p_j = particles[j]

    di = shell_deficit(i, particles, bonds)
    dj = shell_deficit(j, particles, bonds)
    if di <= 0 and dj <= 0:
        return

    dx = p_j.x - p_i.x
    dy = p_j.y - p_i.y
    dx, dy = wrap(dx, dy)
    r2 = dx*dx + dy*dy

    if r2 > CAPTURE_RADIUS * CAPTURE_RADIUS:
        return

    deficit_factor = min(1.0, (di + dj) / 5.0)
    temp_factor = 1.0 if temperature <= T_GEOM_MAX else 0.6
    capture = 0.25 * deficit_factor * temp_factor

    rvx = p_j.vx - p_i.vx
    rvy = p_j.vy - p_i.vy
    p_i.vx += rvx * capture * 0.5
    p_i.vy += rvy * capture * 0.5
    p_j.vx -= rvx * capture * 0.5
    p_j.vy -= rvy * capture * 0.5

def try_form_bond(i, j, particles, bonds, temperature):
    p_i = particles[i]
    p_j = particles[j]

    if len(p_i.bonds) >= dynamic_max_bonds(i, particles, bonds):
        return False
    if len(p_j.bonds) >= dynamic_max_bonds(j, particles, bonds):
        return False

    di = shell_deficit(i, particles, bonds)
    dj = shell_deficit(j, particles, bonds)
    if di <= 0 and dj <= 0:
        return False

    dx = p_j.x - p_i.x
    dy = p_j.y - p_i.y
    dx, dy = wrap(dx, dy)
    r2 = dx*dx + dy*dy
    if r2 > BOND_FORM_RADIUS * BOND_FORM_RADIUS:
        return False

    rvx = p_j.vx - p_i.vx
    rvy = p_j.vy - p_i.vy
    relv2 = rvx*rvx + rvy*rvy

    stabil_bonus = REL_V2_STABIL_BONUS * (di + dj) / 3.5
    relv2_limit = REL_V2_BASE * temperature + stabil_bonus
    if relv2 > relv2_limit:
        return False

    pref_i = TYPE_PREF.get(p_i.type, {})
    pref_j = TYPE_PREF.get(p_j.type, {})

    p_form = P_BOND_BASE + P_BOND_DEFICIT_GAIN * (di + dj)

    if temperature <= T_GEOM_MAX:
        p_form += P_BOND_LOW_T_BOOST
    elif temperature >= T_IONIC_MIN:
        p_form += P_BOND_HIGH_T_BOOST

    if pref_i.get("prefer_branch") or pref_j.get("prefer_branch"):
        p_form += 0.08

    if random.random() > min(0.97, p_form):
        return False

    bond_type = classify_bond_type_temp(p_i, p_j, temperature)

    order = 1
    prefer_multi = pref_i.get("prefer_multiple") or pref_j.get("prefer_multiple")
    if bond_type == "covalent":
        if temperature <= T_GEOM_MAX and (di + dj) >= 3:
            base_prob2 = 0.35
            if prefer_multi:
                base_prob2 += 0.20
            if random.random() < base_prob2:
                order = 2

        if temperature <= T_GEOM_MAX and (di + dj) >= 5 and order == 2:
            base_prob3 = 0.18 if prefer_multi else 0.08
            if random.random() < base_prob3:
                order = 3

    r = math.sqrt(r2)
    if bond_type == "covalent":
        length = max(12.5, min(18.0, r))
    elif bond_type == "ionic":
        length = max(14.0, min(21.0, r))
    else:
        length = max(15.0, min(22.0, r))

    b = Bond(i, j, bond_type=bond_type, order=order, length=length)
    b.dipole = dipole_from_EN(p_i, p_j)

    if bond_type == "ionic":
        if ELECTRONEG[p_i.type] > ELECTRONEG[p_j.type]:
            b.acceptor, b.donor = i, j
        else:
            b.acceptor, b.donor = j, i

    bonds.append(b)
    p_i.bonds.append(j)
    p_j.bonds.append(i)
    return True

def try_break_bond(bond, particles, temperature):
    p_i = particles[bond.i]
    p_j = particles[bond.j]

    dx = p_j.x - p_i.x
    dy = p_j.y - p_i.y
    dx, dy = wrap(dx, dy)
    r2 = dx*dx + dy*dy
    r = math.sqrt(r2)

    if r > bond.length * BOND_BREAK_THRESHOLD * (1.0 + 0.2 * (bond.order - 1)):
        return True

    rvx = p_j.vx - p_i.vx
    rvy = p_j.vy - p_i.vy
    relv2 = rvx*rvx + rvy*rvy

    base = 18.0 if bond.type == "covalent" else 12.0
    if temperature >= T_IONIC_MIN and bond.type != "covalent":
        base *= 0.7

    cluster_factor = 1.0 + 0.2 * min(len(p_i.bonds), len(p_j.bonds))
    base *= cluster_factor

    if relv2 > base * temperature * bond.order:
        return True

    return False

def reaction_addition(i, j, bonds, particles, temperature):
    p_i = particles[i]
    p_j = particles[j]
    if (len(p_i.bonds) > 0) ^ (len(p_j.bonds) > 0):
        return try_form_bond(i, j, particles, bonds, temperature)
    return False


def reaction_substitution(i, j, bonds, particles, temperature):
    if temperature < 1.05:
        return False

    p_i = particles[i]
    p_j = particles[j]
    if len(p_i.bonds) == 0:
        return False

    e_j = ELECTRONEG[p_j.type]

    weakest_k = None
    weakest_en = 999
    for k in p_i.bonds:
        e_k = ELECTRONEG[particles[k].type]
        if e_k < weakest_en:
            weakest_en = e_k
            weakest_k = k

    if weakest_k is None:
        return False

    if e_j > weakest_en + 0.10:
        remove_bond_between(i, weakest_k, bonds, particles)
        return try_form_bond(i, j, particles, bonds, temperature)

    return False


def reaction_dissociation(i, bonds, particles, temperature):
    p_i = particles[i]
    if len(p_i.bonds) <= 1:
        return False

    prob = P_DISSOC_BASE
    if temperature >= T_IONIC_MIN:
        prob += P_DISSOC_HIGH_T * (temperature - T_IONIC_MIN)

    # cuanto más enlaces tenga ese átomo, más protegido (más difícil disociar)
    prob /= (1.0 + 0.6 * (len(p_i.bonds) - 1))

    if random.random() < min(0.5, prob):
        k = random.choice(p_i.bonds)
        remove_bond_between(i, k, bonds, particles)
        return True
    return False


def remove_bond_between(i, j, bonds, particles):
    for idx, b in enumerate(bonds):
        if (b.i == i and b.j == j) or (b.i == j and b.j == i):
            bonds.pop(idx)
            break

    if j in particles[i].bonds:
        particles[i].bonds.remove(j)
    if i in particles[j].bonds:
        particles[j].bonds.remove(i)
