import pygame
from physics import WORLD_WIDTH, WORLD_HEIGHT, DIPOLE_VISUAL_SCALE, wrap

TYPE_COLORS = [
    (60, 160, 255),
    (255, 90, 90),
    (120, 255, 120),
    (255, 240, 90),
    (180, 120, 255),
]

BOND_COLORS = {
    "covalent": (200, 200, 255),
    "ionic": (255, 180, 120),
    "metallic": (180, 180, 180),
}


class Visualizer:
    def __init__(self, sim, show_visual=True):
        self.sim = sim
        self.show_visual = show_visual

        pygame.init()
        if show_visual:
            self.screen = pygame.display.set_mode((WORLD_WIDTH, WORLD_HEIGHT))
        else:
            self.screen = pygame.Surface((WORLD_WIDTH, WORLD_HEIGHT))

        pygame.display.set_caption("Simulación Molecular Avanzada")
        self.clock = pygame.time.Clock()
        self.font = pygame.font.SysFont("Arial", 16)

    def toggle(self):
        self.show_visual = not self.show_visual
        print("Visualización:", "ON" if self.show_visual else "OFF")

    def _draw_particles(self, screen, particles):
        for p in particles:
            pygame.draw.circle(
                screen,
                TYPE_COLORS[p.type],
                (int(p.x), int(p.y)),
                4
            )

    def _draw_dipole_arrow(self, screen, x1, y1, x2, y2, bond):
        midx = (x1 + x2) / 2
        midy = (y1 + y2) / 2

        dpx, dpy = bond.dipole
        ax = midx + dpx * DIPOLE_VISUAL_SCALE
        ay = midy + dpy * DIPOLE_VISUAL_SCALE

        color = (255, 255, 255)
        pygame.draw.line(screen, color, (midx, midy), (ax, ay), 2)

        left_x = ax - dpx * 4 - dpy * 3
        left_y = ay - dpy * 4 + dpx * 3
        right_x = ax - dpx * 4 + dpy * 3
        right_y = ay - dpy * 4 - dpx * 3

        pygame.draw.line(screen, color, (ax, ay), (left_x, left_y), 2)
        pygame.draw.line(screen, color, (ax, ay), (right_x, right_y), 2)

    def _draw_bonds(self, screen, sim):
        particles = sim.get_particles()
        for bond in sim.get_bonds():
            p_i = particles[bond.i]
            p_j = particles[bond.j]

            dx = p_j.x - p_i.x
            dy = p_j.y - p_i.y
            dx, dy = wrap(dx, dy)

            x1, y1 = p_i.x, p_i.y
            x2, y2 = p_i.x + dx, p_i.y + dy

            width = max(1, bond.order * 2)
            color = BOND_COLORS.get(bond.type, (255, 255, 255))

            pygame.draw.line(screen, color, (x1, y1), (x2, y2), width)
            self._draw_dipole_arrow(screen, x1, y1, x2, y2, bond)

    def _draw_hud(self, screen, sim, fps):
        t = sim.temperature
        nb = len(sim.get_bonds())

        screen.blit(self.font.render(f"T = {t:.2f}", True, (255, 255, 255)), (10, 10))
        screen.blit(self.font.render(f"Enlaces = {nb}", True, (255, 255, 255)), (10, 30))
        screen.blit(self.font.render(f"FPS = {fps:.1f}", True, (255, 255, 255)), (10, 50))

    def update(self):
        if not self.show_visual:
            return

        screen = self.screen
        sim = self.sim
        screen.fill((0, 0, 0))

        self._draw_bonds(screen, sim)
        self._draw_particles(screen, sim.get_particles())

        fps = self.clock.get_fps()
        self._draw_hud(screen, sim, fps)

        pygame.display.flip()

    def tick(self, fps=60):
        return self.clock.tick(fps)

