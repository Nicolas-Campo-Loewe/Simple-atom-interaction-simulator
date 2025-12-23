import pygame
import sim_core as sc
from visual import Visualizer

INITIAL_TEMPERATURE = 0.75      # temperatura inicial
PARTICLES_PER_TYPE = 45         # partículas por tipo (5 tipos → 5 * 45)

def main():
    pygame.init()

    print(">>> creando Simulation")
    sim = sc.Simulation(
        n_types=5,
        particles_per_type=PARTICLES_PER_TYPE,
        temperature=INITIAL_TEMPERATURE,
    )

    print(">>> creando Visualizer")
    vis = Visualizer(sim, show_visual=True)

    running = True
    print("Controles: M toggle visual, ↑/↓ temperatura, ESC salir")

    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    running = False
                if event.key == pygame.K_m:
                    vis.toggle()
                if event.key == pygame.K_UP:
                    sim.increase_temperature(0.1)
                    print(f"T = {sim.temperature:.2f}")
                if event.key == pygame.K_DOWN:
                    sim.decrease_temperature(0.1)
                    print(f"T = {sim.temperature:.2f}")

        sim.step()
        vis.update()
        vis.tick(60)

    pygame.quit()

if __name__ == "__main__":
    main()
