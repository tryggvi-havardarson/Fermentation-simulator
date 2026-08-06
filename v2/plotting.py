import matplotlib.pyplot as plt
import numpy as np


def draw_fermentation_graph(dt, biomass_concentration, sugar_concentration, ethanol_concentration) -> tuple:

    time_array = np.arange(len(biomass_concentration)) * dt

    fig, ax = plt.subplots(figsize=(11, 6))

    ax.plot(
        time_array,
        biomass_concentration,
        linewidth=2.2,
        label="Biomass Concentration",
    )
    ax.plot(
        time_array,
        sugar_concentration,
        linewidth=2.2,
        label="Sugar Concentration",
    )
    ax.plot(
        time_array,
        ethanol_concentration,
        linewidth=2.2,
        label="Ethanol Concentration",
    )

    ax.set_title("Fermentation Process Profile", fontsize=15, fontweight="bold")
    ax.set_xlabel("Fermentation time (h)")
    ax.set_ylabel("Concentration (g/L)")

    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend()
    ax.margins(x=0)
    fig.tight_layout()

    plt.show()

    return fig, ax