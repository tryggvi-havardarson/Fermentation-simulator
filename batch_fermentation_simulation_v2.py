import numpy as np
import matplotlib.pyplot as plt

yeast_database = {
    "CBS 8066": {
        "mu_ref": 0.31,  # h^-1
        "Ks": 0.099,  # g/L assuming glucose
        "T_ref": 25,  # °C
        "Ea": 50000,  # J/mol
        "Y_xs": 0.10,  # g biomass / g glucose
    }
}

chemicals = {
    "glucose": {"chemical_formula": "C6H12O6", "molar_mass": 180.16},
    "ethanol": {"chemical_formula": "C2H5OH", "molar_mass": 46.07},
    "carbon_dioxide": {"chemical_formula": "CO2", "molar_mass": 44.01},
    "water": {"chemical_formula": "H2O", "molar_mass": 18.02},
}


def celsius_to_kelvin(temperature_in_celsius):
    return temperature_in_celsius + 273.15


class Reactor:
    def __init__(self, volume: float, T_set: float, ph: float) -> None:
        self.volume = volume
        self.T_set = T_set
        self.ph = ph


class Yeast:
    def __init__(self, name: str) -> None:

        if name not in yeast_database:
            raise ValueError(f"Unknown yeast strain: {name}")

        data = yeast_database[name]

        self.name = name
        self.mu_ref = data["mu_ref"]
        self.Ks = data["Ks"]
        self.T_ref = data["T_ref"]
        self.Ea = data["Ea"]
        self.Y_xs = data["Y_xs"]

    def summary(self):
        print(f"""
name: {self.name}
mu_ref: {self.mu_ref}
Ks: {self.Ks}
T_red: {self.T_ref}
Ea: {self.Ea}
Y_xs: {self.Y_xs}
        """)


class FermentationSimulator:
    def __init__(
        self, reactor: Reactor, yeast: Yeast, sugar: float, biomass: float, time: float
    ) -> None:
        self.reactor = reactor
        self.yeast = yeast
        self.sugar = sugar
        self.biomass = biomass
        self.time = time
        self.dt = float(0.001)
        self.duration_time = 0

    def mass_to_concentration(self, mass: float) -> float:
        return mass / self.reactor.volume

    def mu_max_function(self) -> float:
        R = 8.314
        T_ref_kelvin = celsius_to_kelvin(self.yeast.T_ref)
        T_set_kelvin = celsius_to_kelvin(self.reactor.T_set)

        return self.yeast.mu_ref * (
            np.exp((self.yeast.Ea / R) * ((1 / T_ref_kelvin) - (1 / T_set_kelvin)))
        )

    def euler(
        self,
        Xn,
        Sn,
        mu_max,
    ):
        Ks = self.yeast.Ks
        Y_XS = self.yeast.Y_xs
        dt = self.dt
        time = self.time
        En = 0
        S0 = Sn
        count = 0

        biomass = [Xn]
        sugar_concentration = [Sn]
        ethanol_concentration = [En]
        while count < (time / dt) and Sn > 1e-9:
            X = Xn + (dt * (Xn * mu_max * Sn)) / (Ks + Sn)
            S = Sn - (dt * (mu_max * Xn * Sn)) / (Y_XS * (Ks + Sn))

            Xn = X
            Sn = max(S, 0)

            glucose_to_ethanol = (S0 - Sn) * (1 - Y_XS)

            En = (
                2
                * glucose_to_ethanol
                * (
                    chemicals["ethanol"]["molar_mass"]
                    / chemicals["glucose"]["molar_mass"]
                )
            )
            biomass.append(Xn)
            sugar_concentration.append(Sn)
            ethanol_concentration.append(En)

            count += 1

        return (
            np.array(biomass),
            np.array(sugar_concentration),
            np.array(ethanol_concentration),
            count * dt,
        )

    def run(self) -> None:

        print("Simulation is now running")

        sugar_concentration = self.mass_to_concentration(self.sugar)
        biomass_concentration = self.mass_to_concentration(self.biomass)
        self.mu_max = self.mu_max_function()

        (
            self.biomass_concentration,
            self.sugar_concentration,
            self.ethanol_concentration,
            self.duration_time,
        ) = self.euler(biomass_concentration, sugar_concentration, self.mu_max)

    def draw_fermentation_graph(self):

        biomass_concentration = self.biomass_concentration
        sugar_concentration = self.sugar_concentration
        ethanol_concentration = self.ethanol_concentration
        dt = self.dt

        time = np.arange(len(biomass_concentration)) * dt
        biomass_concentration = np.asarray(biomass_concentration)
        sugar_concentration = np.asarray(sugar_concentration)
        ethanol_concentration = np.asarray(ethanol_concentration)

        fig, ax = plt.subplots(figsize=(11, 6))

        ax.plot(
            time,
            biomass_concentration,
            linewidth=2.2,
            label="Biomass Concentration",
        )
        ax.plot(
            time,
            sugar_concentration,
            linewidth=2.2,
            label="Sugar Concentration",
        )
        ax.plot(
            time,
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

    def print_status(self) -> None:
        print(f"""----- Fermentation Simulator -----

Reactor
Volume: {self.reactor.volume} L
Temperature: {self.reactor.T_set} °C
pH: {self.reactor.ph}

Yeast
Name: {self.yeast.name}
μmax: {self.mu_max} h⁻¹
Ks: {self.yeast.Ks} g/L

Initial sugar: {self.sugar} g
Initial biomass: {self.biomass} g

Duration time: {self.duration_time} hours
        """)


yeast1 = Yeast("CBS 8066")

reactor1 = Reactor(20, 25, 7)

fermentation1 = FermentationSimulator(reactor1, yeast1, 1000, 50, 5)

fermentation1.run()
fermentation1.draw_fermentation_graph()
fermentation1.print_status()
