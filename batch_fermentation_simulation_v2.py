import numpy as np
import matplotlib.pyplot as plt

yeast_strains = {
    "CBS 8066": {
        "mu_reference": 0.31,  # h^-1
        "temperature_reference": 25,  # °C
        "activation_energy": 50000,  # J/mol
        "Ks": 0.099,  # g/L assuming glucose
        "Y_X/S": 0.10,  # g biomass / g glucose
        "maximum_temperature": 40,  # °C
    }
}

chemicals = {
    "glucose": {"chemical_formula": "C6H12O6", "molar_mass": 180.16},
    "ethanol": {"chemical_formula": "C2H5OH", "molar_mass": 46.07},
    "carbon_dioxide": {"chemical_formula": "CO2", "molar_mass": 44.01},
    "water": {"chemical_formula": "H2O", "molar_mass": 18.02},
}


def draw_fermentation_graph(
    biomass_concentration,
    sugar_concentration,
    ethanol_concentration,
    dt,
):

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


def celsius_to_kelvin(temperature_in_celsius):
    return temperature_in_celsius + 273.15


class Reactor:
    def __init__(self, volume: float, T_set: float, ph: float) -> None:
        self.volume = volume
        self.T_set = T_set
        self.ph = ph


class Yeast:
    def __init__(
        self,
        name: str,
        mu_ref: float,
        Ks: float,
        optimum_temperature: float,
        T_ref: float,
        Ea: float,
        Y_xs: float,
    ) -> None:
        self.name = name
        self.mu_ref = mu_ref
        self.Ks = Ks
        self.optimum_temperature = optimum_temperature
        self.T_ref = T_ref
        self.Ea = Ea
        self.Y_xs = Y_xs

    def set_optimum_temperature(self, new_temperature) -> None:
        self.optimum_temperature = new_temperature

        print(f"Temperature has been set to {self.optimum_temperature}")


class Fermentation:
    def __init__(
        self, reactor: Reactor, yeast: Yeast, sugar: float, biomass: float, time: float
    ) -> None:
        self.reactor = reactor
        self.yeast = yeast
        self.sugar = sugar
        self.biomass = biomass
        self.time = time
        self.dt = float(0.001)

    def mass_to_concentration(self, mass) -> float:
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
        mu_max = self.mu_max_function()

        biomass_concentration, sugar_concentration, ethanol_concentration, self.time = (
            self.euler(biomass_concentration, sugar_concentration, mu_max)
        )

        draw_fermentation_graph(
            biomass_concentration, sugar_concentration, ethanol_concentration, self.dt
        )

    def print_status(self) -> None:
        print(f"""----- Fermentation Simulator -----

Reactor
Volume: {self.reactor.volume} L
Temperature: {self.reactor.T_set} °C
pH: {self.reactor.ph}

Yeast
Name: {self.yeast.name}
μmax: {self.yeast.mu_ref} h⁻¹
Ks: {self.yeast.Ks} g/L
Optimum temperature: {self.yeast.optimum_temperature} °C

Initial sugar: {self.sugar} g
Initial biomass: {self.biomass} g

Time: {self.time} hours
        """)


yeast1 = Yeast("Saccharomyces cerevisiae", 0.45, 0.1, 30, 25, 50000, 0.1)

reactor1 = Reactor(20, 25, 7)

fermentation1 = Fermentation(reactor1, yeast1, 100, 5, 5)

fermentation1.run()
fermentation1.print_status()
