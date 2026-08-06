import matplotlib.pyplot as plt
import numpy as np

yeast_database = {
    "yeast_proxy": {
        "Ks": 0.18,  # g/L glucose
        "Y_xs": 0.097,  # g dry biomass / g glucose
        "T_min": 3.08,  # °C
        "T_opt": 30.03,  # °C
        "T_max": 41.21,  # °C
        "mu_opt": 0.368,  # h^-1
    }
}

chemicals = {
    "glucose": {"chemical_formula": "C6H12O6", "molar_mass": 180.16},
    "ethanol": {"chemical_formula": "C2H5OH", "molar_mass": 46.07},
    "carbon_dioxide": {"chemical_formula": "CO2", "molar_mass": 44.01},
    "water": {"chemical_formula": "H2O", "molar_mass": 18.02},
}


class Reactor:
    def __init__(self, volume: float, T_set: float) -> None:
        if volume <= 0:
            raise ValueError("Reactor volume must be greater than zero.")

        self.volume = volume
        self.T_set = T_set


class Yeast:
    def __init__(self, name: str) -> None:

        if name not in yeast_database:
            raise ValueError(f"Unknown yeast strain: {name}")

        data = yeast_database[name]

        self.name = name

        self.Ks = data["Ks"]
        self.Y_xs = data["Y_xs"]
        self.T_min = data["T_min"]
        self.T_opt = data["T_opt"]
        self.T_max = data["T_max"]
        self.mu_opt = data["mu_opt"]


class FermentationSimulator:
    def __init__(
        self,
        reactor: Reactor,
        yeast: Yeast,
        sugar_mass: float,
        biomass_mass: float,
        simulation_time: float,
    ) -> None:
        if sugar_mass < 0:
            raise ValueError("Sugar mass cannot be negative.")

        if biomass_mass < 0:
            raise ValueError("Biomass cannot be negative.")

        if simulation_time <= 0:
            raise ValueError("Simulation time must be greater than zero.")

        self.reactor = reactor
        self.yeast = yeast
        self.sugar_mass = sugar_mass
        self.biomass_mass = biomass_mass
        self.simulation_time = simulation_time

        self.dt = 0.001
        self.fermentation_time = 0

        self.mu_max = None
        self.biomass_concentration = None
        self.sugar_concentration = None
        self.ethanol_concentration = None

    def mass_to_concentration(self, mass: float) -> float:
        return mass / self.reactor.volume

    def rosso_cardinal(self) -> float:
        if (
            self.reactor.T_set < self.yeast.T_min
            or self.reactor.T_set > self.yeast.T_max
        ):
            mu_max = 0
        else:
            mu_max = (
                (self.yeast.mu_opt * (self.reactor.T_set - self.yeast.T_max))
                * (self.reactor.T_set - self.yeast.T_min) ** 2
            ) / (
                (self.yeast.T_opt - self.yeast.T_min)
                * (
                    (self.yeast.T_opt - self.yeast.T_min)
                    * (self.reactor.T_set - self.yeast.T_opt)
                    - (self.yeast.T_opt - self.yeast.T_max)
                    * (self.yeast.T_opt + self.yeast.T_min - 2 * self.reactor.T_set)
                )
            )

        return mu_max

    def euler(self, Xn: float, Sn: float, mu_max: float) -> tuple:
        En = 0
        S0 = Sn
        count = 0

        biomass_mass = [Xn]
        sugar_concentration = [Sn]
        ethanol_concentration = [En]

        while count < (self.simulation_time / self.dt) and Sn > 1e-9:
            X = Xn + (self.dt * (Xn * mu_max * Sn)) / (self.yeast.Ks + Sn)
            S = Sn - (self.dt * (mu_max * Xn * Sn)) / (
                self.yeast.Y_xs * (self.yeast.Ks + Sn)
            )

            Xn = X
            Sn = max(S, 0)

            glucose_to_ethanol = (S0 - Sn) * (1 - self.yeast.Y_xs)

            En = (
                2
                * glucose_to_ethanol
                * (
                    chemicals["ethanol"]["molar_mass"]
                    / chemicals["glucose"]["molar_mass"]
                )
            )

            biomass_mass.append(Xn)
            sugar_concentration.append(Sn)
            ethanol_concentration.append(En)

            count += 1

        return (
            np.array(biomass_mass),
            np.array(sugar_concentration),
            np.array(ethanol_concentration),
            count * self.dt,
        )

    def prepare(self) -> None:

        sugar_concentration = self.mass_to_concentration(self.sugar_mass)
        biomass_concentration = self.mass_to_concentration(self.biomass_mass)
        self.mu_max = self.rosso_cardinal()

        (
            self.biomass_concentration,
            self.sugar_concentration,
            self.ethanol_concentration,
            self.fermentation_time,
        ) = self.euler(biomass_concentration, sugar_concentration, self.mu_max)

    def draw_fermentation_graph(self) -> tuple:

        time_array = np.arange(len(self.biomass_concentration)) * self.dt

        fig, ax = plt.subplots(figsize=(11, 6))

        ax.plot(
            time_array,
            self.biomass_concentration,
            linewidth=2.2,
            label="Biomass Concentration",
        )
        ax.plot(
            time_array,
            self.sugar_concentration,
            linewidth=2.2,
            label="Sugar Concentration",
        )
        ax.plot(
            time_array,
            self.ethanol_concentration,
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
        print(f"""
        ======================================================================
                                FERMENTATION SIMULATION
        ======================================================================

        Reactor
        ----------------------------------------------------------------------
        Volume                : {self.reactor.volume} L
        Temperature           : {self.reactor.T_set} °C


        Yeast
        ----------------------------------------------------------------------
        Strain                : {self.yeast.name}
        Maximum growth (μmax) : {self.mu_max:.3f} h⁻¹
        Ks                    : {self.yeast.Ks} g/L


        Simulation
        ----------------------------------------------------------------------
        Requested time        : {self.simulation_time:.2f} h
        Fermentation time     : {self.fermentation_time:.2f} h
        Sugar conversion      : {((self.sugar_mass - self.sugar_concentration[-1]) / self.sugar_mass) * 100:.1f} %


        Component Mass Balance
        ----------------------------------------------------------------------
        Component              Initial (g)        Final (g)
        ----------------------------------------------------------------------
        Biomass             {self.biomass_mass:12.2f}   {self.reactor.volume * self.biomass_concentration[-1]:12.2f}
        Sugar               {self.sugar_mass:12.2f}   {self.reactor.volume * self.sugar_concentration[-1]:12.2f}
        Ethanol             {0:12.2f}   {self.reactor.volume * self.ethanol_concentration[-1]:12.2f}
        Carbon dioxide      {0:12.2f}   {((self.reactor.volume * self.ethanol_concentration[-1]) / (chemicals["ethanol"]["molar_mass"])) * chemicals["carbon_dioxide"]["molar_mass"]:12.2f}
        ----------------------------------------------------------------------
        """)

    def run(self) -> None:
        print("Simulation is now running")

        self.prepare()
        self.draw_fermentation_graph()
        self.print_status()


yeast1 = Yeast("yeast_proxy")
reactor1 = Reactor(20, 32)

fermentation1 = FermentationSimulator(reactor1, yeast1, 1000, 50, 5)

fermentation1.run()
