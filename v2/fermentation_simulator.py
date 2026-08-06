import database
import kinetics
import numpy as np
import plotting
from reactor import Reactor
from yeast import Yeast

chemicals = database.chemicals


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
        self.mu_max = kinetics.rosso_cardinal(
            self.reactor.T_set,
            self.yeast.T_min,
            self.yeast.T_max,
            self.yeast.T_opt,
            self.yeast.mu_opt,
        )

        (
            self.biomass_concentration,
            self.sugar_concentration,
            self.ethanol_concentration,
            self.fermentation_time,
        ) = self.euler(biomass_concentration, sugar_concentration, self.mu_max)


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
        plotting.draw_fermentation_graph(self.dt,self.biomass_concentration,self.sugar_concentration,self.ethanol_concentration)
        self.print_status()
