class Reactor:
    def __init__(self, volume: float, temperature: float, ph: float) -> None:
        self.volume = volume
        self.temperature = temperature
        self.ph = ph


class Yeast:
    def __init__(
        self, name: str, mu_max: float, ks: float, optimum_temperature: float
    ) -> None:
        self.name = name
        self.mu_max = mu_max
        self.ks = ks
        self.optimum_temperature = optimum_temperature

    def set_optimum_temperature(self, new_temperature) -> None:
        self.optimum_temperature = new_temperature

        print(f"Temperature has been set to {self.optimum_temperature}")


class Fermentation:
    def __init__(
        self, reactor: Reactor, yeast: Yeast, sugar: float, biomass: float
    ) -> None:
        self.reactor = reactor
        self.yeast = yeast
        self.sugar = sugar
        self.biomass = biomass
        self.running: bool = False

    def run(self) -> None:
        if not self.running:
            self.running = True
            print("Simulation is now running")
        else:
            print("Simulation is already running")

    def stop(self) -> None:
        if self.running:
            self.running = False
            print("Simulation has stopped")
        else:
            print("Simulation was already stopped")

    def print_status(self) -> None:
        print(f"""----- Fermentation Simulator -----

Reactor
Volume: {self.reactor.volume} L
Temperature: {self.reactor.temperature} °C
pH: {self.reactor.ph}

Yeast
Name: {self.yeast.name}
μmax: {self.yeast.mu_max} h⁻¹
Ks: {self.yeast.ks} g/L
Optimum temperature: {self.yeast.optimum_temperature} °C

Initial sugar: {self.sugar} g
Initial biomass: {self.biomass} g
        """)


yeast1 = Yeast("Saccharomyces cerevisiae", 0.45, 0.1, 30)
yeast2 = Yeast("Kluyveromyces marxianus", 0.7, 0.2, 42)

reactor1 = Reactor(20, 25, 7)
reactor2 = Reactor(10, 30, 10)
reactor3 = Reactor(5, 40, 3)

fermentation1 = Fermentation(reactor1, yeast1, 100, 5)
fermentation2 = Fermentation(reactor2, yeast1, 200, 5)
fermentation3 = Fermentation(reactor3, yeast2, 50, 2)

fermentation1.stop()
fermentation1.run()
fermentation1.run()
fermentation1.stop()

fermentation1.print_status()
