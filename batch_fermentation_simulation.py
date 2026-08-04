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


def mass_to_concentration(mass, volume):
    concentration = mass / volume
    return concentration


def celsius_to_kelvin(temperature_in_celsius):
    temperature_in_kelvin = temperature_in_celsius + 273.15
    return temperature_in_kelvin


def mu_max_function(mu_ref, T_ref, Ea, T_set):
    R = 8.314
    mu_max = mu_ref * (np.exp((Ea / R) * ((1 / T_ref) - (1 / T_set))))

    return mu_max


def euler(yeast_strain, time, Xn, Sn, mu_max, dt):
    yeast_strain = str(yeast_strain)
    Ks = yeast_strains[yeast_strain]["Ks"]
    Y_XS = yeast_strains[yeast_strain]["Y_X/S"]
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
            * (chemicals["ethanol"]["molar_mass"] / chemicals["glucose"]["molar_mass"])
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


def main():
    # Initial conditions
    dt = 0.001  # h
    initial_glucose_mass = 5000  # g
    initial_biomass = 50  # g
    volume_in_liters = 15  # L
    yeast_strain = "CBS 8066"
    time_in_hours = 72  # h
    temperature_in_celsius = 31  # °C

    if yeast_strain not in yeast_strains:
        print(f"No data available for strain: {yeast_strain}")
        return

    initial_glucose_concentration = mass_to_concentration(
        initial_glucose_mass, volume_in_liters
    )
    initial_biomass_concentration = mass_to_concentration(
        initial_biomass, volume_in_liters
    )

    mu_max = mu_max_function(
        yeast_strains[yeast_strain]["mu_reference"],
        celsius_to_kelvin(yeast_strains[yeast_strain]["temperature_reference"]),
        yeast_strains[yeast_strain]["activation_energy"],
        celsius_to_kelvin(temperature_in_celsius),
    )

    if temperature_in_celsius < yeast_strains[yeast_strain]["maximum_temperature"]:
        biomass_concentration, sugar_concentration, ethanol_concentration, time = euler(
            yeast_strain,
            time_in_hours,
            initial_biomass_concentration,
            initial_glucose_concentration,
            mu_max,
            dt,
        )

        final_biomass = biomass_concentration[-1] * volume_in_liters
        final_glucose_mass = sugar_concentration[-1] * volume_in_liters
        final_ethanol_mass = ethanol_concentration[-1] * volume_in_liters
        final_carbondioxide_mass = (
            final_ethanol_mass / (chemicals["ethanol"]["molar_mass"])
        ) * chemicals["carbon_dioxide"]["molar_mass"]

        draw_fermentation_graph(
            biomass_concentration, sugar_concentration, ethanol_concentration, dt
        )

        print(f"""
        ============================================================
                            Simulation Summary
        ============================================================

        Yeast strain:           {yeast_strain}
        Fermentation time:      {time} h
        Temperature:            {temperature_in_celsius:.1f} °C

        Maximum growth rate:    {mu_max:.3f} h⁻¹
        Sugar conversion:       {(1 - sugar_concentration[-1] / initial_glucose_concentration) * 100:.1f} %

                           Component Mass Balance
        ------------------------------------------------------------
        Component          | Initial (g) | Final (g)
        ------------------------------------------------------------
        Biomass            | {initial_biomass:11.2f} | {final_biomass:11.2f}
        Sugar              | {initial_glucose_mass:11.2f} | {final_glucose_mass:11.2f}
        Ethanol            | {0:11.2f} | {final_ethanol_mass:11.2f}
        Carbon dioxide     | {0:11.2f} | {final_carbondioxide_mass:11.2f}
        ------------------------------------------------------------
        """)
    else:
        print("Temperature is too high")
        print(
            f"Maximum temperature for strain is: {yeast_strains[yeast_strain]['maximum_temperature']}"
        )


if __name__ == "__main__":
    main()
