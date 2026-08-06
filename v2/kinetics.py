

def rosso_cardinal(T_set,T_min,T_max,T_opt,mu_opt) -> float:
    if T_set < T_min or T_set > T_max:
        mu_max = 0
    else:
        mu_max = (
            (mu_opt * (T_set - T_max))
            * (T_set - T_min) ** 2
        ) / (
            (T_opt - T_min)
            * (
                (T_opt - T_min)
                * (T_set - T_opt)
                - (T_opt - T_max)
                * (T_opt + T_min - 2 * T_set)
            )
        )

    return mu_max