import database

yeast_database = database.yeast_database

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