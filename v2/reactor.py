class Reactor:
    def __init__(self, volume: float, T_set: float) -> None:
        if volume <= 0:
            raise ValueError("Reactor volume must be greater than zero.")

        self.volume = volume
        self.T_set = T_set