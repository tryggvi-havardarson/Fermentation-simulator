from reactor import Reactor
from yeast import Yeast
from fermentation_simulator import FermentationSimulator


yeast1 = Yeast("yeast_proxy")
reactor1 = Reactor(20, 32)

fermentation1 = FermentationSimulator(reactor1, yeast1, 1000, 50, 5)

fermentation1.run()
