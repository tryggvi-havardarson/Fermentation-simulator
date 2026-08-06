from fermentation_simulator import FermentationSimulator
from reactor import Reactor
from yeast import Yeast

yeast1 = Yeast("yeast_proxy")
reactor1 = Reactor(20, 30)

fermentation1 = FermentationSimulator(reactor1, yeast1, 1000, 50, 5)

fermentation1.run()
