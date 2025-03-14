from rxcropmaturity import RXCROPMATURITYSHARED


class RXCROPMATURITYINST(RXCROPMATURITYSHARED):
    def run_phase(self):
        self._run_phase(h1_inst=True)
