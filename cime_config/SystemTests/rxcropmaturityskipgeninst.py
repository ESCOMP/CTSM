from rxcropmaturity import RXCROPMATURITYSHARED


class RXCROPMATURITYSKIPGENINST(RXCROPMATURITYSHARED):
    def run_phase(self):
        self._run_phase(skip_gen=True, h1_inst=True)
