from rxcropmaturity import RXCROPMATURITYSHARED


class RXCROPMATURITYSKIPGEN(RXCROPMATURITYSHARED):

    def run_phase(self):
        self._run_phase(skip_run=False, skip_gen=True)
