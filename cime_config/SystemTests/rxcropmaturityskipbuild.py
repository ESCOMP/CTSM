from rxcropmaturityshared import RXCROPMATURITYSHARED

class RXCROPMATURITYSKIPBUILD(RXCROPMATURITYSHARED):
    def build_indv(self, sharedlib_only=False, model_only=False):
        self._case.set_value("BUILD_COMPLETE", "TRUE")


    def run_phase(self):
        self._run_phase(skip_run=False)