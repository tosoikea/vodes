from pymbolic.primitives import Variable

class Noise(Variable):
    glb_index = 0

    def __init__(self,index=None):
        # TODO : Find better way for unique index -> overflow, race condition etc.
        if index is None:
            index = Noise.glb_index
            Noise.glb_index += 1

        self.index = index
        super().__init__(f"e_{index}")
        
    def __lt__(self, other):
        if other > 1:
            return True
        elif other <-1:
            return False
        else:
            return NotImplemented