class Empty_Mangle(object):
    def __init__(self):
        pass

    def weight(self, ra, dec):
        return 0.0

    def polyid(self, ra, dec):
        return [-1] * len(ra)

    def polyid_and_weight(self, ra, dec):
        return tuple([[-1, 0.0]] * len(ra))
