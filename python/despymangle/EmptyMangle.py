class Empty_Mangle:
    """ Class to mimic a Mangle object when there are 0 input polygons

        The class takes no arguments

    """
    def __init__(self):
        pass

    def weight(self, ra, dec):
        """ Method to return the default weight (0.0)

            Parameters
            ----------
            ra - float
                The RA of the point in degrees.

            dec - float
                The DEC of the point in degrees.

            Returns
            -------
            0.0

        """
        return 0.0

    def polyid(self, ra, dec):
        """ Method to return the poly ids an array of -1's

            Parameters
            ----------
            ra - list of floats
                A list of RA points

            dec - list of floats
                A list of DEC points

            Returns
            -------
            A list of the same length as the input lists, filled with -1

        """
        return [-1] * len(ra)

    def polyid_and_weight(self, ra, dec):
        """ Method to return a list of poly ids and weights (-1 and 0.0)

            Parameters
            ----------

            ra - list of floats
                A list of RA points

            dec - list of floats
                A list of DEC points

            Returns
            -------
            A tuple of the same length as the input lists, filled with [-1, 0.0]

        """
        return tuple([[-1, 0.0]] * len(ra))
