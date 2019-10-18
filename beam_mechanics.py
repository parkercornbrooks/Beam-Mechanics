import numpy as np
import matplotlib.pyplot as plt

class Beam:
    """
    A class used to represent a simple beam under loading

    ...

    Atrributes
    ----------
    length : float
        the length of the beam
    units : str
        the units of the length (default "meter")
    loads : list 
        list of loads applied to beam
    supports : list
        list of beam supports
    support_reactions : list
        list of support reactions
    
    Methods
    -------
    add_pointload(pos, force)
        adds a point load to the beam of given force at given position
    add_uniform_distload(start, end, force)
        adds a uniform distributed load from a start to end position of given force
    add_distLoad(start, end, startForce, endForce)
        add a distributed load with a given start/end position/force (units in force/length)
    add_moment(pos, moment)
        add a point moment at the given position
    add_fixed_support(location)
        add a fixed support at "left" or "right" side of beam
    add_pin_roller_support(pinPos, rollerPos)
        add a pin and a roller suport at the defined positions
    solve(position = None, N = 100)
        solve for the shear and moment values at the given positions
    """

    def __init__(self, length, units = "meter"):
        """
        Parameters
        ----------
        length : float
            the length of the beam
        units : str, optional
            the units of length (default is "meter")
        """
        self.length = length
        self.units = units
        self.loads = []
        self.supports = []
        self.support_reactions = []
    
    def __str__(self):
        return f"Beam with length {self.length}{self.units}"
    
    def add_pointLoad(self, pos, force):
        """Creates a point load and adds it to the loads list

        Parameters
        ----------
        pos : float
            the position on the beam where the point load is applied
        force : float
            the force of the point load
        """
        _check_location(self.length, pos)
        self.loads.append(PointLoad(pos, force))

    def add_uniform_distLoad(self, start, end, force):
        """Creates a uniform distributed load and adds it to the loads list

        Parameters
        ----------
        start : float
            the position on the beam where the distload starts
        end : float
            the position on the beam where the distload ends
        force : float
            the force, in force per unit length
        """
        _check_location(self.length, start, end)
        self.loads.append(DistLoad(start, end, force, force))
    
    def add_distLoad(self, start, end, startForce, endForce):
        """Creates a uniform distributed load and adds it to the loads list

        Parameters
        ----------
        start : float
            the position on the beam where the distload starts
        end : float
            the position on the beam where the distload ends
        startForce : float
            the force, in force per unit length, at the start position
        endForce : float
            the force, in force per unit length, at the end position
        """
        _check_location(self.length, start, end)
        self.loads.append(DistLoad(start, end, startForce, endForce))
    
    def add_moment(self, pos, moment):
        """Creates a point moment and adds it to the loads list

        Parameters
        ----------
        pos : float
            the position on the beam where the moment is applied
        moment : float
            the magnitude of the moment
        """
        _check_location(self.length, pos)
        self.loads.append(PointMoment(pos, moment))

    def add_fixed_support(self, location = "left"):
        """Creates a fixed support and adds it to the supports list

        Parameters
        ----------
        location : string
            the position of the fixed support (default "left")

        Returns
        -------
        False if beam already has supports
        True if supports are added
        """

        if len(self.supports) > 0:
            print("Beam already contains supports")
            return False
        elif location == "left":
            pos = 0
        elif location == "right":
            pos = self.length
        else:
            raise ValueError("Fixed support must be located at 'left' or 'right' of beam")
        self.supports.append(pos)
    
    def add_pin_roller_support(self, pinPos, rollerPos):
        """Creates a pin/roller support and adds it to the supports list

        Parameters
        ----------
        pinPos : float
            the position of the pin support
        rollerPos : float
            the position of the roller support

        Returns
        -------
        False if beam already has supports
        True if supports are added
        """

        if len(self.supports) > 0:
            print("Beam already contains supports")
            return False
        else:
            _check_location(self.length, pinPos, rollerPos)
            self.supports.extend([pinPos, rollerPos])
            return True

    def solve(self, positions = None, N = 100, pointsofinterest = True):
        """Solve for the shear and moment values at the given positions

        Parameters
        ----------
        positions : int/float, list of int/float, or numpy array of int/float, optional
            the position or positions to solve for shear and moment
            default is a linear spaced array in [0, beam length] with N points
        N : int, optional
            the number of positions to solve within the linear spaced array (default 100)
            if positions are given, this value is not used
        pointsofinterest : bool, optional
            whether to include more specificity at loading and support points
            only added if no positions are provided, default is True
        Returns
        -------
        False if there are no supports, or more than 2
        dictionary
            keys are "positions", "moment", "shear"
            values are numpy arrays
        """

        if positions == None: # create array of positions
            positions = np.linspace(0, self.length, N)
            if pointsofinterest:
                specificity = .00001
                pointsOI = list(self.supports)
                for load in self.loads:
                    if isinstance(load, DistLoad):
                        pointsOI.extend([load.start, load.end])
                    else:
                        pointsOI.append(load.pos)
                for point in pointsOI:
                    to_add = [point + specificity if np.isin(point, positions) else point, point + specificity]
                    idx = positions.searchsorted(point + specificity)
                    positions = np.concatenate((positions[:idx], to_add, positions[idx:]))

        else:  # convert list to numpy array of floats
            positions = np.array(positions, dtype=float)
        
        self.support_reactions = []

        # Solve for support reactions
        if len(self.supports) == 0:  # no supports
            print("Beam must be supported before it can be solved")
            return False
        elif len(self.supports) == 1:  # fixed support
            pos = self.supports[0]
            F_1 = -sum([load.force for load in self.loads])
            self.support_reactions.append(PointLoad(pos, F_1))

            M_1 = sum([load.force * (pos - load.pos) - load.moment for load in self.loads])
            self.support_reactions.append(PointMoment(pos, M_1))
        elif len(self.supports) == 2:  # pin and roller supports
            pinPos, rollerPos = self.supports
            F_roller = sum([load.force * (pinPos - load.pos) - load.moment for load in self.loads]) / (rollerPos - pinPos)
            self.support_reactions.append(PointLoad(rollerPos, F_roller))

            F_pin = -sum([load.force for load in self.loads]) - F_roller
            self.support_reactions.append(PointLoad(pinPos, F_pin))
        else:  # more than 2 supports??
            print("There are more than 2 supports?? Start over.")
            return False

        all_loads = self.loads + self.support_reactions

        shear = sum([load.calc_shear(positions) for load in all_loads])
        moment = sum([load.calc_moment(positions) for load in all_loads])
        return {"positions" : positions,
                "moment" : moment,
                "shear" : shear}

class PointLoad:
    """
    A class used to represent a Point Load

    ...

    Atrributes
    ----------
    pos : float
        the position of the load
    force : float
        the force of the load
    
    Methods
    -------
    calc_shear(locarray)
        returns the value of the shear due to load at all positions in the locarray
    calc_moment(locarray)
        returns the value of the moment due to load at all positions in the locarray
    """
    def __init__(self, pos, force):
        """
        Parameters
        ----------
        pos : float
            the position of the load
        force : float
            the magnitude and direction of the force
        """
        self.pos = pos
        self.force = force
        self.moment = 0

    def __str__(self):
        return f"Point Load of force {self.force} at {self.pos}m"

    def calc_shear(self, locarray):
        """ Calculates shear at defined points on the beam

        Parameters
        ----------
        locarray : numpy array
            array of positions on the beam

        Returns
        -------
        Array of shear due to this load at each position given
        """

        return np.where(locarray > self.pos, (self.force), 0)

    def calc_moment(self, locarray):
        """ Calculates bending moment at defined points on the beam

        Parameters
        ----------
        locarray : numpy array
            array of positions on the beam

        Returns
        -------
        Array of moments due to this load at each position given
        """

        return np.where(locarray > self.pos, (locarray-self.pos) * self.force, 0)

class DistLoad:
    """
    A class used to represent a Distributed Load

    ...

    Atrributes
    ----------
    start : float
        the position of start of the distributed load
    end : float
        the position of the end of the distributed load
    startForce : float
        the force at the start of the distributed load
    endForce : float
        the force at the end of the distributed load
    
    Methods
    -------
    _check_loc(start, end)
        ensures that start < end and they are not equal
    calc_shear(locarray)
        returns the value of the shear due to load at all positions in the locarray
    calc_moment(locarray)
        returns the value of the moment due to load at all positions in the locarray
    """
    def __init__(self, start, end, startForce, endForce):
        """
        Parameters
        ----------
        start : float
            the starting position of the load
        end : float
            the starting position of the load
        startForce : float
            the magnitude and direction of the force at the start of the load
        endForce : float
            the magnitude and direction of the force at the end of the load
        """
        self.start, self.end = self._check_loc(start, end)
        self.startForce = startForce
        self.endForce = endForce
        self.moment = 0
        (self.force, self.pos) = resultant(self.start, self.end, self.startForce, self.endForce)
    
    def __str__(self):
        return f"Distributed Load from {self.start}-{self.end}m, force ({self.startForce})-({self.endForce})"

    def _check_loc(self, start, end):
        """ Check the given start and end locations to ensure correct format

        Parameters
        ----------
        start : float
            given start position when initialized
        end : float
            given end position when initialized

        Returns
        -------
        tuple of properly ordered start and end values

        Raises
        ------
        ValueError if start and end values are equal
        """
        if start == end:
            raise ValueError("Start and end positions must be different")
        if start > end:
            start, end = end, start
        return start, end
    
    def calc_shear(self, locarray):
        """ Calculates shear at defined points on the beam

        Parameters
        ----------
        locarray : numpy array
            array of positions on the beam

        Returns
        -------
        Array of shear due to this load at each position given
        """

        sh = np.zeros_like(locarray)
        
        # interpolate value of distforce at all locations in locarray
        forceval = np.interp(locarray, [self.start, self.end],[self.startForce, self.endForce])
        vresultant = np.vectorize(resultant)
        (mag, dist) = vresultant(self.start, locarray, self.startForce, forceval)
        # find resultant at 
        sh = np.where(np.logical_and(locarray > self.start, locarray <= self.end), mag, sh)
        sh = np.where(locarray > self.end, self.force, sh)
        return sh

    def calc_moment(self, locarray):
        """ Calculates bending moment at defined points on the beam

        Parameters
        ----------
        locarray : numpy array
            array of positions on the beam

        Returns
        -------
        Array of moments due to this load at each position given
        """

        m = np.zeros_like(locarray)
        
        # interpolate value of distforce at all locations in locarray
        forceval = np.interp(locarray, [self.start, self.end],[self.startForce, self.endForce])
        vresultant = np.vectorize(resultant)
        (mag, dist) = vresultant(self.start, locarray, self.startForce, forceval)
        # find resultant at 
        m = np.where(np.logical_and(locarray > self.start, locarray <= self.end), mag * (locarray-dist), m)
        m = np.where(locarray > self.end, self.force * (locarray-self.pos), m)
        return m

class PointMoment:
    """
    A class used to represent a Point Moment

    ...

    Atrributes
    ----------
    pos : float
        the position of the moment
    moment : float
        the magnitude and direction of the moment
    
    Methods
    -------
    calc_shear(locarray)
        returns the value of the shear due to moment at all positions in the locarray
    calc_moment(locarray)
        returns the value of the moment due to this moment at all positions in the locarray
    """
    def __init__(self, pos, moment):
        """
        Parameters
        ----------
        pos : float
            the position of the moment
        moment : float
            the magnitude and direction of the moment
        """
        self.pos = pos
        self.force = 0
        self.moment = moment
    
    def __str__(self):
        return f"Point Moment of {self.moment} at {self.pos}m"

    def calc_shear(self, locarray):
        """ Calculates shear at defined points on the beam

        Parameters
        ----------
        locarray : numpy array
            array of positions on the beam

        Returns
        -------
        Array of shear due to this load at each position given
        """

        return np.zeros_like(locarray)

    def calc_moment(self, locarray):
        """ Calculates bending moment at defined points on the beam

        Parameters
        ----------
        locarray : numpy array
            array of positions on the beam

        Returns
        -------
        Array of moments due to this load at each position given
        """

        return np.where(locarray > self.pos, -self.moment , 0)


def resultant(x1, x2, y1, y2):
    ''' Calculates resultant force of a distributed load
    Can be used for uniform loads (y1 = y2), triangle (y1 or y2 = 0), or trapezoidal
    
    Parameters
    ----------
    x1 : float
        the start position of the load
    x2 : float
        the end position of the load
    y1 : float
        the start force of the load
    y2 : float
        the end force of the load

    Returns
    -------
    tuple 
        a tuple of resultant force and location (from position 0)
    '''

    ymin = min(y1, y2)
    b = (x2 - x1)
        
    F1 = ymin * b
    d1 = (x1 + x2) / 2
        
    h , d = (y1 - y2), 1/3
    if ymin < y2:
        h , d = (y2 - y1) , 2/3
        
    F2 = b * h / 2
    d2 = x1 + d * b

    R = F1 + F2
    dr = -1. if R == 0 else ((F1 * d1) + (F2 * d2)) / R
    return (R, dr)

def _check_location(length, *points):
    """ Checks all points to see if they fall within [0,length]
    
    Parameters
    ----------
    length : float
        the length of the beam
    *points : float
        variable number of points to be checked
    
    Raises
    ------
    Value Error if points fall outside of beam span
    """

    for point in points:
        if point < 0 or point > length:
            raise ValueError(f"Point {point} is not located on the beam of length 0-{length}")

def plot(solved_dict):
    """ Produce a plot of the shear and moment along the span of the beam

    Parameters
    ----------
    solved_dict : dict
        the return of the beam.solve() method
    """
    linecolor = 'blue'
    x, y1, y2 = solved_dict["positions"], solved_dict["shear"], solved_dict["moment"]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,12))
    
    ax1.plot(x, y1, color=linecolor, linewidth=2.5)
    ax1.fill_between(x, 0, y1, facecolor="gray", alpha=0.3)
    ax1.set_title("Shear", fontsize=18)
    ax1.set_ylabel("Shear (N)", fontsize=14)
    ax1.axhline(y=0, color='black')
    
    ax2.plot(x, y2, color=linecolor, linewidth=2.5)
    ax2.fill_between(x, 0, y2, facecolor="gray", alpha=0.3)
    ax2.set_title("Moment", fontsize=18)
    ax2.set_ylabel("Moment (N-m)", fontsize=14)
    ax2.set_xlabel("Position (m)", fontsize=14)
    ax2.axhline(y=0, color='black')
    plt.show()


def main():
    pass

if __name__ == "__main__":
    main()