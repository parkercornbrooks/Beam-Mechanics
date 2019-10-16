import numpy as np

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
    add_support(style, pos)
        add a support of given style at the given position
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
    
    def __str__(self):
        s = f"Beam with length {self.length}{self.units}"
        #s += f"\n\t Supports: "
        #s += "None" if len(self.supports) == 0 else str([str(support) for support in self.supports])
        return s
    
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

    def add_support(self, style, pos):
        """Creates a support and adds it to the supports list

        Parameters
        ----------
        style : str
            the style of support (pin, roller, fixed)
        pos : float
            the position on the beam where the moment is applied
        """
        _check_location(self.length, pos)
        check = True
        if len(self.supports) == 2:
            print("Adding more than 2 supports is not allowed")
            check = False
        elif len(self.supports) == 1:
            if self.supports[0].unknowns == 3:
                print("Beam is already fixed. Adding more supports will cause the problem to be indeterminate.")
                check = False
            elif style == "fixed":
                print("Cannot add a fixed support to a beam that already has a support.")
                check = False
            elif self.supports[0].pos == pos:
                print(f"There is already a support at {pos}.")
                check = False
        
        if check:
            if style == "roller":
                self.supports.append(RollerSupport(pos))
                print(f"Roller support added at {pos}.")
            elif style == "pin":
                self.supports.append(PinSupport(pos))
                print(f"Pin support added at {pos}.")
            elif style == "fixed":
                if pos not in [0, self.length]:
                    print(f"Fixed supports can only be added at the end of a beam (0 or {self.length})")
                else:
                    self.supports.append(FixedSupport(pos))
                    print(f"Fixed support added at {pos}.")

    def solve(self, positions = None, N = 100):
        """Solve for the shear and moment values at the given positions

        Parameters
        ----------
        positions : int/float, list of int/float, or numpy array of int/float, optional
            the position or positions to solve for shear and moment
            default is a linear spaced array in [0, beam length] with N points
        N : int, optional
            the number of positions to solve within the linear spaced array (default 100)
            if positions are given, this value is not used

        Returns
        -------
        dictionary
            keys are "positions", "moment", "shear"
            values are numpy arrays
        """

        if positions == None:
            positions = np.linspace(0, self.length, N)
        else:  # convert list to numpy array of floats
            positions = np.array(positions, dtype=float)
        
        # ToDo: solve for support reactions
        
        shear = sum([load.calc_shear(positions) for load in self.loads])
        moment = sum([load.calc_moment(positions) for load in self.loads])
        return {"positions" : positions,
                "moment" : moment,
                "shear" : shear}

class PointLoad:
    def __init__(self, pos, force):
        self.pos = pos
        self.force = force

    def calc_shear(self, locarray):
        '''takes in a numpy array of locations
        returns array of shear due to load at each (from L to R)'''

        return np.where(locarray > self.pos, (self.pos), 0)

    def calc_moment(self, locarray):
        '''takes in a numpy array of locations
        returns array of moments due to load at each (from L to R)'''

        return np.where(locarray > self.pos, (locarray-self.pos) * self.force, 0)

class DistLoad:
    def __init__(self, start, end, startForce, endForce):
        self.start, self.end = self._check_loc(start, end)
        self.startForce = startForce
        self.endForce = endForce
        (self.force, self.pos) = self.resultant()
    
    def _check_loc(self, start, end):
        if start == end:
            raise ValueError("Start and end positions must be different")
        if start > end:
            start, end = end, start
        return start, end

    def resultant(self):
        return resultant(self.start, self.end, self.startForce, self.endForce)
    
    def calc_shear(self, locarray):
        '''takes in a numpy array of locations
        returns array of shear due to load at each (from L to R)'''

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
        '''takes in a numpy array of locations
        returns array of moments due to load at each (from L to R)'''

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
    def __init__(self, pos, moment):
        self.pos = pos
        self.moment = moment

    def calc_shear(self, locarray):
        '''takes in a numpy array of locations
        returns array of shear due to load at each (from L to R)'''

        return np.zeros_like(locarray)

    def calc_moment(self, locarray):
        '''takes in a numpy array of locations
        returns array of moments due to load at each (from L to R)'''

        return np.where(locarray > self.pos, self.moment , 0)

class Support:
    def __init__(self, style, pos, unknowns):
        self.pos = pos
        self.unknowns = unknowns
        self.style = style
    
    def __str__(self):
        return f"{self.style} support @ {self.pos}"
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.pos == other.pos
        return False

class FixedSupport(Support):
    def __init__(self, pos):
        super().__init__(style = "Fixed", pos = pos, unknowns = 3)

class RollerSupport(Support):
    def __init__(self, pos):
        super().__init__(style = "Roller", pos = pos, unknowns = 1)

class PinSupport(Support):
    def __init__(self, pos):
        super().__init__(style = "Pin", pos = pos, unknowns = 2)


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
    
    Raises ValueError if any point falls out of bounds
    """
    for point in points:
        if point < 0 or point > length:
            raise ValueError(f"Point {point} is not located on the beam of length 0-{length}")



def main():
    pass

if __name__ == "__main__":
    main()