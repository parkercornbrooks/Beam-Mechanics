import numpy as np

class Beam:
    def __init__(self, length, units = "m"):
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
        if pos < 0 or pos > self.length:
            print(f"Force must be within the span of the beam (0-{self.length}m)")
        else:
            self.loads.append(PointLoad(pos, force))

    def add_uniform_distLoad(self, start, end, force):
        self.loads.append(DistLoad(start, end, force, force))
    
    def add_distLoad(self, start, end, startForce, endForce):
        self.loads.append(DistLoad(start, end, startForce, endForce))
    
    def add_moment(self, pos, moment):
        self.loads.append(PointMoment(pos, moment))

    def add_support(self, style, pos):
        check = True
        if len(self.supports) == 2:
            print("Adding more than 2 supports is not allowed")
            check = False
        elif pos < 0 or pos > self.length:
            print(f"Support must be within the span of the beam (0-{self.length}m)")
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
        
        if positions == None:
            positions = np.linspace(0, self.length, N)
        else:
            positions = np.array(positions, dtype=float)
        # ToDo: solve for support reactions
        
        shear = sum([load.calc_shear(positions) for load in self.loads])
        moment = sum([load.calc_moment(positions) for load in self.loads])
        return {"positions" : positions,
                "moment" : moment,
                "shear" : shear}

def is_indeterminate(beam):
    unknowns = sum([support.unknowns for support in beam.supports])
    return True if unknowns > 3 else False

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
        self.start = start
        self.end = end
        self.startForce = startForce
        self.endForce = endForce
        (self.force, self.pos) = self.resultant()

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
    ''' returns a tuple with resultant force and location
    holds true for uniform loads, triangle, and trapezoid'''

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

def main():
    #resultant_test()
    solve_test()
    
    
def resultant_test():
    L = DistLoad(0, 5, 5, 0)
    x = np.linspace(0, 8, 20)
    print(x)
    print(L.calc_shear(x))
    print(L.calc_moment(x))

def solve_test():
    bm = Beam(10)
    bm.add_pointLoad(5,5)
    d = bm.solve()
    print(d["moment"])

if __name__ == "__main__":
    main()