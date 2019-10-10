class Beam:
    def __init__(self, span, units = "m"):
        self.span = span
        self.units = units
        self.point_loads = []
        self.dist_loads = []
        self.point_moments = []
        self.supports = []
    
    def __str__(self):
        s = f"Beam with length {self.span}{self.units}"
        #s += f"\n\t Supports: "
        #s += "None" if len(self.supports) == 0 else str([str(support) for support in self.supports])
        return s
    
    def add_pointLoad(self, pos, force):
        self.point_loads.append(PointLoad(pos, force))

    def add_distLoad(self, start, end, force):
        self.dist_loads.append(DistLoad(start, end, force))
    
    def add_moment(self, pos, moment):
        self.point_moments.append(PointMoment(pos, moment))

    def add_support(self, type, pos):
        if type == "roller":
            self.supports.append(RollerSupport(pos))
        elif type == "pin":
            self.supports.append(PinSupport(pos))
        elif type == "fixed":
            self.supports.append(FixedSupport(pos))
    


def is_indeterminate(beam):
    unknowns = sum([support.unknowns for support in beam.supports])
    return True if unknowns > 3 else False

class PointLoad:
    def __init__(self, pos, force):
        self.pos = pos
        self.force = force

class DistLoad:
    def __init__(self, start, end, startForce, endForce):
        self.start = start
        self.end = end
        self.startForce = startForce
        self.endForce = endForce

class PointMoment:
    def __init__(self, pos, moment):
        self.pos = pos
        self.moment = moment

class Support:
    def __init__(self, style, pos, unknowns):
        self.pos = pos
        self.unknowns = unknowns
        self.style = style
    
    def __str__(self):
        return f"{self.style} support @ {self.pos}"

class FixedSupport(Support):
    def __init__(self, pos):
        super().__init__(style = "Fixed", pos = pos, unknowns = 3)

class RollerSupport(Support):
    def __init__(self, pos):
        super().__init__(style = "Roller", pos = pos, unknowns = 1)

class PinSupport(Support):
    def __init__(self, pos):
        super().__init__(style = "Pin", pos = pos, unknowns = 2)









def main():
    bm = Beam(10)
    bm.add_support("fixed", 0)
    bm.add_support("roller", 5)
    print(bm)
    for support in bm.supports:
        print(support)
    

if __name__ == "__main__":
    main()