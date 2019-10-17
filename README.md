# Beam-Mechanics
Python code for examining beams under load

Allows the user to add
- A beam with defined length
- Supports
    - Fixed (wall)
    - Pin and Roller combination
- Loads
    - Point Load
    - Distributed Load - either uniform, triangle, or trapezoidal
    - Point Moment


# Examples of Use

- Instantiate a beam

        beam = Beam(10)

- Add a fixed support

        beam.add_fixed_support(location = "right")

- Add loads

        beam.add_pointLoad(7, -8)
        beam.add_distLoad(0, 4, -2, -5)

- Solve, returns a dict with keys "positions", "moment", "shear" and arrays as values

        solution = beam.solve(N = 50)


# TODO
- Add plotting tools

