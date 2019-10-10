# Beam-Mechanics
Python code for examining beams under load

Allows the user to add
- A beam with defined length
- Supports
    - Fixed (wall)
    - Pin
    - Roller
    - Note: only allows one fixed or 2 other (fx calculations are ignored)
- Loads
    - Point Load
    - Distributed Load - either uniform, triangle, or trapezoidal
    - Point Moment


TODO
- checks for adding loads
    - ensure they are within beam span
    - make sure distload start and endforce are of same sign
- update support check to allow for only either fixed or pin/roller combo
- calculation of support reactions
- solve for shear/moment diagrams
