# relativitysim

This program plots Lorentz transforms in 2+1 dimentions with an animated transition between rest frames.

Press 1, 2, 3 or 4 to move to the rest frame of each of 4 spaceships. Press q to quit.
(Don't press a number until the animation stops moving or program will freeze)

To add more spaceships or locations, just insert a new line right after the others are defined with the same format:
variableName = Object([X, Y, T], [vel X, vel Y], "color", "name for label on plot")
X, Y, T are coordinates as viewed by a stationary observer at the origin.
vel X and vel Y are relative to c=1, so magnitude of [vel X, vel Y] must be less than 1.
To make a new object have a button for its animation, copy the format near the end of the program.

