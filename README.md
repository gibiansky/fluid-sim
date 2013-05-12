Fluid Simulations
=========

This is a [Smoothed Particle Hydrodynamics](http://en.wikipedia.org/wiki/Smoothed_particle_hydrodynamics) fluid simulation currently implemented in [Go](http://golang.org).

In order to run and compile it, execute

    ./run

Parts of the simulation:
  * OpenGL code for the simulator environment is in the 'simulator' package (src/simulator)
  * Vector operations are in the 'vector' package (src/vector)
  * Main simulation is 'main.go', with supporting files

Requirements:
  * Go programming language installed (with command-line tools, of course)
  * [OpenGL for Go](https://github.com/go-gl): gltext, gl, glfw, glu (possibly others)
