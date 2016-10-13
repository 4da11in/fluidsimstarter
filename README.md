# fluidsimstarter

Fluidsimstarter is a simple 2D marker-and-cell (MAC) fluid simulator intended for those wanting to jump into fluid simulation implementation. I have intentionally removed the sections of code that perform the pressure solve in order to give the programmer a chance to code them up themselves. This implementation is based on the fluid flow for the rest of us paper. It is not necessarily bug free or very efficient, so feel free to improve it as needed.

To use it, you will need to change the paths found in Simulator.cpp to the directories that you would like to serialize your particles and grids out to. To visualize your grids and particles, I have provided sim_visualizer.hipnc. Again, you will need to change the paths in the Houdini scene to your serialized data.

