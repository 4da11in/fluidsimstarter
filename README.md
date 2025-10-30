# MPM stuff
I told my professor for a project class that I wanted to implement FLIP fluid simulation. He suggested it would be even more impressive if I implemented that, then implemented MPM (Material Point Method), which is often used to simulate snow, mud, and other breakable things. I don't think I knew what I was getting myself into when I agreed, but here we are (and it works!)
My advisor and I originally reasoned that because MPM is a generalization of flip, the code/logic might be similar. As I delved into the paper by Disney (and an MPM course published at SIGGRAPH) I realized how mistaken I was. Unfortunately I had already started working on it as a branch of my fluid sim project, so I decided to stick with it. As a result, some of these files aren't actually necessary for the simulation. Sorry about that. :)
This was an awesome opportunity to learn more about how MPM works (I was pleasantly surprised when I had to break out linear algebra again). Special thanks to Victoria Kala who was generously willing to help me when I had questions about my code, and to Seth Holladay for encouraging me to take on this project.

# fluidsimstarter (original readme for fluid sim starter code)

Fluidsimstarter is a simple 2D marker-and-cell (MAC) fluid simulator intended for those wanting to jump into fluid simulation implementation. I have intentionally removed the sections of code that perform the pressure solve in order to give the programmer a chance to code them up themselves. This implementation is based on the fluid flow for the rest of us paper. It is not necessarily bug free or very efficient, so feel free to improve it as needed.

To use it, you will need to change the paths found in Simulator.cpp to the directories that you would like to serialize your particles and grids out to. To visualize your grids and particles, I have provided sim_visualizer.hipnc. Again, you will need to change the paths in the Houdini scene to your serialized data.

