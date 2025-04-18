# CSCI520 Assignment 1: Simulating a Jello Cube

**Assignment Link:** [Simulating a Jello Cube](https://viterbi-web.usc.edu/~jbarbic/cs520-s25/assign1/)  
**Demo Video:** [Youtube Link](https://youtu.be/eozCOgNBC7s)

---

## Introduction

This project simulates a jello cube using a 3D mass-spring system. The simulation confines the cube within a 4×4×4 bounding box. The cube is modeled as an 8×8×8 grid of mass points, connected by structural, shear, and bend springs. The system includes collision detection with the bounding box, an inclined plane, and a sphere. Additional forces such as gravity and an external force field are applied, and the simulation can be run using Euler, RK4, or an implicit integrator.

---

## How to Run

In the project directory, run:

```bash
make
./jello <worldfile>
```
For example:

```bash
./jello world/gravity.w
./jello world/jello.w
./jello world/moveLeft.w
./jello world/rotate.w
./jello world/skewedCorner.w
```

## Keyboard Commands

- **ESC:** Exit the application
- **Camera Control:**
  - **e:** Reset the camera to the original position
  - **z:** Zoom in
  - **x:** Zoom out
- **View Mode:**
  - **h:** Toggle shear springs
  - **b:** Toggle bend springs
  - **s:** Toggle structural springs
  - **v:** Switch viewing mode
- **Jello Control:**
  - **p:** Pause/unpause the simulation
- **Screenshot:**
  - **Space:** Save a screenshot (saved in the `screenshots` folder)
- **Extra Features:**
  - **g:** Toggle gravity on the jello
  - **w:** Toggle a steady collision sphere at the top corner (enables collision detection with the sphere)

## Implementation Details

### Basic Simulation

- **Mass-Spring Model:**  
  The jello cube is represented as a grid of 512 mass points connected by structural, shear, and bend springs. The simulation solves Newton’s second law using Euler or RK4.

- **Collision Detection:**  
  Collisions are detected with the bounding box, using penalty method to generate collision forces with different elasticity and damping parameters to the springs.

- **Force Field:**  
  An external, non-homogeneous force field (if specified in the world file) is provided as a 3D grid. Trilinear interpolation is used to compute the force at an arbitrary position, and the resulting force is added to the net force on each mass point.

### Visual Enhancements

- **Background and Lighting:**  
  The background color and lighting parameters have been modified for improved visual quality.

- **Texture Mapping:**  
  A 1280×1280 texture image (`texture/texture1.jpg`) is applied to all six faces of the jello cube when view mode is on, used stb_image for texture loading.

- **Semi-Transparent Inclined Plane:**  
  When present (as indicated by the world file), the inclined plane is rendered semi-transparently.

- **On-Screen Display:**  
  The simulation displays:
  - Frames Per Second (FPS) in the top-left corner.
  - The last pressed keyboard key in the bottom-left corner.
  - The name of the currently used integrator in the bottom-right corner.
  - 3D coordinate axes with clearly defined colors.

### Small Features

- **Gravity Control:**  
  Press 'g' to toggle a constant gravitational force on the jello.

- **Inclined Plane:**  
  The simulation supports an inclined plane (specified via the `a, b, c, d` fields) for collision detection and response.

- **Collision Sphere:**  
  Press 'w' to display a steady sphere at the top corner. The sphere is integrated into the collision detection system.

- **Pause Functionality:**  
  The pause feature (triggered by 'p') is implemented to stop the simulation temporarily (the assignment instruction mentioned this functionality but did not implement it in the starter code).

### Large Features

- **Multiple Integrators:**  
  The simulation supports Euler, RK4, and IPC. You can choose which integrator to use by specify it at the world file, here are some custom world file to try for different integrator:
  ```bash
  ./jello world/custom-IPC.w
  ./jello world/custom-RK4.w
  ./jello world/custom-Euler.w
  ```

### IPC Details:
My IPC implementation incorporates the **inertia term**, **mass‐spring potential energy**, **gravity energy**, and **barrier energy** with full gradients and Hessians, along with **filtered line search** and **continuous collision detection**; however, it requires relatively high elasticity coefficients, lacks DOF elimination and self-collision detection, and does not support initial interpenetration between objects.

The implicit integrator minimizes the energy:
```math
E(x) = \frac{1}{2}||x-(x^n+hv^n)||^2_M+h^2P(x)
```
with the following components:

- Inertia Term
```math
E_I(x) = \frac{1}{2}||x-\tilde{x}^n||^2_M
```
```math
\nabla E_I(x) = M(x-\tilde{x}^n)
```
```math
\nabla^2 E_I(x) = M
```
- Mass-Spring Potential Energy
```math
P_e(x) = l^2\frac{1}{2}k(\frac{||x_1-x_2||^2}{l^2}-1)^2
```
```math
\frac{\partial P_e}{\partial x_1}(x) 
= - \frac{\partial P_e}{\partial x_2}(x)
= 2k(\frac{||x_1-x_2||^2}{l^2}-1)(x_1-x_2)
```
```math
\frac{\partial^2 P_e}{\partial x_1^2}(x) 
= \frac{\partial^2 P_e}{\partial x_2^2}(x)
= -\frac{\partial^2 P_e}{\partial x_1x_2}(x)
= -\frac{\partial^2 P_e}{\partial x_2x_1}(x) \\
= \frac{4k}{l^2}(x_1 - x_2)(x_1-x_2)^T+2k(\frac{||x_1-x_2||^2}{l^2}-1)I \\
= \frac{2k}{l^2}(2(x_1-x_2)(x_1-x_2)^T+(||x_1-x_2||^2-l^2)I)
```
- Gravity Energy
```math
P(x)=-x^TMg
```
```math
\nabla P(x)=-Mg
```
```math
\nabla^2 P(x)=0
```
- Barrier Energy
```math
d(x)=x_y-y_0
```
```math
\nabla d(x)=
\begin{bmatrix}
0 \\
1
\end{bmatrix}
```
```math
\nabla^2d(x)=0
```

## Sample Images
![sample 1](<sample-images/sample1.jpg>)
![sample 2](<sample-images/sample2.jpg>)
![sample 3](<sample-images/sample3.jpg>)