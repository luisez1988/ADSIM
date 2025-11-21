# Matrices and vectors in ADSIM

This is documentation to create the matrices and vectors variables for ADSIM. 

## Lumped mass system

The formulation of ADSIM is based on fully explicit lumped mass system formulation. Thus the lumped mass matrix can be stored in variable of the form `M_L(Nnodes)`. This mass matrix is the same for all gases.

## Concentration flows

There are five concentration flows vectors that need to be calculated, namely:
1. Advection flow
2. Gravitational flow
3. Diffusion flow
4. Boundary flows
5. Sink/source flows

Each flow is particular for each gas species. Thus the flows should be specified as `q(Nnodes, Ngases)`. The only flow that is only defined for one gas is the sink/source term. Which will be only activated for a unique gas at this moment, and that is CO2.

## Initialization of boundary flows

For boundary conditions where flow is assigned, these must be initialized from data written in the GID interface.

Note that the GID interfae only permits the assignment of unirform steady normal flows. Thus, the parameter assigned to the nodes, must consider the length of influence associated with that node as in the example below.

.------ --------.--------
 h/2       h/2     h/2

## Finding flows in pressure nodes

In the future, it would be important to determine the flow that passes trought areas defined by the user.

# Writting vtk output

ADSIM uses vtk templates to write the output for visualization. 

The output consists of mesh data with nodal and element information. 

At this stage, only the advection and diffusion flows could be measured at the Gauss points, but this will be implemented later.

Currently all data are mapped to the nodes including:
1. Concentrations of all gases (scalar)
2. Total concentration (scalar)
3. Absolute pressure (scalar)
4. Rate of concentrations (scalar)
5. Rate of reactions (scalar)
6. Lime concentration (scalar)
7. CO2 concentration (scalar)
8. CaCO3 concentration (scalar)
9. Degree of carbonation ($DoC$)
10. Volumetric binder content (scalar)
11. Gas seepage velocity (vector)
12. Temperature (scalar)
13. Rate of temperature (scalar)

The output files should be correctly serialized so that Paraview can understand the 