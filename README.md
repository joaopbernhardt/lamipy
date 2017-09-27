# lamipy
Laminated composites engineering simulations in Python.

<img src="https://upload.wikimedia.org/wikipedia/commons/1/13/Composite_3d.png" data-canonical-src="https://upload.wikimedia.org/wikipedia/commons/1/13/Composite_3d.png" width="300" height="200" />

## Brief description:

This project's purpose is to provide simple & accurate computations for engineering simulations of **laminated composites**.

The implementation uses the **Classical Laminate Theory (CLT)** for computations. The summary of this theory can be found in: [Prof. W. Stein's document](http://wstein.org/edu/2010/480b/projects/05-lamination_theory/A%20summary%20of%20Classical%20Lamination%20Theory.pdf).

Currently lamipy is *not ready for use* due to the lack of verification and validation tests, which will be produced soon.

## Project goals:

- Provide a simple interface for laminate testing;
- Display accurate information of results, including graphical analysis;

## Technical features:
- [x] Calculates CLT stresses & strains (lamina & laminate coord. systems);
- [x] Failure tests for individual laminas (Tsai-Wu, Max. Stress, Max Strain and Hashin criteria);
- [x] Progressive failure analysis using *Ply Discount* method;
- [ ] Plotting of progressive failure analysis;
- [ ] Easy entry of laminate data;
- [ ] Puck failure criteria;
- [ ] Thermal & moisture effects on CLT calculations;
- [ ] Laminate materials simple database
- [ ] ...

## How to use:

Currently, *lamipy* is **not ready for general usage**. If you want to test the code, here are the steps:

*Note: using lamipy requires a python environment for executing the Run_FailureTest.py file. It is also required **numpy** and **matplotlib** libraries (both of which are easily encountered in many python distributions).*
1. Clone this repository;
1. Edit Run_FailureTest.py to input the laminate data;
1. Execute.

## Understanding lamipy:

*To be written.*
