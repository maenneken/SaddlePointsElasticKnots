# Minimum Energy Paths between Multistable Elastic Knots

## Overview

This project investigates transition pathways between distinct equilibrium configurations of multistable elastic knots.

Building upon the ElasticKnots framework, the approach combines sampling-based robotic motion planning with energy-based path optimization to compute Minimum Energy Paths (MEPs) between configurations of the same knot topology.

A central challenge is that conventional linear interpolation between configurations may introduce self-intersections and violate the topology of the knot.

To address this, sampling-based tree planners are used to construct feasible initial paths, which are subsequently optimized using the Nudged Elastic Band (NEB) method.

## Research Contributions

- Implemented RRT-EST-based motion planning in C++.
- Developed a topology-preserving approach for connecting distinct equilibrium configurations of elastic knots.
- Computed initial transition paths using sampling-based motion planning.
- Optimized transition pathways using the Nudged Elastic Band (NEB) method.
- Computed and visualized Minimum Energy Paths between multistable equilibrium states.
- Evaluated the resulting paths using quantitative energy and planning metrics.

## Methodology

### 1. Configuration Space Exploration

Sampling-based motion planning is employed to explore the configuration space and identify feasible transitions between equilibrium configurations.

### 2. Initial Path Construction

RRT-EST is utilized to construct an initial path while avoiding self-intersections that would violate the knot topology.

### 3. Minimum Energy Path Optimization

The initial path is refined using the Nudged Elastic Band method to obtain an energy-minimizing transition pathway.

### 4. Transition State Identification

The maximum-energy configuration along the optimized path is used to identify a candidate transition state between two equilibrium configurations.

## Results

The computed transition pathways provide insight into the energy landscape and the transitions between multistable equilibrium states of elastic knots.

[Insert energy profile along the MEP]

[Insert visualization of the RRT-EST exploration]

[Insert animation of the optimized transition pathway]

## Technical Skills

- C++
- Sampling-based Motion Planning
- RRT-EST
- Nudged Elastic Band (NEB)
- Minimum Energy Path Computation
- Computational Mechanics
- Elastic Rod Simulation
- Algorithm Design and Quantitative Evaluation



## Acknowledgements and Framework

This project is based on the **ElasticKnots** framework, developed for the simulation of knotted elastic rods with self-contacts.

The framework is associated with the SIGGRAPH 2023 paper:

*Computational Exploration of Multistable Elastic Knots.*

The original framework was used as the basis for the motion planning and Minimum Energy Path computations presented in this project.

- **Original repository:** [EPFL-LGG/ElasticKnots](https://github.com/EPFL-LGG/ElasticKnots)
- **Research paper:** [Computational Exploration of Multistable Elastic Knots](https://infoscience.epfl.ch/record/302408)
- **Project webpage:** [Elastic Knots – EPFL](https://www.epfl.ch/labs/gcm/research-projects/elastic-knots/)

The original ElasticKnots framework and its underlying simulation capabilities are credited to the respective authors. The motion planning implementation and Minimum Energy Path computations described in this project were developed as part of my Master's thesis at TU Berlin.

# Getting Started

## C++ Code Dependencies
To install the required dependencies, please follow the instructions from [https://github.com/jpanetta/ElasticRods](https://github.com/EPFL-LGG/ElasticKnots)

## Obtaining the Code
Clone this repository *recursively* so that its submodules are also downloaded:

```
git clone git@github.com:maenneken/SaddlePointsElasticKnots.git --recurse-submodules
```
## Note from Simon
you might need to delete pybind11 and clone the newest yourself and add to attr.h "#include <cstdint>". (path: 3rdparty/ElasticKnots/3rdparty/ElasticRods/3rdparty/MeshFEM/3rdparty/pybind11/include/pybind11)

## Building
Build the C++ code and its Python bindings:

```
python setup.py
```


# Debug for me
Build with setup_debug.py

perf logging:

/usr/lib/linux-tools/6.8.0-107-generic/perf record ./my_program
/usr/lib/linux-tools/6.8.0-107-generic/perf report
