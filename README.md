# MyGeoDock

MyGeoDock is a Python prototype for geometric molecular expansion in molecular docking. This project is inspired by the paper ["Tunable Approximations to Control Time-to-Solution in
an HPC Molecular Docking Mini-App"](https://arxiv.org/pdf/1901.06363.pdf).

## Origin and Publication

This code serves as the basis for the research published in:

["Quantum Molecular Unfolding"](https://iopscience.iop.org/article/10.1088/2058-9565/ac73af/meta)

Please note that the original code of the publication is now the intellectual property of Dompe Pharmaceutics.

## Features

### Software Knobs for Tunable Approximations

MyGeoDock introduces "software knobs," which are adjustable parameters in the code that control the trade-off between speed and accuracy. These knobs can be fine-tuned for:

- Precision of calculations
- Number of iterations
- Granularity of the search space in molecular docking

These knobs are implemented as configurable parameters that can be adjusted before or during runtime.

### Mini-App

MyGeoDock includes a "Mini-App," which is a simplified version of a full-scale molecular docking application. This Mini-App is used for:

- Testing and validating the approach
- Easy experimentation and benchmarking

The Mini-App is designed with modularity in mind, allowing for easy tweaking of knobs and integration with the autotuning system.

### High-Performance Computing (HPC) Integration

MyGeoDock is optimized for High-Performance Computing systems, featuring:

- Parallel execution optimization
- Efficient memory usage
- Scalability across multiple processors

The code incorporates:

- Parallelization techniques (e.g., using MPI or OpenMP)
- Optimized data structures
- Efficient I/O operations
- MPI4PY, Multiprocessing
