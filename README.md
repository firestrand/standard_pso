# Standard PSO (SPSO) Language Implementations

This repository contains multiple language implementations of the **Standard Particle Swarm Optimization (SPSO)** algorithm. It serves both as an archive of the original SPSO2011 and as a collaborative development ground for the upcoming SPSO2025 standard.

## Overview

The Standard PSO (SPSO) algorithm is a foundational variant of Particle Swarm Optimization, providing a clear baseline for research, benchmarking, and algorithmic comparisons. This repository preserves implementations across several programming languages while working actively toward defining the next generation standard, **SPSO2025**.

### Available Implementations

- **C Version (Original)**: Direct derivative of the SPSO2011 reference implementation by Maurice Clerc and collaborators. It maintains the integrity and original specifications of the 2011 standard.

- **C++ Version (Port)**: Modernized C++ port of SPSO2011 featuring object-oriented design and enhanced maintainability, faithfully preserving algorithmic behavior.

- **Python Version (Reimplementation)**: Pythonic recreation of SPSO2011 emphasizing readability, ease of use, and suitability for research environments.

- **Rust Version (Coming Soon)**: Planned future implementation prioritizing safety, efficiency, and performance.

## Repository Goals

1. **Preservation**: Document and maintain historical accuracy of the original SPSO2011 implementation.
2. **Consistency**: Offer reliable, consistent implementations across multiple programming languages.
3. **Innovation**: Develop the SPSO2025 standard, incorporating insights and best practices learned over more than a decade of research.
4. **Research Facilitation**: Provide robust benchmarks, clearly defined parameters, and reproducible results to support comparative PSO research.

## Historical Context

The SPSO2011 standard was established in response to the need for a stable, universally accepted baseline PSO variant. Key contributions of SPSO2011 include:

- **Coordinate Independence (Rotation Invariance)**: Particle updates are independent of coordinate systems, improving fairness across benchmarks.
- **Clear Particle Update Rules**: Well-defined particle behavior and parameters ensure reproducibility.
- **Standard Parameterization**: Defaults provide a reliable baseline for comparisons.
- **Benchmark Functions**: Comprehensive standard benchmarks enable meaningful comparative studies.

Detailed historical notes and algorithmic specifics are documented within the [original C implementation directory](c-original/ReadMe.txt).

## SPSO2025 Development

The upcoming SPSO2025 standard aims to:

- Modernize the SPSO algorithm while preserving fundamental principles.
- Integrate valuable research insights gained since 2011.
- Define clear versioning and ensure backward compatibility where possible.
- Establish enhanced and more rigorous benchmark procedures.

Contributions and discussions are encouraged as we shape this new standard.

## Original Contributors (SPSO2011)

The original SPSO2011 algorithm benefited from significant input by numerous researchers:

- Auger, Anne
- Blackwell, Tim
- Bratton, Dan
- Clerc, Maurice
- Croussette, Sylvain
- Dattasharma, Abhi
- Eberhart, Russel
- Hansen, Nikolaus
- Helwig, Sabine
- Keko, Hrvoje
- Kennedy, James
- Krohling, Renato
- Langdon, William
- Li, Wentao
- Liu, Hongbo
- Miranda, Vladimiro
- Poli, Riccardo
- Serra, Pablo
- Silvers, Travis
- Spears, William
- Stickel, Manfred
- Yue, Shuai

## Current Maintainer

**Travis Silvers (2024–present)**

## Building and Running

Each language implementation has its own build and run instructions within respective subdirectories. The original C version utilizes CMake for ease of compilation:

```bash
./build.sh
./run_tests.sh
```

## License Information

- **Original SPSO2011 C Implementation**: Provided under Academic/Research License.
- **New Implementations and Modifications**: Released under the MIT License.
- **Mersenne Twister**: BSD License (Copyright 2004, Makoto Matsumoto and Takuji Nishimura).
- **Wyhash**: The Unlicense (Author: Wang Yi).

For specific licensing queries, contact the current maintainer.

## Contributing

Contributions are welcome in several areas:

- **Bug fixes and enhancements** to existing implementations.
- **Documentation improvements** and example expansions.
- **New language implementations**.
- **Ideas and proposals** for the SPSO2025 specification.

Please direct contributions via GitHub pull requests or contact the current maintainer directly.

## Contact

For questions, suggestions, or contributions, please reach out to the current maintainer.

## References

For a detailed understanding of SPSO2011's motivation, algorithm specifics, equations, and historical development, consult the documentation within the [original C implementation directory](c-version/ReadMe.txt).

