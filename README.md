# Noise induced excitability: bloom bust and extirpation in autotoxic population dynamics

## Code structure
* **M_declarations.f90**: Define public variables $x$, $y$, $t$, $D$, $r$. $x$, $y$, and $t$ are used in functions and subroutines to characterize the state of the system. $D$ and $r$ represent, respectively, the parameters $D$ and $\rho$ in our text (see License for reference).
* **M_subroutines.f90**: Subroutines like the Milstein method to sample realizations of the process.
* **M_functions.f90**: Mathematical functions used in the main code.
* **dranxor.f90**: Pseudo random number generator.
* **main_Draw_trajectories.f90**: example of code to generate and save $x$, $y$, and $\theta_W$ for different values of $t$.

## License
This project is shared for **academic and research purposes**. It is free to use, redistribute, modify, and share for research purposes, provided that proper credit is given to the authors through citation of: 

Moreno-Spiegelberg, P., & Aguilar, J. (2026). arXiv:2601.20670.

*Noise-induced excitability: bloom, bust and extirpation in autotoxic population dynamics.*

