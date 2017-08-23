```@meta
DocTestSetup  = quote
    using Schrodinger, PyPlot
end
```

# DRAG

This example shows how to implement a NOT gate on a qubit using a Gaussian pulse. We then extend this method to a three-level slightly anharmonic energy spectrum with nearest level coupling. As we will see, the gate error increases due to leakage into the third level. To remedy this, we implement Derivative Removal by Adiabatic Gate (DRAG) which offers better gate fidelity than an ordinary Gaussian pulse \[[1]].

## NOT Gate

First, we will apply a simple NOT gate to a qubit in the ground state. The hamiltonian...

## References

\[[1]] F. Motzoi, J.M. Gambetta, P. Rebentrost, and F. K. Wilhelm, "Simple pulses for elimination of leakage in weakly nonlinear qubits," Phys. Rev. Lett.

[1]: http://dx.doi.org/10.1103/PhysRevLett.103.110501
