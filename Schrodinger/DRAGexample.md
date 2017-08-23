```@meta
DocTestSetup  = quote
    using Schrodinger
end
```

# DRAG

This example shows how to implement a NOT gate on a qubit using a Gaussian pulse. We then extend this method to a three-level slightly anharmonic energy spectrum with nearest level coupling. As we will see, the gate error increases due to leakage into the third level. To remedy this, we implement Derivative Removal by Adiabatic Gate (DRAG) which offers better gate fidelity than an ordinary Gaussian pulse \[[1]].

## NOT Gate

First, we will apply a simple NOT gate to a qubit in the ground state. The hamiltonian for our qubit in the lab frame can be written as $$ħ(ω|1⟩⟨1|+ℇ(t)σ<sub>x</sub>)$$ where $$ħω$$ is the transition energy, $$σ<sub>x</sub>$$ is the Pauli-X operator, and $$ℇ(t)=ℇ<sup>x</sup>(t)cos(ω<sub>d</sub>t)$$ represents our control of the system using a drive frequency $$ω<sub>d</sub>$$. Any control $$ℇ<sup>x</sup>(t)$$ such that the integral of $$ℇ<sup>x</sup>$$ over the total gate time equals $$π$$ will result in a complete inversion. It is common to use a Gaussian shaped π-pulse to implement a NOT gate. For our purposes, we will find it more convenient to work in the rotating frame of our drive frequency $$ω<sub>d</sub>$$. When this frequency is resonant with the qubit frequency $$ω$$, the Hamiltonian is give by $$ħℇ<sup>x</sup>(t)σ<sub>x</sub>/2$$.

```jldoctest example1
H = σx/2 # Hamiltonian
tg = 6e-9 # gate time 6ns
σ = 0.5tg # standard deviation of gaussian pulse
ψ0 = basis(2,0) # begin in ground state
tspan = (-tg/2,tg/2) # pulse centred at t=0

res1 = sesolve([qzero(2),(H,rotgaussianpulse,[tg,σ,π])],ψ0,tspan,saveat=linspace(-tg/2,tg/2,100))

plot(res1.times*1e9,levelprobs(res1.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "Excited State"]); grid();
```

## References

\[[1]] F. Motzoi, J.M. Gambetta, P. Rebentrost, and F. K. Wilhelm, "Simple pulses for elimination of leakage in weakly nonlinear qubits," Phys. Rev. Lett.

[1]: http://dx.doi.org/10.1103/PhysRevLett.103.110501
