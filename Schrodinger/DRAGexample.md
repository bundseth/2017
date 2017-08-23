```@meta
DocTestSetup  = quote
    using Schrodinger
end
```

# DRAG

This example shows how to implement a NOT gate on a qubit using a Gaussian pulse. We then extend this method to a three-level slightly anharmonic energy spectrum with nearest level coupling. As we will see, the gate error increases due to leakage into the third level. To remedy this, we implement Derivative Removal by Adiabatic Gate (DRAG) which offers better gate fidelity than an ordinary Gaussian pulse \[[1]].

## A simple NOT Gate

First, we will apply a simple NOT gate to a qubit in the ground state. The hamiltonian for our qubit in the lab frame can be written as $$ħ(ω|1⟩⟨1|+ℇ(t)σ<sub>x</sub>)$$ where $$ħω$$ is the transition energy, $$σ<sub>x</sub>$$ is the Pauli-X operator, and $$ℇ(t)=ℇ<sup>x</sup>(t)cos(ω<sub>d</sub>t)$$ represents our control of the system using a drive frequency $$ω<sub>d</sub>$$. Any control $$ℇ<sup>x</sup>(t)$$ such that the integral of $$ℇ<sup>x</sup>$$ over the total gate time equals $$π$$ will result in a complete inversion. It is common to use a Gaussian shaped π-pulse to implement a NOT gate. For our purposes, we will find it more convenient to work in the rotating frame with respect to our drive frequency $$ω<sub>d</sub>$$. When this frequency is resonant with the qubit frequency $$ω$$, the Hamiltonian is given by $$ħℇ<sup>x</sup>(t)σ<sub>x</sub>/2$$. If we are going to solve the time dynamics for this system, we also have to define our pulse in the rotating frame.
```jldoctest example1
function rotgaussianpulse(t::Real,p::Vector)
    # normalized pulse centered on t=0, begins and ends at 0
    tg  = p[1] # gate time
    σ   = p[2] # standard dev (0.5tg)
    A   = p[3] # amplitude
    B   = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ)) # normalize
    Ɛˣ = A*B*(gaussian(t,σ)-gaussian(0.5tg,σ))
    return Ɛˣ
end

H = σx/2 # Hamiltonian, using natural units where ħ=1
# output
2×2 Schrodinger.Operator{SparseMatrixCSC{Float64,Int64},1} with space dimensions 2:
 0.0  0.5
 0.5  0.0
```

We are now ready to solve for the time dynamics and plot our results. We'll use parameters that are typical for Transmon qubits.

```@setup plot1
using Schrodinger, PyPlot
function rotgaussianpulse(t::Real,p::Vector)
    # normalized pulse centered on t=0, begins and ends at 0
    tg  = p[1] # gate time
    σ   = p[2] # standard dev (0.5tg)
    A   = p[3] # amplitude
    B   = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ)) # normalize
    Ɛˣ = A*B*(gaussian(t,σ)-gaussian(0.5tg,σ))
    return Ɛˣ
end
H = σx/2 # Hamiltonian
tg = 6e-9 # gate time 6ns
σ = 0.5tg # standard deviation of gaussian pulse
g = basis(2,0) # begin in ground state
tspan = (-tg/2,tg/2) # pulse centred at t=0
res1 = sesolve([qzero(2),(H,rotgaussianpulse,[tg,σ,π])],g,tspan,saveat=linspace(-tg/2,tg/2,100))
!isdir("img") && mkdir("img")
plot(res1.times*1e9,levelprobs(res1.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "Excited State"]); grid();
savefig(joinpath("img","qubitNOT.svg"))
```
```jldoctest example1
tg = 6e-9 # gate time 6ns
σ = 0.5tg # standard deviation of gaussian pulse
g = basis(2,0) # begin in ground state
tspan = (-tg/2,tg/2) # pulse centred at t=0

res1 = sesolve([qzero(2),(H,rotgaussianpulse,[tg,σ,π])],g,tspan,saveat=linspace(-tg/2,tg/2,100))
plot(res1.times*1e9,levelprobs(res1.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "Excited State"]); grid();
```
![qubit NOT gate](img/qubitNOT.svg)

As is expected, the system moves from the ground state to the excited state.

## Three-level system

In reality, one must worry about leakage into other states of the system, especially when dealing with short gate times. Let's see how our Gaussian pulse performs on a three-level anharmonic system. We will again work in the rotating frame at resonance with the qubit frequency. This system will have an anharmonicity $$Δ$$, which is the detuning of the 2nd excited state with respect to the drive frequency, and an additional parameter $$λ$$ describing the relative strength of the 1-2 transition compared to the 0-1 transition (see \[[1]] for more details). Let's create the Hamiltonian and perform the same time evolution.

```@setup plot2
using Schrodinger, PyPlot
function rotgaussianpulse(t::Real,p::Vector)
    # normalized pulse centered on t=0, begins and ends at 0
    tg  = p[1] # gate time
    σ   = p[2] # standard dev (0.5tg)
    A   = p[3] # amplitude
    B   = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ)) # normalize
    Ɛˣ = A*B*(gaussian(t,σ)-gaussian(0.5tg,σ))
    return Ɛˣ
end
Δ = 2π*(-400e6) # anharmonicity
λ = √2 # relative transition strength
Π₂ = basis(3,2) * basis(3,2)' # projector for the 2nd level
Hc = Δ*Π₂ # constant Hamiltonian
Hd = create(3)/2 + destroy(3)/2 # affected by Ɛˣ
g = basis(3,0)
res2 = sesolve([Hc,(Hd,rotgaussianpulse,[tg,σ,π])],g,tspan)
!isdir("img") && mkdir("img")
plot(res2.times*1e9,levelprobs(res2.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "1st Excited State", "2nd Excited State"]); grid()
savefig(joinpath("img","3levelNOT.svg"))
```
```jldoctest example1
Δ = 2π*(-400e6) # anharmonicity
λ = √2 # relative transition strength
Π₂ = basis(3,2) * basis(3,2)' # projector for the 2nd level
Hc = Δ*Π₂ # constant Hamiltonian
Hd = create(3)/2 + destroy(3)/2 # affected by Ɛˣ
g = basis(3,0)

res2 = sesolve([Hc,(Hd,rotgaussianpulse,[tg,σ,π])],g,tspan)
figure()
plot(res2.times*1e9,levelprobs(res2.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "1st Excited State", "2nd Excited State"]); grid()
savefig(joinpath("img","3levelNOT.svg"))
```
![3-level NOT gate](img/3levelNOT.svg)

Instead of working perfectly, our system leaks into the 2nd energy level (we'll quantify this [later](#fidelity)). This becomes problematic once we're trying 

## Fidelity

TODO

## References

\[[1]] F. Motzoi, J.M. Gambetta, P. Rebentrost, and F. K. Wilhelm, "Simple pulses for elimination of leakage in weakly nonlinear qubits," Phys. Rev. Lett.

[1]: http://dx.doi.org/10.1103/PhysRevLett.103.110501
