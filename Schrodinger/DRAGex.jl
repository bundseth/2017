using Schrodinger, PyPlot

function rotgaussianpulse(t::Real,p::Vector)
    # normalized pulse centered on t=0, begins and ends at 0
    σ   = p[1] # standard dev (0.5tg)
    tg  = p[2] # gate time (6e-9)
    ω10 = p[3] # 0 - 1 transition angular freq (0)
    ϕ   = p[4] # phase (0)
    A   = p[5] # amplitude
    B   = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ)) # normalize with: integral from -g/2 to g/2 (exp(-0.5(t/s)^2) - exp(-0.5(0.5g/s)^2))dt
    Ɛˣ = A*B*(gaussian(t,σ)-gaussian(0.5tg,σ)) # normalized pulse that is 0 at t=-tg/2,tg/2
    return Ɛˣ # assumes resonance and no phase
end

# create 2-level Hamiltonian

# lab frame
# N = 2 # number of levels
# Π₁ = basis(N,1) * basis(N,1)'

# rotating frame
N = 2
δ₁ = 0
Hc = qzero(N)
Hdx = create(N)/2 + destroy(N)/2
Hdy = im*create(N)/2 - im*destroy(N)/2

# time-evo setup
ψ0 = basis(N,0)
tg = 6e-9
σ = 0.5tg
tspan = (-tg/2,tg/2)

res1 = sesolve([Hc,(Hdx,rotgaussianpulse,[σ,tg,0,0,π])],ψ0,tspan,saveat=linspace(-tg/2,tg/2,101))
plot(res1.times*1e9,levelprobs(res1.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "Excited State"]); grid();
