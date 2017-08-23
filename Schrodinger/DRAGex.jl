using Schrodinger, PyPlot

function rotgaussianpulse(t::Real,p::Vector)
    # normalized pulse centered on t=0, begins and ends at 0
    tg  = p[1] # gate time
    σ   = p[2] # standard dev (0.5tg)
    A   = p[3] # amplitude
    B   = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ))
    Ɛˣ = A*B*(gaussian(t,σ)-gaussian(0.5tg,σ))
    return Ɛˣ
end

function DRAGx(t::Real,p::Vector)
    # Ɛˣ for fifth order DRAG
    tg  = p[1] # gate time
    σ   = p[2] # standard deviation
    A   = p[3] # amplitude
    Δ   = p[4] # anharmonicity
    λ   = p[5] # relative strength of transitions
    Ɛˣ = rotgaussianpulse(t,[tg,σ,A]) + (λ^2-4)*rotgaussianpulse(t,[tg,σ,A])^3/(8*Δ^2) - (13λ^4-76λ^2+112)*rotgaussianpulse(t,[tg,σ,A])^5/(128Δ^4)
    return Ɛˣ
end

function DRAGy(t::Real,p::Vector)
    # Ɛʸ for fifth order DRAG
    tg  = p[1] # gate time
    σ   = p[2] # standard deviation
    A   = p[3] # amplitude
    Δ   = p[4] # anharmonicity
    λ   = p[5] # relative strength of transitions
    B   = inv(√(2π)*σ*erf(tg/(√(8)*σ))-tg*gaussian(0.5tg,σ))
    Ɛˣ′ = A*B*Schrodinger.gaussianprime(t,σ)
    Ɛʸ = -Ɛˣ′/Δ + 33*(λ^2-2)*rotgaussianpulse(t,[tg,σ,A])^2*Ɛˣ′/(24*Δ^3)
    return Ɛʸ
end

function detuning(t::Real,p::Vector)
    # dynamical detuning for fifth order DRAG
    tg  = p[1] # gate time
    σ   = p[2] # standard deviation
    A   = p[3] # amplitude
    Δ   = p[4] # anharmonicity
    λ   = p[5] # relative strength of transitions
    δ₁ = (λ^2-4)*rotgaussianpulse(t,[tg,σ,A])^2/(4*Δ) - (λ^4-7λ^2+12)*rotgaussianpulse(t,[tg,σ,A])^4/(16Δ^3)
    return δ₁
end

H = σx/2 # same as create(2)/2 + destroy(2)/2

g = basis(2,0)
tg = 6e-9
σ = 0.5tg
tspan = (-tg/2,tg/2)

res1 = sesolve([qzero(2),(H,rotgaussianpulse,[tg,σ,π])],g,tspan,saveat=linspace(-tg/2,tg/2,100))
!isdir("img") && mkdir("img")
plot(res1.times*1e9,levelprobs(res1.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "Excited State"]); grid()
savefig(joinpath("img","qubitNOT.svg"))
clf()

Δ = 2π*(-400e6) # anharmonicity
λ = √2 # relative transition strength
Π₂ = basis(3,2) * basis(3,2)' # projector for the 2nd level
Hc = Δ*Π₂ # constant Hamiltonian
Hd = create(3)/2 + destroy(3)/2 # affected by Ɛˣ
g = basis(3,0)

res2 = sesolve([Hc,(Hd,rotgaussianpulse,[tg,σ,π])],g,tspan)
plot(res2.times*1e9,levelprobs(res2.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "1st Excited State", "2nd Excited State"]); grid()
savefig(joinpath("img","3levelNOT.svg"))
clf()

Π₁ = basis(3,1) * basis(3,1)' # projector for 1st level
Hdet = Π₁ # affected by dynamical detuning
Hdx = create(3)/2 + destroy(3)/2 # affected by Ɛˣ
Hdy = im*create(3)/2 - im*destroy(3)/2 # affected by Ɛʸ

res3 = sesolve([Hc,(Hdet,detuning,[tg,σ,π,Δ,λ]),(Hdx,DRAGx,[tg,σ,π,Δ,λ]),(Hdy,DRAGy,[tg,σ,π,Δ,λ])],g,tspan)
plot(res3.times*1e9,levelprobs(res3.states)); xlabel("Time (ns)"); ylabel("Level Probabilities"); legend(["Ground State", "1st Excited State", "2nd Excited State"]); grid()
savefig(joinpath("img","3levelDRAG.svg"))
