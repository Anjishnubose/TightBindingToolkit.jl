module SpecResponse

using LinearAlgebra

using ..TightBindingToolkit.BZone: BZ, GetQIndex
using ..TightBindingToolkit.Hams: Hamiltonian
using ..TightBindingToolkit.Useful: DistFunction

export spectral, response



function spectral(omega::Float64, eta::Float64, ham::Hamiltonian)

    G = Ref((omega+im*eta)*I) .- ham.H
    G = inv.(G)

    return (G - adjoint.(G)) / (2*im)
end

function spectral(eta::Float64, ham::Hamiltonian ; n_w::Int = 101)

    omegas = collect(range(ham.bandwidth..., n_w))
    return spectral.(omegas, eta, Ref(ham))
end


function response(omega::Float64, eta::Float64, QIndex::Vector{Int64}, Aws, bandwidth::Tuple{Float64, Float64} ;
    n_w::Int = 101, mu::Float64 = 0.0, T::Float64 = 0.001)

    es = collect(range(bandwidth..., n_w))
    response = zeros(Complex{Float64}, size(Aws[1][1, 1]) .^ 2)
    for (i1, e1) in enumerate(es)
        for (i2, e2) in enumerate(es)
            Ak1 = Aws[i1]
            AQk2 = circshift(Aws[i2], -QIndex)
            Aksum = sum(kron.(AQk2, Ak1))
            nfs = (DistFunction(e1 ; T=T, mu=mu, stat=-1) - DistFunction(e2 ; T=T, mu=mu, stat=1))/(omega +im*eta + e1 - e2)
            response += nfs * Aksum

        end
    end
println("response done")
    return response
end


function response(omegas::Vector{Float64}, eta::Float64, Qs::Vector{Vector{Float64}}, ham::Hamiltonian, bz::BZ ;
    n_w::Int = 101, mu::Float64 = 0.0, T::Float64 = 0.001)

    QIndices = GetQIndex.(Qs, Ref(bz)) .- Ref(GetQIndex(zeros(length(bz.basis)), bz))
    Aws = spectral(eta, ham; n_w=n_w)
    println("Aws done")
    return response.(omegas', Ref(eta), QIndices, Ref(Aws), Ref(ham.bandwidth); n_w=n_w, mu=mu, T=T)
end































end
