using LinearAlgebra

function SpinMats(spin::Rational{Int64}) :: Vector{Matrix{ComplexF64}}

    δ(x,y) = isapprox(x, y, atol=1e-5, rtol=1e-5)

    dims    =   Int(2*spin+1)
    Sx      =   zeros(ComplexF64, dims, dims)
    Sy      =   zeros(ComplexF64, dims, dims)
    Sz      =   zeros(ComplexF64, dims, dims)
    S4      =   Matrix{ComplexF64}(1.0I, dims, dims)

    spins   =   Float64.([spin - (sz - 1) for sz in 1:dims])
    Sz      =   diagm(spins)

    for i in spins
        for j in spins
            Sx[Int(1+spin-i), Int(1+spin-j)]  =         (1/2) * (δ(i, j+1) + δ(i, j-1)) * sqrt(spin * (spin+1) - (i * j))
            Sy[Int(1+spin-i), Int(1+spin-j)]  =   -im * (1/2) * (δ(i, j+1) - δ(i, j-1)) * sqrt(spin * (spin+1) - (i * j))
        end
    end

    return [Sx, Sy, Sz, S4]

end