module SpinMatrices
    export SpinMats, HermitianBasis

    using LinearAlgebra

@doc """
```julia
SpinMats(spin::Rational{Int64}) --> Vector{Matrix{ComplexF64}}
```
Returns the 3 generators of SU(2) in the given spin-representation, along with the Identity matrix of d = 2*spin+1, in the format of [Sx, Sy, Sz, Id]

"""
    function SpinMats(spin::Rational{Int64}) :: Vector{Matrix{ComplexF64}}

        δ(x,y)  = isapprox(x, y, atol=1e-5, rtol=1e-5)

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


@doc """
```julia
HermitianBasis(localDim::Int64) --> Vector{Matrix{ComplexF64}}
```
Returns a vector basis of Traceless, localDim x localDim, Hermitian matrices, normalized such that `` Tr(T^{a †} . T^{b}) = (1/2)δ_{ab} ``, along with the identity. 

"""
    function HermitianBasis(localDim::Int64) :: Vector{Matrix{ComplexF64}}

        basis   =   [zeros(ComplexF64, localDim, localDim) for i in 1:localDim^2]
        count   =   1

        for i in 2:localDim
            for j in 1:i-1
                for phase in 1:2

                    basis[count][i, j]  =   Bool(mod(count-1, 2)) ? 0.5 * im : 0.5 + 0.0 * im
                    basis[count][j, i]  =   conj(basis[count][i, j])

                    count   +=  1
                end
            end
        end

        for i in 1:localDim - 1
            diag            =   vcat( ones(ComplexF64, i), -i * ones(ComplexF64, 1), zeros(ComplexF64, localDim - i - 1))
            basis[count]    =   diagm(diag) / sqrt(2 * i * (i+1))

            count   +=  1
        end

        basis[count]    =   Matrix{ComplexF64}(I, localDim, localDim)

        return basis
    end

end