module Gauge

    export LandauGauge, RadialGauge

    using LinearAlgebra


    function LandauGauge(r::Vector{Float64} ; direction::Vector{Float64} = [1.0, 0.0], normal::Vector{Float64} = [0.0, 1.0], B::Float64 = 1.0) :: Vector{Float64}

        distance    =   dot(r, normal) / norm(normal)
        return direction .* B .* distance / norm(direction)
    end


    function RadialGauge(r::Vector{Float64} ; B::Vector{Float64}) :: Vector{Float64}

        return -(1/2) * cross(r, B)
    end


end