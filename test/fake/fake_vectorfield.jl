divergence!(div_u::VectorField, u::VectorField) = ([div_u[i] .= 5im.*u[i] for i in eachindex(u)]; return div_u)
laplacian!(Δu::VectorField, u::VectorField) = ([Δu[i] .= 5im.*u[i] for i in eachindex(u)]; return Δu)
