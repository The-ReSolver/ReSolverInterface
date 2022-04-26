# The ReSolver Interface
The generic interface that should be defined for a custom implementation of the functionality requried to perform the Resolvent optimisation.

## Fields
Below is a list of the required and optional methods that need to be defined for a concrete implementation of the abstract fields type. The absolute bare minimum is included here, as such there are quite a few options related to the abstract array interface that are not listed here. Details about such things can be found in the [interfaces documentation](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array).

| Methods to implement                                          |                                        | Brief description |
|:------------------------------------------------------------- |:-------------------------------------- |:---------------------------------------------- |
| `size(A)`                                                     |                                        | Returns a tuple containing the dimensions of `A` |
| `getindex(A, i::Int)`                                         |                                        | (if `IndexLinear`) Linear scalar indexing |
| `setindex!(A, v, i::Int)`                                     |                                        | (if `IndexLinear`) Scalar indexed assignment |
| `ddy!(dxdy::AbstractArray{T, N}, x::AbstractArray{T, N})`     |                                        | Derivative with respect to wall-normal direction |
| `ddz!(dxdz::AbstractArray{T, N}, x::AbstractArray{T, N})`     |                                        | Derivative with respect to spanwise direction |
| `ddt!(dxdt::AbstractArray{T, N}, x::AbstractArray{T, N})`     |                                        | Derivative with respect to temporal direction |
| **Optional methods**                                          | **Default definition**                 | **Brief description** |
| `d2dy2!(d2xdy2::AbstractArray{T, N}, x::AbstractArray{T, N})` | `ddy!(d2xdy2, ddy!(similar(x), x))`    | Second derivative with respect to wall-normal direction |
| `d2dz2!(d2zdy2::AbstractArray{T, N}, x::AbstractArray{T, N})` | `ddz!(d2xdz2, ddz!(similar(x), x))`    | Second derivative with respect to spanwise direction |
| `d2dt2!(d2tdy2::AbstractArray{T, N}, x::AbstractArray{T, N})` | `ddt!(d2xdt2, ddt!(similar(x), x))`    | Second derivative with respect to temporal direction |

TO CHECK:
 - ensure `IndexStyle(::Type{<:AbstractArray})` can be defined without breaking anything, otherwise it needs to implemented explicitly by user.
 - `similar` method required?
 - how to get physical grid size?
 - how to initialise Laplacian?
 - how to solve Poisson equation?
 - how to "hide" convolution? (probably just define method for "*").
