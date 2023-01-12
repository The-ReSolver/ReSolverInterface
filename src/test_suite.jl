# This file provides the function necessary to ensure that any concrete#
# implementation of the abstract types defined in this module has the correct
# behviour from the base methods that need to implemented.


"""
    test_interface(
                grid::AbstractGrid,
                field::AbstractScalarField,
                laplace::AbstractLaplace
                ) -> nothing

Test that the required types for the interface behave as expected by the
module.
"""
function test_interface(grid::AbstractGrid, field::AbstractScalarField, laplace::AbstractLaplace)
    test_grid(grid)
    test_scalarfield(field)
    test_laplace(laplace)
end

"""
    test_grid(grid::AbstractGrid) -> nothing

Test that an `AbstractGrid` follows the expected behaviour for the interface.
"""
function test_grid(grid::AbstractGrid)

end

"""
    test_scalarfield(field::AbstractScalarField) -> nothing

Test that an `AbstractScalarField` follows the expected behaviour for the
interface.
"""
function test_scalarfield(field::AbstractScalarField)

end

"""
    test_laplace(laplace::AbstractLaplace) -> nothing

Test than an `AbstractLaplace` follows the expected behaviour for the
interface.
"""
function test_laplace(laplace::AbstractLaplace)

end
