# Tests the solvers on simple problems in various types
@testset "Solve simple problems" begin
  # All solvers here should accept (x, A, b, c)
  solvers = [
    primalbasic!
  ]

  # All problems here should generate (A, b, c, sol)
  # All problems here should succeed
  problems = [
    luenberger1,
    luenberger2,
    luenberger3
  ]

  Base.rationalize(x :: Rational) = x

  for s in solvers, p in problems
    A, b, c, sol = p()
    m, n = size(A)
    for (T,Tc) in [
      (Float64, Float64),
      (Rational{Int}, rationalize),
      (Float16, Float16),
      (BigFloat, BigFloat)
    ]
      x = zeros(T, n)
      status, _ = s(x, Tc.(A), Tc.(b), Tc.(c))
      @test eltype(x) == T
      @test status == :solved
      @test x ≈ Tc.(sol)
    end
  end
end