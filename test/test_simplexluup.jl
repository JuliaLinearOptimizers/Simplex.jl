include("generators/transport_problem.jl")

function test_simplexluup()
  #examples 3.1 and 3.2 (luenberger)
  @testset "Luenberger examples" begin
    c = [4.; 1; 1]
    A = [2. 1 2; 3 3 1]
    b = [4.; 3]
    sol = [0.; 2/5; 9/5]
    x, z, status = simplexluup(c, sparse(A), b)
    @test x ≈ sol
    @test z ≈ dot(sol, c)
    @test status == :Optimal
    c = [4.; 1; 1]
    A = [-2. -1 -2; -3 -3 -1]
    b = [-4.; -3]
    sol = [0; 2/5; 9/5]
    x, z, status = simplexluup(c, sparse(A), b)
    @test x ≈ sol
    @test z ≈ dot(sol, c)
    @test status == :Optimal
    c = [-1.; 1; 0; 0; 0]
    A = [-2. 1 1 0 0; 1 -3 0 1 0; 1 1 0 0 1]
    b = [2.; 2; 4]
    sol = [7/2; 1/2; 17/2; 0; 0]
    x, z, status = simplexluup(c, sparse(A), b)
    @test x ≈ sol
    @test z ≈ dot(sol, c)
    @test status == :Optimal
  end
  @testset "Optimal solution without sufficient condition" begin
   c = [1.; 2]
   A = [1. 1]
   b = [0.]
   sol = [0.; 0]
   basis = [2]
   x, z, status = simplexluup(c, sparse(A), b, basis)
   @test x ≈ sol
   @test z ≈ dot(sol, c)
   @test status == :Optimal
  end
  @testset "Infeasibility" begin
    c = [1.; 1]
    A = [1. 1; 1 1]
    b = [1.; 2]
    x, z, status = simplexluup(c, sparse(A), b)
    @test status == :Infeasible
  end
  @testset "Bertsimas example - modified" begin
    # 3.8 - has redundant equations
    c = [1.; 1; 1; 0]
    A = [0. 4 9 0; 1 2 3 0; -1 2 6 0; 0 4 9 0; 0 0 3 4; 0 0 6 8; 0 4 9 0]
    b = [5.; 3; 2; 5; 1; 2; 5]
    sol = [0.5; 1.25; 0; 0.25]
    x, z, status = simplexluup(c, sparse(A), b)
    @test x ≈ sol
    @test z ≈ dot(sol, c)
    @test status == :Optimal
  end
end
test_simplexluup()

function test_transpluup()
  @testset "Transport tests" begin
    @testset "scale $scale_m × $scale_n" for scale_m = 10 .^(0:1), scale_n = 10 .^(0:1)
      for t = collect(1:10)
        m, n = scale_m * rand(2:4), scale_n * rand(2:4)
        A, b, c = transport_instance(m, n)

        model = Model(with_optimizer(GLPK.Optimizer))
        @variable(model, x[1:m*n] >= 0)
        @objective(model, Min, dot(x, c))
        @constraint(model, A * x .== b)
        optimize!(model)
        status = termination_status(model)
        xj = Vector{Float64}(undef,m*n)
        for i in 1:m*n
          xj[i] = value(x[i])
        end
        zj = objective_value(model)

        x, z, status = simplexluup(c, A, b, max_iter = 1000 * scale_m * scale_n)
        Δz = dot(c, x) - zj
        # @test x ≈ xj
        @test dot(c, x) ≈ zj atol = 1e-3
        @test status == :Optimal
      end
    end
  end
end
test_transpluup()
