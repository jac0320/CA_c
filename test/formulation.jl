@testset "Deterministic Solver" begin
    @testset "14-Bus Case Network" begin
        resultProb = adcc(PROBLEM="ieee14",
                            MODEL="network",
                            ALG="regular",
                            S=10,
                            EPS=0.0,
                            STOCHMODE="evolving",
                            FEATURES=["sample-based-risk", "surge-load-shed"])
        @test isapprox(getobjectivevalue(resultProb.model), 609.6159; atol = 1)
    end

    @testset "14-Bus Case DC" begin
        resultProb = adcc(PROBLEM="ieee14",
                            MODEL="dc",
                            ALG="regular",
                            S=10,
                            EPS=0.0,
                            STOCHMODE="evolving",
                            FEATURES=["sample-based-risk", "surge-load-shed"])
        @test isapprox(getobjectivevalue(resultProb.model), 610.7427; atol = 1)
    end

    @testset "14-Bus Case Network EPS0.1" begin
        resultProb = adcc(PROBLEM="ieee14",
                            MODEL="network",
                            ALG="regular",
                            S=10,
                            EPS=0.1,
                            STOCHMODE="evolving",
                            FEATURES=["sample-based-risk", "surge-load-shed"])
        @test isapprox(getobjectivevalue(resultProb.model), 385.4993; atol = 1)
    end

    @testset "14-Bus Case DC EPS0.1" begin
        resultProb = adcc(PROBLEM="ieee14",
                            MODEL="dc",
                            ALG="regular",
                            S=10,
                            EPS=0.1,
                            STOCHMODE="evolving",
                            FEATURES=["sample-based-risk", "surge-load-shed"])
        @test isapprox(getobjectivevalue(resultProb.model), 385.4993; atol = 1)
    end
end
