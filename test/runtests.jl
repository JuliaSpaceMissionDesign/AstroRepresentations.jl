
using Test

const GROUP = get(ENV, "GROUP", "All")
const TEST_GROUPS = ("Representations", )

const TEST_TIME = Dict{String, Float64}()

for tg in TEST_GROUPS
    if GROUP == tg || GROUP == "All"
        TEST_TIME[tg] = @elapsed begin 
            @testset "$tg" verbose=true begin 
                include("$tg.jl")
            end
        end
    end
end

println("\033[1mTest Execution Time Symmary:\033[0m")

for pair in TEST_TIME 
    println("Group: $(pair[1]), elapsed time $(pair[2]) s.")
end