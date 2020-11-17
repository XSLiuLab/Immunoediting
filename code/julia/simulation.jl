##the following code is modified from the paper provided : https://doi.org/10.1038/s41588-020-0687-1

using Distributions, StatsBase, DataFrames, GLM, CSV, DelimitedFiles

function getFitness(n)
    (1 + s*n)
end

mutable struct cancercell
    mutations::Array{Int64,1}
    mut_neop::Array{Bool,1}
    fitness::Float64
    neonumber::Float64
    escaped::Array{Bool,1}
end

function newmutations(cancercell, mutID, p, pesc)
    cancercell.mutations = append!(cancercell.mutations, mutID)
    mutID = mutID + 1
    
    neoep = rand()<p
    if neoep
        cancercell.neonumber = cancercell.neonumber + 1
        cancercell.fitness = getFitness(cancercell.neonumber) #fitness is affected by the number of mutations
        cancercell.mut_neop = append!(cancercell.mut_neop,true)
    else
        cancercell.mut_neop = append!(cancercell.mut_neop,false)
    end

    mut_escape = rand() < pesc
    if mut_escape 
        cancercell.escaped = append!(cancercell.escaped,true)
    else
        cancercell.escaped = append!(cancercell.escaped,false)
    end

    return cancercell, mutID
end

function copycell(cancercellold::cancercell)
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.mut_neop), copy(cancercellold.fitness), copy(cancercellold.neonumber), copy(cancercellold.escaped))
end

function start_population(p, pesc, initial_mut)
    mutID = 1
    N = 1
    cells = cancercell[]
    #muts = Dict{Int64, Float64}()
    push!(cells,cancercell([],[],1,0,[]))
    for i=1:initial_mut
        cells[1],mutID = newmutations(cells[1],mutID, p, pesc)
        #muts[mutID-1] = cells[1].mut_neop[1]
    end

    return cells, mutID,  N

end

function birthdeath_neoep(b0, d0, Nmax, p, initial_mut, mu, pesc)

    dmax = d0 #dmax is updated throughout, starts from d0

    #initialize arrays and parameters
    cells, mutID, N = start_population(p, pesc, initial_mut )
    Nvec = Int64[]
    push!(Nvec,N)
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)

    while (N < Nmax) & (t < 300) #set so we can exit simulation where there is a lot of death

        #pick a random cell
        randcell = rand(1:N)
        Nt = N
        
        #a cell's immunogenicity depends on its fitness, i.e. the summed antigenicity of neoepitopes
        d = max(0, (d0 - b0)*cells[randcell].fitness + b0)

        if true in cells[randcell].escaped # discard effect of fitness if the cell escape
            d = d0
        end

        if (d > dmax) #update dmax to keep track of the highest death rate in the whole population
            dmax = d
        end

        Rmax = b0 + dmax

        r = rand(Uniform(0,Rmax)) #Pick which reaction should happen to cell     

        # If r < birthrate, a birth event happens: a new cell is created and randcell updated as a new one
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))

            #add new mutations to both new cells, the number of mutations is Poisson distributed
            for i=1:(rand(Poisson(mu)))
                cells[randcell],mutID = newmutations(cells[randcell],mutID, p, pesc)
                #muts[mutID-1] = cells[randcell].mut_neop[1]
            end
            for i=1:(rand(Poisson(mu)))
                cells[end],mutID = newmutations(cells[end],mutID, p, pesc)
                #muts[mutID-1] = cells[end].mut_neop[1]
            end

            push!(Nvec, N)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
            
        end

        #if r has neither birth or death (only possible if it is a non-immunogenic cell), nothing happens
        if  (b0+d)<= r
          push!(Nvec, N)
          Δt =  1/(Rmax * Nt) .* - log(rand())
          t = t + Δt
          push!(tvec,t)
        end

        #death event if r > b but < d
        if b0 <= r < (b0+d)

            #population decreases by 1, overall fitness score also decreases if it was non-zero
            N = N - 1

            #remove deleted cell
            deleteat!(cells,randcell)
            push!(Nvec,N)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
        end

        #if every cell dies, restart simulation from a single cell again
        if (N == 0)
            cells, mutID, N = start_population(p, pesc, initial_mut)
            push!(Nvec,N)
            push!(tvec,t)
        end

    end
    
    return Nvec, tvec, mutID, cells
end

function process_mutations(cells)
    mutVec = []
    mut_is_neo = []
    mut_is_escape = []
    for i=1:length(cells)
        append!(mutVec, cells[i].mutations)
        append!(mut_is_neo, cells[i].mut_neop)
        append!(mut_is_escape, cells[i].escaped)

    end

    #detMutDict = filter((k, v) -> v > detLim, countmap(mutVec))
    detMutDict = countmap(mutVec)

    println("Mutations processed for ", length(cells), " cells.")
    return detMutDict, mutVec, mut_is_neo, mut_is_escape

end

b0 = 1
d0 = 0.1
popSize = 1e5
p = 0.1
initial_mut = 0
mu = 1
pesc = 1e-6

s=<s>
for i=1:200
    Nvec, tvec, mutID, cells = birthdeath_neoep(b0, d0, popSize, p, initial_mut, mu,pesc);
    outNDFsim = DataFrame(t=tvec, N=Nvec)
    CSV.write("/public/slst/home/wutao2/julia_simulation/out/minus_"*string(-s)*"/preIT_"*string(i)*".txt", outNDFsim)

    detMutDict, mutVec, mut_is_neo, mut_is_escape = process_mutations(cells)
    writedlm("/public/slst/home/wutao2/julia_simulation/out/minus_"*string(-s)*"/vaf_preIT_"*string(i)*".txt",detMutDict)

    mut_summ = DataFrame(mut=mutVec,mut_neo=mut_is_neo,mut_escaped=mut_is_escape)
    CSV.write("/public/slst/home/wutao2/julia_simulation/out/minus_"*string(-s)*"/mutsumm_"*string(i)*".txt", mut_summ)

end
