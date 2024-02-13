# this program reads a file with partial inbreeding coefficients 
# in the format  ancestor individual, value 
# then reads into a sparse array (T)
# then it reads ped from file and builds I-P
# computes K=T(I-P)
# finally writes that to a file
# andres.legarra@uscdcb.com, 7 Jun 2023
using LinearAlgebra
using StatsBase
using DelimitedFiles
using SparseArrays
using Printf

# for SparseArrays: https://docs.julialang.org/en/v1/stdlib/SparseArrays/

if length(ARGS)<3
    println(Base.stderr,"usage: Julia pedfile PartialInbreedingFile maxval")
    println(Base.stderr," e.g. julia TtimesImP.jl ../pedMTR.ren ../PartialInbreedingMTR.txt 0.01")
    println("exiting")
    println()
    exit()
end
pedfile=ARGS[1]
partialFfile=ARGS[2]
minval=parse(Float64,ARGS[3])
ped=Int.(readdlm(pedfile)[:,1:3])
nanim=size(ped,1)
println(Base.stderr,"animals read: ",nanim)

# read T
# cheap option in programming 
@time temp=readdlm(partialFfile,skipstart=1)
println(Base.stderr,"coefficients in T read: ",size(temp,1))
anc=Int.(temp[:,1])
id=Int.(temp[:,2])
val=temp[:,3]
# alternatively
#= anc=Array{Int64,1}()
id=Int64[]
val=Float64[]
f=open(partialFfile,"r");
@time for (i,line) in enumerate(eachline(f))
    if(i>1)
        a,b,c,d=split(line)
        #println(a,b,c)
        anc1=parse(Int,a)
        id1=parse(Int,b)
        val1=parse(Float64,c)
        push!(anc,anc1)
        push!(id,id1)
        push!(val,val1)
    end
end    
close(f)
println(Base.stderr,"coefficients in T read: ",size(anc,1))
 =#

# add an extra zero element at the last animal so that the matrix has 
# the right dimensions
anc=[anc;nanim]
id=[id;nanim]
val=[val;0]

# id indicates rows; anc indicates cols
@time T=sparse(id,anc,val)
println(Base.stderr,"sparse T built: ",nnz(T))

# build I-P

ImP=spzeros(nanim,nanim)+I
@time for i in 1:nanim
    for j in 2:3
        if ped[i,j]!=0
            ImP[i,ped[i,j]] = -0.5
        end
    end
end
println(Base.stderr,"I-P computed: ")

# compute
@time K=T*ImP
dropzeros!(K)
println(Base.stderr,"K computed: ")

fout=open("K.txt","w")

# now we need to extract the non-zero coefficients and write them
@time for i in 1:nanim
    # extract row
    Ki=K[i,:]
    # how many elements?
    nonZero=nnz(Ki)
    if nonZero>0
        # split index and values
        ancOfI,coefKofI=findnz(Ki)
        if any(coefKofI .> minval)
            # write down how many non-zero coeffs (nonZero) 
            # and how many above minval (count)
            @printf(fout,"%9d%9d%9d",i,nonZero,count(abs.(coefKofI) .> minval))
            for j in 1:nonZero
                if abs(coefKofI[j]) > minval
                    @printf(fout,"%20.10f%9d",coefKofI[j],ancOfI[j])
                end
            end
            # end of line
            println(fout)
        end
    end
end
close(fout)

println(Base.stderr,"K.txt written")
