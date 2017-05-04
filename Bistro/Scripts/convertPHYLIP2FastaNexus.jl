## julia script to convert Phylip file to Fasta/Nexus file
## Claudia May 2017

file = "043.phy"
f = open(file)
p = split(file,".")
fasta = string(p[1],".fasta")
nexus = string(p[1],".nex")
verbose = false

firstLine = true
taxtotal = 0
taxnum = 1
taxa = String[]
seq = String[]
readNumTaxa = 0
readSeqLength = 0
for line in eachline(f)
    verbose && @show line
    if(firstLine) ## first line
        firstLine = false
        v = split(line)
        readNumTaxa = parse(Int,v[1])
        readSeqLength = parse(Int,v[2])
        verbose && @show readNumTaxa
        verbose && @show readSeqLength
        continue
    end
    v = convert(Array{String}, split(line,"  "))
    verbose && @show v
    v = filter(x->x!="",v)
    if(length(v) == 0)
        continue
    elseif(length(v) == 2) ## the first lines that have the taxon name
        push!(taxa,strip(string(v[1])))
        taxtotal += 1
        s = replace(strip(string(v[2]))," ","")
        push!(seq,s)
    elseif(length(v) == 1)
        verbose && @show taxnum
        s = replace(strip(string(v[1]))," ","")
        seq[taxnum] = string(seq[taxnum],s)
        taxnum = mod(taxnum+1,taxtotal)
        taxnum = (taxnum == 0) ? taxtotal : taxnum
        verbose && @show taxnum
    else
        error("error in line read $v")
    end
end
close(f)

taxtotal == readNumTaxa || warn("Counted $taxtotal taxa, but in file it says there are $readNumTaxa")
length(seq) == taxtotal || error("Something wrong when reading the PHYLIP file. I found $(length(seq)) sequences, and there should be $taxtotal")
length(seq[1]) == readSeqLength || warn("Counted $(length(seq[1])) sites, but in file it says there are $readSeqLength")

## writing FASTA file
g = open(fasta,"w")
for i in 1:taxtotal
    write(g,">$(taxa[i])\n")
    write(g,"$(seq[i])\n")
end
close(g)

## writing NEXUS file
g = open(nexus,"w")
write(g,"#NEXUS\n\nbegin taxa;\ndimensions ntax=$taxtotal;\ntaxlabels\n")
for i in 1:taxtotal
    write(g,"$(taxa[i])\n")
end
write(g,";\nend;\n\nbegin characters;\ndimensions nchar=$readSeqLength;\nformat datatype=dna gap=- missing=?;\nmatrix\n")
for i in 1:taxtotal
    write(g,"$(taxa[i])\t\t")
    write(g,"$(seq[i])\n")
end
write(g,";\nend;\nbegin mrbayes;\n      set autoclose=yes nowarn=yes;\n      lset nst=6;\n      outgroup 1;\n\n      mcmc ngen=1100000 printfreq=10000 samplefreq=50;\n      sumt burnin=100000;\n      sump burnin=100000;\nend;")
close(g)
