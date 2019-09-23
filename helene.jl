using DataFrames
using Distributions
using Printf

function bayesPR(genoTrain, phenoTrain, snpInfo, chrs, fixedRegSize, varGenotypic, varResidual, chainLength, burnIn, outputFreq, onScreen)
    SNPgroups = prepRegionData(snpInfo, chrs, genoTrain, fixedRegSize)
    these2Keep = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations
    nRegions    = length(SNPgroups)
    println("number of regions: ", nRegions)
    dfEffectVar = 4.0
    dfRes       = 4.0
    X           = convert(Array{Float64}, genoTrain[2:end])
    println("X is this size", size(X))
    y           = convert(Array{Float64}, phenoTrain)
    println("y is this size", size(y))
    nTraits, nRecords , nMarkers   = size(y,2), size(y,1), size(X,2)
    fileControlSt(fixedRegSize)
    p           = mean(X,dims=1)./2.0
    sum2pq      = sum(2*(1 .- p).*p)
    if varGenotypic==0.0
        varBeta      = fill(0.0005, nRegions)
        else varBeta = fill(varGenotypic/sum2pq, nRegions)
    end
    if varResidual==0.0
        varResidual  = 0.0005
    end
    scaleVar        = varBeta[1]*(dfEffectVar-2.0)/dfEffectVar
    νS_β            = scaleVar*dfEffectVar
    df_β            = dfEffectVar
    scaleRes        = varResidual*(dfRes-2.0)/dfRes
    νS_e            = scaleRes*dfRes
    df_e            = dfRes
    tempBetaVec     = zeros(Float64,nMarkers) #initial values as "0"
    μ               = mean(y)
    X              .-= ones(Float64,nRecords)*2p
    xpx             = diag(X'X)
    ycorr           = y .- μ
    #MCMC starts here
    for iter in 1:chainLength
        #sample residual variance
        varE = sampleVarE(νS_e,ycorr,df_e,nRecords)
        #sample intercept
        ycorr    .+= μ
        rhs      = sum(ycorr)
        invLhs   = 1.0/nRecords
        meanMu   = rhs*invLhs
        μ        = rand(Normal(meanMu,sqrt(invLhs*varE)))
        ycorr    .-= μ
        for r in 1:nRegions
            theseLoci = SNPgroups[r]
            regionSize = length(theseLoci)
            λ_r = varE/varBeta[r]
            for l in theseLoci::UnitRange{Int64}
                BLAS.axpy!(tempBetaVec[l], view(X,:,l), ycorr)
                rhs = view(X,:,l)'*ycorr
                lhs = xpx[l] + λ_r
                meanBeta = lhs\rhs
                tempBetaVec[l] = sampleBeta(meanBeta, lhs, varE)
                BLAS.axpy!(-1*tempBetaVec[l], view(X,:,l), ycorr)
            end
            varBeta[r] = sampleVarBeta(νS_β,tempBetaVec[theseLoci],df_β,regionSize)
        end
        outputControlSt(onScreen,iter,these2Keep,X,tempBetaVec,μ,varBeta,varE,fixedRegSize)
    end
#    betaFromFile =  readcsv(pwd()"/betaOut",header=false)
#    print("read $(size(betaFromFile,1)) samples for $(size(betaFromFile,2)) markers from betaOut \n")
#    bayesRegOut = mean(betaFromFile,1)
#    @printf("acc %.6f \n", cor(y,X*bayesRegOut')[1])
end

function prepRegionData(snpInfo,chrs,genoTrain,fixedRegSize)
    accRegion = 0
    accRegionVec = [0]
    SNPgroups = []
#    mapData = readtable(pwd()"/$mapFile", header=false)
    ###only for our map file
    mapData = readtable("$snpInfo", header=true, separator=' ')
    if size(mapData,2)<4
        mapData = hcat(collect(1:size(mapData,1)),mapData,makeunique=true)
    end
    headMap = [:snpID, :chrID, :pos, :snpOrder]
    rename!(mapData, names(mapData), headMap)
    print(mapData[1:4,:])
    mapData[:snpID] = ["M$i" for i in 1:size(mapData,1)] #to convert original IDs like "HAPMAP43437-BTA-101873"
    print(mapData[1:10,:])
    ###
    mapData = mapData[mapData[:chrID] .<= chrs,:]
    # if first col in genoTrain is ID
    # I find cols that are in mapData (<chrs), and select those
#    usedLoci = intersect(names(genoTrain),Symbol.(mapData[:snpID]))
#    mapData = mapData[[find(usedLoci[i].==Symbol.(mapData[:snpID]))[] for i in 1:length(usedLoci)],:] #trim map data
#    genoX = genoTrain[vcat(Symbol("ID"),usedLoci)]    #trim genoData
#     genoX = genoTrain[[1; [find(i -> i == j, names(genoTrain))[] for j in [Symbol(mapData[:snpID][i]) for i in 1:size(mapData,1)]]]]
    #genoX = genoTrain[[find(i -> i == j, names(genoTrain))[] for j in [Symbol(mapData[:snpID][i]) for i in 1:size(mapData,1)]]]
    #genoX = genoTrain
    totLoci = size(genoTrain[2:end],2) # first col is ID
    snpInfoFinal = DataFrame(Any, 0, 3)
    if fixedRegSize==99
        println("fixedRedSize $fixedRegSize")
        snpInfoFinal = mapData[:,[:snpID,:snpOrder,:chrID]]
        accRegion    = length(unique(mapData[:chrID]))
        elseif fixedRegSize==9999
            snpInfoFinal = mapData[:,[:snpID,:snpOrder,:chrID]]
            snpInfoFinal[:,:chrID]  = 1 #was ".=1"
            accRegion    = 1
        else
        for c in 1:chrs
            thisChr = mapData[mapData[:chrID] .== c,:]
            totLociChr = size(thisChr,1)
            TotRegions = ceil(Int,totLociChr/fixedRegSize)
            accRegion += TotRegions
            push!(accRegionVec, accRegion)
            tempGroups = sort(repeat(collect(accRegionVec[c]+1:accRegionVec[c+1]),fixedRegSize))
            snpInfo = DataFrame(Any, length(tempGroups), 3)
            snpInfo[1:totLociChr,1] = collect(1:totLociChr)
            snpInfo[1:totLociChr,2] = thisChr[:snpID]
            snpInfo[:,3] = tempGroups
            dropmissing!(snpInfo)
            snpInfoFinal = vcat(snpInfoFinal,snpInfo)
            @printf("chr %.0f has %.0f groups \n", c, TotRegions)
            println(by(snpInfo, :x3, nrow)[:,2])
        end
        end  #ends if control flow
#    print(snpInfoFinal)
    writecsv("snpInfo",convert(Array,snpInfoFinal))
    for g in 1:accRegion
        push!(SNPgroups,searchsorted(snpInfoFinal[:,3], g))
    end
    return SNPgroups #, genoX
end

function outputControlSt(onScreen,iter,these2Keep,X,tempBetaVec,μ,varBeta,varE,fixedRegSize)
    if iter in these2Keep
        out0 = open(pwd()*"/muOutST$fixedRegSize", "a")
        writecsv(out0, μ)
        close(out0) 
        out1 = open(pwd()*"/betaOutST$fixedRegSize", "a")
        writecsv(out1, tempBetaVec')
        close(out1)
        out2 = open(pwd()*"/varBetaOutST$fixedRegSize", "a")
        writecsv(out2, varBeta')
        close(out2)
        out3 = open(pwd()*"/varEOutST$fixedRegSize", "a")
        writecsv(out3, varE)
        close(out3)
        varUhat = var(X*tempBetaVec)
        out4 = open(pwd()*"/varUhatOutST$fixedRegSize", "a")
        writecsv(out4, varUhat)
        close(out4)
        if onScreen==true
#            varU = var(X*tempBetaVec)
            @printf("iter %s varUhat %.2f varE %.2f\n", iter, varUhat, varE)
        elseif onScreen==false
             @printf("iter %s\n", iter)
        end
    end
end

function fileControlSt(fixedRegSize)
    for f in ["muOutST$fixedRegSize" "betaOutST$fixedRegSize" "varBetaOutST$fixedRegSize" "varEOutST$fixedRegSize" "varUhatOutST$fixedRegSize"]
        if isfile(f)==true
            rm(f)
            println("$f removed")
        end
    end
end

function outputControl(nTraits,onScreen,iter,these2Keep,X,tempBetaMat,μ,covBeta,Rmat,fixedRegSize,nRegions)
    if iter in these2Keep
        out0 = open(pwd()*"/muOutMTPR$fixedRegSize", "a")
        writecsv(out0, μ)
        close(out0)
        for t in 1:nTraits
            out1 = open(pwd()*"/beta"*"$t"*"OutMTPR$fixedRegSize", "a")
            writecsv(out1, tempBetaMat[t,:]')
            close(out1)
            out2 = open(pwd()*"/varBeta"*"$t"*"OutMTPR$fixedRegSize", "a")
            printThis = [vcat(covBeta[r]...)[t^2] for r in 1:nRegions]'
            writecsv(out2, printThis) #works only for bivariate
            close(out2)
        end
        outCov = open(pwd()*"/covBetaOutMTPR$fixedRegSize", "a")
        printThis = [vcat(covBeta[r]...)[3] for r in 1:nRegions]'
        writecsv(outCov, printThis) #works only for bivariate
        close(outCov)
        out3 = open(pwd()*"/varEOutMTPR$fixedRegSize", "a")
        writecsv(out3, vec(Rmat)')
        close(out3)
        coVarUhat = cov(X*tempBetaMat')
        out4 = open(pwd()*"/coVarUhatOutMTPR$fixedRegSize", "a")
        writecsv(out4, vec(coVarUhat)')
        close(out4)    
        if onScreen==true
#            varU = cov(X*tempBetaMat')
            println("iter $iter \ncoVarUhat: $coVarUhat \nvarE: $Rmat \n")
        elseif onScreen==false
             @printf("iter %s\n", iter)
        end
    end
end

function fileControl(nTraits,fixedRegSize)
    files2Remove = ["muOutMTPR$fixedRegSize", "varEOutMTPR$fixedRegSize", "covBetaOutMTPR$fixedRegSize", "coVarUhatOutMTPR$fixedRegSize"]
    for t in 1:nTraits
        push!(files2Remove,"beta"*"$t"*"OutMTPR$fixedRegSize")
        push!(files2Remove,"varBeta"*"$t"*"OutMTPR$fixedRegSize")
    end
    for f in files2Remove
        if isfile(f)==true
            rm(f)
            println("$f removed")
        end
    end
end

function sampleBeta(meanBeta, lhs, varE)
    return rand(Normal(meanBeta,sqrt(lhs\varE)))
end

function sampleVarBeta(νS_β,whichLoci,df_β,regionSize)
    return((νS_β + dot(whichLoci,whichLoci))/rand(Chisq(df_β + regionSize)))
end
function sampleVarE(νS_e,yCorVec,df_e,nRecords)
    return((νS_e + dot(yCorVec,yCorVec))/rand(Chisq(df_e + nRecords)))
end
function sampleCovBeta(dfβ, regionSize, Vb , tempBetaMat, theseLoci)
    Sb = tempBetaMat[:,theseLoci]*tempBetaMat[:,theseLoci]'
    return rand(InverseWishart(dfβ + regionSize, Vb + Sb))
end
function sampleCovarE(dfR, nRecords, VR, ycorr1, ycorr2)
     Sr = [ycorr1 ycorr2]'* [ycorr1 ycorr2]
    return rand(InverseWishart(dfR + nRecords,VR + Sr))
end

function mmeRunFast(xp,Ri,locus,xpx,ycorr1,ycorr2,invB)
#    r1 = xp*Ri[1]
#    r2 = xp*Ri[2]
#    r3 = xp*Ri[3]
#    r4 = xp*Ri[4]
    
#    rhs    = [r1 r2;r3 r4]*[ycorr1;ycorr2]

    rhs    = [xp*ycorr1*Ri[1]+xp*ycorr2*Ri[2];xp*ycorr1*Ri[3]+xp*ycorr2*Ri[4]]
    invLhs = fastInv(xpx[locus].*Ri .+ invB)
    
    meanBeta = invLhs*rhs    
    return rand(MvNormal(meanBeta,convert(Array,Symmetric(invLhs))))
end

function fastInv(a::Matrix)
    c = copy(a)

    detv = a[1] * a[4] - a[2] * a[3]
    inv_d = 1 / detv

    c[1] = a[4] * inv_d
    c[2] = -a[2] * inv_d
    c[3] = -a[3] * inv_d
    c[4] = a[1] * inv_d
    return c
end
