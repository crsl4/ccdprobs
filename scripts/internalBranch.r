## r script to study the likelihood of internal branch
## in q 4taxon tree, given all other branch lengths

## computes alpha_x = Px1(d1x)Px2(d2x)
## identical to beta_y, so will only code one
alpha = function(x,site1, site2, d1x,d2x){
    if(site1 == x){
        if(site2 == x){
            p = (0.25+0.75*exp(-4*d1x/3)) * (0.25+0.75*exp(-4*d2x/3))
        } else {
            p = (0.25+0.75*exp(-4*d1x/3)) * (0.25-0.25*exp(-4*d2x/3))
        }
    } else {
        if(site2 == x){
            p = (0.25-0.25*exp(-4*d1x/3)) * (0.25+0.75*exp(-4*d2x/3))
        } else {
            p = (0.25-0.25*exp(-4*d1x/3)) * (0.25-0.25*exp(-4*d2x/3))
        }
    }
    return ( p )
}

## computes theta_k = alpha_A*beta_A+alpha_C*beta_C+alpha_G*beta_G+alpha_T*beta_T
theta = function(site1,site2,site3,site4,d1x,d2x,d3y,d4y){
    t = alpha('a',site1,site2,d1x,d2x) * alpha('a',site3,site4,d3y,d4y) +
        alpha('c',site1,site2,d1x,d2x) * alpha('c',site3,site4,d3y,d4y) +
            alpha('g',site1,site2,d1x,d2x) * alpha('g',site3,site4,d3y,d4y) +
                alpha('t',site1,site2,d1x,d2x) * alpha('t',site3,site4,d3y,d4y)
    return ( t )
}

## computes gamma_k = sum_x sum_y alpha_x beta_y (x!=y)
gamma = function(site1,site2,site3,site4,d1x,d2x,d3y,d4y){
    suma = 0
    for(x in c('a','c','g','t')){
        for(y in c('a','c','g','t')){
            if(x != y){
                suma = suma + alpha(x,site1,site2,d1x,d2x)*alpha(y,site3,site4,d3y,d4y)
            }
        }
    }
    return ( suma )
}

## computes f(y), which is l'(t)=0, and f'(y)
## for Newton-Raphson
## y \in (0,1)
fy = function(y,seq1,seq2,seq3,seq4, d1x,d2x,d3y,d4y){
    nsites = length(seq1)
    if(length(seq2) != nsites || length(seq3) != nsites || length(seq4) != nsites)
        stop("wrong number of sites in one sequence")

    suma = 0
    suma2 = 0
    for(i in 1:nsites){
        t = theta(seq1[i], seq2[i], seq3[i], seq4[i], d1x,d2x,d3y,d4y)
        g = gamma(seq1[i], seq2[i], seq3[i], seq4[i], d1x,d2x,d3y,d4y)
        suma = suma + (0.75*t-0.25*g)/(0.25*(t+g+(3*t-g)*y))
        suma2 = suma2 - (0.75*t-0.25*g)^2/(0.25*(t+g+(3*t-g)*y))^2
    }
    return ( list(fy=suma, fyprime=suma2) )
}


## function that will find the MLE for the internal branch length
## given the sequences and the other branch lengths
findMLE = function(seq1,seq2,seq3,seq4, d1x,d2x,d3y,d4y){
    y0 = 0.8 # chosen arbitrarily, true value is 0.9659
    it = 100 #number of iterations
    y = rep(0, it)
    y[1] = y0

    for(i in 2:it){
        f = fy(y[i-1],seq1,seq2,seq3,seq4, d1x,d2x,d3y,d4y)
        y[i] = y[i-1] - f$fy/f$fyprime
    }
    return ( y )
}

