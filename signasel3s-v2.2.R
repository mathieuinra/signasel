####################
DO NOT USE THIS FILE, IT IS ERROR PRONE
####################

















################################################################
## signasel3s-v2.1.R:
## Reprise plus propre de l'utilisation de i1 S1 i2 S2: les calculs de
## proba et de matrice de recurrence se font en utilisant la frequence
## p et la taille N uniquement. C'est au niveau du calcul de la
## vraisemblance qu'on utilise i1 S1 etc. Permet de faire des choix
## sur p=i1/S1 ou somme sur les p, etc., cf testne et slatkin...
################################################################
## signasel3s-v1.1.R: prise en compte explicite de S1 S2
################################################################
## PARADIGME:
##### A REVOIR !! ###
################################################################
## i1 et i2 sont TOUJOURS relatifs a S1 et S2 n1 ne sert a rien mais
## on le garde pour compatibilite et coherence, cf ci-dessous.
## Les donnees sont de la forme:
## For each sample (row) the columns are: g, i, S, N:
## g1 i1 S1 N1
## g2 i2 S2 N2
## g3 i3 S3 N3
## etc.
## with g: generation, i: number of allele copies, S: sample size,
## N: (effective) population size
## De g1 a g2:
## depart en g1 frequence p=i1/S1 puis binomiale B(N2,p) sauf en g2 ou
## B(S2,p). On pourrait faire B(N2,p) en g2 puis tirer un echantillon
## de taille S2 mais autant faire comme si la pop etait de taille S2
## (ne change rien).
## De g2 a g3:
## depart en g2 frequence p=i2/S2 puis binomiale B(N3,p) etc.
################################################################
## TWO ALLELES ONLY
## DIPLOID VERSION
################################################################
## Version 3:
################################################################
## Objectifs:
## 1/ Information sur generations intermediaires
## 2/ Combiner estimation Ne ET selection (sur =/= marqueurs)
## 3/ Souplesse recursion a/ exacte b/ diffusion c/ simul
################################################################
## This program is aimed at detecting selection by computing the
## likelihood of SNP allele frequencies.
##
## GENERAL NOTATIONS:
##
## '1' and '2' refer to the 'first' (initial) and 'last' (final)
## populations, respectively, regardless of the number of generations
## between them.
## ng: total number of generations from population 1 to population 2.
## So, ng = 1 if pop 2 is the direct offspring of pop 1. More generally,
## ng = t2 - t1 if pop 1 lives at time t1 and pop 2 at time t2.
## n1, n2: sizes of populations 1 and 2, resp. Warning: if ng > 1, the
## size of population(s) at intermediate(s) generation(s) is n2, so:
## generation 1: pop size n1->n2
## generations 2 to (ng-1): pop size n2->n2
## generation ng: pop size n2->S2 (if ng>1) or n1->S2 (if ng=1)
## S1, S2: sizes of the samples taken from populations 1 and 2, resp.
## NB: Computations are made as if the total size of population 2 was S2
## (not n2). This is strictly equivalent to first considering n2
## individuals then drawing S2 individuals among those n2, and obviously
## saves time. Important: as a generalization S2 > n2 is allowed.
## i1, i2: number of copies of the allele in samples S1 and S2, resp.
## s: coefficient of selection, such that the fitnesses are 1, 1+s,
## 1+2s for genotypes aa, aA, AA (in the case where A is the
## positively selected allele)
##
## IMPORTANT: population sizes (n1, n2) and sample sizes (S1, S2) are
## expressed in number of (diploid) individuals, so allele frequency
## reads p = i / 2n
##
################################################################
VecProb <- function(px, N, s) {
################################################################
    ## Gives the vector of probabilities of all possible allele
    ## numbers (0 to 2N) in a population of size N, when the frequency
    ## at previous generation was px and selection coefficient is s,
    ## according to Wrigt-Fisher model. This vector will serve as a
    ## column of the recursion matrix, so f(i,j) is element a(j,i) of
    ## the matrix.
    ## px: allele frequency at generations t.
    ## N: population sizes at generations (t+1)
    ## s: coefficient of selection
    ## NB: px est ici une variable, fixee par la fonction qui cree la
    ## matrice comme px=k/2/Nt (Nt taille de pop a la generation t,
    ## qui n'a pas besoin d'etre connue ici.
################################################################
    ## checkparam(px, N, s)
    ## Fitnesses of the 3 genotypes A1A1, A1A2, A2A2
    ## (A1 is the selected allele)
    w11 <- 1 + 2 * s
    w12 <- 1 + s
    w22 <- 1
    a <- w11 - w12
    b <- w12 - w22
    c <- a - b
    ## Allele frequency after selection at generations t.
    pxprime <- 1. * (px * (px * a + w12)) / (2. * px * b + w22 + px * px * c)
    ## Probabilities of all allele numbers at generations (t+1)
    prob <- dbinom(0:(2*N), 2*N, pxprime)
    return(prob)
}
################################################################
MatProb <- function(N, s, ng){
################################################################
  ## Create the Wright-Fisher recursion matrix over ng generations
  ## s: coefficient of selection
  ## N: population size (Ne), for generations 1 to ng-1
  ## ng: number of generations (iterations)
################################################################
    ## checkparam(ng)
    mat <- NULL
    ## Make the (2N)*(2N) recursion matrix from one generation to the
    ## next (pop size is always N)
    for (k in 0:(2*N)) {
        px <- 1.*k/2./N
        mat <- cbind(mat, VecProb(px, N, s))
    }
    ## iterate the matrix over generations
    ## nb: ng is the number of loops
    for (g in 2:ng) {                   # check 2:ng
        mat <- mat %*% mat
    }
    # print("Check!")
    return(mat)
}
################################################################
Like0 <- function(matp0, i1, S1, i2, S2, N, ng, s) {
################################################################
    ## Calculates the Likelihood, ie simply multiplies the recursion
    ## matrix by the vector of initial frequencies
################################################################
    ## p0 is the supposedly known allele frequency in the pop at
    ## generation 1, from which the (known) sample (i1,S1) was drawn,
    ## and from which the iteration starts. It is the responsability
    ## of the calling function to choose p0.
    ## p0 must be passed as a matrix of the form [[prob,value]]
    ## rappel/slatkin:
    ## L=Pr(p1|N)Pr(i1|p1,S1)Pr(i2|i1,p1,S1,S2,N)
    ## ici on ne s'occupe pas de Pr(p1|N)
################################################################
        ## non, il faut traiter p0 ici sinon ca fait recalculer la matrice a chaque fois....??
    ## en fait il faut donner toute la distribution supposee de p0 d'un coup ici (matrice ?)
################################################################
    ## NB: au lieu de fixer p0 on pourrait donner i0, le 'vrai' nombre
    ## d'alleles dans N correspondant à p0...? OUI !! ??
################################################################
    ## 'negative' selection is done by symetry, ie positive selection
    ## for the other allele
    # print("place?")
    if (s<0) {
        s <- -s
        i1 <- 2 * S1 - i1
        i2 <- 2 * S2 - i2
    }
    stopifnot(s>=0)
    ## c'est ici qu'on traite explicitement i1 S1 i2 S2...
    ## Q1: choix de p?
    ## Q2: echantillonnage en t1?
    ## Q3: echantillonnage en t2?
    ## ...
    mat <- MatProb(N, s, ng)
    ## placer boucle sur p0 ici !!
    probp0 <- matp0[,1]
    valp0 <- matp0[,2]
    ## normalization
    probp0 <- probp0 / sum(probp0)
    stopifnot(sum(probp0) == 1)
    L <- 0
    for (kp0 in seq_along(probp0)) {
        p0 <- valp0[kp0]
        ## probability of sampling at t1
        ps1 = dbinom(i1, 2*S1, p0)
        ## probability pop at t1 (N) -> pop at t2 (N) (vector)
        v1 <- rep(0, 2*N+1)
        i0 <- round(p0 * 2 * N)
        v1[i0+1] <- 1 ## i+1 here because i is in 0:(2n) but the first
                      ## element of a vector is 1 in R
        v2 <- mat %*% v1
        ## probability of sampling at t2
        ps2 <- 0
        for (k in 0:(2*N)) {
            ps2 <- ps2 + dbinom(i2, 2*S2, v2[k+1])
        }
        L <- L + probp0[kp0] * ps1 * ps2
    }
    return(L)
}
################################################################

  ## ############################################################
  ## enchainer les like = produit des like, mais avec le meme s!!
  ## OK!
  ## ############################################################

## Big question: multimarqueurs: maximiser N par marqueur ou par intervalle de temps ??
## depend si on cherche N ou s ??

################################################################
WFLike2S <- function(s, data) {
################################################################
  ## Calculates the Likelihood for multiple time samples
  ## Data is assumed to have the format of a matrix with nrow=number
  ## of samples and for each sample the columns: g, i, S, N:
  ## g1 i1 S1 N1
  ## g2 i2 S2 N2
  ## g3 i3 S3 N3
  ## etc.
  ## with g: generation, i: number of allele copies, S: sample size,
  ## N: (effective) population size
  ## We simply call WFLike1S from g1 to g2 (ng=g2-g1), then from
  ## g2 to g3, etc. and multiply the resulting probabilities. s is the
  ## same for all samples and is maximized over all generations.
################################################################
  g <- data[,1]
  i <- data[,2]
  S <- data[,3]
  N <- data[,4]
  p2 <- 1
  for (k in 1:(length(g)-1)) {
      p2 <- p2 *
          WFLike1S(i[k], S[k], i[k+1], S[k+1], N[k+1],g[k+1]-g[k], s)
  }
  return(p2)
}
################################################################

## ############################################################
  ## ce n'est pas bon parce qu'on recalcule mat a chaque fois, il
  ## faudrait le faire avant ??
  ## ==> on verra plus tard, utiliser tel quel d'abord
  ## ############################################################

################################################################
WFMaxiLike2S <- function(data, maxs) {
################################################################
  ## Finds the maximum likelihood on s, returns L(smax), smax
  ## Uses the 'optimize' function of R
  ## WARNING: this function does not find the max in some
  ## cases. Must be checked
  ##
  ## IMPORTANT: maxs gives min and max possible values for s 
  ## must be realistic (??)  
################################################################
  res <- optimize(WFLike2S, interval=c(-maxs, maxs), 
                  data,
                  maximum = TRUE)
  warn <- 0
  ## Issue warning message if maxs is reached
  if ((res$maximum - maxs)**2 < 1e-6) {
    print(paste("Warning: Maximum value smax = ", maxs,
                  " was reached."))
    warn <- 1
    ## it would be possible to launch something here (exhaustive search?)
    ## however i doubt it improves, because likelihood is continuous
    ## and has a single maximum (this must be checked!) hence we are
    ## at the extreme s value. Still it is BAD (unrealistic?)
    ## return value 
  }
  return(c(res$objective, res$maximum, warn))
}
################################################################
checkparam <- function(data) {
################################################################
  ## Data is assumed to have the format of a matrix with nrow=number
  ## of samples and for each sample the columns: g, i, S, N:
  ## g1 i1 S1 N1
  ## g2 i2 S2 N2
  ## g3 i3 S3 N3
  ## etc.
  ## with g: generation, i: number of allele copies, S: sample size,
  ## N: (effective) population size
################################################################
  g <- data[,1]
  i <- data[,2]
  S <- data[,3]
  n <- data[,4]
  stopifnot(all(data[,c(1,3,4)]>0))
  stopifnot(all(data[,2]>=0))
  stopifnot(all(S<=n))
  stopifnot(all(i<=(2*S)))
}
################################################################
signaseltest <- function(data, maxs = 1) {
################################################################
  ## Compute the test statistics to detect selection
  ## Data is a matrix with nrow = number of samples
  ## For each sample (row) the columns are: g, i, S, N:
  ## g1 i1 S1 N1
  ## g2 i2 S2 N2
  ## g3 i3 S3 N3
  ## etc.
  ## with g: generation, i: number of allele copies, S: sample size,
  ## N: (effective) population size
################################################################
  ## verify parameter values
  checkparam(data)
  ## likelihood of the null (s=0)
  L0 <- WFLike2S(s=0, data)
  ## maximum likelihood for the alternative
  x <- WFMaxiLike2S(data, maxs)
  Lmax <- x[1]
  smax <- x[2]
  warn <- x[3]
  ## likelihood ratio test
  LRT <- -2 * log(L0 / Lmax)
  ## pvalue assuming LRT follows Chi-square with 1 df
  pvalue <- -log10(1 - pchisq(LRT, 1))
  res <- matrix(c(L0, Lmax, smax, LRT, pvalue, warn), nrow=1)
  colnames(res) <- c('L0', 'Lmax', 'smax', 'LRT', '-log10pvalue', 'warn')
  rownames(res) <- ""
  return(res)
}
################################################################


