### internal function
stable_gap_decideK <- function(indata,SD,maxK,Gapplot){

    set.seed(SD)
    b <- clusGap(indata,kmeans,maxK+20)
    GapK <- b$Tab[,3]

    Knum = 0
    GapCUTOFF <- mean(GapK[maxK:(maxK+20)]) - 2*sd(GapK[maxK:(maxK+20)])
    for(i in seq(maxK-1)){
        if(GapK[i] > GapCUTOFF & GapK[i] > GapK[i+1]){
            Knum = i
            break}
    }
    pdf(file=Gapplot)
    plot(b)
    abline(v=Knum,lwd=2,col="blue")
    abline(h=GapCUTOFF,lwd=2,col="darkgreen")
    legend("bottomright",legend=c("selected K value","cutoff for stable gap score"),col=c("blue","darkgreen"),bty="n",lwd=2)
    dev.off()

    return(Knum)
}
maxSE_decideK <- function(indata,SD,maxK,Gapplot){
    set.seed(SD)
    b <- clusGap(indata,kmeans,maxK)
    Knum <- maxSE(b$Tab[,3],b$Tab[,4],method="Tibs2001SEmax",SE.factor=0.01)
    pdf(file=Gapplot)
    plot(b)
    abline(v=Knum,lwd=2,col="blue")
    legend("topleft",legend=c("selected K value"),col=c("blue"),bty="n",lwd=2)
    dev.off()
    return(Knum)
}


givenK_kmeans <- function(indata,SD,Knum,clusterplot){
    set.seed(SD)
    pdf(file=clusterplot)
    km <- kmeans(indata,Knum)
    plot(indata,pch=16,xlab="t-SNE 1",ylab="t-SNE 2",main=paste("k-means (k=", Knum, ")",sep=""))
    rain <- rainbow(length(km$size))
    for( i in 1:length(km$size)){
	    points(indata[which(km$cluster==i),],col=rain[i],pch=16)	
    }
    text(km$centers,labels=seq(nrow(km$centers)))
    dev.off()
    cluster_result <- cbind(indata,km$cluster)
    return(cluster_result)
}


givenE_dbscan <- function(indata,EPS,clusterplot){
    pdf(file=clusterplot)
    ds <- dbscan(indata,eps=EPS)
    kmsize <- sort(cluster_result[,3],decreasing=T)[1]
    plot(indata,pch=16,xlab="t-SNE 1",ylab="t-SNE 2",main=paste("DBSCAN (eps=", EPS, ")",sep=""))
    rain <- rainbow(kmsize)
    for(i in 1:kmsize){
        points(cluster_result[which(cluster_result[, 3]==i),1:2],col=rain[i],pch=16)
    }
    dev.off()
    cluster_result <- cbind(indata,ds$cluster)
    return(cluster_result)
}

getTPM <- function(indata){
  return(log10(indata*1e4/sum(indata) + 1))
}
getcoverGnum <- function(indata){
  return(length(which(indata > 0)))
}


selct_high_var_gene <- function(inputdata,highvarZ){
Gave <- apply(inputdata,1,mean)
Gvar <- apply(inputdata,1,var)
Gdp <- Gvar/Gave
Gname <- rownames(inputdata)

each <- floor(length(Gname)/20)

g1 <- Gname[order(Gave)][(1+each*0):(each*1)]
g2 <- Gname[order(Gave)[(1+each*1):(each*2)]]
g3 <- Gname[order(Gave)[(1+each*2):(each*3)]]
g4 <- Gname[order(Gave)[(1+each*3):(each*4)]]
g5 <- Gname[order(Gave)[(1+each*4):(each*5)]]
g6 <- Gname[order(Gave)[(1+each*5):(each*6)]]
g7 <- Gname[order(Gave)[(1+each*6):(each*7)]]
g8 <- Gname[order(Gave)[(1+each*7):(each*8)]]
g9 <- Gname[order(Gave)[(1+each*8):(each*9)]]
g10 <- Gname[order(Gave)[(1+each*9):(each*10)]]
g11 <- Gname[order(Gave)[(1+each*10):(each*11)]]
g12 <- Gname[order(Gave)[(1+each*11):(each*12)]]
g13 <- Gname[order(Gave)[(1+each*12):(each*13)]]
g14 <- Gname[order(Gave)[(1+each*13):(each*14)]]
g15 <- Gname[order(Gave)[(1+each*14):(each*15)]]
g16 <- Gname[order(Gave)[(1+each*15):(each*16)]]
g17 <- Gname[order(Gave)[(1+each*16):(each*17)]]
g18 <- Gname[order(Gave)[(1+each*17):(each*18)]]
g19 <- Gname[order(Gave)[(1+each*18):(each*19)]]
g20 <- Gname[order(Gave)[(1+each*19):(length(Gname))]]

highvargene <- c(g2[which((scale(Gdp[g2]))>highvarZ)],g3[which((scale(Gdp[g3]))>highvarZ)],g4[which((scale(Gdp[g4]))>highvarZ)],g5[which((scale(Gdp[g5]))>highvarZ)],g6[which((scale(Gdp[g6]))>highvarZ)],g7[which((scale(Gdp[g7]))>highvarZ)],g8[which((scale(Gdp[g8]))>highvarZ)],g9[which((scale(Gdp[g9]))>highvarZ)],g10[which((scale(Gdp[g10]))>highvarZ)],g11[which((scale(Gdp[g11]))>highvarZ)],g12[which((scale(Gdp[g12]))>highvarZ)],g13[which((scale(Gdp[g13]))>highvarZ)],g14[which((scale(Gdp[g14]))>highvarZ)],g15[which((scale(Gdp[g15]))>highvarZ)],g16[which((scale(Gdp[g16]))>highvarZ)],g17[which((scale(Gdp[g17]))>highvarZ)],g18[which((scale(Gdp[g18]))>highvarZ)],g19[which((scale(Gdp[g19]))>highvarZ)],g20[which((scale(Gdp[g20]))>highvarZ)])
return(highvargene)
}

ref_tsne <-
function(X, initial_config = NULL, k=2, initial_dims=30, perplexity=30, max_iter = 1000, min_cost=0, epoch_callback=NULL,whiten=TRUE, epoch=100 ){
	#### original code from http://lvdmaaten.github.io/tsne/
	#### L.J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE. 2008
	if (class(X) == 'dist') { 
		n = attr(X,'Size')
		}
	else 	{
		X = as.matrix(X)
		X = X - min(X)
		X = X/max(X)
		initial_dims = min(initial_dims,ncol(X))
		if (whiten) X<-.whiten(as.matrix(X),n.comp=initial_dims)
		n = nrow(X)
	}

	momentum = .5
	final_momentum = .8
	mom_switch_iter = 250

	epsilon = 500
	min_gain = .01
	initial_P_gain = 4

	eps = 2^(-52) # typical machine precision

	if (!is.null(initial_config) && is.matrix(initial_config)) { 		
		if (nrow(initial_config) != n | ncol(initial_config) != k){
			stop('initial_config argument does not match necessary configuration for X')
		}
		ydata = initial_config
		initial_P_gain = 1
		
	} else {
		ydata = matrix(rnorm(k * n),n)
	}
	
	P = .x2p(X,perplexity, 1e-5)$P
	# P[is.nan(P)]<-eps
	P = .5 * (P + t(P))

	P[P < eps]<-eps
	P = P/sum(P)
	

	
	P = P * initial_P_gain
	grads =  matrix(0,nrow(ydata),ncol(ydata))
	incs =  matrix(0,nrow(ydata),ncol(ydata))
	gains = matrix(1,nrow(ydata),ncol(ydata))

	
	for (iter in 1:max_iter){
		if (iter %% epoch == 0) { # epoch
			cost =  sum(apply(P * log((P+eps)/(Q+eps)),1,sum))
			message("Epoch: Iteration #",iter," error is: ",cost)
			if (cost < min_cost) break
			if (!is.null(epoch_callback)) epoch_callback(ydata)

		}


		sum_ydata = apply(ydata^2, 1, sum)
		num =  1/(1 + sum_ydata +    sweep(-2 * ydata %*% t(ydata),2, -t(sum_ydata))) 
		diag(num)=0
		Q = num / sum(num)
		if (any(is.nan(num))) message ('NaN in grad. descent')
		Q[Q < eps] = eps
		stiffnesses = 4 * (P-Q) * num
		for (i in 1:n){
			grads[i,] = apply(sweep(-ydata, 2, -ydata[i,]) * stiffnesses[,i],2,sum)
		}
		
		gains = (gains + .2) * abs(sign(grads) != sign(incs)) + gains * .8 * abs(sign(grads) == sign(incs))		
		gains[gains < min_gain] = min_gain

		incs = momentum * incs - epsilon * (gains * grads)
		ydata = ydata + incs
		ydata = sweep(ydata,2,apply(ydata,2,mean))
		if (iter == mom_switch_iter) momentum = final_momentum
		
		if (iter == 100 && is.null(initial_config)) P = P/4
		

	
		
	}
	ydata
}

.Hbeta <-
function(D, beta){
	P = exp(-D * beta)
	sumP = sum(P)
	if (sumP == 0){
		H = 0
		P = D * 0
	} else {
		H = log(sumP) + beta * sum(D %*% P) /sumP
		P = P/sumP
	}
	r = {}
	r$H = H
	r$P = P
	r
}

.x2p <-
function(X,perplexity = 15,tol = 1e-5){
	if (class(X) == 'dist') {
		D = X
		n = attr(D,'Size')
	} else{
		D = dist(X)
		n = attr(D,'Size')
	}

	D = as.matrix(D)
	P = matrix(0, n, n )		
	beta = rep(1, n)
	logU = log(perplexity)
	
	for (i in 1:n){
		betamin = -Inf
		betamax = Inf
		Di = D[i, -i]
		hbeta = .Hbeta(Di, beta[i])
		H = hbeta$H; 
		thisP = hbeta$P
		Hdiff = H - logU;
		tries = 0;

		while(abs(Hdiff) > tol && tries < 50){
			if (Hdiff > 0){
				betamin = beta[i]
				if (is.infinite(betamax)) beta[i] = beta[i] * 2
				else beta[i] = (beta[i] + betamax)/2
			} else{
				betamax = beta[i]
				if (is.infinite(betamin))  beta[i] = beta[i]/ 2
				else beta[i] = ( beta[i] + betamin) / 2
			}
			
			hbeta = .Hbeta(Di, beta[i])
			H = hbeta$H
			thisP = hbeta$P
			Hdiff = H - logU
			tries = tries + 1
		}	
			P[i,-i]  = thisP	
	}	
	
	r = {}
	r$P = P
	r$beta = beta
	sigma = sqrt(1/beta)
	
	message('sigma summary: ', paste(names(summary(sigma)),':',summary(sigma),'|',collapse=''))

	r 
}

.whiten <-
function(X, row.norm=FALSE, verbose=FALSE, n.comp=ncol(X))
{  
	n.comp; # forces an eval/save of n.comp
	if (verbose) message("Centering")
   n = nrow(X)
	p = ncol(X)
	X <- scale(X, scale = FALSE)
   X <- if (row.norm) 
       t(scale(X, scale = row.norm))
   else t(X)

   if (verbose) message("Whitening")
   V <- X %*% t(X)/n
   s <- La.svd(V)
   D <- diag(c(1/sqrt(s$d)))
   K <- D %*% t(s$u)
   K <- matrix(K[1:n.comp, ], n.comp, p)
   X = t(K %*% X)
	X
}


### data pre-process
# read parameter

a<-commandArgs(T)

inmatrix <- a[1]
outname <- a[2]
hvZ <- as.numeric(a[3])
selectPCcutoff <- as.numeric(a[4])
RDnumber <- as.numeric(a[5])
maxKnum <- as.numeric(a[6])
pcatableY <- as.numeric(a[7])
cortableY <- as.numeric(a[8])
clustering_method <- as.numeric(a[9])
custom_k <- as.numeric(a[10])
custom_d <- as.numeric(a[11])

Rdata <- read.table(inmatrix,row.names=1,header=T)

### transform to TPM, transcript per million reads(log10 scale)
Ndata <- apply(Rdata,2,getTPM)

traindata <- Ndata

### select high variance gene
highvargene <- selct_high_var_gene(traindata,hvZ)


### normalize data and conduct PCA
indata <- traindata[highvargene , ]
PCAdata <- t(indata)
cmean <- apply(PCAdata,2,mean)
csd <- apply(PCAdata,2,sd)
scPCAdata <- t((t(PCAdata)-cmean)/csd)
PCAresult <- prcomp(scPCAdata)


if(pcatableY == 1){
    pcatable <- as.matrix(scPCAdata) %*% as.matrix(PCAresult[2]$rotation)[,1:2]
    write.table(pcatable,file=paste(outname,'_pctable.txt',sep=""),row.names=T,col.names=T,sep="\t",quote=F)
}
if(cortableY == 1){
    cortable <- cor(traindata)
    write.table(cortable,file=paste(outname,'_correlation_table.txt',sep=""),row.names=T,col.names=T,sep="\t",quote=F)
}

selectPC <- max(10,length(which(summary(PCAresult)[[6]][3,]< selectPCcutoff)))
compdata <- as.matrix(scPCAdata) %*% as.matrix(PCAresult[2]$rotation)[,1:selectPC]
set.seed(RDnumber)
tsne_result <- ref_tsne(compdata,max_iter=200)
#write.table(tsne_result,file="a.txt",row.names=F,col.names=F,sep="\t",quote=F)
#tsne_result <- read.table("a.txt")

gap_plot <- paste(outname,'_GapStat.pdf',sep="")
cluster_plot <- paste(outname,'_cluster.pdf',sep="")
if (clustering_method == 4){
    ### if user choose dbscan
    library(fpc)
    final_result <- givenE_dbscan(tsne_result,custom_d,cluster_plot)
}else{
    library(cluster)
    ### if user choose kmeans
    ## different method to decide k
    if(clustering_method == 1){
        
        KNUM <- stable_gap_decideK(tsne_result,RDnumber,maxKnum,gap_plot)
    }else if (clustering_method == 2){
        KNUM <- maxSE_decideK(tsne_result,RDnumber,maxKnum,gap_plot)     
    }else{
        KNUM <- custom_k
    }

    final_result <- givenK_kmeans(tsne_result,RDnumber,KNUM,cluster_plot)
}
row.names(final_result) <- row.names(PCAdata)
if (clustering_method == 4){
    colnames(final_result) <- c("t-SNE_d1","t-SNE_d2","dbscan_cluster")
}else{
    colnames(final_result) <- c("t-SNE_d1","t-SNE_d2","kmeans_cluster")
}
write.table(final_result,file=paste(outname,'_cluster.txt',sep=""),row.names=T,col.names=T,sep="\t",quote=F)