# dataset = generateDataset()
# X = dataset$X
# real_cluster = dataset$col
# plot.curve.dataset(X,col=real_cluster,ratio=0.3)
#
# # Haar wavelet
# result = funHDDCwavelet(X,K=3,minD=1,maxD=1:4,wavelet.family="DaubExPhase",wavelet.filter.number=2,viz=F,minIter=20)
# clusters = apply(result$tm,1,which.max)
# adjustedRandIndex(clusters,real_cluster)
# plot.curve.dataset(X,col=clusters,ratio=0.3)
# plot.mean.curve.dataset(X,col=clusters)
