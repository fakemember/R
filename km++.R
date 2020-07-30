x=mtcars[3:10]
x=as.matrix(x)
c1x=sample(1:dim(x)[1],1)
c1=x[c1x,]

d1=c()
for (i in 1:dim(x)[1]){d1=c(d1,sum((x[i,]-c1)^2))}
p1=d1/sum(d1)
c2x=sample(1:dim(x)[1],size = 1,prob = d1)
c2=x[c2x,]



d2=c()
for (i in 1:dim(x)[1]){d2=c(d2,sum((x[i,]-c2)^2))}
p2=d2/sum(d2)
c3x=sample(1:dim(x)[1],size = 1,prob = pmax(d1,d2))

c3=x[c3x,]


s1=kmeans(x,centers = x[c(c1x,c2x,c3x),],nstart = 20)
str(s1)
kmeans(x,3,nstart = 20)






abalone=read.csv('abalone.csv',header = F)
x=abalone[,2:9]
x=scale(x)


kmeans(x,3)


c1x=sample(1:dim(x)[1],1)
c1=x[c1x,]

d1=c()
for (i in 1:dim(x)[1]){d1=c(d1,sum((x[i,]-c1)^2))}
p1=d1/sum(d1)
c2x=sample(1:dim(x)[1],size = 1,prob = p1)
c2=x[c2x,]



d2=c()
for (i in 1:dim(x)[1]){d2=c(d2,sum((x[i,]-c2)^2))}
p2=d2/sum(d2)
c3x=sample(1:dim(x)[1],size = 1,prob = p2)

c3=x[c3x,]
kmeans(x,centers = x[c(1:3),])





x=matrix(rnorm(200),100,2)
xmean=matrix(rnorm(8,sd = 4),4,2)
whic=sample(1:4,100,replace = T)
x=x+xmean[whic]
plot(x,col=whic)


k1=kmeans(x,4,nstart = 25)
k1
plot(x,col=k1$cluster,cex=2,pch=1,lwd=2)
points(x,col=whic,pch=19)

