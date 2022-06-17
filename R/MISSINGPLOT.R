#' Missing Plot in Randomized Block Design(RBD)
#' @export
#' @param m a matrix containing values in a RBD where row of the matrix denotes the treatments and the column of the matrix denotes  block. In this matrix, we will replace the missing value with 0.
#' @param r the index no. of row/Treatment containing the missing value.
#' @param c the index no. of column/Block- containing the missing value.
#' @author  Arnab Roy , Debarghya Baul.
#' @importFrom stats aov qf
#' @description This function analyses RBD when there is one missing observation.
#' @details In RBD setup , if there is one missing observation we can use this function to estimate the missing observation along with Sum of Squares for testing the differential effect of the treatments. Here we estimate the missing observation twice by minimizing the SSE of the design.
#' @return x.hat : the least sqaure estimate of the missing observation.
#' @return  SSE.x.hat : Sum of Squares of Error of x.hat.
#' @return x.double.hat : the least square estimate of the missing observation under the null hypothesis , H0.
#' @return SSE.x.double.hat : Sum of Squares of Error of x.double.hat.
#' @return F.stat : Observed value of the Test Statistic.
#' @return F.crit.value : Critical value of the Test Statistic.
#' @examples  p=matrix(c(12,15,16,18,16,21,0,27,29,30,35,36),nrow=4,ncol=3,byrow=TRUE )
#' @examples  Missing.RBD(p,3,1)




Missing.RBD=function(m,r,c){
  B=apply(m,2,sum)[c]
  T=apply(m,1,sum)[r]

  G=sum(m)

  b=ncol(m)
  v=nrow(m)
  x.hat=(b*B+v*T-G)/((b-1)*(v-1))
  m[r,c]=x.hat

  y=as.vector(t(m))



  bl=as.factor(rep(1:b,times=v))
  trt=as.factor(rep(1:v,each=b))
  s=summary(aov(y~bl+trt))
  s1=s[[1]][2]$`Sum Sq`
  SSE.x.hat=s1[3]
  x.double.hat=B/(v-1)

  m[r,c]=x.double.hat
  y1=as.vector(t(m))

  s2=summary(aov(y1~bl))
  s3=s2[[1]][2]$`Sum Sq`
  SSE.x.double.hat=s3[2]

  num=(SSE.x.double.hat-SSE.x.hat)/(v-1)
  dem=SSE.x.hat/((b-1)*(v-1)-1)
  Fstat=num/dem

  F.crit.value=qf(0.95,v-1,((b-1)*(v-1)-1))
  list(x.hat=x.hat,SSE.x.hat=SSE.x.hat,
       x.double.hat=x.double.hat,SSE.x.double.hat=SSE.x.double.hat,
       F.stat=Fstat,F.crit.value=F.crit.value)
}


