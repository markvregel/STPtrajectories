#' @title alibi
#' @description alibi query
#' @param t1 t1
#' @param x1 x1
#' @param y1 y1
#' @param t2 t2
#' @param x2 x2
#' @param y2 y2
#' @param v1 v1
#' @param t3 t3
#' @param x3 x4
#' @param y3 y3
#' @param t4 t4
#' @param x4 x4
#' @param y4 y4
#' @param v2 v2
#' @return True or False for the alibi query
#' @author Sytze de Bruin and Mark ten Vregelaar
#' @importFrom polynom polynomial
#' @export
#' @keywords internal
alibi<-function(t1, x1, y1, t2, x2, y2, v1, t3, x3, y3, t4, x4, y4,v2){

  if(!(IsBeadNonEmpty(t1, x1, y1, t2, x2, y2, v1) &&
       IsBeadNonEmpty(t3, x3, y3, t4, x4, y4, v2))){
    stop("STP is not valid")
  }


  if(Case1(t1, x1, y1, t2, x2, y2, v1, t3, x3, y3, t4, x4, y4, v2)){
    return(1)
  }

  else if(Case2(t1, x1, y1, t2, x2, y2, v1, t3, x3, y3, t4, x4, y4, v2)){
    return(2)
  }
  else if(Case3(t1, x1, y1, t2, x2, y2, v1, t3, x3, y3, t4, x4, y4, v2)){
    return(3)
  }
  else{
    return(FALSE)
  }}

# ++++ lower level FUNCTIONS ++++
IsBeadNonEmpty <- function(t1,x1,y1,t2,x2,y2,v)
  t1 <= t2 && (x2-x1)^2 + (y2-y1)^2 <= v^2*(t2-t1)^2

InBead <- function(t,x,y,t1,x1,y1,t2,x2,y2,v)
  (x-x1)^2 + (y-y1)^2 <= v^2*(t-t1)^2 && (x-x2)^2 + (y-y2)^2 <=
  v^2*(t-t2)^2 && t1 <= t && t <= t2

Case1 <- function(t1,x1,y1,t2,x2,y2,v1,t3,x3,y3,t4,x4,y4,v2)
  InBead(t3,x3,y3,t1,x1,y1,t2,x2,y2,v1) ||
  InBead(t4,x4,y4,t1,x1,y1,t2,x2,y2,v1) ||
  InBead(t1,x1,y1,t3,x3,y3,t4,x4,y4,v2) ||
  InBead(t2,x2,y2,t3,x3,y3,t4,x4,y4,v2)

CheckICCondition <- function(t1,x1,y1,v1,t3,x3,y3,v2, lc)
  !((x3 - x1)^2 + (y3 - y1)^2 <= v1^2*(t3 - t1)^2 && ifelse(lc, t3 >= t1, t3 <= t1)) &&
  !((x3 - x1)^2 + (y3 - y1)^2 <= v2^2*(t3 - t1)^2 && ifelse(lc, t3 <= t1, t3 >= t1)) &&
  !(v1 == 0 && v2 == 0)

CheckHalfSpace <- function(t,x,y,t1,x1,y1,t2,x2,y2,v1,lc)
  ifelse(lc, (t1 <= t && t <= t2) &&  2*x*(x1 - x2) + x2^2 - x1^2 + 2*y*(y1 - y2) + y2^2 -
           y1^2 <= v1^2*(2*t*(t1 - t2) + t2^2 - t1^2), (t1 <= t && t <= t2) && 2*x*(x1 - x2) +
           x2^2 - x1^2 + 2*y*(y1 - y2) + y2^2 - y1^2 >= v1^2*(2*t*(t1 - t2) + t2^2 - t1^2))

Case2 <- function(t1,x1,y1,t2,x2,y2,v1,t3,x3,y3,t4,x4,y4,v2)  # this is case 3 in paper
{
  if(CheckICCondition(t1, x1, y1, v1, t3, x3, y3, v2, T))
  {
    t0 <- (t1*v1 + t3*v2 + sqrt((x1 - x3)^2 + (y1 - y3)^2))/(v1 + v2)
    x0 <- x1 + v1*(t0 - t1)*(x3 - x1)/sqrt((x1 - x3)^2 + (y1 - y3)^2)
    y0 <- y1 + v1*(t0 - t1)*(y3 - y1)/sqrt((x1 - x3)^2 + (y1 - y3)^2)
    lcs <- CheckHalfSpace(t0, x0, y0, t1, x1, y1, t2, x2, y2, v1, T) &&
      CheckHalfSpace(t0, x0, y0, t3, x3, y3, t4, x4, y4, v2, T)
  }
  else lcs <- F
  if(!lcs)
  {
    if(CheckICCondition(t2, x2, y2, v1, t4, x4, y4, v2, F))
    {
      t0 <- (t2*v1 + t4*v2 + sqrt((x2 - x4)^2 + (y2 - y4)^2))/(v1 + v2)
      x0 <- x2 + v1*(t0 - t2)*(x4 - x2)/sqrt((x2 - x4)^2 + (y2 - y4)^2)
      y0 <- y2 + v1*(t0 - t2)*(y4 - y2)/sqrt((x2 - x4)^2 + (y2 - y4)^2)
      lcs <- CheckHalfSpace(t0, x0, y0, t1, x1, y1, t2, x2, y2, v1, F) &&
        CheckHalfSpace(t0, x0, y0, t3, x3, y3, t4, x4, y4, v2, F)
    }
    else lcs <- F}
  return(lcs)
}

CoordinateChange <- function(t1,x1,y1,t2,x2,y2,tp1,xp1,yp1,tp2,
                             xp2,yp2,tp3,xp3,yp3,tp4,xp4,yp4)
{
  at1 = sqrt((x2 - x1)^2 + (y2 - y1)^2)
  a11 = x2 - x1; a12 = y2 - y1
  a21 = y1 - y2; a22 = x2 - x1
  if(y1 == y2) f <- function(t, x, y) c(t - t1, x - x1, y - y1)
  else         f <- function(t, x, y) c(at1*(t - t1), (x - x1)*a11 + a12*(y - y1),
                                        a21*(x - x1) + a22*(y - y1))
  return(c(f(tp1, xp1, yp1), f(tp2, xp2, yp2), f(tp3, xp3, yp3), f(tp4, xp4, yp4)))
}

ParamY <- function(x, t2, x2, v1, t3, x3, y3, v2)
  ifelse(y3 != 0, (((2*x*x2 - x2^2 + v1^2*t2^2)/(2*v1^2*t2) - t3)^2 *
                     v2^2 - (x - x3)^2 - v1^2 * ((2*x*x2 - x2^2 + v1^2*t2^2)/(2*v1^2*t2))^2 -
                     x^2 - y3^2)/(2*y3), sqrt((v1*(2*x*x2 - x2^2 + v1^2*t2^2)/(2*v1^2*t2))^2 - x^2))

ParamYRoot <- function(x, t2, x2, v1)
  (v1*(2*x*x2 - x2^2 + v1^2*t2^2)/(2*v1^2*t2))^2 - x^2 >= 0

ParamT <- function(x, t2, x2, v1)
  (2*x*x2 - x2^2 + v1^2*t2^2)/(2*v1^2*t2)

ComputeRoots <- function(t2, x2, v1, t3, x3, y3, v2){
  # solves the polynomial of case 3
  p1 <- polynomial(c(-x2^2 + v1^2*t2^2, 2*x2))
  a <- 2*v1^2*t2
  p2 <- polynomial(c(0,0,1))  # x-square
  p3 <- (2*y3)^2 * (p1^2 * v1^2/a^2 - p2)
  p4 <- ((p1/a - t3)^2 * v2^2 -
           (polynomial(c(-x3, 1)))^2 - y3^2 - (p1^2*v1^2/a^2 - p2))^2
  solve(p3-p4)
}

Case3 <- function(t1,x1,y1,t2,x2,y2,v1,t3,x3,y3,t4,x4,y4,v2){  # this is case 2 in paper
  TINY <- 1e-6
  if(v1*v2 == 0) return(F)
  tmp_crds <- CoordinateChange(t1, x1, y1, t2, x2, y2, t1, x1, y1, t2, x2, y2, t3, x3, y3, t4, x4, y4)
  t1s <- tmp_crds[1]; x1s <- tmp_crds[2]; y1s <- tmp_crds[3]; t2s <- tmp_crds[4]; x2s <- tmp_crds[5];
  y2s <- tmp_crds[6]; t3s <- tmp_crds[7]; x3s <- tmp_crds[8]; y3s <- tmp_crds[9]; t4s <- tmp_crds[10];
  x4s <- tmp_crds[11]; y4s <- tmp_crds[12]
  lroots <- ComputeRoots(t2s, x2s, v1, t3s, x3s, y3s, v2)
  Found <- F
  if(length(lroots)>0) for(i in 1:length(lroots))
  {
    if(abs(Im(lroots[i])) < TINY){
      Found <- ParamYRoot(Re(lroots[i]), t2s, x2s, v1) &&
        (0 <= ParamT(Re(lroots[i]), t2s, x2s, v1) &&
           ParamT(Re(lroots[i]), t2s, x2s, v1) <= t2s) && (CheckHalfSpace(
             ParamT(Re(lroots[i]), t2s, x2s, v1), Re(lroots[i]),
             ParamY(Re(lroots[i]), t2s, x2s, v1, t3s, x3s, y3s, v2),
             t3s, x3s, y3s, t4s, x4s, y4s, v2, T) ||
               ifelse(abs(y3s) == 0,
                      0 <= ParamT(Re(lroots[i]), t2s, x2s, v1) && ParamT(Re(lroots[i]), t2s, x2s, v1) <= t2s &&
                        CheckHalfSpace(ParamT(Re(lroots[i]), t2s, x2s, v1),
                                       Re(lroots[i]), -ParamY(Re(lroots[i]), t2s, x2s, v1, t3s, x3s,
                                                              y3s, v2), t3s, x3s, y3s, t4s, x4s, y4s, v2, T), F))
    }
    if(Found)
      break

  }

  if(!Found){                                              # step 2
    lroots = ComputeRoots(t2s, x2s, v1, t4s, x4s, y4s, v2)
    if(length(lroots)>0) for(i in 1:length(lroots))
    {
      if(abs(Im(lroots[i])) < TINY){
        Found <- ParamYRoot(Re(lroots[i]), t2s, x2s, v1) &&
          (0 <= ParamT(Re(lroots[i]), t2s, x2s, v1) &&
             ParamT(Re(lroots[i]), t2s, x2s, v1) <= t2s) && (CheckHalfSpace(
               ParamT(Re(lroots[i]), t2s, x2s, v1), Re(lroots[i]),
               ParamY(Re(lroots[i]), t2s, x2s, v1, t3s, x3s, y3s, v2),
               t3s, x3s, y3s, t4s, x4s, y4s, v2, F) ||
                 ifelse(y4s == 0,
                        CheckHalfSpace(ParamT(Re(lroots[i]), t2s, x2s, v1),
                                       Re(lroots[i]), -ParamY(Re(lroots[i]), t2s, x2s, v1, t3s, x3s,
                                                              y3s, v2), t3s, x3s, y3s, t4s, x4s, y4s, v2, F), F))
      }
      if(Found)
        break
    }
  }

  if(!Found){                                              # step 3
    tmp_crds <- CoordinateChange(t3, x3, y3, t4, x4, y4, t3, x3, y3, t4, x4, y4, t1, x1, y1, t2, x2, y2)
    t3s<- tmp_crds[1]; x3s<- tmp_crds[2]; y3s<- tmp_crds[3]; t4s<- tmp_crds[4]; x4s<- tmp_crds[5];
    y4s<- tmp_crds[6]; t1s<- tmp_crds[7]; x1s<- tmp_crds[8]; y1s<- tmp_crds[9]; t2s<- tmp_crds[10];
    x2s<- tmp_crds[11]; y2s <- tmp_crds[12]
    lroots <- ComputeRoots(t4s, x4s, v2, t1s, x1s, y1s, v1)
    if(length(lroots)>0) for(i in 1:length(lroots))
    {
      if(abs(Im(lroots[i])) < TINY){
        Found <- ParamYRoot(Re(lroots[i]), t4s, x4s, v2) &&
          (0 <= ParamT(Re(lroots[i]), t4s, x4s, v2) &&
             ParamT(Re(lroots[i]), t4s, x4s, v2) <= t4s) && (CheckHalfSpace(
               ParamT(Re(lroots[i]), t4s, x4s, v2), Re(lroots[i]),
               ParamY(Re(lroots[i]), t4s, x4s, v2, t1s, x1s, y1s, v1),
               t1s, x1s, y1s, t2s, x2s, y2s, v1, T) ||
                 ifelse(y1s == 0,
                        CheckHalfSpace(ParamT(Re(lroots[i]), t4s, x4s, v2),
                                       Re(lroots[i]), -ParamY(Re(lroots[i]), t4s, x4s, v2, t1s, x1s,
                                                              y1s, v1), t1s, x1s, y1s, t2s, x2s, y2s, v1, T), F))
      }
      if(Found)
        break
    }
  }

  if(!Found){                                              # step 4
    lroots <- ComputeRoots(t4s, x4s, v2, t2s, x2s, y2s, v1)
    if(length(lroots)>0) for(i in 1:length(lroots))
    {
      if(abs(Im(lroots[i])) < TINY){
        Found <- ParamYRoot(Re(lroots[i]), t4s, x4s, v2) &&
          (0 <= ParamT(Re(lroots[i]), t4s, x4s, v2) &&
             ParamT(Re(lroots[i]), t4s, x4s, v2) <= t4s) && (CheckHalfSpace(
               ParamT(Re(lroots[i]), t4s, x4s, v2), Re(lroots[i]),
               ParamY(Re(lroots[i]), t4s, x4s, v2, t2s, x2s, y2s, v1),
               t1s, x1s, y1s, t2s, x2s, y2s, v1, F) ||
                 ifelse(y2s == 0,
                        CheckHalfSpace(ParamT(Re(lroots[i]), t4s, x4s, v2),
                                       Re(lroots[i]), -ParamY(Re(lroots[i]), t4s, x4s, v2, t2s, x2s,
                                                              y2s, v1), t1s, x1s, y1s, t2s, x2s, y2s, v1, F), F))
      }
      if(Found)
        break
    }
  }
  return(Found)
}











