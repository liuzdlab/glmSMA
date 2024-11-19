#Input: st data, spots coordiantes
library('ape')
MoranI = function(st.data,st.coor,TotalNumberOfGenes,Distance_CUTOFF){

st.dists = as.matrix(dist(cbind(st.coor$x_coor,st.coor$y_coor)))

#distance cutoff
CUTOFF = Distance_CUTOFF
st.dists.bin = (st.dists > 0 & st.dists <= CUTOFF)

#  #s of genes
NumbersOfGenes = TotalNumberOfGenes

Moran.I.score.st = rep(0,NumbersOfGenes)
Moran.I.score.st.pval = rep(0,NumbersOfGenes)

i = 1
while(i <= NumbersOfGenes){
  if (sum(st.data[,i]) != 0){
    Moran.I.score.temp = Moran.I(st.data[,i],st.dists.bin)
    Moran.I.score.st[i] = Moran.I.score.temp$observed
    Moran.I.score.st.pval[i] = Moran.I.score.temp$p.value
  }
  i = i + 1
  #print(i)
}


return(Moran.I.score.st)
}
