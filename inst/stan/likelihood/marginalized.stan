for(j in 1:sum_D){ Ky += KX[j];}
target += multi_normal_lpdf(y | mu, Ky);
