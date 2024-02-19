library(cmdstanr)
library(posterior)
library(tidyverse)
library(ape)
library(V.PhyloMaker) 

# generate the phylogenetic distance
splist <-read.csv("splist.csv")
tree.a <-phylo.maker(sp.list = splist, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S1")
tree.a$scenario.1$edge.length[tree.a$scenario.1$edge.length<0]<-0 
tree_cop_zero_GBOTB <- ape::cophenetic.phylo(tree.a$scenario.1)
treedata<-tree_cop_zero_GBOTB


#Specifying Analysis Data
tree_dist<-as.matrix(treedata)
colnames(tree_dist)<-sub("^X","",colnames(tree_dist))
spdata<-read.csv("spdata.csv")

#Scaling genetic distance up to 2
tree_dist<-tree_dist/max(tree_dist)*2

nsp<-nrow(tree_dist)

#Rearrange the genetic distance matrix (so that the closest ones are next to each other)
cmdres<-cmdscale(tree_dist)
tree_dist_ord<-tree_dist[order(cmdres[,1]), order(cmdres[,1])]

#spdata renumbering
replace<-data.frame(species_id=as.numeric(colnames(tree_dist_ord)),num=1:nsp)
spdata<-spdata%>%left_join(replace,by="species_id")

####build neighbor index
M<-5
distM<-tree_dist_ord
dist2nn<-function(distM,M){
  n<-nrow(distM)
  #NN_dist & NN_ind
  NN_dist<-NN_ind<-matrix(0,n-1,M)
  for(i in 2:n){
    tempdist<-distM[i,1:(i-1)]
    ntemp<-length(tempdist)
    nimp<-pmin(ntemp,M)
    NN_dist[i-1,1:nimp]<-sort(tempdist)[1:nimp]
    NN_ind[i-1,1:nimp]<-(1:(n-1))[order(tempdist)][1:nimp]
  }
  
  #NN_distM
  NN_distM<-matrix(0,n-1,M*(M-1)/2)
  for(i in 3:n){
    ind<-NN_ind[i-1,];	ind<-ind[ind!=0]
    distM_ind<-distM[ind,ind,drop=F]
    distM_lower<-distM_ind[lower.tri(distM_ind)]
    nlower<-length(distM_lower)
    NN_distM[i-1,1:nlower]<-distM_lower
  }
  res<-list(NN_ind=NN_ind,NN_dist=NN_dist,NN_distM=NN_distM)
  return(res)
}

NN.matrix<-dist2nn(distM,M)
NN_ind<-NN.matrix$NN_ind
NN_dist<-NN.matrix$NN_dist
NN_distM<-NN.matrix$NN_distM

#######data
#response variable
Y<-ifelse(spdata$RL2==1,0,1)
N<-length(Y)
#Explanatory variables
ncat_RL1<-length(unique(spdata$RL1))
contr_RL1<-contr.poly(ncat_RL1)
Xmat_RL1<-contr_RL1[spdata$RL1,]

RL1_L<-Xmat_RL1[,1]
RL1_Q<-Xmat_RL1[,2]
RL1_C<-Xmat_RL1[,3]

Temperature<-spdata$Temperature
Artificial_land<-spdata$Artificial.land
Agricultural_land<-spdata$Agricultural.land
RiverLake<-spdata$RiverLake
Seashore<-spdata$Seashore
Wasteland<-spdata$Wasteland
Precipitation<-spdata$Precipitation
Volcanic_area<-spdata$Volcanic.area
Protected_area<-spdata$Protected.area

X<-cbind(rep(1,N),RL1_L,RL1_Q,RL1_C,Temperature,Artificial_land,Agricultural_land,RiverLake,Seashore,Wasteland,Precipitation,Volcanic_area,Protected_area)
Xname<-c("RL1_L","RL1_Q","RL1_C","Temperature","Artificial_land","Agricultural_land","RiverLake","Seashore","Wasteland","Precipitation","Volcanic_area","Protected_area")
K<-ncol(X)


spnum<-spdata$num

Nsp<-length(unique(spnum))

#prior setting
#rho0<-0.3;	alpha1<-0.05	#p(rho<0.3)=0.05
#sigma0<-3;	alpha2<-0.05	#p(sigma>3)=0.05
#d<-1	#dimension=1
nu<-2
scale<-1	#scale parameter for Cauchy prior

alpha<-10
beta<-10
#When not TRUE, it is better to set alpha and beta larger
quantile(tree_dist[lower.tri(tree_dist)],0.025)<min(1/rgamma(100000,alpha,beta))

grainsize<-1

data_list<-list(N=N, M=M, K=K, Nsp=Nsp, 
                Y=Y, spnum=spnum, X=X,
                NN_ind=NN_ind, NN_dist=NN_dist, NN_distM=NN_distM,
                alpha=alpha,beta=beta,
                scale=scale,
                grainsize=grainsize)


# NNGP based models in Stan
# Code c 2018, Lu Zhang, licensed under BSD(3-clause)


stancode<-"
	functions{
		real nngp_lpdf(vector Y, vector X_beta, real sigmasq, real tausq,
                     real phi, matrix NN_dist, matrix NN_distM, array[,] int NN_ind,
                     int N, int M){

          vector[N] V;
          vector[N] YXb = Y - X_beta;
          vector[N] U = YXb;
          real kappa_p_1 = tausq / sigmasq + 1;
          int dim;
          int h;

          for (i in 2:N) {
              matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M]
              iNNdistM;
              matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M]
              iNNCholL;
              vector[ i < (M + 1) ? (i - 1) : M] iNNcorr;
              vector[ i < (M + 1) ? (i - 1) : M] v;
              row_vector[i < (M + 1) ? (i - 1) : M] v2;
              dim = (i < (M + 1))? (i - 1) : M;

              if(dim == 1){iNNdistM[1, 1] = kappa_p_1;}
              else{
                  h = 0;
                  for (j in 1:(dim - 1)){
                      for (k in (j + 1):dim){
                          h = h + 1;
                          iNNdistM[j, k] = exp(- phi * NN_distM[(i - 1), h]);
                          iNNdistM[k, j] = iNNdistM[j, k];
                      }
                  }
                  for(j in 1:dim){
                      iNNdistM[j, j] = kappa_p_1;
                  }
              }

              iNNCholL = cholesky_decompose(iNNdistM);
              iNNcorr = to_vector(exp(- phi * NN_dist[(i - 1), 1: dim]));

             v = mdivide_left_tri_low(iNNCholL, iNNcorr);

             V[i] = kappa_p_1 - dot_self(v);

             v2 = mdivide_right_tri_low(v', iNNCholL);

             U[i] = U[i] - v2 * YXb[NN_ind[(i - 1), 1:dim]];
          }
          V[1] = kappa_p_1;
          return - 0.5 * ( 1 / sigmasq * dot_product(U, (U ./ V)) +
                          sum(log(V)) + N * log(sigmasq));
		}
		
		/* integer sequence of values
		* Args: 
		*   start: starting integer
		*   end: ending integer
		* Returns: 
		*   an integer sequence from start to end
		*/ 
		array[] int sequence(int start, int end) { 
			array[end - start + 1] int seq;
			for (n in 1:num_elements(seq)) {
				seq[n] = n + start - 1;
			}
			return seq; 
		} 
		
		// compute partial sums of the log-likelihood
		real partial_log_lik_lpmf(array[] int seq, int start, int end, data array[] int Y, data matrix Xc, vector b, real Intercept, data array[] int J_1,
			data vector Z_1_0, data vector Z_1_1, data vector Z_1_2, data vector Z_1_3, data vector Z_1_4, data vector Z_1_5, data vector Z_1_6, data vector Z_1_7, data vector Z_1_8, data vector Z_1_9, data vector Z_1_10, data vector Z_1_11, data vector Z_1_12,
			vector r_1_0, vector r_1_1, vector r_1_2, vector r_1_3, vector r_1_4, vector r_1_5, vector r_1_6, vector r_1_7, vector r_1_8, vector r_1_9, vector r_1_10, vector r_1_11, vector r_1_12) {
			real ptarget = 0;
			int N = end - start + 1;
			// initialize linear predictor term
			vector[N] mu = Intercept + rep_vector(0.0, N);
			for (n in 1:N) {
				// add more terms to the linear predictor
				int nn = n + start - 1;
				mu[n] += r_1_0[J_1[nn]] * Z_1_0[nn] + r_1_1[J_1[nn]] * Z_1_1[nn] + r_1_2[J_1[nn]] * Z_1_2[nn] + r_1_3[J_1[nn]] * Z_1_3[nn] + r_1_4[J_1[nn]] * Z_1_4[nn] + r_1_5[J_1[nn]] * Z_1_5[nn] + r_1_6[J_1[nn]] * Z_1_6[nn] + r_1_7[J_1[nn]] * Z_1_7[nn] + r_1_8[J_1[nn]] * Z_1_8[nn] + r_1_9[J_1[nn]] * Z_1_9[nn] + r_1_10[J_1[nn]] * Z_1_10[nn] + r_1_11[J_1[nn]] * Z_1_11[nn] + r_1_12[J_1[nn]] * Z_1_12[nn];
			}
			ptarget += bernoulli_logit_glm_lpmf(Y[start:end] | Xc[start:end], mu, b);
			return ptarget;
		}
		
	}

	data{
		int<lower=0> N;
		int<lower=0> M;
		int<lower=0> K;
		int<lower=0> Nsp;
		array[N] int<lower=0> Y;
		array[N] int<lower=0> spnum;
		matrix[N,K] X;
		array[Nsp - 1, M] int NN_ind;
		matrix[Nsp - 1, M] NN_dist;
		matrix[Nsp - 1, ((M * (M - 1)) %/% 2)] NN_distM;
		real alpha;
		real beta;
		real scale;
		int grainsize;
	}
	
	transformed data{
		int Kc = K - 1;
		matrix[N, Kc] Xc;  // X without an intercept
		array[N] int seq = sequence(1, N);
		for (i in 2:K) {
			Xc[, i - 1] = X[, i];
		}
		vector[N] Z_1_0 = col(X,1);
		vector[N] Z_1_1 = col(X,2);
		vector[N] Z_1_2 = col(X,3);
		vector[N] Z_1_3 = col(X,4);
		vector[N] Z_1_4 = col(X,5);
		vector[N] Z_1_5 = col(X,6);
		vector[N] Z_1_6 = col(X,7);
		vector[N] Z_1_7 = col(X,8);
		vector[N] Z_1_8 = col(X,9);
		vector[N] Z_1_9 = col(X,10);
		vector[N] Z_1_10 = col(X,11);
		vector[N] Z_1_11 = col(X,12);
		vector[N] Z_1_12 = col(X,13);
		vector[Nsp] zerovec;
		for(i in 1:Nsp){
			zerovec[i] = 0.0;
		}
	}
	
	parameters{
		vector<lower=0,upper=1>[K-1] b_unif;
		real<lower=0,upper=1> Intercept_unif;
		//vector<lower=0>[3] sppar0;
		//matrix<lower=0>[3,K-1] sppar;
		real<lower=0> rho;
		real<lower=0.5,upper=1> sigma0_unif;
		real<lower=0.5,upper=1> tau0_unif;
		vector<lower=0.5,upper=1>[K-1] sigma_unif;
		vector<lower=0.5,upper=1>[K-1] tau_unif;
		
		vector[Nsp] s_1_0;
		vector[Nsp] s_1_1;
		vector[Nsp] s_1_2;
		vector[Nsp] s_1_3;
		vector[Nsp] s_1_4;
		vector[Nsp] s_1_5;
		vector[Nsp] s_1_6;
		vector[Nsp] s_1_7;
		vector[Nsp] s_1_8;
		vector[Nsp] s_1_9;
		vector[Nsp] s_1_10;
		vector[Nsp] s_1_11;
		vector[Nsp] s_1_12;
		
	}
	
	transformed parameters{
		real phi = sqrt(8.0)/rho;
		real sigma0 = tan(pi()*(sigma0_unif-1.0/2));
		real tau0 = inv_Phi(tau0_unif);
		array[K-1] real sigma;
		array[K-1] real tau;
		for(i in 1:(K-1)){
			sigma[i] = tan(pi()*(sigma_unif[i]-1.0/2));
			tau[i] = inv_Phi(tau_unif[i]);
		}
		real Intercept = scale*inv_Phi(Intercept_unif);	//scale*tan(pi()*(Intercept_unif-1/2.0))
		vector[K-1] b;
		for(i in 1:(K-1)){
			b[i] = scale*inv_Phi(b_unif[i]);	//scale*tan(pi()*(b_unif[i]-1/2.0))
		}
		vector[Nsp] r_1_0 = tau0*s_1_0;
		vector[Nsp] r_1_1 = tau[1]*s_1_1;
		vector[Nsp] r_1_2 = tau[2]*s_1_2;
		vector[Nsp] r_1_3 = tau[3]*s_1_3;
		vector[Nsp] r_1_4 = tau[4]*s_1_4;
		vector[Nsp] r_1_5 = tau[5]*s_1_5;
		vector[Nsp] r_1_6 = tau[6]*s_1_6;
		vector[Nsp] r_1_7 = tau[7]*s_1_7;
		vector[Nsp] r_1_8 = tau[8]*s_1_8;
		vector[Nsp] r_1_9 = tau[9]*s_1_9;
		vector[Nsp] r_1_10 = tau[10]*s_1_10;
		vector[Nsp] r_1_11 = tau[11]*s_1_11;
		vector[Nsp] r_1_12 = tau[12]*s_1_12;
	}
		
	model{
		s_1_0~nngp(zerovec, sigma0^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);
		s_1_1~nngp(zerovec, sigma[1]^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);
		s_1_2~nngp(zerovec, sigma[2]^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);
		s_1_3~nngp(zerovec, sigma[3]^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);
		s_1_4~nngp(zerovec, sigma[4]^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);
		s_1_5~nngp(zerovec, sigma[5]^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);
		s_1_6~nngp(zerovec, sigma[6]^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);
		s_1_7~nngp(zerovec, sigma[7]^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);
		s_1_8~nngp(zerovec, sigma[8]^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);
		s_1_9~nngp(zerovec, sigma[9]^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);
		s_1_10~nngp(zerovec, sigma[10]^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);
		s_1_11~nngp(zerovec, sigma[11]^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);
		s_1_12~nngp(zerovec, sigma[12]^2, 1, phi, NN_dist, NN_distM, NN_ind, Nsp, M);

		target += reduce_sum(partial_log_lik_lpmf, seq, grainsize, Y, Xc, b, Intercept, spnum, Z_1_0, Z_1_1, Z_1_2, Z_1_3, Z_1_4, Z_1_5, Z_1_6, Z_1_7, Z_1_8, Z_1_9, Z_1_10, Z_1_11, Z_1_12, r_1_0, r_1_1, r_1_2, r_1_3, r_1_4, r_1_5, r_1_6, r_1_7, r_1_8, r_1_9, r_1_10, r_1_11, r_1_12);
		target += inv_gamma_lpdf(rho | alpha, beta);
//		for(i in 1:(K-1)){
//			target += inv_gamma_lpdf(rho[i] | alpha, beta);
//		}
	}
	generated quantities {
		vector[N] mu = Intercept + rep_vector(0.0, N);
		vector[N] log_lik;
		for (i in 1:N) {
			// add more terms to the linear predictor
			mu[i] += r_1_0[spnum[i]] * Z_1_0[i] + r_1_1[spnum[i]] * Z_1_1[i] + r_1_2[spnum[i]] * Z_1_2[i] + r_1_3[spnum[i]] * Z_1_3[i] + r_1_4[spnum[i]] * Z_1_4[i] + r_1_5[spnum[i]] * Z_1_5[i] + r_1_6[spnum[i]] * Z_1_6[i] + r_1_7[spnum[i]] * Z_1_7[i] + r_1_8[spnum[i]] * Z_1_8[i] + r_1_9[spnum[i]] * Z_1_9[i] + r_1_10[spnum[i]] * Z_1_10[i] + r_1_11[spnum[i]] * Z_1_11[i] + r_1_12[spnum[i]] * Z_1_12[i];
			log_lik[i] = bernoulli_logit_glm_lpmf(Y[i] | to_matrix(Xc[i,]), mu[i], b);
		}
	}
"

#export stan code
filename<-"RLplant_nngp.stan"
sink(filename);cat(stancode);sink()

#compile
mod <- cmdstan_model(filename)


init_b_unif<-rep(0.5,K-1)
init_sigma0_unif<-0.55
init_tau0_unif<-0.9
init_sigma_unif<-rep(0.55,K-1)
init_tau_unif<-rep(0.9,K-1)

init_Intercept_unif<-0.1
#init_sppar0<-c(fit$draws("sppar0")[100,1,])

init<-list(list(b_unif=init_b_unif,tau0_unif=init_tau0_unif,tau_unif=init_tau_unif,sigma0_unif=init_sigma0_unif,sigma_unif=init_sigma_unif,Intercept_unif=init_Intercept_unif),
           list(b_unif=init_b_unif,tau0_unif=init_tau0_unif,tau_unif=init_tau_unif,sigma0_unif=init_sigma0_unif,sigma_unif=init_sigma_unif,Intercept_unif=init_Intercept_unif),
           list(b_unif=init_b_unif,tau0_unif=init_tau0_unif,tau_unif=init_tau_unif,sigma0_unif=init_sigma0_unif,sigma_unif=init_sigma_unif,Intercept_unif=init_Intercept_unif),
           list(b_unif=init_b_unif,tau0_unif=init_tau0_unif,tau_unif=init_tau_unif,sigma0_unif=init_sigma0_unif,sigma_unif=init_sigma_unif,Intercept_unif=init_Intercept_unif))

fit1 <- mod$sample(
  data = data_list, 
  init=init,
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  iter_warmup=3000,
  iter_sampling=10000,
  thin=10,
  max_treedepth=15,
  adapt_delta=0.99,
  refresh = 10
)

fit1$save_output_files(dir=getwd(),"RLnngp_")
fit1$save_object(file="RLnngp.RDS")
save.image("RLnngp.Rdata")



