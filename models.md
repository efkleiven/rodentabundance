## Models used so far

### Notations

Following previous work by Simon and Eivind, we consider 
- $i$ the species
- $t$ the primary occasion 
- $k$ the secondary occasion
- $m$ the site (or camera)
- $r$ the region. Could also be a block, the spatial unit containing multiple sampling sites.

Please mind that Farr et al. Cons Biol 2022 used $j$ for site and $k$ for repeat visits (secondary occasions). While Kleiven et al. JABES 2023 used $j$ for secondary occasions and $k$ for sites, and $b$ for blocks. 
The equivalent of blocks in that code might be regions. 


### Initial model used

For the moment we consider one species so we drop index $i$ and one region / pop so we drop index $r$ (the code considers $r \in \{1,..., R \}$ though. The detection submodel writes

$$ Y_{m,t,k} \sim \mathcal{B}(p_{m,t,k}), \ p_{m,t,k}  = 1 - (1-\theta_{m,t})^{N_t}$$ 

while the biological process model writes

$$ N_{t} = S_{t-1} + G_{t-1}$$  
$$ S_{t-1} \sim \mathcal{B}(N_{t-1},\omega)$$
$$ G_{t-1} \sim \mathcal{P}(\gamma_{t-1})$$

Because N is defined over all sampling sites (within a region or block) it is assumed that all individuals have equal access to all cameras traps. In reality they likely do not. 
This might lead to the true N being larger than the one estimated here, a little underestimation of true pop size. It does solve the problem that individuals might move from one sampling site to the next between primary occasions though.

### Model of Farr et al. Cons Biol 2022

This model assumes that all sampling sites have independent abundances and that the total population size can be obtained by summing over them. The detection submodel writes

$$ Y_{m,t,k} \sim \mathcal{B}(p_{m,t,k}), \ p_{m,t,k}  = 1 - (1-\theta_{m,t,k})^{N_{t,m}}$$ 

while the biological process model writes

$$ N_{t,m} = S_{t-1,m} + G_{t-1,m} $$
$$ S_{t-1,m} \sim \mathcal{B}(N_{t-1,m},\omega_{m,t-1}) $$     
$$ G_{t-1,m} \sim \mathcal{P}(\gamma_{m,t-1}) $$

Because $N_m$ is sampling site $m$ (camera)-specific and a sum has to be taken so that $N_t = \sum_{m=1}^{M} N_{m,t}$, it is assumed that individuals do not flow from one sampling site to the next, within a secondary occasion. In reality they might from time to time. This might lead to a true N lower than the one being estimated, 
hence some overestimation of true pop size. They might also move from site to site, from one primary occasion to the next, which the biological process model does not take into account (assumes independent dynamics per site). 

### Consensus model

Let's index the region by $r$--which can be a block of cameras in a given habitat or a larger spatial extent if need be

$$ Y_{m,t,k} \sim \mathcal{B}(p_{m,t,k}), \ p_{m,t,k}  = 1 - (1-\theta_{m,t,k})^{N_{t,r[m]}}$$ 

while the biological process model writes

$$ N_{t,r} = S_{t-1,r} + G_{t-1,r} $$
$$ S_{t-1,r} \sim \mathcal{B}(N_{t-1,r},\omega) $$     
$$ G_{t-1,r} \sim \mathcal{P}(\gamma_{t-1}) $$

Here we assume that recruitment and survival are synchronized throughout the whole area, which in our case of rodent dynamics may be a reasonable assumption. 

### How could we introduce heterogeneity between individuals? 

(Should we want to do that). Let's drop indices and consider a simplified detection model at a given site, primary occasion, and secondary occasion. We just sketch some ideas

$$ Y \sim \mathcal{B}(p), \ p  = 1 - (1-\theta)^N $$ 

This assumes individuals have all the same $\theta$. But let's assume there is variation between thetas, for instance through a beta distribution. If $\theta \sim Beta(a,b)$ then $1-\theta \sim Beta(b,a)$. Although it is possible to compute $\mathbb{E}((1-\theta)^N)$ (moment of the Beta distribution) the expression is a product of $N$ terms rather than a function of N, which does not help us a lot. A more fruitful way to proceed might be to consider a finite mixture of e.g. two different probabilities so that $(1-p) = (1-\theta_1)^{\alphaN}(1-\theta_2)^{(1-\alpha)N}$ although this model might not be identifiable (useful to think). Maybe some colleagues already working with the Royle-Nichols model have found simpler ways to do this. 
