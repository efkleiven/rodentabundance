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

Because N is defined over all sampling sites (within a region or a block) it is assumed that all individuals that compose to pop size $N_t$ have equal access to all cameras traps. In reality they likely do not. While this may work in case of strong synchrony, on other cases this might lead to a true N being larger than the one estimated, hence some underestimation of true pop size. Having a global N does solve, however, isssues of movement of individuals from one sampling site to the next between primary occasions though. That said, CMR data suggests in Porsanger that individuals do not move between sampling sites. 

### Model of Farr et al. Cons Biol 2022

This model assumes that all sampling sites have independent abundances and that the total population size can be obtained by summing over them. The detection submodel writes

$$ Y_{m,t,k} \sim \mathcal{B}(p_{m,t,k}), \ p_{m,t,k}  = 1 - (1-\theta_{m,t,k})^{N_{t,m}}$$ 

while the biological process model writes

$$ N_{t,m} = S_{t-1,m} + G_{t-1,m} $$
$$ S_{t-1,m} \sim \mathcal{B}(N_{t-1,m},\omega_{m,t-1}) $$     
$$ G_{t-1,m} \sim \mathcal{P}(\gamma_{m,t-1}) $$

Because $N_m$ is sampling site $m$ (camera)-specific and a sum has to be taken so that $N_t = \sum_{m=1}^{M} N_{m,t}$, it is assumed that individuals do not flow from one sampling site to the next, not only within secondary occasions but also between primary ones (assumes independent dynamics per site). In reality, depending on the distance between sites considered, they might from time to time. This might lead to a true N lower than the one being estimated, hence some overestimation of true pop size. 

### Consensus model

Let's index the region by $r$--which can be a block of cameras in a given habitat or a larger spatial extent if need be

$$ Y_{m,t,k} \sim \mathcal{B}(p_{m,t,k}), \ p_{m,t,k}  = 1 - (1-\theta_{m,t,k})^{N_{t,r[m]}}$$ 

while the biological process model writes

$$ N_{t,r} = S_{t-1,r} + G_{t-1,r} $$
$$ S_{t-1,r} \sim \mathcal{B}(N_{t-1,r},\omega_t) $$     
$$ G_{t-1,r} \sim \mathcal{P}(\gamma_{t-1}) $$

we then use covariates on those rates to be able to estimate season-specific and year-specific effects, and predict for a given year x season combination. We could have structured by phase, but with 5 or 6 years of data it really makes less sense than for a very long time series. 

$$ \text{logit}(\omega_t) = \mu_{\omega} + \beta_{\omega,1,\text{year}[t]} + \beta_{\omega,2,\text{season}[t]} $$
$$ \log(\gamma_t) = \mu_{\gamma} + \beta_{\gamma,1,\text{year}[t]} + \beta_{\gamma,2,\text{season}[t]}$$

Here we assume that recruitment and survival are synchronized throughout the whole area, which in our case of rodent dynamics in Porsanger (and Håkøya) may be a reasonable assumption. Then as we move to Varanger we may want to have differing recruitment and survival rates. If there are many blocks nested within larger regions this could potentially be made as well with a (hierarchical) random effects structure. 

$\theta_{m,t,k}$ can be either constant or influenced by covariables (but not freely estimated at that scale), and covariates must be not too heavily correlated with $N_{t,r[m]}$.

### How could we introduce heterogeneity between individuals? 

(Should we want to do that). Let's drop indices and consider a simplified detection model at a given site, primary occasion, and secondary occasion. We just sketch some ideas

$$ Y \sim \mathcal{B}(p), \ p  = 1 - (1-\theta)^N $$ 

This assumes individuals have all the same $\theta$. But let's assume there is variation between thetas, for instance through a beta distribution. If $\theta \sim Beta(a,b)$ then $1-\theta \sim Beta(b,a)$. Although it is possible to compute $\mathbb{E}((1-\theta)^N)$ (moment of the Beta distribution) the expression is a product of $N$ terms rather than a function of N, which does not help us a lot. A more fruitful way to proceed might be to consider a finite mixture of e.g. two different probabilities so that $(1-p) = (1-\theta_1)^{\alpha N}(1-\theta_2)^{(1-\alpha)N}$ although this model might not be identifiable (useful to think). Maybe some colleagues already working with the Royle-Nichols model have found simpler ways to do this. 
