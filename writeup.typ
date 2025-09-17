#set text(font: "New Computer Modern", size: 12pt)
#set page(margin: 
            ( left: 2.5cm, 
              right: 2.5cm, 
              top: 2.5cm, 
              bottom: 2.5cm),
          numbering: "1",
          number-align: center)

#set heading(numbering: "1.1.1.") 

#show heading: it => [
  #v(0.5em)
  #set align(left)
  #set text(13pt, weight: "regular", font: "New Computer Modern")

  #block([
    // Only display number if numbering is enabled for this heading
    #(if it.numbering != none { counter(heading).display(it.numbering) + " " })
    #smallcaps(it.body)
  ])
  #v(0.5em)
]

#let appendix(body) = {
  set heading(numbering: "A", supplement: [Appendix])
  counter(heading).update(0)
  body
}

#set par( justify: true,
          leading: 1.4em,
          spacing: 2.6em)

/*
#pagebreak()
#pagebreak()
*/

#set heading(numbering: "1.")
#outline()

#pagebreak()
= Introduction

Proteins are biomolecules that are fundamental to nearly all biological processes. Their diverse roles include transporting nutrients, catalyzing chemical reactions, providing structural support, and more. The function of a protein is determined by the composition of its primary sequence created from amino acids. A single protein sequence can vary greatly in length and order of its amino acids, leading to a very large number of possible configurations. 

The size of the protein sequence space makes exhaustive experimental exploration infeasible. However, the rapid growth of biological databases, driven by advances in sequencing technologies, has transformed biology into a data-rich discipline. Large repositories such as UniProt, Pfam, and the Protein Data Bank store millions of sequences and structures, providing crucial resources for computational approaches@weigt2020.

This wealth of data has enabled the development of advanced statistical and machine learning models capable of simulating protein sequence evolution, prediciting structural conformations, and generating novel sequences with desired properties. In addition, breakthroughs in protein structure prediction--most notably the Nobel Prize winning AlphaFold@alphafold\-- show that computational methods can rival experimental accuracy, in a much cheaper manner.

In this work, we will implement an efficient autoregressive network for protein sequence generation leveraging the Direct Coupling Analysis (DCA) in Python. In Section 2, the biological background will be laid out. Section 3 will introduce the mathematical foundations. Section 4 will detail previous iterations of similar models, while Section 5 will explore our implementation.

= Biological Background

== Proteins and their structure 
Proteins are essential biological molecules responsible for a wide range of functions in living organisms. Despite their functional diversity, all proteins are polymers of the same set of standard building blocks: the 20 canonical amino acids arranged in different assortments. Amino acids share a common core structure consisting of a central carbon atom bonded to a hydrogen atom, an amino group, a carboxyl group, and a variable side chain. This side chain is the defining feature of each amino acid, giving rise to differences in size, shape, chemical reactivity, and polarity. Proteins are formed by peptide bonds of distinct amino acids, where the term polypeptides comes from. In this process, certain atoms are lost, thus, within a protein sequence, individual amino acids are typically referred to as residues. Generally, protein sequences are made up of between 50 and 2000 amino acids. The ordering of the amino acids dictates how a protein folds into its three-dimensional structure, known as its conformation. Although each conformation is unique, two common folding patterns occur in many proteins: the $alpha$ helix and the $beta$ sheet @alberts2002. 

Each amino acid is chemically distinct and can occur at any position in a protein chain, giving rise to $20^n$ possible polypeptide sequences of length $n$. For a typical protein of 300 amino acids, the number of possible sequences is astronomically large ($approx 2 dot 10^390$). However, only a small fraction of these sequences are capable of folding into a stable three-dimensional conformation. 

== Protein families and evolutionary information
Proteins do not evolve in isolation; they often belong to protein families, groups of proteins that share a common evolutionary origin, therefore exhibit similar sequence features and functional properties @ebi_protein_families. Over evolutionary timescales, mutations accumulate in these families: some are beneficial, altering the protein activity in ways that give rise to new functions, while many others are neutral and have no effect on stability or activity. Harmful mutations, by contrast, disrupt protein folding or its function. These changes are  eliminated through natural selection. The result is a collection of homologous proteins that retain overall structural and functional characteristics, but also display sequence variability that encodes the evolutionary history of the family @alberts2002. 

The study of evolutionary history and relationships among biological entities is referred to as phylogenetics. Protein families provide the domain for phylogenetic analysis, as examining families provides insight that cannot be obtained from a single sequence@EMBL-EBI_Phylogenetics. From homologous proteins, we can detect important correlated mutations between amino acid positions which represents constraints enforced to maintain protein integrity by evolution. These statistical patterns are crucial for computational approaches as they attempt to learn them to predict three-dimensional structure and understand protein function.

== Multiple sequence alignments
To extract important statistical information from protein families, the homologous sequences need to be organized in a systematic way. This is done through a multiple sequence alignment (MSA), in which three or more sequences are arranged so that homologous sites are placed in the same column @WILTGEN201938. Alignment gaps are introduced throughout the MSA to maximise the positional correspondence of sequences and enable them to have the same fixed length. 

Aligning sequences in this manner, patterns of conservation and variation are revealed across the family. Conserved positions, those which are unchanged in multiple sequences, represent sites that are critical for maintaining structure or function, while variable positions indicate sites that can tolerate mutations without disruption to the protein. Beyond conservation, MSAs also capture covariation: pairs of positions that mutate in a correlated way across sequences. These covariation signals reflect couplings, where a mutation at one site requires compensation at another to maintain protein integrity @biom14121531.


== Protein contacts
One of the earliest challenges for researchers in computational biology--historically encouraged by the CASP (Critical Assessment of Structure Prediction) competition--has been to understand how the linear sequence of amino acids folds into its conformation. Researchers explored spatially close pairs of residues in a protein's folded three-dimensional structure, called contacts. As the structures can be represented in Cartesian coordinates $(x,y,z)$,  contacts are defined using distance thresholds. Two residues are generally considered to be in contact if the distance between the selected atoms is below a set threshold, usually 8 Ångströms (Å). 

Residues that are far apart in the linear sequence, may be close together in the 3D structure. The contact distances are further categorized into short, medium, and long range predictions. Most computational approaches evaluate the long range contacts separately, as they are the most important for accurate predictions and unsurprisingly, the hardest to predict @adhikari2016. 

== The problem of inference

A central challenge in protein sequence analysis is disentangling the true structural and functional constraints of protein families from spurious correlations introduced by other factors. Correlations can arise indirectly, for example, if residue A correlates with B, and B with C, then A and C may appear correlated without a direct interaction. Moreover, shared evolutionary history (phylogeny) and biases in sequence databases can further obscure the true couplings that underlie protein folding and function @dietler2023.

= Mathematical Foundations

== Proteins as statistical systems
Protein sequences can be thought of as random variables produced by a certain distribution.  Each sequence of length $L$ can be written as:
#set math.equation(numbering: "(1)", supplement: "Eq.")
$
bold(sigma) = (sigma_1, sigma_2, ..., sigma_L), quad sigma_i in cal(A),
$ where $cal(A)$ is the alphabet of size $q = 21$ (20 amino acids and the alignment gap) and $sigma_i$ are the residue sites. We can organize these sequences into multiple sequence alignments, a table ${bold(sigma)^((m))}^M_(m=1)$ of $M$ empirical samples. These samples are aligned to have a common length $L$. Each row in the MSA represents a protein, and each column a position, or site, in the sequence. From the alignments, we can define single and pairwise frequency counts for the columns. The single-site frequency for MSA column $i$ can be computed as:
$
f_i (A) = 1/M sum^M_(m=1) delta(sigma_i^((m)), A), quad A in cal(A),quad  
"where" delta "is the Kronecker delta" #footnote[$delta(x, y) := cases( 0 "if" x != y, 1 "if" x=y) $]
$
while the pairwise frequency of MSA columns $i, j$ is computed as:
$
  f_(i j) (A, B) = 1/M sum^M_(m=1) delta(sigma_i^((m)), A) delta(sigma_j^((m)), B)
$

The empirical frequencies will serve as constraints for the distribution we want to infer.

== Maximum Entropy Principle
To find the probability distribution $P(bold(sigma))$ that will satisfy constraints, in other words reproduce the empirical marginals $f_i (A)$ and $f_(i j) (A, B)$, we can use the maximum entropy principle@Jaynes1957a @Jaynes1957b. The first step to set up the MEP is extracting information from the system. Usually, this information is given in the form of averages of functions $angle.l f_k (x) angle.r$. For example, in a physical system, one could compute average magnetization or energy of the observed system. We also need to define a probability of occupancy of states, $p(x)$, which runs over all the possible states. This distribution has the usual property of mapping each state to a value within 0 and 1 and adding up to 1 when considering all states @MITinfent. Our uncertainty on the system is expressed quantitatively through Shannon entropy $S$@Shannon1948:
$
  S = - sum_x p(x) ln p(x) 
$ <entropy>
Generally, the distribution which maximizes entropy is the uniform distribution. However, since constraints affect the probabilities of states, the uniform is not suitable. The constraints on our system are:
$
  sum_x p(x) = 1, quad "and" quad sum_x p(x) f_k (x) = angle.l f_k angle.r quad (k = 1, ..., m).
$ <constraints>

Selecting Shannon entropy as our measure of information, we can introduce Lagrange Multipliers $lambda_0, lambda_1, ... lambda_m$ to maximize @entropy subject to the constraints @constraints, yielding:
$
  p(x) = exp(-lambda_0 - sum_(k=1)^m lambda_k f_k (x)).
$
Define the partition function $Z$ as:
$
  Z(lambda) = sum_x exp(-sum_(k=1)^m lambda_k f_k (x)).
$
With normalization $lambda_0 = ln Z$, we can write the moments as:
$
  angle.l f_k angle.r = - partial / (partial lambda_k) ln Z(lambda)
$

The entropy of the distribution  then reduces to#footnote[The full derivation can be found in @app1.]:
$
  S_max = lambda_0 + sum_(k=1)^m lambda_k angle.l f_k angle.r
$
The maximum-entropy (MaxEnt) distribution for this system is exactly the Potts model.

== Connection to statistical mechanics
The Potts model is rooted in statistical mechanics, where it was introduced as a generalization of the Ising spin model. Both models were developed to study ferromagnetism and phase transitions.

=== Ising Model
The Ising model @ising1925 describes a system of spins $sigma_i$ arranged on a lattice. Each spin can take one of two possible values:
$
  sigma_i in {+1, -1}
$
The Hamiltonian function, which represents the energy, of a spin configuration ${sigma_i}$ is 
$
  H_("Ising")({sigma}) = -J sum_(angle.l i, j angle.r) sigma_i sigma_j - h sum_i sigma_i,
$
where $J$ is the coupling constant calculated over nearest-neighbor pairs (represented by $angle.l i, j angle.r$) that defines the strength between paired interactions. The $h$ term is the external magnetic field acting on the spins.

At thermal equilibrium, the probability of a configuration is given by the Boltzmann distribution:
$
  P({sigma}) = 1 / Z exp(-beta H_"Ising" ({sigma})),
$
with $beta = 1 \/ (k_b T)$ and partition function
$
  Z = sum_({sigma}) exp(-beta H_"Ising" ({sigma})).
$

=== Potts Model
The Potts model @potts1952 generalizes the Ising model by allowing each spin to take on $q$ possible states. The $q$-state Potts model states take values in ${1, 2, ... q}$. The Hamiltonian is written as
$
  H_("Potts") ({sigma}) = -J sum_(angle.l i, j angle.r) delta(sigma_i, sigma_j).
$
$J$ again represents the ferromagnetic coupling that encourages neighboring spins to align in the same state. The original Potts was defined without the local fields, but they can be added if needed. The partition function becomes
$
  Z = sum_({sigma}) exp(-beta H_"Potts" ({sigma})).
$
The Potts model simplifies to an Ising model when $q=2$. Putting this together, the probability distribution from the Potts model is written as:
$
  P(bold(sigma)) = 1 / Z exp(-beta H_"Potts" ({sigma})).
$

It is important to note--anticipating further exploration--that the Potts model is over-parametrized: many different $(h, J)$ sets define exactly the same distribution. Without a gauge choice, the parameters will not be uniquely identifiable and the norms of $J$ can be misleading, hindering the process of optimization. Some common gauges, which will be explored in detail further, are reference states, zero-sum, and defining it implicitly via regularization.

== Direct Coupling Analysis
The statistical dependencies observed in MSAs encode important information about protein structure and function. However, naive measures such as covariance cannot distinguish between direct correlations, which arise from true physical or functional contacts, and indirect correlations, propagated via other residues. Direct Coupling Analysis (DCA) addresses this problem by introducing a method to separate direct correlations.

DCA assumes that sequences in the MSA are independent realizations of a Potts model distribution:
$
  P(bold(sigma)) = 1 / Z exp(sum_i h_i (sigma_i) + sum_(i<j) J_(i j) (sigma_i, sigma_j))
$ <potts>
where $h_i$ are site-specific fields and $J_(i j)$ are coupling parameters between residue pairs. The central goal of DCA is to infer the interaction parameters $J_(i j)$ that best explain the observed single- and pairwise frequencies. In practice, direct inference of these parameters is computationally challenging. Evaluating the partition function $Z$ is intractable for realistic proteins, as it requires summing over $q^L$ possible sequences. Thus, distinct implementations of DCA introduce approximations to bypass this computation. The original message passing implementation @weigt2009, directly attempted to solve for couplings through a slowly converging iterative scheme, making it inapplicable in practice. In this paper, we will explore some important algorithmic innovations in DCA following its conception.

= Evolution of DCA Methods

== Mean-field DCA (2011)
The mean-field Direct Coupling Analysis (mfDCA) algorithm, introduced by Morcos et al. @mfDCA, provided the first computationally feasible approximation of the Potts model. The idea is to approximate a small-coupling expansion, valid when correlations between sites are weak which reduces the inference problem to the inversion of a correlation matrix.

=== Method
To begin, the raw frequency counts suffer from sampling bias, so the weight of highly similar sequences is reduced. Each sequence $sigma^((a))$ is assigned a weight
$
  w_m = 1 / lr(|{b in{1, ..., M}:"sim"(sigma^((a)),sigma^((b))) >= x}|),
$ <weight>
where M is the number of sequences in the MSA and sim is their similarity #footnote[The original paper used "seqid" to represent percentage identity. We chose "sim" as the future methods adopted this notation.]. The $x$ represents the sequence identity threshold, which was empirically found to be best at $x=0.8$. The effective weight of a sequence $m$ is $w_m$, and the total effective number of sequences is
$
  M_("eff") = sum_(m=1)^M w_m.
$ <M_eff>
From these weights, the new regularized empirical single and pairwise frequencies are computed as
$
  f_i (A) &= 1 / (M_("eff") + lambda) (lambda / q + sum_(m=1)^M w_m delta(A, sigma_i^((m)))) \
  f_(i j)(A, B) &= 1 / (M_("eff") + lambda) (lambda / q^2 + sum_(m=1)^M w_m delta(A, sigma_i^((m))) delta(B, sigma_j^((m))))
$
where $q=21$ is the amino acid alphabet with the gap symbol, and $lambda$ is a pseudocount parameter used for regularization. The gauge choice employed in mfDCA is defining a reference state to remove redundancy among parameters and ensure identifiability, while preserving all physically relevant couplings. Typically, it is set as the last amino acid $A=q$. This gives us
$
  forall i,j: J_(i j)(a, q) = J_(i j)(q, a) = h_i (q) = 0
$
*Finding the Couplings* \
To circumvent the intractable $Z$ computation, mfDCA assumes weak correlations between sites, expanding the exponential in @potts by the Taylor series to first order@plefka1982 @georges1991. This results in the relation
$
  C_(i j)(A,B)=f_(i j) (A,B) - f_i (A)f_j (B).
$ <corr>
between the couplings and the connected correlation matrix. The couplings are then approximated as through naive mean-field inversion
$
  e_(i j) (A, B) = - (C^(-1))_(i j)(A,B),
$
where $C$ is treated as a $((q-1)L) times ((q-1)L)$ matrix, and the pair $(i, A)$ is regarded as a single index. The full derivation can be found in @app2. In practice, the correlation matrix is often singular without regularization. mfDCA opts to use a strong pseudocount ($lambda approx M_"eff"$) to stabilize the inversion and prevent spurious large couplings. 

*Pair Scoring* \
Once the estimate of the pair couplings $e_(i j) (A,B)$ is inferred, we need a way to rank residue pairs by their interaction strength. The $(q-1) times (q-1)$-dimensional coupling matrices need to map to a single scalar parameter. We can do this through the direct information (DI) @weigt2009. To begin, for each pair $(i, j)$ a two-site model is constructed:
$
  P_(i j)^"(dir)" (A,B) = 1 / Z_(i j) exp( e_(i j) (A,B) + tilde(h)_i (A) + tilde(h)_j (B)),
$
with auxiliary fields $tilde(h)_i, tilde(h)_j$ chosen such that
$
  sum_(B=1)^q P_(i j)^"(dir)" (A,B) = f_i (A), quad sum_(A=1)^q P_(i j)^"(dir)" (A,B) = f_j (B).
$
The direct information is the mutual information associated to this distribution:
$
  cal(S)^"DI"_(i j) = sum_(A,B=1)^q P_(i j)^"(dir)" (A,B)ln (P_(i j)^"(dir)" (A,B) )/ (f_i (A) f_j (B))
$
The top-scoring pairs are predicted to be structural contacts.

=== Limitations
While mfDCA represented a breakthrough in usability of DCA, the simplifying approximations impose limitations on the model:
- Weak coupling assumptions: the small-coupling expasion assumes nearly linear correlations, which can underestimate strong epistatic effects in proteins @Haldane_2019.
- Computational scaling: the inversion of the correlation matrix scales as $cal(O)(L^3)$, which is costly for very large MSAs.
- Pseudocount dependence: due to the algorithm's vast parameter size (in the order of $400N^2$) strong regularization is required. This makes the choice of the pseudocount $lambda$ significantly affect performance.

== Pseudo-likelihood Maximization DCA (2013)
While mfDCA provided a fast approximation, it relied on strong assumptions. Pseudolikelihood maximization (plmDCA) relaxed these by fitting sitewise conditional probabilities directly. Concretely, the inverse Potts problem on an alignment of length $L$ becomes $L$ coupled multinomial logistic regressions rather than a single optimization over the global partition function.

=== Method
We retain the Potts parametrization of fields and couplings introduced in @potts. Similar to mfDCA, sequences are reweighted using $w_m$ defined in @weight and the effective sample size $M_"eff"$ defined in @M_eff.

*Method Setup*

The pivotal idea is to optimize pseudolikelihoods. For site $r$, the conditional distribution given all other sites $sigma_(\\ r)$ is
$
  P(sigma_r = A | sigma_(\\ r)) = (exp(h_r (A)+sum_(i!=r)J_(r i)(A, sigma_i))) / (sum_(B=1)^q exp(h_r (B)+sum_(i!=r)J_(r i)(B, sigma_i)))).
$
The weighted sitewise negative log-pseudolikelihood is
$
  g_r (h_r, J_r) = - 1 / M_"eff" sum_(m=1)^M w_m log P(sigma_r^((m))| sigma_(\\r)^((m))).
$
and the global objective aggregates these points:
$
  cal(L)_"pseudo" (h, J) = sum_(r=1)^L g_r (h_r, J_r).
$
To curtail overfitting, we add convex $l_2$ penalties,
$
  R_(cal(l)_2) = lambda_h sum_(r=1)^L lr(norm(h_r)_2^2) + lambda_J sum_(1 <=i < j <=L) norm(J_(i j))_2^2
$
and minimize
$
  {h^"PLM", J^"PLM"} = arg min_(h, J) {cal(l)_"pseudo" (h, J) + R_(cal(l)_2) (h, J)}
$
The $cal(l)_2$ penalty fixes the gauge implicitly by selecting a unique representative among gauge-equivalent parameters. Each $g_r$ is precisely the loss of a multinomial logistic regression (softmax): classes are the $q$ states of $sigma_r$; features are one-hot encodings of ${sigma_i}_(i!=r)$ with $(L-1)(q-1)$ degrees of freedom after dropping a reference state. Hence (asymmetric) plmDCA is implemented as $L$ independent, weighted softmax problems with $cal(L)_2$ penalty (e.g. solved with L-BFGS or mini-batch SGD). Two variants are used in practice: asymmetric PLM, which fits each independently and then symmetrizes averaging:
$
  hat(J)_(i j) <- 1/2(J_(i j)^((i))+J_(i j)^((j)))
$
and symmetric (joint) PLM, which minimizes $cal(L)_"pseudo"$ over all parameters at once. Results are typically similar, but optimization in the symmetric variant can be heavier.

*Pair Scoring* \
For pair scoring, plmDCA avoids the previously defined Direct Information (DI) as it would introduce a regularization parameter for the pseudocounts, on top of the existing $lambda_h "and" lambda_J$. Instead it considers the Frobenius Norm (FN)
$
  norm(J_(i j))_2 = sqrt(sum_(k,l=1)^q J_(i j)(k,l)^2).
$
Unlike the DI score, the FN is not independent of gauge choice so a direct decision must be made. As noted in @weigt2009, the zero-sum gauge minimizes the FN, making it the most appropriate gauge choice available. The procedure is hence 
(i) convert to zero-sum gauge $J'_(i j)$: 
$
  J'_(i j) (k, l) = J_(i j) (k, l) - J_(i j) (dot, l) - J_(i j) (k, dot) + J_(i j) (dot, dot),
$ <jprime>
where "$dot$" denotes a simple average over the $q$ states at that position; (ii) compute the Frobenius norm 
$
  cal(S)_(i j)^"FN" = norm(J'_(i j))_2 = sqrt(sum_(k,l=1)^q (J'_(i j)(k,l))^2);
$
(iii) apply Average Product Correction (APC) adjusted from its use in @jones2012psicov to reduce phylogenetic bias,
$
  cal(S)_(i j)^"CN" = cal(S)_(i j)^"FN" - (cal(S)_(i dot)^"FN" cal(S)_(dot j)^"FN") / cal(S)_(dot dot)^"FN",
$
where $cal(S)_(i dot)^"FN"$ and $cal(S)_(dot j)^"FN"$ are row/column means and $cal(S)_(dot dot)^"FN"$ is the grand mean. 
Residue pairs $(i, j)$ are ranked by $cal(S)_(i j)^"CN"$ to predict structural contacts. 
=== Limitations
Despite its practical impact, plmDCA remains sensitive to sampling: accurate contact recovery still requires large $M_"eff"$, and sparse or biased MSAs degrade estimates. Phylogenetic and positional bias persist (reweighting and APC help but do not eliminate them), which can inflate false positives.


=== Fast plmDCA (2014)
To make plmDCA deployable at scale, two of the original authors and a third collaborator revisited the optimization choices of the method @fastplmDCA. In this second version, the asymmetric variant is used: $L$ independent, weighted softmax regressions are computed, one per site, and the final couplings are symmetrized by averaging, 
$
  J_(i j) = 1/2 (J_(i j)^((i))+J_(i j)^((j))).
$
This decomposition reduces per-solve dimensionality and, crucially, enables trivial parallelization across CPU cores or nodes, which is the primary source of the runtime gain. Since each $J_(i j)$ is regularized twice (once in the regression for $i$ and once for $j$), the coupling penalty parameter must be halved relative to the symmetric formulation. Before averaging, the two independently inferred coupling blocks are each shifted into the zero-sum gauge to ensure consistency. For scoring, the method adopts the Corrected Frobenius Norm (CFN): first compute the Frobenius norm of $hat(J)_(i j)$ excluding the gap state,
$
  "FN"_(i j) = sqrt(sum_(k,l != "gap") hat(J)_(i j) (k, l)^2),
$
then apply the Average Product Correction
$
  cal(S)^"CFN"_(i j) = "FN"_(i j) - ("FN"_(i :) "FN"_(: j)) / "FN"_(: :)
$

The combination of asymmetric regression and symmetrization yields similar contact accuracy as the original plmDCA, but at a fraction of the runtime, making large protein families and long sequences tractable in practice.

== Boltzmann Machine DCA (2018)
In mfDCA accuracy was traded for speed via the small-coupling inversion; in plmDCA, a product of sitewise conditionals was optimized, avoiding the global partition function. bmDCA@bmDCA instead proposes a solution closer to the statistical ideal: fitting the full Potts model by maximizing the true likelihood, using Monte Carlo Markov Chains to estimate intractable expectations. Improving on this, adabmDCA@adabmDCA introduces an adaptive procedure to keep sampling reliable and efficient. The result is a generative model that (by construction) matches the empirical one- and two-site statistics of the reweighted MSA.

=== Method
bmDCA aims to infer the full Potts model defined as @potts that reproduces the single- and pair-wise frequencies of the MSA. \ 
The training procedure follows the classical Boltzmann machine learning algorithm:
+ Sampling step. For current parameters {J, h}, estimate model frequencies $p_i, p_(i j)$ by MCMC
+ Update step. Adjust parameters when the estimated frequencies deviate from empirical ones: $Delta h_i prop f_i - p_i , quad  Delta J_(i j) prop f_(i j) - p_(i j)$
+ Iterate until empirical and model moments match within tolerance.

Regularization is imposed to stabilize the fit and avoid overfitting.

=== Limitations
- Computational cost: Direct Boltzmann machine learning with MCMC is extremely slow for realistic protein lengths ($10^7-10^9$ parameters).
- Sampling difficulty: Accurate moment estimation requires long chains and careful equilibration, making the approach impractical without approximations.
- Scaling: While bmDCA is theoretically optimal (full likelihood optimization), it is limited in practice to small systems where large-scale computation is feasible.

=== adabmDCA (2021)
adabmDCA revisits the same goal of fitting the full Potts by likelihood maximization, but introduces an adaptive MCMC procedure to make the approach more practical. Instead of relying on fixed sampling parameters, it dynamically tunes equilibration, waiting times, and chain persistance. adabmDCA pushes the bmDCA approach closer to practical usability, at the cost of a more complex implementation and still significant compute requirements. 

*Method* \
Let the MSA have length L, and M sequences $bold(sigma) = (sigma_1, ..., sigma_L)$. The maximum-entropy distribution that reproduces chosen moments is the Potts model defined in @potts. To reduce phylogenetic redundancy, assign each sequence a weight $w_m$ as in @weight. Empirical frequencies are reweighted and smoothed with pseudocount $lambda$ as
$
  f_i (A) = (1 - lambda) f_i^"data" (A) + lambda / q, quad f_(i j)(A,B) &= (1- lambda)f_(i j)^"data" (A, B)+ lambda/q^2
\
  f_i^"data" (A) = 1/M_"eff" sum_m w_m delta(A, sigma_i^((m))), quad f_(i j)^"data" (A, B) &= 1/M_"eff" sum_m w_m delta(A, sigma_i^((m))) (B, sigma_j^((m)))
$

A practical starting point for the system is the profile model, an independent-site Potts model where the first empirical moments are matched by means of the fields,
$
  h_i^"prof" (A) = log f_i (A) + "const",
$ 
but all couplings are set to zero $ J equiv 0$. In addition, zero or user-provided initial parameters can also work.

*Likelihood and moment matching*

The average log-likelihood of the MSA under $P(dot|J,h)$ is
$
  cal(L)(J, h) = 1/M sum_(m=1)^M [sum_i h_i (sigma_i^((m))) + sum_(i<j) J_(i j) (sigma_i^((m)), sigma_j^((m)))]-log Z(J, h).
$
As an exponential-family model, $cal(L)$ is concave in the natural parameters, so gradient ascent converges to the unique optimum. The gradients are moment gaps:
$
  (partial cal(L)) / (partial h_i (A)) = f_i (A) - p_i (A), quad (partial cal(L)) / (partial J_(i j) (A,B)) = f_(i j) (A,B) - p_(i j) (A, B),
$
where $p_i$ and $p_(i j)$ are model marginals under current $(J, h)$. Hence the update
$
  h_i^(t+1)(A) <- h_i^t (A) + eta_h [f_i (A) - p_i^((t))(A)], \
  J_(i j)^(t+1)(A, B) <- J_(i j)^t (A, B) + eta_J [f_(i j) (A, B) - p_(i j)^((t))(A, B)].
$
drives the model toward exact moment matching f = p.
The obstacle is that $p_i$ and $p_(i j)$ are not analytically computable at scale; adabmDCA estimates them by MCMC at each epoch.

*Estimating model expectations via adaptive MCMC*

At training epoch $t$, run $N_s$ independent Markov chains using Metropolis-Hastings @Metropolis1953 @mh-algorithm (Gibbs sampling strategy may also be used @gibbs-algorithm), each producing $N_c$ samples after an equilibration period $T_"eq"$ and with an inter-sample waiting time $T_"wait"$. The Monte Carlo estimators are
$
  p_i^((t)) (A) = 1/(N_s N_c) sum_(m=1)^(N_s N_c) delta(A, sigma_i^((m))(t)), \ p_(i j)^((t)) (A, B) = 1/(N_s N_c) sum_(nu=1)^(N_s N_c) delta(A, sigma_i^((m))(t)) delta(B, sigma_j^((m))(t)).
$
Chains may be transient, reinitialized every epoch, or persistent, initialized only at the first epoch. Equilibration is often sped up through persistence of chains. For more details on the adaptive scheme, refer to @app3.

*Convergence and quality control*

A practical convergence proxy is the difference between the empirical and the model two-site connected correlations:
$
  epsilon_c = max_(i,j,a,B) abs(c_(i j)^"model" (A,B)-c_(i j)^"emp" (A,B)) quad"where" \ quad c_(i j)^"model" = p_(i j)-p_i p_j, quad c_(i j)^"emp" = f_(i j)-f_i f_j
$
with a target $epsilon_c approx 10^-2$ In addition to this, some other commonly used diagnostics is the Pearson correlation between $c^"model"$ and $c^"emp"$, one- and two-site fitting errors, and optionally, a third connected correlation on a subset of triples is used to assess generative fidelity beyond pairwise constraints.

*Priors and sparsity*

A fully connected Potts model has $~ (L(L-1))/2 q^2 + L q$, meaning a number in the order of $10^7$-$10^9$ parameters for a realistic $L$ @adabmDCA. Due to the finite sample size of MSAs, they rarely contain enough independent information to estimate all of them robustly. Not controlling this uncertainty could lead to overfitting, high variance or instability, and bad conditioning in the model. To address these issues, adabmDCA employs two complementary strategies.
First, we can place a prior on $P(J, h)$ and maximize the posterior, equivalent to adding penalties to the objective. The two standard choices are the $cal(l)_1$ and $cal(l)_2$ priors:
$
  R_(cal(l)_1) (J, h) = theta_(1, h) sum_i norm(h_i)_1 + theta_(1, J) sum_(i<j) norm(J_(i j))_1 \
  R_(cal(l)_2) (J, h) = theta_(2, h) sum_i norm(h_i)_2^2 + theta_(2, J) sum_(i<j)norm(J_(i j))_2^2
$
Under $cal(l)_2$, the gradients are shrunk toward zero:
$
  partial / (partial h_i (A)) : f_i (A) - p_i (A) - theta_(2,h) h_i (A), \ partial / (partial J_(i j) (A,B)) : f_(i j) (A, B) - p_(i j) (A, B) - theta_(2,h) J_(i j) (A,B). 
$
Under $cal(l)_1$, they include subgradient terms that promote exact zeros.
$
  partial / (partial h_i (A)) : f_i (A) - p_i (A) - theta_(1,h) "sign"(h_i (A)), \ partial / (partial J_(i j) (A,B)) : f_(i j) (A, B) - p_(i j) (A, B) - theta_(1,h) "sign"(J_(i j) (A,B)). 
$
$cal(l)_2$ reduces variance and improves conditioning by smoothly shrinking all parameters. It selects a unique gauge and tends to preserve the relative ordering of strong couplings. On the other hand, $cal(l)_1$ induces sparsity by zeroing weak parameters, also reducing overfitting, though at a cost of biasing small effects downward. Generally, in stochastic settings, elastic-net is used (a combination of both parameters), that stabilizes training near zero. In addition, a separate parameter is used for fields and couplings, typically regularizing $J$ more strongly than $h$.

The second method is introducing sparsity via pruning or decimation. The reason for this is that true contact maps in nature are indeed sparse; most residue pairs are not in direct physical contact. Encoding this structural prior can reduces variance and speed up learning. There are two approaches:
+ A priori topology. Reduce the number of parameters by starting from a restricted edge set, for example pairs with high MI, and learn only those $J_(i j)$, omitting the rest.
+ Information-based decimation. In this approach, start dense and iteratively remove the least informative couplings until target sparsity. This can be done by comparing the KL divergence for a candidate element. This directly controls overfitting by only keeping the parameters that actually affect the model's predictions.
In short, pruning and decimation prevent overfitting by removing parameters that don't materially alter the model, and have the added benefit of aligning with the biological prior of contact sparsity.


== Autoregressive Network DCA
ArDCA was built to explore the ability of generative models for protein design coming from sequence-data. It was not the first generative model, DeepSequence, ... ArDCA emerges as a capable model for the extraction of structural and functional protein information which encoded in rapidly growing protein databases.

=== Method
In arDCA, the exponential-family MaxEnt distribution is replaced with a conditional probability model, where each residue is predicted from the previous ones. This arises from the chain rule decomposition of the joint probability distributions:
$
  P(bold(sigma)) = product_(i=1)^L P(sigma_i|sigma_1,...,sigma_(i-1))
$
*Parameter inference* \
The inference of the parameters is done through likelihood maximization. Following a Bayesian setting with a uniform prior, the optimal parameters are those that maximize the probability of the data #footnote[Here $a$ is used in place of A for sake of space]:
$
  {J^*, h^*} &= arg max_{J, h} P(cal(M)|{J, h}) \
  &= arg max_{J, h} log P(cal(M)|{J, h}) \
  &= arg max_{J, h} sum_(m=1)^M log product_(i=1)^L P(a_i^m|a_(i-1)^m,...,a_1^m) \
  &= arg max_{J, h} sum_(m=1)^M sum_(i=1)^L log P(a_i^m|a_(i-1)^m,...,a_1^m)
$ <lhmax>
Each $h_i (a)$ and $J_(i j)(a,b)$ is present in only one conditional probability $P(a_i|a_(i-1),...,a_1)$, thus we can maximize each conditional probability independently in @lhmax:
$
  {J_(i j)^*, h_i^*} = arg max_{J_(i j), h_i} sum_(m=1)^M [h_i (a_i^m) + sum_(j=1)^(i-1) J_(i j) (a_i^m, a_j^m) - log z_i (a_(i-1)^m, ..., a_1^m)]
$
where
$
  z_i (a_(i-1), .., a_1) = sum_(a_i) exp(h_i (a_i)+sum_(j=1)^(i-1) J_(i j)(a_i, a_j))
$ <cpnorm>
is the normalization factor of the conditional probability of $a_i$. \
Taking the derivative with respect to $h_i(a)$ or $J_(i j) (a,b)$, with $j=1, ..., i-1$, we get:
$
  0 &= 1/M sum_(m=1)^M [ delta(a, a_i^m) - (partial log z_i (a_(i-1)^m, .., a_1^m)) / (partial h_i (a))] \ 
  0 &= 1/M sum_(m=1)^M [ delta(a, a_i^m) delta(b, a_j^m) - (partial log z_i (a_(i-1)^m, .., a_1^m)) / (partial J_(i j) (a,b))].
$
Using @cpnorm, we find
$
  (partial log z_i (a_(i-1)^m, .., a_1^m)) / (partial h_i (a)) &= P(a_i = a |a_(i-1)^m,...,a_1^m) \ 
  (partial log z_i (a_(i-1)^m, .., a_1^m)) / (partial J_(i j) (a,b)) &= P(a_i = a |a_(i-1)^m,...,a_1^m) delta(a_j^m, b).
$
The set of equations reduces to a simple form:
$
  f_i (a) &= angle.l P(a_i = a |a_(i-1)^m,...,a_1^m) angle.r_cal(M), \
  f_(i j) (a, b) &= angle.l P(a_i = a |a_(i-1)^m,...,a_1^m) delta(a_j^m, b) angle.r_cal(M),
$
where $angle.l dot angle.r_cal(M) = 1/M sum_(m=1)^M dot^m$ denotes the empirical average. The first variable, $i=1$, is unconditioned, therefore the coupling equation is $J equiv 0$
and the equation for the field is the profile model also used in bmDCA, $h_1(a) = log f_1(a) + "const"$.

Unlike bmDCA, the equations do not enforce exact matching between model marginals and empirical frequencies. The ability to reproduce the frequencies is thus a good proxy for the generative properties of the model, on top of the fitting quality of current parameters. In practice, the parameters are updated using gradient descent on the likelihood. As we are able to take the expectations directly over the MSA, not needed MCMC, the gradients are exact. For regularization, arDCA makes use of $ell_2$ penalties, with different strengths depending on whether the task is generating sequences or just contact prediction. It was noted that small regularization values improved the generative quality while larger were beneficial for contact prediction.

=== Key Contributions

Unlike previous bmDCA, expectations are taken directly over the data rather than approximated by MCMC. This means the gradients will be exact which in turn will speed up training.Because of this, ArDCA is a lightweight approach that performs at a similar, or better, accuracy of previous iterations, but at a substantially lower computational cost (a factor between $10^2$ and $10^3$) @Trinquier2021. In addition, it allows for the calculation of exact sequence probabilities, which is crucial in homology detection.

= Implementation and Exploration
The original ArDCA model was developed in Julia, a language designed for high-performance numerical and scientific computing. Julia's strengths lie in just-in-time (JIT) compilation, seamless handling of linear algebra, and built-in support for parallelism, making it well-suited for implementing large-scale statistical models.

In contrast, the re-implementation was carried out in Python, a widely adopted language for machine learning and data science. While Python itself is generally slower due to its interpreted nature, performance-critical operations are offloaded to highly optimized libraries (e.g., NumPy, SciPy, PyTorch), which use underlying C/C++ or Fortran code. This ecosystem makes Python particularly attractive for model development, as it provides a wide range of frameworks, pre-built optimization routines, and strong community support.

*Important Differences Between Julia and Python*

There are some features of the languages that makes the porting not a simple one-to-one translation. To begin, Julia uses column-major order, which akes column-wise operations very fast and easily cached. This is in complete contrast to Python which is structured with row-major order where row-wise operations are more optimized. The impact of this difference is that loops and reshaping operations need to be reconsidered for memory efficiency. In addition, Julia is 1-indexed whereas Python is 0-indexed. This has no effect on the performance, but it changes the bounds on loops, the indexing of objects, and the range functions.

Furthermore, loops have different behaviors in the languages. Julia loops are compiled down to efficient machine code with JIT. There exist a macro that removes the bounds checking, making loops fast without the need for vectorization. On the other hand, pure Python loops are slow. To reach the same level of efficiency vectorization is needed in Python, implemented via NumPy or broadcasting.

The final important difference is in compilation versus interpretation. Julia JIT compiles functions to machine code, which has an initial compilation cost, but allows future calls to run much faster. On the other hand, Python is interpreted, thus heavy computations require external libraries which are precompiled in other, more efficient, languages.

== Implementation details
We implement the autoregressive network within the PyTorch ecosystem, using the $mono("nn.Module")$ interface to leverage its training and optimization capabilities. The model was trained on the CPU with an Apple M1 chip. The central operation of the model is the computation of the autoregressive logits, which define the conditional distribution at each sequence position. 

*Logit Computation* \
For site $i$, the logit vector is given by
$
  z[m,i,a] = h[i,a] + sum_(j<i) sum_b J[i,j,a,b] X[m,j,b]
$ 
where h and J are the local biases and pairwise couplings, and  $X$ is a tensor holding the one-hot encoding of the symbol $b$ observed at site $j$ in sequence $m$.

To compute the autoregressive logits, our initial implementation employed a masked matrix multiplication, performed in a single step using the `einsum` operator. While correct, and more efficient than naive loops, this aproach proved computationally inefficient: the full interaction matrix was calculated, even though the contributions from the upper-triangular part were masked to zero.
To address this, we devised a more efficient formulation that explicitly exploits the block-sparse structure of the coupling matrix. Specifically, we collect all lower-triangular index pairs $(i, j)$ with $j<i$, extract the corresponding coupling blocks extracted $J[i, j]$, and accumulate their contributions through an `einsum` operation. The resulting implementation only computes terms required by the autoregressive factorization, removing unnecessary computation.

```Python
def compute_ar_logits(self, X_oh: torch.Tensor):
        """z[m,i,a] = h[i,a] + sum_{j<i} sum_b J[i,j,a,b] * X[m,j,b]"""
        X_oh = X_oh.to(self.J.dtype)
        M, L, q = X_oh.shape
        logits = self.h_pos.unsqueeze(0).expand(M, -1, -1).clone()  # (M,L,q)

        # collect lower-triangular index pairs
        i_idx, j_idx = torch.tril_indices(L, L, offset=-1)
        J_blocks = self.J[i_idx, j_idx]   # (n_pairs, q, q)
        X_blocks = X_oh[:, j_idx]         # (M, n_pairs, q)

        # compute pairwaise contributions and accumulate into logits
        contrib = torch.einsum("mpq,pqr->mpr", X_blocks, J_blocks) 
        logits = logits.index_add(1, i_idx, contrib)
        return logits
```
With this formulation, the computational complexity scales with the number of lower-triangular pairs rather than the full $L times L$ coupling matrix. For long sequences, this is especially beneficial. In addition, by avoiding the construction of the full masked interaction matrix, the memory footprint is roughly halved, from $L^2$ couplings to $L(L-1) / 2$.

Throughout the model, an explicit binary mask `J_mask` is maintained that encodes the lower-triangular structure, to ensure consistency in the autoregressive aspect. After each optimization step, unused entries of $J$ are clamped to zero.

*Training details* \
The model is initialized with the local bias of the first position, $h[0]$, being defined from the empirical frequencies. The rest of the parameters are drawn from small random distributions, an exploration will be done on initial conditions in the following section.
Training also incorporates sequence reweighting to correct for redundancy in MSAs. The weights for each sequence are computed at the start, with parameter `theta` controlling the sequence similarity threshold.
The model is set up to be optimized with both L-BFGS, suited for energy-based models, as well as AdamW, offering a more scalable and robust option for large datasets. Regularization is applied separately to the field and coulings via $ell_2$-norm penalties with hyperparameters $lambda_h$ and $lambda_J$.

*Evaluation* \
For evaluation, the model reports the negative log-likelihood (NLL), perplexity, and the effective sample size, to make the results more easily comparable across different datasets. Finally, an ancestral sampling procedure is provided, which generates new sequences position by position using the learned conditional distributions. This procedure is inherently slower than parallel inference, but it guarantees that samples are consistent with the autoregressive structure.


== Experiments

=== Effect of Gaps in Multiple Sequence Alignments
When constructing MSAs, gaps are intro introduced to better align homologous positions across sequences. These gaps serve two main purposes: (i) ensuring residues align correctly at homologous sites, and (ii) extending sequences to a specific length for models like ArDCA, that require fixed length input.

Since the presence of gaps can influence the statistical properties of the alignment, we examined how restricting the fraction of gaps affects results. Specifically, we introduced two filtering parameters:
- `max_gap_fraction`: sets the maximum allowed fraction of gaps per sequence. 
- `max_col_gap_fraction`: sets the maximum allowed fraction of gaps in a column.

For the first experiment, we varied `max_gap_fraction`, while holding the other fixed, with three thresholds: $0.05, 0.1, "and" 1$. 1 does not exclude any sequences, therefore is used as a baseline for comparison. The difference in time of the models was also recorded, but given the small nature of the change, no conclusive results were observed. 

#align(center)[
#block(width: 80%,)[
  #figure(
    table(
      columns: (auto, 1fr, 1fr, auto),
      table.vline(x: 1, start: 1),
      table.header[`max_gap_fraction`][M][$M_"eff"$][Perplexity],
      [0.05 (-6%)], [12784], [3876.5], [1.72],
      [0.1 (-1.6%)], [13376], [4145.4], [1.63],
      [1], [13600], [4248.8], [1.77]
    ),
    caption: [Results for `max_gap_fraction` experiment]
  )
]]
Restricting the allowed fraction of gaps resulted in the removal of sequences and a corresponding decrease in the effective number of sequences $M_"eff"$. The strictest filterint criterion ($0.03$) removed approximately 8.7% of sequences relative to the baseline, yielding a smaller dataset and slightly worse perplexity. The lowest perplexity ($1.63$) and fastest time ($344 s$) was observed when 10% of gaps were excluded. Although the full dataset had a higher $M_"eff"$, the higher perplexity suggests that allowing too many gaps may negatively impact the quality fo the alignment.

A similar exploration was performed for the column gaps with values  $0.025, 0.1, 1$. Again, 1 represents the baseline with no columns removed.
#align(center)[
#block(width: 80%,)[
  #figure(
    table(
      columns: (auto, 1fr, 1fr, auto),
      table.vline(x: 1, start: 1),
      table.header[`max_col_gap_fraction`][L][$M_"eff"$][Perplexity],
      [0.025 (-7.5%)], [49], [4363.0], [1.79],
      [0.1 (-3.7%)], [51], [4105.9], [1.62],
      [1], [53], [4248.8], [1.77]
    ),
    caption: [Results for `max_col_gap_fraction` experiment]
  )
]]

Both sequence-level and site-level filtering yield a better result when using moderate thresholds, where exclude just the data with the most gaps improves model performance. Overly aggressive filtering discards useful information and hurts perplexity, while using the full dataset also doesn't reach the best performance.

=== Datasets 
Different protein families were analyzed and explored in the evaluation of the model.

#align(center)[
#block(width: 70%,)[
  #table(
    columns: (auto, 1fr, 1fr, 1fr),
    table.vline(x: 1, start: 1),
    table.header[$bold("Protein Family")$][M][L][$M_"eff"$ #footnote[Calculated with sequence identity threshold 0.8]],
    [PF00014], [13600], [53], [4208.5],
    [PF00076], [137605], [70], [],
    [PF00072], [823798], [112], [],
    [PF13354], [7515], [202], [],
  )
]]

=== Evaluation Metrics



= Results and Discussion

== Training Behavior
ArDCA benefits from exact gradients, meaning the loss curves are smooth and the convergence is direct. In all of the trials, there was minimal overfitting, the validation loss decreased just as the training loss did.

== Generative Model Quality
The final test for the models is their ability to generate accurate sequences. For this we can use AlphaFold's plDDT @AlphaFold_pLDDT, 


== Limitations and Insights
The Python implementation, with its lower-triangular masked multiplication is efficient for Python standards, but due to Julia's improved the memory management and computation power, still doesn't compare.

= Conclusion
In this work, I reimplemented the autoregressive DCA model in Python and demonstrated its viability as a practical tool for protein sequence analysis. The PyTorch-based implementation introduced an efficient block-sparse formulation for computing autoregressive logits, reducing the computational cost and memory usage. Empirical validation across several protein families confirmed the correctness and robustness of the model, with smooth training dynamics and consistent improvements in negative log-likelihood, perplexity, and structural plausibility as measured by pLDDT. Additional experiments on gap filtering highlighted the importance of sequence quality control, while tests on larger Pfam families showed that the method remains tractable at scale. In the future, this work could benefit from GPU acceleration, to improve its speed, but also allow it to more easily work with large datasets. In addition, 

#pagebreak()
#bibliography("references.bib")

#pagebreak()
#show: appendix

= Full Langrage Multipliers Calculation for Maximum Entropy Principle <app1>
$
  p_i = e^( - lambda -mu f(x_i)) \
  sum_i p_i = 1 \
  sum_i e^( - lambda -mu f(x_i)) = 1", factor out " e^(-lambda) \
  e^(-lambda) sum_i e^( -mu f(x_i)) = 1,
$
Define the partition function $Z(mu):$
$
  Z(mu) = sum_i e^( -mu f(x_i))", thus" \
  e^(-lambda)Z(mu) = 1 => lambda = ln Z(mu)
$

#pagebreak()
= Small Coupling Mean-Field Derivation <app2>
#counter(math.equation).update(0)
Start with the perturbed Hamiltonian
$
  cal(H)(alpha) = -alpha sum_(1 <= i < j <= L) e_(i j) (A_i, A_j) - sum_(i=1)^L h_i (A_i)
$
which allows interpolation between independent couplings, $alpha = 0$, and the original model, $alpha = 1$. We also define the Gibbs potential
$
  -cal(G)(alpha) = ln[sum_({A_i|i=1,...,L})e^(-cal(H)(alpha))]-sum_(i=1)^L sum_(B=1)^(q-1)h_i (B) P_i (B)
$ <gibbspotential>
as the Legendre transform of the free energy $cal(F) = -ln Z$. The fields can be found via 
$
  h_i (A) = (partial cal(G) (alpha)) / ( partial P_i (A) ),
$
and
$
  (C^(-1))_(i j) (A, B) = (partial h_i (A)) / (partial P_j (B)) = (partial^2 cal(G)(alpha)) / (partial P_i (A) partial P_j (B)).
$
Our aim is to expand the Gibbs potential up to first order around the independent-site case $alpha=0$,
$
  cal(G)(alpha) = cal(G)(0) + lr((d cal(G)(alpha)) / (d alpha) |)_(alpha = 0) alpha + cal(O)(alpha^2)
$ <gibbsalpha>
*Independent-site approximation* \
To start, let us consider the Gibbs potential in $alpha=0$. In this case, the Gibbs potential equals the negative entropy of an ensemble of $L$ uncoupled Potts spins,
$
  cal(G)(0) &= sum_(i=1)^L sum_(A=1)^q P_i (A) ln P_i (a) \
  &= sum_(i=1)^L sum_(A=1)^(q-1) P_i (A) ln P_i (a) + sum_(i=1)^L [1-sum_(A=1)^(q-1) P_i (A)]ln[1- sum_(A=1)^(q-1) P_i (A)]
$ <g0>

*Mean-field approximation* \
To get the first order in @gibbsalpha, we have to determine the derivative at $alpha=0$. Recalling the definition of Gibbs potential in @gibbspotential, 
$
  (d cal(G)(alpha)) / (d alpha) &= (-d) / (d alpha) ln Z (alpha) - sum_(i=1)^L sum_(A=1)^(q-1) (d h_i (A)) / (d alpha) P_i (A) \ 
  &= -sum_({A_i})[sum_(i<j)e_(i j) (A_i, A_j )+sum_i (d h_i (A)) / (d alpha)](e^(-cal(H)(alpha)))/Z(alpha) - sum_(i=1)^L sum_(A=1)^(q-1) (d h_i (A)) / (d alpha) P_i (A) \ 
  &= - lr(angle.l sum_(i<j) e_(i j) (A_i, A_j)angle.r)_alpha
$
The first derivative of the Gibbs potential with respect to $alpha$ this is the average of the coupling term in the Hamiltonian. At $alpha=0$, this average can be done easily due to the joint distribution of all variables becoming factorized over the single sites,
$
  lr((d cal(G)(alpha)) / (d alpha) |)_(alpha = 0) = -sum_(i<j) sum_(A,B)e_(i j) (A_i, A_j)P_i (A) P_j (B).
$
Plugging this and @g0 into @gibbsalpha, we find the first-order approximation of the Gibbs potential. The first and second partial derivatives with respect to $P_i (A)$ provide self-consistent equations for the local fields,
$
  (P_i (A)) / (P_i (q)) = exp( h_i (A) + sum_({j|j!=i}) sum_(B=1)^(q-1) e_(i j) (A, B) P_j (B))
$
and the inverse of the connected correlation matrix,
$
  lr((C^(-1))_(i j) (A, B)|)_(alpha=0) = cases(-e_(i j) (A, B) "for" i != j, (delta(A,B))/(P_i (A)) + 1/(P_i (q)) "for" i=j ) quad .
$
This equation allows us to solve the original inference in the mean-field approximation in a single step. To determine the marginals of the empirical freqiencies, we just need to determine the empirical connected correlation matrix and invert it to get the couplings $e_(i j)$. 

#pagebreak()
= Adaptive MCMC Sampling in adabmDCA <app3>
To ensure accurate estimation of the marginals, adabmDCA monitors both equilibration and decorrelation of the Markov chains through sequence overlaps. For two sampled sequences $s^i_n, s^k_m in cal(A)^L$, the normalized overlap is defined as
$
  O(s_n^i,s_m^k) = 1/L sum_(j=1)^L delta(s_n^i (j),s_m^k (j))
$
Three types of overlaps are used:
- External overlap ($O_"ext"$): between samples from different chains at the same time $n$.
- Internal-1 overlap ($O_"int1"$): between consecutive samples of the same chain at times $n$ and $n+1$.
- Internal-2 overlap ($O_"int2"$): between samples of the same chain at times $n$ and $n+2$.
At each epoch, adabmDCA computes averages $mu_alpha$ and standard errors $sigma_alpha$ ($alpha in {"ext","int1","int2"}$). In equilibrium, all three overlaps should agree within statistical error.
*Adaptive update of $T_"wait"$*\
To control correlation between successive samples, $T_"wait"$ is adjusted dynamically: \
Increase $T_"wait"$ (double it) if:
$
  abs(mu_"ext" - mu_"int2") > 5 sqrt(sigma_"ext"^2+sigma_"int2"^2)
$
i.e. samples at lag $2T_"wait"$ are still too correlated.\
Reduce $T_"wait"$ (average with last pre-increase value) if:
$
  abs(mu_"ext" - mu_"int1") < 5 sqrt(sigma_"ext"^2+sigma_"int1"^2)
$
i.e. decorrelation already occurs at lag $T_"wait"$. Equilibration time is then set to $T_"eq" = 2 T_"wait"$.

This adaptive scheme ensures that chain samples are equilibrated, initial bias is removed after $T_"eq"$ steps, and decorrelated, independence between samples at lag $T_"wait" - 2T_"wait"$. As a result, the estimated marginals match those of the equilibrium Potts distribution, guaranteeing stable gradient estimates and convergence of the Boltzmann machine.