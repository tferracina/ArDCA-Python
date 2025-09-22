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

#set page(numbering: none)
#pagebreak()
#pagebreak()
#counter(page).update(1)
#set heading(numbering: none)
== Acknowledgements
\
Vorrei innanzitutto ringraziare i miei nonni per l'amore e supporto durante questi miei tre anni a Milano. \

Grazie anche a i miei genitori per le opportunità che mi hanno sempre dato, per il sostegno e la fiducia. Un grazie speciale alle mie sorelle Carlotta e Martina le cui chiamate mi hanno permesso di fare pause divertenti durante la stesura di questa tesi. Grazie di cuore a Federica che mi ha supportato e sopportato, grazie a Tebe il mio collega e consigliere e a Stefano, George, Federico, Ulysse, per avere letto, o almeno fatto finta di leggere, le mie bozze.

Vorrei infine dedicare questo lavoro alla mia Nonna Nadia che ha sempre fatto l'impossibile per me. Nonna, ti voglio tanto bene.


#pagebreak()
#pagebreak()
#counter(page).update(1)
#set page(numbering: "1")
#set heading(numbering: "1.")
#outline()

#pagebreak()
= Introduction

Proteins are biomolecules that are fundamental to nearly all biological processes. Their diverse roles include transporting nutrients, catalyzing chemical reactions, providing structural support, and more. The function of a protein is determined by the composition of its primary sequence created from amino acids. A single protein sequence can vary greatly in length and order of its amino acids, leading to a very large number of possible configurations. The massive size of the protein sequence space makes exhaustive experimental exploration infeasible. 

Advances in sequencing technologies has ushered in a rapid growth of biological databases, transforming biology into a data-rich discipline. Large repositories such as UniProt, Pfam, and the Protein Data Bank store millions of sequences and structures, providing crucial resources for computational approaches @weigt2020. This wealth of data has enabled the development of advanced statistical and machine learning models capable of simulating protein sequence evolution, predicting structural conformations, and generating novel sequences with desired properties. In addition, breakthroughs in protein structure prediction--most notably the Nobel Prize winning AlphaFold @alphafold\-- show that computational methods can rival experimental accuracy, in a much cheaper manner.

In this work, we will implement arDCA, an efficient autoregressive network for protein sequence generation leveraging the Direct Coupling Analysis (DCA) method in Python. The biological background will be introduced in  Section 2. The mathematical foundations will be laid out in Section 3. Section 4 will describe previous iterations of Direct Coupling Analysis methods and arDCA, while Section 5 will explore our implementation.

= Biological Background

== Proteins and their structure 
Proteins are essential biological molecules responsible for a wide range of functions in living organisms. Despite their functional diversity, all proteins are polymers of the same set of standard building blocks: the 20 canonical amino acids. Amino acids share a common core structure consisting of a central carbon atom bonded to a hydrogen atom, an amino group, a carboxyl group, and a variable side chain. The side chain is the defining feature of an amino acid, giving rise to differences in size, shape, chemical reactivity, and polarity. Proteins, or polypeptides, are formed through peptide bonds of distinct amino acids. Certain atoms are lost in this process, thus within a protein sequence, individual amino acids are typically referred to as residues. Generally, protein sequences are made up of between 50 and 2000 amino acids. The ordering of the amino acids dictates how a protein folds into its three-dimensional structure, known as its conformation. Although each conformation is unique, two common folding patterns occur in many proteins: the $alpha$ helix and the $beta$ sheet @alberts2002. 

#figure(
  image("code/out/PF13354protein.png", width: 40%),
  caption: [Protein from PF13354 folded using `ColabFold` @Mirdita2022_ColabFold]
)

Each amino acid is chemically distinct and can occur at any position in a protein chain, giving rise to $20^n$ possible polypeptide sequences of length $n$. For a typical protein of 300 amino acids, the number of possible sequences is astronomically large ($approx 2 dot 10^390$). However, only a small fraction of these sequences are capable of folding into a stable three-dimensional conformation. 

== Protein families and evolutionary information
Proteins do not evolve in isolation; they often belong to protein families, groups of proteins that share a common evolutionary origin, therefore exhibit similar sequence features and functional properties @ebi_protein_families. Over evolutionary timescales, mutations accumulate in these families: some are beneficial, altering the protein activity in ways that give rise to new functions, while many others are neutral and have no effect on stability or activity. Harmful mutations, by contrast, disrupt protein folding or its function. These destructive changes are eliminated through natural selection. The result is a collection of homologous proteins that retain overall structural and functional characteristics, but also display sequence variability that encodes the evolutionary history of the family @alberts2002. 

The study of evolutionary history and relationships among biological entities is referred to as phylogenetics. Protein families provide the domain for phylogenetic analysis, as examining families provides insight that cannot be obtained from a single sequence @EMBL-EBI_Phylogenetics. From homologous proteins, we can detect important correlated mutations between amino acid positions which represents constraints enforced to maintain protein integrity by evolution. These statistical patterns are crucial for computational approaches which learn them to predict three-dimensional structure and understand protein function.

== Multiple sequence alignments
To extract important statistical information from protein families, the homologous sequences need to be organized in a systematic way. This is done through a multiple sequence alignment (MSA), in which three or more sequences are arranged so that homologous sites are placed in the same column @WILTGEN201938. Alignment gaps are introduced throughout the MSA to maximize the positional correspondence of sequences and enable them to have the same fixed length. 
#figure(
  image("code/out/MSA.png"),
  caption: [MSA of first 10 proteins in PF00014, created using `pyMSAviz` @moshi4_pyMSAviz_2024]
)

Aligning sequences in this manner, patterns of conservation and variation are revealed across the family. Conserved positions, those which are unchanged in multiple sequences, represent sites that are critical for maintaining structure or function, while variable positions indicate sites that can tolerate mutations without disruption to the protein. Beyond conservation, MSAs also capture covariation: pairs of positions that mutate in a correlated way across sequences. These covariation signals reflect couplings, where a mutation at one site requires compensation at another to maintain protein integrity @biom14121531.


== Protein contacts
One of the earliest challenges for researchers in computational biology--historically encouraged by the CASP (Critical Assessment of Structure Prediction) @PredictionCenter competition--was understanding how the linear sequence of amino acids folds into its conformation. A key insight was the role of contacts, pairs of residues that, while possibly far apart in the linear sequence, end up spatially close in the folded three-dimensional structure. Since structures can be represented in Cartesian coordinates $(x,y,z)$, contacts are usually defined by distance thresholds, with two residues considered in contact if the distance between the selected atoms is below a set threshold, usually 8 Ångströms (Å). These contacts are further categorized into short, medium, and long range, with long-range contacts playing a disproportionately important role in determining the fold. Importantly, this illustrates that correlations in protein sequences are not limited to neighboring sites along the chain, many of the most informative dependencies span large separations @adhikari2016. 

== The problem of inference

The central challenge in protein sequence analysis is disentangling the true structural and functional constraints of protein families from spurious correlations introduced by other factors. Correlations can arise indirectly, for example, if residue A correlates with B, and B with C, then A and C may appear correlated without a direct interaction. Moreover, shared evolutionary history (phylogeny) and biases in sequence databases can further obscure the true couplings that underlie protein folding and function @dietler2023. Addressing this problem requires a mathematical framework capable of distinguishing direct from indirect correlations in the data.

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
To find the probability distribution $P(bold(sigma))$ that will satisfy constraints, in other words reproduce the empirical marginals $f_i (A)$ and $f_(i j) (A, B)$, we can use the maximum entropy principle @Jaynes1957a @Jaynes1957b. The first step to set up the MEP is extracting information from the system. Usually, this information is given in the form of averages of functions $angle.l f_k (x) angle.r$. For example, in a physical system, one could compute average magnetization or energy of the observed system. We also need to define a probability of occupancy of states, $p(x)$, which runs over all the possible states. This distribution satisfies the standard properties of a probability distribution: for all states $x$, $0<=p(x)<=1$, and the distribution is normalized, $sum_x p(x) =1$ @MITinfent. Our uncertainty on the system is expressed quantitatively through Shannon entropy $S$ @Shannon1948:
$
  S = - sum_x p(x) ln p(x) 
$ <entropy>
In general, the distribution that maximizes entropy under no constraints is the uniform distribution. However, once constraints are imposed, the uniform is no longer suitable. Our system is constrained by the probability distribution's normalization condition and the expectation constraints:
$
  sum_x p(x) = 1, quad "and" quad sum_x p(x) f_k (x) = angle.l f_k angle.r quad (k = 1, ..., m),
$ <constraints>
where each $f_k (x)$ represents a feature, and $angle.l f_k (x) angle.r$ its empirical average. By selecting Shannon entropy as our measure of information, we can introduce Lagrange Multipliers $lambda_0, lambda_1, ... lambda_m$ to maximize @entropy subject to the constraints defined in @constraints, yielding:
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

The maximum entropy of the distribution  then reduces to#footnote[The full derivation can be found in @app1.]:
$
  S_max = lambda_0 + sum_(k=1)^m lambda_k angle.l f_k angle.r
$
When the constraints are chosen to match the empirical single-site and pairwise frequencies of the MSA, the maximum-entropy (MaxEnt) distribution reduces to the Potts model.

== Connection to statistical mechanics
The Potts model is rooted in statistical mechanics, where it was introduced as a generalization of the Ising spin model. Both models were developed to study ferromagnetism and phase transitions.

=== Ising Model
The Ising model @ising1925 describes a system of spins $sigma_i$ arranged on a lattice. Each spin can take one of two possible values:
$
  sigma_i in {+1, -1}.
$
The Hamiltonian function, which represents the energy, of a spin configuration ${sigma}$ is 
$
  H_("Ising")({sigma}) = -J sum_(angle.l i, j angle.r) sigma_i sigma_j - h sum_i sigma_i,
$
where $J$ is the coupling constant calculated over nearest-neighbor pairs (represented by $angle.l i, j angle.r$) that defines the strength between paired interactions. The $h$ term is the external magnetic field acting on the spins. At thermal equilibrium, the probability of a configuration is given by the Boltzmann distribution:
$
  P({sigma}) = 1 / Z exp(-beta H_"Ising" ({sigma})),
$
with $beta = 1 \/ (k_b T)$ and partition function
$
  Z = sum_({sigma}) exp(-beta H_"Ising" ({sigma})).
$

=== Potts Model
The Potts model @potts1952 generalizes the Ising model by allowing each spin to take on $q$ possible states. Its Hamiltonian is 
$
  H_("Potts") ({sigma}) = -J sum_(angle.l i, j angle.r) delta(sigma_i, sigma_j).
$
where $J$ is the ferromagnetic coupling that encourages neighboring spins to align in the same state. While the original formulation did not include external fields, they can be incorporated into the Hamiltonian to allow for site-dependent biases. The corresponding partition function is
$
  Z = sum_({sigma}) exp(-beta H_"Potts" ({sigma})).
$
For $q=2$, the Potts model reduces to the Ising model. The probability distribution over spin configurations is then given by
$
  P(sigma) = 1 / Z exp(-beta H_"Potts" ({sigma})).
$

It is important to note that the Potts model is over-parametrized: multiple choices of $(h, J)$ yield the same probability distribution. As a result, the parameters are not uniquely identifiable without fixing a gauge. This ambiguity can make the norms of $J$ misleading and complicate optimization. Several gauge conventions are commonly used--such as the reference-state gauge, the zero-sum gauge, or implicit definitions through regularization--which will be discussed in detail in the following sections. In the context of protein sequence analysis, the Potts model 

== Direct Coupling Analysis
The statistical dependencies observed in MSAs encode important information about protein structure and function. However, naive measures such as covariance cannot distinguish between direct correlations, which arise from true physical or functional contacts, and indirect correlations, propagated via other residues. Direct Coupling Analysis (DCA) addresses this problem by introducing a method to separate direct correlations.

DCA assumes that sequences in the MSA are independent realizations of a Potts model distribution:
$
  P(bold(sigma)) = 1 / Z exp(sum_i h_i (sigma_i) + sum_(i<j) J_(i j) (sigma_i, sigma_j))
$ <potts>
where $h_i$ are site-specific fields and $J_(i j)$ are coupling parameters between residue pairs. The central goal of DCA is to infer the interaction parameters $J_(i j)$ that best explain the observed single- and pairwise frequencies. In practice, direct inference of these parameters is computationally challenging. Evaluating the partition function $Z$ is intractable for realistic proteins, as it requires summing over $q^L$ possible sequences. Thus, distinct implementations of DCA introduce approximations to bypass this computation. The original message passing implementation @weigt2009, directly attempted to solve for couplings through a slowly converging iterative scheme, making it inapplicable in practice. In this paper, we will explore some important algorithmic innovations in DCA following its conception.

= Evolution of DCA Methods

== Mean-field DCA
The mean-field Direct Coupling Analysis (mfDCA) algorithm, introduced by Morcos et al. @mfDCA, provided the first computationally feasible approximation of the Potts model. The idea is to approximate a small-coupling expansion, valid when correlations between sites are weak which reduces the inference problem to the inversion of a correlation matrix.

=== Method
To begin, raw frequency counts in an MSA are affected by sampling bias, since closely related sequences contribute redundant information. To mitigate this, a sequence reweighting procedure is applied. Each sequence $sigma^((a))$ is assigned a weight
$
  w_m = 1 / lr(|{b in{1, ..., M}:"sim"(sigma^((a)),sigma^((b))) >= x}|),
$ <weight>
where M is the number of sequences in the MSA, and $"sim"(dot,dot)$ denotes the pairwise sequence similarity #footnote[The original paper used "seqid" for percentage identity. We adopt the notation "sim" as the future methods used this notation.]. Intuitively, if a sequence is similar to many others, with a similarity greater than the threshold $x$, its contribution is down-weighted; if it is unique, it retains full weight. Through empirical tests $x=0.8$ was found to be most effective. The weight of a sequence $m$ is $w_m$, and the effective number of  sequences is
$
  M_("eff") = sum_(m=1)^M w_m.
$ <M_eff>
From these weights, the new regularized empirical single and pairwise frequencies are computed as
$
  f_i (A) &= 1 / (M_("eff") + lambda) (lambda / q + sum_(m=1)^M w_m delta(A, sigma_i^((m)))) \
  f_(i j)(A, B) &= 1 / (M_("eff") + lambda) (lambda / q^2 + sum_(m=1)^M w_m delta(A, sigma_i^((m))) delta(B, sigma_j^((m))))
$
where $q=21$ is the amino acid alphabet with the gap symbol, and $lambda$ is a pseudocount parameter used for regularization. To address the overparametrization of the Potts, mfDCA employs a reference-state gauge, which removes redundancy among parameters while preserving all physically relevant couplings. Typically, it is set as the last amino acid $A=q$. Implementing this gauge gives us
$
  forall i,j: J_(i j)(a, q) = J_(i j)(q, a) = h_i (q) = 0
$
*Coupling Inference: Mean-field Approximation* \
To circumvent the intractable computation of the partition function $Z$, mfDCA assumes that the correlations between sites are weak. Under this approximation, the exponential of the Potts defined in @potts can be expanded to first order in a Taylor series  @plefka1982 @georges1991. This yields a relation between the couplings and the connected correlation matrix:
$
  C_(i j)(A,B)=f_(i j) (A,B) - f_i (A)f_j (B).
$ <corr>
In this framework, the couplings can be approximated through a naive mean-field inversion,
$
  e_(i j) (A, B) = - (C^(-1))_(i j)(A,B),
$
where $C$ is treated as a $((q-1)L) times ((q-1)L)$ matrix, and the pair $(i, A)$ is regarded as a single index #footnote[A detailed derivation is provided in @app2]. In practice, however, the correlation matrix is often singular or ill-conditioned, making the inversion unstable. To address this, mfDCA opts to use a strong pseudocount ($lambda approx M_"eff"$) ensuring numerical stability and preventing spurious large couplings. 

*Pair Scoring: Direct Information* \
Once the pairwise couplings $e_(i j) (A,B)$ are inferred, residue pairs must be ranked by their interaction strength. The $(q-1) times (q-1)$-dimensional coupling matrices need to be mapped to a single scalar parameter. A natural starting point is mutual information (MI), which measures statistical dependence between two MSA columns. However, MI is not able to disentangle direct correlations. Therefore a new measure, direct information (DI) @weigt2009 was introduced. To begin, for each pair $(i, j)$ a two-site model is constructed:
$
  P_(i j)^"(dir)" (A,B) = 1 / Z_(i j) exp( e_(i j) (A,B) + tilde(h)_i (A) + tilde(h)_j (B)),
$
with auxiliary fields $tilde(h)_i, tilde(h)_j$ chosen such that the single-site marginals are preserved:
$
  sum_(B=1)^q P_(i j)^"(dir)" (A,B) = f_i (A), quad sum_(A=1)^q P_(i j)^"(dir)" (A,B) = f_j (B).
$
The direct information is then defined as the mutual information of this two-site model:
$
  cal(S)^"DI"_(i j) = sum_(A,B=1)^q P_(i j)^"(dir)" (A,B)ln (P_(i j)^"(dir)" (A,B) )/ (f_i (A) f_j (B))
$
Residue pairs with the highest DI scores are predicted to be structural contacts.

=== Limitations
While mfDCA represented a breakthrough in usability of DCA, the simplifying approximations impose limitations on the model:
- Weak coupling assumptions: the small-coupling expasion assumes nearly linear correlations, which can underestimate strong epistatic effects in proteins @Haldane_2019.
- Computational scaling: the inversion of the correlation matrix scales as $cal(O)(L^3)$, which is costly for very large MSAs.
- Pseudocount dependence: due to the algorithm's vast parameter size (in the order of $400N^2$) strong regularization is required. The choice of pseudocount $lambda$ can significantly affect performance.

== Pseudo-likelihood Maximization DCA
While mfDCA provided a fast approximation, it relied on strong assumptions. Pseudolikelihood maximization (plmDCA) relaxed these by fitting sitewise conditional probabilities directly. Concretely, the inverse Potts problem on an alignment of length $L$ becomes $L$ coupled multinomial logistic regressions rather than a single optimization over the global partition function.

=== Method
For plmDCA, the Potts parametrization of fields and couplings introduced in @potts is retained. Similar to mfDCA, sequences are reweighted using $w_m$ defined in @weight and the effective sample size $M_"eff"$ defined in @M_eff.

*Coupling Inference: Pseudolikelihood Maximization* \
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
The $ell_2$ penalty fixes the gauge implicitly by selecting a unique representative among gauge-equivalent parameters. Each $g_r$ is precisely the loss of a multinomial logistic regression (softmax): classes are the $q$ states of $sigma_r$; features are one-hot encodings of ${sigma_i}_(i!=r)$ with $(L-1)(q-1)$ degrees of freedom after dropping a reference state. Hence (asymmetric) plmDCA is implemented as $L$ independent, weighted softmax problems with $ell_2$ penalty. Two variants are used in practice: asymmetric PLM, which fits each independently and then symmetrizes averaging:
$
  hat(J)_(i j) <- 1/2(J_(i j)^((i))+J_(i j)^((j)))
$
and symmetric (joint) PLM, which minimizes $cal(L)_"pseudo"$ over all parameters at once. Results are typically similar, but optimization in the symmetric variant can be heavier.

*Pair Scoring: Corrected Frobenius Norm* \
For pair scoring, plmDCA avoids the previously defined Direct Information (DI) as it would introduce a pseudocount regularization parameter, on top of the existing $lambda_h "and" lambda_J$. Instead, a new score is defined starting from the Frobenius Norm (FN)
$
  norm(J_(i j))_2 = sqrt(sum_(A,B=1)^q J_(i j)(A,B)^2).
$
Unlike the DI score, the FN is not independent of gauge choice so a direct decision must be made. As noted in @weigt2009, the zero-sum gauge minimizes the FN, making it the most appropriate gauge choice available. The procedure is hence 
(i) convert to zero-sum gauge $J'_(i j)$: 
$
  J'_(i j) (A, B) = J_(i j) (A, B) - J_(i j) (dot, B) - J_(i j) (A, dot) + J_(i j) (dot, dot),
$ <jprime>
where "$dot$" denotes a simple average over the $q$ states at that position; (ii) compute the Frobenius norm 
$
  cal(S)_(i j)^"FN" = norm(J'_(i j))_2 = sqrt(sum_(A,B=1)^q (J'_(i j)(A,B))^2);
$
(iii) apply Average Product Correction (APC) adjusted from its use in @jones2012psicov to reduce phylogenetic bias,
$
  cal(S)_(i j)^"CFN" = cal(S)_(i j)^"FN" - (cal(S)_(i dot)^"FN" cal(S)_(dot j)^"FN") / cal(S)_(dot dot)^"FN",
$
where $cal(S)_(i dot)^"FN"$ and $cal(S)_(dot j)^"FN"$ are row/column means and $cal(S)_(dot dot)^"FN"$ is the grand mean. 
Residue pairs $(i, j)$ are ranked by $cal(S)_(i j)^"CN"$ to predict structural contacts. 
=== Limitations
Despite its impact for practical applications, plmDCA remains sensitive to sampling: accurate contact recovery still requires large $M_"eff"$, and sparse or biased MSAs degrade estimates. Phylogenetic and positional bias persist (reweighting and APC help but do not eliminate them), which can inflate false positives.


=== Fast plmDCA
To make plmDCA deployable at scale, two of the original authors and a collaborator revisited its optimization strategy @fastplmDCA. The revised method adopts an asymmetric variant, instead of optimizing a single global pseudolikelihood, it solves $L$ independent, weighted softmax regressions--one per site. The resulting couplings are then symmetrized by averaging, 
$
  J_(i j) = 1/2 (J_(i j)^((i))+J_(i j)^((j))).
$
This decomposition reduces the dimensionality of each optimization problem and, crucially, enables trivial parallelization across CPU cores or nodes--the primary source of runtime gain. Since each $J_(i j)$ is regularized twice (once in the regression for $i$ and once for $j$), the coupling penalty parameter must be halved relative to the symmetric formulation. Before averaging, the two coupling blocks are shifted into the zero-sum gauge to ensure consistency. For contact prediction, the method adopts the Corrected Frobenius Norm (CFN) scoring function defined above. The combination of asymmetric regression and symmetrization achieves accuracy comparable to the original plmDCA, but at a fraction of the runtime, making large protein families and long sequences tractable in practice.

== Boltzmann Machine DCA
In mfDCA accuracy was traded for speed via the small-coupling inversion; in plmDCA, a product of sitewise conditionals was optimized, avoiding the global partition function. bmDCA @bmDCA instead proposes a solution closer to the statistical ideal: fitting the full Potts model by maximizing the true likelihood, using Monte Carlo Markov Chains to estimate intractable expectations. The result is a generative model that (by construction) matches the empirical one- and two-site statistics of the reweighted MSA.

=== Method
bmDCA aims to infer the full Potts model defined as @potts that reproduces the single- and pair-wise frequencies of the MSA. The training procedure follows the classical Boltzmann machine learning algorithm:
+ Sampling step. For current parameters {J, h}, estimate model frequencies $p_i, p_(i j)$ by MCMC
+ Update step. Adjust parameters when the estimated frequencies deviate from empirical ones: $Delta h_i prop f_i - p_i , quad  Delta J_(i j) prop f_(i j) - p_(i j)$
+ Iterate until empirical and model moments match within tolerance.

Regularization is imposed to stabilize the fit and avoid overfitting.

=== Limitations
- Computational cost: Direct Boltzmann machine learning with MCMC is extremely slow for realistic protein lengths ($10^7-10^9$ parameters).
- Sampling difficulty: Accurate moment estimation requires long chains and careful equilibration, making the approach impractical without approximations.
- Scaling: While bmDCA is theoretically optimal (full likelihood optimization), it is limited in practice to small systems where large-scale computation is feasible.

=== adabmDCA
adabmDCA @adabmDCA the same goal of fitting the full Potts by likelihood maximization, but introduces an adaptive MCMC procedure to make the approach more practical. Instead of relying on fixed sampling parameters, it dynamically tunes equilibration, waiting times, and chain persistence. adabmDCA pushes the bmDCA approach closer to practical usability, at the cost of a more complex implementation and still significant compute requirements. 

*Coupling Infernce: Full Likelihood and Moment Matching* \
The maximum-entropy distribution that reproduces chosen moments is the Potts model defined in @potts. To reduce phylogenetic redundancy, each sequence is assigned a weight $w_m$ as in @weight. Reweighted empirical frequencies, smoothed pseudocount $lambda$ are defined as
$
  f_i (A) = (1 - lambda) f_i^"data" (A) + lambda / q, quad f_(i j)(A,B) &= (1- lambda)f_(i j)^"data" (A, B)+ lambda/q^2
\
  f_i^"data" (A) = 1/M_"eff" sum_m w_m delta(A, sigma_i^((m))), quad f_(i j)^"data" (A, B) &= 1/M_"eff" sum_m w_m delta(A, sigma_i^((m))) (B, sigma_j^((m)))
$
The average log-likelihood of the MSA under $P(dot|J,h)$ is
$
  cal(L)(J, h) = 1/M sum_(m=1)^M [sum_i h_i (sigma_i^((m))) + sum_(i<j) J_(i j) (sigma_i^((m)), sigma_j^((m)))]-log Z(J, h).
$
As an exponential-family model, $cal(L)$ is concave in the natural parameters, so gradient ascent converges to the unique optimum. The gradients are moment gaps:
$
  (partial cal(L)) / (partial h_i (A)) = f_i (A) - p_i (A), quad (partial cal(L)) / (partial J_(i j) (A,B)) = f_(i j) (A,B) - p_(i j) (A, B),
$
where $p_i$ and $p_(i j)$ are model marginals under current $(J, h)$. Hence the updates
$
  h_i^(t+1)(A) <- h_i^t (A) + eta_h [f_i (A) - p_i^((t))(A)], \
  J_(i j)^(t+1)(A, B) <- J_(i j)^t (A, B) + eta_J [f_(i j) (A, B) - p_(i j)^((t))(A, B)].
$
drive the model toward the maximum likelihood-solution where $f=p$.
The obstacle is that the marginals $p_i$ and $p_(i j)$ are not analytically computable at scale; adabmDCA estimates them by MCMC at each epoch.

*Estimating Model Expectations via Adaptive MCMC* \
At training epoch $t$, run $N_s$ independent Markov chains using Metropolis-Hastings @Metropolis1953 @mh-algorithm (Gibbs sampling strategy may also be used @gibbs-algorithm), each producing $N_c$ samples after an equilibration period $T_"eq"$ and with an inter-sample waiting time $T_"wait"$. The Monte Carlo estimators are
$
  p_i^((t)) (A) = 1/(N_s N_c) sum_(m=1)^(N_s N_c) delta(A, sigma_i^((m))(t)), \ p_(i j)^((t)) (A, B) = 1/(N_s N_c) sum_(nu=1)^(N_s N_c) delta(A, sigma_i^((m))(t)) delta(B, sigma_j^((m))(t)).
$
Chains may be transient, reinitialized every epoch, or persistent, initialized only at the first epoch. Equilibration is often sped up through persistence of chains. For more details on the adaptive scheme, refer to @app3.

A practical convergence proxy is the difference between the empirical and the model two-site connected correlations:
$
  epsilon_c = max_(i,j,a,B) abs(c_(i j)^"model" (A,B)-c_(i j)^"emp" (A,B)) quad"where" \ quad c_(i j)^"model" = p_(i j)-p_i p_j, quad c_(i j)^"emp" = f_(i j)-f_i f_j
$
with a target $epsilon_c approx 10^-2$ In addition to this, some other commonly used diagnostics is the Pearson correlation between $c^"model"$ and $c^"emp"$, one- and two-site fitting errors, and optionally, a third connected correlation on a subset of triples is used to assess generative fidelity beyond pairwise constraints. 

*Regularization* \
A fully connected Potts model contains $~ (L(L-1))/2 q^2 + L q$ parameters, in the order of $10^7$-$10^9$ for realistic $L$. Due to the finite sample size of MSAs, they rarely contain enough independent information to estimate all of the parameters robustly. Not controlling this uncertainty could lead to overfitting, high variance or instability, and bad conditioning in the model. adabmDCA addresses this with two strategies: regularization and sparsity priors.
First, parameter penalties are added via $ell_1$ and $ell_2$ norms:
$
  R_(cal(l)_1) (J, h) = theta_(1, h) sum_i norm(h_i)_1 + theta_(1, J) sum_(i<j) norm(J_(i j))_1 \
  R_(cal(l)_2) (J, h) = theta_(2, h) sum_i norm(h_i)_2^2 + theta_(2, J) sum_(i<j)norm(J_(i j))_2^2
$
$ell_2$ reduces variance and improves conditioning by smoothly shrinking all parameters. It selects a unique gauge and tends to preserve the relative ordering of strong couplings. On the other hand, $ell_1$ induces sparsity by zeroing weak parameters, also reducing overfitting, though at a cost of biasing small effects downward. Generally, in stochastic settings, elastic-net is used (a combination of both parameters), that stabilizes training near zero. In addition, a separate parameter is used for fields and couplings, typically regularizing $J$ more strongly than $h$.

Second, structural sparsity is imposed, reflecting that only a small fraction of residue pairs form contacts. Encoding a structural prior can reduce variance and speed up learning. There are two approaches:
+ A priori topology. Reduce the number of parameters by starting from a restricted edge set, for example pairs with high MI, and learn only those $J_(i j)$, omitting the rest.
+ Information-based decimation. Starting dense and iteratively removing the least informative couplings until target sparsity. This directly controls overfitting by only keeping the parameters that actually affect the model's predictions.
In short, pruning and decimation prevent overfitting by removing parameters that don't materially alter the model, and have the added benefit of aligning with the biological prior of contact sparsity.

== Autoregressive Network DCA
Building on the pseudolikelihood maximization of plmDCA, and Boltzmann-machine parametrization of bmDCA, arDCA was introduced as a faster and more efficient alternative @Trinquier2021. While plmDCA is limited in scope and bmDCA suffers from slow MCMC sampling, reformulating the problem as an autoregressive network provides some key improvements. The idea is to decompose the joint sequence probability distribution into its conditionals, thereby turning the problem of inference into a supervised learning task. In this way, arDCA emerges as a powerful model for extracting structural and functional information from the rapidly growing protein databases.

=== Method
In arDCA, the exponential-family MaxEnt distribution underlying previous DCA methods is factorized into conditional probabilities, predicting each residue from its predecessors. This follows directly from the chain rule of probability:
$
  P(bold(sigma)) = product_(i=1)^L P(sigma_i|sigma_1,...,sigma_(i-1))
$ <ar>
Inspired by approaches in classical @wu2019solving and quantum @sharir2020deep statistical mechanics, the conditional distribution is parametrized as:
$
  P(sigma_i | sigma_(i-1), ..., sigma_1) = exp(h_i (sigma_i)+sum_(j=1)^(i-1) J_(i j)(sigma_i, sigma_j)) / (sum_sigma_i exp(h_i (sigma_i)+sum_(j=1)^(i-1) J_(i j)(sigma_i, sigma_j)))
$<ardca>
This is a multiclass softmax regression @hastie2009elements, generalizing logistic regression to multiple residue states. The model is defined by site-specific fields $h_i (A)$ and directed couplings $J_(i j) (A,B)$, and was therefore termed arDCA. For comparison, the authors also introduced a simpler profile model (independent-site model), which only includes field terms. Its joint probability factorizes across sites,
$
P(sigma_1,...,sigma_L) = product_(i=1...L)f_i (sigma_i)
$.
making it a useful baseline.

Although arDCA shares the same number of parameters as standard DCA, their interpretation differs. In the Potts formulation used by standard DCA, couplings are symmetric, $J_(i j)(a, b) = J_(j i)(b, a)$. In contrast, arDCA employs directed couplings, describing the influence of site $j$ on $i$ only for $j<i$. Thus, only the lower-triangular part of the coupling matrix is populated. Inference is carried out via pseudo-likelihood maximization @balakrishnan2011learning, allowing exact gradients to be computed from data--unlike bmDCA, which requires costly MCMC sampling. Moreover, plmDCA conditions each residue on all others in the sequence and then symmetrizes the couplings to align with Potts models. This symmetrization shifts parameters away from their maximum-likelihood values and reduces generative accuracy. arDCA avoids this step, preserving likelihood consistency.

The chain rule decomposition is valid under any site ordering, but once the conditional parametrization of @ardca is adopted, the likelihood becomes order dependent. Optimizing over $L!$ permutations is infeasible in practice, so the original authors proposed entropic ordering as a principled heuristic, ranking sites from lowest entropy (most conserved position) to the highest entropy (most variable). This choice reflects a natural interpretation: conserved positions provide little additional information if occupied by their dominant residue, whereas deviations must be compensated downstream by more variable sites. By placing these variable positions later in the order, the model can capture compensatory effects in a structured way. Although we do not optimize or reorder sites in our implementation, entropic ordering illustrates how  positional dependencies can be learned by arDCA. 

*Coupling Inference: Likelihood Maximization* \
The inference of the parameters is done through likelihood maximization. Following a Bayesian setting with a uniform prior, the optimal parameters are those that maximize the probability of the data:
$
  {J^*, h^*} &= arg max_{J, h} P(cal(M)|{J, h}) \
  &= arg max_{J, h} log P(cal(M)|{J, h}) \
  &= arg max_{J, h} sum_(m=1)^M log product_(i=1)^L P(sigma_i^m|sigma_(i-1)^m,...,sigma_1^m) \
  &= arg max_{J, h} sum_(m=1)^M sum_(i=1)^L log P(sigma_i^m|sigma_(i-1)^m,...,sigma_1^m)
$ <lhmax>
Each $h_i (A)$ and $J_(i j)(A,B)$ is present in only one conditional probability $P(sigma_i|sigma_(i-1),...,sigma_1)$, thus we can maximize each conditional probability independently in @lhmax:
$
  {J_(i j)^*, h_i^*} = arg max_{J_(i j), h_i} sum_(m=1)^M [h_i (sigma_i^m) + sum_(j=1)^(i-1) J_(i j) (sigma_i^m, sigma_j^m) - log z_i (a_(i-1)^m, ..., a_1^m)]
$
where
$
  z_i (sigma_(i-1), .., sigma_1) = sum_(sigma_i) exp(h_i (sigma_i)+sum_(j=1)^(i-1) J_(i j)(sigma_i, sigma_j))
$ <cpnorm>
is the normalization factor of the conditional probability of $sigma_i$. Taking the derivative with respect to $h_i (A)$ or $J_(i j) (A,B)$, with $j=1, ..., i-1$, we get:
$
  0 &= 1/M sum_(m=1)^M [ delta(A, sigma_i^m) - (partial log z_i (sigma_(i-1)^m, .., sigma_1^m)) / (partial h_i (A))] \ 
  0 &= 1/M sum_(m=1)^M [ delta(A, sigma_i^m) delta(B, sigma_j^m) - (partial log z_i (sigma_(i-1)^m, .., sigma_1^m)) / (partial J_(i j) (A,B))].
$
Using @cpnorm, we find
$
  (partial log z_i (sigma_(i-1)^m, .., sigma_1^m)) / (partial h_i (A)) &= P(sigma_i = A |sigma_(i-1)^m,...,sigma_1^m) \ 
  (partial log z_i (sigma_(i-1)^m, .., sigma_1^m)) / (partial J_(i j) (A,B)) &= P(sigma_i = A |sigma_(i-1)^m,...,sigma_1^m) delta(sigma_j^m, B).
$
The set of equations reduces to a simple form:
$
  f_i (A) &= angle.l P(sigma_i = A |sigma_(i-1)^m,...,sigma_1^m) angle.r_cal(M), \
  f_(i j) (A, B) &= angle.l P(sigma_i = A |sigma_(i-1)^m,...,sigma_1^m) delta(sigma_j^m, b) angle.r_cal(M),
$
where $angle.l dot angle.r_cal(M) = 1/M sum_(m=1)^M dot^m$ denotes the empirical average. The first variable, $i=1$, is unconditioned, therefore the coupling equation is $J equiv 0$
and the equation for the field is the profile model, $h_1(a) = log f_1(a) + "const"$.

Unlike bmDCA, the equations do not enforce exact matching between model marginals and empirical frequencies. Thus the ability to reproduce the frequencies is a good proxy for the generative properties of the model, on top of the fitting quality of current parameters. In practice, the parameters are updated using gradient descent on the likelihood. These gradients are exact as we can take the expectations directly over the MSA. For regularization, arDCA makes use of $ell_2$ penalties, with different strengths depending on whether the task is generating sequences or contact prediction. It was noted that small regularization values improved the generative quality while larger were beneficial for contact prediction.

=== Key Contributions

arDCA represents a major advance over previous DCA methods by enabling exact gradient computations from the data. This eliminates the need for costly MCMC sampling, leading to training speeds that are two to three orders of magnitude faster while maintaining, or even improving, predictive accuracy @Trinquier2021. The resulting model is therefore lightweight and computationally scalable, making it feasible for application on large protein families. For the largest families tested, although arDCA was most time efficient, the deep learning model DeepSequence @deepsequence prevailed in prediction accuracy.

A second key contribution is that arDCA allows for the exact calculation of sequence probabilities, a task intractable for previous models. In bmDCA, only unnormalized sequence weights could be computed, with the partition function requiring expensive thermodynamic integration. By contrast, each conditional probability in arDCA is normalized locally. This reduces the computational burden from summing over $q^L$ possible sequences to only $L$ sums over q residue states. This feature enables direct sequence-level comparisons across models, essential for applications such as homology detection, protein family classification, and model-based sequence evaluation.

= Implementation and Exploration
The original arDCA model was developed in Julia, a high-performance language designed for numerical and scientific computing. Julia's just-in-time (JIT) compilation, efficient handling of linear algebra, and native support for parallelism, make it well-suited for implementing large-scale statistical models. For broader accessibility and integration with standard machine-learning workflows, we re-implemented the model in Python.

*Important Differences Between Julia and Python* \
Several structural differences complicate a direct translation between Julia and Python. Julia is column-major and 1-indexed, while Python is row-major and 0-indexed; these distinctions require careful adjustment of array reshaping, indexing, and loop boundaries to maintain correctness and efficiency. Julia loops are compiled to fast machine code and can avoid bounds checking with the `@inbounds` macro, whereas Python loops are inherently slow and require vectorization through NumPy or broadcasting for comparable performance. Finally, Julia's JIT compilation introduces an initial overhead but accelerates subsequent calls, while Python relies on interpretation and delegates heavy computation to optimized C/C++ or Fortran backends.

Julia favors low-level efficiency with straightforward loop-based implementations while for Python, a vectorized style is encouraged, relying on well-developed numerical libraries. Our Python re-implementation therefore required significant restructuring of the original Julia code to align with Python's capabilities without sacrificing performance.

== Implementation details
We implement the autoregressive network within the PyTorch ecosystem, using the $mono("nn.Module")$ interface to leverage its training and optimization capabilities. The central operation of the model is the computation of the autoregressive logits, which define the conditional distribution at each sequence position. 

*Logit Computation* \
For site $i$, the logit vector is given by
$
  z[m,i,a] = h[i,a] + sum_(j<i) sum_b J[i,j,a,b] X[m,j,b]
$ 
where $h$ are the local biases and $J$ the pairwise couplings, and  $X$ is a tensor holding the one-hot encoding of the symbol $b$ observed at site $j$ in sequence $m$.

To compute the autoregressive logits, our initial implementation employed a masked matrix multiplication, performed in a single step using the `einsum` operator. While correct, and more efficient than naive loops, calculating the full interaction matrix is computationally inefficient because the upper-triangular entries are all masked to zero. To address this, we devised a more efficient formulation that exploits the block-sparse structure of the coupling matrix. Specifically, we collect all lower-triangular index pairs $(i, j)$ with $j<i$, extract the corresponding coupling blocks $J[i, j]$, and accumulate their contributions through an `einsum` operation. The resulting implementation only computes terms required by the autoregressive factorization, removing unnecessary computation.

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

  # compute pairwise contributions and accumulate into logits
  contrib = torch.einsum("mpq,pqr->mpr", X_blocks, J_blocks)
  contrib = torch.einsum("mpb,pab->mpa", X_blocks, J_blocks) 
  logits = logits.index_add(1, i_idx, contrib)
  return logits
```
With this formulation, the computational complexity scales with the number of lower-triangular pairs rather than the full $L times L$ coupling matrix. For long sequences, this is especially beneficial. In addition, by avoiding the construction of the full masked interaction matrix, the memory footprint is roughly halved, from $L^2$ couplings to $L(L-1) / 2$.

Throughout the model, an explicit binary mask `J_mask` is maintained that encodes the lower-triangular structure, to ensure consistency in the autoregressive aspect. After each optimization step, unused entries of $J$ are clamped to zero.

*Ancestral Sampling*\
Once model parameters $(h, J)$ are inferred, sequences can be generated by ancestral sampling. The autoregressive factorization in @ar allows residues to be sampled sequentially, one site at a time. The procedure is:
+ Draw the first residue from $P(sigma_1)$
+ For each subsequent site $i$, compute the conditional distribution
  $
  P(sigma_i|sigma_1,...,sigma_i-1)
  $
  given the residues already sampled.
+ Sample $sigma_i$ from this distribution.
+ Repeat until all $L$ positions are generated.
Since each conditional distribution only has $q=21$ states, the sampling step is fast. The entire sequence generation scales linearly with sequence length $L$.
```Python
def sample_sequences(self, n_samples: int = 1) -> torch.Tensor:
  self.eval()
  L, q = self.L, self.q
  samples = torch.zeros((n_samples, L), dtype=torch.long)
  for pos in range(L):
      if pos == 0:
          logits = self.h_pos[0].unsqueeze(0).expand(n_samples, -1)
      else:
          partial_oh = torch.zeros((n_samples, L, q), device=device)
          for i in range(pos):
              partial_oh[torch.arange(n_samples), i, samples[:, i]] = 1.0
          full_logits = self.compute_ar_logits(partial_oh) # (n_samples, L, q)
          logits = full_logits[:, pos, :]
      probs = F.softmax(logits, dim=-1)
      samples[:, pos] = torch.multinomial(probs, 1).squeeze(-1)
  return samples
```
The python implementation above iteratively constructs one-hot encodings of previously sampled positions, computes logits for the current site via the previously defined autoregressive logits function, converts them to probabilities through a softmax, and samples the next residue via multinomial sampling.

*Training and Evaluation Details* \
The model is initialized with the local field of the first position, $h[0]$, estimated from the empirical frequencies observed at site 0. The rest of the parameters are drawn from small random distributions. Training incorporates sequence reweighting to correct for redundancy in MSAs. The weights for each sequence are computed before model training, with parameter `theta` controlling the sequence similarity threshold. The model can be optimized with either L-BFGS @LiuNocedal1989, suited for energy-based models, or AdamW @AdamW, offering a more scalable and robust option for large datasets. Regularization is applied separately to the field and couplings via $ell_2$-norm penalties with hyperparameters $lambda_h$ and $lambda_J$.

For evaluation, the model reports the negative log-likelihood (NLL) per position,
$
  "nll"_(m,i) = -w_m log p_theta (X_(m,i)|X_(m,<i))
$
with the average NLL per position 
$
  overline("NLL") = 1 / (M L) sum_(m=1)^M sum_(i=1)^L "nll"_(m,i).
$
Another metric used for evaluation is the perplexity, 
$
  "Perplexity" = exp(overline("NLL")).
$
Perplexity is widely used in Natural Language Processing and represents the average number of equally likely choices the model is making per residue @perplexity. Lower perplexity means the model is more confident and better fits the data. The effective sample size is also reported, to make the results more easily comparable across different datasets.

== Experiment on the effect of gaps in multiple sequence alignments
When constructing MSAs, gaps are introduced to better align homologous positions across sequences. These gaps serve two main purposes: (i) ensuring residues align correctly at homologous sites, and (ii) extending sequences to a specific length for models like ArDCA, that require fixed length input. Since the presence of gaps can influence the statistical properties of the alignment, we examined how restricting the fraction of gaps affects results. Specifically, we introduced two filtering parameters:
- `max_gap_fraction`: sets the maximum allowed fraction of gaps per sequence. 
- `max_col_gap_fraction`: sets the maximum allowed fraction of gaps in a column.
We will observe how the $M_"eff"$ is affected by removing rows and columns, as well as some key metrics like average per-site Perplexity, and average per-site train and validation NLL.  For the first experiment, we varied `max_gap_fraction`, while holding the other fixed, with three thresholds: $0.05, 0.1, "and" 1$. The last threshold does not exclude any sequences, therefore is used as a baseline for comparison. The difference in training time of the models was also recorded, but given the small nature of the change, no conclusive results were observed. 

#align(center)[
#block(width: 100%,)[
  #figure(
    table(
      columns: (auto, auto, auto, auto, auto, auto, auto),
      table.vline(x: 1, start: 1),
      table.header[version][`max_gap_fraction`][M][$M_"eff"$][Perplexity][Train NLL][Val NLL],
      [0], [1], [13600], [*4248.8*], [1.72], [0.32], [0.54],
      [1], [0.1 (-1.6%)], [13376], [4145.4], [*1.58*], [0.32], [*0.46*],
      [2], [0.05 (-6%)], [12784], [3876.5], [1.68], [0.3], [0.51]  
    ),
    caption: [Perplexity and NLL improve when moderate gap filtering (0.1) is applied]
  )
]]
Restricting the allowed fraction of gaps removed sequences therefore reduced the effective number of sequences $M_"eff"$. The strictest filtering criterion ($0.05$), removed approximately 6% of sequences relative to the baseline, yielding a smaller dataset and just slightly better perplexity value. The lowest perplexity was observed when 10% of gaps were excluded. A perplexity of $1.58$ means the model is, on average, choosing between fewer than 2 equally likely residues. Although the full dataset had the highest $M_"eff"$, higher perplexity suggests that allowing too many gaps may negatively impact the quality of the alignment. A similar exploration was performed for the column gaps with values  $0.025, 0.1, "and" 1$. Again, 1 represents the baseline with no columns removed.
#align(center)[
#block(width: 100%,)[
  #figure(
    table(
      columns: (auto, auto, auto, auto, auto, auto, auto),
      table.vline(x: 1, start: 1),
      table.header[version][`max_col_gap_fraction`][L][$M_"eff"$][Perplexity][Train NLL][Val NLL],
      [0], [1], [53], [4248.8], [1.72], [0.32], [0.54],
      [3], [0.1 (-3.7%)], [51], [4208.5], [*1.71*], [0.32], [*0.53*],
      [4], [0.025 (-7.5%)], [49], [*4363.0*], [1.75], [0.34], [0.56]      
    ),
    caption: [Similar perplexity and NLL throughout the experiments]
  )
]]
In this experiment, the change in perplexity and $"M_eff"$ wasn't as drastic, though again we observe the best perplexity value in the second threshold. It is interesting that the effective sample size increased when removing the four columns of version 4. This occurred because of the change in the weight computation of the sequences. Both sequence-level and site-level filtering yield a better result when using moderate thresholds, where excluding only the data with the most gaps improves model performance. Overly aggressive filtering discards useful information and hurts perplexity, while using the full dataset also doesn't reach the best performance.

= Results and Discussion

== Training Behavior
We trained arDCA with L-BFGS, using exact gradients of the negative log-likelihood (NLL). A 10% validation split was held out, and early stopping prevented overfitting. Across all settings, the training NLL decreased smoothly and the validation NLL tracked it closely, indicating stable optimization and good generalization. A few values of the $lambda_J$ and $lambda_h$ parameters were tested but the final choices aligned with the papers recommendations of `1e-4` and `1e-6` respectively. @nll shows the runs in the fraction gap study, the NLL drops from around 1.7 to under 1.1. In practice. we observed little sensitivity to random seed for the reported settings, consistent with convex, per-site softmax problems and well-behaved $ell_2$ regularization. Exact gradients and the autoregressive factorization yield smooth, monotonic learning curves and minimal overfitting under regularization.

#figure( 
 image("code/out/train_nll.png", width: 80%),
 caption: [Training NLL plotted for baseline and gap fraction experiments]
) <nll>


== Generative Model Quality
To assess the structural plausibility of sequences generated by the model, we employed AlphaFold's predicted Local Distance Difference Test (pLDDT) score @AlphaFold_pLDDT. The pLDDT provides a residue-level confidence measure ranging from 0 to 100, reflecting the predicted accuracy of local environments. High pLDDT values (>90) correspond to very reliable backbone and side-chain placement, values between 70-90 indicate generally well-resolved domain-level structures, while values below 70 are associated with uncertain or potentially disordered regions. By evaluating generated sequences with this metric, we obtain an intrinsic estimate of whether the model captures the structural constraints encoded in the protein family. pLDDT represents a quantitative proxy for generative quality, complementing the evaluation metrics used previously such as negative log-likelihood and perplexity. It serves to provide a bridge between sequence statistics and three-dimensional structure accuracy.

#figure(grid(columns: 2, row-gutter: 2mm, column-gutter: 1mm,

  image("code/out/sample1plddt.png"), image("code/out/sample2plddt.png"), 

  "a) Sample 1", "b) Sample 2"),

  caption: "pLDDT per site for model-generated samples"

)<plDDTplot>
@plDDTplot shows the pLDDT profiles for two representative model-generated samples. Both sequences display consistently high pLDDT values across most positions, indicating that the model generates foldable and structurally coherent proteins. The slight drop in the end of a sequence is likely as these sites often correspond to poorly aligned or frequently gapped positions in the alignment. These results complement our previously observed statistical metrics by linking generative performance to structural integrity.

== Limitations and Future Work
The Python implementation, while optimized through a lower-triangular masked multiplication, remains less efficient than the original Julia version due to the differences in performance between the languages. GPU acceleration could substantially improve speed and scalability, enabling the model to handle longer, or more populated, protein families. From an architectural perspective, arDCA is restricted to fixed length sequences, limiting its ability to model families with variable domain lengths. Furthermore, our study used the natural site ordering exclusively; prior work has demonstrated that entropy-based permutations can yield superior likelihood estimates that improve the quality of the model. Exploring alternative orderings could yield more expressive models with better structural accuracy. Another important limitation is the strong dependence of arDCA on the quality of the multiple sequence alignment. Families with limited diversity or strong phylogenetic bias may provide insufficient or misleading signals. In such cases, the model risks poor generalization, generating sequences that fail to capture the true structural and functional constraints of the protein family. Careful sequence reweighting and gap filtering could mitigate these effects. Our analysis was constrained by evaluation throughput. The pLDDT metric from AlphaFold provides valuable structural insight but running it for each sample is computationally expensive, limiting us to a small number of generated samples. A more rigorous statistical analysis would require folding a larger set of generated sequences and comparing them to natural ones.

Future work should thus focus on shifting the heavy computations to be done on the GPU, investigating alternative site orderings, such as the entropic order heuristic, and exploring more expressive generative architectures such as transformer-based models @Caredda2024.02.06.579080. Finally, the study could benefit from a large-scale structural benchmarking through AlphaFold's pLDDT. Collectively, these improvements would make a more scalable, flexible, and biologically accurate model, further bridging the gap between statistical models and natural protein sequences.

= Conclusion
In this work, we demonstrated that arDCA, originally implemented in Julia, can be successfully adapted in Python as a practical tool for protein sequence analysis. By leveraging PyTorch and introducing an efficient block-sparse formulation for computing autoregressive logits, our implementation reduces computational cost and memory usage. Empirical evaluation shows smooth optimization dynamics, improvements in negative log-likelihood and perplexity, and structural plausibility as measured by pLDDT. Our experiments on gap filtering highlight the importance of sequence quality control, where moderate gap filtering improves model confidence. In the future, this work could benefit from GPU acceleration to improve its speed and allow it to work more efficiently with large datasets. Overall, this work adapted arDCA beyond its Julia origin and produced a powerful protein-generating model, opening the doors for downstream tasks like novel protein generation and structural modeling.


*Code availability*
The full Python code is available at https://shorturl.at/5ufHk.

#pagebreak()
#bibliography("references.bib")

#pagebreak()
#show: appendix

= Langrage Multipliers Calculation for Maximum Entropy Principle <app1>
We seek the probability distribution ${p_i}$ over stattes $x_i$ that maximizes the Shannon entropy subject to the normalization and expectation value constraints. Introduce multipliers $lambda$ and $mu$:
$
  Phi = - sum_i p_i ln p_i - lambda (sum_i p_i-1)-mu (sum_i p_i f(x_i)-angle.l f angle.r)
$
From here, we can set the stationarity condition
$
  (partial Phi) / (partial p_i) = -(ln p_i +1)-lambda-mu f(x_i) =0,
$
rearranging:
$
  ln p_i = -1 -lambda - mu f(x_i).
$
We then solve for $p_i$, absorbing the -1 constant in $lambda$,
$
  p_i = exp( - lambda -mu f(x_i)) \
$
Impose the normalization condition:
$
  1 = sum_i exp( - lambda -mu f(x_i)), quad "and factor our " e^(-lambda) \
  
  1= e^(-lambda) sum_i e^( -mu f(x_i)).
$
Define the partition function:
$
  Z(mu) = sum_i e^( -mu f(x_i))", thus" \
  e^(-lambda)Z(mu) = 1 quad => quad lambda = ln Z(mu)
$
Thus, the maximum-entropy distribution is
$
  p_i = 1 / (Z(mu)) exp(-mu f(x_i)).
$
The multiplier $mu$ is determined by the constraint:
$
  angle.l f angle.r = sum_i p_i f(x_i) = - partial / (partial mu) ln Z(mu).
$
Finally, inserting this into the entropy formula,
$
  S_max = ln Z(mu) + mu angle.l f angle.r 
$
giving us the canonical exponential-family distribution.
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
This equation allows us to solve the original inference in the mean-field approximation in a single step. To determine the marginals of the empirical frequencies, we just need to determine the empirical connected correlation matrix and invert it to get the couplings $e_(i j)$. 

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