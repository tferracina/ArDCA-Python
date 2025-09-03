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

#outline(target: heading.where(supplement: [Appendix]), title: [Appendix])


= Introduction

Proteins are biomolecules that are fundamental to nearly all biological processes. Their diverse roles include transporting nutrients, catalyzing chemical reactions, providing structural support, and more. The function of a protein is determined by the composition of its primary sequence, composed of 20 amino acids. A single protein sequence can vary greatly in length and order of its amino acids, leading to a very large number of possible configurations. 

The size of the protein sequence space makes exhaustive experimental exploration infeasible. However, the rapid growth of biological databases, driven by advances in sequencing technologies, has transformed biology into a data-rich discipline. Large repositories such as UniProt, Pfam, and the Protein Data Bank store millions of sequences and structures, providing crucial resources for computational approaches@weigt2020.

This wealth of data has enabled the development of advanced statistical and machine learning models capable of simulating protein sequence evolution, prediciting structural conformations, and generating novel sequences with desired properties. In addition, breakthroughs in protein structure prediction--most notably the Nobel Prize winning AlphaFold 2--demonstrated that computational methods can rival experimental accuracy, in a much cheaper manner.

In this work, I explore the implementation of an autoregressive network for protein sequence generation leveraging the Direct Coupling Analysis (DCA) method to model residue-residue dependencies. In Section 2, the biological background will be laid out. Section 3 will introduce the mathematical foundations. Section 4 will detail previous iterations of similar models, while Section 5 will explore my implementation (and improvements).

= Biological Background

== Proteins and their structure 
Proteins are essential biological molecules responsible for a wide range of functions in living organisms. Despite their functional diversity, all proteins are polymers of the same set of standard building blocks: the 20 canonical amino acids arranged in different assortments.

Amino acids share a common core structure consisting of a central carbon atom bonded to a hydrogen atom, an amino group, a carboxyl group, and a variable side chain. This side chain is the defining feature of each amino acid, giving rise to differences in size, shape, chemical reactivity, and polarity. The distinct amino acids are bonded together through peptide bonds to form proteins, also known as polypeptides. In this process, certain atoms are lost, and what remains of each amino acid is called a residue. Thus, within a protein sequence, individual amino acids are typically referred to as residues. Generally, protein sequences are made up of between 50 and 2000 amino acids. The ordering of the amino acids dictates how a protein folds into its three-dimensional structure, known as its conformation. A protein's conformation defines its function. Although each conformation is unique, two common folding patterns occur in many proteins: the $alpha$ helix and the $beta$ sheet. @alberts2002

Protein structure is broken down into four levels of organization. The amino acid sequence is known as the primary structure. The regularly repeating local structures stabilized by chemical bonds, such as the previously mentioned helices and sheets, make up the secondary structure. The overall three-dimensional shape of a single protein molecule constitutes the tertiary structure. Finally, if a protein molecule is formed as a complex of more than one polypeptide chain, the complete structure is known as the quaternary structure @alberts2002.

Each amino acid is chemically distinct and can occur at any position in a protein chain, giving rise to $20^n$ possible polypeptide sequences of length $n$. For a typical protein of 300 amino acids, the number of possible sequences is astronomically large. However, only a small fraction of these sequences are capable of folding into a stable three-dimensional conformation. Natural selection has enabled living organisms to explore sequence space, favoring those sequences that reliably fold into stable structures.

== Protein families and evolutionary information
Proteins do not evolve in isolation; they often belong to protein families, groups of proteins that share a common evolutionary origin, therefore exhibiting related sequence features and functional properties @ebi_protein_families. Over evolutionary timescales, mutations accumulate in these families: some are beneficial, altering the protein activity in ways that give rise to new functions, while many others are neutral and have no effect on stability or activity. Harmful changes, by contrast, disrupt folding or function and are eliminated by natural selection. The result is a collection of homologous proteins that retain overall structural and functional characteristics, but also display sequence variability that encodes the evolutionary history of the family @alberts2002. 

The study of evolutionary history and relationships among biological entities is referred to as phylogenetics. Protein families provide the domain for phylogenetic analysis, as examining families provides insight that cannot be obtained from a single sequence. Analyzing homologous proteins across diverse organisms allows the detection of correlated mutations between amino acid positions, which in turn represent structural or functional constraints enforced by evolution. These statistical patterns are exploited by computational approaches to predict three-dimensional structure and understand protein function.

== Multiple sequence alignments
To extract important evolutionary clues from protein families, the homologous sequences need to be organized in a systematic way. This is done through a multiple sequence alignment (MSA), in which three or more sequences are arranged so that homologous sites are placed in the same column @WILTGEN201938. To maximise the positional correspondence of sequences with varied length, alignment gaps are introduced when necessary. 

Multiple sequence alignments reveal patterns of conservation and variation across the family. Conserved positions, those which are unchanged in multiple sequences, represent sites that are critical for maintaining structure or function, while variable positions indicate sites that can tolerate mutations without disruption to the protein. Beyond conservation, MSAs also capture covariation: pairs of positions that mutate in a correlated way across sequences. These covariation signals reflect couplings, where a mutation at one site requires compensation at another to maintain protein integrity. @adhikari2016


== Protein contacts
One of the earliest challenges for researchers in computational biology--historically encouraged by the CASP (Critical Assessment of Structure Prediction) competition--has been to understand how the linear sequence of amino acids folds into its conformation. One of the ways researchers did this was by exploring pairs of residues that are spatially close in a protein's folded three-dimensional structure, called contacts. As the structures can be represented in Cartesian coordinates $(x,y,z)$, the contacts are defined using distance thresholds. Two residues are generally considered to be in contact if the distance between selected atoms is below a set threshold, usually 8 Ångströms (Å). Residues that are far apart in the linear sequence, may be close together in the 3D structure @adhikari2016.

Contacts are further broken down into short, medium, and long range predictions. Most computational approaches evaluate the long range contacts separately, as they are the most important for accurate predictions and unsurprisingly, the hardest to predict.

== The problem of inference
The specific challenge in computational biology that will be explored in this paper is the prediction and generation of protein sequences given a multiple sequence alignment. As previously mentioned, the protein families contain correlations between residues and we want to build a model to take advantage of that structure. The  covariance is confounded by phylogeny and sampling bias @dietler2023, but most importantly, the inherent challenge of this problem is to distinguish between the direct and indirect correlations. 

= Mathematical Foundations

== Proteins as statistical systems
Protein sequences can be thought of as random variables produced by a certain distribution.  Each sequence of length $L$ can be written as:
#set math.equation(numbering: "(1)", supplement: "Eq.")
$
bold("S") = (s_1, s_2, ..., s_L), quad s_i in cal(A),
$ where $cal(A)$ is the alphabet of size $q = 21$ (20 amino acids and the alignment gap) and $s_i$ are the residue sites. When considering multiple sequence alignments, we have $M$ empirical samples: ${s^((1)) , s^((2)), ..., s^((M))}$.
The sites $s_i$, columns in the MSA, are random variables and the whole sequence is therefore a joint random variable $bold("S")$ with a distribution over $cal(A)^L$, where $|cal(A)| = 21$.

From the alignments, we can define single and pairwise frequency counts for the columns. The single-site frequency for MSA column $i$ can be computed as:
$
f_i (A) = 1/M sum^M_(m=1) delta(s_i^((m)), A), quad A in cal(A),quad  
"where" delta "is the Kronecker delta" #footnote[$delta(x, y) := cases( 0 "if" x != y, 1 "if" x=y) $]
$
While the pairwise frequencies of MSA columns $i, j$ are computed as:
$
  f_(i j) (A, B) = 1/M sum^M_(m=1) delta(s_i^((m)), A) delta(s_j^((m)), B)
$

These empirical frequencies will serve as constraints for the model we want to infer.

== Maximum Entropy Principle
To find the probability distribution $P(bold(S))$ that will satisfy constraints, in other words reproduce the empirical marginals $f_i (A)$ and $f_(i j) (A, B)$, we can use the maximum entropy principle@MITinfent.

The first step in setting up the MEP, is to extract information from the system. Usually, this information is given in the form of averages of functions $angle.l f(x) angle.r$. For example, in a physical system, one could compute average magnetization or energy of the observed system. We also need to define a probability of occupancy of states, $p(x)$, which runs over all the possible states. This distribution has the usual property of mapping each state to a value within 0 and 1 and adding up to 1 when considering all states. 

Our uncertainty on the system is expressed quantitatively through Shannon entropy $S$@Shannon1948:
$
  S = - sum_x p(x) ln p(x) 
$ <entropy>
Generally, the distribution which maximizes entropy is the uniform distribution. In this situation, constraints affect the probabilities of states, so the uniform is not suitable. The constraints on our system are:
$
  sum_x p(x) = 1, quad "and" quad sum_x p(x) f_k (x) = angle.l f_k angle.r quad (k = 1, ..., m).
$ <constraints>

By selecting the Shannon entropy as our measure of information, it allows us to use Lagrange Multipliers for our maximization problem@Jaynes1957. To maximize @entropy subject to the constraints @constraints, we introduce Langrange multipliers $lambda_0, lambda_1, ... lambda_m$ which yield:
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
The MaxEnt distribution with single and pairwise constraints is exactly the Potts model.

== Connection to statistical mechanics
The Potts model is defined in the context of statistical mechanics, where it was introduced as a generalization of the Ising spin model. Both models were defined as a way to explore ferromagnetism and phase transitions.

=== Ising Model
The Ising model describes a system of spins $sigma_i$ arranged on a lattice. Each spin can take one of two possible values:
$
  sigma_i in {+1, -1}
$
The Hamiltonian function, which represents the energy, of a spin configuration ${sigma_i}$ is 
$
  H_("Ising")({sigma}) = -J sum_(angle.l i, j angle.r) sigma_i sigma_j - h sum_i sigma_i,
$
where $J$ is the coupling constant that defines the strength between paired interactions, where the sum is calculated over nearest-neighbor pairs on the lattice (represented by $angle.l i, j angle.r$). The $h$ term is the external magnetic field acting on the spins.

At thermal equilibrium, the probability of a configuration is given by the Boltzmann distribution:
$
  P({sigma}) = 1 / Z e^(-beta H({sigma})),
$
with $beta = 1 \/ (k_b T)$ and partition function
$
  Z = sum_({sigma}) e^(-beta H({sigma})).
$

=== Potts Model
The Potts model generalizes the Ising model by allowing each spin to take on $q$ possible states. The $q$-state Potts model spin can take values in ${1, 2, ... q}$. The Hamiltonian is written as
$
  H_("Potts") ({sigma}) = -J sum_(angle.l i, j angle.r) delta(sigma_i, sigma_j).
$
$J$ again represents the ferromagnetic coupling that encourages neighboring spins to align in the same state. The partition function becomes
$
  Z = sum_({sigma}) exp(beta J sum_(angle.l i, j angle.r) delta(sigma_i, sigma_j)).
$

The Potts model simplifies to an Ising model when $q=2$.


Putting this together, the probability distribution from the Potts model is written as:
$
  P(bold(S)) = 1 / Z exp(sum_(i=1)^L h_i (s_i) + sum_(1 <= i < j <= L) J_(i j)(s_i, s_j)).
$

== Direct Coupling Analysis
Naively, correlations in the alignment can be capture by covariance, but this simple approach will not be able to separate the direct correlations, arising from structural or functional contacts, and the indirect correlations, which are propagated via other residues. Taking the empirical frequencies as before we can define the covariance matrix:
$
  C_(i j)(A,B)=f_(i j) (A,B) - f_i (A)f_j (B).
$ <corr>
A positive $C_(i j)(A, B)$ means that $A$ at $i$ and $B$ at $j$ co-occur more often than expected by chance. The reverse is true for a negative value, where the residues occur less often than expected. The covariance matrix is computed across the whole alignment, and it contains crucial information about pairwise correlations between residues.

From here, all classical DCA methods use the empirical frequency counts as constraints under which the Maximum Entropy Principle will find the most accurate and least biased distribution. As previously alluded, the solution to this system is the Potts model:

$
  P(bold(S)) prop exp(sum_(i) h_i (s_i) + sum_(i < j) J_(i j)(s_i, s_j))
$
where the local field represents the single-column, and the couplings represent the column pairs.  

The problem with this approach is that learning all of the $h_i, J_(i j)$ directly is difficult. To find them, an evaluation of the partition function $Z$ is required, which entails summing over the potentionally astronomical $q^L$ possible sequences. The distinct implementations of DCA introduce approximations to bypass this complex computation.

The first implementation of DCA was done through a message passing algorithm. This algorithm was computationally costly as it was based on a slowly converging iterative scheme@weigt2009. In this paper, we will explore the innovations in models using the DCA framework following its conception.

= Evolution of DCA Methods

== Mean-field DCA (2011)
The mean-field Direct Coupling Analysis (mfDCA) algorithm, introduced by Morcos et al. (2011), provides a computationally feasible approximation of the Potts model used to disentangle direct from indirect correlations in multiple sequence alignments. mfDCA's method includes reweighing the sequences, using maximum entropy formulation of the distribution, and a small-coupling expansion to reduce the inference problem to the inversion of the correlation matrix@morcos2011.

=== Method
The authors found the raw frequency counts to be suffering from sampling bias, so the weight of highly similar sequences was reduced. Each sequence $A^a$ is assigned a weight
$
  m^a = |{b in{1, ..., M}|"seqid"(A^a,A^b) > 80%}|,
$
where M is the number of sequences in the MSA and seqid is their percentage identity. The effective weight of a sequence $a$ is $1 \/ m^a$, and the total effective number of sequences is
$
  M_("eff") = sum_(a=1)^M 1 \/ m^a.
$
From these weights, the regularized empirical single-site frequencies are computed as
$
  f_i(A) = 1 / (M_("eff") + lambda) (lambda / q + sum_(a=1)^M 1 / m^a delta(A, A_i^a)) 
$
and the pairwise frequencies as
$
  f_(i j)(A, B) = 1 / (M_("eff") + lambda) (lambda / q^2 + sum_(a=1)^M 1 / m^a delta(A, A_i^a) delta(B, A_j^a))
$
where $q=21$ is the amino acid alphabet with the gap symbol, and $lambda$ is a pseudocount parameter used for regularization.

The maximum entropy principle is applied to reproduce the empirical single- and pairwise frequencies. This yields the previously derived Potts model distribution

$
  P(A_1, ... A_L) = 1 / Z exp(sum_( 1<= i < j <= L) e_(i j) (A_i, A_j) + sum_(i=1) h_i (A_i))
$ <mfequation>
with local fields $h_i(A)$ and pairwise couplings $E_(i j) (A, B)$, and partition function
$
  Z = sum_(A_1, ..., A_L) exp(sum_( 1<= i < j <= L)(A_i, A_j) + sum_(i=1) h_i (A_i))
$ <mfpartition>

Since the number of free parameters in @mfequation exceeds the number of constraints, multiple equivalent solutions can be found for the fitting. Therefore, the couplings and fields are measured relative to the reference state, typically the last amino acid $A=q$:
$
  forall i,j: e_(i j)(A, q) = e_(i j)(q, A) = 0, quad h_i (q) =0
$

The partition function $Z$ cannot be computed exactly because it requires summing over $q^L$ sequences. To make the problem tractable, mfDCA assumes weak correlations between sites, expanding the exponential  in @mfequation by the Taylor series to first order. This results in the relation found in @corr between the couplings and the connected correlation matrix.

The couplings are then approximated as
$
  e_(i j) (A, B) = - (C^(-1))_(i j)(A,B),
$
where $C$ is treated as a $((q-1)L) crossmark ((q-1)L)$ matrix, and the pair $(i, A)$ is regarded as a single index.

In practice, the correlation matrix is often singular without regularization. mfDCA opts to use a strong pseudocount ($lambda approx M_"eff"$) to stabilize the inversion and prevent spurious large couplings.

Once the couplings are inferred, the next step is to rank residue pairs by their likelihood of physical contact. For each pair $(i, j)$, a two-site model is constructed:
$
  P_(i j)^"(dir)" (A,B) = 1 / Z_(i j) exp( e_(i j) (A,B) + tilde(h_i) (A) + tilde(h_j) (B)),
$
with auxiliary fields $tilde(h_i), tilde(h_j)$ chose such that
$
  sum_B P_(i j)^"(dir)" (A,B) = f_i(A), quad sum_A P_(i j)^"(dir)" (A,B) = f_j(B).
$
The Direct Information (DI) between sites $i, j$ is then defined as the mutual information of this two-site distribution:
$
  "DI"_(i j) = sum_(A B)P_(i j)^"(dir)" (A,B)ln (P_(i j)^"(dir)" (A,B) )/ (f_i (A) f_j (B))
$
Residue pairs are ranked by $"DI"_(i j)$, and the top-scoring pairs are predicted to be structural contacts.
=== Limitations
While mfDCA represented a breakthrough in usability of these models, it relies on approximations that impose limitations:
- Weak coupling assumptions: the small-coupling expasion assumes nearly linear correlations, which can underestimate strong epistatic effects in proteins. {add in-text}
- Computational scaling: the inversion of the correlation matrix scales as $cal(O)(((q-1)L)^3)$, which is costly for very large MSAs.
- Pseudocount dependence: the algorithm requires strong pseudocount regularization, and the choice of $lambda$ significantly affects performance.

== Pseudo-likelihood Maximization DCA

=== Method
In the plmDCA, the full likelihood is replaced with conditional likelihoods:
$
  log P(bold(S)) arrow sum_i log P(s_i | bold(S_( \\ i)))
$
Each conditional distribution does not require the difficult to compute partition function $Z$

This model is widely used in practice due to its feasible and scalable optimization.

== Boltzmann Machine DCA

=== Method
Use Monte Carlo sampling from the Potts model to adjust the parameters until the model reproduces the observed single- and pairwise frequencies.

This is the most faithful approach to the original inference problem. Because of this, it is computationally expensive as each gradient step requires expensive sampling. The cost of this model makes it unsuitable for long protein sequences.

== Autoregressive Network DCA

ArDCA was built to explore the ability of generative models for protein design coming from sequence-data. ArDCA emerges as a capable model for the extraction of structural and functional protein information which encoded in rapidly growing protein databases.

=== Method
In arDCA, the exponential-family MaxEnt distribution is replaced with a conditional probability model, where each residue is predicted from the previous ones. This arises from the chain rule decomposition of the join probability distributions:
$
  P(bold(S)) = product_(i=1)^L P(s_i|s_1,...,s_(i-1))
$

The parameters are learned by predicting each residue given the previous ones, a similar approach to that taken in NLP methods.

This method has the advantage of being tractable, able to generate new sequences from the given learned parameters, and be scalable.

=== Key Contributions

This lightweight approach performs at a similar accuracy of previous iterations, but at a substantially lower computational cost (a factor between $10^2$ and $10^3$) @weigt2020. It presents an important innovation also due to its mathematical advantages, which will be explored further, leading to improved applicability in sequence generation and evaluation.

== Attention DCA - (MAYBE)
- Attention-Potts Model: factored self-attention -> potts model

= Implementation and Extension

== Implementation details
- Julia -> Python re-implementation (vectorization, frameworks, missing functions)
== Computational challenges

- Benchmarks

== (MAYBE) Improvements:
- Allowing arbitrary sequence length (GPT-style transformers)
- Incorporating attention mechanism (Potts with attention)

Evaluation with Boltz-2 / AlphaFold3

= Results and Discussion

#bibliography("references.bib")

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

