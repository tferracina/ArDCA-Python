#set text(font: "Arial", size: 12pt)
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
  #set text(13pt, weight: "regular", font: "Arial CE MT")

  #block([
    // Only display number if numbering is enabled for this heading
    #(if it.numbering != none { counter(heading).display(it.numbering) + " " })
    #smallcaps(it.body)
  ])
  #v(0.5em)
]

#set par( justify: true,
          leading: 1.4em,
          spacing: 2.6em)

/*
#pagebreak()
#pagebreak()
*/

#set heading(numbering: "1.")
#outline()


= Introduction

== Introduction

Proteins are biomolecules that are fundamental to nearly all biological processes. Their diverse roles include transporting nutrients, catalyzing chemical reactions, providing structural support, and more. In alignment-based analyses, the function of a protein is determined by its primary sequence, a chain composed of 20 amino acids and an additional symbol representing an alignment gap. A single protein sequence can vary greatly in length and order of its amino acids, leading to a very large number of possible configurations. 

The size of the protein sequence space makes exhaustive experimental exploration infeasible. However, the rapid growth of biological databases, driven by advances in sequencing technologies, has transformed biology into a data-rich discipline. Large repositories such as UniProt, Pfam, and the Protein Data Bank store millions of sequences and structures, providing crucial resources for computational approaches.

This wealth of data has enabled the development of advanced statistical and machine learning models capable of simulating protein sequence evolution, prediciting structural conformations, and generating novel seuqneces with desired properties. Breakthroughs in protein structure prediction, most notably the Nobel Prize winning AlphaFold, demonstrated that computational methods can rival experimental accuracy, in a much cheaper manner.

In this work, I explore the implementation of an autoregressive neural network for protein sequence generation leveraging the Direct Coupling Analysis (DCA) framework to model residue-residue dependencies as defined in [REFER PAPER]. I then evaluate the generated sequences for structural plausibility and functional relevance using AlphaFold3 and Boltz-2.

= Biological Background

== Proteins and their structure 
Proteins are essential biological molecules responsible for a wide range of functions in living organisms. Despite their functional diversity, all proteins are polymers of the same set of standard building blocks: the 20 canonical amino acids arranged in different assortments. These amino acids are linked though peptide bonds, which is why proteins are also known as polypeptides.

Amino acids share a common core structure consisting of a central carbon atom bonded to a hydrogen atom, an amino group, a carboxyl group, and a variable side chain. This side chain is the defining feature of each amino acid, giving rise to differences in size, shape, chemical reactivity, and polarity. [+LINKING SENTENCE HIGHLIGHTING THE FACT THAT DIFFERENT AAs PRODUCE DIFFERENT PURPOSE PROTEINS- atm unclear why amino acid chemical properties are important]. Generally the sequences are made up of between 50 and 2000 amino acids. The ordering of the amino acids dictates how a protein folds into its three-dimensional structure, known as its conformation. Although each conformation is unique, two common folding patterns occur in many proteins: the $alpha$ helix and the $beta$ sheet.

Protein structure is broken down into four levels of organization. The amino acid sequence is known as the primary structure. The regularly repeating local structures stabilized by chemical bonds, such as helices and sheets, are referred to as the secondary structure. The overall three-dimensional shape of a single protein molecule constitutes the tertiary structure. Finally, if a protein molecule is formed as a complex of more than one polypeptide chain, the complete structure is the quaternary structure @alberts2002.

Each amino acid is chemically distinct and can occur at any position in a protein chain, giving rise to $20^n$ possible polypeptide sequences of length $n$. For a typical protein of 300 amino acids, the number of possible sequences is astronomically large. However, only a small fraction of these sequences are capable of folding into a stable three-dimensional conformation. Natural selection has enabled living organisms to explore sequence space, favoring those sequences that reliably fold into stable structures.


== Protein families and evolutionary information
Proteins do not evolve in isolation; instead, they belong to protein families, groups of proteins that share a common evolutionary origin and therefore exhibit related structures and functions @ebi_protein_families. 


Define: proteins with shared evolutionary origin, similar structure/function.

Point out that protein sequences vary across species but retain conserved features critical for function.

Explain why families are important: by comparing homologous sequences across organisms, we can extract information not obvious from a single sequence.


== Multiple Sequence Alignments
Studying proteins at the sequence level can be challenging because of the immense diveristy and complexity of possible configurations. However, comparitive analyses across related proteins, grouped into protein families, allow researchers to identify conserved positions and co-evolving residue pairs, offering insights into structural and functional constraints. This is the basis of approaches such as Direct Coupling Analysis, which seeks to uncover the statistical dependencies between amino acid positions that reflect physical contacts in the folded structure.


== Protein Residues and Contacts
One of the earliest challenges for researchers in computational biology, encouraged by the CASP (Critical Assessment of Structure Prediction) competition, has been understanding how a linear sequence of amino acids folds into a three-dimensional protein structure. Since a protein's 3D shape largely determines its biological function, predicting structure from sequence is a fundamental goal. 
Two important concepts in this context are residues and contacts:
- *Residues*: An individual amino acid within a protein sequence. During peptide bond formations, the chemically linked amino acids usually lose certain atoms, and what remains is referred to as the "residue" of the original amino acid. In structural biology, “residue” also refers to a specific position in a protein sequence.
- *Contacts*: A pair of residues that are spatially close in a protein's folded three-dimensional structure. These structures can be represented in Cartesian coordinates $(x,y,z)$, hence the contacts are defined using distance thresholds. Two residues are generally considered to be in contact if the distance between selected atoms is below a set threshold, usually 8 Ångströms (Å). Residues that are far apart in the linear sequence, may be close together in the 3D structure. may be far apart in @adhikari2016

important to add: there is short, mid and long range + most models evalute long range separately- it is the most important and hardest to predicting

This leads us to: 

== The problem of inference
computational biology challenge: given MSA, oobserve correlations between residues (MI, co-variation)
but correlations can be direct vs. indirect

lead into DCA 

AlphaFold / Boltz-2 ... will be used for  structural integrity check.

= Mathematical Foundations

== Proteins as Statistical Systems

== Direct Coupling Analysis

= Evolution of DCA Methods

== Mean-field DCA (mfDCA)

- Simplified inference method; first successful applications to contact prediction

- Strengths/weaknesses

== Pseudo-likelihood Maximization DCA (plmDCA)

- More accurate & scalable inference

- Widely used in practice

== Boltzmann Machine DCA (bmDCA)

- Boltzmann learning to directly fit the Potts model

- Shown to generate functional protein variants (e.g., chorismate mutase)

== Deep Generative and Hybrid Models

- DeepSequence (variational autoencoder on MSAs)

- arDCA (autoregressive formulation, avoids MCMC)

ArDCA is an autoregressive network built on the basis of the Direct Coupling Analyis framework defined by Trinquier et. al. It will serve as the basis of this exploration. ArDCA was built to explore the ability of generative models for protein design coming from sequence-data. ArDCA emerges as a capable model for the extraction of structural and functional protein information which encoded in rapidbly growing protein databases. The paper compares the generative model's performance to existing solutions involving Boltzmann machines and deep generative models. It finds that this lightweight approach not only performs at a similar accuracy, but at a substantially lower computational cost (a factor between $10^2$ and $10^3$). It presents an important innovation also due to its mathematical advantages, which will be explored further, leading to improved applicability in sequence generation and evaluation.

- Attention-based DCA (bridging DCA and transformer attention)

== Transformer-based models (could be combined above^)
- MSA Transformer evoformer
- Attention-Potts Model: factored self-attention -> potts model

= Implementation and Extension
- Julia -> Python re-implementation (vectorization, frameworks, missing functions)
- Benchmarks
(MAYBE) Improvements:
- Allowing arbitrary sequence length (GPT-style transformers)
- Incorporating attention mechanism (Potts with attention)

Evaluation with Boltz-2 / AlphaFold3

= Results and Discussion

#bibliography("references.bib")
