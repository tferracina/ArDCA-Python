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

#pagebreak()
#pagebreak()


= Introduction

Proteins are biomolecules that are fundamental to nearly all biological processes. Their diverse roles include transporting nutrients, catalyzing chemical reactions, providing structural support, and more. In alignment-based analyses, the function of a protein is determined by its primary sequence, a chain composed of 20 amino acids and an additional symbol representing an alignment gap. A single protein sequence can vary greatly in length and order of its amino acids, leading to a very large number of possible configurations. 

The size of the protein sequence space makes exhaustive experimental exploration infeasible. However, the rapid growth of biological databases, driven by advances in sequencing technologies, has transformed biology into a data-rich discipline. Large repositories such as UniProt, Pfam, and the Protein Data Bank store millions of sequences and structures, providing crucial resources for computational approaches.

This wealth of data has enabled the development of advanced statistical and machine learning models capable of simulating protein sequence evolution, prediciting structural conformations, and generating novel seuqneces with desired properties. Breakthroughs in protein structure prediction, most notably the Nobel Prize winning AlphaFold, demonstrated that computational methods can rival experimental accuracy, in a much cheaper manner.

In this work, I explore the implementation of an autoregressive neural network for protein sequence generation leveraging the Direct Coupling Analysis (DCA) framework to model residue-residue dependencies as defined in [REFER PAPER]. I then evaluate the generated sequences for structural plausibility and functional relevance using AlphaFold3 and Boltz-2.

== Biological Background

=== Proteins
Proteins are essential biological molecules responsible for a wide range of functions in living organisms. Despite their functional diversity, all proteins are polymers of the same set of 20 standard building blocks, the 20 canonical amino acids arranged in different assortments. Each amino acid shares a common core structure consisting of a central carbon atom bonded to a hydrogen atom, an amino group, a carboxyl group, and a variable side chain. This side chain is the defining feature of each amino acid, giving rise to differences in size, shape, chemical reactivity, and polarity. The chemical composition of amino acids involves carbon, hydrogen, oxygen, and nitrogen atoms, with some also containing sulfur. Based on the properties of their side chains, amino acids can be broadly classified as hydrophobic (nonpolar), hydrophilic (polar), and charged. The specific sequence of amino acids dictates how a protein folds into its three-dimensional structure, which in turn dictates its function. Even small changes in sequence can dramatically affect stability, activity, and interaction patterns.

=== Multiple Sequence Alignments
Studying proteins at the sequence level can be challenging because of the immense diveristy and complexity of possible configurations. However, comparitive analyses across related proteins, grouped into protein families, allow researchers to identify conserved positions and co-evolving residue pairs, offering insights into structural and functional constraints. This is the basis of approaches such as Direct Coupling Analysis, which seeks to uncover the statistical dependencies between amino acid positions that reflect physical contacts in the folded structure.

Protein families
@adhikari2016

Here define the LENGTH L, NUMBER M of sequences, --- this will introduce the importance of picking the correct generative model in the next section


=== Protein Residues and Contacts
One of the earliest challenges for researchers in computational biology, encouraged by the CASP (Critical Assessment of Structure Prediction) competition, has been understanding how a linear sequence of amino acids folds into a three-dimensional protein structure. Since a protein's 3D shape largely determines its biological function, predicting structure from sequence is a fundamental goal. 
Two important concepts in this context are residues and contacts:
- *Residues*: An individual amino acid within a protein sequence. During peptide bond formations, the chemically linked amino acids usually lose certain atoms, and what remains is referred to as the "residue" of the original amino acid. In structural biology, “residue” also refers to a specific position in a protein sequence.
- *Contacts*: A pair of residues that are spatially close in a protein's folded three-dimensional structure. These structures can be represented in Cartesian coordinates $(x,y,z)$, hence the contacts are defined using distance thresholds. Two residues are generally considered to be in contact if the distance between selected atoms is below a set threshold, usually 8 Ångströms (Å). Residues that are far apart in the linear sequence, may be close together in the 3D structure. may be far apart in @adhikari2016

important to add: there is short, mid and long range + most models evalute long range separately- it is the most important and hardest to predicting

This leads us to: 

=== Succesful protein folding models

AlphaFold / Boltz-2 ... will be used for  structural integrity check.


== Brief Review of the paper

ArDCA is an autoregressive network built on the basis of the Direct Coupling Analyis framework defined by Trinquier et. al. It will serve as the basis of this exploration. ArDCA was built to explore the ability of generative models for protein design coming from sequence-data. ArDCA emerges as a capable model for the extraction of structural and functional protein information which encoded in rapidbly growing protein databases. The paper compares the generative model's performance to existing solutions involving Boltzmann machines and deep generative models. It finds that this lightweight approach not only performs at a similar accuracy, but at a substantially lower computational cost (a factor between $10^2$ and $10^3$). It presents an important innovation also due to its mathematical advantages, which will be explored further, leading to improved applicability in sequence generation and evaluation.
- what it is used for
- why it is important
- what sets it apart

= Literature Review
== Deep Learning for Proteins - Historical Context & SOTA
=== Early DCA approaches
- mean-field DCA
- bmDCA
- plmDCA

=== Deep Generative models
- deep sequence
- ardca

=== Transformer-based models 
- MSA Transformer evoformer
- Attention-Potts Model: factored self-attention -> potts model

== Contextualization & Key advances
- arrange chronologically to trace evolution from statistical physics to deep interpretable neural approaches
- highlight SOTA for different tasks

= Preliminary Methods

== Direct Coupling Analysis

== arDCA - technical review

= Implementation

= Results and Conclusions

#bibliography("references.bib")
