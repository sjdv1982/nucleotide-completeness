# nucleotide-completeness

Analyses of the completeness of RNA at the dinucleotide and trinucleotide level.

## Instructions

This repo contains code in submodules. To check them out, do:

`git submodule update --init --recursive --progress`

## Synopsis

We are interested in the question whether every possible RNA fragment conformation has been seen in the PDB. Here we present evidence that this is so at the dinucleotide level, but not at the trinucleotide level.

We built a reproducible pipeline that analyzes all protein-bound RNAs at the level of dinucleotide/trinucleotide fragments. For each non-redundant fragment, we determined the most similar fragment among all other RNAs. For 93.2 % of the trinucleotides and 97.3 % of the dinucleotides, the most similar fragment was within 1 A. The remaining fragments were in many cases either terminal fragments or originated from a limited number of RNAs. After eliminating these cases, 99.0 % of all dinucleotide fragments had a similar fragment within 1 A among other RNAs. This shows that essentially all known RNA fragments have been seen at least twice. We performed additional analysis that shows that this level of completeness mostly holds among different protein and RNA families. Finally, we show that fragments can be clustered almost perfectly, in the sense that essentially all of the above 99 % can be organized into 1A clusters that come from multiple RNAs.

We have made an online analysis tool that can analyze the completeness at the family level (Pfam or Rfam). Using this tool, we found that completeness levels are not much lower at the family level than at the single PDB level. Concretely, one can exclude from a fragment library all fragments that are implicated in binding to the RRM family, and find that nearly all those fragments are still very well modeled by (at least one fragment of) the library.

Finally, we analyze the relation between dinucleotide and trinucleotide fragments. At 0.5 A precision, the PDB is 77 % complete; in other words, the primary 0.5A trinucleotide library covers 77 % of the fragments within 0.5. An additional 16 % of the fragments is within 1A, i.e. a total of 93 %. It turns out that this primary library can be described very well (99 % with both fits <0.5A) as *pairs* of primary (i.e. seen in more than one PDB) 0.5A *dinucleotide* fragments. These dinucleotide pairs have an excellent compatibility RMSD (97 % < 0.5A). In conclusion, *primary trinucleotides (that cover 90 % of the cases) can be synthetically generated from superimposed pairs of primary dinucleotides*. It is then observed that the pairing is *very sparse*: among all potential pairs of primary dinucleotides, only about 1 in 1000 are actually observed. Yet, many potential pairs are *plausible*: in terms of compatibility RMSD, 1 out of 18 pairings are under 0.5A. One could ascribe this sparsity to a simple lack of observations, but statistical analysis shows otherwise: a number of pairings between dinucleotide fragments that are common and plausible are nevertheless not observed, and this is shown to be highly statistically significant. We provide a small "anti-fragment library" of trinucleotides where the probability of absence is <1 % chance to be due to sparse observations. This probability of absence was calculated after removing redundancies at the pentanucleotide level, which reduced the number of observations by about half.  
