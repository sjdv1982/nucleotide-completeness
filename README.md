# nucleotide-completeness

Analyses of the completeness of RNA at the dinucleotide and trinucleotide level.

## Synopsis

We are interested in the question whether every possible RNA fragment conformation has been seen in the PDB. Here we present evidence that this is so at the dinucleotide level, but not at the trinucleotide level.

We built a reproducible pipeline that analyzes all protein-bound RNAs at the level of dinucleotide/trinucleotide fragments. For each non-redundant fragment, we determined the most similar fragment among all other RNAs. For 93.2 % of the trinucleotides and 97.3 % of the dinucleotides, the most similar fragment was within 1 A. The remaining fragments were in many cases either terminal fragments or originated from a limited number of RNAs. After eliminating these cases, 99.0 % of all dinucleotide fragments had a similar fragment within 1 A among other RNAs. This shows that essentially all known RNA fragments have been seen at least twice. We performed additional analysis that shows that this level of completeness mostly holds among different protein and RNA families. Finally, we show that fragments can be clustered almost perfectly, in the sense that essentially all of the above 99 % can be organized into 1A clusters that come from multiple RNAs.
