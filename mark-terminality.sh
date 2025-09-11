for motif in AA AC CA CC; do
    outfile=output/dinuc-$motif-terminality.txt
    python3 mark-terminality.py dinuc $motif > $outfile
done

for motif in AAA AAC ACA ACC CAA CAC CCA CCC; do
    outfile=output/trinuc-$motif-terminality.txt
    python3 mark-terminality.py trinuc $motif > $outfile
done
