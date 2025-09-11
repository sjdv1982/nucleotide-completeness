for lib in dinuc trinuc; do
    python3 failure-analysis.py $lib > failure-analysis-$lib.out
done
