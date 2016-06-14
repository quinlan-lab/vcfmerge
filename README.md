```
python vcfmerge.py --remove-ref $GATK_VCF $SV_VCF | bgzip -c > $MERGED_VCF
```

if the input is sorted, the output will also be sorted.
