```
python vcfmerge.py --remove-ref $GATK_VCF $SV_VCF | bgzip -c > $MERGED_VCF
```

the input files must both be sorted
