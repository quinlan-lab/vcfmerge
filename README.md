```
python vcfmerge.py --remove-ref $GATK_VCF $SV_VCF | bgzip -c > $MERGED_VCF
```

note that output is unsorted. It is up to the caller to sort by position if needed.
