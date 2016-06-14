echo "\
    672 1
      6 2
      2 3"

python vcfmerge.py test.vcf test.vcf | grep -v ^# | cut -f 1 | uniq -c
