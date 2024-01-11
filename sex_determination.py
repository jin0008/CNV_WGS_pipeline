# coding: utf-8

# Author: Julien Buratti
#!/usr/bin/env python3

import os
import glob
import pysam
import sys

sample_l = []
sex_d = {}
bams = glob.glob("*.dedup.bam")

print("\n************************************")
print("Sex determination script openning.")
print("************************************\n")

for i in bams:
    
    sample = os.path.basename(i.split(".")[0])
    sample_l.append(sample)
    bamfile = pysam.AlignmentFile(i, "rb", check_sq=False)
    sry_count = bamfile.count(contig='chrY', start=2786989, stop=2787603, until_eof=False, read_callback='all')
    print(sry_count)


    if sry_count >= 20:
        sex = "M"
    elif sry_count <= 10:
        sex = "F"
    else:
        sex = "?"
    
    sex_d[sample] = sex

if os.path.isfile("samples.txt"):
    print("Sex determination was already done. No changes.")
else: 
    with open('samples.txt', 'w') as samp:
        samp.write("sample\tsex\n")
        
        for s in sample_l:
            samp.write(s + "\t" + sex_d[s] + "\n")
        
        print("Sex determination done for {} samples.".format(str(len(sample_l))))
        print("samples.txt generated.")

print("\nSex determination script job done!\n")
