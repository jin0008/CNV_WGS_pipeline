# coding: utf-8

# Author: Julien Buratti

import os
import glob
import pysam

sample_l = []
sex_d = {}
bams = glob.glob("*.analysisready.bam")

print("\nSex determination script openning.\n")

for i in bams:
    sample = os.path.basename(i.split(".")[0])
    sample_l.append(sample)
    bamfile = pysam.AlignmentFile(i, "rb")
    sry_count = bamfile.count(contig='Y', start=2786989, stop=2787603, until_eof=False, read_callback='all')
    if sry_count >= 50:
        sex = "M"
    elif sry_count <= 10:
        sex = "F"
    else:
        sex = "?"
        sys.exit("/!\\ Ambigous Sex Determination, please check!")    
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
