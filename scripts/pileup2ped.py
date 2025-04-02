#!/usr/bin/python3
__author__ = 'Daniel Fernandes'

import os
import sys
import random
import re

cwd = os.getcwd()

# Parse arguments
if len(sys.argv) < 2:
    print("Usage: pileup2ped.py <excludeTerminals> [--sample_id=<sample_id>]")
    sys.exit(1)

exclTerminals = sys.argv[1]
sample_id = None

# Check for the --sample_id flag
for arg in sys.argv[2:]:
    if arg.startswith("--sample_id="):
        sample_id = arg.split("=")[1]

# Ensure a sample_id is provided
if not sample_id:
    print("Error: --sample_id flag is required to specify the sample.")
    sys.exit(1)

# Locate the target file for the specified sample
target_file = f"{sample_id}.pileupsamtools.gwv.txt"
if not os.path.exists(target_file):
    print(f"Error: File {target_file} not found for sample_id {sample_id}.")
    sys.exit(1)

# Process the specified sample
pedout = open(sample_id + '.ped', 'w')
pileupin = open(target_file, 'r')
snpcount, discard = re.split(' ', os.popen('wc -l ' + target_file).read())
snpcount = int(snpcount)
pedout.write(sample_id + ' ' + sample_id + ' 0 0 0 1 ')
counter = 1
countbads = 0

for line in pileupin:
    if counter <= snpcount:
        line = line.rstrip()
        if line.split("\t")[3] == "0":
            pedout.write('0 0 ')
            countbads += 1
        else:
            chro, coord, refb, cov, rbase, qual = line.split('\t')
            counter += 1
            rbase = rbase.replace(".", refb.capitalize())
            rbase = rbase.replace(",", refb.capitalize())
            if int(cov) == 1:
                if r"*" in rbase or r"+" in rbase or r"-" in rbase:
                    pedout.write('0 0 ')
                elif r"$" in rbase:
                    rbase = rbase[0]
                    if exclTerminals == 'TRUE': 
                        pedout.write('0 0 ')
                    else:
                        pedout.write(rbase.capitalize() + ' ' + rbase.capitalize() + ' ')
                elif r"^" in rbase:
                    rbase = rbase[2]
                    if exclTerminals == 'TRUE':
                        pedout.write('0 0 ')
                    else:
                        pedout.write(rbase.capitalize() + ' ' + rbase.capitalize() + ' ')
                else:
                    pedout.write(rbase.capitalize() + ' ' + rbase.capitalize() + ' ')
            elif int(cov) > 1:
                listBases = list(rbase)
                while r"-" in listBases:
                    for i in listBases:
                        if i == r"-":
                            indexi = listBases.index(r"-")
                            delL = listBases[indexi+1]
                            del listBases[indexi:(indexi+int(delL)+2)]
                            del listBases[indexi-1]
                while r"+" in listBases:
                    for i in listBases:
                        if i == r"+":
                            indexi = listBases.index(r"+")
                            delL = listBases[indexi+1]
                            del listBases[indexi:(indexi+int(delL)+2)]
                            del listBases[indexi-1]
                while r"^" in listBases:
                    for i in listBases:
                        if i == r"^":
                            indexi = listBases.index(r"^")
                            del listBases[indexi]
                            del listBases[indexi]
                            if exclTerminals == 'TRUE' and len(listBases) > 0:
                                del listBases[indexi-1]
                while r"$" in listBases:
                    for i in listBases:
                        if i == r"$":
                            indexi = listBases.index(r"$")
                            del listBases[indexi]
                            if exclTerminals == 'TRUE' and len(listBases) > 0:
                                del listBases[indexi-1]
                if len(listBases) != 0:
                    chosenbase = random.choice(listBases)
                    pedout.write(chosenbase.capitalize() + ' ' + chosenbase.capitalize() + ' ')
                else:
                    countbads += 1
                    pedout.write('0 0 ')
    else:
        pass
pedout.close()
print(f"[INFO] Processing completed for sample: {sample_id}")

