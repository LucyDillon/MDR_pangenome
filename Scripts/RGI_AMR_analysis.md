RGI AMR analysis:

# E. coli
To get the number of unique ARGs I run this command: 
directory: 
 ```
 $ cut -f9 *.faa.txt | grep -v 'Best_Hit_ARO' | sort | uniq | wc -l
 $ 465
 ```
How many drug classes does this cover? 32
```
$ cut -f15 *.faa.txt | grep -v 'Drug Class' | sed 's/; /\n/g' | sort | uniq
acridine dye
aminocoumarin antibiotic
aminoglycoside antibiotic
benzalkonium chloride
carbapenem
cephalosporin
cephamycin
diaminopyrimidine antibiotic
elfamycin antibiotic
fluoroquinolone antibiotic
fosfomycin
glycopeptide antibiotic
glycylcycline
isoniazid
lincosamide antibiotic
macrolide antibiotic
monobactam
nitrofuran antibiotic
nitroimidazole antibiotic
nucleoside antibiotic
oxazolidinone antibiotic
penam
penem
peptide antibiotic
phenicol antibiotic
pleuromutilin antibiotic
rhodamine
rifamycin antibiotic
streptogramin antibiotic
sulfonamide antibiotic
tetracycline antibiotic
triclosan
```
How many ARGs did we find all together? 559038
```
cut -f9 *.faa.txt | grep -v 'Best_Hit_ARO' | wc -l
```


Now lets work out how the genes were distributed across the pangenome (or all genomes):
- first lets get a uniq list of genes:
```
cut -f9 *.faa.txt | grep -v 'Best_Hit_ARO' | sort | uniq > uniq_ARGs.txt
```

- secondly, lets search for the ARGs in each file to see how many of the files contain the ARGs and return the file name of those with matches, we can then count these 
test:
```
grep -ol -m 1 "AAC(2')-Id" *.faa.txt | wc -l
```

- get a column just of the ARGs in a new file
```
for i in *.faa.txt; do cut -f9 $i > ${i%%.txt}.ARGs.txt; done
```

Now lets iterate over each .ARGs.txt file:
```
while read i; do   count=$(grep -l -m 1 -- "$i" *.ARGs.txt | wc -l);   echo -e "$i\t$count"; done < uniq_ARGs.txt > ARG_and_no_of_genomes_present_in.txt
```

to check this is correct, search for a gene in a low amount to go through the results physically and count instead of relying on the code


# P. aeruginosa:

How many unique ARGs? 416
```
cut -f9 *.faa.txt | grep -v 'Best_Hit_ARO' | sort | uniq | wc -l
```


How many drug classes does this cover? 30
```
$ cut -f15 *.faa.txt | grep -v 'Drug Class' | sed 's/; /\n/g' | sort | uniq
acridine dye
aminocoumarin antibiotic
aminoglycoside antibiotic
benzalkonium chloride
bicyclomycin
carbapenem
cephalosporin
cephamycin
diaminopyrimidine antibiotic
elfamycin antibiotic
fluoroquinolone antibiotic
fosfomycin
glycopeptide antibiotic
glycylcycline
lincosamide antibiotic
macrolide antibiotic
monobactam
mupirocin
nitrofuran antibiotic
oxazolidinone antibiotic
penam
penem
peptide antibiotic
phenicol antibiotic
pleuromutilin antibiotic
rifamycin antibiotic
streptogramin antibiotic
sulfonamide antibiotic
tetracycline antibiotic
triclosan
```


How many ARGs did we find all together? 426164
```
cut -f9 *.faa.txt | grep -v 'Best_Hit_ARO' | wc -l
```


Now lets work out how the genes were distributed across the pangenome (or all genomes):
- first lets get a uniq list of genes:
```
cut -f9 *.faa.txt | grep -v 'Best_Hit_ARO' | sort | uniq > uniq_ARGs.txt
```

- secondly lets create a new file with just the ARGs in (this is just to speed up the processing):
```
for i in *.faa.txt; do cut -f9 $i > ${i%%.txt}.ARGs.txt; done
```


Now lets iterate over each .ARGs.txt file:
```
while read i; do   count=$(grep -l -m 1 -- "$i" *.ARGs.txt | wc -l);   echo -e "$i\t$count"; done < uniq_ARGs.txt > ARG_and_no_of_genomes_present_in.txt
```

to check this is correct, search for a gene in a low amount to go through the results physically and count instead of relying on the code

