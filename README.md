# Internship Notebook

## Day 1: 10/04/2017
###	Objectives > Status:
- Install software > Completed
- Read papers > Pending
- GitHub account > Completed
 https://education.github.com/
- Read about Markdown > Completed

###	Notes on all Objectives:
- Install software: I now have a running bash shell which I used to make a directory called 'Internship'. This is going to be the root folder for all files relevant for this internship. 
- Read Papers: I was able to read 3 of the 6 recommended articles. 1 was particularly difficult to understand ( Finotello F, Di Camillo B. Measuring differential gene expression with RNA-seq: challenges and strategies for data analysis. Brief Funct Genomics. 2015 Mar;14(2):130-42. doi:10.1093/bfgp/elu035. Review. PubMed PMID: 2524000) because it had many new terms and programming jargons.
-  GitHub account: I opened a github account under the name theNNvader and also registered under the github education platform.
- Read about Markdown: This was a little confusing at first but I got help from Sung-Hyan and Richa . I need to learn the syntax but I understand the concept.

[theNNvader]: https://github.com/theNNvader 

# Day 2: 11/04/2017
### objectives > Status
- Read Papers > Completed
- Create Markdown File for Internship Notebook > Completed
- Download the Genome of Salmonella Thyphimurium SL1344 > Completed
- Count the occurance of all feautures (genes, CDS, rRNA, tRNA) in the file. Just use standard libraries. > Completed
- Calculate the average gene length > Pending
- Calculate the average gene legth of each feature class > Pending
- Clean up your code - write functions > Pending
- use main function and further function to structure you code > Pending
- Rewrite the code using object oriented programming > Pending
- Rewrite the code and use pandas > Pending
- Generate a histogram of the length distribution using matplotlib > Pending
- Generate a histogram of the length distribution using bokeh > Pending

### Codes
#### Count the occurance of all feautures (genes, CDS, rRNA, tRNA) in the file
- first try:*This code worked but gave inflated values because it scanned the entire document.
```
file = 'GCF_000210855.2_ASM21085v2_genomic.gff'
f_o = open(file)
cds_count = 0
gene_count = 0
rRNA_count = 0
tRNA_count = 0
for lines in f_o:
    if lines.startswith('#'): continue
    attr = lines.split()
    for ele in attr:
        if ele == 'CDS':
            cds_count = cds_count + 1
        if ele == 'gene':
            gene_count = gene_count + 1
        if ele == 'rRNA':
            rRNA_count = rRNA_count + 1
        if ele == 'tRNA':
            tRNA_count =tRNA_count + 1
print(cds_count, gene_count, rRNA_count, tRNA_count)
```
- with a dictionary:*This code did the same thing as the one above
```
file  = 'GCF_000210855.2_ASM21085v2_genomic.gff'
gff_dic = dict()
o_f = open(file)
for lines in o_f:
    if lines.startswith('#'):continue
    words = lines.split()
    for key in words:
        if key not in gff_dic:
            gff_dic[key] = 1
        else:
            gff_dic[key] = gff_dic[key] + 1
print(gff_dic['gene'])
print(gff_dic['CDS'])
print(gff_dic['tRNA'])
print(gff_dic['rRNA'])
```
- Richa Improved:*This code only scanned the keys from the second column hence gave lower values.
```
file  = 'GCF_000210855.2_ASM21085v2_genomic.gff'
gff_dic = dict()
o_f = open(file)
for lines in o_f:
    if lines.startswith('#'):continue
    words = lines.split("\t")
    gff_key = words[2]
    if gff_key not in gff_dic:
        gff_dic[gff_key] = 1
    else:
        gff_dic[gff_key] = gff_dic[gff_key] + 1
        
print(gff_dic['tRNA']) 
```
# Day 3: 12/04/2017
### Objective > Status
- Calculate the average gene length > Completed
- Calculate the average gene legth of each feature class > Completed
- Install R studio > Completed
- Clean up your code - write functions > Pending
- use main function and further function to structure you code > Pending
- Rewrite the code using object oriented programming > Pending
- Rewrite the code and use pandas > Pending
- Generate a histogram of the length distribution using matplotlib > Pending
- Generate a histogram of the length distribution using bokeh > Pending

### Codes
#### Calculate the average gene length
```
file  = 'GCF_000210855.2_ASM21085v2_genomic.gff'
gff_dic = dict()
o_f = open(file)
l_g_l = []
total_gene_length = 0
for lines in o_f:
    if lines.startswith('#'):continue
    words = lines.split("\t")
    start_length = int(words[3])
    end_length = int(words[4])
    gene_length = end_length - start_length
    l_g_l.append(gene_length)
total_gene_length = sum(l_g_l)
num_gene_length = len(l_g_l)
average_gene_length = total_gene_length/num_gene_length
print(average_gene_length) 
```
#### Calculate the average gene length of each feature class
```
file  = 'GCF_000210855.2_ASM21085v2_genomic.gff'
gff_dic = dict()
gff_dic_count = dict()
feature_dic = dict()
averag_gene_length = 0
o_f = open(file)
for lines in o_f:
    if lines.startswith('#'):continue
    words = lines.split("\t")
    gff_key = words[2]
    start_length = int(words[3])
    end_length = int(words[4])
    gene_length = end_length - start_length
    if gff_key not in gff_dic:
        gff_dic[gff_key] = gene_length
        gff_dic_count[gff_key] = 1
    else:
        gff_dic[gff_key] = gff_dic[gff_key] + gene_length
        gff_dic_count[gff_key] = gff_dic_count[gff_key] + 1
        
for values in gff_dic:
    for value in gff_dic_count:
        if values == value:
            average = gff_dic[values]/gff_dic_count[value]
        feature_dic[values] = average
print(feature_dic)
```
