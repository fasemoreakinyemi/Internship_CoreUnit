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
- Calculate the average gene length
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
- Calculate the average gene length of each feature class
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
# Day 4: 13/04/2017
### Objective > Status
- Rewrite the code with pandas > Completed

### Codes
- Count the occurance of all feautures (genes, CDS, rRNA, tRNA) in the file
```
import pandas as pd
pd.set_option('display.mpl_style', 'default')
pd.set_option('display.line_width', 5000) 
pd.set_option('display.max_columns', 60)
file = pd.read_csv('GCF_000210855.2_ASM21085v2_genomic.gff', sep = "\t", comment = "#")
file['region'].value_counts()
```
- Calculate the average gene length 
```
import pandas as pd
file = pd.read_csv('GCF_000210855.2_ASM21085v2_genomic.gff', sep = "\t", comment = "#")
sum(file[file.columns[4]]-file[file.columns[3]])/len(file)
```
- Calculate the average gene length of each feature class
```
file = pd.read_csv('GCF_000210855.2_ASM21085v2_genomic.gff', sep = "\t", comment = "#")
gff_dic = dict()
avg_l = 0
fet_list = file['region'].unique()
for row in file['region']:
    if row not in gff_dic:
        v_row = file[file['region'] == row]
        avg_l = sum(v_row[v_row.columns[4]]-v_row[v_row.columns[3]])/len(v_row)
        gff_dic[row] = avg_l
print(gff_dic)
```
# Day 5: 18/04/2017
### Activties
- Introduction to R
    - Data Types, basic functions and commands, installing packages
    - Solved simple problems like creating vectors and dataframe slicing
- Installed bokeh 
```
conda install bokeh
```
- Set the enviroment variable for R.exe on command prompt
# Day 6: 19/04/2017
### Activities
- Installed DESeq2 package on R
```
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```
- Calculated the average gene lenght from the gff file with R
```
gff = read.csv(file = file.choose(), sep = '\t', comment = '#')
sum(gff[,5] - gff[,4])/(length(gff[,4]))
```
- Attempted differential gene expression analysis with data and instructions from
https://export.uppmax.uu.se/b2013006/courses/RNAseq201410/build/html/courseSource/diffexp-lab.html#deseq
# Day 7: 20/04/2017
### Activities
- Visualization using matplotlib.
    - Bar graph of feature counts in gff file
    ```
    file = pd.read_csv('GCF_000210855.2_ASM21085v2_genomic.gff', sep = "\t", comment = "#")
    r = file['region'].value_counts()
    r.plot(kind='bar')
    ```
    - Bar graph of average length of each feature
    ```
    file = pd.read_csv('GCF_000210855.2_ASM21085v2_genomic.gff', sep = "\t", comment = "#")
    gff_dic = dict()
    avg_l = 0
    xas = []
    yas = []
    fet_list = file['region'].unique()
    for row in file['region']:
        if row not in gff_dic:
            v_row = file[file['region'] == row]
            avg_l = sum(v_row[v_row.columns[4]]-v_row[v_row.columns[3]])/len(v_row)
            gff_dic[row] = avg_l
    for i in gff_dic:
        xas.append(i)
        yas.append(gff_dic[i])
    fr = pd.Series(yas, index = xas)
    fr.plot(kind='bar')
    ```
# Day 8: 21/4/2017
### Activities
- plotted stacked bar of relative mean of bacteria in abundance.summary file.
```
abs = pd.read_csv('abundance.summary', sep = '\t')
abs.index = abs["taxon"]
lvl2 = abs[abs['taxlevel'] == 2]
taxn = lvl2['taxon']
rel_dic = dict()
absco = lvl2[lvl2.columns[5:]]
for i in absco:
    total = sum(absco[i])
    if total < 0: continue
    abnd = (absco[i]/total)*100.0
    rel_dic[i] = abnd
df = pd.DataFrame.from_dict(data = rel_dic, orient = 'index' )
df.plot(kind = 'bar', stacked = True)
```
- Visualization with bokeh
    - bar chart of feature counts from gff file
    ```
    from bokeh.charts import Bar
    file = pd.read_csv('GCF_000210855.2_ASM21085v2_genomic.gff', sep = "\t", comment = "#")
    r = file['region'].value_counts()
    df = pd.DataFrame(data = r )
    p = Bar(df, legend = 'top_right'
    show(p)
    ```
# Day 9: 24/04/2017
Activities:
- In search for the perfect text editor.
    - Downloaded: Emacs, Pycharm and Vim.
# Day 10: 25/05/2017
Activities:
- Wrote an Rscript for counting gene features in a gff file.
```
args = commandArgs (trailingOnly = TRUE)
if (length(args)==0) {
	stop("Please provide a file name saved in this working directory", call. = FALSE)
}
file_open = read.csv(args[1], sep = ('\t'), comment = '#')
fet_c = table(file_open['region'])
write.table(fet_c, file = args[2], row.names = FALSE)
```
- Modified the script to also calculate the average gene length.
```
args = commandArgs (trailingOnly = TRUE)
if (length(args)==0) {
	stop("Please provide a file name saved in this working directory", call. = FALSE)
}
file_open = read.csv(args[1], sep = ('\t'), comment = '#')
feat_count = table(file_open['region'])
avg_len = sum(file_open[,5] - file_open[,4])/(length(file_open[,4]))
write.table(feat_count, file = args[2], append = FALSE, row.names = FALSE)
write(avg_len, file = args[2], append = TRUE)
```
    
  
   

