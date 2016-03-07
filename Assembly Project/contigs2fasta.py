l = []
with open("k25_sample_fasta_contigs.txt",'r') as f:
    counter = 0
    for row in f:
        counter += 1
        text = ">contig_%s\n" % counter
        text += row
        print text
        l.append(text)

with open("k25_sample_fasta_annotated_contigs.txt","w") as g:
    for item in l:
        g.write(item)
