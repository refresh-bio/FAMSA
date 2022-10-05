from Bio import SeqIO
from Bio.Seq import Seq
import ete3


for i in [1,2]:

    file = open(f'upgma.part{i}.dnd', 'r')
    tree_str = file.read()
    tree_str = tree_str[:-1].replace(';','___') + ';'
    tree = ete3.Tree(tree_str)
    file.close()

    nodes = { node.name.replace('___',';') for node in tree.traverse("postorder")}

    records = []
    for seq in SeqIO.parse('upgma.no_refine.fasta', 'fasta'):
        if seq.id in nodes:
            records.append(seq)

    length = len(records[0].seq)
    to_remove = set()

    for c in range(length):
        gaps_only = True
        for r in records:
            if r.seq[c] != '-':
                gaps_only = False
                break

        if gaps_only == True:
            to_remove.add(c)

    for ir in range(len(records)):
        seq = records[ir].seq
        temp = Seq(''.join([ seq[c] for c in range(length) if c not in to_remove ]))
        records[ir].seq = temp


    SeqIO.write(records, f'upgma.part{i}.fasta', 'fasta')


