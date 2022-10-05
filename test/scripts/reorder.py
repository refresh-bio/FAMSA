from Bio import SeqIO
from Bio.Seq import Seq


# reorder part
inputs = [seq for seq in SeqIO.parse('adeno_fiber', 'fasta')]
outputs = {seq.id : seq for seq in SeqIO.parse('upgma.pp.fasta', 'fasta')}

ordered_outputs = []
for r in inputs:
    ordered_outputs.append(outputs[r.id])

SeqIO.write(ordered_outputs, f'upgma.pp.ordered.fasta', 'fasta')

