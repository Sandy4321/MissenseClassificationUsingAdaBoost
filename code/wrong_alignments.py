
import os

proteins = os.listdir('../alignments')
gg = [k for k in proteins if k!='realigned']

wrongs = []

for prot in gg:
    with open('../alignments/%s'%prot) as fp:
        line = fp.readline()
        if len(line.strip().split(' ')) == 2:
	    wrongs.append(prot)

print(len(wrongs))
