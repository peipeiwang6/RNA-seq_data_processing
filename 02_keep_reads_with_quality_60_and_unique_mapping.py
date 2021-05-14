import sys,os
sam = sys.argv[1]
out = open(sam.split('.sam')[0] + '_quality_60_unique.sam','w')
inp = open(sam,'r')
inl = inp.readline()
read = []
while inl:
	if inl.startswith('@'):
		out.write(inl)
	else:
		tem = inl.split('\t')
		if len(read)==0:
			read.append(tem[0])
			if inl.split('\t')[4] == '60':
				read.append(inl)
		else:
			if tem[0] == read[0]:
				if tem[4] == '60':
					read.append(inl)
			else:
				if len(read)==3:
					out.write(read[1])
					out.write(read[2])
					out.flush()
				read = []
				read.append(tem[0])
				if inl.split('\t')[4] == '60':
					read.append(inl)
	inl = inp.readline()

out.close()



