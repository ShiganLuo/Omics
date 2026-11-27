#! /usr/bin/env python3


import subprocess, os, re
from Bio import AlignIO
from collections import defaultdict, Counter
from itertools import groupby
from operator import itemgetter

stkname = ''

def aligntRNA(tRNAseqs, out, threads=1):
	global stkname
	stkname = tRNAseqs.split(".fa")[0] + '_align.stk'
	cmfile = os.path.dirname(os.path.realpath(__file__)) + '/data/tRNAmatureseq.cm'

	if threads > 1:
		cpus = threads
	else:
		cpus = 0

	cmcommand = ['cmalign', '-o', stkname, '--nonbanded', '--cpu', str(cpus), '-g', cmfile, tRNAseqs]
	subprocess.check_call(cmcommand, stdout = open(out + 'cm.log', 'a'))

def extraCCA():
	extra_cca = list()
	stk = AlignIO.read(stkname, "stockholm")
	for record in stk:
		if record.seq[-3:] == 'cca':
			extra_cca.append(record.name)

	os.remove(stkname)

	return(extra_cca)

def tRNAclassifier(ungapped = False):

	struct_dict = structureParser()
	stk = AlignIO.read(stkname, "stockholm")

	ss_cons = stk.column_annotations['secondary_structure']
	if ungapped:
		ss_cons_orig = ss_cons
		ss_cons = ss_cons.replace(".", "")
	cons_pos = 0
	cons_pos_list = list()
	cons_pos_dict = defaultdict()
	openstem_count = 0
	closestem_count = 0
	enum = 0
	for pos, char in enumerate(ss_cons):
		if not ss_cons[pos] == ".":
			if cons_pos < 46:
				if (not cons_pos == 17) and (not cons_pos == 20):
					cons_pos_list.append(str(cons_pos))
					cons_pos_dict[pos+1] = str(cons_pos)
					cons_pos += 1
				elif (cons_pos == 17) and not ('17' in cons_pos_list):
					cons_pos_dict[pos+1] = '17'
					cons_pos_list.append('17')
				elif (cons_pos == 17) and ('17' in cons_pos_list):
					cons_pos_dict[pos+1] = '17a'
					cons_pos_list.append('17a')
					cons_pos += 1

				elif cons_pos == 20: 
					if not '20' in cons_pos_list:
						cons_pos_dict[pos+1] = '20'
						cons_pos_list.append('20')
					elif not '20a' in cons_pos_list:
						cons_pos_dict[pos+1] = '20a'
						cons_pos_list.append('20a')
					elif not '20b' in cons_pos_list:
						cons_pos_dict[pos+1] = '20b'
						cons_pos_list.append('20b')
						cons_pos += 1

			elif cons_pos == 46:
				if (not closestem_count == openstem_count) or (closestem_count == 0 or openstem_count == 0):
					if ss_cons[pos] == "<":
						openstem_count += 1
						enum += 1
						current_e = "e" + str(enum)
						cons_pos_dict[pos+1] = current_e
						cons_pos_list.append(current_e)
					elif ss_cons[pos] == ">":
						closestem_count += 1
						enum += 1
						current_e = "e" + str(enum)
						cons_pos_dict[pos+1] = current_e
						cons_pos_list.append(current_e)
					elif ss_cons[pos] == "_":
						enum += 1
						current_e = "e" + str(enum)
						cons_pos_dict[pos+1] = current_e
						cons_pos_list.append(current_e)

				elif (closestem_count == openstem_count) and (not closestem_count == 0 or not openstem_count == 0):
					cons_pos_dict[pos+1] = str(cons_pos)
					cons_pos_list.append(str(cons_pos))
					cons_pos += 1

			elif cons_pos > 46:
				cons_pos_dict[pos+1] = str(cons_pos)
				cons_pos_list.append(str(cons_pos))
				cons_pos += 1

		elif ss_cons[pos] == ".":
			cons_pos_dict[pos+1] = '-'

	cons_pos_list = "_".join(cons_pos_list)

	tRNA_struct = defaultdict(dict)
	tRNA_ungap2canon = defaultdict(dict)

	for record in stk:
		tRNA = record.id
		seq = record.seq
		ungapped_pos = 0
		bases = ["A", "C", "G", "U"]

		if not ungapped:
			for i, letter in enumerate(seq, 1):
				if letter.upper() in bases:
					tRNA_ungap2canon[tRNA][ungapped_pos] = cons_pos_dict[i]
					ungapped_pos += 1
					tRNA_struct[tRNA][i] = struct_dict[i]
				else:
					tRNA_struct[tRNA][i] = 'gap'
		elif ungapped:
			n = 0
			for i, letter in enumerate(seq,1):
				if letter.upper() in bases:
					n += 1
					try:
						tRNA_ungap2canon[tRNA][ungapped_pos] = cons_pos_dict[n]
					except KeyError:
						break
					tRNA_struct[tRNA][n] = struct_dict[i]
					ungapped_pos += 1
				elif letter == "-" and ss_cons_orig[i-1] != ".":
					n += 1
					tRNA_struct[tRNA][n] = 'gap'

	return(tRNA_struct, tRNA_ungap2canon, cons_pos_list, cons_pos_dict)

def tRNAclassifier_nogaps(oneBased = False):

	struct_dict = structureParser()
	tRNA_struct = defaultdict(dict)

	stk = AlignIO.read(stkname, "stockholm")
	for record in stk:
		tRNA = record.id
		seq = record.seq
		pos = 0
		bases = ["A", "C", "G", "U"]

		for i, letter in enumerate(seq):
			if letter.upper() in bases:
				if oneBased:
					tRNA_struct[tRNA][pos+1] = struct_dict[i+1]
				else:
					tRNA_struct[tRNA][pos] = struct_dict[i+1]
				pos += 1

	return(tRNA_struct)

def getAnticodon():

	anticodon = list()
	stk = AlignIO.read(stkname, "stockholm")
	rf_cons = stk.column_annotations['reference_annotation']
	for pos, char in enumerate(rf_cons):
		if char == "*":
			anticodon.append(pos)
	
	return(anticodon)

def clusterAnticodon(cons_anticodon, cluster):

	bases = ["A", "C", "G", "U"]
	stk = AlignIO.read(stkname, "stockholm")
	cluster_anticodon = list()
	for record in stk:
		if record.id == cluster:
			for pos in cons_anticodon:
				gapcount = 0
				for char in record.seq[:pos]:
					if char.upper() not in bases:
						gapcount += 1
				cluster_anticodon.append(pos - gapcount)

	return(cluster_anticodon)

def modContext(out, unsplitCluster_lookup):

	cons_pos_list, cons_pos_dict = tRNAclassifier()[2:4]

	sites_dict = defaultdict()
	mod_sites = ['9', '20', '26', '32', '34', '37', '58']
	
	for mod in mod_sites:
		sites_dict[mod] = list(cons_pos_dict.keys())[list(cons_pos_dict.values()).index(mod)]

	upstream_dict = defaultdict(lambda: defaultdict(list))

	stk = AlignIO.read(stkname, "stockholm") 
	for record in stk:
		gene = record.id
		seq = record.seq
		for pos in sites_dict.values():
			identity = seq[pos-1]
			if identity in ['A','C','G','U','T']:
				up = pos - 2
				down = pos
				while seq[up].upper() not in ['A','C','G','U','T']:
					up -= 1
				while seq[down].upper() not in ['A','C','G','U','T']:
					down += 1
				canon_pos = cons_pos_dict[pos]
				upstream_dict[gene][canon_pos].append(identity) 
				upstream_dict[gene][canon_pos].append(seq[up])
				upstream_dict[gene][canon_pos].append(seq[down])

	try:
		os.mkdir(out + "mods")
	except FileExistsError:
		pass

	with open(out + "mods/modContext.txt", 'w') as outfile:
		outfile.write("cluster\tcanon_pos\tidentity\tupstream\tdownstream\n")
		for cluster, data in upstream_dict.items():
			if cluster in unsplitCluster_lookup.keys():
				cluster = unsplitCluster_lookup[cluster]
			for pos, base in data.items():
				outfile.write(cluster + "\t" + str(pos) + "\t" + base[0] + "\t" + base[1] + "\t" + base[2] + "\n")

	mod_sites = str("_".join(str(e) for e in mod_sites))

	return(mod_sites, cons_pos_list)

def structureParser():
	
	struct_dict = dict()
	ss_cons = "".join([line.split()[-1] for line in open(stkname) if line.startswith("#=GC SS_cons")])

	acc = defaultdict()

	term = defaultdict()
	term_type = "5'"

	bulge_list = []
	bulge_items = []
	bulge_count = 0

	stemloops = defaultdict()
	stemloops_count = 0
	stemloops_type = ['D stem-loop','Anticodon stem-loop','Variable loop','T stem-loop']
	open_count = 0
	close_count = 0

	for pos, char in enumerate(ss_cons):
		if char == ':':
			if pos < 10:
				term[pos+1] = term_type
			else:
				term_type = "3'"
				term[pos+1] = term_type

		if char == '(':
			acc[pos+1] = "Acceptor stem 5'"

		if char == ')':
			acc[pos+1] = "Acceptor stem 3'"

		if char == "<":
			open_count += 1
			stemloops[pos+1] = stemloops_type[stemloops_count]
		if char == ">":
			close_count +=1
			stemloops[pos+1] = stemloops_type[stemloops_count]
			if close_count == open_count:
				stemloops_count += 1
				open_count = 0
				close_count = 0

	for i in ["Acceptor stem 5'","Acceptor stem 3'"]:
		pos_list = [k for k, v in acc.items() if v == i]
		start = min(pos_list)
		stop = max(pos_list)
		acc.update([(n, i) for n in range(start, stop +1)])

	for i in stemloops_type:
		pos_list = [k for k, v in stemloops.items() if v == i]
		start = min(pos_list)
		stop = max(pos_list)
		stemloops.update([(n, i) for n in range(start,stop+1)])

	struct_dict = {**term, **acc}
	struct_dict.update(stemloops)

	for pos, char in enumerate(ss_cons):
		if (pos+1) not in [x for x in struct_dict.keys()]:
			bulge_list.append(pos+1)

	for k, g in groupby(enumerate(bulge_list), lambda x:x[0] - x[1]):
		group = map(itemgetter(1), g)
		group = list(map(int,group))
		bulge_items.append(group)

	for bulge in bulge_items:
		bulge_count += 1
		if bulge_count == 3 or bulge_count == 4:
			for pos in bulge:
				struct_dict[pos] = 'Variable loop'
		elif bulge_count == 5:
			for pos in bulge:
				struct_dict[pos] = 'bulge3'
		else:
			for pos in bulge:
				struct_dict[pos] = 'bulge' + str(bulge_count)

	return(struct_dict)
