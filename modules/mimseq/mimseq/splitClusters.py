#! /usr/bin/env python3

from __future__ import absolute_import
import logging
from collections import defaultdict
from itertools import chain, combinations as comb
from numpy.core.fromnumeric import shape
from .ssAlign import aligntRNA, tRNAclassifier
from .getCoverage import getBamList
import re, copy
import pybedtools
from multiprocessing import Pool
from functools import partial
import subprocess


log = logging.getLogger(__name__)

tRNA_ungap2canon = ''

def dd_set():
	return(defaultdict(set))

def dd():
	return(defaultdict(list))

def GetPowerset(s):
    
    return list(chain(*map(lambda x: comb(s, x), range(0, len(s)+1))))

def reformatInDelDict (dict, newDict, type):

    for cluster, data in dict.items():
        for pos, members in data.items():
            for member in members:
                isodecoder = "-".join(member.split("-")[:4]) if not "chr" in member else member
                posIdentity = str(pos) + type
                newDict[cluster][isodecoder].add(posIdentity)

    return(newDict)

def updateMismatchPosDict (mismatchPosDict, inputDict):
    for cluster, data in inputDict.items():
        for isodecoder, posInfo in data.items():
            mismatchPosDict[cluster][isodecoder].update(posInfo)

    return(mismatchPosDict)

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def natural_keys_list(text):
    l = list()
    for a in text:
        for c in re.split(r'(\d+)', a):
            l.append(atoi(c))
    
    return(l)

def findUniqueSubset (inputDict, outputDict, tRNA_dict):

    tRNA_ungap2canon = tRNAclassifier()[1]
    mods = {'9':(['G', 'A'],[]), 
            '20':(['C', 'T'],[]), 
            '20a':(['T'],[]), 
            '20b':(['T'],[]), 
            '26':(['G'],[]), 
            '32':(['C'],[]), 
            '37':(['G', 'A'],[]), 
            '58':(['A'],[]),
            'e9':(['C'],['Leu-CAG','Ser-AGA','Ser-CGA','Ser-GCT','Ser-TGA'])}

    notSplit_mods_clusterInfo = defaultdict(list) 
    notSplit_mods_posInfo = defaultdict(set)

    for cluster, data in inputDict.items():
        powersets  = [GetPowerset(data[s]) for s in data]
        temp_uss_dict = defaultdict(list)
        temp_fullset_dict = defaultdict(list)
        
        modsFilter = {i:v[0] for i, v in mods.items() if (i in tRNA_ungap2canon[cluster].values()) and (any([re.search(x, cluster) for x in v[1]] if v[1] else ['empty']))}
        clusterMods_pos = {list(tRNA_ungap2canon[cluster].keys())[list(tRNA_ungap2canon[cluster].values()).index(i)]:v for i, v in modsFilter.items()}
        clusterMods_list = [str(i) + str(x) for i, nucs in clusterMods_pos.items() for x in nucs]

        for i in range(0, len(data)):
            isodecoder, mismatch_powerset = list(data)[i], powersets[i]
            temp_fullset_dict[isodecoder] = tuple(sorted((max(mismatch_powerset, key=len))))
            ops = [x for ind, x in enumerate(powersets) if ind != i]
            uss = list(set(mismatch_powerset)-set([x for y in ops for x in y if x in mismatch_powerset]))
            if not uss:
                uss = [(max(mismatch_powerset, key=len))]
            uss = list(filter(None, uss))
            uss = sorted(uss, key = natural_keys_list, reverse = True)

            temp_modPos = list()
            for i, s in enumerate(uss):
                uss[i] = tuple(sorted(s, key=natural_keys))
                ussPosOnly = [int(re.search("([0-9]*)[A-Z]?", x).group(1)) for x in uss[i] if not re.search("Ins|Del", x)]
                ussClusterParent = [str(x) + tRNA_dict[cluster]['sequence'][x] for x in ussPosOnly]
                if any(x in clusterMods_list for x in uss[i]) or any(y in clusterMods_list for y in ussClusterParent):
                    temp_modPos.append(s)
                    del uss[i]

            if temp_modPos:
                temp_modPos = [set(x) for x in temp_modPos]
                modPos = set.intersection(*temp_modPos)
                modPos =  {int(re.search("([0-9]*)[A-Z]?", x).group(1)) for x in modPos}
                canonmodPos = {tRNA_ungap2canon[cluster][x] for x in modPos}
                canonmodPos = canonmodPos.intersection(set(mods.keys()))

            if not uss:
                notSplit_mods_clusterInfo[cluster].append(isodecoder)
                notSplit_mods_posInfo[cluster] = canonmodPos
                continue

            temp_uss_dict[isodecoder] = uss
            min_uss = tuple(sorted(min(uss, key = len)))
            index = 0
            while min_uss in outputDict[cluster].keys():
                index += 1
                try:
                    min_uss = uss[index]
                except IndexError:
                    conflictingIso = outputDict[cluster][min_uss][0]
                    conflictingIso_uss = None
                    iso_ss = temp_uss_dict[conflictingIso]
                    for ss in iso_ss:
                        if not ss in outputDict[cluster].keys():
                            conflictingIso_uss = ss
                            outputDict[cluster][conflictingIso_uss].append(conflictingIso)
                            outputDict[cluster][conflictingIso_uss].append(temp_fullset_dict[conflictingIso])
                            del outputDict[cluster][min_uss]
                            break
                    if not conflictingIso_uss:
                        conflictingIso_uss = temp_fullset_dict[isodecoder]
                        outputDict[cluster][conflictingIso_uss].append(conflictingIso)
                        outputDict[cluster][conflictingIso_uss].append(conflictingIso_uss)
                        del outputDict[cluster][min_uss]
        
            outputDict[cluster][min_uss].append(isodecoder)
            outputDict[cluster][min_uss].append(temp_fullset_dict[isodecoder])

    for cluster, data in outputDict.items():
        for uss, isodecoder in data.items():
            isodecoder_list = [x for x in tRNA_dict.keys() if isodecoder[0] in x and not "chr" in isodecoder[0]]
            iso_min = min([x.split("-")[-1] for x in isodecoder_list])
            newIso = isodecoder[0] + "-" + str(iso_min) if not "chr" in isodecoder[0] else isodecoder[0]
            newData = [newIso, isodecoder[1]]
            outputDict[cluster][uss] = newData
    for cluster, data in notSplit_mods_clusterInfo.items():
        for i, isodecoder in enumerate(data):
            isodecoder_list = [x for x in tRNA_dict.keys() if isodecoder in x and not "chr" in isodecoder]
            iso_min = min([x.split("-")[-1] for x in isodecoder_list])
            newIso = isodecoder + "-" + str(iso_min) if not "chr" in isodecoder else isodecoder
            notSplit_mods_clusterInfo[cluster][i] = newIso

    return(outputDict, notSplit_mods_clusterInfo, notSplit_mods_posInfo)

def splitIsodecoder(cluster_perPos_mismatchMembers, insert_dict, del_dict, tRNA_dict, cluster_dict, out_dir, experiment_name):
    log.info("\n+------------------------------------------------------------------------------+\
        \n| Characterizing cluster mismatches for read splitting by unique tRNA sequence |\
       \n+------------------------------------------------------------------------------+")

    log.info("** Assessing mismatches between cluster members and parent... **")

    cluster_MemberMismatchPos = defaultdict(dd_set)

    for cluster, data in cluster_perPos_mismatchMembers.items():
        for pos, members in data.items():
            for member in members:
                member_seq = tRNA_dict[member]['sequence']
                ins_num = len(set([ins for ins in insert_dict[cluster] if member in insert_dict[cluster][ins] and ins < pos]))
                del_num = len(set([deletion for deletion in del_dict[cluster] if member in del_dict[cluster][deletion] and deletion <= pos]))
                identity = member_seq[pos-ins_num+del_num]

                isodecoder = "-".join(member.split("-")[:-1]) if not "chr" in member else member
                posIdentity = str(pos) + identity
                cluster_MemberMismatchPos[cluster][isodecoder].add(posIdentity)
                
    cluster_MemberInsertPos = defaultdict(dd_set)
    cluster_MemberInsertPos = reformatInDelDict(insert_dict, cluster_MemberInsertPos, "Ins")
    cluster_MemberMismatchPos =updateMismatchPosDict(cluster_MemberMismatchPos, cluster_MemberInsertPos)

    cluster_MemberDeletePos = defaultdict(dd_set)
    cluster_MemberDeletePos = reformatInDelDict(del_dict, cluster_MemberDeletePos, "Del")
    cluster_MemberMismatchPos = updateMismatchPosDict(cluster_MemberMismatchPos, cluster_MemberDeletePos)
    
    unique_isodecoderMMs = defaultdict(dd)
    unique_isodecoderMMs, notSplit_mods_clusterInfo, notSplit_mods_posInfo = findUniqueSubset(cluster_MemberMismatchPos, unique_isodecoderMMs, tRNA_dict)
    
    singleSeq_clusters_num = 0
    multiSeq_size = 0
    multiSeq_names = set()
    for cluster, members in cluster_dict.items():
        child_iso_set = {tRNA_dict[member]['sequence'].upper() for member in members}
        if len(child_iso_set) == 1:
            singleSeq_clusters_num += 1
        else:
            multiSeq_size += len(child_iso_set)
            multiSeq_names.update({"-".join(member.split("-")[:-1]) for member in members})

    deconv_sequences_num = 0
    deconv_names = set()
    for cluster, data in unique_isodecoderMMs.items():
        unique_deconv = {member[0] for member in data.values()}
        deconv_sequences_num += len(unique_deconv)
        if len(unique_deconv) == len(cluster_MemberMismatchPos[cluster]):
            deconv_sequences_num += 1
            deconv_names.update({"-".join(cluster.split("-")[:-1])})
        deconv_names.update({"-".join(member[0].split("-")[:-1]) for member in data.values()})

    notSplit_mods_clusterInfo = defaultdict(set, ((k, set(v)) for k, v in notSplit_mods_clusterInfo.items()))
    
    log.info("Total unique sequences: {}".format(singleSeq_clusters_num + multiSeq_size))
    log.info("Single sequence clusters: {}".format(singleSeq_clusters_num))
    log.info("Deconvoluted sequences: {}".format(deconv_sequences_num))
    log.info("Single transcript resolution for {} sequences".format(singleSeq_clusters_num + deconv_sequences_num))
    log.info("Deconvolution not possible for {} sequences".format((multiSeq_size - deconv_sequences_num)))

 
    nondeconv_isos = list(multiSeq_names - deconv_names)
    if ((singleSeq_clusters_num + multiSeq_size) != (singleSeq_clusters_num + deconv_sequences_num)) or (notSplit_mods_clusterInfo):
        log.warning("\n*** Deconvolution for some unique transcripts not possible ***\n\
*** due to overlap of distinguishing mismatches with canonical, potentially modified sites ***\n\
This might also be due to errors in input sequence naming in reference.\n\
Ensure all transcripts with identical isodecoder numbers are indeed unique in mature sequence (see *_tRNATranscripts.fa).")
        log.warning("Continuing with analysis - clusters with unsplit transcripts will be renamed accordingly in output.")

    splitBool = notSplit_mods_clusterInfo

    return(unique_isodecoderMMs, splitBool, notSplit_mods_posInfo)

def covCheck_mp(bedTool, unique_isodecoderMMs, tRNA_ungap2canon, splitBool, covDiff, input):
    unsplit = set()
    unsplit_isosOnly = set()
    notSplit_cov_posInfo = defaultdict(list)

    log.info("Calculating nucleotide coverage for {}".format(input))
    bam = pybedtools.BedTool(input)
    temp_chrom = input + "_chrom.txt"
    cmd = "samtools view -H " + input + " | grep @SQ | sed 's/@SQ\tSN:\|@sq\tSN:\|LN://g' > " + temp_chrom
    subprocess.call(cmd, shell = True)
    temp_chrom_l = list()
    with open(temp_chrom, "r") as temp_chrom_f:
        for line in temp_chrom_f:
            line = line.strip()
            temp_chrom_l.append(line.split("\t")[0])
    bedTool = bedTool.filter(lambda b: b.chrom in temp_chrom_l)
    bedTool = bedTool.sort(faidx=temp_chrom)
    cov = bedTool.coverage(bam, s = True, d = True, sorted = True, g = temp_chrom)
    cmd = "rm " + temp_chrom
    subprocess.call(cmd, shell = True)
    cov_df = cov.to_dataframe()

    for cluster in unique_isodecoderMMs.keys():
        if cov_df.name.str.contains(cluster).any():
            for mismatch in unique_isodecoderMMs[cluster]:
                mismatchNums = [int(re.search('(\d+\.|\d+)+', pos).group(0)) for pos in mismatch]
                minMismatch = min(mismatchNums) + 1

                end = max(cov_df.loc[cov_df.name == cluster,'thickStart'])
                if minMismatch < end:
                    endCov = int(cov_df.loc[(cov_df.name == cluster) & (cov_df.thickStart == end - 5), 'thickEnd'])
                    ratio = cov_df.loc[(cov_df.name == cluster) & (cov_df.thickStart == minMismatch), 'thickEnd'].astype(float) / endCov

                    if float(ratio) < covDiff:
                        unsplit.add(cluster)
                        unsplit.add(unique_isodecoderMMs[cluster][mismatch][0])
                        unsplit_isosOnly.add(unique_isodecoderMMs[cluster][mismatch][0])
                        splitBool[cluster].add(unique_isodecoderMMs[cluster][mismatch][0])
                        notSplit_cov_posInfo[cluster].append(minMismatch)

    notSplit_cov_posInfo = {k:[tRNA_ungap2canon[k][x]] for k, v in notSplit_cov_posInfo.items() for x in v}

    return(unsplit, unsplit_isosOnly, splitBool, notSplit_cov_posInfo)

def unsplitClustersCov(coverageData, coverageBed, unique_isodecoderMMs, splitBool, threads, map_round, covDiff = 0.5):
   
    log.info("\n+---------------------------------------------------------------------------------+\
        \n| Determining non-deconvoluted clusters due to insufficient coverage at mismatches |\
        \n+---------------------------------------------------------------------------------+")

    log.info("*** Un-deconvoluted sequences being defined based on coverage difference more than {:.2%}".format(covDiff))
    baminfo, bamlist = getBamList(coverageData)
    
    if len(baminfo) > threads:
        multi = threads
    else:
        multi = len(baminfo)

    log.info("Determining unsplittable sequences...")
    pool = Pool(multi)
    bed = pybedtools.BedTool(coverageBed)

    if map_round == 1:
        global tRNA_ungap2canon
        tRNA_ungap2canon = tRNAclassifier()[1]
    func = partial(covCheck_mp, bed, unique_isodecoderMMs, tRNA_ungap2canon, splitBool, covDiff)
    unsplit, unsplit_isosOnly, splitBool, notSplit_cov_posInfo = zip(*pool.map(func, bamlist))
    pool.close()
    pool.join()

    newSplitBool = {cluster: set().union(data) for proc in splitBool for cluster, data in proc.items()}
    newSplitBool = defaultdict(set, newSplitBool)

    unique_isodecoderMMs_new = copy.deepcopy(unique_isodecoderMMs)
    for cluster, isos in newSplitBool.items():
        if cluster in unique_isodecoderMMs.keys():
            for mismatch, data in unique_isodecoderMMs[cluster].items():
                if data[0] in isos:
                    del unique_isodecoderMMs_new[cluster][mismatch]
    
    unsplit_all = list(set().union(*unsplit))
    unsplit_isosAll = len(list(set().union(*unsplit_isosOnly)))
    log.info("{} unique sequences not split from cluster parent".format(unsplit_isosAll))
    log.info("{} parents and isododecoders not deconvoluted due to insufficient coverage at mismatches".format(len(unsplit_all)))

    notSplit_cov_posInfo_all = {k:set(v).union(set(v)) for tup in notSplit_cov_posInfo for k, v in tup.items()}

    return(newSplitBool, unique_isodecoderMMs_new, notSplit_cov_posInfo_all)

def getDeconvSizes(splitBool, tRNA_dict, cluster_dict, unique_isodecoderMMs):

    isodecoder_sizes = defaultdict(int)

    unsplitIsos = [x for k in splitBool.values() for x in k]
    unsplitClusters = [x for x in splitBool.keys()]

    for cluster, members in cluster_dict.items():
        child_iso_list = [tRNA_dict[member]['sequence'].upper() for member in members]
        if len(set(child_iso_list)) == 1:
            isodecoder_sizes[cluster] = len(child_iso_list)

    unsplitCluster_lookup = defaultdict()
    for cluster, isos in splitBool.items():
        if "tRNA" in cluster:
            member_IsoNums = tuple(int(iso.split("-")[-2]) if "tRNA" in iso else "tRX" + iso.split("-")[-2] for iso in isos)
        elif "tRX" in cluster:
            member_IsoNums = tuple(int(iso.split("-")[-2]) if "tRX" in iso else "tRNA" + iso.split("-")[-2] for iso in isos)
        member_IsoNums = sorted(tuple(x for x in member_IsoNums if isinstance(x, int))) + [x for x in member_IsoNums if isinstance(x, str)]
        member_IsoString = "/" + "/".join([str(iso) for iso in member_IsoNums])
        newClusterName = "-".join(cluster.split("-")[:-1]) + member_IsoString
        shortClusterName = "-".join(cluster.split("-")[:-1])
        unsplitCluster_lookup[shortClusterName] = newClusterName

    for cluster, isos in splitBool.items():
        shortClusterName = "-".join(cluster.split("-")[:-1])
        clusterSize = len([info['sequence'].upper() for info in tRNA_dict.values() if info['sequence'].upper() == tRNA_dict[cluster]['sequence'].upper()])
        totalIsoSize = 0
        for isodecoder in isos:
            totalIsoSize += len([info['sequence'].upper() for info in tRNA_dict.values() if info['sequence'].upper() == tRNA_dict[isodecoder]['sequence'].upper()])
        isodecoder_sizes[unsplitCluster_lookup[shortClusterName]] = clusterSize + totalIsoSize

    for cluster, data in unique_isodecoderMMs.items():
        if not cluster in unsplitClusters:
            isodecoder_sizes[cluster] = len([info['sequence'].upper() for info in tRNA_dict.values() if info['sequence'].upper() == tRNA_dict[cluster]['sequence'].upper()])
        for isodecoder in data.values():
            isodecoder_sizes[isodecoder[0]] = len([info['sequence'].upper() for info in tRNA_dict.values() if info['sequence'].upper() == tRNA_dict[isodecoder[0]]['sequence'].upper()])

    return(isodecoder_sizes, unsplitCluster_lookup)

def writeDeconvTranscripts(out_dir, experiment_name, tRNA_dict, isodecoder_sizes, clustering):

    with open(out_dir + experiment_name + '_isodecoderTranscripts.fa', 'w') as tempSeqs:
        for seq in isodecoder_sizes.keys():
            if "/" in seq:
                tempSeq = seq.split("/")[0] + "-"
                isodecoder_list = [x for x in tRNA_dict.keys() if tempSeq in x and not "chr" in x]
                iso_min = min([x.split("-")[-1] for x in isodecoder_list])
                newSeq = tempSeq + str(iso_min)
            else:
                newSeq = seq
            if clustering:
                shortname = "-".join(seq.split("-")[:-1]) if not re.search("chr|/", seq) else seq.split("/")[0]
            else:
                shortname = seq
            tempSeqs.write(">" + shortname + "\n" + tRNA_dict[newSeq]['sequence'] + "\n")
    aligntRNA(tempSeqs.name, out_dir)

def writeIsodecoderInfo(out_dir, experiment_name, isodecoder_sizes, readRef_unsplit_newNames, tRNA_dict):

    for name in readRef_unsplit_newNames:
        base = "-".join(name.split("-")[:-1])
        count = 0
        for num in name.split("-")[-1].split("/"):
            count += 1
            if "tRNA" in base:
                if not "tRX" in num:
                    iso = base + "-" + num
                else:
                    num = num.replace("tRX", "")
                    iso = base.replace("tRNA", "tRX") + "-" + num
            elif "tRX" in base:
                if not "tRNA" in num:
                    iso = base + "-" + num
                else:
                    num = num.replace("tRX", "")
                    iso = base.replace("tRNA", "tRX") + "-" + num
            isodecoder_list = [x for x in tRNA_dict.keys() if iso == "-".join(x.split("-")[:-1]) and not "chr" in x]
            iso_min = min([x.split("-")[-1] for x in isodecoder_list])
            gene = iso + "-" + str(iso_min)
            try:
                del isodecoder_sizes[gene]
            except KeyError:
                while not gene in isodecoder_sizes:
                    iso_min = int(iso_min) + 1
                    gene = iso + "-" + str(iso_min)
                del isodecoder_sizes[gene]
        isodecoder_sizes[name] = count

    with open(out_dir + experiment_name + "isodecoderInfo.txt", "w") as isodecoderInfo:
        isodecoderInfo.write("Isodecoder\tsize\n")
        for isodecoder, size in isodecoder_sizes.items():
            isodecoder = "-".join(isodecoder.split("-")[:-1]) if not re.search("chr|/", isodecoder) else isodecoder
            isodecoderInfo.write(isodecoder + "\t" + str(size) + "\n")

    return(isodecoder_sizes)

def writeSplitInfo(out, name, splitBool, notSplit_mods_posInfo, notSplit_cov_posInfo):

    with open(out + name + "_unsplitClusterInfo.txt", "w") as unsplitOut:
        unsplitOut.write("Parent\tSize (total unsplit transcripts)\tUnsplit transcripts\tReason\n")
        for cluster, members in splitBool.items():
            size = len(members) + 1
            members = ["-".join(x.split("-")[:-1]) for x  in members]
            member_print = ", ".join(members)
            if cluster in notSplit_mods_posInfo:
                reason = "potential mod at mismatch {}".format(", ".join(notSplit_mods_posInfo[cluster]))
            elif cluster in notSplit_cov_posInfo:
                reason = "insufficient coverage at mismatch {}".format(", ".join(notSplit_cov_posInfo[cluster]))
            else:
                reason = ">10% parent assigned reads do not match parent sequence"
            unsplitOut.write(cluster + "\t" + str(size) + "\t" + member_print + "\t" + reason + "\n")
    log.info("Info on all clusters not deconvoluted written to {}".format(out + "annotation/" + name + "_unsplitClusterInfo.txt"))

def writeIsodecoderTranscripts(out_dir, experiment_name, cluster_dict, tRNA_dict):

	with open(out_dir + experiment_name + '_isodecoderTranscripts.fa', 'w') as tempSeqs:
		for seq in cluster_dict.keys():
			shortname = "-".join(seq.split("-")[:-1]) if not "chr" in seq else seq
			tempSeqs.write(">" + shortname + "\n" + tRNA_dict[seq]['sequence'] + "\n")
	aligntRNA(tempSeqs.name, out_dir)

def getIsodecoderSizes(out_dir, experiment_name, tRNAdict):

	isodecoder_sizes = defaultdict(int)
	already_added = set()
	for tRNA in tRNAdict:
		if tRNA not in already_added:
			sameSeq = [tRNAs for tRNAs, data in tRNAdict.items() if data['sequence'] == tRNAdict[tRNA]['sequence']]
			already_added.update(sameSeq)
			isodecoder_sizes[tRNA] = len(sameSeq)

	with open(out_dir + experiment_name + "isodecoderInfo.txt", "w") as isodecoderInfo:
		isodecoderInfo.write("Isodecoder\tsize\n")
		for isodecoder, size in isodecoder_sizes.items():
			isodecoderInfo.write(isodecoder + "\t" + str(size) + "\n")

	return(isodecoder_sizes)