    # Fix scTE bugs: M<->MT mitochondrial naming mismatch

    # 1. Fix splitAllChrs in base.py (BAM processing)
    SCTE_BASE=$(find /opt/conda/envs/scTE -path "*/site-packages/scTE/base.py" | head -1) && \
        python3 -c "
p='$SCTE_BASE'
with open(p) as f: src=f.read()
old='''            # Force chrMT -> chrM
            if chrom == 'MT':
                chrom = 'M'
            else:
                continue'''
new='''            # Handle M<->MT mitochondrial naming mismatch (UCSC chrM vs Ensembl MT)
            if chrom == 'M' and 'MT' in chromosome_list:
                chrom = 'MT'
            elif chrom == 'MT' and 'M' in chromosome_list:
                chrom = 'M'
            else:
                continue'''
if old in src:
    src=src.replace(old,new)
    with open(p,'w') as f: f.write(src)
    print('patched base.py:',p)
else:
    print('base.py: already patched or different version')
"

    # 2. Fix scTE_build: chr_list and M<->MT normalization
    SCTE_BUILD=$(find /opt/conda/envs/scTE -path "*/bin/scTE_build" | head -1) && \
        python3 -c "
p='$SCTE_BUILD'
with open(p) as f: src=f.read()
changed = False

# Fix chr_list to include both M and MT
old1 = \"chr_list = [str(k) for k in list(range(1,50))] + ['X','Y','M']\"
new1 = \"chr_list = [str(k) for k in list(range(1,50))] + ['X','Y','M', 'MT']\"
if old1 in src:
    src = src.replace(old1, new1)
    print('fixed chr_list')
    changed = True

# Fix readGtf: normalize M<->MT
old2 = \"            if chrom.replace('chr','') not in chr_list:\\n                continue\"
new2 = \"            _chr = chrom.replace('chr','')\\n            if _chr not in chr_list:\\n                if _chr == 'M' and 'MT' in chr_list:\\n                    pass\\n                elif _chr == 'MT' and 'M' in chr_list:\\n                    pass\\n                else:\\n                    continue\"
if old2 in src:
    src = src.replace(old2, new2)
    print('fixed readGtf chr check')
    changed = True

# Fix TE processing: normalize M<->MT
old3 = \"            if chr not in chr_list:\\n                continue\"
new3 = \"            if chr not in chr_list:\\n                if chr == 'M' and 'MT' in chr_list:\\n                    pass\\n                elif chr == 'MT' and 'M' in chr_list:\\n                    pass\\n                else:\\n                    continue\"
if old3 in src:
    src = src.replace(old3, new3)
    print('fixed TE chr check')
    changed = True

if changed:
    with open(p,'w') as f: f.write(src)
    print('saved', p)
else:
    print('scTE_build: already patched or no matching text')
"
