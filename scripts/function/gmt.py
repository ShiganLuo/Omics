def read_gmt(gmt:str) -> dict:
    gene_sets = {}
    try:
        with open(gmt,'r',encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >=3:
                    set_name = parts[0]
                    genes = parts[2:]
                    gene_sets[set_name] = [g for g in genes if g]
        return gene_sets
    except FileNotFoundError:
        print(f"file not found: {gmt}")
        return None
    except Exception as e:
        print(f"some unexpected wrong with reading file: {e} ")
    return None
    
    
if __name__ == "__main__":
    read_gmt()