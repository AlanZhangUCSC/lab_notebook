import sys
from read_fasta import read_record

def compare_sequences(seq1, seq2):
    """
    Analyze two sequences for matches, SNPs, and gaps.
    
    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        
    Returns:
        tuple: (match_string, num_errors, num_snps, num_gaps, num_gaps_from_ends)
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")
        
    match_string = ''
    canonical = {
        'a': True, 't': True, 'c': True, 'g': True,
        'A': True, 'T': True, 'C': True, 'G': True
    }
    
    num_errors = 0
    num_snps = 0
    num_gaps = 0
    num_gaps_from_ends = 0
    snp_positions = []
    gap_positions = []
    left_nongap_start = -1
    right_nongap_end = -1
    
    # Process sequence differences
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            match_string += ':'
            if seq1[i] == '-' or seq2[i] == '-':
                num_errors += 1
                num_gaps += 1
                gap_positions.append(i)
            elif seq1[i] in canonical and seq2[i] in canonical:
                num_snps += 1
                num_errors += 1
                snp_positions.append(i)
        else:
            match_string += '|'
            
        if i % 1000 == 0:
            print(f'\r{i} positions processed', end='', file=sys.stderr)
    
    # Count gaps from ends
    for i in range(len(seq1)):
        if seq1[i] == '-' or seq2[i] == '-':
            num_gaps_from_ends += 1
        else:
            left_nongap_start = i
            break
    
    for i in range(len(seq1)-1, -1, -1):
        if seq1[i] == '-' or seq2[i] == '-':
            num_gaps_from_ends += 1
        else:
            right_nongap_end = i
            break
    
    return match_string, num_errors, num_snps, num_gaps, num_gaps_from_ends, snp_positions, gap_positions, left_nongap_start, right_nongap_end

def main(mafft_path):
    import sys
    from read_fasta import read_record
    
    seqs = []
    for idn, seq_str in read_record(mafft_path):
        seqs.append((idn, seq_str))
    
    idn1, seq1 = seqs[0]
    idn2, seq2 = seqs[1]
    
    print(f"Sequence lengths are {len(seq1)} and {len(seq2)}", file=sys.stderr)
    
    match_string, num_errors, num_snps, num_gaps, num_gaps_from_ends, snp_positions, gap_positions, left_nongap_start, right_nongap_end = compare_sequences(seq1, seq2)
    
    # Print results
    line_size = 80
    for i in range(0, len(match_string), line_size):
        print(seq1[i:i+line_size])
        print(match_string[i:i+line_size])
        print(seq2[i:i+line_size])
        print()
    

    print(f'{num_errors} errors out of {len(seq1)} positions')
    print(f'{num_snps} snps', '.' if not snp_positions else ','.join(map(str, snp_positions)))
    print(f'{num_gaps} gaps', left_nongap_start, right_nongap_end, '.' if not gap_positions else ','.join(map(str, gap_positions)))
    print(f'{num_gaps - num_gaps_from_ends} gaps edge corrected')

if __name__ == "__main__":
    main(sys.argv[1])