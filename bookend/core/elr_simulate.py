import sys
import os
import bookend.core.cython_utils._rnaseq_utils as ru
if __name__ == '__main__':
    sys.path.append('../../bookend')
    from argument_parsers import simulate_parser as parser

class ELRsimulator:
    def __init__(self, args):
        """Simulates a set of ELR reads based on a reference transcriptome and abundances."""
        self.output = args['OUT']
        self.reference = args['REFERENCE']
        self.abundance = args['ABUNDANCE']
        self.count = args['COUNT']
        self.read_length = args['READ_LENGTH']
        self.min_length = args['MIN_LENGTH']
        self.adapter_5p = args['ADAPTER_5P']
        self.adapter_3p = args['ADAPTER_3P']
        self.paired = args['PAIRED']
        self.fragment_mean = args['FRAGMENT_MEAN']
        self.fragment_sd = args['FRAGMENT_SD']
        self.var_5p = args['VAR_5P']
        self.var_3p = args['VAR_3P']
        self.percent_intact = args['PERCENT_INTACT']
        self.percent_sense = args['PERCENT_SENSE']
        self.label_noise = args['LABEL_NOISE']
        self.seed = args['SEED']
        self.bed = args['BED']
        
        self.percent_intact = float(1) if self.percent_intact < 1 else (float(100) if self.percent_intact > 100 else self.percent_intact)
        self.percent_sense = float(0) if self.percent_sense < 0 else (float(100) if self.percent_sense > 100 else self.percent_sense)
        random.seed(args.SEED)
        
    def run(self):
        if self.output != 'stdout':
            print(self.display_options())
        
        self.import_annotation()
        self.generate_reads()
        if self.output != 'stdout':
            print(self.display_summary())
    
    def import_annotation():
        self.dataset = ru.AnnotationDataset(
            annotation_files=[],
            reference=self.reference
        )
    
    def generate_reads():
        anno_dict = self.dataset.reference
        prob_dict = self.parse_abundance_file(args.TPM)
        chroms_dict = {}
        reverse_chroms_dict = {}
        chrom_counter = 0
        names_dict = {}
        reverse_names_dict = {}
        name_counter = 0    
        randomized_transcripts = self.choose_transcripts(prob_dict, self.count)
        for i in range(self.count):
            selected_transcript = randomized_transcripts[i]
            read_object = anno_dict[selected_transcript]
            r = self.simulate_read(read_object, verbose=False)
            if r is not None:
                if r[0] not in reverse_chroms_dict:
                    chroms_dict[chrom_counter] = r[0]
                    if not args.BED:
                        print('#G {} {}'.format(chrom_counter,chroms_dict[chrom_counter]))
                    
                    reverse_chroms_dict[r[0]] = chrom_counter
                    chrom_counter += 1
                
                if r[4] not in reverse_names_dict:
                    names_dict[name_counter] = r[4]
                    if not args.BED:
                        print('#N {} {}'.format(name_counter,names_dict[name_counter]))
                    
                    reverse_names_dict[r[4]] = name_counter
                    name_counter += 1
                
                r[0] = reverse_chroms_dict[r[0]]
                r[4] = reverse_names_dict[r[4]]
                elr_line = '\t'.join([str(i) for i in r])
                if args.BED:
                    print(eu.elr_line_to_bed(elr_line,chroms_dict,names_dict))
                else:
                    print(elr_line)
    
    def parse_abundance_file(self, filename):
        """Parses a 2-column TSV file with relative transcript 
        abundances to generate a dict of proportions per transcript."""
        abundance = {}
        abundance_table = open(filename)
        for line in abundance_table:
            if line[0] == '#':
                continue
            
            l = line.rstrip().split('\t')
            abundance[l[0]] = float(l[1])
        
        abundance_table.close()
        total_abundance = sum(abundance.values())
        proportions = {}
        for k,v in abundance.items():
            proportions[k] = v/total_abundance
        
        return proportions
    
    def choose_transcripts(self, prob_dict, read_number):
        """Selects one transcript with a weighted distribution
        based on TPM"""
        return list(
            random.choice(
                list(prob_dict.keys()),
                read_number,
                p=list(prob_dict.values())
            )
        )
    
    def display_options(self):
        """Returns a string describing all input args"""
        options_string = "\n/| bookend simulate |\\\n¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯\n"
        options_string += "  Reference file:        {}\n".format(self.reference)
        options_string += "  Output file:           {}\n".format(self.output)
        options_string += "  *** Parameters ***\n"
        options_string += "  Read length:           {}\n".format(self.read_length)
        options_string += "  Paired-end?            {}\n".format(self.paired)
        options_string += "  Min length:            {}\n".format(self.min_length)
        options_string += "  Number of reads:       {}\n".format(self.count)
        options_string += "  Fragment length mean:  {}\n".format(self.fragment_mean)
        options_string += "  Fragment length sd:    {}\n".format(self.fragment_sd)
        options_string += "  5' adapter length:     {}\n".format(self.adapter_5p)
        options_string += "  3' adapter length:     {}\n".format(self.adapter_3p)
        options_string += "  5' end pos variance:   {}\n".format(self.var_5p)
        options_string += "  3' end pos variance:   {}\n".format(self.var_3p)
        options_string += "  RNA intact [0-100]%:   {}\n".format(self.percent_intact)
        options_string += "  cDNA sense [0-100]%:   {}\n".format(self.percent_sense)
        options_string += "  Label noise [0-100]%:  {}\n".format(self.label_noise)
        options_string += "  RNG seed:              {}\n".format(self.seed)
        return options_string
    
    def display_summary(self):
        summary = '\n'
        summary += "Wrote {} simulated reads to {}.\n".format(self.count, self.output)
        return summary
    
    def simulate_read(self, read_object, verbose=False):
        """Randomly fragment a ReadObject and write ELR
        lines simulating sequencing reads from the object"""
        base_length = read_object.get_length()
        base_positions = read_object.explode()
        if read_object.strand == '-':
            base_positions = base_positions[::-1]
        
        base_5p = base_positions[0]
        base_3p = base_positions[-1]
        var_3p = abs(var_3p)
        var_5p = abs(var_5p)
        
        if var_5p > 0:
            extension_5p = int(round(random.laplace(0,var_5p)))
            if extension_5p > 0:
                # If extension, make an array of extended positions
                if read_object.strand == '+':
                    ext_pos_5p = list(range(base_5p-extension_5p,base_5p))
                else:
                    ext_pos_5p = list(range(base_5p+extension_5p,base_5p,-1))
            
                offset_5p = 0
            else:
                ext_pos_5p = []
                offset_5p = -extension_5p
        else:
            offset_5p = 0
            ext_pos_5p = []

        if var_3p > 0:
            extension_3p = int(round(random.laplace(0,var_3p)))
            if extension_3p > 0:
                # If extension, make an array of extended positions
                if read_object.strand == '+':
                    ext_pos_3p = list(range(base_3p+1,base_3p+extension_3p+1))
                else:
                    ext_pos_3p = list(range(base_3p-1,base_3p-extension_3p-1,-1))
                
                offset_3p = 0
            else:
                ext_pos_3p = []
                offset_3p = extension_3p
        else:
            offset_3p = 0
            ext_pos_3p = []
        
        if offset_3p < 0:
            adjusted_positions = ext_pos_5p+base_positions[offset_5p:offset_3p]+ext_pos_3p
        else:
            adjusted_positions = ext_pos_5p+base_positions[offset_5p:]+ext_pos_3p
        
        adapter_seq_5p = ['S']*adapter_5p
        adapter_seq_3p = ['E']*adapter_3p
        
        if random.randint(0,100) >= percent_intact:
            break_pos = random.randint(0,len(adjusted_positions))
        else:
            break_pos = 0
        
        cdna_positions = adapter_seq_5p + adjusted_positions[break_pos:] + adapter_seq_3p
        fragment_length = int(round(random.normal(fragment_mean,fragment_sd)))
        
        if fragment_length < minlen:
            return
        
        if fragment_length >= len(cdna_positions):
            cdna_fragment = cdna_positions
            fragment_start = 0
        else:
            fragment_start = random.randint(0,len(cdna_positions)-fragment_length)
            cdna_fragment = cdna_positions[fragment_start:(fragment_start+fragment_length)]
        
        if random.randint(0,100) >= percent_sense:
            # cDNA fragment is antisense to the read_object
            cdna_fragment = cdna_fragment[::-1]
            if read_object.strand == '+':
                fragment_strand = '-'
            else:
                fragment_strand = '+'
        else:
            fragment_strand = read_object.strand
        
        if verbose:
            print(cdna_positions)
            print(fragment_strand)
        
        read_mate1 = cdna_fragment[:read_length]
        read_mate2 = []
        if paired:
            read_mate2 = cdna_fragment[-read_length:]
        
        if read_mate1[0] in ['S','E']:
            first_label = read_mate1[0]
        else:
            first_label = '.'
        
        if paired:
            if read_mate2[-1] in ['S','E']:
                last_label = read_mate2[-1]
            else:
                last_label = '.'
        else:
            if read_mate1[-1] in ['S','E']:
                last_label = read_mate1[-1]
            else:
                last_label = '.'
        
        if verbose:
            print(read_mate1)
            print(read_mate2)
        
        mate1_blocks = list(eu.array_to_blocks(read_mate1))
        mate2_blocks = list(eu.array_to_blocks(read_mate2))
        
        mate1_length = sum([b-a+1 for a,b in mate1_blocks])
        mate2_length = sum([b-a+1 for a,b in mate2_blocks])
        if not mate1_length >= minlen and not mate2_length >= minlen:
            # The fragment has less than minlen transcript-matching nucleotides
            return
        
        if first_label in ['S','E'] or last_label in ['S','E'] or len(mate1_blocks) > 1 or len(mate2_blocks) > 1:
            # Flip the fragment's strand back to the correct orientation
            fragment_strand = read_object.strand
        else:
            fragment_strand = '.'

        if fragment_strand == '-':
            # Flip 'first' and 'last' to reflect 'leftmost' and 'rightmost' labels
            first_label,last_label = (last_label,first_label)

        if len(mate2_blocks) > 0 and len(mate1_blocks) > 0:
            combined_blocks = sorted(mate1_blocks + mate2_blocks)
            leftmost = min(su.flatten(combined_blocks))
            # Check if the combined blocks are strictly increasing with no overlap
            if all([sum(su.overlap_type(range_a,range_b))==4 for range_a,range_b in zip(combined_blocks[1:],combined_blocks[:-1])]):
                # A gap exists between mate1 and mate2.
                # Print each ELCIGAR separately, combined by '.GAP.'
                if mate1_blocks[0][0] < mate2_blocks[0][0]:
                    gap_length = mate2_blocks[0][0] - mate1_blocks[-1][-1]
                    LEFT_CIGAR = eu.blocks_to_ELCIGAR(mate1_blocks,fragment_strand,first_label=first_label,gap_type='J')
                    RIGHT_CIGAR = eu.blocks_to_ELCIGAR(mate2_blocks,fragment_strand,last_label=last_label,gap_type='J')
                else:
                    gap_length = mate1_blocks[0][0] - mate2_blocks[-1][-1]
                    LEFT_CIGAR = eu.blocks_to_ELCIGAR(mate2_blocks,fragment_strand,first_label=first_label,gap_type='J')
                    RIGHT_CIGAR = eu.blocks_to_ELCIGAR(mate1_blocks,fragment_strand,last_label=last_label,gap_type='J')
                
                ELCIGAR = LEFT_CIGAR + str(gap_length+1) + RIGHT_CIGAR
            else:
                # Mate1 and mate2 overlap to some degree.
                # Find a single merged set of blocks to convert to ELCIGAR
                new_blocks = eu.collapse_blocks(combined_blocks)
                leftmost = min(su.flatten(new_blocks))
                ELCIGAR = eu.blocks_to_ELCIGAR(new_blocks,fragment_strand,first_label=first_label,last_label=last_label,gap_type='J')
        else:
            # Either mate1 or mate2 are missing
            if len(mate1_blocks) > 0:
                leftmost = min(su.flatten(mate1_blocks))
                ELCIGAR = eu.blocks_to_ELCIGAR(
                    mate1_blocks,
                    fragment_strand,
                    first_label,
                    last_label,
                    gap_type='J'
                )
            elif len(mate2_blocks) > 0:
                leftmost = min(su.flatten(mate2_blocks))
                ELCIGAR = eu.blocks_to_ELCIGAR(
                    mate2_blocks,
                    fragment_strand,
                    first_label,
                    last_label,
                    gap_type='J'
                )
            else:
                return
        
        # Add noise to end labels
        if random.randint(0,100) < label_noise:
            ELCIGAR = random.choice(['S','E','.'],1)[0] + ELCIGAR[1:]
        
        if random.randint(0,100) < label_noise:
            ELCIGAR = ELCIGAR[:-1] + random.choice(['S','E','.'],1)[0]
        
        # Write the simulated read
        return [
            read_object.chrom, # chromosome
            leftmost, # chromStart
            fragment_strand, # strand
            ELCIGAR,
            read_object.readname,
            '1'
        ]

    
if __name__ == "__main__":
    args = vars(parser.parse_args())
    obj = ELRsimulator(args)
    sys.exit(obj.run())
