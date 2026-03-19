import random
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.homopolymer_checker import HomopolymerChecker
from genedesign.checkers.codon_checker import CodonChecker

class TranscriptDesigner:
    def __init__(self):
        self.aminoAcidToCodons = {}
        self.codon_weights = {}
        self.rbsChooser = None
        self.forbidden_checker = None
        self.promoter_checker = None
        self.homopolymer_checker = None
        self.codon_checker = None
        
        # Scoring Weights
        self.PENALTY_FORBIDDEN = -1000
        self.PENALTY_PROMOTER = -750
        self.PENALTY_HOMOPOLYMER = -200
        self.PENALTY_CODON_METRICS = -500  # Penalty for failing CodonChecker thresholds
        self.PENALTY_HAIRPIN = -100
        self.WEIGHT_CAI = 100  # Multiplier for the 0.0-1.0 codon weight

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser(); self.rbsChooser.initiate()
        self.forbidden_checker = ForbiddenSequenceChecker(); self.forbidden_checker.initiate()
        self.promoter_checker = PromoterChecker(); self.promoter_checker.initiate()
        self.homopolymer_checker = HomopolymerChecker(max_run_length=5); self.homopolymer_checker.initiate()
        self.codon_checker = CodonChecker(); self.codon_checker.initiate()

        self.codon_data = {
        'A': [("GCG", 0.35), ("GCC", 0.27), ("GCA", 0.22), ("GCT", 0.16)], # [cite: 3]
        'C': [("TGC", 0.55), ("TGT", 0.45)], # 
        'D': [("GAT", 0.63), ("GAC", 0.37)], # [cite: 3]
        'E': [("GAA", 0.68), ("GAG", 0.32)], # [cite: 3]
        'F': [("TTT", 0.58), ("TTC", 0.42)], # 
        'G': [("GGC", 0.39), ("GGT", 0.33), ("GGG", 0.16), ("GGA", 0.12)], # [cite: 3]
        'H': [("CAT", 0.58), ("CAC", 0.42)], # [cite: 1, 2]
        'I': [("ATT", 0.50), ("ATC", 0.40), ("ATA", 0.09)], # 
        'K': [("AAA", 0.76), ("AAG", 0.24)], # 
        'L': [("CTG", 0.49), ("TTA", 0.13), ("TTG", 0.13), ("CTT", 0.11), ("CTC", 0.10), ("CTA", 0.04)], # [cite: 1, 2]
        'M': [("ATG", 1.00)], # 
        'N': [("AAC", 0.53), ("AAT", 0.47)], # 
        'P': [("CCG", 0.51), ("CCA", 0.19), ("CCT", 0.17), ("CCC", 0.13)], # [cite: 1, 2]
        'Q': [("CAG", 0.67), ("CAA", 0.33)], # 
        'R': [("CGC", 0.37), ("CGT", 0.36), ("CGG", 0.11), ("CGA", 0.07), ("AGA", 0.05), ("AGG", 0.03)], # [cite: 1, 2]
        'S': [("AGC", 0.27), ("AGT", 0.16), ("TCC", 0.15), ("TCG", 0.15), ("TCT", 0.14), ("TCA", 0.14)], # [cite: 1, 2]
        'T': [("ACC", 0.41), ("ACG", 0.27), ("ACT", 0.17), ("ACA", 0.15)], # 
        'V': [("GTG", 0.37), ("GTT", 0.26), ("GTC", 0.21), ("GTA", 0.16)], # [cite: 3]
        'W': [("TGG", 1.00)], # 
        'Y': [("TAT", 0.57), ("TAC", 0.43)], # 
        '*': [("TAA", 0.59), ("TGA", 0.33), ("TAG", 0.08)] # 
    }
        
        for aa, data in self.codon_data.items():
            self.aminoAcidToCodons[aa] = [d[0] for d in data]
            for codon, weight in data:
                self.codon_weights[codon] = weight
    
    def get_codons(dna_sequence: str) -> list:
        dna = dna_sequence.upper()
        codons = [dna[i : i + 3] for i in range(0, len(dna), 3)]
        if len(codons) > 0 and len(codons[-1]) < 3:
            print(f"Warning: The last fragment '{codons[-1]}' is not a full codon.")
            
        return codons

    
    def _score_candidate(self, sequence_window: str, variant_codons: list, committed_codons: list) -> tuple[bool, str]:
        """
        Acts as a guard clause. Returns True only if all sequence constraints are met.
        """

        # 1. Guard: Internal Promoters
        passed_prom, _ = self.promoter_checker.run(sequence_window)
        if not passed_prom:
            return False

        # 2. Guard: Forbidden Sequences
        passed_forb, _ = self.forbidden_checker.run(sequence_window)
        if not passed_forb:
            return False
        
        # 3. Guard: Homopolymers
        passed_homo, _ = self.homopolymer_checker.run(sequence_window)
        if not passed_homo:
            return False
        
        # 4. Guard: Hairpins (Checked against the full DNA for global safety)
        passed_hairpin, _ = hairpin_checker(sequence_window)
        if not passed_hairpin:
            return False
        
        return True
    
    def select_random_codon(self, codon_sequence: list):
        codon_list = []
        for aa in codon_sequence:
            choices = self.aminoAcidToCodons[aa]
            weights = []
            for c in choices:
                weights.append(self.codon_weights[c])
            selected_codon = random.choices(choices, weights=weights, k=1)[0] # randomly selects based on weights
            codon_list.append(selected_codon)
        return codon_list

    def _design_transcript(self, peptide: str, ignores: set) -> Transcript:
        full_peptide = peptide if peptide.endswith('*') else peptide + '*'
        selectedRBS = self.rbsChooser.run("", ignores)
        utr = selectedRBS.utr.upper()
        
        committed_codons = []
        window_size = 9

        idx = 0
        while idx < len(full_peptide):
            window_pep = full_peptide[idx : idx + window_size]
            
            context = (utr + "".join(committed_codons))[-6:]
            best_legal_variant = None
            
            for i in range(1_000_001):
                variant = self.select_random_codon(window_pep)
                
                test_seq = context + "".join(variant)
                
                # Use checkers as guard clauses
                success = self._score_candidate(test_seq, variant, committed_codons)

                # If the candidate (original or repaired) passes all guards, we take it!
                if success:
                    best_legal_variant = variant
                    break
                if i == 1_000_000:
                    print("Attempt failed: restarting...")
                    return self.run(peptide)
            
            committed_codons.append(best_legal_variant[0])
            idx += 1
            

        return Transcript(selectedRBS, full_peptide, committed_codons)
    
    def run(self, peptide, ignores=set(), max_retries=10):
        """
        Supervisor function: Calls the design logic and validates the 
        final product against the global checkers.
        """
        attempt = 0
        while (True):
            # 1. Generate a candidate using the original logic
            # (Renamed from run to _design_transcript)
            transcript = self._design_transcript(peptide, ignores)
            
            # 2. Combine RBS and Codons for a full sequence audit
            # The benchmark fails often at this junction!
            full_dna = (transcript.rbs.utr + "".join(transcript.codons)).upper()
            
            # 3. Explicitly unpack the tuples (Boolean, Message)
            # Using [0] ensures we get the True/False value, not the 'Truthy' tuple
            forbidden_pass = self.forbidden_checker.run(full_dna)[0]
            promoter_pass = self.promoter_checker.run(full_dna)[0]
            hairpin_pass = hairpin_checker(full_dna)[0]
            
            if forbidden_pass and promoter_pass and hairpin_pass:
                return transcript

            print(f"Attempt {attempt + 1} failed global checks. Retrying...")
            attempt += 1

        raise Exception(f"Could not produce a valid transcript for {peptide[:10]}... after {max_retries} retries.")
        


if __name__ == "__main__":
    # 1. Initialize the designer
    designer = TranscriptDesigner()
    designer.initiate()

    # 2. Define a test amino acid sequence
    # 'A' (Alanine) has 4 codons with different weights
    test_peptide = "MAAAAAAAAAAA" 
    for i in range(1):
        codons = designer.select_random_codon(test_peptide)
        print(codons)

