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
            'A': [("GCG", 0.36), ("GCC", 0.26), ("GCA", 0.21), ("GCT", 0.17)],
            'C': [("TGC", 0.55), ("TGT", 0.45)],
            'D': [("GAT", 0.63), ("GAC", 0.37)],
            'E': [("GAA", 0.68), ("GAG", 0.32)],
            'F': [("TTC", 0.57), ("TTT", 0.43)],
            'G': [("GGT", 0.34), ("GGC", 0.40), ("GGA", 0.11), ("GGG", 0.15)],
            'H': [("CAC", 0.43), ("CAT", 0.57)],
            'I': [("ATC", 0.48), ("ATT", 0.49), ("ATA", 0.03)],
            'K': [("AAA", 0.76), ("AAG", 0.24)],
            'L': [("CTG", 0.50), ("TTA", 0.13), ("TTG", 0.13), ("CTT", 0.11), ("CTC", 0.10), ("CTA", 0.04)],
            'M': [("ATG", 1.0)],
            'N': [("AAC", 0.55), ("AAT", 0.45)],
            'P': [("CCG", 0.52), ("CCA", 0.19), ("CCT", 0.16), ("CCC", 0.12)],
            'Q': [("CAG", 0.65), ("CAA", 0.35)],
            'R': [("CGT", 0.38), ("CGC", 0.40), ("CGA", 0.07), ("CGG", 0.11), ("AGA", 0.02), ("AGG", 0.01)],
            'S': [("AGC", 0.28), ("TCT", 0.15), ("TCC", 0.15), ("TCA", 0.13), ("TCG", 0.15), ("AGT", 0.15)],
            'T': [("ACC", 0.44), ("ACA", 0.13), ("ACT", 0.17), ("ACG", 0.26)],
            'V': [("GTG", 0.38), ("GTT", 0.26), ("GTC", 0.22), ("GTA", 0.15)],
            'W': [("TGG", 1.0)],
            'Y': [("TAT", 0.57), ("TAC", 0.43)],
            '*': [("TAA", 0.61), ("TAG", 0.09), ("TGA", 0.30)]
        }
        
        for aa, data in self.codon_data.items():
            self.aminoAcidToCodons[aa] = [d[0] for d in data]
            for codon, weight in data:
                self.codon_weights[codon] = weight

    def _exhaustive_optimize(self, utr: str, codons: list[str], amino_acids: str) -> list[str]:
        optimized_codons = list(codons)
        
        for i in range(len(optimized_codons)):
            aa = amino_acids[i]
            synonymous_choices = self.aminoAcidToCodons[aa]
            
            if len(synonymous_choices) <= 1:
                continue
                
            # 1. Prepare context: Previous codons + current + future
            # History is everything before index i
            committed_before = optimized_codons[:i]
            # Context string for the checkers (50bp history)
            context_str = (utr + "".join(committed_before))[-50:]
            
            # 2. Score current codon at this position
            # We must join the variant list into a string for the first argument!
            current_variant = [optimized_codons[i]]
            test_seq = context_str + "".join(current_variant)
            
            current_best_score, _, _ = self._score_candidate(test_seq, current_variant, committed_before, utr)
            best_codon_for_pos = optimized_codons[i]
            
            for alternative in synonymous_choices:
                if alternative == optimized_codons[i]:
                    continue
                
                test_variant = [alternative]
                # JOIN into a string to avoid the 'list' object has no attribute 'upper' error
                test_seq = context_str + "".join(test_variant)
                
                # Pass the joined string as the first argument
                test_score, _, _ = self._score_candidate(test_seq, test_variant, committed_before, utr)
                
                if test_score > current_best_score:
                    current_best_score = test_score
                    best_codon_for_pos = alternative
            
            optimized_codons[i] = best_codon_for_pos
            
        return optimized_codons
    
    def _break_hairpin(self, utr: str, committed_codons: list[str], variant_codons: list[str], bad_seq: str, window_pep: str) -> tuple[list[str], float, bool]:
        """
        Targeted Mutation: Identifies the specific codons overlapping the hairpin 
        and swaps them for the best synonymous alternative that breaks the stem.
        """
        best_variant = list(variant_codons)
        context_str = (utr + "".join(committed_codons))[-50:]
        
        # Get baseline
        current_test_seq = context_str + "".join(variant_codons)
        best_score, best_passed, _ = self._score_candidate(current_test_seq, variant_codons, committed_codons, utr)

        # Identify which codons in the variant overlap with the bad_seq
        # We check each codon's position in the joined variant string
        variant_dna = "".join(variant_codons)
        
        for i in range(len(variant_codons)):
            codon_start = i * 3
            codon_seq = variant_codons[i]
            
            # If any part of this codon is inside the bad_seq
            # We use a simple string check, but could be more precise with indices
            if codon_seq in bad_seq:
                aa = window_pep[i]
                choices = self.aminoAcidToCodons[aa]
                
                if len(choices) > 1:
                    # Test every alternative synonymous codon for this specific position
                    for alt_codon in choices:
                        if alt_codon == variant_codons[i]:
                            continue
                        
                        test_variant = list(variant_codons)
                        test_variant[i] = alt_codon
                        
                        test_seq = context_str + "".join(test_variant)
                        test_score, test_passed, _ = self._score_candidate(test_seq, test_variant, committed_codons, utr)
                        
                        # Optimization: If this alternative breaks the hairpin AND 
                        # maintains a better score than our current best, keep it.
                        if test_score > best_score:
                            best_score = test_score
                            best_variant = test_variant
                            best_passed = test_passed
                            
                            # If we fixed the hairpin and passed all other checks, stop
                            if best_passed:
                                return best_variant, best_score, True

        return best_variant, best_score, best_passed
        
    def hairpin_counter(self, dna, min_stem=3, min_loop=4, max_loop=9):
        """ Counts the total number of distinct hairpin structures in a DNA sequence. """
        def get_reverse_complement(seq):
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
            return "".join(complement.get(base, base) for base in reversed(seq))

        hairpin_count = 0
        first_hairpin_found = None
        n = len(dna)
        for i in range(n):
            for stem_len in range(min_stem, (n - max_loop) // 2 + 1):
                stem_a = dna[i : i + stem_len]
                if len(stem_a) < min_stem: break
                for loop_len in range(min_loop, max_loop + 1):
                    start_b = i + stem_len + loop_len
                    end_b = start_b + stem_len
                    if end_b > n: break
                    stem_b = dna[start_b : end_b]
                    if stem_a == get_reverse_complement(stem_b):
                        hairpin_count += 1
                        if not first_hairpin_found: first_hairpin_found = dna[i : end_b]
                        break 
        return hairpin_count, first_hairpin_found

    def _score_candidate(self, sequence_window: str, variant_codons: list, committed_codons: list, utr) -> tuple[float, bool]:
        score = 0.0
        all_passed = True
        bad_seq = None
        
        # 1. FASTEST: Internal Promoters        
        passed_prom, _ = self.promoter_checker.run(sequence_window)
        if not passed_prom:
            score += self.PENALTY_PROMOTER
            all_passed = False
        
        # 2. Fast: Forbidden Sequences
        pass_forb, _ = self.forbidden_checker.run(sequence_window)
        if not pass_forb:
            score += self.PENALTY_FORBIDDEN
            all_passed = False

        
        """
        # 2. HEAVY: Hairpins (Only check if promoters passed to save time)
        full_dna = utr + "".join(committed_codons) + "".join(variant_codons)
        if all_passed:
            passed_hairpin, bad_seq = hairpin_checker(full_dna) # Use the imported checker
            if not passed_hairpin:
                score += self.PENALTY_HAIRPIN
                all_passed = False
        """
        
        # 3. GLOBAL: Codon Metrics (Diversity, Rare Codons, Global CAI)
        full_potential_list = committed_codons + variant_codons
        codons_above_board, _, _, _ = self.codon_checker.run(full_potential_list)
        
        # Update all_passed based on codon metrics
        # We only enforce this after 40 codons so diversity has a chance to build up
        if len(full_potential_list) > 40 and not codons_above_board:
            score += self.PENALTY_CODON_METRICS
            all_passed = False # CRITICAL: Prevents Early Exit if diversity is too low

        # 4. SCORING: CAI and Diversity Bonus (Always calculated for ranking)
        avg_cai = sum(self.codon_weights.get(c, 0) for c in variant_codons) / len(variant_codons)
        score += avg_cai * self.WEIGHT_CAI

        history_set = set(committed_codons)
        new_codons = [c for c in variant_codons if c not in history_set]
        score += len(set(new_codons)) * 50 
        
        return score, all_passed, bad_seq
        

    def run(self, peptide: str, ignores: set) -> Transcript:
        full_peptide = peptide if peptide.endswith('*') else peptide + '*'
        
        # Choose RBS
        selectedRBS = self.rbsChooser.run("", ignores)
        utr = selectedRBS.utr.upper()
        
        committed_codons = []
        window_size = 9 # Number of amino acids to look ahead
        
        idx = 0
        while idx < len(full_peptide):
            remaining_len = len(full_peptide) - idx
            num_to_commit = min(3, remaining_len)
            window_pep = full_peptide[idx : idx + window_size]
            best_candidate = None
            context = (utr + "".join(committed_codons))[-50:]
            
            for _ in range(10): # tries 10 sequences, uses only the best
                variant = []
                for aa in window_pep:
                    choices = self.aminoAcidToCodons[aa]
                    variant.append(random.choices(choices, weights=[self.codon_weights[c] for c in choices])[0])
                test_seq = context + "".join(variant)
                
                # Modified scorer returns (score, all_passed_boolean)
                score, all_passed, bad_seq = self._score_candidate(test_seq, variant, committed_codons, utr)
                
                # tries replacing a codon to break a hairpin
                if not all_passed and bad_seq:
                    variant, score, all_passed = self._break_hairpin(utr, committed_codons, variant, bad_seq, window_pep)

                # Update the best candidate for this position
                current_candidate = (variant, score)
                if best_candidate is None or score > best_candidate[1]:
                    best_candidate = current_candidate
                
                # FIRST-FIT LOGIC: Exit early if this candidate is perfect
                if all_passed:
                    best_candidate = current_candidate
                    break 
            
            for i in range(num_to_commit):
                committed_codons.append(best_candidate[0][i])
            idx += num_to_commit

        #final_codons = self._exhaustive_optimize(utr, committed_codons, full_peptide) 
        return Transcript(selectedRBS, full_peptide, committed_codons)