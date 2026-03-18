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
        self.PENALTY_PROMOTER = -500
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
        """
        Refines the sequence by aggressively minimizing hairpins and removing hard motifs.
        Uses local windowing to avoid expensive string operations.
        """
        current_codons = list(codons)
        
        for i in range(len(current_codons)):
            aa = amino_acids[i]
            candidates = self.aminoAcidToCodons.get(aa, [])
            if len(candidates) < 2:
                continue
                
            # Pre-calculate the full string once per position to extract context
            full_seq_str = utr + "".join(current_codons)
            pos = len(utr) + (i * 3)
            
            # Extract +/- 100bp context to handle the maximum checker window
            left_context = full_seq_str[max(0, pos - 100) : pos]
            right_context = full_seq_str[pos + 3 : pos + 103]
            
            scored_candidates = []
            for candidate in candidates:
                # Construct local strings for testing
                motif_win = left_context[-60:] + candidate + right_context[:60]
                hair_win = left_context[-100:] + candidate + right_context[:100]
                
                # 1. Check Hard Constraints (Forbidden/Promoters)
                passed_forb, _ = self.forbidden_checker.run(motif_win)
                passed_prom, _ = self.promoter_checker.run(motif_win)
                passed_homo, _ = self.homopolymer_checker.run(motif_win)
                
                # 2. Count Hairpins (The Primary Optimization Goal)
                _, hair_report = hairpin_checker(hair_win)
                h_count = hair_report.count("Hairpin") if isinstance(hair_report, str) else 0
                
                # 3. Check Global Metrics (Diversity/Rare Codons)
                # We update the list temporarily to check global diversity
                original_at_i = current_codons[i]
                current_codons[i] = candidate
                passed_metrics, _, _, _ = self.codon_checker.run(current_codons)
                current_codons[i] = original_at_i

                # Assign Tier: 
                # Tier 0: Safe (No Promoters/Forbidden) AND Clean (No Hairpins)
                # Tier 1: Safe (No Promoters/Forbidden) BUT has Hairpins
                # Tier 2: Unsafe (Contains Promoter or Forbidden sequence)
                if passed_forb and passed_prom and passed_homo:
                    tier = 0 if h_count == 0 else 1
                else:
                    tier = 2
                
                scored_candidates.append({
                    'codon': candidate,
                    'tier': tier,
                    'h_count': h_count,
                    'weight': self.codon_weights.get(candidate, 0)
                })
            
            # SORTING LOGIC:
            # 1. Tier (Safety First)
            # 2. h_count (Lowest number of hairpins in window)
            # 3. weight (Best CAI as the final tie-breaker)
            scored_candidates.sort(key=lambda x: (x['tier'], x['h_count'], -x['weight']))
            
            # Update the codon at this position with the winner
            current_codons[i] = scored_candidates[0]['codon']
            
        return current_codons
    
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

    def _score_candidate(self, sequence_window: str, variant_codons: list, committed_codons: list) -> tuple[float, bool]:
        score = 0.0
        all_passed = True
        
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
        if all_passed:
            h_count, _ = self.hairpin_counter(sequence_window)
            if h_count > 0:
                score += (h_count * self.PENALTY_HAIRPIN)
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
        
        return score, all_passed
        

    def run(self, peptide: str, ignores: set) -> Transcript:
        full_peptide = peptide if peptide.endswith('*') else peptide + '*'
        
        # Choose RBS
        selectedRBS = self.rbsChooser.run("", ignores)
        utr = selectedRBS.utr.upper()
        
        committed_codons = []
        window_size = 3 # Number of amino acids to look ahead
        
        idx = 0
        while idx < len(full_peptide):
            window_pep = full_peptide[idx : idx + window_size]
            best_candidate = None
            
            for _ in range(10): # tries 10 sequences, uses only the best
                variant = []
                for aa in window_pep:
                    choices = self.aminoAcidToCodons[aa]
                    variant.append(random.choices(choices, weights=[self.codon_weights[c] for c in choices])[0])
                
                context = (utr + "".join(committed_codons))[-50:]
                test_seq = context + "".join(variant)
                
                # Modified scorer returns (score, all_passed_boolean)
                score, all_passed = self._score_candidate(test_seq, variant, committed_codons)
                
                current_candidate = (variant, score)
                
                # Keep track of the 'best' so far in case none pass perfectly
                if best_candidate is None or score > best_candidate[1]:
                    best_candidate = current_candidate
                
                # FIRST-FIT LOGIC: Exit early if this candidate is perfect
                if all_passed:
                    best_candidate = current_candidate
                    break 
            
            committed_codons.append(best_candidate[0][0])
            idx += 1

        #final_codons = self._exhaustive_optimize(utr, committed_codons, full_peptide) 

        return Transcript(selectedRBS, full_peptide, committed_codons)