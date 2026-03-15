import random
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript

# Import all unmodified checkers
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.homopolymer_checker import HomopolymerChecker
from genedesign.checkers.codon_checker import CodonChecker

class TranscriptDesigner:
    def __init__(self):
        self.aminoAcidToCodons = {}
        self.codon_weights = {} # For scoring
        self.rbsChooser = None
        self.forbidden_checker = None
        self.promoter_checker = None
        self.homopolymer_checker = None
        self.codon_checker = None

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()
        
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.forbidden_checker.initiate()
        
        self.promoter_checker = PromoterChecker()
        self.promoter_checker.initiate()
        
        self.homopolymer_checker = HomopolymerChecker(max_run_length=5)
        self.homopolymer_checker.initiate()
        
        self.codon_checker = CodonChecker()
        self.codon_checker.initiate()

        # Codon table with approximate E. coli frequency weights (0.0 to 1.0)
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

    def _get_weighted_shuffled_codons(self, aa: str) -> list[str]:
        codons = self.aminoAcidToCodons.get(aa, []).copy()
        # Multiplying by random.random() creates a weighted lottery. 
        # This keeps CAI high while introducing enough variation to pass the Diversity >= 0.5 threshold.
        codons.sort(key=lambda c: self.codon_weights.get(c, 0) * random.random(), reverse=True)
        return codons

    # current best implementation
    def _optimize_window(self, global_history: str, local_context: str, peptide: str, strict=True) -> list[str]:
        best_codons = []
        eval_counter = [0]
        max_evals = 50

        def dfs(depth, current_local_context, current_codons):
            if eval_counter[0] > max_evals: return False
            if depth == len(peptide):
                eval_counter[0] += 1
                proposed_seq = "".join(current_codons)
                test_seq = global_history + proposed_seq
                
                if strict:
                    passed_hairpin, _ = hairpin_checker(test_seq)
                    if not passed_hairpin: return False
                    
                    passed_promoter, _ = self.promoter_checker.run(test_seq)
                    if not passed_promoter: return False
                
                best_codons.extend(current_codons)
                return True

            amino_acid = peptide[depth]
            available_codons = self._get_weighted_shuffled_codons(amino_acid)
            
            for codon in available_codons:
                if eval_counter[0] > max_evals: break
                test_local_context = current_local_context + codon
                
                # FIX 1: Pass the whole local_context (which we will expand to 50bp) 
                # instead of slicing to [-12:] inside the loop
                passed_forbidden, _ = self.forbidden_checker.run(test_local_context)
                passed_homopolymer, _ = self.homopolymer_checker.run(test_local_context)
                
                if passed_forbidden and passed_homopolymer:
                    if dfs(depth + 1, test_local_context, current_codons + [codon]):
                        return True
            return False

        dfs(0, local_context, [])
        return best_codons

    def run(self, peptide: str, ignores: set) -> Transcript:
        if not peptide.startswith('M'):
            peptide = 'M' + peptide
        full_peptide = peptide if peptide.endswith('*') else peptide + '*'
        
        max_attempts = 15
        local_ignores = set(ignores)
        attempt_history = {}

        for attempt in range(max_attempts):
            try:
                selectedRBS = self.rbsChooser.run("", local_ignores)
                utr = selectedRBS.utr.upper()
                committed_codons = []
                
                for i in range(len(full_peptide)):
                    target_peptide = full_peptide[i : i + 3]
                    full_context = utr + ''.join(committed_codons)
                    
                    # FIX 2: Expand structure history to 150bp, and local forbidden context to 50bp
                    window_codons = self._optimize_window(full_context[-150:], full_context[-50:], target_peptide, strict=True)
                    
                    if not window_codons:
                        raise ValueError("Trap")
                    
                    committed_codons.append(window_codons[0])
                
                codons_above_board, codon_div, rare_count, cai_val = self.codon_checker.run(committed_codons)
                if not codons_above_board:
                    raise ValueError(f"Codon metrics failed: CAI={cai_val:.2f}, Div={codon_div:.2f}")
                
                # FIX 3: Return full_peptide to satisfy the completeness checker
                return Transcript(selectedRBS, full_peptide, committed_codons)

            except ValueError:
                attempt_history[selectedRBS] = len(committed_codons)
                local_ignores.add(selectedRBS)
            except Exception:
                break

        # Emergency Fallback
        best_rbs = max(attempt_history, key=attempt_history.get) if attempt_history else self.rbsChooser.run("", ignores)
        utr = best_rbs.utr.upper()
        committed_codons = []
        
        for i in range(len(full_peptide)):
            target_peptide = full_peptide[i : i + 1]
            full_context = utr + ''.join(committed_codons)
            
            # Apply expanded lookbacks here as well
            window_codons = self._optimize_window(full_context[-150:], full_context[-50:], target_peptide, strict=False)
            
            if not window_codons:
                aa = full_peptide[i]
                committed_codons.append(self._get_weighted_shuffled_codons(aa)[0])
            else:
                committed_codons.append(window_codons[0])
                
        # FIX 3: Return full_peptide here as well
        return Transcript(best_rbs, full_peptide, committed_codons)
     # the best approach so far