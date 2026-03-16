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
        codons.sort(key=lambda c: self.codon_weights.get(c, 0) * random.random(), reverse=True)
        return codons

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
    
    def minimize_hairpins(self, utr: str, codons: list[str], amino_acids: str, passes: int = 3) -> list[str]:
        """
        Targetedly removes hairpins by performing multiple sweeps of the sequence.
        Prioritizes hairpin reduction over codon usage (CAI).
        """
        current_codons = list(codons)
        
        for p in range(passes):
            changes_made = 0
            for i in range(len(current_codons)):
                # 1. Setup Context with a large 500bp window
                full_seq_str = utr + "".join(current_codons)
                pos = len(utr) + (i * 3)
                win_start = max(0, pos - 250)
                win_end = pos + 253 
                
                current_win_seq = full_seq_str[win_start : win_end]
                has_hairpin, report = hairpin_checker(current_win_seq)
                
                if not has_hairpin:
                    continue
                
                current_h_count = report.count("Hairpin") if isinstance(report, str) else 1
                
                # 2. Evaluate Synonymous Alternatives
                aa = amino_acids[i]
                candidates = self.aminoAcidToCodons.get(aa, [])
                if len(candidates) < 2: continue
                
                left_side = full_seq_str[win_start : pos]
                right_side = full_seq_str[pos + 3 : win_end]
                
                best_candidate = current_codons[i]
                min_h_count = current_h_count
                best_weight = self.codon_weights.get(best_candidate, 0)
                
                for candidate in candidates:
                    if candidate == current_codons[i]: continue
                    test_win = left_side + candidate + right_side
                    
                    # Safety Check: No new forbidden sequences or promoters
                    motif_win = test_win[max(0, (pos-win_start)-50) : (pos-win_start)+53]
                    pass_f, _ = self.forbidden_checker.run(motif_win)
                    pass_p, _ = self.promoter_checker.run(motif_win)
                    
                    if pass_f and pass_p:
                        _, test_report = hairpin_checker(test_win)
                        test_h_count = test_report.count("Hairpin") if isinstance(test_report, str) else 0
                        test_weight = self.codon_weights.get(candidate, 0)
                        
                        # PRIORITY: Lower Hairpin Count > Higher CAI
                        if test_h_count < min_h_count:
                            min_h_count = test_h_count
                            best_candidate = candidate
                            best_weight = test_weight
                        elif test_h_count == min_h_count and test_weight > best_weight:
                            best_candidate = candidate
                            best_weight = test_weight
                
                if best_candidate != current_codons[i]:
                    current_codons[i] = best_candidate
                    changes_made += 1
            
            if changes_made == 0: break
                
        return current_codons

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
                    window_codons = self._optimize_window(full_context[-150:], full_context[-50:], target_peptide, strict=True)
                    if not window_codons:
                        raise ValueError("Trap")
                    committed_codons.append(window_codons[0])
                
                # Perform Exhaustive Ranking Optimization
                final_codons = self._exhaustive_optimize(utr, committed_codons, full_peptide)
                
                codons_above_board, codon_div, _, cai_val = self.codon_checker.run(final_codons)
                if not codons_above_board:
                    raise ValueError(f"Codon metrics failed: CAI={cai_val:.2f}, Div={codon_div:.2f}")
                
                return Transcript(selectedRBS, full_peptide, final_codons)

            except ValueError:
                attempt_history[selectedRBS] = len(committed_codons)
                local_ignores.add(selectedRBS)
            except Exception:
                break

        # Emergency Fallback
        best_rbs = max(attempt_history, key=attempt_history.get) if attempt_history else self.rbsChooser.run("", ignores)
        utr, committed_codons = best_rbs.utr.upper(), []
        for i in range(len(full_peptide)):
            target_peptide = full_peptide[i : i + 1]
            full_context = utr + ''.join(committed_codons)
            window_codons = self._optimize_window(full_context[-150:], full_context[-50:], target_peptide, strict=False)
            committed_codons.append(window_codons[0] if window_codons else self._get_weighted_shuffled_codons(full_peptide[i])[0])
            
        minized_hairpins = self.minimize_hairpins(utr, committed_codons, full_peptide)
        final_codons = self._exhaustive_optimize(utr, minized_hairpins, full_peptide)
        return Transcript(best_rbs, full_peptide, final_codons)
        # best implementation so far