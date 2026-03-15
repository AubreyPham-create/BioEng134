import random
import math
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript

# Import all biological checkers
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

        self.codon_data = {
            'A': [('GCG', 0.36), ('GCC', 0.26), ('GCA', 0.21), ('GCT', 0.17)],
            'C': [('TGC', 0.55), ('TGT', 0.45)],
            'D': [('GAT', 0.63), ('GAC', 0.37)],
            'E': [('GAA', 0.68), ('GAG', 0.32)],
            'F': [('TTC', 0.57), ('TTT', 0.43)],
            'G': [('GGT', 0.34), ('GGC', 0.40), ('GGA', 0.11), ('GGG', 0.15)],
            'H': [('CAC', 0.43), ('CAT', 0.57)],
            'I': [('ATC', 0.48), ('ATT', 0.49), ('ATA', 0.03)],
            'K': [('AAA', 0.76), ('AAG', 0.24)],
            'L': [('CTG', 0.50), ('TTA', 0.13), ('TTG', 0.13), ('CTT', 0.11), ('CTC', 0.10), ('CTA', 0.04)],
            'M': [('ATG', 1.0)],
            'N': [('AAC', 0.55), ('AAT', 0.45)],
            'P': [('CCG', 0.52), ('CCA', 0.19), ('CCT', 0.16), ('CCC', 0.12)],
            'Q': [('CAG', 0.65), ('CAA', 0.35)],
            'R': [('CGT', 0.38), ('CGC', 0.40), ('CGA', 0.07), ('CGG', 0.11), ('AGA', 0.02), ('AGG', 0.01)],
            'S': [('AGC', 0.28), ('TCT', 0.15), ('TCC', 0.15), ('TCA', 0.13), ('TCG', 0.15), ('AGT', 0.15)],
            'T': [('ACC', 0.44), ('ACA', 0.13), ('ACT', 0.17), ('ACG', 0.26)],
            'V': [('GTG', 0.38), ('GTT', 0.26), ('GTC', 0.22), ('GTA', 0.15)],
            'W': [('TGG', 1.0)],
            'Y': [('TAT', 0.57), ('TAC', 0.43)],
            '*': [('TAA', 0.61), ('TAG', 0.09), ('TGA', 0.30)]
        }
        for aa, data in self.codon_data.items():
            self.aminoAcidToCodons[aa] = data
            for codon, weight in data:
                self.codon_weights[codon] = weight

    def _get_candidate_codons(self, aa: str, count: int = 2) -> list[str]:
        options = self.codon_data.get(aa, [])
        filtered = [c for c, w in options if w >= 0.1]
        if not filtered:
            filtered = [max(options, key=lambda x: x[1])[0]]
        filtered.sort(key=lambda c: self.codon_weights[c], reverse=True)
        return filtered[:count]

    def run(self, peptide: str, ignores: set) -> Transcript:
        if not peptide.startswith('M'):
            peptide = 'M' + peptide
        full_peptide = peptide if peptide.endswith('*') else peptide + '*'
        
        selectedRBS = self.rbsChooser.run('', ignores)
        utr = selectedRBS.utr.upper()
        committed_codons = []
        
        for i in range(0, len(full_peptide), 3):
            window_peptide = full_peptide[i : i + 9]
            commit_size = min(3, len(full_peptide) - i)
            
            best_window_codons = []
            best_score = -float('inf')
            
            current_dna = ''.join(committed_codons)
            context_dna = (utr + current_dna)[-50:]
            
            # Monte Carlo sampling: Try 50 variants to find structural safety
            for _ in range(50):
                trial_codons = []
                trial_weights = []
                for aa in window_peptide:
                    c_options = self._get_candidate_codons(aa, count=3)
                    chosen = random.choice(c_options)
                    trial_codons.append(chosen)
                    trial_weights.append(self.codon_weights[chosen])
                
                trial_dna = ''.join(trial_codons)
                full_test_dna = context_dna + trial_dna
                
                cai_score = sum([math.log(w) for w in trial_weights]) / len(trial_weights) if trial_weights else 0
                penalty = 0
                
                if not self.forbidden_checker.run(full_test_dna)[0]:
                    penalty += 500.0
                if not self.homopolymer_checker.run(full_test_dna)[0]:
                    penalty += 500.0
                
                # Evaluation for hairpins
                passed_hairpin, _ = hairpin_checker(full_test_dna)
                if not passed_hairpin:
                    penalty += 300.0
                    
                passed_promoter, _ = self.promoter_checker.run(full_test_dna)
                if not passed_promoter:
                    penalty += 50.0
                
                score = cai_score - penalty
                if score > best_score:
                    best_score = score
                    best_window_codons = trial_codons
            
            committed_codons.extend(best_window_codons[:commit_size])
            
        return Transcript(selectedRBS, full_peptide, committed_codons)