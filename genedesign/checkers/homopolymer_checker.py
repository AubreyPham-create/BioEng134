import re

class HomopolymerChecker:
    """
    Scans a DNA sequence for homopolymers (long runs of the same nucleotide).
    """
    def __init__(self, max_run_length: int = 5):
        """
        Args:
            max_run_length: The maximum allowed consecutive identical nucleotides.
                            If set to 5, 'AAAAA' is legal, but 'AAAAAA' (6) fails.
        """
        self.max_run_length = max_run_length
        self.pattern = None

    def initiate(self) -> None:
        """
        Compiles the regex pattern for maximum execution speed.
        The regex \1 matches whatever the first capture group ([ATCG]) found.
        """
        # {max_run_length,} means "repeat the captured character this many times or more"
        self.pattern = re.compile(rf'([ATCG])\1{{{self.max_run_length},}}', re.IGNORECASE)

    def run(self, sequence: str) -> tuple[bool, str]:
        """
        Runs the homopolymer check.
        
        Returns:
            tuple: (Passed (bool), Failed Sequence (str)). Returns (True, "") if safe.
        """
        match = self.pattern.search(sequence)
        
        if match:
            # Match was found, meaning a homopolymer exceeded the limit
            return False, match.group(0)
            
        return True, ""