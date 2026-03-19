import pytest
from genedesign.checkers.homopolymer_checker import HomopolymerChecker

def test_default_limit_pass():
    """Checks that a sequence within the default limit (5) passes."""
    checker = HomopolymerChecker(max_run_length=5)
    checker.initiate()
    
    # Exactly 5 'A's should pass
    passed, match = checker.run("ATCGAAAAATC")
    assert passed is True
    assert match == ""

def test_default_limit_fail():
    """Checks that a sequence exceeding the default limit (5) fails."""
    checker = HomopolymerChecker(max_run_length=5)
    checker.initiate()
    
    # 6 'G's should fail
    passed, match = checker.run("GGGGGGATCG")
    assert passed is False
    assert match == "GGGGGG"

def test_custom_limit():
    """Verifies the checker respects a smaller custom limit."""
    # Limit set to 2: 'AAA' should fail
    checker = HomopolymerChecker(max_run_length=2)
    checker.initiate()
    
    passed, match = checker.run("AAATCG")
    assert passed is False
    assert match == "AAA"

def test_case_insensitivity():
    """Ensures 'aaaaaa' fails just like 'AAAAAA'."""
    checker = HomopolymerChecker(max_run_length=3)
    checker.initiate()
    
    passed, match = checker.run("atcgtttta")
    assert passed is False
    assert match.lower() == "tttt"

@pytest.mark.parametrize("sequence, expected_pass", [
    ("ATCGATCG", True),        # No repeats
    ("AATTAATT", True),        # Short repeats
    ("AAAAAAAAAA", False),     # Massive repeat
    ("", True),                # Empty string
    ("GCGCGCGC", True),        # Alternating
])
def test_various_sequences(sequence, expected_pass):
    """Parametrized test for multiple sequence types."""
    checker = HomopolymerChecker(max_run_length=5)
    checker.initiate()
    passed, _ = checker.run(sequence)
    assert passed == expected_pass