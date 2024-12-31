
import pytest
from seq import calculate_statistics

def test_a_seq():
    seq = "ACCTGXXCXXGTTACTGGGCXTTGTXX"
    stats = calculate_statistics(seq)  # Only one return value
    assert stats == {"A": 2, "C": 5, "G": 6, "T": 7, "Unknown": 7, "Total": 27}

def test_b_seq():
    seq = "ACCGGGTTTT"
    stats = calculate_statistics(seq)  # Only one return value
    assert stats == {"A": 1, "C": 2, "G": 3, "T": 4, "Unknown": 0, "Total": 10}
