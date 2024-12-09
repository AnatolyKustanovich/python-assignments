
import pytest
from seq import calculate_statistics

def test_a_seq():
    seq = "ACCTGXXCXXGTTACTGGGCXTTGTXX"
    stats, percentages, total = calculate_statistics(seq)
    assert stats == {"A": 2, "C": 5, "G": 6, "T": 7, "Unknown": 7}
    assert percentages == {"A": 2, "C": 5, "G": 6, "T": 7, "Unknown": 7}
    assert total == 27

def test_b_seq():
    seq = "ACCGGGTTTT"
    stats, percentages, total = calculate_statistics(seq)
    assert stats == {"A": 1, "C": 2, "G": 3, "T": 4, "Unknown": 0}
    assert percentages == {"A": 2, "C": 5, "G": 6, "T": 7, "Unknown": 7}
    assert total == 10