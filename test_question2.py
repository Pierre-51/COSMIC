"""
Question 2.3 – Pytest tests for q2_analysis.py functions
Run with:  pytest test_q2_analysis.py -v
"""

import pandas as pd
import pytest
from question2_analysis import get_allele_mutation_counts, get_sample_mutation_extremes

# Fixtures

@pytest.fixture
def simple_df():
    """Minimal DataFrame that mimics the real file's relevant columns."""
    data = {
        "icgc_mutation_id": [
            "MUT_1", "MUT_1",   # same mutation, two transcripts  → count once
            "MUT_2",
            "MUT_3", "MUT_3",   # same mutation, two transcripts  → count once
            "MUT_4",
        ],
        "mutated_from_allele": ["A", "A",  "C", "G", "G",  "T"],
        "mutated_to_allele":   ["T", "T",  "G", "A", "A",  "C"],
        "icgc_sample_id": [
            "SA_001", "SA_001",
            "SA_001",
            "SA_002", "SA_002",
            "SA_002",
        ],
    }
    return pd.DataFrame(data)


@pytest.fixture
def multi_sample_df():
    """DataFrame with three samples having different mutation counts."""
    data = {
        "icgc_mutation_id": [
            "MUT_A", "MUT_B", "MUT_C",         # SA_HIGH  → 3
            "MUT_D",                             # SA_LOW   → 1
            "MUT_E", "MUT_F",                   # SA_MID   → 2
        ],
        "mutated_from_allele": ["A", "C", "G", "T", "A", "C"],
        "mutated_to_allele":   ["T", "G", "A", "C", "G", "T"],
        "icgc_sample_id": [
            "SA_HIGH", "SA_HIGH", "SA_HIGH",
            "SA_LOW",
            "SA_MID",  "SA_MID",
        ],
    }
    return pd.DataFrame(data)

# Tests for get_allele_mutation_counts

class TestGetAlleleMutationCounts:

    def test_returns_dataframe(self, simple_df):
        result = get_allele_mutation_counts(simple_df)
        assert isinstance(result, pd.DataFrame)

    def test_expected_columns(self, simple_df):
        result = get_allele_mutation_counts(simple_df)
        assert set(result.columns) == {
            "mutated_from_allele",
            "mutated_to_allele",
            "count_unique_icgc_mutation_id",
        }

    def test_duplicate_mutations_counted_once(self, simple_df):
        """MUT_1 appears twice (two transcripts) – must be counted as 1."""
        result = get_allele_mutation_counts(simple_df)
        a_to_t = result[
            (result["mutated_from_allele"] == "A") &
            (result["mutated_to_allele"]   == "T")
        ]
        assert not a_to_t.empty, "Expected A→T row to exist"
        assert a_to_t["count_unique_icgc_mutation_id"].iloc[0] == 1

    def test_correct_count_for_each_pattern(self, simple_df):
        """Check all four allele patterns return count = 1."""
        result = get_allele_mutation_counts(simple_df)
        # Every pattern in our fixture has exactly 1 unique mutation
        assert (result["count_unique_icgc_mutation_id"] == 1).all()

    def test_number_of_patterns(self, simple_df):
        """Fixture contains 4 distinct allele-change patterns."""
        result = get_allele_mutation_counts(simple_df)
        assert len(result) == 4

    def test_empty_dataframe(self):
        """Empty input should return empty result without error."""
        empty = pd.DataFrame(columns=[
            "icgc_mutation_id",
            "mutated_from_allele",
            "mutated_to_allele",
            "icgc_sample_id",
        ])
        result = get_allele_mutation_counts(empty)
        assert len(result) == 0

    def test_all_same_pattern(self):
        """All rows share the same allele pattern."""
        df = pd.DataFrame({
            "icgc_mutation_id":    ["MUT_1", "MUT_2", "MUT_3"],
            "mutated_from_allele": ["A", "A", "A"],
            "mutated_to_allele":   ["T", "T", "T"],
            "icgc_sample_id":      ["SA_1", "SA_1", "SA_1"],
        })
        result = get_allele_mutation_counts(df)
        assert len(result) == 1
        assert result["count_unique_icgc_mutation_id"].iloc[0] == 3


# ──────────────────────────────────────────────
# Tests for get_sample_mutation_extremes
# ──────────────────────────────────────────────

class TestGetSampleMutationExtremes:

    def test_returns_dict_with_correct_keys(self, multi_sample_df):
        result = get_sample_mutation_extremes(multi_sample_df)
        assert "highest" in result
        assert "lowest"  in result
        for key in ("highest", "lowest"):
            assert "sample_id" in result[key]
            assert "count"     in result[key]

    def test_highest_sample(self, multi_sample_df):
        result = get_sample_mutation_extremes(multi_sample_df)
        assert result["highest"]["sample_id"] == "SA_HIGH"
        assert result["highest"]["count"]     == 3

    def test_lowest_sample(self, multi_sample_df):
        result = get_sample_mutation_extremes(multi_sample_df)
        assert result["lowest"]["sample_id"] == "SA_LOW"
        assert result["lowest"]["count"]     == 1

    def test_count_is_integer(self, multi_sample_df):
        result = get_sample_mutation_extremes(multi_sample_df)
        assert isinstance(result["highest"]["count"], int)
        assert isinstance(result["lowest"]["count"],  int)

    def test_duplicate_mutation_sample_pairs_counted_once(self):
        """Same (sample, mutation) pair appearing twice must count as 1."""
        df = pd.DataFrame({
            "icgc_mutation_id": ["MUT_X", "MUT_X", "MUT_Y"],
            "mutated_from_allele": ["A", "A", "G"],
            "mutated_to_allele":   ["T", "T", "C"],
            "icgc_sample_id": ["SA_1", "SA_1", "SA_2"],
        })
        result = get_sample_mutation_extremes(df)
        # SA_1 has 1 unique mutation (MUT_X counted once), SA_2 has 1 → tie
        assert result["highest"]["count"] == 1
        assert result["lowest"]["count"]  == 1

    def test_single_sample(self):
        """Only one sample → it is both highest and lowest."""
        df = pd.DataFrame({
            "icgc_mutation_id":    ["MUT_1", "MUT_2"],
            "mutated_from_allele": ["A", "C"],
            "mutated_to_allele":   ["T", "G"],
            "icgc_sample_id":      ["SA_ONLY", "SA_ONLY"],
        })
        result = get_sample_mutation_extremes(df)
        assert result["highest"]["sample_id"] == "SA_ONLY"
        assert result["lowest"]["sample_id"]  == "SA_ONLY"
        assert result["highest"]["count"] == 2