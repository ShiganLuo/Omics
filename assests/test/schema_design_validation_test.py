#!/usr/bin/env python3
"""Schema validation cases for sample_params.design (scRNAseq).

P0 acceptance for the comparison-design unification plan:
docs(../../.hermes)scrnaseq-compare-unification-plan.md §3.5. Exercises
SchemaValidator.validate() against the scRNAseq schema, including the
additionalProperties recursion added for dynamic-key maps (tissue / 组名 /
对比名).

Run:  python3 assests/test/schema_design_validation_test.py
Exit 0 = all cases pass.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "src", "common", "util"))

from SchemaValidatorUtil import SchemaValidator  # noqa: E402

SCHEMA = os.path.join(ROOT, "config", "scRNAseq.schema.json")


def make_design(design):
    return {"sample_params": {"design": design}}


VALID_DESIGN = {
    "uterus": {
        "groups": {
            "Youth": ["s1", "s2"],
            "Middle": ["s3", "s4"],
            "Aged": ["s5", "s6"],
        },
        "comparisons": {
            "Aging_series": {
                "groups": ["Youth", "Middle", "Aged"],
                "ref": "Youth",
                "coding": "consecutive",
            },
            "Aged_as_ref": {
                "groups": ["Youth", "Middle", "Aged"],
                "ref": "Aged",
                "coding": "vs_ref",
            },
            "Pair_only": {
                "groups": ["Youth", "Aged"],
            },
        },
        "tests": ["pairwise", "trend"],
    }
}

CASES = [
    # (name, design, expect_ok, expect_substring)
    ("valid_full", VALID_DESIGN, True, None),
    ("valid_empty", {}, True, None),
    ("valid_defaults", {
        "t1": {
            "groups": {"A": ["s1"], "B": ["s2"]},
            "comparisons": {"A_vs_B": {"groups": ["A", "B"]}},
        }
    }, True, None),
    ("bad_coding_enum", {
        "t1": {
            "groups": {"A": ["s1"], "B": ["s2"]},
            "comparisons": {"X": {"groups": ["A", "B"], "coding": "successive"}},
        }
    }, False, "coding"),
    ("missing_groups_in_comparison", {
        "t1": {
            "groups": {"A": ["s1"], "B": ["s2"]},
            "comparisons": {"X": {"ref": "A"}},
        }
    }, False, "groups"),
    ("missing_comparisons_in_tissue", {
        "t1": {"groups": {"A": ["s1"], "B": ["s2"]}},
    }, False, "comparisons"),
    ("missing_groups_in_tissue", {
        "t1": {"comparisons": {"X": {"groups": ["A", "B"]}}},
    }, False, "groups"),
    ("group_value_not_list", {
        "t1": {
            "groups": {"A": "s1", "B": ["s2"]},
            "comparisons": {"X": {"groups": ["A", "B"]}},
        }
    }, False, "t1.groups.A"),
    ("group_item_not_str", {
        "t1": {
            "groups": {"A": ["s1", 2], "B": ["s2"]},
            "comparisons": {"X": {"groups": ["A", "B"]}},
        }
    }, False, "t1.groups.A"),
    ("comparison_groups_not_list", {
        "t1": {
            "groups": {"A": ["s1"], "B": ["s2"]},
            "comparisons": {"X": {"groups": "A"}},
        }
    }, False, "comparisons.X.groups"),
    ("bad_tests_enum", {
        "t1": {
            "groups": {"A": ["s1"], "B": ["s2"]},
            "comparisons": {"X": {"groups": ["A", "B"]}},
            "tests": ["pairwise", "anova"],
        }
    }, False, "tests"),
    ("ref_not_str", {
        "t1": {
            "groups": {"A": ["s1"], "B": ["s2"]},
            "comparisons": {"X": {"groups": ["A", "B"], "ref": 3}},
        }
    }, False, "ref"),
]


def main():
    sv = SchemaValidator()
    sv.load(SCHEMA)
    failed = 0
    for name, design, expect_ok, needle in CASES:
        errors = sv.validate(make_design(design))
        ok = not errors if expect_ok else any(needle in e for e in errors)
        status = "PASS" if ok else "FAIL"
        if not ok:
            failed += 1
        print(f"[{status}] {name}")
        for e in errors:
            print(f"        {e}")
    print(f"\n{len(CASES) - failed}/{len(CASES)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
