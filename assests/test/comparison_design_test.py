#!/usr/bin/env python3
"""Unit cases for common.smk comparison_design() (P1 acceptance).

Exercises the four codings, both views, derived-key naming and the semantic
validations of modules/common/common.smk::comparison_design. The function is
extracted from the .smk source via ast (common.smk is a snakemake include,
not an importable module).

Run:  python3 assests/test/comparison_design_test.py
Exit 0 = all cases pass.
"""
import ast
import os
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
COMMON_SMK = os.path.join(ROOT, "modules", "common", "common.smk")

WANTED = {"comparison_design", "write_compare_design", "_COMPARISON_CODINGS"}


def load_functions():
    src = open(COMMON_SMK, encoding="utf-8").read()
    tree = ast.parse(src)
    ns = {"os": os}
    for node in tree.body:
        keep = (
            (isinstance(node, ast.FunctionDef) and node.name in WANTED)
            or (isinstance(node, ast.Assign)
                and any(getattr(t, "id", None) in WANTED for t in node.targets))
        )
        if keep:
            exec(compile(ast.Module([node], []), COMMON_SMK, "exec"), ns)
    missing = WANTED - set(ns)
    assert not missing, f"missing definitions in common.smk: {missing}"
    return ns


NS = load_functions()
comparison_design = NS["comparison_design"]

DESIGN = {
    "uterus": {
        "groups": {
            "Youth": ["s1", "s2"],
            "Middle": ["s3", "s4"],
            "Aged": ["s5", "s6"],
            "Control": ["s7"],
            "Drug": ["s8"],
        },
        "comparisons": {},
        "tests": ["pairwise"],
    }
}
TISSUE_SAMPLES = {"uterus": [f"s{i}" for i in range(1, 9)]}


def entry(cname, groups, **kw):
    e = {"groups": groups}
    e.update(kw)
    return {"uterus": dict(DESIGN["uterus"], comparisons={cname: e})}


CASES = []


def case(fn):
    CASES.append(fn)
    return fn


@case
def consecutive_three():
    got = comparison_design(
        entry("A", ["Youth", "Middle", "Aged"]), "uterus", "pairwise")
    assert sorted(got) == ["A.Middle_vs_Aged", "A.Youth_vs_Middle"], got
    e = got["A.Youth_vs_Middle"]
    assert e["control_group_name"] == "Youth"
    assert e["experimental_group_name"] == "Middle"
    assert e["control_samples"] == ["s1", "s2"]
    assert e["experimental_samples"] == ["s3", "s4"]


@case
def vs_ref_default_first():
    got = comparison_design(
        entry("A", ["Youth", "Middle", "Aged"], coding="vs_ref"),
        "uterus", "pairwise")
    assert sorted(got) == ["A.Youth_vs_Aged", "A.Youth_vs_Middle"], got


@case
def vs_ref_explicit_nonfirst():
    got = comparison_design(
        entry("Aged_as_ref", ["Youth", "Middle", "Aged"],
              coding="vs_ref", ref="Aged"),
        "uterus", "pairwise")
    assert sorted(got) == ["Aged_as_ref.Aged_vs_Middle",
                           "Aged_as_ref.Aged_vs_Youth"], got
    e = got["Aged_as_ref.Aged_vs_Youth"]
    assert e["control_group_name"] == "Aged"
    assert e["control_samples"] == ["s5", "s6"]


@case
def all_pairs_three():
    got = comparison_design(
        entry("A", ["Youth", "Middle", "Aged"], coding="all_pairs"),
        "uterus", "pairwise")
    assert sorted(got) == ["A.Middle_vs_Aged", "A.Youth_vs_Aged",
                           "A.Youth_vs_Middle"], got


@case
def helmert_three():
    got = comparison_design(
        entry("A", ["Youth", "Middle", "Aged"], coding="helmert"),
        "uterus", "pairwise")
    assert sorted(got) == ["A.Aged_vs_prev_mean", "A.Middle_vs_prev_mean"], got
    e = got["A.Aged_vs_prev_mean"]
    assert e["control_group_name"] == "prev_mean"
    assert e["control_samples"] == ["s1", "s2", "s3", "s4"]
    assert e["experimental_group_name"] == "Aged"
    assert e["experimental_samples"] == ["s5", "s6"]


@case
def two_group_pairwise():
    got = comparison_design(
        entry("DrugCheck", ["Control", "Drug"], coding="vs_ref"),
        "uterus", "pairwise")
    assert list(got) == ["DrugCheck.Control_vs_Drug"], got
    e = got["DrugCheck.Control_vs_Drug"]
    assert e["control_samples"] == ["s7"]
    assert e["experimental_samples"] == ["s8"]


@case
def native_multigroup():
    got = comparison_design(
        entry("Aging_series", ["Youth", "Middle", "Aged"], ref="Youth"),
        "uterus", "native")
    assert list(got) == ["Aging_series"], got
    groups = got["Aging_series"]["groups"]
    assert [g["name"] for g in groups] == ["Youth", "Middle", "Aged"]
    assert groups[1]["samples"] == ["s3", "s4"]


@case
def native_pairwise_group_pairs_shape():
    got = comparison_design(
        entry("DrugCheck", ["Control", "Drug"]), "uterus", "native")
    e = got["DrugCheck"]
    assert e["control_group_name"] == "Control"
    assert e["experimental_group_name"] == "Drug"


@case
def native_multigroup_ref_must_be_first():
    try:
        comparison_design(
            entry("A", ["Youth", "Middle", "Aged"], ref="Aged"),
            "uterus", "native")
    except ValueError as e:
        assert "首柱即基线" in str(e)
        return
    raise AssertionError("expected ValueError")


@case
def error_unknown_tissue():
    try:
        comparison_design(DESIGN, "ovaries", "pairwise")
    except ValueError as e:
        assert "ovaries" in str(e)
        return
    raise AssertionError("expected ValueError")


@case
def error_unknown_group():
    try:
        comparison_design(entry("A", ["Youth", "Ghost"]), "uterus", "pairwise")
    except ValueError as e:
        assert "Ghost" in str(e)
        return
    raise AssertionError("expected ValueError")


@case
def error_sample_in_two_groups():
    d = entry("A", ["Youth", "Middle"])
    d["uterus"]["groups"]["Middle"] = ["s1", "s3"]  # s1 也属 Youth
    try:
        comparison_design(d, "uterus", "pairwise")
    except ValueError as e:
        assert "s1" in str(e)
        return
    raise AssertionError("expected ValueError")


@case
def error_cross_tissue_sample():
    d = entry("A", ["Youth", "Middle"])
    d["uterus"]["groups"]["Middle"] = ["s3", "other_tissue_sample"]
    try:
        comparison_design(d, "uterus", "pairwise",
                          tissue_samples=TISSUE_SAMPLES)
    except ValueError as e:
        assert "other_tissue_sample" in str(e)
        return
    raise AssertionError("expected ValueError")


@case
def error_derived_key_collision():
    # 复合键 {对比名}.{派生名} 下唯一能构造的碰撞:组名带点。
    # comp "X" 派生 "X.A.B_vs_C";comp "X.A" 派生 "X.A.B_vs_C"。
    d = {
        "uterus": {
            "groups": {"A.B": ["s1"], "B": ["s2"], "C": ["s3"]},
            "comparisons": {
                "X": {"groups": ["A.B", "C"], "coding": "consecutive"},
                "X.A": {"groups": ["B", "C"], "coding": "consecutive"},
            },
        }
    }
    try:
        comparison_design(d, "uterus", "pairwise")
    except ValueError as e:
        assert "collides" in str(e)
        return
    raise AssertionError("expected ValueError")


@case
def error_bad_coding():
    try:
        comparison_design(entry("A", ["Youth", "Middle"], coding="successive"),
                          "uterus", "pairwise")
    except ValueError as e:
        assert "successive" in str(e)
        return
    raise AssertionError("expected ValueError")


def main():
    failed = 0
    for fn in CASES:
        try:
            fn()
            print(f"[PASS] {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"[FAIL] {fn.__name__}: {e}")
    print(f"\n{len(CASES) - failed}/{len(CASES)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
