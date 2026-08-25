"""Integration tests for leafcutter2.

These tests run the full pipeline (or sub-steps) against the example data
bundled in the example/ directory of the repository.

Regression tests (test_scenario*) compare output against reference files
generated from the pre-refactoring code (commit 55c3009) and stored in
tests/reference_outputs/.
"""

import gzip
import subprocess
import sys
from pathlib import Path

import pytest

# Paths to example data (relative to repo root)
REPO_ROOT = Path(__file__).parent.parent
EXAMPLE_DIR = REPO_ROOT / "example"
ANNOT_DIR = EXAMPLE_DIR / "annotation"
REF_DIR = Path(__file__).parent / "reference_outputs"
ANNOTATION_GTF = ANNOT_DIR / "chr10.gtf.gz"
ANNOTATION_FA = ANNOT_DIR / "chr10.fa.gz"
JUNCTION_FILES_TXT = EXAMPLE_DIR / "junction_files.txt"


def _make_abs_juncfiles(tmp_path: Path) -> Path:
    """Write a junction_files.txt with absolute paths into tmp_path."""
    abs_junc_txt = tmp_path / "junction_files.txt"
    lines = []
    for line in JUNCTION_FILES_TXT.read_text().splitlines():
        line = line.strip()
        if line:
            lines.append(str(EXAMPLE_DIR / line))
    abs_junc_txt.write_text("\n".join(lines) + "\n")
    return abs_junc_txt


def _make_renamed_juncfiles(tmp_path: Path) -> Path:
    """Write junction_files.txt with absolute paths; first two samples renamed A and B."""
    out = tmp_path / "junction_files_renamed.txt"
    lines = []
    for i, raw in enumerate(JUNCTION_FILES_TXT.read_text().splitlines()):
        raw = raw.strip()
        if not raw:
            continue
        abs_path = str(EXAMPLE_DIR / raw)
        if i == 0:
            lines.append(f"{abs_path}\tA")
        elif i == 1:
            lines.append(f"{abs_path}\tB")
        else:
            lines.append(abs_path)
    out.write_text("\n".join(lines) + "\n")
    return out


def _compare_gz_to_ref(gz_path: Path, ref_txt: Path):
    """Assert decompressed gz content matches reference text file."""
    actual = gzip.open(gz_path, "rt").read()
    expected = ref_txt.read_text()
    assert actual == expected, (
        f"Output {gz_path.name} differs from reference {ref_txt}.\n"
        f"First differing line: "
        + next(
            (f"actual={a!r}  expected={e!r}" for a, e in zip(actual.splitlines(), expected.splitlines()) if a != e),
            "(files have different line counts)",
        )
    )


# ── Session fixture: run scenario 1 once, shared with scenario 2 ─────────────

@pytest.fixture(scope="session")
def scenario1_rundir(tmp_path_factory):
    """Run scenario 1 once per test session; shared with scenario 2 for pickle reuse."""
    rundir = tmp_path_factory.mktemp("scenario1")
    juncfiles = _make_abs_juncfiles(rundir)
    result = subprocess.run(
        [
            "leafcutter2",
            "-j", str(juncfiles),
            "-A", str(ANNOTATION_GTF),
            "-G", str(ANNOTATION_FA),
            "-P", "-L",
            "-r", str(rundir),
        ],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, (
        f"scenario 1 pipeline failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    return rundir


# ── Basic sanity tests ────────────────────────────────────────────────────────

def test_package_imports():
    """The package and all sub-modules import without error."""
    import importlib
    for module in [
        "leafcutter2",
        "leafcutter2.cli",
        "leafcutter2.star_utils",
    ]:
        importlib.import_module(module)


def test_entry_points_help():
    """All registered entry-point commands respond to --help."""
    for cmd in ["leafcutter2", "leafcutter2-make-clusters", "leafcutter2-star2junc", "leafcutter2-transcript-tools"]:
        result = subprocess.run(
            [cmd, "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, (
            f"`{cmd} --help` exited with code {result.returncode}.\n"
            f"stderr: {result.stderr}"
        )


def test_make_clusters(tmp_path):
    """leafcutter2-make-clusters runs on the example data and produces a clusters file."""
    abs_junc_txt = _make_abs_juncfiles(tmp_path)
    result = subprocess.run(
        [
            "leafcutter2-make-clusters",
            "-j", str(abs_junc_txt),
            "-r", str(tmp_path),
            "-o", "leafcutter2",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"leafcutter2-make-clusters failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    clusters_file = tmp_path / "clustering" / "leafcutter2_clusters"
    assert clusters_file.exists(), f"Expected clusters file not found: {clusters_file}"
    assert clusters_file.stat().st_size > 0, "Clusters file is empty"


# ── Regression tests: output must match pre-refactoring reference ─────────────

@pytest.mark.skipif(
    not ANNOTATION_GTF.exists() or not ANNOTATION_FA.exists(),
    reason="Example annotation files not found",
)
def test_scenario1_standard_with_pickle_and_leafcutter1(scenario1_rundir):
    """Standard run with -P and -L; output must match pre-refactoring reference."""
    _compare_gz_to_ref(
        scenario1_rundir / "leafcutter2.cluster_ratios.gz",
        REF_DIR / "test1/cluster_ratios.txt",
    )
    _compare_gz_to_ref(
        scenario1_rundir / "leafcutter2.junction_counts.gz",
        REF_DIR / "test1/junction_counts.txt",
    )


@pytest.mark.skipif(
    not ANNOTATION_GTF.exists() or not ANNOTATION_FA.exists(),
    reason="Example annotation files not found",
)
def test_scenario2_renamed_samples(scenario1_rundir, tmp_path):
    """First two samples renamed A/B; pickle from scenario 1 is reused."""
    juncfiles = _make_renamed_juncfiles(tmp_path)
    result = subprocess.run(
        [
            "leafcutter2",
            "-j", str(juncfiles),
            "-A", str(ANNOTATION_GTF),
            "-G", str(ANNOTATION_FA),
            "-P",
            "-r", str(scenario1_rundir),
            "-o", "leafcutter2_renamed",
        ],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, (
        f"scenario 2 pipeline failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    _compare_gz_to_ref(
        scenario1_rundir / "leafcutter2_renamed.cluster_ratios.gz",
        REF_DIR / "test2/cluster_ratios.txt",
    )
    _compare_gz_to_ref(
        scenario1_rundir / "leafcutter2_renamed.junction_counts.gz",
        REF_DIR / "test2/junction_counts.txt",
    )


@pytest.mark.skipif(
    not (ANNOT_DIR / "chr10_minimal_w_CDS.gtf.gz").exists() or not ANNOTATION_FA.exists(),
    reason="Example annotation files not found",
)
def test_scenario3_minimal_cds_gtf(tmp_path):
    """Minimal CDS GTF annotation; skip clustering via pre-built counts file."""
    counts_file = REF_DIR / "test1/leafcutter1_files/leafcutter2_perind.counts.gz"
    result = subprocess.run(
        [
            "leafcutter2",
            "--leafcutter1-counts-file", str(counts_file),
            "-A", str(ANNOT_DIR / "chr10_minimal_w_CDS.gtf.gz"),
            "-G", str(ANNOTATION_FA),
            "-r", str(tmp_path),
        ],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, (
        f"scenario 3 pipeline failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    _compare_gz_to_ref(
        tmp_path / "leafcutter2.cluster_ratios.gz",
        REF_DIR / "test3/cluster_ratios.txt",
    )
    _compare_gz_to_ref(
        tmp_path / "leafcutter2.junction_counts.gz",
        REF_DIR / "test3/junction_counts.txt",
    )


@pytest.mark.skipif(
    not (ANNOT_DIR / "chr10_minimal.gtf.gz").exists() or not ANNOTATION_FA.exists(),
    reason="Example annotation files not found",
)
def test_scenario4_minimal_gtf_no_cds(tmp_path):
    """Minimal GTF without CDS annotation; skip clustering via pre-built counts file."""
    counts_file = REF_DIR / "test1/leafcutter1_files/leafcutter2_perind.counts.gz"
    result = subprocess.run(
        [
            "leafcutter2",
            "--leafcutter1-counts-file", str(counts_file),
            "-A", str(ANNOT_DIR / "chr10_minimal.gtf.gz"),
            "-G", str(ANNOTATION_FA),
            "-r", str(tmp_path),
        ],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, (
        f"scenario 4 pipeline failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    _compare_gz_to_ref(
        tmp_path / "leafcutter2.cluster_ratios.gz",
        REF_DIR / "test4/cluster_ratios.txt",
    )
    _compare_gz_to_ref(
        tmp_path / "leafcutter2.junction_counts.gz",
        REF_DIR / "test4/junction_counts.txt",
    )


# ── leafcutter2-transcript-tools ─────────────────────────────────────────────
#
# The tests above only ever invoke transcript-tools via `--help`, which returns
# before any input is read. Its GTF-input path shells out to `bedparse gtf2bed`,
# so a broken or missing bedparse went unnoticed until it failed at runtime for
# users -- see the setuptools<81 constraint on bedparse's pkg_resources import.

def test_transcript_tools_gtf_to_bed12(tmp_path):
    """transcript-tools converts the example GTF to bed12 (exercises bedparse)."""
    out_bed = tmp_path / "transcripts.bed.gz"
    result = subprocess.run(
        [
            "leafcutter2-transcript-tools",
            "-i", str(ANNOTATION_GTF),
            "-fa", str(ANNOTATION_FA),
            "-bed12_out", str(out_bed),
        ],
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert result.returncode == 0, (
        f"transcript-tools failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    assert out_bed.exists(), "no bed12 output was written"

    lines = gzip.open(out_bed, "rt").read().splitlines()
    assert lines, "bed12 output is empty"

    header = lines[0].lstrip("#").split("\t")
    for expected in ("chr", "start", "end", "name", "gene_id", "transcript_id",
                     "transcript_type", "NMDFinderB"):
        assert expected in header, f"missing column {expected!r} in {header}"

    assert len(lines) > 1, "bed12 output has a header but no transcripts"
    # Every record must have as many fields as the header.
    for line in lines[1:]:
        assert len(line.split("\t")) == len(header), (
            f"field count mismatch:\n{line}"
        )


def test_transcript_tools_gtf_to_gtf(tmp_path):
    """transcript-tools round-trips the example GTF to an annotated GTF."""
    out_gtf = tmp_path / "annotated.gtf"
    result = subprocess.run(
        [
            "leafcutter2-transcript-tools",
            "-i", str(ANNOTATION_GTF),
            "-fa", str(ANNOTATION_FA),
            "-o", str(out_gtf),
        ],
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert result.returncode == 0, (
        f"transcript-tools failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    assert out_gtf.exists() and out_gtf.stat().st_size > 0, "no GTF output was written"
    assert any(
        "\ttranscript\t" in line for line in out_gtf.read_text().splitlines()[:2000]
    ), "output GTF has no transcript features"
