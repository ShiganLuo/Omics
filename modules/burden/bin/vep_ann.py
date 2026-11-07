#!/usr/bin/env python
"""Cache-based Ensembl VEP annotation wrapper for mutation-burden inputs.

Wraps the ``vep`` command-line tool: offline cache annotation of a VCF
(``--everything --pick --per_gene``), with automatic cache installation via
``vep_install`` when the species cache is missing. Used by ``burden.py`` to
annotate filtered somatic VCFs before TMB counting.

Ported from ``genome-stability/index/burden/vep/vep_ann.py`` and reworked
for repository house style (parameterized CLI, English documentation,
shared LogUtil logging).

Usage:
    python vep_ann.py -i in.vcf -o out.vep.vcf -c ~/.vep \\
        -s mus_musculus -y GRCm39
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

_BIN_DIR = Path(__file__).resolve().parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))
from LogUtil import setup_logger  # noqa: E402

logger = setup_logger("VepAnn", stream=sys.stdout)


class VEP_Annotation:
    """Cache-based Ensembl VEP annotation of VCF files.

    Parameters
    ----------
    vep_cache_dir : str
        Root directory of the VEP cache (``~`` is expanded).
    species : str
        VEP species name, e.g. ``mus_musculus``.
    assembly : str
        Assembly name as known to the cache, e.g. ``GRCm39``.
    """

    def __init__(self, vep_cache_dir: str = "~/.vep",
                 species: str = "mus_musculus",
                 assembly: str = "GRCm39"):
        self.vep_cache_dir = str(Path(vep_cache_dir).expanduser().resolve())
        self.species = species
        self.assembly = assembly
        os.makedirs(self.vep_cache_dir, exist_ok=True)
        logger.debug(f"VEP cache root: {self.vep_cache_dir}")

    def _run_cmd(self, cmd):
        """Run an external command; log output, raise on failure.

        Raises RuntimeError with the captured stdout/stderr when the
        command exits non-zero, and a clear message when the binary is
        missing from $PATH.

        Args:
            cmd: argv list (cmd[0] is the executable).

        Returns:
            str: the command's stdout.
        """
        logger.info("Running: " + " ".join(cmd))
        if shutil.which(cmd[0]) is None:
            raise RuntimeError(
                f"Command not found: {cmd[0]!r}; install it and put it on $PATH")
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            logger.error(f"Command failed with exit code {exc.returncode}")
            logger.error(f"STDOUT:\n{exc.stdout or '[empty]'}")
            logger.error(f"STDERR:\n{exc.stderr or '[empty]'}")
            raise RuntimeError(f"Command execution failed: {' '.join(cmd)}") from exc
        if result.stdout:
            logger.debug(f"Command output:\n{result.stdout}")
        return result.stdout

    def vep_annotation_install(self):
        """Download and install the configured species/assembly VEP cache."""
        self._run_cmd([
            "vep_install", "-a", "cf",
            "-s", self.species,
            "-y", self.assembly,
            "-c", self.vep_cache_dir,
        ])

    def annotate_sv_vep(self, in_vcf: str, outfile: str,
                        result_format: str = "vcf"):
        """Annotate a VCF with VEP (offline cache; everything + pick per gene).

        Installs the species cache first when it is missing under the cache
        root. Output format follows ``result_format`` (``vcf`` or ``tab``).

        Args:
            in_vcf: input VCF path (plain or gzipped).
            outfile: output path for the annotated result.
            result_format: VEP output format flag (default: vcf).

        Returns:
            str: the output path.
        """
        species_cache = os.path.join(self.vep_cache_dir, self.species)
        if not os.path.exists(species_cache):
            logger.warning(
                f"Cache for {self.species} not found under "
                f"{self.vep_cache_dir}; attempting install")
            self.vep_annotation_install()
        self._run_cmd([
            "vep", "-i", in_vcf, "-o", outfile,
            "--cache", "--dir_cache", self.vep_cache_dir,
            "--species", self.species,
            "--assembly", self.assembly,
            "--format", "vcf", f"--{result_format}", "--force_overwrite",
            "--everything", "--pick", "--per_gene", "--offline",
        ])
        logger.info(f"VEP annotation finished: {outfile}")
        return outfile


def main(argv=None):
    """CLI entry point: annotate one VCF and write the annotated output."""
    parser = argparse.ArgumentParser(
        description="Cache-based Ensembl VEP annotation of a VCF.")
    parser.add_argument("-i", "--vcf", required=True, help="input VCF path")
    parser.add_argument("-o", "--out", required=True,
                        help="annotated output path")
    parser.add_argument("-c", "--vep-cache", default="~/.vep",
                        help="VEP cache root directory (default: ~/.vep)")
    parser.add_argument("-s", "--species", default="mus_musculus",
                        help="VEP species name (default: mus_musculus)")
    parser.add_argument("-y", "--assembly", default="GRCm39",
                        help="assembly name (default: GRCm39)")
    parser.add_argument("-l", "--log", default=None,
                        help="optional log file path (default: stdout only)")
    args = parser.parse_args(argv)

    global logger
    if args.log:
        log_file = Path(args.log)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        logger = setup_logger("VepAnn", log_file=str(log_file), stream=sys.stdout)
    analysis = VEP_Annotation(
        vep_cache_dir=args.vep_cache,
        species=args.species,
        assembly=args.assembly,
    )
    analysis.annotate_sv_vep(args.vcf, args.out)


if __name__ == "__main__":
    main()
