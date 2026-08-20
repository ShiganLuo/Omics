"""Environment utility: generate Apptainer definition files from conda env
YAMLs, build SIF images, and (optionally) the legacy Docker build/save flow.

Each Omics module directory under ``modules/`` contains one or more conda env
YAML files.  This utility:

1. Generates an Apptainer ``.def`` file next to each conda YAML (in the module
   directory), turning ``<name>.yaml`` into ``<name>.def``.
2. Runs ``apptainer build`` to produce a ``.sif`` image.
3. Copies the ``.def`` and ``.yaml`` to the output directory, preserving the
   module-level hierarchy.

The legacy Docker flow (``build_image`` / ``save_image`` / ``process_dockerfile``)
remains for backward compatibility but the primary workflow is now Apptainer.
"""

import re
import os
import shutil
import shlex
import logging
import argparse
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import urlparse
try:
    from LogUtil import setup_logger
except ImportError:
    from .LogUtil import setup_logger

logger = setup_logger(__name__, level=logging.INFO)

# ---------------------------------------------------------------------------
# Regex patterns for parsing Dockerfile comments
# ---------------------------------------------------------------------------
_BUILD_RE = re.compile(r"#\s*Build:\s*docker\s+build\s+-t\s+(\S+)")
_SOURCE_RE = re.compile(r"#\s*Source:\s*(.+)")
_SAVE_RE = re.compile(r"#\s*Save:\s*docker\s+save\s+\S+\s+-o\s+(.+)")
_COPY_RE = re.compile(r"^COPY\s+(\S+)\s+", re.IGNORECASE)

# --------------------------------------------------------------------------- #
# Regex patterns for parsing Apptainer .def comments
# --------------------------------------------------------------------------- #
_DEF_BUILD_RE = re.compile(r"#\s*Build:\s*apptainer\s+build\s+(\S+)\.sif")
_DEF_SOURCE_RE = re.compile(r"#\s*Source YAML:\s*(.+)")


class EnvUtil:
    """Build Apptainer SIF / Docker images from conda env YAMLs.

    Supports two backends with parallel workflows:

    - **Apptainer**: YAML -> ``.def`` -> ``.sif``
    - **Docker**:    YAML -> ``.dockerfile`` -> image -> ``.tar``

    Each backend supports three action levels:

    1. ``gen``:   generate build file (.def / .dockerfile) from YAML
    2. ``build``: build image (.sif / .tar) from existing build file
    3. ``all``:   full pipeline from YAML to image (gen + build + copy)

    Parameters
    ----------
    modules_dir : str
        Root directory containing module sub-directories (e.g. ``.../modules``).
    output_dir : str
        Base output directory for images and copied build files / YAMLs.
        The module-level hierarchy below ``modules_dir`` is preserved.
    """

    def __init__(
        self,
        modules_dir: str,
        output_dir: str,
    ) -> None:
        self.modules_dir = Path(modules_dir).resolve()
        self.output_dir = Path(output_dir).resolve()

        if not self.modules_dir.is_dir():
            raise FileNotFoundError(f"modules_dir not found: {self.modules_dir}")

        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------------
    @staticmethod
    def parse_dockerfile(dockerfile_path: str) -> Dict[str, Optional[str]]:
        """Extract image name, conda YAML path, and save path from a Dockerfile.

        Parameters
        ----------
        dockerfile_path : str
            Path to the Dockerfile.

        Returns
        -------
        Dict[str, Optional[str]]
            Keys: ``image_name``, ``yaml_filename``, ``source_path``,
            ``save_path``.
        """
        dockerfile = Path(dockerfile_path)
        if not dockerfile.is_file():
            raise FileNotFoundError(f"Dockerfile not found: {dockerfile}")

        text = dockerfile.read_text(encoding="utf-8")

        image_name: Optional[str] = None
        source_path: Optional[str] = None
        save_path: Optional[str] = None
        yaml_filename: Optional[str] = None

        for line in text.splitlines():
            m = _BUILD_RE.search(line)
            if m and image_name is None:
                image_name = m.group(1)

            m = _SOURCE_RE.search(line)
            if m and source_path is None:
                source_path = m.group(1).strip()

            m = _SAVE_RE.search(line)
            if m and save_path is None:
                save_path = m.group(1).strip()

            m = _COPY_RE.match(line)
            if m and yaml_filename is None:
                yaml_filename = m.group(1)

        return {
            "image_name": image_name,
            "yaml_filename": yaml_filename,
            "source_path": source_path,
            "save_path": save_path,
        }

    def resolve_image_name(
        self,
        dockerfile_path: str,
        parsed: Optional[Dict[str, Optional[str]]] = None,
    ) -> str:
        """Determine the Docker image name.

        Priority:
        1. ``# Build: docker build -t <name>`` comment in Dockerfile.
        2. ``name:`` field in the conda YAML referenced by COPY.
        3. Module directory name (last component of the parent path).

        Parameters
        ----------
        dockerfile_path : str
            Path to the Dockerfile.
        parsed : Optional[Dict], optional
            Pre-parsed Dockerfile info from :meth:`parse_dockerfile`.

        Returns
        -------
        str
            The resolved image name.
        """
        if parsed is None:
            parsed = self.parse_dockerfile(dockerfile_path)

        # 1. From Dockerfile comment
        if parsed["image_name"]:
            return parsed["image_name"]

        dockerfile = Path(dockerfile_path)
        module_dir = dockerfile.parent

        # 2. From conda YAML name: field
        yaml_name = parsed.get("yaml_filename")
        if yaml_name:
            yaml_path = module_dir / yaml_name
            if yaml_path.is_file():
                name = self._read_yaml_name(yaml_path)
                if name:
                    return name

        # 3. Fallback: directory name
        return module_dir.name

    @staticmethod
    def _read_yaml_name(yaml_path: Path) -> Optional[str]:
        """Read the ``name:`` field from a conda env YAML without full YAML parse."""
        try:
            for line in yaml_path.read_text(encoding="utf-8").splitlines():
                m = re.match(r"^name:\s*(\S+)", line)
                if m:
                    return m.group(1)
        except OSError:
            pass
        return None

    # ------------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------------
    def get_module_relpath(self, dockerfile_path: str) -> Path:
        """Return the module path relative to ``modules_dir``.

        Example::

            modules_dir = .../modules
            dockerfile  = .../modules/openms/decoydatabase/Dockerfile
            -> "openms/decoydatabase"

        Parameters
        ----------
        dockerfile_path : str
            Path to the Dockerfile.

        Returns
        -------
        Path
            Relative path from ``modules_dir`` to the Dockerfile's parent.
        """
        dockerfile = Path(dockerfile_path).resolve()
        try:
            rel = dockerfile.parent.relative_to(self.modules_dir)
        except ValueError:
            # Dockerfile is not under modules_dir — use its parent name
            logger.warning(
                f"Dockerfile {dockerfile} is not under modules_dir {self.modules_dir}, "
                f"using parent directory name only"
            )
            rel = Path(dockerfile.parent.name)
        return rel

    def get_output_subdir(self, dockerfile_path: str) -> Path:
        """Return the output sub-directory for a given Dockerfile.

        Preserves the module-level hierarchy.
        """
        rel = self.get_module_relpath(dockerfile_path)
        return self.output_dir / rel

    # ------------------------------------------------------------------
    # Docker operations
    # ------------------------------------------------------------------
    @staticmethod
    def _run_streaming(
        cmd: List[str],
        cwd: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> int:
        """Run a command and stream merged stdout/stderr line by line.

        Parameters
        ----------
        cmd : List[str]
            Command tokens.
        cwd : Optional[str]
            Working directory for the subprocess.
        timeout : Optional[int]
            Maximum seconds to wait.

        Returns
        -------
        int
            Return code.

        Raises
        ------
        RuntimeError
            If the command exits with a non-zero code.
        """
        cmd_str = " ".join(shlex.quote(c) for c in cmd)
        logger.info(f"Running: {cmd_str}")

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=cwd,
        )
        assert process.stdout is not None
        for line in process.stdout:
            logger.info(line.rstrip())
        process.wait(timeout=timeout)

        if process.returncode != 0:
            logger.error(f"Command failed (rc={process.returncode}): {cmd_str}")
            raise RuntimeError(f"Command failed with code {process.returncode}: {cmd_str}")

        return process.returncode

    # ------------------------------------------------------------------
    # YAML parsing helpers
    # ------------------------------------------------------------------
    @staticmethod
    def find_env_yamls(module_dir: str) -> List[Path]:
        """Find conda env YAML files in a module directory.

        Excludes ``*.schema.yaml`` files which are Snakemake schema definitions,
        not conda environment files.

        Parameters
        ----------
        module_dir : str
            Path to the module directory.

        Returns
        -------
        List[Path]
            Sorted list of conda env YAML paths.
        """
        d = Path(module_dir)
        if not d.is_dir():
            return []
        result: List[Path] = []
        for f in sorted(d.iterdir()):
            if not f.is_file():
                continue
            if f.suffix not in (".yaml", ".yml"):
                continue
            if f.name.endswith(".schema.yaml") or f.name.endswith(".schema.yml"):
                continue
            result.append(f)
        return result

    @staticmethod
    def parse_conda_yaml(yaml_path: str) -> Dict[str, object]:
        """Parse a conda env YAML to extract name, channels, and dependencies.

        Uses a lightweight line-based parser to avoid a hard dependency on
        PyYAML/ruamel.  Only the top-level ``name``, ``channels``, and
        ``dependencies`` keys are extracted.

        Parameters
        ----------
        yaml_path : str
            Path to the conda env YAML.

        Returns
        -------
        Dict[str, object]
            Keys: ``name`` (Optional[str]), ``channels`` (List[str]),
            ``dependencies`` (List[str]).
        """
        p = Path(yaml_path)
        text = p.read_text(encoding="utf-8")

        name: Optional[str] = None
        channels: List[str] = []
        deps: List[str] = []

        current_section: Optional[str] = None
        for raw_line in text.splitlines():
            stripped = raw_line.rstrip()
            if not stripped or stripped.startswith("#"):
                continue

            # Top-level key (no leading whitespace)
            if not raw_line[0].isspace():
                m = re.match(r"^(\w+):\s*(.*)", stripped)
                if not m:
                    current_section = None
                    continue
                key, val = m.group(1), m.group(2)
                if key == "name" and val:
                    name = val.strip()
                    current_section = None
                elif key in ("channels", "dependencies") and val:
                    # Inline value, e.g. "dependencies: [a, b]"
                    current_section = key
                    items = [x.strip() for x in val.strip("[]").split(",") if x.strip()]
                    if key == "channels":
                        channels.extend(items)
                    else:
                        deps.extend(items)
                elif key in ("channels", "dependencies"):
                    current_section = key
                else:
                    current_section = None
            else:
                # Indented list item
                m = re.match(r"^\s+-\s+(.+)", stripped)
                if m and current_section:
                    item = m.group(1).strip()
                    if current_section == "channels":
                        channels.append(item)
                    elif current_section == "dependencies":
                        # Skip nested dicts (pip: section) -- handled separately
                        if ":" in item and not item.startswith("'") and not item.startswith('"'):
                            continue
                        deps.append(item)

        return {
            "name": name,
            "channels": channels,
            "dependencies": deps,
        }

    @staticmethod
    def resolve_env_name(yaml_path: str, fallback: Optional[str] = None) -> str:
        """Determine the conda env / image name from a YAML file.

        Priority:
        1. ``name:`` field in the YAML.
        2. Filename stem (e.g. ``bedtools.yaml`` -> ``bedtools``).
        3. ``fallback`` argument.

        Parameters
        ----------
        yaml_path : str
            Path to the conda env YAML.
        fallback : Optional[str]
            Fallback name if YAML has no ``name:`` and stem is generic.

        Returns
        -------
        str
            The resolved env name.
        """
        parsed = EnvUtil.parse_conda_yaml(yaml_path)
        if parsed["name"]:
            return str(parsed["name"])
        stem = Path(yaml_path).stem
        return stem if stem else (fallback or "env")

    # ------------------------------------------------------------------
    # Apptainer definition file generation
    # ------------------------------------------------------------------
    APPTAINER_BOOTSTRAP = "docker"
    APPTAINER_FROM = "docker.m.daocloud.io/continuumio/miniconda3:latest"

    @classmethod
    def generate_def_file(
        cls,
        yaml_path: str,
        def_path: Optional[str] = None,
        env_name: Optional[str] = None,
    ) -> str:
        """Generate an Apptainer definition file from a conda env YAML.

        The ``.def`` file mirrors the existing Dockerfile pattern:
        ``FROM continuumio/miniconda3:latest``, copy the YAML, run
        ``conda env create``, set PATH.

        Parameters
        ----------
        yaml_path : str
            Path to the source conda env YAML.
        def_path : Optional[str]
            Output ``.def`` file path. If ``None``, placed next to the YAML
            as ``<parent_dir>.def`` (matching the dockerfile convention).
        env_name : Optional[str]
            Override the conda env name. If ``None``, resolved from the YAML.

        Returns
        -------
        str
            Path to the generated ``.def`` file.
        """
        yaml = Path(yaml_path).resolve()
        if not yaml.is_file():
            raise FileNotFoundError(f"Conda YAML not found: {yaml}")

        if env_name is None:
            env_name = cls.resolve_env_name(str(yaml))

        if def_path is None:
            # Default: "<parent_dir>.def" in the YAML's directory.
            # If a .def already exists but references a *different* YAML,
            # use "<env_name>.def" to avoid clobbering it (multi-env modules).
            default_path = yaml.parent / f"{yaml.parent.name}.def"
            if default_path.is_file():
                existing_text = default_path.read_text(encoding="utf-8")
                if yaml.name not in existing_text:
                    def_path = str(yaml.parent / f"{env_name}.def")
                    logger.info(
                        f"Existing .def references a different YAML; "
                        f"using {env_name}.def to avoid collision"
                    )
                else:
                    def_path = str(default_path)
            else:
                def_path = str(default_path)

        yaml_content = yaml.read_text(encoding="utf-8")

        def_text = cls._render_def(
            env_name=env_name,
            yaml_filename=yaml.name,
            yaml_content=yaml_content,
            def_filename=Path(def_path).name,
        )

        out = Path(def_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(def_text, encoding="utf-8")
        logger.info(f"Generated Apptainer def file: {out}")

        return str(out)

    @classmethod
    def _render_def(
        cls,
        env_name: str,
        yaml_filename: str,
        yaml_content: str,
        def_filename: str = "env.def",
    ) -> str:
        """Render the Apptainer def file text from a template.

        Parameters
        ----------
        env_name : str
            Conda environment name (used for PATH).
        yaml_filename : str
            Filename of the conda YAML (for COPY into the container).
        yaml_content : str
            Full text of the conda YAML (embedded as comment for reference).
        def_filename : str
            Filename of the .def file (for the ``# Build:`` comment).

        Returns
        -------
        str
            The rendered ``.def`` file content.
        """
        return f"""# Auto-generated Apptainer definition file for {env_name} conda environment
# Source YAML: {yaml_filename}
# Build: apptainer build {env_name}.sif {def_filename}
# Or:   apptainer build <output.sif> <def_file>

Bootstrap: {cls.APPTAINER_BOOTSTRAP}
From: {cls.APPTAINER_FROM}

%files
    {yaml_filename} /opt/conda/{yaml_filename}

%post
    export PATH="/opt/conda/bin:$PATH"
    export CONDA_NO_PLUGINS=true
    conda config --set solver classic
    conda env create -f /opt/conda/{yaml_filename} && \\
        conda clean -afy && \\
        rm /opt/conda/{yaml_filename}

    echo 'export PATH="/opt/conda/envs/{env_name}/bin:$PATH"' >> /etc/profile.d/conda_{env_name}.sh
    echo 'export CONDA_DEFAULT_ENV={env_name}' >> /etc/profile.d/conda_{env_name}.sh

%environment
    export PATH="/opt/conda/envs/{env_name}/bin:$PATH"
    export CONDA_DEFAULT_ENV={env_name}

%runscript
    exec "$@"

%labels
    Author genomestability
    EnvName {env_name}

%help
    Apptainer image for conda environment '{env_name}'.
    Generated from {yaml_filename}.
"""

    # ------------------------------------------------------------------
    # Dockerfile generation (from conda YAML, same pattern as .def)
    # ------------------------------------------------------------------
    DOCKER_FROM = "docker.m.daocloud.io/continuumio/miniconda3:latest"

    @classmethod
    def generate_dockerfile(
        cls,
        yaml_path: str,
        dockerfile_path: Optional[str] = None,
        image_name: Optional[str] = None,
    ) -> str:
        """Generate a Dockerfile from a conda env YAML.

        Mirrors the existing Dockerfile pattern found across modules:
        ``FROM continuumio/miniconda3:latest``, ``COPY <yaml>``,
        ``conda env create``, set ``ENV PATH``.

        Parameters
        ----------
        yaml_path : str
            Path to the source conda env YAML.
        dockerfile_path : Optional[str]
            Output Dockerfile path. If ``None``, placed next to the YAML
            as ``Dockerfile`` (matching the existing convention).
        image_name : Optional[str]
            Override the image name used in the ``# Build:`` comment.
            If ``None``, resolved from the YAML.

        Returns
        -------
        str
            Path to the generated Dockerfile.
        """
        yaml = Path(yaml_path).resolve()
        if not yaml.is_file():
            raise FileNotFoundError(f"Conda YAML not found: {yaml}")

        if image_name is None:
            image_name = cls.resolve_env_name(str(yaml))

        if dockerfile_path is None:
            # Default: "<parent_dir>.dockerfile" in the YAML's directory.
            # If a dockerfile already exists but references a *different* YAML,
            # use "<env_name>.dockerfile" to avoid clobbering it (multi-env modules).
            default_path = yaml.parent / f"{yaml.parent.name}.dockerfile"
            if default_path.is_file():
                existing_text = default_path.read_text(encoding="utf-8")
                if yaml.name not in existing_text:
                    dockerfile_path = str(yaml.parent / f"{image_name}.dockerfile")
                    logger.info(
                        f"Existing dockerfile references a different YAML; "
                        f"using {image_name}.dockerfile to avoid collision"
                    )
                else:
                    dockerfile_path = str(default_path)
            else:
                dockerfile_path = str(default_path)

        df_text = cls._render_dockerfile(
            image_name=image_name,
            yaml_filename=yaml.name,
            yaml_path=str(yaml),
            dockerfile_filename=Path(dockerfile_path).name,
        )

        out = Path(dockerfile_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(df_text, encoding="utf-8")
        logger.info(f"Generated Dockerfile: {out}")

        return str(out)

    @classmethod
    def _render_dockerfile(
        cls,
        image_name: str,
        yaml_filename: str,
        yaml_path: str,
        dockerfile_filename: str = "Dockerfile",
    ) -> str:
        """Render the Dockerfile text from a template.

        Parameters
        ----------
        image_name : str
            Docker image name (tag).
        yaml_filename : str
            Filename of the conda YAML (for COPY).
        yaml_path : str
            Absolute path to the YAML (for the ``# Source:`` comment).
        dockerfile_filename : str
            Filename of the dockerfile (for the ``# Build:`` comment).

        Returns
        -------
        str
            The rendered Dockerfile content.
        """
        return f"""# Auto-generated Dockerfile for {image_name} conda environment
# Source: {yaml_path}
# Build: docker build -t {image_name} -f {dockerfile_filename} .
# Save:  docker save {image_name} -o /home/luosg/Database/env/{image_name}.tar

FROM {cls.DOCKER_FROM}

COPY {yaml_filename} /tmp/{yaml_filename}

RUN conda env create -f /tmp/{yaml_filename} && \\
    conda clean -afy && \\
    rm /tmp/{yaml_filename}

ENV CONDA_DEFAULT_ENV={image_name}
ENV PATH="/opt/conda/envs/{image_name}/bin:$PATH"

CMD ["bash"]
"""

    # ------------------------------------------------------------------
    # Apptainer build
    # ------------------------------------------------------------------
    def build_sif(
        self,
        def_path: str,
        output_sif: str,
        fakeroot: bool = True,
        force: bool = True,
        timeout: Optional[int] = None,
    ) -> str:
        """Run ``apptainer build`` to create a SIF image from a def file.

        Parameters
        ----------
        def_path : str
            Path to the Apptainer ``.def`` file.
        output_sif : str
            Output ``.sif`` file path.
        fakeroot : bool
            Use ``--fakeroot`` for unprivileged builds (default: True).
        force : bool
            Use ``--force`` to overwrite existing SIF (default: True).
        timeout : Optional[int]
            Build timeout in seconds.

        Returns
        -------
        str
            Path to the built ``.sif`` file.
        """
        df = Path(def_path).resolve()
        if not df.is_file():
            raise FileNotFoundError(f"Def file not found: {df}")

        out = Path(output_sif).resolve()
        out.parent.mkdir(parents=True, exist_ok=True)

        cmd: List[str] = ["apptainer", "build"]
        if fakeroot:
            cmd.append("--fakeroot")
        if force:
            cmd.append("--force")
        cmd.extend([str(out), str(df)])

        logger.info(f"Building SIF: {out}")
        self._run_streaming(cmd, cwd=str(df.parent), timeout=timeout)
        logger.info(f"SIF built: {out}")

        return str(out)

    # ------------------------------------------------------------------
    # Full Apptainer workflow
    # ------------------------------------------------------------------
    def process_yaml(
        self,
        yaml_path: str,
        env_name: Optional[str] = None,
        output_subdir: Optional[str] = None,
        skip_existing: bool = True,
        fakeroot: bool = True,
        force: bool = True,
        build_timeout: Optional[int] = None,
        generate_only: bool = False,
    ) -> Dict[str, str]:
        """Full Apptainer pipeline for a single conda YAML: generate def, build SIF, copy.

        Steps:
        1. Generate ``<parent_dir>.def`` next to the source YAML (in the module dir).
        2. ``apptainer build <env_name>.sif <parent_dir>.def``.
        3. Copy ``.def`` and ``.yaml`` to the output sub-directory.

        Parameters
        ----------
        yaml_path : str
            Path to the conda env YAML.
        env_name : Optional[str]
            Override env / image name.
        output_subdir : Optional[str]
            Override output sub-directory. If ``None``, derived from the
            YAML's module path relative to ``modules_dir``.
        skip_existing : bool
            Skip if the SIF already exists in the output directory.
        fakeroot : bool
            Pass ``--fakeroot`` to ``apptainer build``.
        force : bool
            Pass ``--force`` to ``apptainer build``.
        build_timeout : Optional[int]
            Build timeout in seconds.
        generate_only : bool
            If ``True``, only generate the ``.def`` file without building.

        Returns
        -------
        Dict[str, str]
            Keys: ``env_name``, ``def_path``, ``sif_path``, ``output_dir``,
            ``yaml_dest``, ``def_dest``, ``module_relpath``.
        """
        yaml = Path(yaml_path).resolve()
        if not yaml.is_file():
            raise FileNotFoundError(f"YAML not found: {yaml}")

        if env_name is None:
            env_name = self.resolve_env_name(str(yaml))

        # Generate .def next to the YAML (in the module directory)
        def_path = self.generate_def_file(
            yaml_path=str(yaml),
            env_name=env_name,
        )

        # Determine output sub-directory
        if output_subdir is None:
            output_subdir = str(self.get_output_subdir_for_module(yaml.parent))
        out_dir = Path(output_subdir).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

        sif_path = out_dir / f"{env_name}.sif"

        result: Dict[str, str] = {
            "env_name": env_name,
            "def_path": def_path,
            "sif_path": str(sif_path),
            "output_dir": str(out_dir),
            "yaml_dest": str(out_dir / yaml.name),
            "def_dest": str(out_dir / Path(def_path).name),
            "module_relpath": str(self.get_module_relpath_for_module(yaml.parent)),
        }

        # Copy .def and .yaml to output directory
        shutil.copy2(yaml, out_dir / yaml.name)
        shutil.copy2(def_path, out_dir / Path(def_path).name)
        logger.info(f"Copied {yaml.name} and {Path(def_path).name} to {out_dir}")

        if generate_only:
            logger.info(f"[GENERATE_ONLY] Def file created: {def_path}")
            return result

        if skip_existing and sif_path.exists():
            logger.info(f"[SKIP] SIF already exists: {sif_path}")
            return result

        # Build SIF
        self.build_sif(
            def_path=def_path,
            output_sif=str(sif_path),
            fakeroot=fakeroot,
            force=force,
            timeout=build_timeout,
        )

        logger.info(
            f"[DONE] {yaml.name} -> env={env_name}, "
            f"sif={sif_path}, dir={out_dir}"
        )
        return result

    # ------------------------------------------------------------------
    # Full Docker workflow (from YAML: generate Dockerfile, build, save)
    # ------------------------------------------------------------------
    def process_yaml_docker(
        self,
        yaml_path: str,
        image_name: Optional[str] = None,
        output_subdir: Optional[str] = None,
        skip_existing: bool = True,
        no_cache: bool = False,
        build_timeout: Optional[int] = None,
        save_timeout: Optional[int] = None,
        generate_only: bool = False,
    ) -> Dict[str, str]:
        """Full Docker pipeline for a single conda YAML: generate Dockerfile, build, save.

        Steps:
        1. Generate ``Dockerfile`` next to the source YAML (in the module dir).
        2. ``docker build -t <image> -f Dockerfile .``
        3. ``docker save <image> -o <image>.tar``
        4. Copy ``Dockerfile`` and ``.yaml`` to the output sub-directory.

        Parameters
        ----------
        yaml_path : str
            Path to the conda env YAML.
        image_name : Optional[str]
            Override Docker image name. If ``None``, resolved from the YAML.
        output_subdir : Optional[str]
            Override output sub-directory. If ``None``, derived from the
            YAML's module path relative to ``modules_dir``.
        skip_existing : bool
            Skip if the tar archive already exists in the output directory.
        no_cache : bool
            Pass ``--no-cache`` to ``docker build``.
        build_timeout : Optional[int]
        save_timeout : Optional[int]
        generate_only : bool
            If ``True``, only generate the Dockerfile without building.

        Returns
        -------
        Dict[str, str]
            Keys: ``image_name``, ``dockerfile_path``, ``tar_path``,
            ``output_dir``, ``yaml_dest``, ``dockerfile_dest``,
            ``module_relpath``.
        """
        yaml = Path(yaml_path).resolve()
        if not yaml.is_file():
            raise FileNotFoundError(f"YAML not found: {yaml}")

        if image_name is None:
            image_name = self.resolve_env_name(str(yaml))

        # Generate Dockerfile next to the YAML (in the module directory)
        dockerfile_path = self.generate_dockerfile(
            yaml_path=str(yaml),
            image_name=image_name,
        )

        # Determine output sub-directory
        if output_subdir is None:
            output_subdir = str(self.get_output_subdir_for_module(yaml.parent))
        out_dir = Path(output_subdir).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

        tar_path = out_dir / f"{image_name}.tar"

        result: Dict[str, str] = {
            "image_name": image_name,
            "dockerfile_path": dockerfile_path,
            "tar_path": str(tar_path),
            "output_dir": str(out_dir),
            "yaml_dest": str(out_dir / yaml.name),
            "dockerfile_dest": str(out_dir / Path(dockerfile_path).name),
            "module_relpath": str(self.get_module_relpath_for_module(yaml.parent)),
        }

        # Copy Dockerfile and YAML to output directory
        shutil.copy2(yaml, out_dir / yaml.name)
        shutil.copy2(dockerfile_path, out_dir / Path(dockerfile_path).name)
        logger.info(f"Copied {yaml.name} and {Path(dockerfile_path).name} to {out_dir}")

        if generate_only:
            logger.info(f"[GENERATE_ONLY] Dockerfile created: {dockerfile_path}")
            return result

        if skip_existing and tar_path.exists():
            logger.info(f"[SKIP] Tar already exists: {tar_path}")
            return result

        # Build Docker image
        self.build_image(
            dockerfile_path=dockerfile_path,
            image_name=image_name,
            no_cache=no_cache,
            timeout=build_timeout,
        )

        # Save as tar
        self.save_image(
            image_name=image_name,
            output_tar=str(tar_path),
            timeout=save_timeout,
        )

        logger.info(
            f"[DONE] {yaml.name} -> image={image_name}, "
            f"tar={tar_path}, dir={out_dir}"
        )
        return result

    def process_module(
        self,
        module_dir: str,
        skip_existing: bool = True,
        fakeroot: bool = True,
        force: bool = True,
        build_timeout: Optional[int] = None,
        generate_only: bool = False,
    ) -> List[Dict[str, str]]:
        """Process all conda env YAMLs in a single module directory.

        Parameters
        ----------
        module_dir : str
            Path to the module directory.
        skip_existing : bool
        fakeroot : bool
        force : bool
        build_timeout : Optional[int]
        generate_only : bool

        Returns
        -------
        List[Dict[str, str]]
            Results from :meth:`process_yaml` for each YAML found.
        """
        yamls = self.find_env_yamls(module_dir)
        if not yamls:
            logger.warning(f"No conda env YAML found in: {module_dir}")
            return []

        results: List[Dict[str, str]] = []
        for yaml in yamls:
            logger.info(f"Processing YAML: {yaml}")
            try:
                result = self.process_yaml(
                    yaml_path=str(yaml),
                    skip_existing=skip_existing,
                    fakeroot=fakeroot,
                    force=force,
                    build_timeout=build_timeout,
                    generate_only=generate_only,
                )
                results.append(result)
            except Exception as e:
                logger.error(f"[FAIL] {yaml}: {e}")
        return results

    def process_all_apptainer(
        self,
        skip_existing: bool = True,
        fakeroot: bool = True,
        force: bool = True,
        build_timeout: Optional[int] = None,
        generate_only: bool = False,
    ) -> List[Dict[str, str]]:
        """Process every module under ``modules_dir``: generate defs, build SIFs.

        Parameters
        ----------
        skip_existing : bool
        fakeroot : bool
        force : bool
        build_timeout : Optional[int]
        generate_only : bool

        Returns
        -------
        List[Dict[str, str]]
            Results from :meth:`process_yaml` for each conda YAML found.
        """
        yamls = self.find_all_env_yamls()
        logger.info(f"Found {len(yamls)} conda env YAML(s) under {self.modules_dir}")

        results: List[Dict[str, str]] = []
        failed: List[Tuple[str, str]] = []

        for i, yaml in enumerate(yamls, 1):
            rel = self.get_module_relpath_for_module(yaml.parent)
            logger.info(f"[{i}/{len(yamls)}] Processing: {rel}/{yaml.name}")
            try:
                result = self.process_yaml(
                    yaml_path=str(yaml),
                    skip_existing=skip_existing,
                    fakeroot=fakeroot,
                    force=force,
                    build_timeout=build_timeout,
                    generate_only=generate_only,
                )
                results.append(result)
            except Exception as e:
                logger.error(f"[FAIL] {rel}/{yaml.name}: {e}")
                failed.append((str(yaml), str(e)))

        logger.info(
            f"Summary: {len(results)}/{len(yamls)} succeeded, "
            f"{len(failed)} failed"
        )
        for yaml_path, err in failed:
            logger.error(f"  FAILED: {yaml_path} -> {err}")

        return results

    def find_all_env_yamls(self) -> List[Path]:
        """Find all conda env YAML files under ``modules_dir``.

        Excludes ``*.schema.yaml`` files.

        Returns
        -------
        List[Path]
            Sorted list of YAML paths.
        """
        results: List[Path] = []
        for f in sorted(self.modules_dir.rglob("*.yaml")):
            if not f.is_file():
                continue
            if f.name.endswith(".schema.yaml"):
                continue
            results.append(f)
        return results

    def get_module_relpath_for_module(self, module_dir: Path) -> Path:
        """Return the module path relative to ``modules_dir``.

        Parameters
        ----------
        module_dir : Path
            A module directory under ``modules_dir``.

        Returns
        -------
        Path
            Relative path from ``modules_dir`` to ``module_dir``.
        """
        try:
            return Path(module_dir).resolve().relative_to(self.modules_dir)
        except ValueError:
            logger.warning(
                f"Module dir {module_dir} is not under modules_dir {self.modules_dir}, "
                f"using directory name only"
            )
            return Path(module_dir.name)

    def get_output_subdir_for_module(self, module_dir: Path) -> Path:
        """Return the output sub-directory for a module directory.

        Preserves the module-level hierarchy.
        """
        rel = self.get_module_relpath_for_module(module_dir)
        return self.output_dir / rel

    # ------------------------------------------------------------------
    # Batch Docker workflow (from YAML)
    # ------------------------------------------------------------------
    def process_module_docker(
        self,
        module_dir: str,
        skip_existing: bool = True,
        no_cache: bool = False,
        build_timeout: Optional[int] = None,
        save_timeout: Optional[int] = None,
        generate_only: bool = False,
    ) -> List[Dict[str, str]]:
        """Process all conda env YAMLs in a single module directory (Docker flow).

        Parameters
        ----------
        module_dir : str
            Path to the module directory.
        skip_existing : bool
        no_cache : bool
        build_timeout : Optional[int]
        save_timeout : Optional[int]
        generate_only : bool

        Returns
        -------
        List[Dict[str, str]]
            Results from :meth:`process_yaml_docker` for each YAML found.
        """
        yamls = self.find_env_yamls(module_dir)
        if not yamls:
            logger.warning(f"No conda env YAML found in: {module_dir}")
            return []

        results: List[Dict[str, str]] = []
        for yaml in yamls:
            logger.info(f"Processing YAML (docker): {yaml}")
            try:
                result = self.process_yaml_docker(
                    yaml_path=str(yaml),
                    skip_existing=skip_existing,
                    no_cache=no_cache,
                    build_timeout=build_timeout,
                    save_timeout=save_timeout,
                    generate_only=generate_only,
                )
                results.append(result)
            except Exception as e:
                logger.error(f"[FAIL] {yaml}: {e}")
        return results

    def process_all_docker(
        self,
        skip_existing: bool = True,
        no_cache: bool = False,
        build_timeout: Optional[int] = None,
        save_timeout: Optional[int] = None,
        generate_only: bool = False,
    ) -> List[Dict[str, str]]:
        """Process every module under ``modules_dir``: generate Dockerfiles, build, save.

        Parameters
        ----------
        skip_existing : bool
        no_cache : bool
        build_timeout : Optional[int]
        save_timeout : Optional[int]
        generate_only : bool

        Returns
        -------
        List[Dict[str, str]]
            Results from :meth:`process_yaml_docker` for each conda YAML found.
        """
        yamls = self.find_all_env_yamls()
        logger.info(f"Found {len(yamls)} conda env YAML(s) under {self.modules_dir}")

        results: List[Dict[str, str]] = []
        failed: List[Tuple[str, str]] = []

        for i, yaml in enumerate(yamls, 1):
            rel = self.get_module_relpath_for_module(yaml.parent)
            logger.info(f"[{i}/{len(yamls)}] Processing (docker): {rel}/{yaml.name}")
            try:
                result = self.process_yaml_docker(
                    yaml_path=str(yaml),
                    skip_existing=skip_existing,
                    no_cache=no_cache,
                    build_timeout=build_timeout,
                    save_timeout=save_timeout,
                    generate_only=generate_only,
                )
                results.append(result)
            except Exception as e:
                logger.error(f"[FAIL] {rel}/{yaml.name}: {e}")
                failed.append((str(yaml), str(e)))

        logger.info(
            f"Summary: {len(results)}/{len(yamls)} succeeded, "
            f"{len(failed)} failed"
        )
        for yaml_path, err in failed:
            logger.error(f"  FAILED: {yaml_path} -> {err}")

        return results

    def build_image(
        self,
        dockerfile_path: str,
        image_name: Optional[str] = None,
        no_cache: bool = False,
        timeout: Optional[int] = None,
    ) -> str:
        """Run ``docker build`` for a module Dockerfile.

        Parameters
        ----------
        dockerfile_path : str
            Path to the Dockerfile.
        image_name : Optional[str]
            Override image name. If ``None``, resolved from the Dockerfile.
        no_cache : bool
            Pass ``--no-cache`` to ``docker build``.
        timeout : Optional[int]
            Build timeout in seconds.

        Returns
        -------
        str
            The image name (tag) that was built.
        """
        dockerfile = Path(dockerfile_path).resolve()
        if not dockerfile.is_file():
            raise FileNotFoundError(f"Dockerfile not found: {dockerfile}")

        if image_name is None:
            image_name = self.resolve_image_name(dockerfile_path)

        module_dir = dockerfile.parent

        cmd: List[str] = [
            "docker", "build",
            "-t", image_name,
            "-f", str(dockerfile),
        ]
        if no_cache:
            cmd.append("--no-cache")
        cmd.append(".")

        logger.info(f"Building Docker image '{image_name}' from {dockerfile}")
        self._run_streaming(cmd, cwd=str(module_dir), timeout=timeout)
        logger.info(f"Successfully built image: {image_name}")

        return image_name

    def save_image(
        self,
        image_name: str,
        output_tar: str,
        timeout: Optional[int] = None,
    ) -> str:
        """Run ``docker save`` to export an image as a tar archive.

        Parameters
        ----------
        image_name : str
            Name/tag of the Docker image to save.
        output_tar : str
            Output ``.tar`` file path.
        timeout : Optional[int]
            Save timeout in seconds.

        Returns
        -------
        str
            The output tar path.
        """
        output_path = Path(output_tar).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)

        cmd: List[str] = [
            "docker", "save",
            image_name,
            "-o", str(output_path),
        ]

        logger.info(f"Saving image '{image_name}' to {output_path}")
        self._run_streaming(cmd, timeout=timeout)
        logger.info(f"Image saved to: {output_path}")

        return str(output_path)

    # ------------------------------------------------------------------
    # File copying
    # ------------------------------------------------------------------
    def copy_artefacts(
        self,
        dockerfile_path: str,
        output_subdir: str,
        parsed: Optional[Dict[str, Optional[str]]] = None,
    ) -> List[str]:
        """Copy the Dockerfile and its referenced conda YAML to the output directory.

        Parameters
        ----------
        dockerfile_path : str
            Path to the source Dockerfile.
        output_subdir : str
            Destination directory.
        parsed : Optional[Dict]
            Pre-parsed Dockerfile info.

        Returns
        -------
        List[str]
            List of copied file paths.
        """
        dockerfile = Path(dockerfile_path).resolve()
        dest_dir = Path(output_subdir).resolve()
        dest_dir.mkdir(parents=True, exist_ok=True)

        if parsed is None:
            parsed = self.parse_dockerfile(dockerfile_path)

        copied: List[str] = []

        # Copy Dockerfile
        dest_dockerfile = dest_dir / dockerfile.name
        shutil.copy2(dockerfile, dest_dockerfile)
        copied.append(str(dest_dockerfile))
        logger.info(f"Copied {dockerfile.name} -> {dest_dockerfile}")

        # Copy conda YAML if referenced by COPY directive
        yaml_name = parsed.get("yaml_filename")
        if yaml_name:
            yaml_src = dockerfile.parent / yaml_name
            if yaml_src.is_file():
                yaml_dest = dest_dir / yaml_name
                shutil.copy2(yaml_src, yaml_dest)
                copied.append(str(yaml_dest))
                logger.info(f"Copied conda YAML -> {yaml_dest}")
            else:
                logger.warning(f"Conda YAML not found: {yaml_src}")

        return copied

    # ------------------------------------------------------------------
    # Full workflow
    # ------------------------------------------------------------------
    def process_dockerfile(
        self,
        dockerfile_path: str,
        image_name: Optional[str] = None,
        no_cache: bool = False,
        build_timeout: Optional[int] = None,
        save_timeout: Optional[int] = None,
        skip_existing: bool = False,
    ) -> Dict[str, str]:
        """Full pipeline: build image, save tar, copy Dockerfile + YAML.

        Output layout::

            output_dir/<module_rel_path>/
                Dockerfile          (copied)
                <conda>.yaml        (copied if present)
                <image_name>.tar    (docker save)

        Parameters
        ----------
        dockerfile_path : str
            Path to the Dockerfile.
        image_name : Optional[str]
            Override image name.
        no_cache : bool
            Pass ``--no-cache`` to ``docker build``.
        build_timeout : Optional[int]
        save_timeout : Optional[int]
        skip_existing : bool
            If ``True``, skip when the tar archive already exists.

        Returns
        -------
        Dict[str, str]
            Keys: ``image_name``, ``tar_path``, ``output_dir``,
            ``dockerfile_dest``, ``module_relpath``.
        """
        dockerfile = Path(dockerfile_path).resolve()
        parsed = self.parse_dockerfile(dockerfile)

        if image_name is None:
            image_name = self.resolve_image_name(str(dockerfile), parsed)

        output_subdir = self.get_output_subdir(str(dockerfile))
        output_subdir.mkdir(parents=True, exist_ok=True)

        tar_path = output_subdir / f"{image_name}.tar"

        if skip_existing and tar_path.exists():
            logger.info(f"[SKIP] Tar already exists: {tar_path}")
            return {
                "image_name": image_name,
                "tar_path": str(tar_path),
                "output_dir": str(output_subdir),
                "dockerfile_dest": str(output_subdir / dockerfile.name),
                "module_relpath": str(self.get_module_relpath(str(dockerfile))),
            }

        # Build
        self.build_image(
            dockerfile_path=str(dockerfile),
            image_name=image_name,
            no_cache=no_cache,
            timeout=build_timeout,
        )

        # Save
        self.save_image(
            image_name=image_name,
            output_tar=str(tar_path),
            timeout=save_timeout,
        )

        # Copy Dockerfile + YAML
        copied = self.copy_artefacts(
            dockerfile_path=str(dockerfile),
            output_subdir=str(output_subdir),
            parsed=parsed,
        )
        dockerfile_dest = copied[0] if copied else ""

        logger.info(
            f"[DONE] {dockerfile.name} -> image={image_name}, "
            f"tar={tar_path}, dir={output_subdir}"
        )

        return {
            "image_name": image_name,
            "tar_path": str(tar_path),
            "output_dir": str(output_subdir),
            "dockerfile_dest": dockerfile_dest,
            "module_relpath": str(self.get_module_relpath(str(dockerfile))),
        }

    # ------------------------------------------------------------------
    # Apptainer: parse / find / process existing .def files
    # ------------------------------------------------------------------
    @staticmethod
    def parse_def_file(def_path: str) -> Dict[str, Optional[str]]:
        """Extract env name and YAML filename from an Apptainer ``.def`` file.

        Parameters
        ----------
        def_path : str
            Path to the ``.def`` file.

        Returns
        -------
        Dict[str, Optional[str]]
            Keys: ``env_name``, ``yaml_filename``.
        """
        df = Path(def_path)
        if not df.is_file():
            raise FileNotFoundError(f"Def file not found: {df}")

        text = df.read_text(encoding="utf-8")
        env_name: Optional[str] = None
        yaml_filename: Optional[str] = None

        for line in text.splitlines():
            m = _DEF_BUILD_RE.search(line)
            if m and env_name is None:
                env_name = m.group(1)
            m = _DEF_SOURCE_RE.search(line)
            if m and yaml_filename is None:
                yaml_filename = Path(m.group(1).strip()).name

        return {"env_name": env_name, "yaml_filename": yaml_filename}

    def resolve_env_name_from_def(
        self,
        def_path: str,
        parsed: Optional[Dict[str, Optional[str]]] = None,
    ) -> str:
        """Determine the env / SIF name from a ``.def`` file.

        Priority:
        1. ``# Build: apptainer build <name>.sif`` comment.
        2. Module directory name.

        Parameters
        ----------
        def_path : str
            Path to the ``.def`` file.
        parsed : Optional[Dict]
            Pre-parsed def info from :meth:`parse_def_file`.

        Returns
        -------
        str
            The resolved env name.
        """
        if parsed is None:
            parsed = self.parse_def_file(def_path)
        if parsed["env_name"]:
            return parsed["env_name"]
        return Path(def_path).parent.name

    def find_def_files(self) -> List[str]:
        """Find all Apptainer ``.def`` files under ``modules_dir``.

        Returns
        -------
        List[str]
            Sorted list of ``.def`` file paths.
        """
        results: List[str] = []
        for df in sorted(self.modules_dir.rglob("*.def")):
            if df.is_file():
                results.append(str(df))
        return results

    def process_def_file(
        self,
        def_path: str,
        env_name: Optional[str] = None,
        fakeroot: bool = True,
        force: bool = True,
        skip_existing: bool = True,
        build_timeout: Optional[int] = None,
    ) -> Dict[str, str]:
        """Build a SIF from an existing ``.def`` file and copy artefacts.

        Output layout::

            output_dir/<module_rel_path>/
                <parent>.def        (copied)
                <conda>.yaml        (copied if present)
                <env_name>.sif      (apptainer build)

        Parameters
        ----------
        def_path : str
            Path to the ``.def`` file.
        env_name : Optional[str]
            Override env / SIF name.
        fakeroot : bool
            Pass ``--fakeroot`` to ``apptainer build``.
        force : bool
            Pass ``--force`` to ``apptainer build``.
        skip_existing : bool
            Skip if the SIF already exists.
        build_timeout : Optional[int]
            Build timeout in seconds.

        Returns
        -------
        Dict[str, str]
            Keys: ``env_name``, ``sif_path``, ``output_dir``,
            ``def_dest``, ``module_relpath``.
        """
        df = Path(def_path).resolve()
        parsed = self.parse_def_file(str(df))

        if env_name is None:
            env_name = self.resolve_env_name_from_def(str(df), parsed)

        output_subdir = self.get_output_subdir(str(df))
        output_subdir.mkdir(parents=True, exist_ok=True)
        sif_path = output_subdir / f"{env_name}.sif"

        # Copy .def and YAML to output directory
        def_dest = output_subdir / df.name
        shutil.copy2(df, def_dest)
        yaml_name = parsed.get("yaml_filename")
        if yaml_name:
            yaml_src = df.parent / yaml_name
            if yaml_src.is_file():
                shutil.copy2(yaml_src, output_subdir / yaml_name)
                logger.info(f"Copied {df.name} and {yaml_name} to {output_subdir}")
            else:
                logger.warning(f"Conda YAML not found: {yaml_src}")
        else:
            logger.info(f"Copied {df.name} to {output_subdir}")

        if skip_existing and sif_path.exists():
            logger.info(f"[SKIP] SIF already exists: {sif_path}")
            return {
                "env_name": env_name,
                "sif_path": str(sif_path),
                "output_dir": str(output_subdir),
                "def_dest": str(def_dest),
                "module_relpath": str(self.get_module_relpath(str(df))),
            }

        # Build SIF
        self.build_sif(
            def_path=str(df),
            output_sif=str(sif_path),
            fakeroot=fakeroot,
            force=force,
            timeout=build_timeout,
        )

        logger.info(
            f"[DONE] {df.name} -> env={env_name}, "
            f"sif={sif_path}, dir={output_subdir}"
        )
        return {
            "env_name": env_name,
            "sif_path": str(sif_path),
            "output_dir": str(output_subdir),
            "def_dest": str(def_dest),
            "module_relpath": str(self.get_module_relpath(str(df))),
        }

    def find_dockerfiles(self) -> List[str]:
        """Find all Dockerfiles under ``modules_dir``.

        Recognises both the classic ``Dockerfile`` name and the
        ``<name>.dockerfile`` / ``<name>.Dockerfile`` variants.

        Returns
        -------
        List[str]
            Sorted list of Dockerfile paths.
        """
        seen: set = set()
        results: List[str] = []
        for pattern in ("Dockerfile", "*.dockerfile", "*.Dockerfile"):
            for df in self.modules_dir.rglob(pattern):
                resolved = str(df.resolve())
                if resolved not in seen:
                    seen.add(resolved)
                    results.append(str(df))
        results.sort()
        return results

    def process_all(
        self,
        no_cache: bool = False,
        skip_existing: bool = True,
        build_timeout: Optional[int] = None,
        save_timeout: Optional[int] = None,
    ) -> List[Dict[str, str]]:
        """Process every Dockerfile found under ``modules_dir``.

        Parameters
        ----------
        no_cache : bool
        skip_existing : bool
            Skip modules whose tar archive already exists (default: True).
        build_timeout : Optional[int]
        save_timeout : Optional[int]

        Returns
        -------
        List[Dict[str, str]]
            Results from :meth:`process_dockerfile` for each module.
        """
        dockerfiles = self.find_dockerfiles()
        logger.info(f"Found {len(dockerfiles)} Dockerfile(s) under {self.modules_dir}")

        results: List[Dict[str, str]] = []
        failed: List[Tuple[str, str]] = []

        for i, df in enumerate(dockerfiles, 1):
            rel = self.get_module_relpath(df)
            logger.info(f"[{i}/{len(dockerfiles)}] Processing: {rel}")
            try:
                result = self.process_dockerfile(
                    dockerfile_path=df,
                    no_cache=no_cache,
                    skip_existing=skip_existing,
                    build_timeout=build_timeout,
                    save_timeout=save_timeout,
                )
                results.append(result)
            except Exception as e:
                logger.error(f"[FAIL] {rel}: {e}")
                failed.append((df, str(e)))

        logger.info(
            f"Summary: {len(results)}/{len(dockerfiles)} succeeded, "
            f"{len(failed)} failed"
        )
        for df, err in failed:
            logger.error(f"  FAILED: {df} -> {err}")

        return results

    def process_all_def(
        self,
        fakeroot: bool = True,
        force: bool = True,
        skip_existing: bool = True,
        build_timeout: Optional[int] = None,
    ) -> List[Dict[str, str]]:
        """Build SIF from every ``.def`` file found under ``modules_dir``.

        Parameters
        ----------
        fakeroot : bool
        force : bool
        skip_existing : bool
            Skip modules whose SIF already exists (default: True).
        build_timeout : Optional[int]

        Returns
        -------
        List[Dict[str, str]]
            Results from :meth:`process_def_file` for each ``.def`` file.
        """
        def_files = self.find_def_files()
        logger.info(f"Found {len(def_files)} .def file(s) under {self.modules_dir}")

        results: List[Dict[str, str]] = []
        failed: List[Tuple[str, str]] = []

        for i, df in enumerate(def_files, 1):
            rel = self.get_module_relpath(df)
            logger.info(f"[{i}/{len(def_files)}] Processing: {rel}")
            try:
                result = self.process_def_file(
                    def_path=df,
                    fakeroot=fakeroot,
                    force=force,
                    skip_existing=skip_existing,
                    build_timeout=build_timeout,
                )
                results.append(result)
            except Exception as e:
                logger.error(f"[FAIL] {rel}: {e}")
                failed.append((df, str(e)))

        logger.info(
            f"Summary: {len(results)}/{len(def_files)} succeeded, "
            f"{len(failed)} failed"
        )
        for df, err in failed:
            logger.error(f"  FAILED: {df} -> {err}")

        return results


def is_path_like(value: object) -> bool:
    """Return True if *value* looks like a local filesystem path."""
    if not isinstance(value, str):
        return False

    value = value.strip()
    if not value:
        return False

    # URL is not a local filesystem path.
    parsed = urlparse(value)
    if parsed.scheme.lower() in {
        "http",
        "https",
        "ftp",
        "s3",
        "ssh",
    }:
        return False

    # Obvious shell/pipeline expressions.
    if any(x in value for x in (" -> ", " | ", " && ", " || ")):
        return False

    # Windows absolute paths.
    if re.match(r"^[A-Za-z]:[\\/]", value):
        return True

    # UNC path.
    if value.startswith("\\\\"):
        return True

    # Unix absolute path.
    if value.startswith("/"):
        return True

    # Explicit relative path.
    if value.startswith(("./", "../", "~/")):
        return True

    return False
# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        prog="EnvUtil.py",
        description=(
            "Build Apptainer SIF / Docker images from conda env YAMLs.\n\n"
            "Two backends, each with three actions:\n"
            "  apptainer gen    YAML -> .def\n"
            "  apptainer build  .def  -> .sif\n"
            "  apptainer all    YAML -> .def -> .sif  (build file synced to output)\n"
            "  docker gen       YAML -> .dockerfile\n"
            "  docker build     .dockerfile -> image -> .tar\n"
            "  docker all       YAML -> .dockerfile -> image -> .tar"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="backend", required=True)

    # --- args shared by ALL subcommands ---
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "-m", "--modules-dir",
        type=str,
        default=str(Path(__file__).resolve().parents[3] / "modules"),
        help="Root modules directory (default: .../modules).",
    )
    common.add_argument(
        "-s", "--skip-existing",
        action="store_true",
        default=True,
        help="Skip if the image already exists (default: True).",
    )
    common.add_argument(
        "--no-skip-existing",
        dest="skip_existing",
        action="store_false",
        help="Rebuild even if the image already exists.",
    )
    common.add_argument(
        "-t", "--build-timeout",
        type=int,
        default=None,
        help="Build timeout in seconds.",
    )

    # --- args shared by apptainer subcommands ---
    appt_common = argparse.ArgumentParser(add_help=False)
    appt_common.add_argument(
        "--no-fakeroot",
        action="store_false",
        dest="fakeroot",
        default=True,
        help="Disable --fakeroot in apptainer build.",
    )

    # --- args shared by docker subcommands ---
    docker_common = argparse.ArgumentParser(add_help=False)
    docker_common.add_argument(
        "--no-cache",
        action="store_true",
        help="Pass --no-cache to docker build.",
    )
    docker_common.add_argument(
        "--save-timeout",
        type=int,
        default=None,
        help="docker save timeout in seconds.",
    )

    # --- args for 'gen' actions (YAML -> build file, output next to YAML) ---
    gen_common = argparse.ArgumentParser(add_help=False)
    gen_common.add_argument(
        "-y", "--yaml",
        type=str,
        default=None,
        help="Single conda YAML to process. If omitted, process all under --modules-dir.",
    )
    gen_common.add_argument(
        "-e", "--env-name",
        type=str,
        default=None,
        help="Override env / image name (single YAML mode only).",
    )

    # --- args for 'build'/'all' actions (image output requires -o) ---
    image_common = argparse.ArgumentParser(add_help=False)
    image_common.add_argument(
        "-o", "--output-dir",
        type=str,
        required=True,
        help=(
            "Output directory for images and copied build files / YAMLs. "
            "The module-level hierarchy under --modules-dir is preserved."
        ),
    )

    # --- args for 'build' actions (build file -> image) ---
    build_file_common = argparse.ArgumentParser(add_help=False)
    build_file_common.add_argument(
        "-f", "--file",
        type=str,
        default=None,
        help=(
            "Single build file (.def / .dockerfile) to build from. "
            "If omitted, process all found under --modules-dir."
        ),
    )
    build_file_common.add_argument(
        "-e", "--env-name",
        type=str,
        default=None,
        help="Override env / image name (single file mode only).",
    )

    # --- args for 'all' actions (YAML -> image) ---
    all_common = argparse.ArgumentParser(add_help=False)
    all_common.add_argument(
        "-y", "--yaml",
        type=str,
        default=None,
        help="Single conda YAML to process. If omitted, process all under --modules-dir.",
    )
    all_common.add_argument(
        "-e", "--env-name",
        type=str,
        default=None,
        help="Override env / image name (single YAML mode only).",
    )

    # ===== apptainer subcommands =====
    p_appt = sub.add_parser(
        "apptainer",
        help="Apptainer backend: YAML -> .def -> .sif",
    )
    appt_sub = p_appt.add_subparsers(dest="action", required=True)

    appt_sub.add_parser(
        "gen",
        parents=[common, appt_common, gen_common],
        help="Generate .def files from conda YAMLs (no build).",
    )
    appt_sub.add_parser(
        "build",
        parents=[common, appt_common, image_common, build_file_common],
        help="Build SIF images from existing .def files.",
    )
    appt_sub.add_parser(
        "all",
        parents=[common, appt_common, image_common, all_common],
        help="Full pipeline: YAML -> .def -> SIF (build file synced to output).",
    )

    # ===== docker subcommands =====
    p_docker = sub.add_parser(
        "docker",
        help="Docker backend: YAML -> .dockerfile -> image -> .tar",
    )
    docker_sub = p_docker.add_subparsers(dest="action", required=True)

    docker_sub.add_parser(
        "gen",
        parents=[common, docker_common, gen_common],
        help="Generate .dockerfile files from conda YAMLs (no build).",
    )
    docker_sub.add_parser(
        "build",
        parents=[common, docker_common, image_common, build_file_common],
        help="Build Docker images from existing .dockerfile files and save as tar.",
    )
    docker_sub.add_parser(
        "all",
        parents=[common, docker_common, image_common, all_common],
        help="Full pipeline: YAML -> .dockerfile -> image -> tar.",
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Validate -y / -f file arguments
    # ------------------------------------------------------------------
    yaml_arg = getattr(args, "yaml", None)
    if yaml_arg:
        p = Path(yaml_arg)
        if p.suffix not in (".yaml", ".yml"):
            parser.error(
                f"-y/--yaml expects a conda environment YAML file (.yaml/.yml), "
                f"got: {p.name} (suffix={p.suffix or 'none'}). "
                f"Use -f/--file for build files (.def/.dockerfile)."
            )
        if not p.is_file():
            parser.error(f"-y/--yaml file not found: {p}")
        # Lightweight content check: must have channels: or dependencies:
        text = p.read_text(encoding="utf-8")
        if "channels:" not in text and "dependencies:" not in text:
            parser.error(
                f"-y/--yaml does not look like a conda environment file: {p}\n"
                f"  Expected top-level 'channels:' or 'dependencies:' key."
            )

    file_arg = getattr(args, "file", None)
    if file_arg:
        p = Path(file_arg)
        expected = {".def"} if args.backend == "apptainer" else {".dockerfile"}
        if p.suffix not in expected:
            parser.error(
                f"-f/--file expects a {'/'.join(sorted(expected))} file for "
                f"{args.backend} backend, got: {p.name} (suffix={p.suffix or 'none'})."
            )
        if not p.is_file():
            parser.error(f"-f/--file file not found: {p}")

    # gen: no output-dir needed (build files go next to source YAMLs)
    # build/all: output-dir is required (image_common enforces this)
    output_dir = getattr(args, "output_dir", None) or args.modules_dir

    util = EnvUtil(
        modules_dir=args.modules_dir,
        output_dir=output_dir,
    )

    # ---- apptainer ----
    if args.backend == "apptainer":
        if args.action == "gen":
            if args.yaml:
                def_path = util.generate_def_file(
                    yaml_path=args.yaml,
                    env_name=args.env_name,
                )
                logger.info(f"Generated: {def_path}")
            else:
                yamls = util.find_all_env_yamls()
                logger.info(f"Found {len(yamls)} conda env YAML(s)")
                for y in yamls:
                    try:
                        util.generate_def_file(yaml_path=str(y))
                    except Exception as e:
                        logger.error(f"[FAIL] {y}: {e}")

        elif args.action == "build":
            if args.file:
                result = util.process_def_file(
                    def_path=args.file,
                    env_name=args.env_name,
                    fakeroot=args.fakeroot,
                    skip_existing=args.skip_existing,
                    build_timeout=args.build_timeout,
                )
                logger.info(f"Result: {result}")
            else:
                results = util.process_all_def(
                    fakeroot=args.fakeroot,
                    skip_existing=args.skip_existing,
                    build_timeout=args.build_timeout,
                )
                logger.info(f"Processed {len(results)} .def file(s)")

        elif args.action == "all":
            if args.yaml:
                result = util.process_yaml(
                    yaml_path=args.yaml,
                    env_name=args.env_name,
                    skip_existing=args.skip_existing,
                    fakeroot=args.fakeroot,
                    build_timeout=args.build_timeout,
                )
                logger.info(f"Result: {result}")
            else:
                results = util.process_all_apptainer(
                    skip_existing=args.skip_existing,
                    fakeroot=args.fakeroot,
                    build_timeout=args.build_timeout,
                )
                logger.info(f"Processed {len(results)} YAML(s)")

    # ---- docker ----
    elif args.backend == "docker":
        if args.action == "gen":
            if args.yaml:
                df_path = util.generate_dockerfile(
                    yaml_path=args.yaml,
                    image_name=args.env_name,
                )
                logger.info(f"Generated: {df_path}")
            else:
                yamls = util.find_all_env_yamls()
                logger.info(f"Found {len(yamls)} conda env YAML(s)")
                for y in yamls:
                    try:
                        util.generate_dockerfile(yaml_path=str(y))
                    except Exception as e:
                        logger.error(f"[FAIL] {y}: {e}")

        elif args.action == "build":
            if args.file:
                result = util.process_dockerfile(
                    dockerfile_path=args.file,
                    image_name=args.env_name,
                    no_cache=args.no_cache,
                    skip_existing=args.skip_existing,
                    build_timeout=args.build_timeout,
                    save_timeout=args.save_timeout,
                )
                logger.info(f"Result: {result}")
            else:
                results = util.process_all(
                    no_cache=args.no_cache,
                    skip_existing=args.skip_existing,
                    build_timeout=args.build_timeout,
                    save_timeout=args.save_timeout,
                )
                logger.info(f"Processed {len(results)} dockerfile(s)")

        elif args.action == "all":
            if args.yaml:
                result = util.process_yaml_docker(
                    yaml_path=args.yaml,
                    image_name=args.env_name,
                    skip_existing=args.skip_existing,
                    no_cache=args.no_cache,
                    build_timeout=args.build_timeout,
                    save_timeout=args.save_timeout,
                )
                logger.info(f"Result: {result}")
            else:
                results = util.process_all_docker(
                    skip_existing=args.skip_existing,
                    no_cache=args.no_cache,
                    build_timeout=args.build_timeout,
                    save_timeout=args.save_timeout,
                )
                logger.info(f"Processed {len(results)} YAML(s)")


if __name__ == "__main__":
    main()
