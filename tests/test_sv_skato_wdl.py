"""Check workflow-level expressions and render both commands with localized paths."""

import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import WDL


ROOT = Path(__file__).resolve().parents[1]


def descendants(node):
    yield node
    for child in node.children:
        yield from descendants(child)


class LocalCommandStdLib(WDL.StdLib.Base):
    def _virtualize_filename(self, filename):
        return filename


def render_and_capture(task, inputs, paths, directory):
    values = WDL.values_from_json(inputs, task.available_inputs, task.required_inputs)
    localized = WDL.Value.rewrite_env_paths(values, lambda value: str(paths[value.value]))
    stdlib = LocalCommandStdLib("1.0", write_dir=str(directory))
    rendered = task.command.eval(localized, stdlib).value
    capture = directory / "captured.json"
    fake_rscript = directory / "Rscript"
    fake_rscript.write_text(
        f"#!{sys.executable}\n"
        "import json, os, sys\n"
        "with open(os.environ['CAPTURE'], 'w') as handle:\n"
        "    json.dump(sys.argv[1:], handle)\n",
        encoding="utf-8",
    )
    fake_rscript.chmod(0o755)
    env = dict(os.environ, CAPTURE=str(capture), PATH=str(directory) + os.pathsep + os.environ["PATH"])
    subprocess.run(["bash", "-c", rendered], cwd=directory, env=env, check=True)
    return json.loads(capture.read_text(encoding="utf-8"))


def check_commands(document, directory):
    tasks = {task.name: task for task in document.tasks}
    runtime = dict(docker_image="test/image:1", task_memory="2 GB", disk_gb=10)
    uris = {name: f"gs://test/{name}" for name in
            ("pgen", "pvar", "psam", "molecular_traits", "gtf", "lof", "cov", "chr1", "chr2")}
    paths = {}
    for name, uri in uris.items():
        paths[uri] = directory / f"local {name}'s $quoted file"
        paths[uri].touch()

    scan = dict(
        runtime,
        files={**{key: uris[key] for key in ("pgen", "pvar", "psam", "molecular_traits")},
               "chromosome": "chr1"},
        cds_gtf=uris["gtf"],
        lof_table=uris["lof"],
        covariates=uris["cov"],
        covariate_columns="age,PC1",
        lof_missing_as_noncarrier=False,
        max_alt_af=0.01,
        max_missing=0.05,
    )
    args = render_and_capture(tasks["test_chromosome"], scan, paths, directory)
    assert args[:2] == ["/opt/sv-skato/sv_pgen_workflow.R", "scan"]
    for flag, key in (("pgen", "pgen"), ("pvar", "pvar"), ("psam", "psam"),
                      ("traits", "molecular_traits"), ("gtf", "gtf"), ("lof", "lof"),
                      ("covariates", "cov")):
        assert args[args.index("--" + flag) + 1] == str(paths[uris[key]])
    assert "--unknown-unlisted-lof" in args
    args = render_and_capture(tasks["test_chromosome"],
                              dict(scan, covariates=None, covariate_columns="", lof_missing_as_noncarrier=True),
                              paths, directory)
    assert all(flag not in args for flag in ("--covariates", "--covariate-cols", "--unknown-unlisted-lof"))

    results = [uris["chr1"], uris["chr2"]]
    args = render_and_capture(tasks["gather_results"], dict(runtime, masks=results, fdr_method="BH"),
                              paths, directory)
    assert args[:2] == ["/opt/sv-skato/sv_pgen_workflow.R", "gather"]
    assert args[-2:] == ["--fdr-method", "BH"]
    contents = Path(args[args.index("--mask-list") + 1]).read_text(encoding="utf-8")
    assert contents.splitlines() == [str(paths[uri]) for uri in results]


def check_manifest(document, directory):
    manifest = directory / "manifest.tsv"
    header = "chromosome\tpgen\tpvar\tpsam\tmolecular_traits\n"
    row = "chr1\tgs://test/chr1.pgen\tgs://test/chr1.pvar\tgs://test/chr1.psam\tgs://test/rna.tsv\n"
    manifest.write_text(header + row, encoding="utf-8")
    output = directory / "inputs.json"
    command = [
        sys.executable, str(ROOT / "scripts/prepare_sv_skato_inputs.py"),
        "--manifest", str(manifest), "--cds-gtf", "gs://test/cds.gtf",
        "--lof-table", "gs://test/lof.tsv", "--docker-image", "test/image:1", "--output", str(output),
    ]
    subprocess.run(command, check=True)
    values = WDL.values_from_json(json.loads(output.read_text(encoding="utf-8")),
                                  document.workflow.available_inputs, document.workflow.required_inputs,
                                  namespace="sv_skato")
    entry = values["chromosome_manifest"].value[0]
    assert isinstance(entry.value["molecular_traits"], WDL.Value.File)
    assert entry.value["molecular_traits"].value == "gs://test/rna.tsv"
    manifest.write_text(header + row + row.replace("chr1\t", "1\t", 1), encoding="utf-8")
    duplicate = subprocess.run(command, capture_output=True, text=True)
    assert duplicate.returncode and "repeated chromosome" in duplicate.stderr


def main():
    document = WDL.load(str(ROOT / "wdl/sv_skato.wdl"))
    assert document.wdl_version == "1.0"
    assert not any(isinstance(node, WDL.Expr.Apply) and node.function_name.startswith("write_")
                   for node in descendants(document.workflow))
    with TemporaryDirectory() as temp:
        check_commands(document, Path(temp))
        check_manifest(document, Path(temp))
    print("PASS: workflow-scope file creation, typed manifest, and both task localization paths.")


if __name__ == "__main__":
    main()
