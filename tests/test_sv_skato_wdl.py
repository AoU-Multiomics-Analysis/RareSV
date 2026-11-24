"""Check workflow-level expressions and render commands with localized paths."""

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


class CloudCommandStdLib(WDL.StdLib.Base):
    def _virtualize_filename(self, filename):
        return f"gs://test/task-commands/{Path(filename).name}"


def render_command(task, inputs, paths, directory, stdlib_class=LocalCommandStdLib):
    values = WDL.values_from_json(inputs, task.available_inputs, task.required_inputs)
    localized = WDL.Value.rewrite_env_paths(values, lambda value: str(paths[value.value]))
    stdlib = stdlib_class("1.0", write_dir=str(directory))
    return task.command.eval(localized, stdlib).value


def render_and_capture(task, inputs, paths, directory, stdlib_class=LocalCommandStdLib):
    rendered = render_command(task, inputs, paths, directory, stdlib_class)
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


def check_commands(document, directory, files):
    tasks = {task.name: task for task in document.tasks}
    runtime = dict(docker_image="test/image:1", task_memory="2 GB", disk_gb=10)
    uris = {name: f"gs://test/{name}" for name in
            ("gtf", "lof", "cov", "chr1", "chr2", "parsed_manifest", "significant_genes")}
    uris.update({key: uri for key, uri in files.items() if key != "chromosome"})
    paths = {}
    for name, uri in uris.items():
        paths[uri] = directory / f"local {name}'s $quoted file"
        paths[uri].touch()

    scan = dict(
        runtime,
        **files,
        cds_gtf=uris["gtf"],
        lof_table=uris["lof"],
        significant_genes=None,
        covariates=uris["cov"],
        covariate_columns="age,PC1",
        lof_missing_as_noncarrier=False,
        missing_pgen_as_reference=True,
        window_padding_bp=250,
        max_alt_af=0.01,
        max_missing=0.05,
    )
    args = render_and_capture(tasks["test_chromosome"], scan, paths, directory)
    assert args[:2] == ["/opt/sv-skato/sv_pgen_workflow.R", "scan"]
    for flag, key in (("pgen", "pgen"), ("pvar", "pvar"), ("psam", "psam"),
                      ("phenotype-bed", "phenotype_bed"), ("gtf", "gtf"), ("lof", "lof"),
                      ("covariates", "cov")):
        assert args[args.index("--" + flag) + 1] == str(paths[uris[key]])
    assert "--unknown-unlisted-lof" in args
    assert "--missing-pgen-as-reference" in args
    assert "--significant-genes" not in args
    assert args[args.index("--window-padding-bp") + 1] == "250"
    follow_up = dict(scan, significant_genes=uris["significant_genes"])
    args = render_and_capture(tasks["test_chromosome"], follow_up, paths, directory)
    assert args[args.index("--significant-genes") + 1] == str(paths[uris["significant_genes"]])
    args = render_and_capture(tasks["test_chromosome"],
                              dict(scan, cds_gtf=None, covariates=None, covariate_columns="",
                                   lof_missing_as_noncarrier=True, missing_pgen_as_reference=False),
                              paths, directory)
    assert all(flag not in args for flag in ("--gtf", "--covariates", "--covariate-cols",
                                            "--unknown-unlisted-lof", "--missing-pgen-as-reference"))

    results = [uris["chr1"], uris["chr2"]]
    paths[results[0]] = directory / "mask's file $(touch unexpected-command).tsv"
    paths[results[0]].touch()
    args = render_and_capture(tasks["gather_results"],
                              dict(runtime, masks=results, manifest=uris["parsed_manifest"], fdr_method="BH"),
                              paths, directory, CloudCommandStdLib)
    assert args[:2] == ["/opt/sv-skato/sv_pgen_workflow.R", "gather"]
    assert args[-2:] == ["--fdr-method", "BH"]
    assert args[args.index("--manifest") + 1] == str(paths[uris["parsed_manifest"]])
    mask_list = args[args.index("--mask-list") + 1]
    assert "://" not in mask_list
    contents = (directory / mask_list).read_text(encoding="utf-8")
    assert contents.splitlines() == [str(paths[uri]) for uri in results]
    assert not (directory / "unexpected-command").exists()
    for files in (results, []):
        args = render_and_capture(tasks["gather_follow_up"], dict(runtime, files=files),
                                  paths, directory, CloudCommandStdLib)
        assert args[:2] == ["/opt/sv-skato/sv_pgen_workflow.R", "gather-follow-up"]
        follow_up_list = args[args.index("--follow-up-list") + 1]
        assert "://" not in follow_up_list
        assert (directory / follow_up_list).read_text(encoding="utf-8").splitlines() == [
            str(paths[uri]) for uri in files
        ]


def check_follow_up_rows(document, directory):
    task = next(task for task in document.tasks if task.name == "gather_results")
    output = next(output for output in task.outputs if output.name == "follow_up_rows")

    class OutputStdLib(WDL.StdLib.Base):
        def _devirtualize_filename(self, filename):
            return str(directory / filename)

    stdlib = OutputStdLib("1.0")
    manifest = directory / "follow_up_manifest.tsv"
    manifest.touch()
    assert output.expr.eval(WDL.Env.Bindings(), stdlib).json == []
    row = ["chr21", "gs://test/21.pgen", "gs://test/21.pvar", "gs://test/21.psam", "gs://test/21.bed"]
    manifest.write_text("\t".join(row) + "\n", encoding="utf-8")
    assert output.expr.eval(WDL.Env.Bindings(), stdlib).json == [row]


def check_manifest(document, directory):
    manifest = directory / "chr21 manifest's $quoted path.tsv"
    uri = "gs://test/chr21_manifest.tsv"
    output = directory / "inputs.json"
    command = [
        sys.executable, str(ROOT / "scripts/prepare_sv_skato_inputs.py"),
        "--manifest", uri,
        "--lof-table", "gs://test/lof.tsv", "--output", str(output),
    ]
    subprocess.run(command, check=True)
    inputs = json.loads(output.read_text(encoding="utf-8"))
    assert "sv_skato.cds_gtf" not in inputs
    assert "sv_skato.window_padding_bp" not in inputs
    assert "sv_skato.docker_image" not in inputs
    assert "sv_skato.missing_pgen_as_reference" not in inputs
    assert "sv_skato.follow_up" not in inputs
    values = WDL.values_from_json(inputs,
                                  document.workflow.available_inputs, document.workflow.required_inputs,
                                  namespace="sv_skato")
    assert isinstance(values["chromosome_manifest"], WDL.Value.File)
    assert values["chromosome_manifest"].value == uri
    subprocess.run(command + ["--cds-gtf", "gs://test/cds.gtf", "--window-padding-bp", "0",
                              "--docker-image", "test/image:2", "--missing-pgen-as-reference",
                              "--no-follow-up"], check=True)
    values = WDL.values_from_json(json.loads(output.read_text(encoding="utf-8")),
                                  document.workflow.available_inputs, document.workflow.required_inputs,
                                  namespace="sv_skato")
    assert isinstance(values["cds_gtf"], WDL.Value.File)
    assert values["cds_gtf"].value == "gs://test/cds.gtf"
    assert values["window_padding_bp"].value == 0
    assert values["docker_image"].value == "test/image:2"
    assert values["missing_pgen_as_reference"].value is True
    assert values["follow_up"].value is False
    task = next(task for task in document.tasks if task.name == "parse_manifest")
    inputs = {"manifest": uri, "docker_image": "test/image:1"}
    rendered = render_command(task, inputs, {uri: manifest}, directory)
    header = "phenotype_bed\tpsam\tchromosome\tpgen\tpvar\n"
    row = "gs://test/rna's $file.bed.gz\tgs://test/chr21.psam\tchr21\tgs://test/chr21.pgen\tgs://test/chr21.pvar\n"
    manifest.write_text("\ufeff" + (header + row).replace("\n", "\r\n"), encoding="utf-8")
    subprocess.run(["bash", "-c", rendered], cwd=directory, check=True)
    fields = (directory / "manifest_rows.tsv").read_text(encoding="utf-8").rstrip("\n").split("\t")
    assert fields == ["chr21", "gs://test/chr21.pgen", "gs://test/chr21.pvar",
                      "gs://test/chr21.psam", "gs://test/rna's $file.bed.gz"]
    scatter = next(node for node in document.workflow.body if isinstance(node, WDL.Tree.Scatter))
    row_value = WDL.Value.Array(WDL.Type.String(), [WDL.Value.String(field) for field in fields])
    env = WDL.Env.Bindings().bind(scatter.variable, row_value)
    stdlib = LocalCommandStdLib("1.0", write_dir=str(directory))
    for declaration in scatter.body:
        if isinstance(declaration, WDL.Tree.Decl):
            env = env.bind(declaration.name, declaration.expr.eval(env, stdlib).coerce(declaration.type))
    scan_task = next(task for task in document.tasks if task.name == "test_chromosome")
    scan_call = next(node for node in scatter.body if isinstance(node, WDL.Tree.Call)
                     and node.name == "test_chromosome")
    scan_inputs = {"chromosome": scan_call.inputs["chromosome"].eval(env, stdlib).json}
    for key, uri_value in zip(("pgen", "pvar", "psam", "phenotype_bed"), fields[1:]):
        assert isinstance(scan_task.required_inputs[key].type, WDL.Type.File)
        value = scan_call.inputs[key].eval(env, stdlib)
        assert isinstance(value, WDL.Value.File)
        assert value.value == uri_value
        scan_inputs[key] = value.json

    manifest.write_text(header + row + row.replace("\tchr21\t", "\t21\t"), encoding="utf-8")
    duplicate = subprocess.run(["bash", "-c", rendered], cwd=directory, capture_output=True, text=True)
    assert duplicate.returncode and "repeated chromosome" in duplicate.stderr
    manifest.write_text(header + row.replace("gs://test/chr21.psam", "/home/rstudio/chr21.psam"), encoding="utf-8")
    local_uri = subprocess.run(["bash", "-c", rendered], cwd=directory, capture_output=True, text=True)
    assert local_uri.returncode and "expected a gs:// path for psam" in local_uri.stderr
    manifest.write_text(header, encoding="utf-8")
    empty = subprocess.run(["bash", "-c", rendered], cwd=directory, capture_output=True, text=True)
    assert empty.returncode and "expected at least one chromosome row" in empty.stderr
    unresolved = render_command(task, inputs, {uri: uri}, directory)
    result = subprocess.run(["bash", "-c", unresolved], cwd=directory, capture_output=True, text=True)
    assert result.returncode and "Localization error for chromosome manifest" in result.stderr
    return scan_inputs


def main():
    document = WDL.load(str(ROOT / "wdl/sv_skato.wdl"))
    assert document.wdl_version == "1.0"
    assert not any(isinstance(node, WDL.Expr.Apply) and node.function_name.startswith("write_")
                   for node in descendants(document.workflow))
    with TemporaryDirectory() as temp:
        files = check_manifest(document, Path(temp))
        check_commands(document, Path(temp), files)
        check_follow_up_rows(document, Path(temp))
    print("PASS: workflow-scope file creation, TSV validation, and cloud-to-task file localization.")


if __name__ == "__main__":
    main()
