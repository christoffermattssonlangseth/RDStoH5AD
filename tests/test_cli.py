from __future__ import annotations

import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLI = ["python3", "-m", "rdstoh5ad.cli"]
ENV = {**os.environ, "PYTHONPATH": str(ROOT / "src"), "OMP_NUM_THREADS": "1"}


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [*CLI, *args],
        cwd=ROOT,
        env=ENV,
        text=True,
        capture_output=True,
        check=False,
    )


class CliSmokeTests(unittest.TestCase):
    def test_help_lists_subcommands(self) -> None:
        result = run_cli("--help")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("inspect,convert,validate", result.stdout)

    def test_inspect_plain_list_rds(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "toy.rds"
            r_command = (
                "x <- list("
                "obs=data.frame(cluster=c(\"a\",\"b\"), row.names=c(\"cell1\",\"cell2\")),"
                "expression=matrix(c(1,0,2,3), nrow=2, dimnames=list(c(\"gene1\",\"gene2\"), c(\"cell1\",\"cell2\"))),"
                "umap=matrix(c(0,1,2,3), ncol=2, byrow=TRUE, dimnames=list(c(\"cell1\",\"cell2\"), c(\"UMAP_1\",\"UMAP_2\"))),"
                "coordinates=matrix(c(10,20,30,40), ncol=2, byrow=TRUE, dimnames=list(c(\"cell1\",\"cell2\"), c(\"x\",\"y\"))));"
                f"saveRDS(x, \"{input_path}\")"
            )
            created = subprocess.run(
                ["Rscript", "-e", r_command],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(created.returncode, 0, created.stderr)

            result = run_cli("inspect", str(input_path))
            self.assertEqual(result.returncode, 0, result.stderr)

            payload = json.loads(result.stdout)
            self.assertEqual(payload["object_type"], "list")
            self.assertEqual(payload["cells"], 2)
            self.assertEqual(payload["genes"], 2)
            self.assertTrue(payload["has_spatial"])

    def test_validate_missing_file_fails_cleanly(self) -> None:
        result = run_cli("validate", "/tmp/does-not-exist.h5ad")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Validation failed.", result.stderr)


if __name__ == "__main__":
    unittest.main()
