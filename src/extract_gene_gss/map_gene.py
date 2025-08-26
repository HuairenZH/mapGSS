import subprocess
from pathlib import Path

def run(gene: str):
    """
    调用 R 脚本 map_gene.R，并传入参数 --gene
    """
    script_path = Path(__file__).resolve().parent / "map_gene.R"

    print(f"🔎 Running map_gene pipeline for {gene} (via R script {script_path})")

    # 调用 Rscript
    result = subprocess.run(
        ["Rscript", str(script_path), "--gene", gene],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print("❌ R script failed!")
        print(result.stderr)
    else:
        print("✅ R script finished successfully")
        print(result.stdout)