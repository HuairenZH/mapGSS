#!/usr/bin/env python3
import argparse
import os
import shutil
from pathlib import Path
import pandas as pd
from pathlib import Path

# ===== 固定路径（按你当前代码保持不变） =====
# 找到项目根目录（假设 cli.py 在 src/extract_gene_gss/ 里）
PROJECT_ROOT = Path(__file__).resolve().parents[2]   # 回到 extract-gene-gss/

DATA_DIR = PROJECT_ROOT / "data"

EMBRYO_FEATHER = DATA_DIR / "ST/E16.5_E1S1/E16.5_E1S1_gene_marker_score.feather"
EMBRYO_META    = DATA_DIR / "ST/E16.5_E1S1/mouse_E16.5_E1S1.txt"

ADULT_FEATHER  = DATA_DIR / "ST/mouse_Adult_Mouse_brain_cell_bin/mouse_Adult_Mouse_brain_cell_bin_gene_marker_score.feather"
ADULT_META     = DATA_DIR / "ST/mouse_Adult_Mouse_brain_cell_bin/mouse_Adult_Mouse_brain_cell_bin.txt"

BANNER = r"""
 __  __             _____    _____   _____   _____ 
|  \/  |     /\    |  __ \  / ____| / ____| / ____|
| \  / |    /  \   | |__) || |  __ | (___  | (___  
| |\/| |   / /\ \  |  ___/ | | |_ | \___ \  \___ \ 
| |  | |  / ____ \ | |     | |__| | ____) | ____) |
|_|  |_| /_/    \_\|_|      \_____| |_____/ |_____/ 



                
"""

def center_banner(text: str) -> str:
    width = shutil.get_terminal_size((80, 20)).columns
    return "\n".join(line.center(width) for line in text.splitlines())




def _extract_and_merge(feather_path: Path, meta_path: Path, gene: str) -> pd.DataFrame:
    """
    从 feather 里取出目标基因一行，转置成两列 [cell_name, {gene}_GSS]，
    与 meta（按 cell_name）内连接，并删除 Lambda/Depression 列。
    """
    if not feather_path.exists():
        raise FileNotFoundError(f"Feather 不存在: {feather_path}")
    if not meta_path.exists():
        raise FileNotFoundError(f"Meta 不存在: {meta_path}")

    # 读 marker matrix（基因名在第 1 列）
    df_marker = pd.read_feather(feather_path)
    df_gene = df_marker[df_marker.iloc[:, 0] == gene]

    if df_gene.empty:
        raise ValueError(f"在 {feather_path} 中找不到基因：{gene}")

    # 转置：第一列成为列名，第一行是原列名，需丢弃
    df_gene_T = df_gene.T.reset_index()
    df_gene_T.columns = ["cell_name", f"{gene}_GSS"]
    df_gene_T = df_gene_T[1:]  # 去掉第一行（原表头）

    # 读 meta（tab 分隔）
    df_meta = pd.read_csv(meta_path, sep="\t")

    # 合并
    df_merged = pd.merge(df_meta, df_gene_T, on="cell_name", how="inner")

    # 删除多余列（若不存在则忽略）
    df_merged = df_merged.drop(columns=["Lambda", "Depression"], errors="ignore")
    return df_merged


def run_pipeline(gene: str) -> Path:
    """
    跑 embryo 与 adult 两个数据集，输出到 ../data/Genes/{gene}/
    返回输出目录路径。
    """
    out_dir = Path(f"./data/Genes/{gene}/")
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) embryo
    df_embryo = _extract_and_merge(EMBRYO_FEATHER, EMBRYO_META, gene)
    (out_dir / f"mouse_E16.5_E1S1_{gene}.csv").write_text(
        df_embryo.to_csv(index=False), encoding="utf-8"
    )

    # 2) adult
    df_adult = _extract_and_merge(ADULT_FEATHER, ADULT_META, gene)
    (out_dir / f"mouse_Adult_Mouse_brain_cell_bin_{gene}.csv").write_text(
        df_adult.to_csv(index=False), encoding="utf-8"
    )

    return out_dir


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Extract per-gene GSS and merge with metadata for two fixed ST datasets"
    )
    p.add_argument(
        "--gene", "-g", required=True, help="基因名（如 SLC32A1）"
    )
    return p


def main(argv=None):
    print(center_banner(BANNER))
    parser = build_parser()
    args = parser.parse_args(argv)

    out_dir = run_pipeline(args.gene)
    
    print(f"✅ Completed {args.gene}, results have been saved to: {out_dir}")


if __name__ == "__main__":
    main()