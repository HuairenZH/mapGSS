#!/usr/bin/env python3
import argparse
from . import map_gene   # 调用 map_gene.py
import shutil

BANNER = r"""
 __  __             _____    _____   _____   _____ 
|  \/  |     /\    |  __ \  / ____| / ____| / ____|
| \  / |    /  \   | |__) || |  __ | (___  | (___  
| |\/| |   / /\ \  |  ___/ | | |_ | \___ \  \___ \ 
| |  | |  / ____ \ | |     | |__| | ____) | ____) |
|_|  |_| /_/    \_\|_|      \_____| |_____/ |_____/ 



                
"""



def build_parser():
    p = argparse.ArgumentParser(description="Run map_gene functionality on a given gene")
    p.add_argument("--gene", "-g", required=True, help="基因名（如 SLC32A1）")
    return p



def center_banner(text: str) -> str:
    width = shutil.get_terminal_size((80, 20)).columns
    return "\n".join(line.center(width) for line in text.splitlines())

def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    # 打印拼凑出来的 ASCII banner
    print(center_banner(BANNER))

    # 调用 map_gene.py 里的函数
    map_gene.run(args.gene)

    print(f"✅ map_gene_GSS finished for {args.gene}")

if __name__ == "__main__":
    main()