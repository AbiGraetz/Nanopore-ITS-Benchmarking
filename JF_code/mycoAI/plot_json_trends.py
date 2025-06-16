#!/usr/bin/env python3
"""
plot_json_trends.py


"""

import argparse
import json
import pathlib
import sys
from typing import Dict, List
import numpy as np

import matplotlib.pyplot as plt

def load_json_files(folder: pathlib.Path, sub_str1:str, sub_str2:str) -> List[pathlib.Path]:
    json_files = []
    for fp in folder.iterdir():
        if fp.suffix.lower() == ".json" and sub_str1 in fp.name and sub_str2 in fp.name:
            json_files.append(fp)
    if not json_files:
        sys.exit("can not find json file")
    return sorted(json_files, key=lambda p: int(p.stem.split('_')[-1]))

def extract_max_values(file_path: pathlib.Path) -> Dict[str, float]:
    with file_path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    if not isinstance(data, dict):
        sys.exit(f"{file_path} not dict")

    result = {}
    for k, sub in data.items():
        if not isinstance(sub, dict) or not sub:
            sys.exit(f"{file_path} key {k} error")
        try:
            result[k] = max(sub.values())
        except TypeError:
            sys.exit(f"{file_path} {k} contains value that not a number")

    return result

def plot_trends(x: List[int], y_dict: Dict[str, List[float]], output_dir:str, fname1:str, fname2:str):
    plt.figure(figsize=(8, 5))
    for key, ys in y_dict.items():
        plt.plot(np.array(x)*5 + 5, ys, marker="o", label=key)

    plt.xlabel("epochs")
    plt.ylabel("correct classified")
    plt.title("classification accuracy by epoches")
    plt.legend()
    
    
    max_x = max(x)
    plt.xticks(np.arange(0, max_x*5 + 1 + 5, 5))



    plt.ylim(1600, 2000)                        
    plt.yticks(np.arange(1600, 2001, 20))    



    plt.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)        
    out_path = output_dir / (fname1 + fname2 + ".png")
    plt.savefig(out_path, dpi=300)                      
    plt.close()                                        

    print(f"saved：{out_path}")

def main():
    parser = argparse.ArgumentParser(description="plot figure")
    parser.add_argument("json_folder", type=pathlib.Path,
                        help="folder cotnains json files")
    parser.add_argument("file_name_sub_string1", type=str,
                        help="which kind of json files you would like to load")
    parser.add_argument("file_name_sub_string2", type=str,
                        help="which kind of json files you would like to load")
    parser.add_argument("out_folder", type=pathlib.Path,
                        help="folder to output the figure")
    args = parser.parse_args()

    if not args.json_folder.is_dir():
        sys.exit(f"not folder：{args.json_folder}")

    json_files = load_json_files(args.json_folder, args.file_name_sub_string1, args.file_name_sub_string2)

    first_vals = extract_max_values(json_files[0])
    keys = list(first_vals.keys())
    if len(keys) != 3:
        sys.exit(f"JSON should have 3 keys， but got {len(keys)}：{keys}")

    trends = {k: [] for k in keys}

    indices = []
    for fp in json_files:
        idx = int(fp.stem.split('_')[-1])       
        indices.append(idx)

        values = extract_max_values(fp)
        if set(values.keys()) != set(keys):
            sys.exit(f"{fp} key does not match")

        for k in keys:
            trends[k].append(values[k])

    plot_trends(indices, trends, args.out_folder, args.file_name_sub_string1, args.file_name_sub_string2)

if __name__ == "__main__":
    main()


