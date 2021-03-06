import glob
import argparse
from PanACoTA.subcommands.pangenome import build_parser
import pandas as pd
import os
from datetime import datetime, timezone
import pytz
from pathlib import Path


class ThrowingArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise Exception(message)


def get_acs(outdir):
    names_map = pd.read_table(Path(outdir)/"names_mapping.tsv", sep="\t", header=None)
    return list(names_map[0])


def read_params_from_log(logpath):
    with open(logpath, "r") as log:
        filtered = filter(lambda line: line.strip().startswith("> PanACoTA pangenome"), log)
        command = next(filtered).strip()
        cmd_parser = ThrowingArgumentParser(add_help=False)
        build_parser(cmd_parser)
        try:
            parsed_args = cmd_parser.parse_args(command.split()[3:])
            parsed_args.mod_time = datetime.fromtimestamp(os.stat(logpath).st_mtime, tz=pytz.timezone("Europe/Moscow"))
            acs = get_acs(parsed_args.outdir)
            parsed_args.acs = ",".join(acs)
            return parsed_args
        except:
            pass


def find_data():
    logfiles = glob.glob("*_pan/PanACoTA-pangenome_*.log")
    arg_map = map(read_params_from_log, logfiles)
    args_valid = filter(lambda x: x is not None, arg_map)
    dict_list = list(map(vars, args_valid))
    available_data = pd.DataFrame(dict_list)
    return available_data.sort_values(by="mod_time", ascending=False, ignore_index=True)


if __name__ == "__main__":
    data_available = find_data()
    data_available.to_csv("available_data.tsv", sep="\t", index=False)
