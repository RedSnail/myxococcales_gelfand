import os
import subprocess as sp
import pandas as pd
from PanACoTA.bin.run_panacota import parse_arguments
from PanACoTA.subcommands.pangenome import build_parser, main_from_parse
import argparse as arp
import concurrent.futures as cf
from pathlib import Path
from bot_utils import HOME, send_text, send_image, send_album
from find_data_available import find_data
import shutil
from functools import reduce


def run_sub(subcommand):
    action, args = parse_arguments(subcommand.split())
    action(args)


def get_pan_args(configagrs, pan_parser):
    pandir = HOME / f"{configagrs.alias}_pan"
    lstinfo_path = pandir / f"LSTINFO-LSTINFO-{configagrs.alias[:4]}.lst"
    pangenome_command = (f"pangenome -l {lstinfo_path} -n {configagrs.alias[:4]} -d {pandir / 'Proteins'} "
                         f"-o {pandir} -m proteinortho --threads 8 -v ")
    if configagrs.eval is not None:
        pangenome_command += f"--eval {configagrs.eval} "

    if configagrs.conn is not None:
        pangenome_command += f"--conn {configagrs.conn} "

    if configagrs.purity is not None:
        pangenome_command += f"--purity {configagrs.purity}"

    if configagrs.minspec is not None:
        pangenome_command += f"--minspecies {configagrs.minspec} "

    args = pan_parser.parse_args(pangenome_command.split()[1:])
    args.threads = 8
    args.argv = pangenome_command.split()

    return args


def parse_task(task, metadata, myxo_parser, pan_parser):
    alias = task.split()[0]
    if task.startswith("WithOut"):
        families = pd.unique(metadata["family"])
        add_args = task.split()[1:]
    elif task.startswith("Bulk"):
        mask = ~genome_metadata.family.eq("Geobacteraceae")
        families = pd.unique(genome_metadata[mask]["family"])
        add_args = task.split()[1:]
    else:
        families = task.split()[1].split(",")
        add_args = task.split()[2:]

    acs = pd.unique(metadata[genome_metadata.family.isin(families)]["AC"])

    config_args = myxo_parser.parse_args(add_args)
    config_args.acs = acs
    config_args.alias = alias

    return config_args, get_pan_args(config_args, pan_parser)


def create_name_map(lstinfo: pd.DataFrame, metadata: pd.DataFrame):
    name_map = lstinfo.copy()
    name_map["AC"] = lstinfo.orig_name.str.extract("(GCF_\d*[.,]?\d)_")
    name_map = pd.merge(name_map, metadata, on="AC")
    return name_map[["AC", "gembase_name", "organism"]]


def run_preparation(taxids, alias):
    tmpdir = HOME / f"{alias}_tmp_pan"

    create_dirs(tmpdir)
    run_sub(f"prepare -t {','.join(str(taxid) for taxid in taxids)} -l complete,chromosome "
            f"-o {tmpdir} -p 4 --min_dist 0 --max_dist 1 --norefseq")
    old_primary_lstinfo = next(tmpdir.glob("LSTINFO-*"))
    new_primary_lstifo = f"{tmpdir}/LSTINFO-{alias[:4]}.txt"
    shutil.move(old_primary_lstinfo, new_primary_lstifo)
    run_sub(f"annotate  --info {new_primary_lstifo} -r {tmpdir} --threads 4 -n {alias[:4]}")
    return tmpdir


def create_dirs(pandir):
    pandir.mkdir(exist_ok=True)
    (pandir / "Proteins").mkdir(exist_ok=True)
    (pandir / "Genes").mkdir(exist_ok=True)


def copy_source(old_path, new_path, new_alias, acs):
    create_dirs(new_path)
    old_lstinfo_path = next(old_path.glob("LSTINFO-LSTINFO-*.lst"))
    old_lstinfo = pd.read_table(old_lstinfo_path, sep="\t")

    new_lstinfo_path = new_path / f"LSTINFO-LSTINFO-{new_alias[:4]}.lst"
    if not new_lstinfo_path.is_file():
        new_lstinfo = pd.DataFrame(columns="gembase_name orig_name to_annotate gsize nb_conts L90".split())
    else:
        new_lstinfo = pd.read_table(new_lstinfo_path)

    existing_acs = set(new_lstinfo.orig_name.str.extract("(GCF_\d*[.,]?\d)_")[0])
    acs = acs - existing_acs
    if len(acs) == 0:
        return

    mask = old_lstinfo.orig_name.str.contains("|".join(acs))
    lstinfo_appendix = old_lstinfo[mask]
    lstinfo_appendix = lstinfo_appendix.reset_index(drop=True)
    old_names = list(lstinfo_appendix.gembase_name)
    name_from_i = lambda i: f"{new_alias[:4]}.0921.{(i + 1):05d}"
    new_names = list(map(name_from_i, range(new_lstinfo.shape[0], new_lstinfo.shape[0] + lstinfo_appendix.shape[0])))
    lstinfo_appendix.gembase_name = pd.Series(new_names)
    lstinfo_appendix.orig_name = lstinfo_appendix.orig_name.str.replace(str(old_path), str(new_path))
    (new_path / "tmp_files").mkdir(exist_ok=True)
    for idx, (old_tmp, new_tmp) in lstinfo_appendix[["to_annotate", "orig_name"]].iterrows():
        shutil.copyfile(old_tmp, new_tmp)
        shutil.copyfile(f"{old_tmp}-prokka.log", f"{new_tmp}-prokka.log")
        shutil.copytree(f"{old_tmp}-prokkaRes", f"{new_tmp}-prokkaRes")

    lstinfo_appendix.to_annotate = lstinfo_appendix.orig_name
    new_lstinfo = pd.concat([new_lstinfo, lstinfo_appendix])
    new_lstinfo.to_csv(new_lstinfo_path, index=False, sep="\t")

    for old_name, new_name in zip(old_names, lstinfo_appendix.gembase_name):
        do_substitution(old_path, old_name, new_path, new_name)


def prepare_data(configargs, available_data, metadata):
    pandir = HOME / f"{configargs.alias}_pan"
    pandir.mkdir(exist_ok=True)
    cur_acs = set(configargs.acs)
    ac_sets = list(map(lambda a: set(a.split(",")), available_data.acs))
    difflist = list(map(lambda dataset_acs: cur_acs - dataset_acs, ac_sets))
    create_dirs(pandir)
    common_diff_set = reduce(lambda a, b: a & b, difflist)
    copied = set()
    for i, acs in enumerate(ac_sets):
        to_copy = (cur_acs & acs) - copied
        copied |= to_copy
        copy_source(Path(available_data.outdir[i]), pandir, configargs.alias, to_copy)
        if copied == cur_acs - common_diff_set:
            print("copied everything drom cache")
            break

    mask = metadata.AC.isin(common_diff_set)
    taxids = pd.unique(metadata[mask].taxid)
    if len(taxids) > 0:
        tmpdir = run_preparation(taxids, configargs.alias)
        tmp_lstinfo = pd.read_table(tmpdir / f"LSTINFO-LSTINFO-{configargs.alias[:4]}.lst", sep="\t")
        name_map = create_name_map(tmp_lstinfo, metadata)
        name_map.to_csv(tmpdir / "names_mapping.tsv", sep='\t', index=False, header=False)
        copy_source(tmpdir, pandir, configargs.alias, common_diff_set)
        shutil.rmtree(tmpdir)

    lstinfo = pd.read_table(pandir / f"LSTINFO-LSTINFO-{configargs.alias[:4]}.lst",
                            sep="\t")
    name_map = create_name_map(lstinfo, metadata)
    name_map.to_csv(pandir / "names_mapping.tsv", sep='\t', index=False, header=False)


def do_substitution(old_dir, old_name, new_dir, new_name):
    dest_prt = new_dir / "Proteins" / f"{new_name}.prt"
    shutil.copyfile(old_dir / "Proteins" / f"{old_name}.prt", dest_prt)
    dest_gene = new_dir / "Genes" / f"{new_name}.gen"
    shutil.copyfile(old_dir / "Genes" / f"{old_name}.gen", dest_gene)
    sp.run(f"sed -i s|{old_name}|{new_name}|g {dest_gene}".split())
    sp.run(f"sed -i s|{old_name}|{new_name}|g {dest_prt}".split())


def subset_blast_graph(old_path, new_path, names_mapping, evalue):
    names_map_path = new_path.parent.parent / "names_for_blast.tsv"
    names_mapping.to_csv(names_map_path, header=None, sep="\t")
    cmd = f"gawk -v eval_th={evalue} -f {HOME / 'subset_blast-graph.awk'} {names_map_path} {old_path}"
    with open(new_path, "w") as outfile:
        sp.run(cmd.split(), stdout=outfile)

    # os.remove(names_map_path)


def inherit_blast_graph(conf, pan, available_data):
    pandir = Path(pan.outdir)
    cur_acs = set(conf.acs)
    has_cached_data = False
    for idx, (lstinfo_path, acs, evalue) in available_data[["lstinfo_file", "acs", "evalue"]].iterrows():
        ac_set = set(acs.split(","))
        if cur_acs <= ac_set and float(conf.eval) <= float(evalue):
            lstinfo_path = Path(lstinfo_path)
            source_pandir = lstinfo_path.parent
            source_po_dir = next(source_pandir.glob("tmp_proteinortho*"))
            has_cached_data = True
            cache_idx = idx
            break

    if has_cached_data:
        source_blast_graph = source_po_dir / f"{available_data.dataset_name[cache_idx]}.blast-graph"
        our_acs = set(conf.acs)
        source_name_map = pd.read_table(source_pandir / "names_mapping.tsv", sep="\t", header=None, index_col=0)
        new_name_map = pd.read_table(pandir / "names_mapping.tsv", sep="\t", index_col=0, header=None)
        old_new_map = pd.DataFrame({"new": list(new_name_map[1][our_acs])},
                                   index=list(source_name_map[1][our_acs]))

        new_po_dir = f"tmp_proteinortho_{conf.alias[:4]}-All_{pan.po_mode}-search"
        if pan.threads > 1:
            new_po_dir = f"{new_po_dir}-th{pan.threads}"

        new_po_dir = pandir / new_po_dir
        new_po_dir.mkdir(exist_ok=True)
        new_blast_graph = new_po_dir / f"{conf.alias[:4]}.blast-graph"
        new_info = new_po_dir / f"{conf.alias[:4]}.info"
        new_info.touch()

        subset_blast_graph(source_blast_graph, new_blast_graph, old_new_map, conf.eval)


def run_plots(conf):
    plots = []
    if conf.venn:
        sp.run(f"{HOME / 'analyse_pan.R'} {conf.alias} venn".split())
        plots.append(f"Venn_plots/{conf.alias}.png")
    else:
        sp.run(f"{HOME / 'analyse_pan.R'} {conf.alias}".split())

    plots.append(f"Error_plots/{conf.alias}.png")
    plots.append(f"Gord_plots/{conf.alias}.png")
    plots.append(f"Fam_distrs/{conf.alias}.png")
    plots.append(f"Copyn_hists/{conf.alias}.png")
    send_album(plots)


def run_tree(conf):
    pandir = HOME / f"{conf.alias}_pan"
    lstinfo_path = pandir / f"LSTINFO-LSTINFO-{conf.alias[:4]}"
    run_sub(f"corepers -p {next(pandir.glob('PanGenome-*.lst'))} -o {pandir}")
    run_sub(f"align -c {next(pandir.glob('PersGenome_*.lst'))} -l {lstinfo_path} "
            f"-d {pandir} -o {pandir} --threads 8 -n {conf.alias[:4]}")


def run_everything(conf, pan):
    prepare_data(conf, available_data, genome_metadata)
    inherit_blast_graph(conf, pan, available_data)
    pan.threads = 8
    run_sub(" ".join(pan.argv))
    run_plots(conf)
    # if not conf.notree:
    #    run_tree(conf)


exe_parser = arp.ArgumentParser(description="Parse execution info.")
exe_parser.add_argument("--threads", type=int, dest="threads", default=1)

with open("families_for_analysis.list") as configfile:
    exe_args = exe_parser.parse_args(configfile.readline().split())
    filtered_lines = list(filter(lambda line: not line.startswith("#"), configfile))

genome_metadata = pd.read_table("myxococcales_table.tsv", sep="\t")

config_parser = arp.ArgumentParser(description="Parse config args")
config_parser.add_argument("--picsonly", dest="pics", action="store_true")
config_parser.add_argument("--notree", dest="notree", action="store_true")
config_parser.add_argument("--venn", action="store_true", dest="venn")
config_parser.add_argument("--evalue", dest="eval", type=str, default=None)
config_parser.add_argument("--conn", dest="conn", type=str, default=None)
config_parser.add_argument("--purity", dest="purity", type=str, default=None)
config_parser.add_argument("--minspec", dest="minspec", type=float, default=None)

pangenome_parser = arp.ArgumentParser(description="Parse PanACoTA pangenome options", add_help=False)
build_parser(pangenome_parser)

argmap = map(lambda line: parse_task(line, genome_metadata, config_parser, pangenome_parser), filtered_lines)
tasklist = list(argmap)
configargs, panargs = tuple(zip(*tasklist))

available_data = find_data()
for conf, pan in zip(configargs, panargs):
    run_everything(conf, pan)



# env = os.getenv("MYXO_ENV", "MUSKRAT")
# if env == "MUSKRAT":
#     with cf.ThreadPoolExecutor(max_workers=exe_args.threads) as executor:
#         executor.map(run_pipeline, filtered_lines)
#
# if env == "MEATGRINDER":
#     run_pipeline(filtered_lines[int(os.environ["SGE_TASK_ID"]) - 1])
