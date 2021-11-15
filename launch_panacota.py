import shutil
import subprocess as sp
from telegram import Bot, InputMediaPhoto
from telegram import ParseMode
import pandas as pd
from PanACoTA.bin.run_panacota import parse_arguments
from pathlib import Path
import argparse as arp
import concurrent.futures as cf


def run_sub(subcommand):
    action, args = parse_arguments(subcommand.split())
    action(args)


bot_token = "1489745929:AAGVfh8pJ8NMrjtBCiuiAaeVsp4tjwA4Z4M"
my_id = 424187748
HOME = Path("/home/oleg/myxococcales_gelfand/")

bot = Bot(token=bot_token)


def send_album(img_paths):
    photos = list(map(lambda path: InputMediaPhoto(open(path, "rb")), img_paths))
    bot.send_media_group(my_id, photos)


def send_text(text):
    bot.sendMessage(chat_id=my_id, text=f"``{text}``", parse_mode=ParseMode.MARKDOWN)


def send_image(path):
    bot.send_photo(chat_id=my_id, photo=(HOME / path).open("rb"))


def perform_pangenome(families, pandir, alias, genome_metadata, notree, eval, conn, purity, minspec):
    acs = pd.unique(genome_metadata[genome_metadata.family.isin(families)]["AC"])
    lstinfo = pd.read_table((HOME / "WithOut_pan" / "LSTINFO-LSTINFO-With.lst"), sep="\t")
    mask = lstinfo.orig_name.str.contains("|".join(acs))
    lstinfo = lstinfo[mask]
    lstinfo.reset_index(drop=True, inplace=True)

    name_map = pd.DataFrame({"old": lstinfo.gembase_name,
                             "new": map(lambda i: f"{alias[:4]}.0921.{i + 1}", range(lstinfo.shape[0]))})

    lstinfo.gembase_name = name_map.new
    lstinfo_path = pandir / f"LSTINFO-LSTINFO-{alias[:4]}.lst"
    lstinfo.to_csv(lstinfo_path, sep='\t', index=False)

    name_map.apply(do_substitution, axis=1, args=(pandir,))
    ac_order = lstinfo.orig_name.str.extract("(" + "|".join(acs) + ")")[0]
    ac_map = pd.DataFrame({"AC": ac_order,
                           "name": name_map.new,
                           "species": genome_metadata.set_index("AC").organism[ac_order].reset_index(drop=True)})
    ac_map.to_csv(pandir / "names_mapping.tsv", sep='\t', index=False, header=False)
    pangenome_command = (f"pangenome -l {lstinfo_path} -n {alias[:4]} -d {pandir / 'Proteins'} "
                         f"-o {pandir} -m proteinortho --threads 8 -v ")

    if eval is not None:
        pangenome_command += f"--eval {eval} "

    if conn is not None:
        pangenome_command += f"--conn {conn} "

    if purity is not None:
        pangenome_command += f"--purity {purity}"

    if minspec is not None:
        pangenome_command += f"--minspecies {minspec} "

    run_sub(pangenome_command)
    send_text("proteinortho done")
    if not notree:
        run_sub(f"corepers -p {next(pandir.glob('PanGenome-*.lst'))} -o {pandir}")
        run_sub(f"align -c {next(pandir.glob('PersGenome_*.lst'))} -l {lstinfo_path} "
                f"-d {pandir} -o {pandir} --threads 8 -n {alias[:4]}")


def run_pipeline(task):
    print(task)
    alias = task.split()[0]
    pandir = HOME / f"{alias}_pan"
    pandir.mkdir(exist_ok=True)
    (pandir / "Proteins").mkdir(exist_ok=True)
    (pandir / "Genes").mkdir(exist_ok=True)
    genome_metadata = pd.read_table("myxococcales_table.tsv", sep="\t")
    no_cow = ~genome_metadata.organism.eq("Pajaroellobacter abortibovis")
    genome_metadata = genome_metadata[no_cow]
    if task.startswith("WithOut"):
        taxids = pd.unique(genome_metadata["taxid"])
        run_sub(f"prepare -t {','.join(str(taxid) for taxid in taxids)} -l complete,chromosome -o "
                f"~/myxococcales_gelfand/WithOut_pan -p 4 --min_dist 0 --max_dist 1")
        run_sub("annotate  --info WithOut_pan/LSTINFO-WithOut.txt -r WithOut_pan --threads 4 -n With")

        families = pd.unique(genome_metadata["family"])

        add_args = task.split()[1:]
    elif task.startswith("Bulk"):
        mask = ~genome_metadata.family.eq("Geobacteraceae")
        families = pd.unique(genome_metadata[mask]["family"])
        add_args = task.split()[1:]
    else:
        families = task.split()[1].split(",")
        add_args = task.split()[2:]

    parser = arp.ArgumentParser(description="Parse config args")
    parser.add_argument("--picsonly", dest="pics", action="store_true")
    parser.add_argument("--notree", dest="notree", action="store_true")
    parser.add_argument("--venn", action="store_true", dest="venn")
    parser.add_argument("--evalue", dest="eval", type=str, default=None)
    parser.add_argument("--conn", dest="conn", type=str, default=None)
    parser.add_argument("--purity", dest="purity", type=str, default=None)
    parser.add_argument("--minspec", dest="minspec", type=float, default=None)
    args = parser.parse_args(add_args)

    if not args.pics:
        perform_pangenome(families, pandir, alias, genome_metadata, args.notree,
                          args.eval, args.conn, args.purity, args.minspec)
    
    plots = []
    if args.venn:
        sp.run(f"{HOME / 'analyse_pan.R'} {alias} venn".split())
        plots.append(f"Venn_plots/{alias}.png")
    else:
        pass
        sp.run(f"{HOME / 'analyse_pan.R'} {alias}".split())

    plots.append(f"Error_plots/{alias}.png")
    plots.append(f"Gord_plots/{alias}.png")
    plots.append(f"Fam_distrs/{alias}.png")
    plots.append(f"Copyn_hists/{alias}.png")

    send_album(plots)


def do_substitution(row, pandir):
    old = row[0]
    new = row[1]
    dest_prt = pandir / "Proteins" / f"{new}.prt"
    shutil.copyfile(HOME / "WithOut_pan" / "Proteins" / f"{old}.prt", dest_prt)
    dest_gene = pandir / "Genes" / f"{new}.gen"
    shutil.copyfile(HOME / "WithOut_pan" / "Genes" / f"{old}.gen", dest_gene)
    sp.run(f"sed -i s|{old}|{new}|g {dest_gene}".split())
    sp.run(f"sed -i s|{old}|{new}|g {dest_prt}".split())


with open("families_for_analysis.list") as configfile:
    exe_parser = arp.ArgumentParser(description="Parse execution info.")
    exe_parser.add_argument("--threads", type=int, dest="threads", default=1)
    exe_args = exe_parser.parse_args(configfile.readline().split())
    active_th = 0
    filtered_lines = filter(lambda line: not line.startswith("#"), configfile)
    with cf.ThreadPoolExecutor(max_workers=exe_args.threads) as executor:
        executor.map(run_pipeline, list(filtered_lines))



