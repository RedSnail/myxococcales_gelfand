import shutil
import subprocess as sp
from telegram import Bot
from telegram import ParseMode
import pandas as pd
from PanACoTA.bin.run_panacota import parse_arguments
from pathlib import Path
import argparse as arp


def run_sub(subcommand):
    action, args = parse_arguments(subcommand.split())
    action(args)


bot_token = "1489745929:AAGVfh8pJ8NMrjtBCiuiAaeVsp4tjwA4Z4M"
my_id = 424187748
HOME = Path("/home/oleg/myxococcales_gelfand/")

bot = Bot(token=bot_token)


def send_text(text):
    bot.sendMessage(chat_id=my_id, text=f"``{text}``", parse_mode=ParseMode.MARKDOWN)


def send_image(path):
    bot.send_photo(chat_id=my_id, photo=(HOME/path).open("rb"))


def perform_pangenome(families, pandir, alias, genome_metadata):
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
    name_map["species"] = genome_metadata.set_index("AC").organism[ac_order].reset_index(drop=True)
    name_map.to_csv(pandir / "names_mapping.tsv", sep='\t', index=False, header=False)
    run_sub(f"pangenome -l {lstinfo_path} -n {alias[:4]} -d {pandir / 'Proteins'} "
            f"-o {pandir} -m proteinortho --threads 8 -v")
    send_text("proteinortho done")
    run_sub(f"corepers -p {next(pandir.glob('PanGenome-*.lst'))} -o {pandir}")
    run_sub(f"align -c {next(pandir.glob('PersGenome_*.lst'))} -l {lstinfo_path} "
            f"-d {pandir} -o {pandir} --threads 8 -n {alias[:4]}")


def run_pipeline(task):
    if task.startswith("#"):
        pass

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
    parser.add_argument("--picsonly", dest="pics", type=bool, default=False)
    args = parser.parse_args(add_args)

    if not args.pics:
        perform_pangenome(families, pandir, alias, genome_metadata)

    sp.run(f"{HOME/'analyse_pan.R'} {alias}".split())

def do_substitution(row, pandir):
    old = row[0]
    new = row[1]
    dest_prt = pandir/"Proteins"/f"{new}.prt"
    shutil.copyfile(HOME/"WithOut_pan"/"Proteins"/f"{old}.prt", dest_prt)
    dest_gene = pandir/"Genes"/f"{new}.gen"
    shutil.copyfile(HOME/"WithOut_pan"/"Genes"/f"{old}.gen", dest_gene)
    sp.run(f"sed -i s|{old}|{new}|g {dest_gene}".split())
    sp.run(f"sed -i s|{old}|{new}|g {dest_prt}".split())


# with open("families_for_analysis.list") as configfile:
#   for line in configfile:
#        run_pipeline(line)

run_pipeline("Testdrive Archangiaceae")




