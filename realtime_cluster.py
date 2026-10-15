#!/usr/bin/env python

"""
Real-time pipeline for Nanopore runs

Each cycle: it checks for new files in fastq_pass/barcodeXX and maps ONLY the new ones.
Then it merges the BAM files with the cumulative one (suffix -N = processed files).
Variant calling > annotates positions > HTML report.

The script stops if '--idle_timeout' seconds (1 hour by default) pass without new files

Usage:
	realtime_prokaion.py -i Run_folder -s barcodes_summary.tsv \
		-r reference.fasta -c annotation.tsv -o Results -t 32
"""


import argparse
import datetime
import glob
import gzip
import html
import json
import logging
import os
import sys
import re
import subprocess
import shutil
import time
import traceback
import pandas as pd
import multiprocessing


logger = logging.getLogger("Real-Time prokaion")
tag = re.compile(r"(Marker|Branch)-([\w.]+)")
min_markers = 10


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
RED = '\033[31m'
GREEN = '\033[32m'
DIM = "\033[2m"


def check_create_dir(path):
    # exists = os.path.isfile(path)
    # exists = os.path.isdir(path)

    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)


def execute_subprocess(cmd, isShell=False, isInfo=False):
    """
    https://crashcourse.housegordon.org/python-subprocess.html
    https://docs.python.org/3/library/subprocess.html 
    Execute and handle errors with subprocess, outputting stderr instead of the subprocess CalledProcessError
    """

    logger.debug("")
    logger.debug(cmd)

    if cmd[0] == "samtools" or cmd[0] == "bwa":
        prog = " ".join(cmd[0:2])
        param = cmd[3:]
    else:
        prog = cmd[0]
        param = cmd[1:]

    try:
        command = subprocess.run(
            cmd, shell=isShell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if command.returncode == 0:
            logger.debug(
                GREEN + DIM + "Program %s successfully executed" % prog + END_FORMATTING)
        else:
            logger.info(RED + BOLD + "Command %s FAILED\n" % prog + END_FORMATTING + BOLD + "with parameters: " + END_FORMATTING + " ".join(
                param) + "\n" + BOLD + "EXIT-CODE: %d\n" % command.returncode + "ERROR:\n" + END_FORMATTING + command.stderr.decode().strip())

        if isInfo:
            logger.info(command.stdout.decode().strip())
        else:
            logger.debug(command.stdout.decode().strip())

        logger.debug(command.stderr.decode().strip())

    except OSError as e:
        sys.exit(RED + BOLD + "Failed to execute program '%s': %s" % (prog,
                                                                      str(e)) + END_FORMATTING)


### Input from the Nanopore run folder ###

def read_samples(path):
	"""barcode	samplename"""
	out = []
	with open(path) as f:
		for line in f:
			if not line.strip() or line.startswith("#"):
				continue
			barcode, sample = line.rstrip("\n").split("\t")[:2]
			barcode = barcode.strip()
			if barcode.isdigit():
				barcode = "barcode{:02d}".format(int(barcode))
			out.append((barcode, sample.strip()))
		return out


def find_barcode_dir(run_dir, barcode):
	hits = glob.glob(os.path.join(run_dir, "**", "fastq_pass", barcode), recursive=True)
	return hits[0] if hits else None


def ready_files(bdir, settle):
	"""Fastq files that are already closed (unchanged for the last settle)"""
	limit = time.time() - settle
	files = [f for f in glob.glob(os. path.join(bdir, "*"))
			 if f.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")) and os.path.getmtime(f) < limit]
	return sorted(files)


def concat(files, out):
    with open(out, "wb") as fo:
        for f in files:
            if f.endswith(".gz"):
                with open(f, "rb") as fi:
                    shutil.copyfileobj(fi, fo)
            else:
                with open(f, "rb") as fi, gzip.GzipFile(fileobj=fo, mode="wb") as gz:
                    shutil.copyfileobj(fi, gz)


### ProkaION preexisting functions ###

def extract_indels(input_vcf):

    input_vcf = os.path.abspath(input_vcf)
    vcf_dir = ("/").join(input_vcf.split("/")[0:-1])
    output_indel_vcf = os.path.join(vcf_dir, "snps.indel.vcf")

    with open(output_indel_vcf, "w+") as fout:
        with open(input_vcf, "r") as f:
            for line in f:
                if "TYPE=ins" in line or "TYPE=del" in line:
                    fout.write(line)


def merge_vcf(snp_vcf, indel_vcf):

    snp_vcf = os.path.abspath(snp_vcf)
    indel_vcf = os.path.abspath(indel_vcf)

    vcf_dir = ("/").join(snp_vcf.split("/")[0:-1])
    output_complete_vcf = os.path.join(vcf_dir, "snps.all.vcf")

    with open(output_complete_vcf, "w+") as fout:
        with open(snp_vcf, "r") as f1:
            for line in f1:
                fout.write(line)
        with open(indel_vcf, "r") as f2:
            for line in f2:
                if not line.startswith("#"):
                    fout.write(line)


def import_VCF42_freebayes_to_tsv(vcf_file, sep="\t"):

    vcf_file = os.path.abspath(vcf_file)
    tsv_file = (".").join(vcf_file.split(".")[:-1]) + ".tsv"

    headers = []
    extra_fields = ["TYPE", "DP", "RO", "AO"]
    with open(tsv_file, "w+") as fout:
        with open(vcf_file, "r") as f:
            next_line = f.readline().strip()
            while next_line.startswith("#"):
                next_line = f.readline().strip()
                if next_line.startswith("#CHROM"):
                    headers = next_line.split("\t")

        headers = headers[:7] + extra_fields + ["OLDVAR"]
        fout.write(("\t").join(headers) + "\n")

        with open(vcf_file, "r") as f:
            for line in f:
                extra_field_list = []
                # and not 'complex' in line and not 'mnp' in line
                if not line.startswith("#"):
                    line_split = line.split(sep)[:8]
                    info = line_split[-1].split(";")
                    for field in extra_fields:
                        extra_field_list.append(
                            [x.split("=")[-1] for x in info if field in x][0])
                    if "OLDVAR" in line:
                        extra_field_list.append(
                            [x.split("=")[-1] for x in info if "OLDVAR" in x][0].split(",")[0])
                    output_line = ("\t").join(
                        line_split[:7] + extra_field_list)
                    fout.write(output_line + "\n")


def import_tsv_freebayes_to_df(tsv_file, sep="\t"):

    tsv_file = os.path.abspath(tsv_file)

    df = pd.read_csv(tsv_file, sep=sep)

    df.rename(columns={"#CHROM": "REGION", "RO": "REF_DP",
                       "DP": "TOTAL_DP", "AO": "ALT_DP", "QUAL": "ALT_QUAL"}, inplace=True)

    df["REF_FREQ"] = df["REF_DP"] / df["TOTAL_DP"]
    df["ALT_FREQ"] = df["ALT_DP"] / df["TOTAL_DP"]

    df = df.sort_values(by=["POS"]).reset_index(drop=True)

    return df[["REGION", "POS", "ID", "REF", "ALT", "ALT_QUAL", "FILTER", "TOTAL_DP", "TYPE", "REF_DP", "ALT_DP", "REF_FREQ", "ALT_FREQ", "OLDVAR"]]


def vcf_to_ivar_tsv(input_vcf, output_tsv):

    input_tsv = (".").join(input_vcf.split(".")[:-1]) + ".tsv"

    import_VCF42_freebayes_to_tsv(input_vcf)

    df = import_tsv_freebayes_to_df(input_tsv)
    df.to_csv(output_tsv, sep="\t", index=False)


def minimap2_mapping(HQ_filename, filename_bam_out, reference):
    """
    https://github.com/lh3/minimap2
        # Oxford Nanopore genomic reads
        minimap2 -ax map-ont ref.fa ont.fq.gz > aln.sam
    http://www.htslib.org/doc/samtools.html
    """

    # -a: Output in the SAM format
    # -x: Preset (always applied before other options; see minimap2.1 for details) []
    #    - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
    #    - map-hifi - PacBio HiFi reads vs reference mapping
    #    - ava-pb/ava-ont - PacBio/Nanopore read overlap
    #    - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    #    - splice/splice:hq - long-read/Pacbio-CCS spliced alignment
    #    - sr - genomic short-read mapping
    # -t: Number of threads

    # -b: Output BAM
    # -S: Ignored (input format is auto-detected)
    # -F: Only include reads with none of the FLAGS in INT present
    # --threads: Number of additional threads to use

    cmd_minimap2 = "minimap2 -ax map-ont {} {} | samtools view -bS -F 4 - | samtools sort -o {}".format(
        reference, HQ_filename, filename_bam_out)
    # print(cmd_minimap2)
    execute_subprocess(cmd_minimap2, isShell=True)

    cmd_indexing = "samtools", "index", filename_bam_out
    # print(cmd_indexing)
    execute_subprocess(cmd_indexing, isShell=False)


def ngmlr_mapping(HQ_filename, filename_bam_out, reference, threads=30):
    """
    https://github.com/philres/ngmlr
    """

    # -r <file>, --reference <file>
    # -q <file>, --query <file>
    # -o <string>, --output <string>
    # -t <int>, --threads <int>
    # -x <pacbio, ont>, --presets <pacbio, ont>
    # -i <0-1>, --min-identity <0-1> - Alignments with an identity lower than this threshold will be discarded [0.65]

    cmd_ngmlr = "ngmlr -t {} -r {} -q {} -x ont | samtools view -bS -F 4 - | samtools sort -o {}".format(
        threads, reference, HQ_filename, filename_bam_out)
    # print(cmd_ngmlr)
    execute_subprocess(cmd_ngmlr, isShell=True)

    cmd_indexing = "samtools", "index", filename_bam_out
    # print(cmd_indexing)
    execute_subprocess(cmd_indexing, isShell=False)


def optimal_chunk_size(reference):
    """
    https://www.bioinformatics.nl/emboss-explorer/
        http://www.sacs.ucsf.edu/Documentation/emboss/infoseq.html
    """

    input_reference = os.path.abspath(reference)
    input_folder = os.path.dirname(reference)
    chunks_size = "chunks_size"
    chunks_size_path = os.path.join(input_folder, "chunks_size")

    cmd_infoseq = ["infoseq", "-auto", "-only", "-length", "-noheading",
                   "-odirectory", input_folder, "-outfile", chunks_size, input_reference]
    # print(cmd_infoseq)
    execute_subprocess(cmd_infoseq, isShell=False)

    df = pd.read_csv(chunks_size_path, header=None)
    size = df.loc[df.index[0]]
    return size


def create_reference_chunks(reference):

    ref_size = optimal_chunk_size(reference)
    min_freebayes_chunk_size = 1000
    chunk_cpu = multiprocessing.cpu_count() - 2
    num_chunks = round(max(min_freebayes_chunk_size, int(ref_size.iloc[0]) / chunk_cpu))

    input_reference = os.path.abspath(reference)
    input_folder = os.path.dirname(reference)
    out_reference_file = os.path.join(
        input_folder, "reference." + str(num_chunks) + ".regions")
    fai_reference = input_reference + ".fai"

    if os.path.isfile(out_reference_file):
        # logger.info(out_reference_file + " already EXISTS")
        pass
    else:
        logger.info(GREEN + "Creating " + out_reference_file + END_FORMATTING)
        cmd_chunks = "fasta_generate_regions.py {} {} > {}".format(
            fai_reference, num_chunks, out_reference_file)
        execute_subprocess(cmd_chunks, isShell=True)

    return out_reference_file


def freebayes_variant(reference, filename_bam_out, output_vcf, threads=36, frequency=0.1, ploidy=1, base_qual=7, map_qual=60):
    """
    https://github.com/freebayes/freebayes
        # Freebayes-parallel
        freebayes-parallel <(fasta_generate_regions.py {fai_reference} {chunks}) {threads} {args} > {output}
    """

    # --region_file: Genome partitioning according to computational threads.
    # --haplotype_length: Allow haplotype calls with contiguous embedded matches of up to this length
    # --use-best-n-alleles: Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores
    # --min-alternate-count: Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position
    # --min-alternate-fraction: Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position
    # -p: Sets the default ploidy for the analysis to N
    # --min-coverage:
    # -q: Exclude alleles from analysis if their supporting base quality is less than Q
    # -m: Exclude alignments from analysis if they have a mapping quality less than Q
    # --strict-vcf: Generate strict VCF format (FORMAT/GQ will be an int)

    region_file = create_reference_chunks(reference)

    cmd_bayes = "freebayes-parallel {} {} -f {} --haplotype-length 0 --use-best-n-alleles 1 --min-alternate-count 0 --min-alternate-fraction {} -p {} --min-coverage 1 -q {} -m {} --strict-vcf {} > {}".format(
        region_file, str(threads), reference, str(frequency), str(ploidy), str(base_qual), str(map_qual), filename_bam_out, output_vcf)
    # print(cmd_bayes)
    execute_subprocess(cmd_bayes, isShell=True)


def freebayes(reference, filename_bam_out, output_vcf, frequency=0.1, ploidy=1, base_qual=7, map_qual=60):

    cmd_bayes = "freebayes -f {} --haplotype-length 0 --use-best-n-alleles 1 --min-alternate-count 0 --min-alternate-fraction {} -p {} --min-coverage 1 -q {} -m {} --strict-vcf {} > {}".format(
        reference, str(frequency), str(ploidy), str(base_qual), str(map_qual), filename_bam_out, output_vcf)
    # print(cmd_bayes)
    execute_subprocess(cmd_bayes, isShell=True)


def bcftool_filter(output_raw_vcf, output_vcf):
    """
    https://samtools.github.io/bcftools/bcftools.html
        # bcftools view: View, subset and filter VCF files by position and filtering expression
        # bcftools annotate: Add or remove annotations
    """

    # --include: include sites for which Expression is true
    # --remove: list of annotations to remove

    # cmd_bcf_view = "bcftools view --include 'QUAL >= 10 && FMT/DP >= 1 && (FMT/AO)/(FMT/DP) >= 0.1' {} -o {}".format(output_raw_vcf, output_vcf)
    # # print(cmd_bcf_view)
    # execute_subprocess(cmd_bcf_view, isShell=True)

    # cmd_bcf_annot = 'bcftools annotate --remove ^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL {} -o {}'.format(output_vcf, output_vcf_filt)
    # # print(cmd_bcf_annot)
    # execute_subprocess(cmd_bcf_annot, isShell=True)

    cmd_bcf = "bcftools view --include '(GT=\"1/1\" || GT=\"0/1\" || GT=\"1\") && QUAL >= 10 && FMT/DP >= 1 && (FMT/AO)/(FMT/DP) >= 0.1' {} | bcftools annotate --remove ^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL -o {}".format(
        output_raw_vcf, output_vcf)
    # print(cmd_bcf)
    execute_subprocess(cmd_bcf, isShell=True)


def snippy_sub(output_vcf_filt, output_vcf_sub):
    """
    https://github.com/tseemann/snippy/tree/master/bin
        # snippy-vcf_extract_subs: Convert MNP,COMPLEX into SNP and ignore INS,DEL
    """

    cmd_snippy_subs = "snippy-vcf_extract_subs {} > {}".format(
        output_vcf_filt, output_vcf_sub)
    # print(cmd_snippy_subs)
    execute_subprocess(cmd_snippy_subs, isShell=True)


### Cluster annotation steps ###

def load_info(path):
	"""File with header: #CHROM POS ID REF ALT QUAL FILTER FORMAT INFO"""
	df = pd.read_csv(path, sep="\t", dtype=str).rename(columns={"#CHROM": "REGION"})
	rows = []
	for r in df.itertuples():
		for kind, cl in tag.findall(str(r.INFO)):
			rows.append((r.REGION, int(r.POS), r.REF, r.ALT, kind, cl))
	return pd.DataFrame(rows, columns=["REGION", "POS", "REF", "ALT", "KIND", "INFO"])


def annotate(tsv, clusters, out_file):
	v = pd.read_csv(tsv, sep="\t")
	tags = (clusters.assign(CLUSTER_INFO=clusters.KIND + "-" + clusters.INFO).groupby(["POS", "REF", "ALT"], as_index=False).CLUSTER_INFO.agg(",".join))
	v.merge(tags, how="left", on=["POS", "REF", "ALT"]).to_csv(out_file, sep="\t", index=False)
	return v


def score(variants, clusters):
	"""Per cluster: X SNPs from Y markers, x SNPs from Z branch markers
	Assigned if x >= min(10, Y): >10 markers quite 10; <=10 require all"""
	hit = clusters.merge(variants[["POS", "REF", "ALT"]].drop_duplicates(), on=["POS", "REF", "ALT"])
	tot = clusters.groupby(["INFO", "KIND"]).size().unstack(fill_value=0)
	got = hit.groupby(["INFO", "KIND"]).size().unstack(fill_value=0)
	res = []
	for c in tot.index:
		Y, Z = int(tot.loc[c].get("Marker", 0)), int(tot.loc[c].get("Branch", 0))
		X = int(got.loc[c].get("Marker", 0)) if c in got.index else 0
		x = int(got.loc[c].get("Branch", 0)) if c in got.index else 0
		need = min(min_markers, Y)
		res.append(dict(cluster=c, X=X, Y=Y, x=x, Z=Z, need=need, assigned=bool(Y and X >= need)))
	res.sort(key=lambda r: (not r["assigned"], -(r["X"] / r["Y"] if r["Y"] else 0), -r["x"]))
	return res


def make_summary(res):
    ok = [r for r in res if r["assigned"]]
    if ok:
        return ", ".join("Cluster {} ({}/{} markers, {}/{} branch markers)".format(
            r["cluster"], r["X"], r["Y"], r["x"], r["Z"]) for r in ok)
    return "Unassigned"


### Analysis ###

def coverage(bam, clusters, assigned, wd, args):
	"""Positions of each assigned cluster with a depth > threshold, read from the BAM file"""
	sub = clusters[clusters["INFO"].isin(assigned)][["REGION", "POS"]].drop_duplicates().sort_values(["REGION", "POS"])
	bed, dp = os.path.join(wd, "cluster_positions.bed"), os.path.join(wd, "cluster_positions.depth")
	sub.assign(start=sub.POS - 1)[["REGION", "start", "POS"]].to_csv(bed, sep="\t", header=False, index=False)
	execute_subprocess("samtools depth -a -b {} -q {} -Q {} {} > {}".format(
		bed, args.min_quality, args.min_mapping, bam, dp), isShell=True)
	if os.path.getsize(dp) == 0:
		raise RuntimeError("samtools depth did not return any positions: Does the REGION name match?")
	d = pd.read_csv(dp, sep="\t", header=None, names=["REGION", "POS", "DEPTH"])
	ok = set(zip(d.REGION[d.DEPTH > args.annot_min_dp], d.POS[d.DEPTH > args.annot_min_dp]))
	out = {}
	for c in assigned:
		s = clusters[clusters["INFO"] == c]
		out[c] = {k: sum((r, p) in ok for r, p in set(zip(g.REGION, g.POS))) for k, g in s.groupby("KIND")}
	return out


def analyse(sample, new, prev, args, clusters):
	"new: new fastq files; prev: files already analysed until now"
	n = prev.get("n", 0) + len(new)
	names = "{}-{}".format(sample, n)
	sdir = os.path.join(args.output, "samples", sample)
	cur, bat, wd = (os.path.join(sdir, d) for d in ("current", "batch", names))
	for d in (cur, bat, wd):
		os.makedirs(d, exist_ok=True)
	bam, fq = os.path.join(cur, names + ".bam"), os.path.join(cur, names + ".fastq.gz")
	p = lambda name: os.path.join(wd, "variants", name)

	if not os.path.isfile(bam + ".bai"):
		nfq, nbam = (os.path.join(bat, f) for f in ("new.fastq.gz", "new.bam"))
		logger.info("[%s] processing %d new files", sample, len(new))
		concat(new, nfq)
		if args.amplicon:
			ngmlr_mapping(nfq, nbam, args.reference, threads=args.threads)
		else:
			minimap2_mapping(nfq, nbam, reference=args.reference)

		if prev.get("n"):
			if not os.path.isfile(prev["bam"]):
				raise RuntimeError("The cumulative BAM " + prev["bam"] + " is missing")
			execute_subprocess(("samtools", "merge", "-f", "-@", str(args.threads), bam, prev["bam"], nbam), isShell=False)
			execute_subprocess(("samtools", "index", bam), isShell=False)

			if os.path.isfile(prev["fq"]):
				os.replace(prev["fq"], fq)
			with open(fq, "ab") as fo, open (nfq, "rb") as fi:
				shutil.copyfileobj(fi, fo)
			for f in (prev["bam"], prev["bam"] + ".bai"):
				if os.path.isfile(f):
					os.remove(f)
		else:
			os.replace(nbam, bam)
			os.replace(nbam + ".bai", bam + ".bai")
			os.replace(nfq, fq)
		shutil.rmtree(bat, ignore_errors=True)

	os.makedirs(os.path.join(wd, "variants"), exist_ok=True)

	if not os.path.isfile(p("snps.raw.vcf")):
		# logger.info("[%s] variant calling", names)
		kw = dict(frequency=args.min_allele_frequency, ploidy=args.ploidy, base_qual=args.min_quality, map_qual=args.min_mapping)
		if args.amplicon:
			freebayes(args.reference, bam, p("snps.raw.vcf"), **kw)
		else:
			freebayes_variant(args.reference, bam, p("snps.raw.vcf"), threads=args.threads, **kw)
	if not os.path.isfile(p("snps.subs.vcf")):
		bcftool_filter(p("snps.raw.vcf"), p("snps.vcf"))
		snippy_sub(p("snps.vcf"), p("snps.subs.vcf"))
	if not os.path.isfile(p("snps.indel.vcf")):
		extract_indels(p("snps.vcf"))
	if not os.path.isfile(p("snps.all.vcf")):
		merge_vcf(p("snps.subs.vcf"), p("snps.indel.vcf"))
	if not os.path.isfile(p("snps.all.ivar.tsv")):
		vcf_to_ivar_tsv(p("snps.all.vcf"), p("snps.all.ivar.tsv"))

	# logger.info("[%s] clusters anotation", names)
	variants = annotate(p("snps.all.ivar.tsv"), clusters, os.path.join(wd, names + ".annotated.tsv"))
	passing = variants[(variants["ALT_FREQ"] >= args.annot_min_freq) & (variants["TOTAL_DP"] > args.annot_min_dp)]
	res = score(variants[(variants["ALT_FREQ"] >= args.annot_min_freq) & (variants["TOTAL_DP"] > args.annot_min_dp)], clusters)
	assigned = [r["cluster"] for r in res if r["assigned"]]

	cov = coverage(bam, clusters, assigned, wd, args) if assigned else {}
	for r in res:
		if r["cluster"] in cov:
			r["Cm"] = cov[r["cluster"]].get("Marker", 0)
			r["Cb"] = cov[r["cluster"]].get("Branch", 0)
			r["mindp"] = args.annot_min_dp

	summary = make_summary(res)
	return dict(n=n, time=datetime.datetime.now().strftime("%H:%M:%S"), res=res, assigned=assigned, summary=summary, bam=bam, fq=fq)



CSS = """
:root{color-scheme:light;--bg:#F2F5F6;--card:#fff;--ink:#16262B;--mute:#5B6E74;--marker:#0E6B73;--branch:#B4691A;--ok:#1F7A4D;--line:#D5DEE1;--track:#E3EAEC}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:860px;margin:0 auto;padding:32px 20px 64px}
h1{font:600 clamp(1.7rem,4.5vw,2.6rem)/1.15 Georgia,"Iowan Old Style",serif;margin:0 0 8px}
.sub{display:block;margin-top:4px;font:500 .55em system-ui,-apple-system,"Segoe UI",sans-serif;color:var(--mute)}
h2{font:600 1.15rem Georgia,serif;margin:36px 0 12px}
.lead{color:var(--mute);margin:0 0 28px;max-width:62ch}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:18px 20px;margin-bottom:14px}
.card.assigned{border-left:6px solid var(--ok)}
.card header{display:flex;justify-content:space-between;align-items:baseline;gap:12px;flex-wrap:wrap}
.card h3{margin:0;font-size:1.1rem}
.badge{font-size:.85rem;color:var(--ok);font-weight:600}.partial .badge{color:var(--mute)}
.row{margin-top:14px}.row p{margin:0 0 6px;font-size:.95rem}
.bar{position:relative;height:14px;background:var(--track);border-radius:3px}
.bar span{display:block;height:100%;border-radius:3px}
.marker span{background:var(--marker)}.branch span{background:var(--branch)}
.tick{position:absolute;top:-5px;bottom:-5px;width:2px;background:var(--ink)}
.hint{color:var(--mute);font-size:.85rem;margin-top:4px}
table{width:100%;border-collapse:collapse;background:var(--card);border:1px solid var(--line);border-radius:10px;overflow:hidden}
th,td{text-align:left;padding:9px 14px;border-bottom:1px solid var(--line);font-size:.95rem}
th{color:var(--mute);font-weight:600}a{color:var(--marker)}
tr.ok td{background:#DDF3E4}
tr.no td{background:#F8DEDC}
@media (max-width:520px){main{padding:20px 14px}}
"""


def esc(x):
	return html.escape(str(x))


def page(title, body):
	return ('<!doctype html><html lang="es"><head><meta charset="utf-8">'
			'<meta name="viewport" content="width=device-width,initial-scale=1">'
			'<meta http-equiv="refresh" content="60"><title>{}</title><style>{}</style></head>'
			'<body><main>{}</main></body></html>').format(esc(title), CSS, body)


def bar(val, total, cls, need=None):
	pct = 100.0 * val / total if total else 0
	tick = '<i class="tick" style="left:{:.1f}%"></i>'.format(100.0 * need / total) if (need and total) else ""
	return '<div class="bar {}"><span style="width:{:.1f}%"></span>{}</div>'.format(cls, pct, tick)


def marker_row(called, covered, total, name, cls):
	"""X called of Y covered (Z in total); bar = called / total"""
	if covered is None:   # iteraciones antiguas guardadas sin cobertura
		txt = "<b>{}</b> total {}: <b>{}</b> markers called".format(total, name, called)
	else:
		txt = "<b>{}</b> total {}: <b>{}</b> markers called of <b>{}</b> covered".format(total, name, called, covered)
	return '<div class="row"><p>{}</p>{}</div>'.format(txt, bar(called, total, cls))


def cov_html(r):
	hint = '<div class="hint">Covered = depth &gt; {}x</div>'.format(r["mindp"]) if "mindp" in r else ""
	return (marker_row(r["X"], r.get("Cm"), r["Y"], "marker SNPs", "marker") +
			marker_row(r["x"], r.get("Cb"), r["Z"], "branch marker SNPs", "branch") + hint)


def card(r):
	badge = "Assigned" + (", with branch" if r["x"] else "")
	return ('<section class="card assigned"><header><h3>Cluster {}</h3>'
			'<span class="badge">{}</span></header>{}</section>'
			).format(esc(r["cluster"]), badge, cov_html(r))


def sample_page(sample, entry, history):
	res = entry["res"]
	assigned = [r for r in res if r["assigned"]]
	partial = [r for r in res if not r["assigned"] and (r["X"] or r["x"])]
	if assigned:
		# head = '{}<br><span class="sub">Cluster {} candidate</span>'.format(esc(sample), esc(", ".join(entry["assigned"])))
		# lead = " ".join("{} SNPs of {} markers, and {} SNPs of {} branch markers".format(r["X"], r["Y"], r["x"], r["Z"]) for r in assigned)
		head = '{}'.format(esc(sample))
		lead = 'Cluster {} candidate'.format(esc(", ".join(entry["assigned"])))
	else:
		head = '{}'.format(esc(sample))
		lead = "There are no new entries into preexisting clusters yet"
	body = "<h1>{}</h1><p class=\"lead\">{}</p><p class=\"hint\">Analysis with {} fastq files, updated at {}.</p>".format(head, lead, entry["n"], entry["time"])
	# body += "".join(card(r) for r in assigned + partial[:5])
	body += "".join(card(r) for r in assigned) 
	rows = "".join("<tr><td>{}</td><td>{}</td><td>{}</td></tr>".format(h["n"], h["time"], esc(h["summary"]))
				   for h in reversed(history))
	body += "<h2>Evolution during the run</h2><table><tr><th>Files</th><th>Time</th><th>Cluster candidate</th></tr>{}</table>".format(rows)
	body += '<p class="hint"><a href="Summary.html">Return to the run summary</a></p>'
	return page("Report " + sample, body)


def index_result(e):
	ok = [r for r in e["res"] if r["assigned"]]
	if not ok:
		return "Unassigned"
	return ", ".join("Cluster {}".format(r["cluster"]) for r in ok)


def index_page(state):
	rows = ""
	for s in sorted(state):
		e = state[s]["history"][-1]
		cls = "ok" if e["assigned"] else "no"
		rows += '<tr class="{0}"><td><a href="{1}.html">{1}</a></td><td>{2}</td><td>{3}</td></tr>'.format(cls, esc(s), esc(index_result(e)), e["time"])
	body = ("<h1>Run summary</h1><p class=\"lead\">{} sequences in analysis. The page refreshes automatically every minute</p>"
			"<table><tr><th>Sequences</th><th>Cluster candidate</th><th>Time</th></tr>{}</table>").format(len(state), rows)
	return page("Run summary", body)
 
 
def write_reports(args, state):
	rdir = os.path.join(args.output, "reports")
	check_create_dir(rdir)
	for s, st in state.items():
		with open(os.path.join(rdir, s + ".html"), "w") as f:
			f.write(sample_page(s, st["history"][-1], st["history"]))
	with open(os.path.join(rdir, "Summary.html"), "w") as f:
		f.write(index_page(state))


def save_state(path, state):
	tmp = path + ".tmp"
	with open(tmp, "w") as f:
		json.dump(state, f)
	os.replace(tmp, path)



def main():
	ap = argparse.ArgumentParser(description="Real-Time analysis of a Nanopore run")
	ap.add_argument("-i", "--input", required=True, help="Run folder (read-onlt)")
	ap.add_argument("-s", "--samples", required=True, help="TSV: barcode<TAB>muestra")
	ap.add_argument("-r", "--reference", required=True)
	ap.add_argument("-c", "--clusters", required=True, help="File with info Marker-X / Branch-X")
	ap.add_argument("-o", "--output", required=True, help="Output folder with every output")
	ap.add_argument("-t", "--threads", type=int, default=8)
	ap.add_argument("--amplicon", action="store_true")
	ap.add_argument("--min_allele_frequency", type=float, default=0.1)
	ap.add_argument("--annot_min_freq", type=float, default=0.8, help="Minimum frequency to annotate (>=)")
	ap.add_argument("--annot_min_dp", type=int, default=10, help="Minimum depth to annotate (>)")
	ap.add_argument("--min_quality", type=int, default=12)
	ap.add_argument("--min_mapping", type=int, default=60)
	ap.add_argument("--ploidy", type=int, default=1)
	ap.add_argument("--poll", type=int, default=60, help="Seconds between inspections")
	ap.add_argument("--settle", type=int, default=30, help="Seconds without modification to mark a file as closed")
	ap.add_argument("--idle_timeout", type=int, default=3600, help="Stop after N sec without any new files")
	args = ap.parse_args()
 
	inp, out = os.path.realpath(args.input), os.path.realpath(args.output)
	if not os.path.isdir(inp):
		ap.error("-i doesn't exists: " + args.input)
	if out == inp or out.startswith(inp + os.sep):
		ap.error("-o cannot be inside -i: input folder must remain intact")
	args.input, args.output = inp, out
 
	logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
	os.makedirs(args.output, exist_ok=True)
	samples = read_samples(args.samples)
	clusters = load_info(args.clusters)
	state_file = os.path.join(args.output, "state.json")
	state = json.load(open(state_file)) if os.path.isfile(state_file) else {}
	if state:
		for st in state.values():
			for h in st["history"]:
				h["summary"] = make_summary(h["res"])
		write_reports(args, state)
	last_new = time.time()
 
	while True:
		new_any = False
		for barcode, sample in samples:
			bdir = find_barcode_dir(args.input, barcode)
			if not bdir:
				logger.warning("There is not %s yet (fastq_pass)", barcode)
				continue
			st = state.get(sample, {"history": [], "files": []})
			seen = set(st["files"])
			new = [f for f in ready_files(bdir, args.settle) if os.path.relpath(f, args.input) not in seen]
			if not new:
				continue
			last = st["history"][-1] if st["history"] else {}
			prev = {"n": len(st["files"]), "bam": last.get("bam"), "fq": last.get("fq")}
			try:
				entry = analyse(sample, new, prev, args, clusters)
			except Exception:
				logger.error("Fail in %s:\n%s", sample, traceback.format_exc())
				continue
			st["history"].append(entry)
			st["files"] += [os.path.relpath(f, args.input) for f in new]
			state[sample] = st
			save_state(state_file, state)
			new_any = True
			COLOR = GREEN + BOLD if entry["assigned"] else RED
			logger.info("[%s] " + COLOR + "%s" + END_FORMATTING, sample, "; ".join("Cluster {} ({}/{} markers, {}/{} branch markers)".format(r["cluster"], r["X"], r["Y"], r["x"], r["Z"]) for r in entry["res"] if r["assigned"]) or entry["summary"])
		if new_any:
			write_reports(args, state)
			last_new = time.time()
		if time.time() - last_new > args.idle_timeout:
			logger.info("No new files in %d s: end of pipeline", args.idle_timeout)
			break
		time.sleep(args.poll)
 
 
if __name__ == "__main__":
	main()