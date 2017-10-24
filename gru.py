#!/usr/bin/python

#
# Description: Gru is a small pipeline for analysing (obolete) 2D Nanopore reads.
# it requires bwa, nanopore tools, samtools and bamstats. This pipeline can be extended by
# an assembler such as minasm or spades (together with Illumina reads).
# unfrotunately the 2D reads a not more supported by Nanopore.
# Usage /gru.py job.yml
#

import os, sys, subprocess
import yaml
import tarfile
import shutil
import time
import base64
from lxml import etree
from pprint import pprint
import re
###################################################################################
##                                                                               ##
##                        USAGE: ./gru.py job.yml                                ##
##                                                                               ##
###################################################################################


# debug switch
# from docutils.nodes import paragraph

debug_enabled = False

# default usage mesaaage
msg_usage = "Gru v0.1\n" \
            "Description: pipleine to map nanopore/minion 2D reads to the reference genomes and calulate statisical data\n" \
            "Usage: gru.py <config-file>\n"

# param vars
defaults = {}
params = {}

# folder vars
project_folder = ""
temp_folder = ""
log_folder = ""


##########################################
#    BEGIN SECTION  MESSAGING            #
##########################################

# abort func to close the programm in case of an error
def abort(message, code, usage=False):
    print("ABORT: " + str(message) + " (ERROR " + str(code) + ")")
    if usage:
        print("\n" + msg_usage)
    sys.exit(code)


# warning func to show warnings
def warning(message):
    print("WARNING: " + str(message))


# debug log func to show debugging messages
def debug(message):
    if bool(params["software_settings"]["gru_debug"]) == True:
        print("DEBUG: " + str(message))


##########################################
#    END SECTION  MESSAGING              #
##########################################

##########################################
#    BEGIN SECTION  CONFIGURATION        #
##########################################

# reads configuration file
def read_config(file):
    global params
    try:
        with open(file, 'r') as config_file:
            # load params in YAML format
            cfg = yaml.load(config_file)
        params = cfg  # assign to global params
        debug(params)
    except:
        abort("Could not read the config file", 1, True)


# checks configuration file for correct settings
def check_config():
    global project_folder

    # check project folder existence
    if os.path.exists(params["project_settings"]["project_folder"]):
        if bool(params["project_settings"]["overwrite_folder"]) == False:
            abort("Project output folder exists. But overwrite_folder is not 'True'", 1, True)
    # check nanopore input files tar.gz or as folder
    if params["nanopore_input"].endswith("tar.gz"):
        if not os.path.isfile(params["nanopore_input"]):
            abort("Nanopore input archive does not exist. Wrong path?", 2, True)
    else:
        if not os.path.exists(params["nanopore_input"]):
            abort("Nanopore input folder does not exist. Wrong path?", 2, True)
    # check illumina reads, if given
    if bool(params["illumina_reads"]["enable_illumina"]) == True:
        if not os.path.exists(params["illumina_reads"]["folder"]):
            abort("Illumina reads folder does not exist, but is enabled. Wrong path?", 3, True)
    # check reference files existence
    if bool(params["references"]["enable_references"]) == True:
        if not os.path.exists(params["references"]["folder"]):
            abort("References folder does not exist, but is enabled. Wrong path?", 4, True)
    # check if reference files for given genomes
    for gen in params["file_mapping"]:
        genome = params["file_mapping"][gen]
        if genome.has_key("reference"):
            genome_reference = params["references"]["folder"] + params["file_mapping"][gen]["reference"]
            if not os.path.isfile(genome_reference):
                abort("The reference file for " + gen + " (" + params["file_mapping"][gen][
                    "reference"] + ") does not exist.", 5, True)
        else:
            abort("Reference genome for " + gen + " is not given, but enabled (enable_references: True)", 35)
        # check for illumina files
        if bool(params["illumina_reads"]["enable_illumina"]) == True:
            if genome.has_key("illumina"):
                genome_illumina = params["file_mapping"][gen]["illumina"]
                if isinstance(genome_illumina, str):
                    if not os.path.isfile(params["illumina_reads"]["folder"] + genome_illumina):
                        abort("The Illumina file for " + gen + " (" + params["file_mapping"][gen][
                            "illumina"] + ") does not exist.", 6, True)
                else:
                    for illumina_file in genome_illumina:
                        if not os.path.isfile(params["illumina_reads"]["folder"] + illumina_file):
                            abort("The illumina file for " + gen + " (" + illumina_file + ") does not exist.", 7, True)

    project_folder = params["project_settings"]["project_folder"]


# checks, if all needed software is installed
def check_software():
    # check for bwa
    try:
        p = subprocess.Popen([params["software_general"]["bwa"]], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        # p.stderr
        # p2 = subprocess.Popen(stdin=p.stderr)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("BWA executable could not be started. Wrong path?", 16)
        else:
            abort("Something went wrong by checking for bwa", 17)

    try:
        out = subprocess.check_output([params["software_general"]["samtools"], "--version"])
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("Samtools executable could not be started. Wrong path?", 18)
        else:
            abort("Something went wrong heng li by checking for Samtools", 19)
            raise
    try:
        out = subprocess.check_output([params["software_general"]["r"], "--version"])
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("R executable could not be started. Wrong path?", 20)
        else:
            abort("Something went wrong by checking for R", 21)
            raise

    try:
        p = subprocess.Popen([params["software_general"]["poretools"]], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        # p.stderr
        # p2 = subprocess.Popen(stdin=p.stderr)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("Poretools executable could not be started. Wrong path?", 22)
        else:
            abort("Something went wrong by checking for poretools", 23)


# cleans project folder if overwrite is enabled and folder exists
def clean_project_folder():
    if os.path.exists(params["project_settings"]["project_folder"]):
        if bool(params["project_settings"]["overwrite_folder"]) == True:
            shutil.rmtree(params["project_settings"]["project_folder"])
        debug("Project output folder cleaned")


# creates project folder if not exists
def create_project_folder():
    if not os.path.exists(params["project_settings"]["project_folder"]):
        try:
            os.makedirs(params["project_settings"]["project_folder"])
        except:
            abort("Could not create the project folder!", 24)


# crates temporary folder
def create_temp_folder():
    global temp_folder
    try:
        os.makedirs(project_folder + params["gru_settings"]["temp_foldername"] + "/")
    except:
        abort("Could not create the nanopore output folder folder!", 26)
    temp_folder = project_folder + params["gru_settings"]["temp_foldername"] + "/"


# creates log folder
def create_log_folder():
    global log_folder
    try:
        os.makedirs(project_folder + "log" + "/")
    except:
        abort("Could not create the log folder!", 34)
    log_folder = project_folder + "log" + "/"


# runs all needed configuration funcs
def run_prerequisites():
    check_config()
    check_software()

    clean_project_folder()
    create_project_folder()
    create_temp_folder()
    create_log_folder()


##########################################
#    END SECTION  CONFIGURATION          #
##########################################

##########################################
#    BEGIN SECTION  TOOLS_EXECUTION      #
##########################################

# runs poretools in case to create fasta sequences from the raw files.
def run_poretools():
    # folder names
    nanopore_reads = params["project_settings"]["project_folder"] + params["gru_settings"][
        "nanopore_reads_foldername"] + "/"
    nanopore_fastq = params["project_settings"]["project_folder"] + params["gru_settings"][
        "nanopore_reads_foldername"] + "/"
    # create reads folder
    try:
        os.makedirs(nanopore_reads)
    except:
        abort("Could not create the nanopore output folder!", 25)
    debug("Extracting fast5")
    fastf_folder = params["nanopore_input"]
    poretools_fasta = params["software_general"]["poretools"] + " fasta "
    poretools_fastq = params["software_general"]["poretools"] + " fastq "
    if params["nanopore_input"].endswith("tar.gz"):
        try:
            os.makedirs(temp_folder + "nanopore_fast5" + "/")
        except:
            abort("Could not create temp folder for nanopore extraction!", 27)
        archive = tarfile.open(params["nanopore_input"])
        for member in archive.getmembers():
            if member.isreg():  # skip if the TarInfo is not files
                member.name = os.path.basename(member.name)  # remove the path by reset it
                archive.extract(member, temp_folder + "nanopore_fast5" + "/")  # extract
        fastf_folder = temp_folder + "nanopore_fast5" + "/"
        debug("Extraction of nanopore input done")
    else:
        if not os.path.exists(params["nanopore_input"]):
            abort("Nanopore input folder does not exist", 36)
    poretools_err = open(log_folder + "poretools_error.log", "wb")
    try:
        #print(["export HDF5_DISABLE_VERSION_CHECK=2; find " + fastf_folder + ' -maxdepth 1 -name "*.fast5" -print0 | xargs -0 -I "{}" ' + poretools_fasta + ' "{}" >> ' + nanopore_reads +
        #          params["gru_settings"]["nanopore_reads_filename"] + ".fasta"])

        p = subprocess.Popen(
            ["export HDF5_DISABLE_VERSION_CHECK=2; find " + fastf_folder + ' -maxdepth 1 -name "*.fast5" -print0 | xargs -0 -I "{}" ' + poretools_fasta +  ' "{}" >> ' + nanopore_reads +
             params["gru_settings"]["nanopore_reads_filename"] + ".fasta"], shell=True, stderr=poretools_err)
        # p.stderr
        # p2 = subprocess.Popen(stdin=p.stderr)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("Could not start poretools executable. Wrong path?", 30)
        else:
            abort("Something went wrong by starting poretools", 31)
    debug("Generated fasta reads")
    try:
        p = subprocess.Popen(["export HDF5_DISABLE_VERSION_CHECK=2; find " + fastf_folder + ' -maxdepth 1 -name "*.fast5" -print0 | xargs -0 -I "{}" ' + poretools_fastq + ' "{}" >> ' + nanopore_reads +
                params["gru_settings"]["nanopore_reads_filename"] + ".fastq"], shell=True, stderr=poretools_err)
        # p.stderr
        # p2 = subprocess.Popen(stdin=p.stderr)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("Could not start poretools executable. Wrong path?", 32)
        else:
            abort("Something went wrong by starting poretools", 33)
    debug("Generated fastq reads")


# builds index files using bwa_index
def run_bwa_index():
    if not bool(params["references"]["enable_references"]) == True:
        warning("Skipped mapping to the reference genome because of disabled reference mapping")
        return
    try:
        os.makedirs(project_folder + params["gru_settings"]["bwa_reference_index_foldername"] + "/")
    except:
        abort("Could not create bwa index output folder!", 37)
    bwa_output = project_folder + params["gru_settings"]["bwa_reference_index_foldername"] + "/"
    allcontigs = open(bwa_output + "_contigs.fasta", "wb")
    # load for each genome reference file and read it. Enumerate and add prefix
    for gen in params["file_mapping"]:
        contig_count = 0
        prefix = params["file_mapping"][gen]["prefix"]
        with open(params["references"]["folder"] + params["file_mapping"][gen]["reference"], "r") as reference_file:
            for line in reference_file:
                if line.startswith('>'):
                    allcontigs.write('>' + prefix + "_contig_" + str(contig_count) + "\n")
                    contig_count = contig_count + 1
                else:
                    allcontigs.write(line)

    allcontigs.close()
    contigfile = bwa_output + "_contigs.fasta"
    bwa_index_err = open(log_folder + "bwa_index_error.log", "wb")
    # print "Run bwa here"

    try:
        p = subprocess.Popen(
            [params["software_general"]["bwa"], "index", contigfile], stderr=bwa_index_err)
        # p.stderr
        # p2 = subprocess.Popen(stdin=p.stderr)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("Could not start bwa indexing.. Wrong path?", 32)
        else:
            abort("Something went wrong by starting bwa indexing", 33)


# mapps reads using bwa_mapper
def run_bwa_mapping():
    # create mapping folder
    try:
        os.makedirs(project_folder + params["gru_settings"]["mapping_foldername"] + "/")
    except:
        abort("Could not create bwa index output folder!", 38)

    index_contigs = project_folder + params["gru_settings"]["bwa_reference_index_foldername"] + "/" + "_contigs.fasta"
    reads_fastq = project_folder + params["gru_settings"]["nanopore_reads_foldername"] + "/" + params["gru_settings"][
        "nanopore_reads_filename"] + ".fastq"
    bwa_command = params["software_general"]["bwa"] + " mem -x ont2d -t " + str(params["software_settings"][
                                                                                    "mapping_threads"]) + " " + index_contigs + " " + reads_fastq + " > " + project_folder + \
                  params["gru_settings"]["mapping_foldername"] + "/" + "mapped.sam"
    debug("running " + bwa_command);

    try:
        p = subprocess.Popen(bwa_command, shell=True)
        # p.stderr
        # p2 = subprocess.Popen(stdin=p.stderr)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("Could not start bwa mapping.. Wrong command?", 32)
        else:
            abort("Something went wrong by starting bwa mapping", 33)


# splits mapping files into single files, takes care about unique identifiers.
def process_mappingfile():
    mapping_folder = project_folder + params["gru_settings"]["mapping_foldername"] + "/"
    splitted_mapping_folder = project_folder + params["gru_settings"]["mapping_foldername"] + "/splitted/"
    mapping_file = project_folder + params["gru_settings"]["mapping_foldername"] + "/" + "mapped.sam"
    stats_folder = project_folder + "stats/"

    try:
        os.makedirs(splitted_mapping_folder)
    except:
        abort("Could not create bwa index output folder!", 38)

    prefix_set = ""

    for gen in params["file_mapping"]:
        prefix_set = prefix_set + params["file_mapping"][gen]["prefix"] + " "

    debug("Prefix set: " + prefix_set)

    split_command = "for i in " + prefix_set + "; do grep $i " + mapping_file + " > " + splitted_mapping_folder + "${i}.sam; done"

    try:
        p = subprocess.Popen(split_command, shell=True)
        # p.stderr
        # p2 = subprocess.Popen(stdin=p.stderr)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("Could not split mappings into single files.. Wrong command?", 32)
        else:
            abort("Something went wrong by splitting the mappings into single files", 33)

    sort_command = "for i in `ls " + splitted_mapping_folder + "*sam | sed 's/.sam//'`; do echo $i; samtools view -buh $i.sam | samtools sort -@ " + str(
        params["software_settings"]["sorting_threads"]) + " -o $i.bam -; done"

    debug("Sorting cmd: " + sort_command)
    try:
        p = subprocess.Popen(sort_command, shell=True)
        # p.stderr
        # p2 = subprocess.Popen(stdin=p.stderr)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("Could not split mappings into single files.. Wrong command?", 32)
        else:
            abort("Something went wrong by splitting the mappings into single files", 33)

    index_command = "for i in " + splitted_mapping_folder + "*bam; do samtools index $i; done"

    try:
        p = subprocess.Popen(index_command, shell=True)
        # p.stderr
        # p2 = subprocess.Popen(stdin=p.stderr)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("Could not split mappings into single files.. Wrong command?", 32)
        else:
            abort("Something went wrong by splitting the mappings into single files", 33)

    stats_command = "for i in " + splitted_mapping_folder + "*bam; do samtools stats $i > $i.stats & done"

    try:
        p = subprocess.Popen(stats_command, shell=True)
        # p.stderr
        # p2 = subprocess.Popen(stdin=p.stderr)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("Could not split mappings into single files.. Wrong command?", 32)
        else:
            abort("Something went wrong by splitting the mappings into single files", 33)

    flagstats_command = "for i in " + splitted_mapping_folder + "*bam; do samtools flagstat $i > $i.flagstat & done"

    try:
        p = subprocess.Popen(flagstats_command, shell=True)
        # p.stderr
        # p2 = subprocess.Popen(stdin=p.stderr)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("Could not split mappings into single files.. Wrong command?", 32)
        else:
            abort("Something went wrong by splitting the mappings into single files", 33)

    try:
        os.makedirs(stats_folder)
    except:
        abort("Could not create bwa index output folder!", 38)

        # try:
        #     p = subprocess.Popen("cp " + splitted_mapping_folder + "*.stats" + stats_folder, shell=True)
        #     # p.stderr
        #     # p2 = subprocess.Popen(stdin=p.stderr)
        #     out, err = p.communicate()
        # except OSError as e:
        #     if e.errno == os.errno.ENOENT:
        #         abort("Could not split mappings into single files.. Wrong command?", 32)
        #     else:
        #         abort("Something went wrong by splitting the mappings into single files", 33)
        #
        # try:
        #     p = subprocess.Popen("cp " + splitted_mapping_folder + "*.flagstat " + stats_folder, shell=True)
        #     # p.stderr
        #     # p2 = subprocess.Popen(stdin=p.stderr)
        #     out, err = p.communicate()
        # except OSError as e:
        #     if e.errno == os.errno.ENOENT:
        #         abort("Could not split mappings into single files.. Wrong command?", 32)
        #     else:
        #         abort("Something went wrong by splitting the mappings into single files", 33)


# run assembler
# TODO
##########################################
#    END SECTION  TOOLS_EXECUTION        #
##########################################

##########################################
#    BEGIN SECTION  STATISTICS           #
##########################################

# creates statistics from bam files using bamstats
def plot_bamstats():
    time.sleep(10)  # take 10 secons to snyc files

    stats_folder = project_folder + "stats/"
    splitted_mapping_folder = project_folder + params["gru_settings"]["mapping_foldername"] + "/splitted/"
    # plot_stats_command = params["software_general"]["plot_bamstats"] + " -p " + stats_folder + " " + stats_folder + "*.stats"

    # plot_stats_command = "for i in " + stats_folder + "*.stats; do  plot-bamstats -p stats_$i/ $i; done"
    plot_stats_command = "for i in `ls " + splitted_mapping_folder + "*.stats`; do base=$(basename $i | sed 's/\///' | sed 's/.bam.stats//'); " + \
                         params["software_general"]["plot_bamstats"] + " -p " + stats_folder + "$base/ $i; done"

    print(plot_stats_command)
    try:
        p = subprocess.Popen(plot_stats_command, shell=True)
        # p.stderr
        # p2 = subprocess.Popen(stdin=p.stderr)
        out, err = p.communicate()
        time.sleep(10)
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            abort("Could not create statisic charts.. Wrong command?", 32)
        else:
            abort("Something went wrong by creating statistic charts", 33)


##########################################
#    BEGIN SECTION  RENDERING            #
##########################################
def render_output():
    output = ""
    output_file = open(params["project_settings"]["project_folder"] + 'gru-output.html', 'w')
    with open('template.html', 'r') as template:
        output = template.read()

    ##produce captions + menu entreis
    gru_menu_general = '<h3 class="gru_menu_general">General</h3>'
    gru_menu_general += '<ul class="nav nav-sidebar"> <li class="active"><a href="#" data-paneclass="gru-general-overview">Input overview</a></li>' \
                        '<li><a href="#" data-paneclass="gru-general-organisms">Organisms</a></li>' \
                        '</ul>'
    output = output.replace("<gru-menu-general/>", gru_menu_general)

    gru_menu_organism_staticstics = '<h3>Organism statistics</h3><ul class="nav nav-sidebar">'
    for gen in params["file_mapping"]:
        prefix = params["file_mapping"][gen]["prefix"] + " "
        gru_menu_organism_staticstics += '<li><a href="#" data-paneclass="gru-stats-' + prefix + '">' + prefix + '</a></li>'  # e.g. gru-stats-nc201

    gru_menu_organism_staticstics += '</ul>'
    output = output.replace("<gru-menu-organism-statistics/>", gru_menu_organism_staticstics)

    gru_menu_assemblies = '<h3>Assemblies</h3>' \
                          '<ul class="nav nav-sidebar">' \
                          '<li><a href="#" data-paneclass="gru-assembler-overview">Assembler overview</a></li> <li><a href="#" data-paneclass="gru-assembler-quast">Quast report</a></li>' \
                          '</ul>'
    output = output.replace("<gru-menu-assemblies/>", gru_menu_assemblies)

    gru_menu_logs = '<h3>Logs</h3>' \
                    '<ul class="nav nav-sidebar">' \
                    '<li><a href="#" data-paneclass="">gru logs</a></li>' \
                    '<li><a href="#" data-paneclass="">3rd party software logs</a></li></ul>'

    output = output.replace("<gru-menu-logs/>", gru_menu_logs)

    # content - general
    ## Input overview
    gru_content_settings_general = "<tr><td>Project folder</td><td>" + params["project_settings"][
        "project_folder"] + "</td></tr>"
    gru_content_settings_general += "<tr><td>Enable overwrite</td><td>" + "Yes" if bool(
        params["project_settings"]["overwrite_folder"]) == True else "No" + "</td></tr>"
    gru_content_settings_general += "<tr><td>Nanopore input</td><td>" + params["nanopore_input"] + "</td></tr>"
    gru_content_settings_general += "<tr><td>Enable references</td><td>" + "Yes" if bool(
        params["references"]["enable_references"]) == True else "No" + "</td></tr>"
    gru_content_settings_general += "<tr><td>References folder</td><td>" + params["references"]["folder"] + "</td></tr>"

    output = output.replace("<gru-content-settings-general/>", gru_content_settings_general)

    ## GRU content file mapping

    gru_content_file_mapping = ""
    for gen in params["file_mapping"]:
        genome = params["file_mapping"][gen]
        prefix = params["file_mapping"][gen]["prefix"]
        reference = params["file_mapping"][gen]["reference"]
        illumina_reads = ""
        if genome.has_key("illumina"):
            genome_illumina = params["file_mapping"][gen]["illumina"]
            if isinstance(genome_illumina, str):
                illumina_reads = genome_illumina
            else:
                for illumina_read in genome_illumina:
                    illumina_reads += illumina_read + "</br>"
        else:
            illumina_reads = "-"

        gru_content_file_mapping += "<tr><td>" + prefix + "</td><td>" + prefix + "</td><td>" + reference + "</td><td>" + illumina_reads + "</td></tr>"
    output = output.replace("<gru-content-file-mapping/>", gru_content_file_mapping)

    ## Software

    gru_content_software_used = ""
    for programm in params["software_general"]:
        gru_content_file_mapping += "<tr><td>" + programm + "</td><td>" + params["software_general"][
            programm] + "</td></tr>"
    output = output.replace("<gru-content-file-mapping/>", gru_content_software_used)

    ##Used assembler
    gru_content_assembler_used = ""
    # TODO


    # Organism details
    # TODO

    # Organism statistics
    gru_content_organisms_stats = ""
    stats_folder = project_folder + "stats/"
    for organism in params["file_mapping"]:
        prefix = params["file_mapping"][organism]["prefix"]
        organism_stat = '<div class="col-sm-9 col-sm-offset-3 col-md-10 col-md-offset-2 main gru-output gru-stats-' + prefix + ' hidden">' \
                                                                                                                               '<h1 class="page-header">Statistics ' + prefix + '</h1>' \
                                                                                                                                                                                '<h2>Overview</h2>'

        # read stats from bamstats html and strip it
        bamstats_file = open(stats_folder + prefix + "/index.html")
        bamstats_s = bamstats_file.read()
        bamstats_html = etree.HTML(bamstats_s)
        reads_stats_a = [stat.strip(' ') for stat in
                         bamstats_html.xpath('//table[@class="nums"]/tr[2]/td/table/tr/td[2]//text()')]
        reads_stats_b = [stat.strip(' ') for stat in
                         bamstats_html.xpath('//table[@class="nums"]/tr[2]/td/table/tr/td[3]//text()')]
        bases_stats_a = [stat.strip(' ') for stat in
                         bamstats_html.xpath('//table[@class="nums"]/tr[4]/td/table/tr/td[2]//text()')]
        bases_stats_b = [stat.strip(' ') for stat in
                         bamstats_html.xpath('//table[@class="nums"]/tr[4]/td/table/tr/td[3]//text()')]

        reads_mapped = reads_stats_b[3].translate(None, '()%').replace(',', '.')
        organism_stat += '<div class="row">' \
                        '<div class="col-sm-8">' \
                        '<div class="progress">' \
                        '<div class="progress-bar progress-bar-success progress-bar-striped" role="progressbar" aria-valuenow="' + reads_mapped +'" aria-valuemin="0" aria-valuemax="100" style="width: ' + reads_mapped + '%" >' \
                        '<span style="color: black;">' + reads_mapped + '% reads mapped</span>' \
                        '</div>' \
                        '</div>' \
                        '</div>' \
                        '</div>'

        organism_stat += '<div class="row">' \
                         '<div class="col-sm-4">' \
                         '<table class="nums">' \
                         '<tbody><tr><th>Reads</th></tr>' \
                         '<tr>' \
                         '<td class="pad"><table>' \
                         '<tbody><tr><td>total: </td><td class="right">' + reads_stats_a[0] + '</td><td class="right"></td></tr>' \
                         '<tr><td>filtered: </td><td class="right">' + reads_stats_a[1] + '</td><td class="right">' + reads_stats_b[0] + '</td></tr>' \
                         '<tr><td>non-primary: </td><td class="right">' + reads_stats_a[2] + '</td><td class="right"> </td></tr>' \
                         '<tr><td>duplicated: </td><td class="right">' + reads_stats_a[3] + '</td><td class="right">' + reads_stats_b[2] + '</td></tr>' \
                         '<tr><td>mapped: </td><td class="right">' + reads_stats_a[4] + '</td><td class="right">' + reads_stats_b[3] + '</td></tr>' \
                         '<tr><td>zero MQ: </td><td class="right">' + reads_stats_a[5] + '</td><td class="right">' + reads_stats_b[4] + '</td></tr>' \
                         '<tr><td>avg read length: </td><td class="right">' + reads_stats_a[6] + '</td><td class="right"></td></tr>' \
                         '</tbody></table></td></tr>' \
                         '<tr><th>Bases</th></tr>' \
                         '<tr>' \
                         '<td class="pad"><table>' \
                         '<tbody><tr><td>total: </td><td class="right">' + bases_stats_a[0] + '</td><td class="right">' + bases_stats_b[0] + '</td></tr>' \
                         '<tr><td>mapped: </td><td class="right">' + bases_stats_a[1] + '</td><td class="right"></td></tr>' \
                         '<tr><td>error rate: </td><td class="right">' + bases_stats_a[2] + '</td><td class="right"></td></tr>' \
                         '</tbody></table></td></tr>' \
                         '</tbody></table>' \
                         '</div>' \
                         '<div class="col-sm-4">' \
                         '<table class="nums">' \
                         '<tbody><tr><th>Futher stats</th></tr>' \
                         '<tr>' \
                         '<td class="pad"><table>' \
                         '<tbody><tr><td>Key</td><td class="right"> value</td><td class="right"></td></tr>' \
                         '</tbody></table></td></tr>' \
                         '</tbody></table>' \
                         '</div>' \
                         '</div>'

        organism_stat += '<div class="row gru-bamstats-charts">'\
                '<h2>Bamstats charts'\
                '<a tabindex="0" class="btn btn-default gru-popover" role="button" data-toggle="popover" data-trigger="focus" title="Information" data-content="These charts were generated by samtools\' plot-bamstats. For further description and meaning of the charts please consider the samtools manual."><span class="glyphicon glyphicon-info-sign" aria-hidden="true"></span></a>'\
                '</h2>'
        stat_graphs = [("gc-content", "GC Content"),
                       ("coverage", "Coverage"),
                       ("quals", "Quality per cycle"),
                       ("quals2", "Quality per cycle"),
                       ("quals3", "Quality per cycle"),
                       ("quals-hm", "Quality per cycle"),
                       ("acgt-cycles", "Per-base sequence content"),
                       ("gc-depth", "Mapped depth vs GC"),
                       ("indel-cycles", "InDels per cycle"),
                       ("indel-dist", "InDel length")]

        for graph, description in stat_graphs:
            graph_file = open(stats_folder + prefix + "/" + graph + ".png")
            print (graph_file.name)
            organism_stat += '<div class="col-sm-6 col-md-6">' \
                                 '<div class="thumbnail">' \
                             '<img src="data:image/png;base64,' + base64.b64encode(graph_file.read()) + '" />' \
                                    '<div class="caption"><h3>' + description + '</h3></div>' \
                                 '</div>' \
                             '</div>'
            graph_file.close()
        organism_stat += '</div>'

        organism_stat += '</div>'

        gru_content_organisms_stats += organism_stat

    output = output.replace("<gru-content-organisms-stats/>", gru_content_organisms_stats)

    output_file.write(output)
    output_file.close()


##########################################
#    END SECTION  RENDERING              #
##########################################


##########################################
#    BEGIN SECTION  MAIN_ROUTINE         #
##########################################
def main(argv):
    read_config(argv[0])
    run_prerequisites()
    run_poretools()
    run_bwa_index()
    run_bwa_mapping()
    process_mappingfile()
    plot_bamstats()
    # TODO more statistics. Compare assemblies!!!

    render_output()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        abort("No argumens given.", 1, True)
    if len(sys.argv) > 2:
        abort("Too much argumens given.", 1, True)
    main(sys.argv[1:])  # start main with arguments

    ##########################################
    #    END SECTION  MAIN_ROUTINE           #
    ##########################################
