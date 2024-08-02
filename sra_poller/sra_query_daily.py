#!/bin/env python3

import os
import xml.parsers.expat
from argparse import ArgumentParser
from datetime import date, timedelta
import pandas as pd
import json


def parse_process_arrays_args(parser: ArgumentParser):
    """Parses the python script arguments from bash and makes sure files/inputs are valid"""
    parser.add_argument(
        "--sra_tracking_dir",
        type=str,
        help="directory to output the results",
        required=True,
    )
    parser.add_argument(
        "--completed_files", type=str, help="File of completed_sra", required=True
    )
    parser.add_argument("--serpent_dir", type=str, help="serpent_dir", required=True)
    parser.add_argument("--config_dir", type=str, help="config_dir", required=True)
    parser.add_argument(
        "--default_config_dir", type=str, help="default_config_dir", required=True
    )


def get_process_arrays_args():
    """Inputs arguments from bash
    Gets the arguments, checks requirements, returns a dictionary of arguments
    Return: args - Arguments as a dictionary
    """
    parser = ArgumentParser()
    parse_process_arrays_args(parser)
    return parser.parse_args()


args = get_process_arrays_args()
# same arguments to a local variable by same name as the argument
sra_tracking_dir = args.sra_tracking_dir
completed_files = args.completed_files
default_config_dir = args.default_config_dir
serpent_dir = args.serpent_dir
config_dir = args.config_dir
downloaded_filepath = os.path.join(sra_tracking_dir, "downloaded_list.txt")
queued_filepath = os.path.join(sra_tracking_dir, "queued_list.txt")
skipped_filepath = os.path.join(sra_tracking_dir, "skipped_list.txt")
todays_date = date.today()
os.makedirs(sra_tracking_dir, exist_ok=True)
year = str(todays_date.year)
month = "{:02d}".format(todays_date.month)
day = "{:02d}".format(todays_date.day)
run_dirname = "SRA_{0}{1}{2}".format(year, month, day)
json_filepath = os.path.join(config_dir, f"{run_dirname}.json")
run_dir = os.path.join(serpent_dir, run_dirname)
print(run_dir)
os.makedirs(run_dir, exist_ok=True)
local_module_in_dir = os.path.join(run_dir, "in_files", "stage_sra")
os.makedirs(local_module_in_dir, exist_ok=True)
print(local_module_in_dir)
# creating the date object of today's date

if True:
    # printing todays date
    print("Current date: ", todays_date)
    # fetching the current year, month and day of today
    start_date = todays_date - timedelta(days=45)
    print("Start year:", start_date.year)
    print("Start month:", start_date.month)
    print("Start day:", start_date.day)

    year_s = str(start_date.year)
    month_s = str("{:02d}".format(start_date.month))
    day_s = str("{:02d}".format(start_date.day))
    search_path = os.path.join(local_module_in_dir, "SRA_Search_Results.html")
    print(
        f"curl -A 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:50.0) Gecko/20100101 Firefox/50.0' 'https://www.ncbi.nlm.nih.gov/sra/?term=sars-cov-2+wastewater' -o {search_path}"
    )
    os.system(
        f"curl -A 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:50.0) Gecko/20100101 Firefox/50.0' 'https://www.ncbi.nlm.nih.gov/sra?term=((%22{year_s}%2F{month_s}%2F{day_s}%22%5BPublication+Date%5D+%3A+%223000%22%5BPublication+Date%5D))+AND+sars-cov-2+wastewater' -o {search_path}"
    )
    MCID = ""
    Key = ""

    with open(search_path, "r") as search_res:
        MCID = search_res.read().split('value="MCID_')[1].split('"')[0]
        search_res.seek(0)
        Key = search_res.read().split("query_key:&quot;")[1].split("&quot")[0]
    xml_out_path = os.path.join(local_module_in_dir, "SRAmetadata.xml")
    print(MCID, Key)
    os.system(
        f"curl 'https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi?rettype=exp&WebEnv=MCID_{MCID}&query_key={Key}' -L -o {xml_out_path}"
    )
    xml_out_path = os.path.join(local_module_in_dir, "SRAmetadata.xml")

    with open(xml_out_path, "r") as xml_fh:
        with open(os.path.join(local_module_in_dir, "SRA_meta.txt"), "w") as txt_fh:
            with open(os.path.join(local_module_in_dir, "SRA_meta.tsv"), "w") as tsv_fh:
                elements = []
                val_dict = {
                    "SRR_acc": "",
                    "BioProject": "",
                    "BioSamp": "",
                    "Submitter": "",
                    "Col_Date": "",
                    "GeoLoc": "",
                    "size": "",
                }
                flag_dict = {
                    "SRR_acc": 0,
                    "BioProject": 0,
                    "BioSamp": 0,
                    "Submitter": 0,
                    "Col_Date": 0,
                    "GeoLoc": 0,
                }
                parse_xml = xml.parsers.expat.ParserCreate()

                def start_element(name, attrs):
                    elements.append(name)
                    if len(elements) > 1:
                        if len(elements) > 2:
                            for ele_ in elements[2:]:
                                txt_fh.write("    ")
                        txt_fh.write(name)
                        txt_fh.write(" :")
                        if attrs:
                            txt_fh.write(" ")
                            txt_fh.write(str(attrs))
                        txt_fh.write("\n")
                        if "BioProject" in str(attrs):
                            flag_dict["BioProject"] = 1
                        if "BioSamp" in str(attrs):
                            flag_dict["BioSamp"] = 1
                        if name == "RUN" and attrs:
                            val_dict["SRR_acc"] = attrs["accession"]
                            if "size" in attrs:
                                val_dict["size"] = attrs["size"]
                            else:
                                val_dict["size"] = "-1"
                        if name == "SUBMISSION" and attrs:
                            try:
                                val_dict["Submitter"] = attrs["center_name"]
                            except:
                                pass
                        if (
                            name == "SUBMITTER_ID"
                            and attrs
                            and not val_dict["Submitter"]
                        ):
                            try:
                                val_dict["Submitter"] = attrs["namespace"]
                            except:
                                pass

                def end_element(name):
                    elements.pop()
                    if elements and len(elements) < 2:
                        tsv_fh.write("\t".join(val_dict.values()))
                        tsv_fh.write("\n")
                        for entry in val_dict:
                            val_dict[entry] = ""
                        txt_fh.write("\n")

                def char_data(data):
                    if data and data.strip():
                        if len(elements) > 1:
                            for ele in elements[1:]:
                                txt_fh.write("    ")
                        txt_fh.write(data)
                        txt_fh.write("\n")
                        if flag_dict["BioProject"] == 1 and data.startswith("PR"):
                            flag_dict["BioProject"] = 0
                            val_dict["BioProject"] = data
                        if flag_dict["BioSamp"] == 1 and data.startswith("SAM"):
                            flag_dict["BioSamp"] = 0
                            val_dict["BioSamp"] = data
                        if flag_dict["Col_Date"] == 1:
                            val_dict["Col_Date"] = data
                            flag_dict["Col_Date"] = 0
                        if data in ("collection_date", "collection date"):
                            flag_dict["Col_Date"] = 1
                        if flag_dict["GeoLoc"] == 1:
                            val_dict["GeoLoc"] += data + ", "
                            flag_dict["GeoLoc"] = 0
                        if (
                            data in ("geo_loc_name", "geo loc name")
                            or "geographic location" in data
                        ):
                            flag_dict["GeoLoc"] = 1

                parse_xml.StartElementHandler = start_element
                parse_xml.EndElementHandler = end_element
                parse_xml.CharacterDataHandler = char_data
                xml_data = xml_fh.read().split(
                    "</EXPERIMENT_PACKAGE_SET>\n<EXPERIMENT_PACKAGE_SET>"
                )
                parse_xml.Parse("\n".join(xml_data))

all_filepath = os.path.join(os.path.join(local_module_in_dir, "SRA_meta.tsv"))
df_completed = pd.read_csv(completed_files)
df_all = pd.read_csv(all_filepath, sep="\t", header=None)
# df_downloaded = pd.read_csv(downloaded_filepath, sep=',', header=None)
df_queued = pd.read_csv(queued_filepath, sep=",", header=None)
df_all = df_all[~df_all[0].isnull()]
df_all = df_all[~df_all[6].isnull()]
df_all[6] = pd.to_numeric(df_all[6])
df_all = df_all[df_all[6] > 1]
all_run_list = set(df_all[0].unique())
queued_list = set(df_queued[0].unique())
all_completed_list = set(df_completed["Run"].unique())
remaining_list = list(all_run_list - all_completed_list - queued_list)
print(len(remaining_list))
if len(remaining_list) < 1:
    print("No new files")

df_all_remain = df_all[df_all[0].isin(remaining_list)]
df_skipped = df_all_remain[df_all_remain[3] == "national institute of biology"]
df_all_remain_filtered = df_all_remain[
    df_all_remain[3] != "national institute of biology"
]
all_filtered_list = set(df_all_remain_filtered[0].unique())
# Output the list to be used in the typical loop
df_all_remain_filtered.to_csv(
    os.path.join(local_module_in_dir, "SRA_meta_run.tsv"),
    sep="\t",
    index=False,
    header=False,
)
# output the list for the skipped with meta data. This may be redundant but that is okay
df_skipped.to_csv(
    os.path.join(local_module_in_dir, "SRA_meta_skipped.tsv"),
    sep="\t",
    index=False,
    header=False,
)
df_skipped_list = set(df_skipped[0].unique())
try:
    df_incoming_skipped = pd.read_csv(skipped_filepath, header=None)
    incoming_skipped_list = set(df_incoming_skipped[0].unique())
except pd.errors.EmptyDataError:
    incoming_skipped_list = set()


df_skipped_added_list = list(df_skipped_list - incoming_skipped_list)
# add to queued list to not re import on the next daily run if it does not finish in time.
if os.path.exists(queued_filepath):
    with open(queued_filepath, "a") as fw:
        for sra_i in all_filtered_list:
            sra_i_s = sra_i.strip()
            fw.write("{0}\n".format(sra_i_s))

# get a skipped list going as a central spot.
if os.path.exists(skipped_filepath):
    with open(skipped_filepath, "a") as fw:
        for sra_i in df_skipped_added_list:
            sra_i_s = sra_i.strip()
            fw.write("{0}\n".format(sra_i_s))
#

json_dict = {
    "all": {"local": {"default_config_dir": default_config_dir}},
    "stage_sra": {
        "local": {
            "mark_as_completed": "False",
            "arguments": {"--node_limit": "30", "--files_per_node": "30"},
        }
    },
}

with open(json_filepath, "w") as j:
    json.dump(json_dict, j)
