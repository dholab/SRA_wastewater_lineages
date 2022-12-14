#!/bin/env python3

import os
import xml.parsers.expat
from argparse import ArgumentParser


def parse_process_arrays_args(parser: ArgumentParser):
    """Parses the python script arguments from bash and makes sure files/inputs are valid"""
    parser.add_argument('--out_dir',
                        type=str,
                        help='directory to output the results',
                        required=True)


def get_process_arrays_args():
    """	Inputs arguments from bash
    Gets the arguments, checks requirements, returns a dictionary of arguments
    Return: args - Arguments as a dictionary
    """
    parser = ArgumentParser()
    parse_process_arrays_args(parser)
    return parser.parse_args()


args = get_process_arrays_args()
# same arguments to a local variable by same name as the argument
out_dir = args.out_dir
search_path = os.path.join(out_dir, "SRA_Search_Results.html")
os.system(
    f"curl -A 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:50.0) Gecko/20100101 Firefox/50.0' 'https://www.ncbi.nlm.nih.gov/sra/?term=sars-cov-2+wastewater' -o {search_path}")
MCID = ""
Key = ""

with open(search_path, "r") as search_res:
    MCID = search_res.read().split('value="MCID_')[1].split('"')[0]
    search_res.seek(0)
    Key = search_res.read().split("query_key:&quot;")[1].split("&quot")[0]
xml_out_path = os.path.join(out_dir, "SRAmetadata.xml")
os.system(
    f"curl 'https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi?rettype=exp&WebEnv=MCID_{MCID}&query_key={Key}' -L -o {xml_out_path}")
xml_out_path = os.path.join(out_dir, "SRAmetadata.xml")
with open(xml_out_path, "r") as xml_fh:
    with open(os.path.join(out_dir, "SRA_meta.txt"), "w") as txt_fh:
        with open(os.path.join(out_dir, "SRA_meta.tsv"), "w") as tsv_fh:

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
                        if 'size' in attrs:
                            val_dict["size"] = attrs["size"]
                        else:
                            val_dict["size"] = '-1'
                    if name == "SUBMISSION" and attrs:
                        try:
                            val_dict["Submitter"] = attrs["center_name"]
                        except:
                            pass
                    if name == "SUBMITTER_ID" and attrs and not val_dict["Submitter"]:
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
                    if data in ("geo_loc_name", "geo loc name") or "geographic location" in data:
                        flag_dict["GeoLoc"] = 1


            parse_xml.StartElementHandler = start_element
            parse_xml.EndElementHandler = end_element
            parse_xml.CharacterDataHandler = char_data
            parse_xml.Parse(xml_fh.read())
