#!/bin/env python3
from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Self

from jsonschema import validate
from jsonschema.exceptions import ValidationError
from loguru import logger


@dataclass
class VariantConfig:
    """
    TODO
    """

    nt_positions: list[str]
    nt_call_variants: list[str]
    nt_multi_variants: list[str]
    aa_variants: list[str]
    multi_var: list[str]
    omis: dict[str, dict[str, list[str]]]
    omi_matches: dict[str, Any]
    nt_outfiles: dict[str, Any]
    nt_mult_outfiles: dict[str, Any]
    aa_outfiles: dict[str, Any]
    multivar_fhs: dict[str, Any]
    nt_call_dict: dict[str, Any]
    nt_SRAs: dict[str, str]
    nt_sample_dict: dict[str, Any]
    nt_multi_sample_dict: dict[str, str]
    ref_sample_dict: dict[str, str]
    pos_sample_dict: dict[str, Any]
    ref_nt_dict: dict[str, str]
    search_omis: int = 0
    search_ref: int = 0

    # TODO
    # this section handles initializing writers for many of the input files that can be accessed
    # for each potential list of input files. We'll need to do the following:
    # 1. move all into a function that handles initializing writers.
    # 2. think about alternative writer container data structures or switching to just hold the path
    # 3. at minimum, rename the lists of files (or whatever) in the JSON-based dataclass to be more intuitive
    # 4. replace + overload string concatenation with f-strings
    # Modifiable source code block for handling NT vall vars
    def init_writers(self) -> Self:
        if variant_config.nt_call_variants:
            i = 2  # why does this start at 2? Is nt_call_variants only initialized when there are at least two? How to validate that if so?
            for variant in variant_config.nt_call_variants:
                out_name = Path(
                    f"{samps}_NT_{i}.tsv",
                )  # why not just store the filename for automatic and more timebound opening of files with context managers?
                variant_config.nt_outfiles[variant] = out_name.open(
                    "w",
                )  # why load a writer into a dictionary for each file? Is there a more appropriate data structure here?
                i += 1

        if (
            variant_config.nt_multi_variants
        ):  # how is this different from nt_call_variants?
            for variant in variant_config.nt_multi_variants:
                out_name = Path(f"{samps}_NTmult_{variant}.tsv")
                variant_config.nt_mult_outfiles[variant] = out_name.open("w")

        if variant_config.aa_variants:
            for variant in variant_config.aa_variants:
                out_name = Path(f"{samps}_AA_{variant}.tsv")
                variant_config.aa_outfiles[variant] = out_name.open(
                    "w",
                )

        if variant_config.multi_var:
            i = 1
            for variant in variant_config.multi_var:
                out_name = Path(f"{samps}_multi{i}.tsv")
                variant_config.multivar_fhs[variant] = out_name.open(
                    "w",
                )
                variant_config.multivar_fhs[variant].write(variant)
                variant_config.multivar_fhs[variant].write("\n")
                i += 1

        return self


# Open both the json and schema which is just the json name with a .schema extension
def validate_json(json_file):
    try:
        with Path(json_file).with_suffix(".schema").open("r") as schema_file:
            schema = json.load(schema_file)
        validate(instance=json.load(json_file), schema=schema)
        logger.info("Validation passed successfully")
    except ValidationError:
        logger.warning("Validation Error with json validation")
        exit(1)
    except FileNotFoundError as e:
        logger.warning(f"File not found: {e.filename}")


def populate_config(json_file: str) -> VariantConfig:
    """
    TODO
    """

    assert (
        Path(json_file).is_file()
    ), f"Please double check that the provided config JSON, {json_file} exists."
    assert validate_json(json_file)

    with Path(json_file).open("r") as variant_extractor_json:
        all_vars_json = json.load(variant_extractor_json)
        return VariantConfig(
            nt_positions=all_vars_json["variant_config.nt_positions"],
            nt_call_variants=all_vars_json["variant_config.nt_call_variants"],
            nt_multi_variants=all_vars_json["variant_config.nt_multi_variants"],
            aa_variants=all_vars_json["variant_config.aa_variants"],
            multi_var=all_vars_json["multi_var"],
            omis=all_vars_json["variant_config.omis"],
            omi_matches=all_vars_json["omi_matches"],
            nt_outfiles=all_vars_json["nt_outfiles"],
            nt_mult_outfiles=all_vars_json["nt_mult_outfiles"],
            aa_outfiles=all_vars_json["aa_outfiles"],
            multivar_fhs=all_vars_json["multivar_fhs"],
            nt_call_dict=all_vars_json["nt_call_dict"],
            nt_SRAs=all_vars_json["nt_SRAs"],
            nt_sample_dict=all_vars_json["nt_sample_dict"],
            nt_multi_sample_dict=all_vars_json["nt_multi_sample_dict"],
            ref_sample_dict=all_vars_json["ref_sample_dict"],
            pos_sample_dict=all_vars_json["pos_sample_dict"],
            ref_nt_dict=all_vars_json["ref_nt_dict"],
            search_omis=all_vars_json["search_variant_config.omis"],
            search_ref=all_vars_json["search_ref"],
        )


# UGGGLY
samps = sys.argv[1]

CONFIG_PATH = "../config/variant_extractor.json"
variant_config = populate_config(CONFIG_PATH)
variant_config.init_writers()

# TODO
files_read = []

testcounter = 0
for subdir, dirs, files in os.walk(os.getcwd()):
    for file in files:
        if file not in files_read and "Assemblies" not in subdir:
            files_read.append(file)
            if (
                (file.endswith("_unique_seqs.tsv"))
                and (
                    variant_config.aa_variants
                    or variant_config.multi_var
                    or variant_config.nt_multi_variants
                )
                and "_AA_" not in file
            ):  # and not 'wgs' in file:  or file.endswith('_reads.tsv') or file.endswith('_covars.tsv')
                in_file = open(os.path.join(subdir, file))
                nt_multi_matches = {}
                aa_matches = {}
                mv_matches = {}
                counts = 0
                ft = 0
                if file.endswith("_reads.tsv"):
                    ft = 1

                for line in in_file:
                    splitline = line.split("\t")
                    try:
                        splitline[1]
                    except:
                        pass
                    else:
                        if splitline[1] != "Count":
                            if variant_config.nt_multi_variants:
                                for variant in variant_config.nt_multi_variants:
                                    mismatch = 0
                                    first = variant.split("and")[0][1:-1]
                                    try:
                                        if int(splitline[ft].split(" ")[0]) > int(
                                            first,
                                        ):
                                            continue
                                    except:
                                        pass

                                    for PM in variant.split("and"):
                                        if "not" in PM:
                                            curPM = PM.strip("not")
                                            if curPM[0] == curPM[-1]:
                                                if PM[:-1] in splitline[ft]:
                                                    # for nt in ['A', 'T', 'C', 'G', 'N', '-']:
                                                    # if PM[:-1]+nt in splitline[0]:
                                                    mismatch += 1
                                            elif curPM in splitline[ft]:
                                                mismatch += 1
                                        elif PM[0] == PM[-1]:
                                            if PM[:-1] in splitline[ft]:
                                                # for nt in ['A', 'T', 'C', 'G', 'N', '-']:
                                                # if PM[:-1]+nt in splitline[0]:
                                                mismatch += 1
                                        elif PM not in splitline[ft]:
                                            mismatch += 1
                                    if mismatch == 0:
                                        try:
                                            nt_multi_matches[variant].append(
                                                line,
                                            )  # str(mismatch) + "\t" + line)
                                            # counts += int(splitline[1])
                                        except:
                                            nt_multi_matches[variant] = [line]
                            if ft == 1 or (int(splitline[1]) >= 4):
                                if variant_config.multi_var:
                                    for variant in variant_config.multi_var:
                                        mismatch = 0
                                        match = 0
                                        for PM in variant.split("and"):
                                            if "not" in PM:
                                                if PM.strip("not") in splitline[ft]:
                                                    mismatch += 1
                                            elif PM in splitline[ft]:
                                                match += 1

                                        if (
                                            ((match - mismatch) > 1)
                                            and len(splitline[ft].split(" ")) < 25
                                            and splitline[ft].count("insert") < 3
                                            and splitline[ft].count("Del") < 4
                                        ):
                                            try:
                                                mv_matches[variant].append(line)
                                            except:
                                                mv_matches[variant] = [line]
                            if variant_config.aa_variants:
                                if int(splitline[1]) > 0 or ft == 1:
                                    for variant in variant_config.aa_variants:
                                        if variant in splitline[ft]:
                                            try:
                                                aa_matches[variant].append(line)
                                            except:
                                                aa_matches[variant] = [line]

                if nt_multi_matches:
                    for variant in nt_multi_matches:
                        variant_config.nt_mult_outfiles[variant].write(file[:-4] + "\t")
                        variant_config.nt_mult_outfiles[variant].write("\n")
                        for line in nt_multi_matches[variant]:
                            variant_config.nt_mult_outfiles[variant].write(line)
                if mv_matches:
                    for variant in mv_matches:
                        variant_config.multivar_fhs[variant].write(
                            subdir + "/" + file[:-4] + "\t",
                        )
                        variant_config.multivar_fhs[variant].write("\n")
                        for line in mv_matches[variant]:
                            variant_config.multivar_fhs[variant].write(line)
                if aa_matches:
                    for variant in aa_matches:
                        if aa_matches[variant]:
                            variant_config.aa_outfiles[variant].write(
                                subdir + "/" + file + "\n",
                            )
                            for line in aa_matches[variant]:
                                variant_config.aa_outfiles[variant].write(line)

                in_file.close()

            if file.endswith("_nt_calls.tsv") and (
                variant_config.nt_call_variants or variant_config.nt_positions
            ):
                in_file = open(os.path.join(subdir, file))
                match_lines = {}
                match_dict = {}
                ref_dict = {}
                for line in in_file:
                    splitline = line.split("\t")
                    try:
                        splitline[1]
                    except:
                        pass
                    else:
                        if splitline[0] == "Position":
                            PosLine = line
                            for i in range(len(splitline)):
                                variant_config.nt_call_dict[splitline[i]] = i
                        if variant_config.nt_positions:
                            if splitline[0] in variant_config.nt_positions:
                                try:
                                    variant_config.pos_sample_dict[
                                        file.split("_nt_calls")[0]
                                    ][splitline[0]] = int(
                                        splitline[
                                            int(
                                                variant_config.nt_call_dict[
                                                    splitline[1]
                                                ],
                                            )
                                        ],
                                    ) / int(
                                        splitline[
                                            int(variant_config.nt_call_dict["Total"])
                                        ],
                                    )
                                except:
                                    variant_config.pos_sample_dict[
                                        file.split("_nt_calls")[0]
                                    ] = {
                                        splitline[0]: (
                                            int(
                                                splitline[
                                                    int(
                                                        variant_config.nt_call_dict[
                                                            splitline[1]
                                                        ],
                                                    )
                                                ],
                                            )
                                            / int(
                                                splitline[
                                                    int(
                                                        variant_config.nt_call_dict[
                                                            "Total"
                                                        ],
                                                    )
                                                ],
                                            )
                                        ),
                                    }
                                try:
                                    variant_config.ref_nt_dict[splitline[0]]
                                except:
                                    variant_config.ref_nt_dict[
                                        splitline[0]
                                    ] = splitline[1]

                        if variant_config.nt_call_variants:
                            for var in variant_config.nt_call_variants:
                                try:
                                    match_lines[var]
                                except:
                                    match_lines[var] = []
                                # if 'and' in var:
                                varPMs = var.split("and")
                                for subvar in varPMs:
                                    if subvar[1:-1] == splitline[0]:
                                        if int(
                                            splitline[
                                                int(
                                                    variant_config.nt_call_dict[
                                                        "Total"
                                                    ],
                                                )
                                            ],
                                        ) > 50 and (
                                            int(
                                                splitline[
                                                    int(
                                                        variant_config.nt_call_dict[
                                                            subvar[-1]
                                                        ],
                                                    )
                                                ],
                                            )
                                            > (
                                                0.2
                                                * int(
                                                    splitline[
                                                        int(
                                                            variant_config.nt_call_dict[
                                                                "Total"
                                                            ],
                                                        )
                                                    ],
                                                )
                                            )
                                        ):
                                            match_lines[var].append(line)
                                            try:
                                                match_dict[var][subvar] = int(
                                                    splitline[
                                                        int(
                                                            variant_config.nt_call_dict[
                                                                subvar[-1]
                                                            ],
                                                        )
                                                    ],
                                                ) / int(
                                                    splitline[
                                                        int(
                                                            variant_config.nt_call_dict[
                                                                "Total"
                                                            ],
                                                        )
                                                    ],
                                                )
                                            except:
                                                match_dict[var] = {
                                                    subvar: (
                                                        int(
                                                            splitline[
                                                                int(
                                                                    variant_config.nt_call_dict[
                                                                        subvar[-1]
                                                                    ],
                                                                )
                                                            ],
                                                        )
                                                        / int(
                                                            splitline[
                                                                int(
                                                                    variant_config.nt_call_dict[
                                                                        "Total"
                                                                    ],
                                                                )
                                                            ],
                                                        )
                                                    ),
                                                }
                                            try:
                                                ref_dict[var][subvar] = int(
                                                    splitline[
                                                        int(
                                                            variant_config.nt_call_dict[
                                                                subvar[0]
                                                            ],
                                                        )
                                                    ],
                                                ) / int(
                                                    splitline[
                                                        int(
                                                            variant_config.nt_call_dict[
                                                                "Total"
                                                            ],
                                                        )
                                                    ],
                                                )
                                            except:
                                                ref_dict[var] = {
                                                    subvar: (
                                                        int(
                                                            splitline[
                                                                int(
                                                                    variant_config.nt_call_dict[
                                                                        subvar[0]
                                                                    ],
                                                                )
                                                            ],
                                                        )
                                                        / int(
                                                            splitline[
                                                                int(
                                                                    variant_config.nt_call_dict[
                                                                        "Total"
                                                                    ],
                                                                )
                                                            ],
                                                        )
                                                    ),
                                                }
                            # elif var[1:-1] == splitline[0]:
                            # if int(splitline[int(nt_call_dict['Total'])]): # > 50 and (int(splitline[int(nt_call_dict[subvar[-1]])]) > (.01 * int(splitline[int(nt_call_dict['Total'])]))):
                            # match_lines[var].append(line)

                if match_lines:
                    for variant in match_lines:
                        if len(match_lines[variant]) == len(variant.split("and")):
                            try:
                                variant_config.nt_sample_dict[variant][
                                    file.split(".")[0]
                                ] = match_dict[variant]
                                variant_config.ref_sample_dict[variant][
                                    file.split(".")[0]
                                ] = ref_dict[variant]
                                variant_config.nt_outfiles[variant].write(
                                    file.split(".")[0] + "\n",
                                )
                                variant_config.nt_outfiles[variant].write(
                                    "\t".join(
                                        PosLine.split("\t")[
                                            : int(variant_config.nt_call_dict["Total"])
                                            + 1
                                        ],
                                    )
                                    + "\tVarient\tAbundance\n",
                                )
                                for line in match_lines[variant]:
                                    splitline = line.split("\t")
                                    variant_config.nt_outfiles[variant].write(
                                        "\t".join(
                                            splitline[
                                                : int(
                                                    variant_config.nt_call_dict[
                                                        "Total"
                                                    ],
                                                )
                                                + 1
                                            ],
                                        ),
                                    )
                                    for subvar in var.split("and"):
                                        if subvar[1:-1] == splitline[0]:
                                            variant_config.nt_outfiles[variant].write(
                                                f"\t{subvar}\t{match_dict[variant][subvar]}",
                                            )
                                    variant_config.nt_outfiles[variant].write("\n")
                                try:
                                    variant_config.nt_SRAs[variant].append(
                                        file.split(".")[0],
                                    )
                                except:
                                    variant_config.nt_SRAs[variant] = [
                                        (file.split(".")[0]),
                                    ]
                            except:
                                pass

                in_file.close()

            if file.endswith("_chim_rm.tsv") or file.endswith("_covar_deconv.tsv"):
                if variant_config.search_omis == 1 or variant_config.search_ref == 1:
                    samp_name = ""
                    if file.endswith("_covar_deconv.tsv"):
                        samp_name = file.strip("_covar_deconv.tsv")
                    else:
                        samp_name = "_".join(file.split("_")[:-3])
                    ref_lines = []
                    amp = ""
                    if "NTD" in file:
                        amp = "NTD"
                    elif "S1S2" in file:
                        amp = "S1S2"
                    elif "RBD" in file:
                        amp = "RBD"
                    if amp:
                        in_file = open(os.path.join(subdir, file))
                        counts = {
                            "1": [0, 0, 0],
                            "2": [0, 0, 0],
                            "3": [0, 0, 0],
                        }  # counts['2'][1] == abundance of seq that had 1 mismatches to BA.2
                        matched = 0
                        for line in in_file:
                            splitline = line.strip("\n\r").split("\t")
                            try:
                                splitline[2]
                            except:
                                pass
                            else:
                                if splitline[1] != "Count":
                                    split_seq = splitline[0].split(" ")
                                    if variant_config.search_omis == 1:
                                        if amp != "S1S2":
                                            for sub in variant_config.omis:
                                                mismatches = 0
                                                for PM in variant_config.omis[sub][amp]:
                                                    if "not" in PM:
                                                        if (
                                                            PM.split(" ")[1]
                                                            in splitline[0]
                                                        ):
                                                            mismatches += 1
                                                    elif PM not in splitline[0]:
                                                        mismatches += 1
                                                if mismatches < 3:
                                                    matched = 1
                                                    counts[sub][mismatches] += float(
                                                        splitline[2],
                                                    )
                                    if variant_config.search_ref == 1:
                                        # if float(splitline[2]) >= .01: # and not '/NY/' in subdir:
                                        if (
                                            amp == "S1S2"
                                            and splitline[0] == "1841G(D614G)"
                                        ) or splitline[0] == "Reference":
                                            ref_lines.append(line)
                    if ref_lines:
                        ref_fh.write(os.path.join(subdir, file) + "\n")
                        for line in ref_lines:
                            ref_fh.write(line)
                        ref_fh.write("\n")

                        # print(counts)
                        if matched == 1:
                            try:
                                variant_config.omi_matches[samp_name]
                            except:
                                variant_config.omi_matches[samp_name] = counts
                            else:
                                for subvar in variant_config.omi_matches[samp_name]:
                                    for i in range(3):
                                        variant_config.omi_matches[samp_name][subvar][
                                            i
                                        ] += counts[subvar][i]

                    in_file.close()

if variant_config.search_ref == 1:
    ref_fh.close()

if variant_config.multi_var:
    for variant in variant_config.multivar_fhs:
        variant_config.multivar_fhs[variant].close()
if variant_config.nt_mult_outfiles:
    for variant in variant_config.nt_mult_outfiles:
        variant_config.nt_mult_outfiles[variant].close()
if variant_config.nt_call_variants:
    for variant in variant_config.nt_call_variants:
        variant_config.nt_outfiles[variant].close()
if variant_config.aa_variants:
    for variant in variant_config.aa_variants:
        variant_config.aa_outfiles[variant].close()
if variant_config.nt_SRAs:
    for variant in variant_config.nt_SRAs:
        print(variant + " " + str(len(variant_config.nt_SRAs[variant])))
        # for SRA in nt_SRAs[variant]:
        # print(SRA)

if variant_config.pos_sample_dict:
    fh_POS_table = open(samps + "_RefNTtable.tsv", "w")
    fh_POS_table.write("Position")
    for sample in variant_config.pos_sample_dict:
        fh_POS_table.write(f"\t{sample}")
    fh_POS_table.write("\n")
    for position in variant_config.nt_positions:
        fh_POS_table.write(variant_config.ref_nt_dict[position] + position)
        for sample in variant_config.pos_sample_dict:
            try:
                fh_POS_table.write(
                    f"\t{variant_config.pos_sample_dict[sample][position]:04f}",
                )
            except:
                fh_POS_table.write("\t")
        fh_POS_table.write("\n")

    fh_POS_table.close()

if variant_config.nt_sample_dict:
    i = 2
    for variant in variant_config.nt_sample_dict:
        fh_NT_table = open(samps + "_" + str(i) + "_NTtable.tsv", "w")
        i += 1
        fh_NT_table.write("mut")
        for sample in variant_config.nt_sample_dict[variant]:
            fh_NT_table.write(f"\t{sample}")
        fh_NT_table.write("\n")
        for subvar in variant.split("and"):
            fh_NT_table.write(subvar)
            for sample in variant_config.nt_sample_dict[variant]:
                try:
                    fh_NT_table.write(
                        f"\t{variant_config.nt_sample_dict[variant][sample][subvar]:04f}",
                    )
                except:
                    fh_NT_table.write("\t")
            fh_NT_table.write("\n")

        fh_NT_table.close()

if variant_config.ref_sample_dict:
    i = 2
    for variant in variant_config.ref_sample_dict:
        fh_Ref_table = open(samps + "_" + str(i) + "_Reftable.tsv", "w")
        i += 1
        fh_Ref_table.write("mut")
        for sample in variant_config.ref_sample_dict[variant]:
            fh_Ref_table.write(f"\t{sample}")
        fh_Ref_table.write("\n")
        for subvar in variant.split("and"):
            fh_Ref_table.write(subvar[:-1] + subvar[0])
            for sample in variant_config.ref_sample_dict[variant]:
                try:
                    fh_Ref_table.write(
                        f"\t{variant_config.ref_sample_dict[variant][sample][subvar]:04f}",
                    )
                except:
                    fh_Ref_table.write("\t")
            fh_Ref_table.write("\n")

        fh_Ref_table.close()


if variant_config.omi_matches:
    Omi_out_fh = open(samps + "_variant_config.omis.tsv", "w")
    sorted_samps = sorted(variant_config.omi_matches.keys())
    for subvar in variant_config.omis:
        for domain in variant_config.omis[subvar]:
            Omi_out_fh.write(
                f"BA.{subvar}\t{domain}\t{', '.join(variant_config.omis[subvar][domain])}",
            )
            Omi_out_fh.write("\n")
    Omi_out_fh.write("\n")
    # Omi_out_fh.write(f"\t\t[Perfect match, 1 mismatch, 2 mismatch]")
    for key in sorted_samps:
        Omi_out_fh.write(key + "\t\tPerfect match\t1 mismatch\t2 mismatch\n")
        for subvar in variant_config.omi_matches[key]:
            Omi_out_fh.write(f"\tBA.{subvar}")
            for val in variant_config.omi_matches[key][subvar]:
                Omi_out_fh.write(f"\t{val}")
            Omi_out_fh.write("\n")
        Omi_out_fh.write("\n")
    Omi_out_fh.close()
