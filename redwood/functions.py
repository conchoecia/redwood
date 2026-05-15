#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Derived from pauvre
# Copyright (c) 2016-2020 Darrin T. Schultz.
# Licensed under the GNU General Public License v3 or later.

import os
import time

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import warnings


def print_images(base=None, image_formats=None, dpi=600, transparent=False, no_timestamp=False, **kwargs):
    """Save the active matplotlib figure in one or more formats."""
    if base is None:
        base = kwargs.pop("base_output_name", None)
    if base is None:
        raise ValueError("print_images requires a base output name")
    if image_formats is None:
        image_formats = ["png"]
    for fmt in image_formats:
        if no_timestamp:
            out_name = "{0}.{1}".format(base, fmt)
        else:
            out_name = "{0}_{1}.{2}".format(base, timestamp(), fmt)
        try:
            if fmt == "png":
                plt.savefig(out_name, dpi=dpi, transparent=transparent)
            else:
                plt.savefig(out_name, format=fmt, transparent=transparent)
        except PermissionError:
            print(
                "You do not have permission to save redwood plots to this "
                "directory. Try changing the directory and running the command again."
            )


class GFFParse:
    """Parse the subset of GFF3 needed for redwood circular annotation tracks."""

    def __init__(self, filename, stop_codons=None, species=None):
        self.filename = filename
        self.samplename = os.path.splitext(os.path.basename(filename))[0]
        self.species = species
        gffnames = [
            "sequence",
            "source",
            "featType",
            "start",
            "stop",
            "dunno1",
            "strand",
            "dunno2",
            "tags",
        ]
        self.features = pd.read_csv(self.filename, comment="#", sep="\t", names=gffnames)
        self.features["name"] = self.features["tags"].apply(self._get_name)
        self.features.drop("dunno1", axis=1, inplace=True)
        self.features.drop("dunno2", axis=1, inplace=True)
        self.features.reset_index(inplace=True, drop=True)
        self._check_triplets()
        self.features.sort_values(by="start", ascending=True, inplace=True)
        if stop_codons:
            strip_codons = ["gene", "CDS"]
            self.features.loc[
                (self.features["featType"].isin(strip_codons)) & (self.features["strand"] == "+"),
                "stop",
            ] = (
                self.features.loc[
                    (self.features["featType"].isin(strip_codons))
                    & (self.features["strand"] == "+"),
                    "stop",
                ]
                - 3
            )
            self.features.loc[
                (self.features["featType"].isin(strip_codons)) & (self.features["strand"] == "-"),
                "start",
            ] = (
                self.features.loc[
                    (self.features["featType"].isin(strip_codons))
                    & (self.features["strand"] == "-"),
                    "start",
                ]
                + 3
            )
        self.features["center"] = self.features["start"] + (
            (self.features["stop"] - self.features["start"]) / 2
        )
        self.features["width"] = abs(self.features["stop"] - self.features["start"]) + 1
        self.features["lmost"] = self.features.apply(self._determine_lmost, axis=1)
        self.features["rmost"] = self.features.apply(self._determine_rmost, axis=1)
        self.features["track"] = 0
        circular_regions = self.features.loc[
            self.features["tags"].fillna("").str.contains("Is_circular=true", regex=False),
            "stop",
        ]
        if len(circular_regions) < 1:
            raise IOError(
                'The GFF file needs a region tag ending in "Is_circular=true", '
                "from 1 to the number of bases in the circular genome."
            )
        self.seqlen = int(circular_regions.iloc[0])
        self.features.reset_index(inplace=True, drop=True)

    def _check_triplets(self):
        genes_cdss = self.features.query("featType == 'CDS' or featType == 'gene'")
        not_trips = genes_cdss.loc[((abs(genes_cdss["stop"] - genes_cdss["start"]) + 1) % 3) > 0]
        if len(not_trips) > 0:
            warnings.warn(
                "There are CDS and gene entries that are not divisible by three\n"
                + str(not_trips),
                SyntaxWarning,
            )

    def _get_name(self, tag_value):
        try:
            if ";" in tag_value:
                return tag_value[5:].split(";")[0]
            return tag_value[5:].split()[0]
        except Exception:
            return tag_value

    def _determine_lmost(self, row):
        return row["start"] if row["start"] < row["stop"] else row["stop"]

    def _determine_rmost(self, row):
        return row["stop"] if row["start"] < row["stop"] else row["start"]


def timestamp():
    """Return the current time in YYYYMMDD_HHMMSS format."""
    return time.strftime("%Y%m%d_%H%M%S")
