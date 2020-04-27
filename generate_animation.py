#!/usr/bin/python
#-*-coding:utf-8-*-

import sys, os
import pandas as pd
import numpy as np
import datetime
import math
from matplotlib.colors import to_rgb
from matplotlib.animation import ArtistAnimation
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import copy

MIN_LENGTH = 29500
MAX_LENGTH = 30000
METADATA = "metadata.csv"
INPUT_NAME = "sequences.fasta"
OUTPUT_NAME = 'nuc_sequence.csv'
ANIMATION_NAME = 'sarscov2_genome.mp4'
base_rgb_map = {"A":to_rgb("b"), "C":to_rgb("r"), "G": to_rgb("g"), "T": to_rgb("y"), "N" : to_rgb("k")}

def valid_sequence(string):
    if len(string) < MIN_LENGTH or len(string) > MAX_LENGTH:
        return False
    valid_strings = {"T", "A", "G", "C"}
    return all([c in valid_strings for c in string])

def create_rgb_array(full_str):
    rows = 200
    columns = 150
    char_array = np.asarray(list(full_str))

    nones = (rows * columns - len(char_array)) * ['N']

    full_array = np.concatenate((char_array, nones)).reshape((rows, columns))

    # n corresponds to "None"
    rgbs = np.empty((rows, columns, 3))
    for i in range(rows):
        for j in range(columns):
            rgbs[i, j, :] = base_rgb_map[full_array[i, j]]
    return rgbs

class Dataset(object):
    def __init__(self, accession, sequence, date, location):
        """
        Args:
            accession (string): id of the accession
            sequence (string): nucleotide bases, represented by one of {A, T, G, C}
            length (int): length of the sequence
            date (datetime): the date of the dataset
            location (string): 
        """
        self.accession = accession
        self.date = date
        self.location = location
        self.length = len(sequence)
        self._sequence = sequence
        self._rgb = None

    @property
    def sequence(self):
        return self._sequence
    
    @sequence.setter
    def sequence(self, new_sequence):
        self._sequence = new_sequence
        self.length = len(sequence)
    
    def frame(self, ax):
        im = plt.imshow(create_rgb_array(self._sequence), animated = True)
        ttl = plt.text(0.5, 1.01, str(self.accession) + " " + self.date.strftime("%m/%d/%y") , horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes, fontsize='large')
        return [im, ttl]



"""
Read the metadata to construct a mapping between Assencion IDs and pieces of data.
"""
# Read the metadata
meta = pd.read_csv(METADATA)

# Construct a mapping between Assencion ID and properties.
data_map = {}
date_number = {} # Keep track of how nay in each date so far.
for i, row in meta.iterrows():
    # row[1] takes the form of 2020-01-29T00:00:00Z, so we
    # convert it to a datetime using the first 9 characters.
    stripped_date = row[1][0:10].replace("-","")
    data_map[row[0]] = dict()
    data_datetime = datetime.datetime.strptime(stripped_date, "%Y%m%d")
    data_map[row[0]]["date"] = data_datetime
    if data_datetime not in date_number.keys():
        date_number[data_datetime] = 0
    data_map[row[0]]["number"] = date_number[data_datetime]
    date_number[data_datetime] += 1
    if not math.isnan(row[3]):
        data_map[row[0]]["length"] = int(row[3])
    else:
        data_map[row[0]]["length"] = None
    data_map[row[0]]["location"] = str(row[4])



"""
Read the data stored in the fasta file.
"""
file = open(INPUT_NAME, 'r')
lines_i = file.readlines()

# datasets is a list
datasets = []
sequence = None

for l in lines_i:
    if ">" in l:
        if not sequence == None:
            if valid_sequence(sequence):
                accession = (l.strip().split(" ")[0]).replace(">", "")
                datasets.append(Dataset(accession, sequence, data_map[accession]["date"], data_map[accession]["location"]))
        sequence = ''
    else:
        sequence += l.strip()

# Filter here
valid_datasets = [dataset for dataset in datasets if "USA" in dataset.location]
# valid_datasets = datasets

"""
Many of the sequences are truncated either at the beginning or the end.
The idea here is to use the longest sequence, find some substring at the beginning, and then
start every single dataset from that substring. Doing so requires us to start a bit further
into the longest sequence.
"""
def start_from_substring(dataset, substring, max_idx = 200):
    idx_start = dataset.sequence.find(starting_substring, 0, max_idx)
    new_sequence = copy.copy(dataset.sequence[idx_start:])
    dataset.sequence = new_sequence
    return dataset

start_idx = 100
substring_len = 10
longest_sequence = max(valid_datasets, key= lambda dataset: dataset.length)
starting_substring = longest_sequence.sequence[start_idx : start_idx + substring_len]
valid_datasets = [start_from_substring(dataset, starting_substring) for dataset in valid_datasets]
valid_datasets = [dataset for dataset in valid_datasets if len(dataset.sequence) > MIN_LENGTH]



# Sort the accessions
sorted_datasets = sorted(valid_datasets, key=lambda data : data.date)


# Create the animation.
fig, ax = plt.subplots()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
frames = [d.frame(ax) for d in sorted_datasets]
base_pairs = ["A", "C", "T", "G"]

patches = [ mpatches.Patch(color=base_rgb_map[bp], label=bp) for bp in base_pairs]
plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
ani = ArtistAnimation(fig, frames, interval=100, blit=False)
ani.save(ANIMATION_NAME)
plt.show()