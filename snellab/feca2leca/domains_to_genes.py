#!/usr/bin/env python3

# Written by John van Dam
# Adapted by Julian Vosseberg

import sys

# File with both tree sequences and representing sequences
ogs_file_location = sys.argv[1]
# Threshold for the overlap percentage between OGs to merge them
overlap_threshold = float(sys.argv[2])

ogs = {}
with open(ogs_file_location) as ogs_file:
    ogs_file.readline()
    for line in ogs_file:
        line = line.rstrip()
        fields = line.split('\t')
        og = fields[0] + '_' + fields[1]
        if len(fields) == 4:
            seq_ids = fields[2].split(',') + fields[3].split(',')
        elif len(fields) == 5:
            seq_ids = fields[3].split(',') + fields[4].split(',')
        else:
            sys.exit('Error: number of tab-separated fields incorrect.')
        seq_ids = set([seq_id if not '_' in seq_id else seq_id[:seq_id.find('_')] for seq_id in seq_ids])
        ogs[og] = seq_ids

old_set = dict()
new_set = ogs.copy()
ogs_merged = 0

while old_set != new_set:  # When both dicts are identical comp() returns 0, which is quite convenient in this case
    old_set = new_set.copy()
    sorted_old_ogs = sorted(old_set, key = lambda x: len(old_set[x]))
    merged = False
    try:
        for index,og_id in enumerate(sorted_old_ogs):
            for larger_og_id in sorted_old_ogs[index:]:  # Not doing index+1 to avoid self matching because at the end this would throw an exception
                if larger_og_id == og_id:
                    continue
                if old_set[og_id].issubset(old_set[larger_og_id]):
                    # og is subset of another, so merge in the new_set
                    sys.stderr.write("{} ({}) is subset of {} ({})\n".format(og_id,len(old_set[og_id]),larger_og_id,len(old_set[larger_og_id])))
                    del new_set[og_id]  # First attempt to delete the old one, so that on KeyError the key does not already get updated.
                    new_set[larger_og_id + '|' + og_id] = new_set.pop(larger_og_id)  # Updating the larger set's key to include the redundant OG id.
                    merged = True
                    break
                elif old_set[og_id] & old_set[larger_og_id]:
                    size_x_set = len(old_set[og_id])
                    size_y_set = len(old_set[larger_og_id])
                    size_overlap = len(old_set[og_id].intersection(old_set[larger_og_id]))
                    pcnt_overlap = size_overlap/size_x_set*100
                    if pcnt_overlap >= overlap_threshold:
                        sys.stderr.write("{} ({}) shares overlap with {} ({}) of size {} ({}%)\n".format(og_id,size_x_set,larger_og_id,size_y_set,size_overlap,pcnt_overlap))
                        # The situation now is different in that I try to add new_set[og_id] to the larger og, so here we will throw the key error before assigning anyway
                        # And we can only delete the old og after we've added it, so here the del statement comes second.
                        new_set[larger_og_id + '|' + og_id] = new_set.pop(larger_og_id).union(new_set[og_id])  # Updating the larger set's key to include the redundant OG id.
                        del new_set[og_id]
                        merged = True
                        break
            if merged:
                ogs_merged += 1
                break
    except KeyError:
        sys.stderr.write("Could not find og_id in updated orthology dictionary, likely because the key was already updated before and removed. Restarting procedure with updated orthology dictionary.\n")

for ogid in sorted(new_set.keys(),key=lambda x: len(new_set[x])):
    print("{}: {}".format(ogid," ".join([str(x) for x in list(new_set[ogid])])))

sys.stderr.write(f"{ogs_merged} OGs merged in total.\n")
