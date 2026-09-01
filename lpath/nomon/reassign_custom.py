import numpy


def reassign_custom(data, pathways, dictionary, assign_file=None):
    """
    Reassign each frame's state_id to its shared-space cluster label.

    The cluster labels were computed jointly across both Eg5 WE runs
    (see ../cluster_shared.py) and stored in each west.h5 as
    'auxdata/labels', then carried through the extract step via
    ``-a labels``.  Here we overwrite the w_assign state_id (col 2) with
    that cluster label so the matching/dendrogram is built on the shared
    cluster space, making the no-monastrol and with-monastrol pathways
    directly comparable.
    """
    # reassign states to be the shared cluster IDs
    for idx, val in enumerate(data):          # each successful pathway
        val_arr = numpy.asarray(val)
        for idx2, val2 in enumerate(val_arr):  # each frame of the pathway
            val2[2] = int(val2[-3])            # renumber state_id with aux label
            pathways[idx, idx2] = val2

    # map each shared cluster id -> single-character state string
    labels = numpy.unique(numpy.concatenate(
        [numpy.asarray(v)[:, -3].astype(int) for v in data]))
    dictionary = {int(l): str(int(l)) for l in labels}
    # last entry is the "unknown" state
    dictionary[max(dictionary) + 1] = "?"

    return dictionary
