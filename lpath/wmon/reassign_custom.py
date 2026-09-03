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

    IMPORTANT: the shared cluster ids that actually occur are non-contiguous
    (e.g. {0, 1, 3, 4}).  lpath's ``calc_dist`` strips "unknown" frames by
    *index*, not value -- ``seq[seq < len(dictionary) - 1]`` -- and treats the
    highest dictionary key as the unknown state.  With non-contiguous ids the
    top real cluster (4) looks like "unknown" and is silently dropped from
    every path string (here cluster 4 is the fully-unbound terminal state,
    visited by nearly every pathway).  To avoid that we remap the occurring
    labels to a CONTIGUOUS 0..K-1 range and place the "unknown" sentinel at K,
    so every real cluster survives the filter.
    """
    # shared-cluster labels that actually occur, mapped to contiguous 0..K-1
    labels = numpy.unique(numpy.concatenate(
        [numpy.asarray(v)[:, -3].astype(int) for v in data]))
    remap = {int(l): i for i, l in enumerate(labels)}

    # renumber state_id (col 2) with the remapped shared-cluster label
    for idx, val in enumerate(data):          # each successful pathway
        val_arr = numpy.asarray(val)
        for idx2, val2 in enumerate(val_arr):  # each frame of the pathway
            val2[2] = remap[int(val2[-3])]
            pathways[idx, idx2] = val2

    # map each contiguous state id -> single-character state string; the
    # character keeps the ORIGINAL shared-cluster id for readability.
    dictionary = {i: str(int(l)) for l, i in remap.items()}
    # last entry is the "unknown" state (now truly the highest key)
    dictionary[len(remap)] = "?"

    return dictionary
