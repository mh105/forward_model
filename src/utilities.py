import numpy as np


def shrink_source_space(scaling_factor=1, src=None):
    # left hemisphere and right hemisphere
    orig_src_0 = src[0]['rr']
    orig_src_1 = src[1]['rr']
    source_combined_location = np.concatenate((orig_src_0, orig_src_1), axis=0)
    coordinate_min, coordinate_max = np.amin(source_combined_location, axis=0), np.amax(source_combined_location, axis=0)
    central_point = np.true_divide(coordinate_min+coordinate_max, 2.0)
    shrinked_source_space_1 = (source_combined_location-central_point)*scaling_factor
    shrinked_source_space = shrinked_source_space_1 + central_point

    src_new_0 = shrinked_source_space[0:orig_src_0.shape[0], :]
    src_new_1 = shrinked_source_space[orig_src_0.shape[0]:, :]

    src[0]['rr'] = src_new_0
    src[1]['rr'] = src_new_1
    return src


def get_vert_atlas_lobe_info(labels, src, atlas_templates):
    # """Generate extract_label_time_course."""
    # if len(src) > 2:
    #     if src[0]['type'] != 'surf' or src[1]['type'] != 'surf':
    #         raise ValueError('The first 2 source spaces have to be surf type')
    #     if any(np.any(s['type'] != 'vol') for s in src[2:]):
    #         raise ValueError('source spaces have to be of vol type')
    #
    #     n_aparc = len(labels)
    #     n_aseg = len(src[2:])
    #     n_labels = n_aparc + n_aseg
    # else:
    #     n_labels = len(labels)

    # get vertices from source space, they have to be the same as in the stcs
    vertno = [s['vertno'] for s in src]
    nvert = [len(vn) for vn in vertno]
#    vert_pos = [s.pos for s in src]

    # do the initialization
    atlas_vertidx = list()
    atlas_names = list()
    atlas_lobe = list()

    for label in labels:
        if label.hemi == 'both':
            # handle BiHemiLabel
            sub_labels = [label.lh, label.rh]
        else:
            sub_labels = [label]
        this_vertidx = list()
        for slabel in sub_labels:
            if slabel.hemi == 'lh':
                this_vertno = np.intersect1d(vertno[0], slabel.vertices)  # find all the vertices in that atlas
                vertidx = np.searchsorted(vertno[0], this_vertno)  # find the correct order.
            elif slabel.hemi == 'rh':
                this_vertno = np.intersect1d(vertno[1], slabel.vertices)
                vertidx = nvert[0] + np.searchsorted(vertno[1], this_vertno)
            else:
                raise ValueError('label %s has invalid hemi' % label.name)
            this_vertidx.append(vertidx)

        # convert it to an array
        this_vertidx = np.concatenate(this_vertidx)

        atlas_vertidx.append(this_vertidx)
        atlas_names.append(label.name)

        atlas_lobe.append(atlas_templates.loc[label.name, 'lobes'])

    atlas_info = dict()
    atlas_info['atlas_vertidx'] = atlas_vertidx
    atlas_info['atlas_names'] = atlas_names
    atlas_info['atlas_lobe'] = atlas_lobe
    return atlas_info


def atlas_matlab_update_idx(atlas_info):
    copy_vertidx = atlas_info['atlas_vertidx']
    for i in range(len(copy_vertidx)):
        copy_vertidx[i] = copy_vertidx[i] + 1
    # atlas_info['atlas_vertidx'] = copy_vertidx  # unnecessary since dictionary is mutable data structure
    return atlas_info

