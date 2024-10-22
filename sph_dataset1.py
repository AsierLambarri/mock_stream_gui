# 
#
# This file contains code generated by A. Lambarri Martinez for the
# Master Thesis in Astrophysics at UCM
#
# Version: Feb 12, 2024 @18:00
#
# SPDX-License-Identifier: GPL-3.0+
#

import yt

def create_sph_dataset(ad, pt, 
                       extra_fields = None, 
                       n_neighbours = 32, 
                       kernel = 'wendland6',
                       use_norm = True
                      ):
    """Creates an sph field using the provided particle field, from the data contained
    in the ds dataset. The smoothing is performed by .add_sph_field. Both the number of
    neighbours and kernel can be selected. The returned field has name 'io'.

    Parameters
    ----------
    ad : yt.data_objects.selection_objects.region.YTRegion
        Original dataset.
    pt : str
        Particle type to be smoothened.
    extra_fields : list, optional
        Fields of the particles that will be added to the smoothened dataset.
        Default is None, which adds only mass, position and velocity of the particles.
        Any extra field passed, will be appended to this list.
    n_neighbours : int, optional
        Number of neighbour particles to include in the smoothing. Default is 32.
    kernel : str, optional
        Kernel used in the smoothing. Default is 'wendland6' (see yt Documentation).
    use_norm : bool, optional
        Value of use_sph_normalization. Default is True.

    Returns
    -------
    ds_sph : yt.dataset
        Dataset containing the sph field. The returned field has name 'io'.
    """
    
    fields = ['particle_mass'] + [f'particle_position_{ax}' for ax in 'xyz'] + [f'particle_velocity_{ax}' for ax in 'xyz']
    if extra_fields is not None:
        for e_field in extra_fields:
            field.append(e_field)
            
    data = {field: ad[pt, field] for field in fields}
    
    ds_sph = yt.load_particles(data, data_source=ad)
    
    ds_sph.use_sph_normalization = use_norm
    
    ds_sph.add_sph_fields(n_neighbors=n_neighbours, kernel=kernel, sph_ptype='io')
    
    return ds_sph








