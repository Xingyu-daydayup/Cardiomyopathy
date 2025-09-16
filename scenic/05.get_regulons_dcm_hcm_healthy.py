import pickle
import numpy as np
from pyscenic.utils import load_motifs
import os
from pyscenic.transform import df2regulons
import operator as op
from cytoolz import compose
def derive_regulons():
    # Load enriched motifs.
    motifs = load_motifs('reg_dcm_hcm_healthy.csv')
    motifs.columns = motifs.columns.droplevel(0)

    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f

    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motifs for the modules that were created using the method 'weight>50.0%' (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.
    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains('weight>50.0%')), motifs.Context), dtype=bool) & \
        # np.fromiter(map(contains('hg19-tss-centered-10kb-10species.mc9nr',
        #                          'hg19-500bp-upstream-10species.mc9nr',
        #                          'hg19-tss-centered-5kb-10species.mc9nr'), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains('activating'), motifs.Context), dtype=bool)]

    # We build regulons only using enriched motifs with a NES of 3.0 or higher; we take only directly annotated TFs or TF annotated
    # for an orthologous gene into account; and we only keep regulons with at least 10 genes.
    regulons = list(filter(lambda r: len(r) >= 10, df2regulons(motifs[(motifs['NES'] >= 3.0)
                                                                      & ((motifs[
                                                                              'Annotation'] == 'gene is directly annotated')
                                                                         | (motifs['Annotation'].str.startswith(
                'gene is orthologous to')
                                                                            & motifs['Annotation'].str.endswith(
                        'which is directly annotated for motif')))
                                                                      ])))

    # Rename regulons, i.e. remove suffix.
    regulons = list(map(lambda r: r.rename(r.transcription_factor), regulons))

    # Pickle these regulons.
    with open('./dcm_hcm_healthy.regulons.dat', 'wb') as f:
        pickle.dump(regulons, f)
derive_regulons()
print('success')
