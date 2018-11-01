import scipy.stats
import multiprocessing, time

import psutil
import numpy as np


def compute_all_scores_multi(spectra_list, gcf_list, strain_list, scoring_function, do_random=True, cpus=8):
    # TODO get num CPUs available from psutil
    spectra_part_list = np.array_split(spectra_list, cpus)

    q = multiprocessing.Queue()
    procs = []
    for i in range(cpus):
        if len(spectra_part_list[i]) == 0:
            continue

        p = multiprocessing.Process(target=compute_all_scores, args=(spectra_part_list[i], gcf_list, strain_list, scoring_function, do_random, q, i))
        procs.append(p)

    for p in procs:
        p.start()

    m_scores = {}

    t = time.time()
    num_finished = 0
    while num_finished < len(procs):
        data = q.get()
        m_scores.update(data)
        num_finished += 1

    print(('Total time: {:.1f}, {:.1f}/s'.format(time.time() - t, len(spectra_list) / (time.time() - t))))
    for p in procs:
        p.join()

    return m_scores

def compute_all_scores(spectra_list,gcf_list,strain_list,scoring_function,do_random = True, q=None, cpu_aff=None):
    m_scores = {}
    best = 0
    best_random = 0

    if cpu_aff is not None:
        psutil.Process().cpu_affinity([cpu_aff])

    for i,spectrum in enumerate(spectra_list):
        m_scores[spectrum] = {}
        if i % 100 == 0:
            print(("Done {} of {}".format(i,len(spectra_list))))
        for gcf in gcf_list:
            s,metadata = scoring_function(spectrum,gcf,strain_list)
            if do_random:
                s_random,_ = scoring_function(spectrum.random_spectrum,gcf.random_gcf,strain_list)
            else:
                s_random = None
            m_scores[spectrum][gcf] = (s,s_random,metadata)


    if q is not None:
        q.put(m_scores)

    return m_scores

def metcalf_scoring(spectral_like,gcf_like,strains,both = 10,met_not_gcf = -10,gcf_not_met = 0,neither = 1):
    cum_score = 0
    shared_strains = set()
    for strain in strains:
        in_spec = spectral_like.has_strain(strain)
        in_gcf = gcf_like.has_strain(strain)
        if in_spec and in_gcf:
            cum_score += both
            shared_strains.add(strain)
        if in_spec and not in_gcf:
            cum_score += met_not_gcf
        if in_gcf and not in_spec:
            cum_score += gcf_not_met
        if not in_gcf and not in_spec:
            cum_score += 1
    return cum_score,shared_strains


def hg_scoring(spectral_like, gcf_like, strains):
    spectral_count = 0
    gcf_count = 0
    overlap_count = 0

    for strain in strains:
        if spectral_like.has_strain(strain):
            spectral_count += 1
        if gcf_like.has_strain(strain):
            gcf_count += 1
        if spectral_like.has_strain(strain) and gcf_like.has_strain(strain):
            overlap_count += 1

    pos_in_sample = overlap_count
    N = spectral_count
    n = gcf_count
    M = len(strains)

    r = scipy.stats.hypergeom.sf(pos_in_sample, M, n, N, 1)
    return r, None



def name_scoring(spectral_like,gcf_like,mibig_map):
	score = 0
	metadata = None
	if len(spectral_like.annotations) == 0:
		print("No annotations")
		return None,None
	mibig_bgcs = gcf_like.get_mibig_bgcs()
	if len(mibig_bgcs) == 0:
		print("no mibig")
		return None,None
	for annotation in spectral_like.annotations:
		for mibig in mibig_bgcs:
			short_mibig = mibig.name.split('.')[0]
			if short_mibig in mibig_map:
				m = match(annotation,mibig_map[short_mibig])
				if m:
					metadata = m
					score += 100
	return (score,metadata)

def match(spectral_annotation,mibig_name):
	metadata = None
	name,source = spectral_annotation
	for m_name in mibig_name:
		if name.lower() == m_name.split()[0].lower():
			print(name,m_name)
			metadata = (name,m_name)
			return metadata
	return False

def knownclusterblast_scoring(spectral_like,gcf_like,mibig_map):
    score = 0
    metadata = None
    if len(spectral_like.annotations) == 0:
        print("No annotations")
        return None,None
    kcb = []
    for bgc in gcf_like.bgc_list:
        these = bgc.known_cluster_blast
        # if hasattr(bgc,'metadata'):
        #     these = bgc.metadata.get('knownclusterblast',None)
        if these:
            for mibig,score in these:
                kcb.append((mibig,score))
    if len(kcb) == 0:
        return None,None
    total_score = 0
    for annotation in spectral_like.annotations:
        for mibig,score in kcb:
            short_mibig = mibig.split('_')[0]
            if short_mibig in mibig_map:
                m = match(annotation,mibig_map[short_mibig])
                if m:
                    metadata = m
                    total_score += int(score)
                    print(m)
    return total_score,metadata


def aa_scoring(spectrum, gcf_like):
    """
    Check for the prescence of AA mass shifts in the spectrum
    """
    tol = 0.01
    from metabolomics import read_aa_losses
    aa_loss_file = 'aa_residues.csv'
    aa_losses = read_aa_losses(aa_loss_file)

    p = 1.0
    for aa, aa_prob in gcf_like.aa_predictions:
        if aa_prob < 0.2 or aa_prob > 0.8:
            if aa in aa_losses:
                mass_iso, mass_avg = aa_losses[aa]
                found_losses = spectrum.has_loss(mass_iso, tol)
                if len(found_losses) > 0:
                    p *= aa_prob
                else:
                    p *= (1 - aa_prob)

    return p

