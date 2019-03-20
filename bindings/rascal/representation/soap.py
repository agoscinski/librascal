import json

from ..utils import get_strict_neighbourlist
from ..lib import RepresentationManager, FeatureManager
from .base import RepresentationFactory

class SOAP(object):

    """
    Computes a SOAP representation, e.g. power spectrum etc.

    Hyperparameters
    ----------
    interaction_cutoff : float
        Maximum pairwise distance for atoms to be considered in
        expansion

    cutoff_smooth_width : float
        The distance over which the the interaction is smoothed to zero

    max_radial : int
        Number of radial basis functions

    max_angular : int
        Highest angular momentum number (l) in the expansion

    n_species : int
        Number of species to be considered separately

    gaussian_sigma_type : str
        How the Gaussian atom sigmas (smearing widths) are allowed to
        vary -- fixed ('Constant'), by species ('PerSpecies'), or by
        distance from the central atom ('Radial').

    gaussian_sigma_constant : float
        Specifies the atomic Gaussian widths, in the case where they're
        fixed.

    soap_type : string
        Specifies the type of representation to be computed
        (power spectrum etc.).

    Methods
    -------
    transform(frames)
        Compute the representation for a list of ase.Atoms object.


    .. [soap] Bartók, Kondor, and Csányi, "On representing chemical
        environments", Phys. Rev. B. 87(18), p. 184115
        http://link.aps.org/doi/10.1103/PhysRevB.87.184115

    """

    def __init__(self, interaction_cutoff, cutoff_smooth_width,
                 max_radial, max_angular, gaussian_sigma_type,
                 gaussian_sigma_constant=0., n_species=1, 
                 soap_type="PowerSpectrum"):
        """Construct a SphericalExpansion representation

        Required arguments are all the hyperparameters named in the
        class documentation
        """
        self.name = 'soap'
        self.hypers = dict()
        self.update_hyperparameters(
            interaction_cutoff=interaction_cutoff,
            cutoff_smooth_width=cutoff_smooth_width,
            max_radial=max_radial, max_angular=max_angular,
            gaussian_sigma_type=gaussian_sigma_type,
            gaussian_sigma_constant=gaussian_sigma_constant,
            n_species=n_species,
            soap_type=soap_type)

        if soap_type == "RadialSpectrum":
            self.update_hyperparameters(max_angular=0)

    def update_hyperparameters(self, **hypers):
        """Store the given dict of hyperparameters

        Also updates the internal json-like representation

        """
        allowed_keys = {'interaction_cutoff', 'cutoff_smooth_width',
                        'max_radial', 'max_angular', 'gaussian_sigma_type',
                        'gaussian_sigma_constant', 'n_species', 'soap_type'}
        hypers_clean = {key: hypers[key] for key in hypers
                        if key in allowed_keys}
        self.hypers.update(hypers_clean)
        return

    def transform(self, frames):
        """Compute the representation.

        Parameters
        ----------
        frames : list(ase.Atoms)
            List of atomic structures.

        Returns
        -------
        FeatureManager.Dense_double
            Object containing the representation

        """
        n_frames = len(frames)
        managers = list(map(
            get_strict_neighbourlist,
            frames, [self.hypers['interaction_cutoff'], ] * n_frames))
        hypers_str = json.dumps(self.hypers)
        n_features = self.get_num_coefficients()
        features = FeatureManager.BlockSparse_double(n_features, hypers_str)
        cms = map(RepresentationFactory(self.name),
                  managers, [hypers_str, ] * n_frames)
        for cm in cms:
            cm.compute()
            features.append(cm)

        return features

    def get_num_coefficients(self):
        """Return the number of coefficients in the representation

        (this is the descriptor size per atomic centre)

        """
        if self.hypers['soap_type'] == 'RadialSpectrum':
            return (self.hypers['n_species']*self.hypers['max_radial'])
        if self.hypers['soap_type'] == 'PowerSpectrum':
            return (self.hypers['n_species']**2*self.hypers['max_radial']**2
                    *(self.hypers['max_angular'] + 1))
        else:
            raise RuntimeError('Only soap_type = RadialSpectrum || '
                               'PowerSpectrum implemented for now')