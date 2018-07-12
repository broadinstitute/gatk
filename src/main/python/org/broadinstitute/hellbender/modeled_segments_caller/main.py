import math
import random
import matplotlib
matplotlib.use('Agg')
import numpy as np
import scipy
from scipy.stats import multivariate_normal
from scipy.optimize import minimize
from sklearn import mixture, cluster
from sklearn.neighbors import KernelDensity
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from sklearn import datasets
from scipy.spatial.distance import pdist
import numpy.ma as ma
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from matplotlib import  ticker
from matplotlib import pyplot as plt
from pymc3 import Normal, Metropolis, sample, MvNormal, Dirichlet, DensityDist, find_MAP, NUTS, Slice
import pymc3 as pm
from pymc3.math import logsumexp
import theano.tensor as tt
from theano.tensor.nlinalg import det
from copy import deepcopy
from typing import List
import logging
import datetime
from statistics import median

RANDOM_SEED = 123


def is_number(s):
    """ Decide if 's' is a number. """
    try:
        float(s)
        return True
    except ValueError:
        return False


class LoadAndSampleCrAndAf:
    """ Class that loads the copy ratio and allele fraction data from the modelFinal.seg file.
        In this data file, both the copy ratio and the allele fraction posterior distribution are
        characterized by the distribution's median, 10th and 90th percentile.
        We attribute a weigth to each segment based on the tightness of their posteriors.
    """
    def __init__(self, filename: str, load_cr: bool = True, load_af: bool=True, do_logging: bool=True,
                 cr_weight_ratio_max=10, af_weight_ratio_max=10, output_log_dir: str= "", output_log_prefix: str= ""):
        """ Inputs:
            - filename: 'modelFinal.seg' file characterizing the posterior distribution of the
              copy ratio and the allele fraction data
            - load_cr: whether to load copy ratio data
            - load_af: whether to load the allele fraction data
            - cr_weight_ratio_max: defines the maximum value that the weights of copy ratio points can take.
              The maximum value is cr_weight_ratio_max * median(copy ratio weights)
            - af_weight_ratio_max: same as the previous one, just using allele fraction instead of copy ratio
        """

        # Start logging
        if do_logging:
            # Set log filename and create log file
            if output_log_dir == "":
                # If no directory is given for logging, we use the input file's directory
                input_filename_ending = filename.split("/")[-1]
                input_dir = filename.split(input_filename_ending)[0]
                self.__output_log_dir=input_dir
            else:
                self.__output_log_dir = output_log_dir
                if output_log_dir[-1] != "/":
                    self.__output_log_dir += "/"
            if output_log_prefix == "":
                # If no prefix is given for the log files, we use the input filename
                log_filename_prefix = filename
                log_filename_extension = log_filename_prefix.split(".")[-1]
                log_filename_prefix = log_filename_prefix.split("." + log_filename_extension)[0]
                log_filename_prefix = log_filename_prefix.split("/")[-1]
                self.__output_log_prefix = log_filename_prefix
            else:
                self.__output_log_prefix = output_log_prefix
            self.__log_filename = self.__output_log_dir + self.__output_log_prefix + ".log"

            logging.basicConfig(level=logging.INFO,
                                datefmt='%Y-%m-%d %H:%M:%S',
                                filename=self.__log_filename,
                                filemode="w")
            self.__logger = logging.getLogger()
            self.__logger.info("* ---- %s ---- *" % str(datetime.datetime.now()))
            self.__logger.info("Initializing class LoadAndSampleCrAndAf")
            self.__logger.info("Input file: %s" % filename)
            self.__logger.info("load_copy_ratio: %s" % load_cr)
            self.__logger.info("load_allele_fraction: %s" % load_af)
        else:
            # In the future, whenever we find that self.__log_filename is an empty string, we don't do logging,
            # otherwise we do.
            self.__log_filename = ""
            self.__output_log_dir = ""
            self.__output_log_prefix = ""
            self.__logger = None

        # Initialize
        np.random.seed(RANDOM_SEED)
        self.__filename = filename
        self.__load_cr = load_cr
        self.__load_af = load_af
        self.__cr_weight_ratio_max = cr_weight_ratio_max
        self.__af_weight_ratio_max = af_weight_ratio_max

        # Load data from file
        if not output_log_prefix == "":
            self.__logger.info("Loading data from file.")

        self.__max_copy_ratio_possible = 40   # we do not load any copy ratio value that is larger than this
        [self.__copy_ratio_median,           # median of the copy ratio posterior for each segment
         self.__copy_ratio_10th_perc,        # 10th percentile ...
         self.__copy_ratio_90th_perc,        # 90th percentile ...
         self.__allele_fraction_median,      # median of the allele fraction posterior
         self.__allele_fraction_10th_perc,   # 10th percentile ...
         self.__allele_fraction_90th_perc,   # 90th percentile ...
         self.__contig,                      # contig (chromosome) of the segment
         self.__segment_start,               # starting position of the segment (bp)
         self.__segment_end,                 # ending position of the segment (bp)
         self.__segment_lengths,             # lengths of the segments
         self.__n_points_cr,                 # number of copy ratio points in each segment
         self.__n_points_af,                 # number of allele fraction in each segment
         self.__n_segments                   # number of segments
         ] = self.__import_copy_ratio_and_allele_fraction()

        # Sample posterior distributions
        # Returns list of copy ratio and allele fraction samples for each segment
        if not output_log_prefix == "":
            self.__logger.info("Determining the weights of the segments.")

        # Determine the weights of each segment based on the posterior
        self.__weights = self.__determine_weights()
        plt.plot(self.__weights)
        plt.show()

        # Auxiliary arrays, used by methods that do not know how to deal with points with weights:
        # we sample the number of points in each segment to be proportional to their weights
        if not output_log_prefix == "":
            self.__logger.info("Sampling data points according to the segment weights.")
        self.__copy_ratio_sampled, self.__allele_fraction_sampled = self.__sample_points()

        if not output_log_prefix == "":
            self.__logger.info("Finished.\n\n\n")

    def get_data(self):
        """ Return data loaded from the file."""
        return [self.__copy_ratio_median,
                self.__copy_ratio_10th_perc,
                self.__copy_ratio_90th_perc,
                self.__allele_fraction_median,
                self.__allele_fraction_10th_perc,
                self.__allele_fraction_90th_perc,
                self.__contig,
                self.__segment_start,
                self.__segment_end,
                self.__segment_lengths,
                self.__n_points_cr,
                self.__n_points_af,
                self.__n_segments,
                self.__max_copy_ratio_possible]

    def get_input_filename(self):
        """Return input filename."""
        return self.__filename

    def get_weights(self):
        """Return the weights of the points."""
        return self.__weights

    def get_sampled_points(self):
        """Return sampled posterior distributions as flattened lists."""
        return [self.__copy_ratio_sampled, self.__allele_fraction_sampled]

    def get_load_cr(self):
        """ Return whether we loaded the copy ratio values. """
        return self.__load_cr

    def get_load_af(self):
        """ Return whether we loaded the allele fraction values. """
        return self.__load_af

    def get_log_filename(self):
        """ Return log file name. """
        return self.__log_filename

    def get_logger(self):
        """Return the handler of the logger associated with the log file. """
        return self.__logger

    def __import_copy_ratio_and_allele_fraction(self):
        """ Read data related to the copy ratio and allele fraction distributions from file. """
        copy_ratio_median = []         # median of the copy ratio posterior distribution for the segment
        copy_ratio_10th_perc = []      # 10th percentile of the copy ratio posterior distribution
        copy_ratio_90th_perc = []      # 90th percentile of the copy ratio posterior distribution
        allele_fraction_median = []    # median of the allele fraction posterior distribution
        allele_fraction_10th_perc = [] # 10th percentile of the allele fraction posterior distribution
        allele_fraction_90th_perc = [] # 90th percentile of the allele fraction posterior distribution
        contig = []                    # contig (which chromosome the segments are on)
        segment_start = []             # starting positions of the segments (base pair)
        segment_end = []               # ending positions of the segments (base pair)
        segment_lengths = []           # length of segments (base pair)
        n_points_cr = []               # number of copy ratio data points taken from the segment (base pair)
        n_points_af = []               # number of allele fraction data points taken from the segment (base pair)

        # Read data from file
        file = open(self.__filename, "r")
        lines = file.readlines()

        # Check the fraction of lines that are NaN for copy ratio and allele fraction, respectively.
        # If we find too many NaNs (more than cr_NaN_threshold or af_NaN_threshold, respectively), then we do not load
        # the corresponding data.
        cr_nan = 0
        af_nan = 0
        all_sites = 0
        for line in lines:
            values = line.strip().split()
            if is_number(values[0]):
                if len(values) >= 5:
                    segment_length = float(values[2]) - float(values[1])
                    if (math.isnan(float(values[5]))
                        or math.isnan(float(values[6]))
                        or math.isnan(float(values[7]))
                        or 2**float(values[5]) > self.__max_copy_ratio_possible
                        or 2**float(values[6]) > self.__max_copy_ratio_possible
                        or 2**float(values[7]) > self.__max_copy_ratio_possible):
                        cr_nan += segment_length
                    if (math.isnan(float(values[8]))
                        or math.isnan(float(values[9]))
                        or math.isnan(float(values[10]))
                        or not (0. <= float(values[8]) <= 0.5)
                        or not (0. <= float(values[9]) <= 0.5)
                        or not (0. <= float(values[10]) <= 0.5)):
                        af_nan += segment_length
                    all_sites += segment_length
        cr_nan_ratio_threshold = 0.80
        af_nan_ratio_threshold = 0.80
        if all_sites == 0:
            cr_nan_ratio = 1
            af_nan_ratio = 1
        else:
            cr_nan_ratio = cr_nan / all_sites
            af_nan_ratio = af_nan / all_sites
        if cr_nan_ratio >= cr_nan_ratio_threshold:
            if self.__load_cr and not self.__log_filename == "":
                self.__logger.info("More than %s%% " % (100 * cr_nan_ratio_threshold) +
                                   "of the lines have NaN copy ratio values. " +
                                   "We will thus not load copy ratio data.")
            self.__load_cr = False
        if af_nan_ratio >= af_nan_ratio_threshold:
            if self.__load_cr and not self.__log_filename == "":
                self.__logger.info("More than %s%% " % (100 * af_nan_ratio_threshold) +
                                   "of the lines have NaN allele fraction values. " +
                                   "We will thus not load copy ratio data.")
            self.__load_af = False
        if (not self.__load_cr) and (not self.__load_af) and (not self.__log_filename == ""):
            self.__logger.info("Error: No copy ratio and no allele fraction data will be loaded.")

        # Load data
        if self.__load_cr and self.__load_af:
            for line in lines:
                values = line.strip().split()
                if is_number(values[0]):
                    if(len(values) >= 5):
                        if (not math.isnan(float(values[5]))
                            and not math.isnan(float(values[6]))
                            and not math.isnan(float(values[7]))
                            and not math.isnan(float(values[8]))
                            and not math.isnan(float(values[9]))
                            and not math.isnan(float(values[10]))
                            ):
                            if (2**float(values[5]) <= self.__max_copy_ratio_possible
                                and 2**float(values[6]) <= self.__max_copy_ratio_possible
                                and 2**float(values[7]) <= self.__max_copy_ratio_possible
                                and 0. <= float(values[8]) <= 0.5
                                and 0. <= float(values[9]) <= 0.5
                                and 0. <= float(values[10]) <= 0.5
                                ):
                                contig.append(int(values[0]))
                                segment_start.append(int(values[1]))
                                segment_end.append(int(values[2]))
                                segment_lengths.append(int(values[2]) - int(values[1]))
                                n_points_cr.append(int(values[3]))
                                n_points_af.append(int(values[4]))
                                copy_ratio_10th_perc.append(2**float(values[5]))
                                copy_ratio_median.append(2**float(values[6]))
                                copy_ratio_90th_perc.append(2**float(values[7]))
                                allele_fraction_10th_perc.append(float(values[8]))
                                allele_fraction_median.append(float(values[9]))
                                allele_fraction_90th_perc.append(float(values[10].rstrip("\\")))
        elif self.__load_cr:
            for line in lines:
                values = line.strip().split()
                if is_number(values[0]):
                    if len(values) >= 5:
                        if (not math.isnan(float(values[5]))
                            and not math.isnan(float(values[6]))
                            and not math.isnan(float(values[7]))
                            ):
                            if (2**float(values[5]) <= self.__max_copy_ratio_possible
                                and 2**float(values[6]) <= self.__max_copy_ratio_possible
                                and 2**float(values[7]) <= self.__max_copy_ratio_possible
                                ):
                                contig.append(int(values[0]))
                                segment_start.append(int(values[1]))
                                segment_end.append(int(values[2]))
                                segment_lengths.append(int(values[2]) - int(values[1]))
                                n_points_cr.append(int(values[3]))
                                n_points_af.append(0)
                                copy_ratio_10th_perc.append(2**float(values[5]))
                                copy_ratio_median.append(2**float(values[6]))
                                copy_ratio_90th_perc.append(2**float(values[7]))
                                allele_fraction_10th_perc.append(float('nan'))
                                allele_fraction_median.append(float('nan'))
                                allele_fraction_90th_perc.append(float('nan'))
        elif self.__load_af:
            for line in lines:
                values = line.strip().split()
                if is_number(values[0]):
                    if len(values) >= 5:
                        if (not math.isnan(float(values[8]))
                            and not math.isnan(float(values[9]))
                            and not math.isnan(float(values[10]))
                            ):
                            if(0. <= float(values[8]) <= 0.5
                               and 0. <= float(values[9]) <= 0.5
                               and 0. <= float(values[10]) <= 0.5
                               ):
                                contig.append(int(values[0]))
                                segment_start.append(int(values[1]))
                                segment_end.append(int(values[2]))
                                segment_lengths.append(int(values[2]) - int(values[1]))
                                n_points_cr.append(0)
                                n_points_af.append(int(values[4]))
                                copy_ratio_median.append(float('nan'))
                                copy_ratio_10th_perc.append(float('nan'))
                                copy_ratio_90th_perc.append(float('nan'))
                                allele_fraction_10th_perc.append(float(values[8]))
                                allele_fraction_median.append(float(values[9]))
                                allele_fraction_90th_perc.append(float(values[10].rstrip("\\")))
        n_segments = len(copy_ratio_median)
        return [copy_ratio_median, copy_ratio_10th_perc,
                copy_ratio_90th_perc, allele_fraction_median,
                allele_fraction_10th_perc, allele_fraction_90th_perc,
                contig, segment_start, segment_end, segment_lengths,
                n_points_cr, n_points_af, n_segments]

    def fit_beta_distribution(self, af_median: float, af_10: float, af_90: float):
        """ Fit the beta distribution to the allele fraction distribution,
            when the 50th (median), 10th and 90th percentiles are known.
        """

        def __find_beta_percentile(a_: float, b_: float, perc_: float):
            # Get the value corresponding to percentile perc from
            # the beta distribution with parameters a and b
            if perc_ < 0.:
                return 0.
            if perc_ > 1.:
                return 1.
            return scipy.special.betaincinv(a_, b_, perc_)

        def __find_initial_ab(af_median_: float, af_10_: float, af_90_: float):
            """ In this method, we set the initial conditions for the
                minimization algorithm that fits the parameters
                of the beta distribution that best describe the
                median, 10th and 90th percentiles.
                This initial condition should be close to the solution.

                Approximations:
                - We approximate the std deviation of the beta distribution with
                  0.5*(90th_percentile - 10th_percentile).
                  This is a rough but not very bad estimate.
                  The variance is given by var = a*b / ((a+b)**2 * (a+b+1)).
                - Since the beta distributions that we need to fit is usually tight
                  (and in the not tight case the minimization just works fine), we
                  approximate the mean by the median. The mean has a simple formula
                  mean = a / (a+b)
                We solve the approximate equations for the mean and the variance
                and get a and b from there.
            """
            std_dev_estimate = 0.5 * (af_90_ - af_10_)
            var_estimate = std_dev_estimate ** 2
            a0_ = (af_median_**2 * (1 - af_median_)) / var_estimate - af_median_
            b0_ = a0_ * (1 - af_median_) / af_median_
            a0_ = min(max(0.001, a0_), 10000)
            b0_ = min(max(0.001, b0_), 10000)
            return [a0_, b0_]

        def __beta_err(ab):
            """ This function defines the fitting error. We require that the median
                and the width of the distribution (difference between the 10th and
                90th percentile) be as accurate as possible. """
            # Note, that we do not take the errors of the 10th and 90th percentile,
            # only the error of determining their difference.
            # The reason for this choice is that the errors in the 10th and 90th
            # percentile can be much larger than those of the median, as these are
            # at the tails of the distribution.
            a_ = ab[0]
            b_ = ab[1]
            af_median_prediction = __find_beta_percentile(a_, b_, 0.5)
            af_10_prediction = __find_beta_percentile(a_, b_, 0.1)
            af_90_prediction = __find_beta_percentile(a_, b_, 0.9)

            objective = (abs((af_90 - af_10) - (af_90_prediction - af_10_prediction))**2
                         + abs(af_median - af_median_prediction)**2)
            return objective

        # We need bounds to make sure that the distribution is well-defined.
        # (It is defined only for a, b > 0.)
        ab_cutoff = [1e-6, None]
        bnds = ((ab_cutoff[0], ab_cutoff[1]), (ab_cutoff[0], ab_cutoff[1]))
        [a0, b0] = __find_initial_ab(af_median, af_10, af_90)
        res = minimize(__beta_err, (a0, b0), tol=1e-6, bounds=bnds, method="L-BFGS-B")
        a = res.x[0]
        b = res.x[1]
        return [a, b]

    def fit_gamma_distribution(self, cr_median: float, cr_10: float, cr_90: float):
        """ Fit the gamma distribution to the copy ratio distribution data, given
            the median, the 10th and the 90th percentile of the distribution.
        """

        def __find_gamma_percentile(a_, b_, perc_):
            # Get the copy ratio value corresponding to the percentile 'perc'
            if perc_ < 0:
                return 0
            if perc_ > 1:
                return 1
            return scipy.special.gammaincinv(a_, perc_) / b_

        def __find_initial_ab(cr_median_: float, cr_10_: float, cr_90_: float):
            """ In this method, we find the initial conditions for the minimization
                algorithm, i.e. we approximate the values of a and b and start
                minimizing from there.
                Approximations:
                - We approximate the std deviation with
                  0.5*(90th_percentile - 10th_percentile).
                  This is a rough but OK estimate.
                  The variance is given by var = a/b**2.
                - Since the beta distributions that we need to fit is usually tight,
                  we approximate the mean with the median.
                  The mean has a simple formula mean = a / b
                We solve the approximate equations for the mean and the
                variance and get a and b from there.
            """
            std_dev_estimate = 0.5 * (cr_90_ - cr_10_)
            var_estimate = std_dev_estimate ** 2
            mean_estimate = 0.3 * cr_median_ + 0.7 * 0.5 * (cr_10_ + cr_90_)

            # The choice of mixing the known parameters is slightly arbitrary,
            # but it approximates the mean relatively well
            b0_ = mean_estimate / var_estimate
            a0_ = mean_estimate * b0_
            return [a0_, b0_]

        def __gamma_err(ab):
            """ This function defines the fitting error. We require that
                the median and the width of the distribution
                (difference between the 10th and 90th percentile)
                be as accurate as possible.
                However, we do not want the 10th and 90th percentiles to be
                accurate, only their difference.
                The reason for this choice is that the errors in the 10th and 90th
                percentile can be large, as these are at the tails of the distribution.
            """
            a_ = ab[0]
            b_ = ab[1]
            cr_10_prediction = __find_gamma_percentile(a_, b_, 0.1)
            cr_90_prediction = __find_gamma_percentile(a_, b_, 0.9)
            cr_median_prediction = __find_gamma_percentile(a_, b_, 0.5)
            objective = (abs((cr_90 - cr_10) - (cr_90_prediction - cr_10_prediction))**2
                         + abs(cr_median - cr_median_prediction)**2)
            return objective

        # We need bounds to make sure that the distribution is well-defined.
        # (It is defined only for a, b > 0.)
        ab_cutoff = [1e-10, None]
        bnds = ((ab_cutoff[0], ab_cutoff[1]), (ab_cutoff[0], ab_cutoff[1]))
        [a0, b0] = __find_initial_ab(cr_median, cr_10, cr_90)
        res = minimize(__gamma_err, (a0, b0), tol=1e-10, bounds=bnds)
        a = res.x[0]
        b = res.x[1]
        return [a, b]

    def __determine_weights(self):
        """ This function determines the weights associated with the segment medians.
            It fits a gamma distribution to the copy ratio and a beta distribution to
            the allele fraction data. Then, it determines the variance of these
            distributions and sets the weight to be equal to the harmonic mean of the
            variances. (It assumes infinite covariance for missing data)
        """
        weights = []
        for i in range(self.__n_segments):
            w = self.__determine_weight_one_segment(self.__copy_ratio_median[i],
                                                    self.__copy_ratio_10th_perc[i],
                                                    self.__copy_ratio_90th_perc[i],
                                                    self.__allele_fraction_median[i],
                                                    self.__allele_fraction_10th_perc[i],
                                                    self.__allele_fraction_90th_perc[i])
            weights.append(w)
        cr_avg_weight = median([weights[i][0] for i in range(len(weights))])
        af_avg_weight = median([weights[i][1] for i in range(len(weights))])
        cr_avg_weight = median([weights[i][0] for i in range(len(weights)) if weights[i][0] > 0.001 * cr_avg_weight])
        af_avg_weight = median([weights[i][1] for i in range(len(weights)) if weights[i][1] > 0.001 * af_avg_weight])
        max_cr_weight = self.__cr_weight_ratio_max * cr_avg_weight
        max_af_weight = self.__af_weight_ratio_max * af_avg_weight
        weights = [[min([weights[i][0], max_cr_weight]), min([weights[i][1], max_af_weight])] for i in range(len(weights))]
        return weights

    def __determine_weight_one_segment(self, cr_median: float, cr_10: float, cr_90: float,
                                       af_median: float, af_10: float, af_90: float):
        """ Sample n_points data points from each segment, according to the gamma
            and beta distributions whose 10th, 50th and 90th percentiles are given.
        """
        if self.__load_cr:
            [cr_a, cr_b] = self.fit_gamma_distribution(cr_median=cr_median, cr_10=cr_10, cr_90=cr_90)
            var_cr = cr_a / cr_b**2
        else:
            var_cr = float("inf")
        if self.__load_af:
            [af_a, af_b] = self.fit_beta_distribution(2 * af_median, 2 * af_10, 2 * af_90)
            var_af = (af_a * af_b) / (af_a + af_b)**2 / (af_a + af_b + 1)
        else:
            var_af = float("inf")
        weight = [1./var_cr, 1./var_af]
        return weight

    def __sample_points(self):
        total_number_of_points = 5000
        avg_number_of_points = total_number_of_points / self.__n_segments
        avg_weight = np.sum(np.ndarray.flatten(np.array(self.__weights))) / self.__n_segments

        copy_ratio_sampled = []
        allele_fraction_sampled = []
        for i in range(self.__n_segments):
            n_pts_cr = int(round(avg_number_of_points * self.__weights[i][0] / avg_weight))
            n_pts_af = int(round(avg_number_of_points * self.__weights[i][1] / avg_weight))
            copy_ratio_sampled += [self.__copy_ratio_median[i]] * n_pts_cr
            allele_fraction_sampled += [self.__allele_fraction_median[i]] * n_pts_af
        return copy_ratio_sampled, allele_fraction_sampled

class ModeledSegmentsCaller:
    """ This class takes an instance of the class LoadAndSampleCrAndAf as input, which contains the
        data from the input file. It then identifies the normal segments. It outputs the plots of
        the segments (and optionally other plots about the fitting procedure) and saves the results
        into files.
    """
    def __init__(self, cr_af_data: LoadAndSampleCrAndAf, interactive: bool=True,
                 output_image_dir: str= "",
                 output_calls_dir: str= "",
                 output_image_prefix:str= "",
                 output_calls_prefix: str= "",
                 output_image_suffix: str=".png",
                 output_calls_suffix: str=".called.seg",
                 interactive_output_calls_image_suffix:str="_classification.png",
                 interactive_output_summary_plot_suffix:str="_summary_plot.png",
                 interactive_output_allele_fraction_plot_suffix:str= "_allele_fraction_CN1_and_CN2_candidate_intervals.png",
                 interactive_output_copy_ratio_suffix:str="_copy_ratio_fit.png",
                 interactive_output_copy_ratio_clustering_suffix:str="_copy_ratio_clusters.png",
                 normal_minor_allele_fraction_threshold: float=0.475,
                 copy_ratio_peak_min_relative_height: float=0.05,
                 copy_ratio_kernel_density_bandwidth: float=None,
                 min_weight_first_cr_peak_cr_data_only: float=0.35,
                 min_fraction_of_points_in_normal_allele_fraction_region: float=0.15,
                 responsibility_threshold_normal: float=0.5
                 ):
        """ On initialization, the caller loads the copy ratio and allele fraction data from
            the LoadAndSampleCrAndAf class. It then identifies normal segments and saves the
            results, including the plots.
        """

        # Initialize random number generator
        np.random.seed(RANDOM_SEED)

        # Start logging, continuing the log file that was written when cr_af_data was created
        self.__log_filename=cr_af_data.get_log_filename()
        self.__logger=cr_af_data.get_logger()
        if not self.__log_filename == "":
            self.__logger.info("* ---- %s ---- *" % str(datetime.datetime.now()))
            self.__logger.info("Initializing class CNVCaller")

        # Copy input data
        self.__cr_af_data = cr_af_data
        self.__load_cr = self.__cr_af_data.get_load_cr()
        self.__load_af = self.__cr_af_data.get_load_af()
        self.__interactive = interactive
        self.__output_image_dir = output_image_dir
        if (output_image_dir != "") and (output_image_dir[-1] != "/"):
            self.__output_image_dir += "/"
        self.__output_calls_dir = output_calls_dir
        if (output_calls_dir != "") and (output_calls_dir[-1] != "/"):
            self.__output_calls_dir += "/"
        self.__output_calls_prefix = output_calls_prefix
        self.__output_image_prefix = output_image_prefix
        self.__output_calls_suffix = output_calls_suffix
        self.__output_image_suffix = output_image_suffix
        self.__interactive_output_calls_image_suffix = interactive_output_calls_image_suffix
        self.__interactive_output_summary_plot_suffix = interactive_output_summary_plot_suffix
        self.__interactive_output_allele_fraction_plot_suffix = interactive_output_allele_fraction_plot_suffix
        self.__interactive_output_copy_ratio_suffix = interactive_output_copy_ratio_suffix
        self.__interactive_output_copy_ratio_clustering_suffix = interactive_output_copy_ratio_clustering_suffix
        self.__normal_minor_allele_fraction_threshold = normal_minor_allele_fraction_threshold
        self.__copy_ratio_peak_min_relative_height = copy_ratio_peak_min_relative_height
        self.__copy_ratio_kernel_density_bandwidth = copy_ratio_kernel_density_bandwidth
        self.__min_fraction_of_points_in_normal_allele_fraction_region = min_fraction_of_points_in_normal_allele_fraction_region
        self.__responsibility_threshold_normal = responsibility_threshold_normal
        self.__min_weight_first_cr_peak_cr_data_only = min_weight_first_cr_peak_cr_data_only

        # Set the maximal value of the PHRED score we allow (since we don't want it to be off the scale on the plots)
        self.__max_PHRED_score = 100.

        # Load data from file
        [self.__copy_ratio_median,
         self.__copy_ratio_10th_perc,
         self.__copy_ratio_90th_perc,
         self.__allele_fraction_median,
         self.__allele_fraction_10th_perc,
         self.__allele_fraction_90th_perc,
         self.__contig,
         self.__segment_start,
         self.__segment_end,
         self.__segment_lengths,
         self.__n_points_cr,
         self.__n_points_af,
         self.__n_segments,
         self.__max_copy_ratio_possible] = self.__cr_af_data.get_data()

        # Get the weights associated with each segment
        self.__weights = self.__cr_af_data.get_weights()

        # Auxiliary arrays: the number of times the values of the segment means appear in these arrays
        # is proportional to the weights associated with those segments
        [self.__copy_ratio_sampled, self.__allele_fraction_sampled] = self.__cr_af_data.get_sampled_points()

        # Set the filenames for the output data
        [self.__output_calls_filename, self.__fig_normal_segments_filename, self.__fig_del_ampl_filename,
         self.__fig_scatter_plot, self.__fig_allele_fraction_cn1_cn2_intervals, self.__fig_copy_ratio_fit,
         self.__fig_copy_ratio_clusters] = self.__set_output_filenames()
        if not self.__log_filename == "":
            self.__logger.info("Setting output filenames:")
            self.__logger.info("   Normal segments image file : %s" % self.__fig_normal_segments_filename)
            self.__logger.info("   Calls file : %s" % self.__output_calls_filename)
            if self.__interactive:
                self.__logger.info("   (interactive mode) Image of deletions and amplifications: %s"
                                   % self.__fig_del_ampl_filename)
                self.__logger.info("   (interactive mode) Scatter plot of segments: %s" % self.__fig_scatter_plot)
                self.__logger.info("   (interactive mode) Plots of allele fraction data with copy ratios in the " +
                                   "copy number 1 and 2 candidate intervals: %s"
                                   % self.__fig_allele_fraction_cn1_cn2_intervals)
                self.__logger.info("   (interactive mode) Gaussian fit to the copy ratio data: %s" % self.__fig_copy_ratio_fit)

        # Find normal segments
        if not self.__log_filename == "":
            self.__logger.info("Determining normal segments.")

        if self.__load_cr and self.__load_af:
            [self.__responsibilities_normal,
             self.__normal_segment_indices,
             self.__gaussian_mixture_fit] = self.__choose_normal_segments__cr_af_data()
        elif self.__load_cr:
            [self.__responsibilities_normal,
             self.__normal_segment_indices] = self.__choose_normal_segments__cr_data_only()
        elif self.__load_af:
            [self.__responsibilities_normal,
             self.__normal_segment_indices] = self.__choose_normal_segments__af_data_only()
        else:
            self.__responsibilities_normal = None
            self.__normal_segment_indices = None

        # Save plots of the segments
        if not self.__log_filename == "":
            self.__logger.info("Plotting and saving segments.")
        self.__plot_and_save_segments()

        # Create auxiliary plots if in interactive mode
        if self.__load_cr and self.__load_af and self.__interactive:
            if not self.__log_filename == "":
                self.__logger.info("Creating auxiliary plots in interactive mode.")
            self.__plot_clustering()

        # Save the results to a file
        if not self.__log_filename == "":
            self.__logger.info("Saving results.")
        self.__save_calls_to_file()

        # Finish logging
        if not self.__log_filename == "":
            self.__logger.info("Finished.\n\n\n")

    def __set_output_filenames(self):
        """ Set file names for all outputs. """
        # Get input directory
        input_filename_ending = self.__cr_af_data.get_input_filename().split("/")[-1]
        input_dir = self.__cr_af_data.get_input_filename().split(input_filename_ending)[0]
        if input_dir[-1] != "/":
            input_dir = input_dir + "/"

        if self.__output_image_dir == "":
            self.__output_image_dir = input_dir

        if self.__output_calls_dir == "":
            self.__output_calls_dir = input_dir

        # Set filenames for the regular output
        if self.__output_calls_prefix == "":
            calls_filename_prefix = self.__cr_af_data.get_input_filename()
            calls_fname_extension = calls_filename_prefix.split(".")[-1]
            calls_filename_prefix = calls_filename_prefix.split("." + calls_fname_extension)[0]
            calls_filename_prefix = calls_filename_prefix.split("/")[-1]
        else:
            calls_filename_prefix = self.__output_calls_prefix
        output_calls_filename = self.__output_calls_dir + calls_filename_prefix + self.__output_calls_suffix

        # Set file names for images produced in interactive mode
        if self.__output_image_prefix == "":
            image_fname_prefix = self.__cr_af_data.get_input_filename()
            image_fname_extension = image_fname_prefix.split(".")[-1]
            image_fname_prefix = image_fname_prefix.split("." + image_fname_extension)[0]
            image_fname_prefix = image_fname_prefix.split("/")[-1]
        else:
            image_fname_prefix = self.__output_image_prefix

        # Set the names of the outputs used in interactive mode
        fig_normal_segments_filename = self.__output_image_dir + image_fname_prefix + self.__output_image_suffix
        fig_del_ampl_filename = self.__output_image_dir +  image_fname_prefix + self.__interactive_output_calls_image_suffix
        fig_scatter_plot = self.__output_image_dir +  image_fname_prefix + self.__interactive_output_summary_plot_suffix
        fig_allele_fraction_cn1_cn2_intervals = self.__output_image_dir +  image_fname_prefix + self.__interactive_output_allele_fraction_plot_suffix
        fig_copy_ratio_fit = self.__output_image_dir +  image_fname_prefix + self.__interactive_output_copy_ratio_suffix
        fig_copy_ratio_clusters = self.__output_image_dir +  image_fname_prefix + self.__interactive_output_copy_ratio_clustering_suffix

        return [output_calls_filename, fig_normal_segments_filename, fig_del_ampl_filename,
                fig_scatter_plot, fig_allele_fraction_cn1_cn2_intervals, fig_copy_ratio_fit,
                fig_copy_ratio_clusters]

    def get_all_output_filenames(self):
        filenames = [self.__output_calls_filename, self.__fig_normal_segments_filename, self.__log_filename]
        if self.__interactive:
            filenames.append(self.__fig_normal_segments_filename)
            filenames.append(self.__fig_del_ampl_filename)
            filenames.append(self.__fig_scatter_plot)
            filenames.append(self.__fig_allele_fraction_cn1_cn2_intervals)
            filenames.append(self.__fig_copy_ratio_fit)
        return filenames

    def get_output_calls_filename(self):
        """Name of the output calls file."""
        return self.__output_calls_filename

    def get_all_output_fig_filenames(self):
        """Name of all filenames output by the tool."""
        if self.__interactive:
            return [self.__fig_normal_segments_filename,
                    self.__fig_del_ampl_filename,
                    self.__fig_scatter_plot,
                    self.__fig_allele_fraction_cn1_cn2_intervals,
                    self.__fig_copy_ratio_fit]
        else:
            return [self.__fig_normal_segments_filename]

    def get_output_fig_filename(self):
        """Name of the output fig file."""
        return self.__fig_normal_segments_filename

    def get_interactive_figure_filenames(self):
        """Get the names of auxiliary files that we save in interactive mode."""
        if not self.__interactive:
            return None
        return [self.__fig_del_ampl_filename,
                self.__fig_scatter_plot,
                self.__fig_allele_fraction_cn1_cn2_intervals,
                self.__fig_copy_ratio_fit]

    def get_normal_segment_indices(self):
        """Output the indices of the normal segments."""
        return self.__normal_segment_indices

    def get_responsibilities_segments_are_normal(self):
        """Output the responsibility associated with each segment being normal."""
        return self.__responsibilities_normal

    def get_interactive(self):
        """Return whether the code was run in interactive mode."""
        return self.__interactive


    def __k_means(self, data, data_weights, max_n_clusters):
        # Determines the optimal number of clusters based on the silhouette coefficient
        silhouette_coeff = []
        for n_cl in range(2, max_n_clusters):
            kmeans_model = KMeans(n_clusters=n_cl, random_state=1).fit(data)
            labels = kmeans_model.labels_
            silhouette_coeff.append([n_cl, metrics.silhouette_score(data, labels, metric='euclidean')])
        i = np.argmax(np.array(silhouette_coeff)[:,1])
        n_clusters = silhouette_coeff[i][0]

        kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(data)
        plt.scatter([data[i][0] for i in range(len(data))],
                    [data[i][1] for i in range(len(data))],
                    c='g', alpha=0.5)
        plt.scatter(np.array(kmeans.cluster_centers_)[:, 0], np.array(kmeans.cluster_centers_)[:, 1], c="r")
        plt.title("Finding cluster centers using k-means")
        plt.show()
        return [kmeans.cluster_centers_, n_clusters]

    def __fit_weighted_Gaussian_mixture_ADVI(self, data_points, data_weights, max_n_Gaussians):
        # Function that fits Gaussians to data with weights. The dimension of the data can be arbitrary.

        if np.array(data_weights).any() == None:
            data_weights = [1, 1] * len(data_points)
        data = np.array([[data_points[i][0], data_points[i][1],
                          data_weights[i][0], data_weights[i][1]]
                         for i in range(len(data_points))])
        data_dim = 2

        # Log likelihood of normal distribution
        def logp_normal(mu, tau, value):
            # log probability of individual samples
            # The first n coordinates of the input 'value' contains the data point coordinates,
            # whereas the last coordinate is the weight of the point.
            k = tau.shape[0]
            delta = lambda mu: value[:, 0:k] - mu
            return (-1 / 2.) * tt.sum([k * tt.log(2 * np.pi) - tt.log(det(tau))
                                       - (tt.log(value[i][-2]) + tt.log(value[i][-1]))
                                       + mu[0]**2 * tau[0][0] * value[i][-2] + mu[1]**2 * tau[1][1] * value[i][-1]
                                       for i in range(len(value))])

        # Log likelihood of Gaussian mixture distribution
        def logp_gmix(mus, pis, taus, n_Gaussians):
            def logp_(value):
                logps = [tt.log(pis[i]) + logp_normal(mus[i], taus[i], value) for i in range(n_Gaussians)]
                return tt.sum(logsumexp(logps), axis=0)
            return logp_

        [mu_estimates, n_Gaussians] = self.__k_means(data=data_points,
                                                     data_weights=data_weights,
                                                     max_n_clusters=max_n_Gaussians)

        pi_init = np.random.uniform(0, 1, n_Gaussians)
        pi_init = pi_init / np.sum(pi_init)
        with pm.Model() as model:
            mus = [MvNormal('mu_%d' % i,
                            mu=pm.floatX(np.zeros(data_dim)),
                            tau=pm.floatX(0.1*np.eye(data_dim)),
                            shape=(data_dim,), testval=mu_estimates[i])
                   for i in range(n_Gaussians)]
            std_dev_data = np.std(data)
            taus =  [tt.nlinalg.alloc_diag([pm.Gamma('tau_%d_%d' % (i, j),
                                                     mu=std_dev_data/3,
                                                     sd=std_dev_data/8,
                                                     testval=std_dev_data/3)
                                            for j in range(data_dim)])
                     for i in range(n_Gaussians)]
            pis = Dirichlet('pis', a=pm.floatX(0.1*np.ones(n_Gaussians)), shape=(n_Gaussians,), testval=pi_init)
            xs = DensityDist('x', logp_gmix(mus, pis, taus, n_Gaussians), observed=data)

        with model:
            approx = pm.fit(n=40000)

        means = approx.bij.rmap(approx.mean.eval())
        cov = approx.cov.eval()
        sds = approx.bij.rmap(np.diag(cov)**.5)

        mu_result = []
        pi_result = means['pis_stickbreaking__']
        cov_result = []
        for i in range(n_Gaussians):
            mu_result.append(means['mu_' + str(i)])
        for i in range(n_Gaussians):
            cov_result.append([np.exp(means["tau_" + str(i) + "_" + str(j) + "_log__"]) for j in range(data_dim)])

        plt.figure(figsize=(5, 5))
        plt.scatter(data[:, 0], data[:, 1], c='g', alpha=0.5)
        colors = ["r", "b", "k", "c", "m"]
        for i in range(n_Gaussians):
            plt.scatter(mu_result[i][0], mu_result[i][1], c=colors[i % len(colors)], s=100)
        plt.xlim(1.5*min(data[:,0]), 1.5*max(data[:,0]))
        plt.ylim(1.5*min(data[:,1]), 1.5*max(data[:,1]))

        return [pi_result, mu_result, cov_result]

    def __choose_normal_segments__cr_af_data(self, n_Gaussians: int=10):
        """ Choose those Gaussians from the fitted distribution that cover the
            normal segments. The Gaussians are fitted using Bayesian variational
            inference in the  two dimensional space of copy ratio and allele fraction axes.
        """
        np.random.seed(RANDOM_SEED)

        # Approximate the standard deviation of the prior for fitting the peaks
        sigma_cr_normal = 0.005
        sigma_af_normal = 0.005

        # Fit mixtures using Gaussian variational inference
        data = [[self.__copy_ratio_median[i], self.__allele_fraction_median[i]] for i in range(self.__n_segments)]
        [pis, mus, covs] = self.__fit_weighted_Gaussian_mixture_ADVI(data, self.__weights, n_Gaussians)

        # We choose those peaks to be normal whose mean's copy ratio value is within the range specified
        # by '__choose_cn2_cr_cluster' and whose allele fraction value is within the range specified by
        # 'normal_range_af'.
        normal_range_af = [self.__normal_minor_allele_fraction_threshold, 0.5]
        normal_range_cr = self.__choose_cn2_cr_cluster()

        normal_peak_indices = []
        for i in range(len(pis)):
            if (normal_range_cr[0] - min(np.sqrt(covs[i][0]), 0.25) <= mus[i][0]
                and mus[i][0] <= normal_range_cr[1] + min(np.sqrt(covs[i][0]), 0.25)
                and (normal_range_af[0] <= mus[i][1] + np.sqrt(covs[i][1]))
                and (mus[i][1] <= normal_range_af[1])
                ):
                normal_peak_indices.append(i)

        responsibilities_normal = self.__responsibility_segments_are_normal(pis, mus, covs, normal_peak_indices)
        normal_segment_indices = []
        for i in range(len(responsibilities_normal)):
            if responsibilities_normal[i] >= self.__responsibility_threshold_normal:
                normal_segment_indices.append(i)
        pis = [abs(p) for p in pis]
        gmix = [pis, mus, covs]
        return [responsibilities_normal, normal_segment_indices, gmix]

    def __choose_normal_segments__af_data_only(self):
        """ This function tries to determine which segments are normal relying only on allele fraction data.
            It simply assumes that the allele fraction simply needs to be above a certain threshold so that
            it is considered normal.
        """
        responsibilities_normal = [0] * len(self.__allele_fraction_median)
        normal_segment_indices = []
        for i in range(len(self.__allele_fraction_median)):
            responsibilities_normal[i] = np.sum([self.__weights[i][1] if a0 > self.__allele_fraction[i] else 0])
            if responsibilities_normal[i] >= self.__responsibility_threshold_normal:
                normal_segment_indices.append(i)
        return [responsibilities_normal, normal_segment_indices]

    def __choose_normal_segments__cr_data_only(self):
        """ This function determines which peak in the 1D copy ratio data is the copy number 2 peak.
        """
        cluster_separators = self.__find_1D_clusters(data=self.__copy_ratio_sampled,
                                                     interactive_image_filename=self.__fig_copy_ratio_fit)

        if len(cluster_separators) == 0:
            cn2_interval = [0, 1.1 * max(self.__copy_ratio_median)]
        else:
            weights_per_cluster = [0] * (len(cluster_separators) + 1)
            ordering = self.__indices_increasing_order(self.__copy_ratio_median)
            j = 0
            total_weight = np.sum(np.array(self.__weights)[:, 0])
            for i in range(len(cluster_separators)):
                while self.__copy_ratio_median[ordering[j]] < cluster_separators[i]:
                    j += 1
                    weights_per_cluster[i] += self.__weights[ordering[j][0]]
                weights_per_cluster[i] /= total_weight
            weights_per_cluster[-1] = 1. - np.sum([weights_per_cluster[i] for i in range(len(weights_per_cluster)-1)])

            if (weights_per_cluster[0] / (weights_per_cluster[0] + weights_per_cluster[1])
                    >= self.__min_weight_first_cr_peak_cr_data_only):
                cn2_interval = [0, cluster_separators[0]]
            else:
                if len(cluster_separators) == 1:
                    cn2_interval = [cluster_separators[0], 1.1 * max(self.__copy_ratio_median)]
                else:
                    cn2_interval = [cluster_separators[0], cluster_separators[1]]

        responsibilities_normal = [0] * len(self.__copy_ratio_median)
        normal_segment_indices = []
        for i in range(len(self.__copy_ratio_median)):
            if cn2_interval[0] <= self.__copy_ratio_median[i] <= cn2_interval[1]:
                responsibilities_normal[i] = 1
                normal_segment_indices.append(i)
            else:
                responsibilities_normal[i] = 0

        return [responsibilities_normal, normal_segment_indices]

    def __responsibility_segments_are_normal(self, pis, means, covariances, normal_peak_indices):
        """ Based on a Gaussian fit to the clusters in copy ratio and allele fraction space,
            figure out the responsibilities of a segment being normal.
        """
        n_Gaussians = len(pis)

        # The Gaussian mixture model determined the Gaussian mixture model that fits the data the best
        # Now, we assign the probability of each segment being Gaussian to each segment.
        # This probability is defined as the 'responsibility' assigned to the segment median
        responsibilities_normal = [0] * self.__n_segments
        for i in range(self.__n_segments):
            normal_pdf = np.sum([pis[j] * scipy.stats.multivariate_normal.pdf(
                [self.__copy_ratio_median[i], self.__allele_fraction_median[i]],
                means[j], covariances[j]) for j in normal_peak_indices])
            total_pdf  = np.sum([pis[j] * scipy.stats.multivariate_normal.pdf(
                [self.__copy_ratio_median[i], self.__allele_fraction_median[i]],
                means[j], covariances[j])
                                 for j in range(n_Gaussians)])
            responsibilities_normal[i] = max(min(normal_pdf / total_pdf, 1), 0)
        return responsibilities_normal

    def __find_1D_clusters(self, data, interactive_image_filename: str):
        """Looks for clusters in 1D data by taking a Gaussian kernel density estimation of the data
           and finding the local minima of the density. The width of the Gaussian is chosen automatically,
           based on the data.
        """

        def __find_extremal_indices(myArray):
            # Find the indices where the data has local minima or maxima
            extremal_indices = []
            for i in range(2, len(myArray) - 1):
                if myArray[i] < myArray[i-1] and myArray[i] < myArray[i+1]:
                    # found local maximum
                    extremal_indices.append(i)
                elif myArray[i-1] < myArray[i] and myArray[i+1] < myArray[i]:
                    # found local maximum
                    extremal_indices.append(i)
            return extremal_indices

        def __find_minimum_indices(myArray):
            # Find the indices where the data has local minima
            minimum_indices = []
            for i in range(2, len(myArray) - 1):
                if myArray[i] < myArray[i-1] and myArray[i] < myArray[i+1]:
                    minimum_indices.append(i)
            return minimum_indices

        # Estimate the bandwidth of the Gaussians used for density estimation: we fit a set of Gaussians
        # to the data and take the standard deviation of the smallest peak as the bandwidth
        n_peaks, _, _, _ = self.__estimate_number_of_cr_clusters(data, max_n_Gaussians = 10,
                                                                 alpha = 0.1, min_std_dev = 0.05)
        n_peaks, w_peaks, mu_peaks, sd_peaks = self.__estimate_number_of_cr_clusters(data,
                                                                                     max_n_Gaussians = n_peaks + 2,
                                                                                     alpha = 0.1, min_std_dev = 0.05)
        ind = self.__indices_increasing_order(list(mu_peaks))
        if self.__copy_ratio_kernel_density_bandwidth == None:
            if len(ind) == 1:
                bandwidth = 0.5 * sd_peaks[0]
            else:
                bandwidth = 0.5 * min([sd_peaks[ind[0]], sd_peaks[ind[1]]])
        else:
            bandwidth = self.__copy_ratio_kernel_density_bandwidth

        bandwidth = 0.05
        print("n_peaks = ", n_peaks)
        print("w_peaks = ", w_peaks)
        print("mu_peaks = ", mu_peaks)
        print("sd_peaks = ", sd_peaks)

        # Convolve the data with a Gaussian kernel ("Gaussian kernel density estimation")
        dataTransposed = np.array(data)[:, np.newaxis]
        x = np.linspace(min(0, min(data)), max(5, max(data)), 1000)[:, np.newaxis]
        kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(dataTransposed)
        kDensity = np.exp(kde.score_samples(x))

        extremal_indices = __find_extremal_indices(kDensity)
        extremal_indices.append(0)
        extremal_indices.append(len(kDensity)-1)
        extremal_indices.sort()

        # Filter out indices corresponding to small fluctuations only
        scale = max(kDensity) - min(kDensity)
        threshold = scale * self.__copy_ratio_peak_min_relative_height
        filtered_extremal_indices = []
        if len(extremal_indices) > 0:
            last_value = kDensity[extremal_indices[0]]
            filtered_extremal_indices = [extremal_indices[0]]
            for i in range(1, len(extremal_indices)):
                if abs(last_value - kDensity[extremal_indices[i]]) >= threshold:
                    filtered_extremal_indices.append(extremal_indices[i])
                    last_value = kDensity[extremal_indices[i]]

        filtered_extremal_values = [kDensity[i] for i in filtered_extremal_indices]
        filered_minimum_indices = __find_minimum_indices(filtered_extremal_values)
        cluster_separators = [x[filtered_extremal_indices[i]] for i in filered_minimum_indices]

        if self.__interactive:
            xPlot = np.linspace(0, 5, 1000)
            fig = plt.figure(3, dpi=400)
            plt.subplot(311)
            cols = ["r", "m", "c", "k", "y"]
            plt.hist(self.__copy_ratio_median, weights=np.array(self.__weights)[:,0], density=1, bins=200)
            for i in range(len(w_peaks)):
                y = [0] * len(xPlot)
                rv = multivariate_normal(mu_peaks[i], sd_peaks[i]**2)
                y = w_peaks[i] * rv.pdf(xPlot)
                plt.plot(xPlot, y, cols[i % len(cols)])
            plt.xlim(0, 5)
            plt.ylabel('Estimate of cluster widths')

            plt.subplot(312)
            plt.plot(x, kDensity, '-')
            plt.plot(cluster_separators, np.zeros(len(cluster_separators)), 'r*')
            plt.xlim(0, 5)
            plt.ylabel('Kernel density')

            plt.subplot(313)
            plt.hist(data, density=1, bins=200)
            plt.plot(cluster_separators, np.zeros(len(cluster_separators)), 'r*')
            plt.xlim(0, 5)
            plt.xlabel('Copy ratio')
            plt.ylabel('Histogram and clusters')
            plt.savefig(interactive_image_filename)
            plt.close(fig)

        return cluster_separators

    def __choose_cn2_cr_cluster(self):
        """ This function determines which peak in the 1D copy ratio data is the copy number 2 peak.
            It also makes use of the allele fraction data to make this decision.
        """
        normal_cr_segment_min = max([0, 0.9 * min(self.__copy_ratio_median)])
        normal_cr_segment_max = min([1.1 * max(self.__copy_ratio_median), 100])

        cluster_separators = self.__find_1D_clusters(self.__copy_ratio_sampled, self.__fig_copy_ratio_fit)

        if len(cluster_separators) == 0:
            cn2_interval = [0, 1.1 * max(self.__copy_ratio_median)]
            return cn2_interval

        cn1_interval_candidate = [0, cluster_separators[0]]
        if len(cluster_separators) == 1:
            cn2_interval_candidate = [cluster_separators[0], normal_cr_segment_max]
        else:
            cn2_interval_candidate = [cluster_separators[0], cluster_separators[1]]

        if self.__interactive:
            fig = plt.figure(1, dpi=300)
            plt.hist(self.__copy_ratio_median, weights=np.array(self.__weights)[:,0], color = "k", bins = 100, density = True)
            plt.plot(cn1_interval_candidate[0], [0], "rv")
            plt.plot(cn1_interval_candidate[1], [0], "rv")
            plt.plot(cn2_interval_candidate[0], [0], "g*")
            plt.plot(cn2_interval_candidate[1], [0], "g*")
            plt.xlim(0, 5)
            plt.xlabel("Copy ratio")
            plt.title("Histogram and segmentation")
            plt.savefig(self.__fig_copy_ratio_clusters)
            plt.close(fig)

        # If at least 'self.__min_fraction_of_points_in_normal_allele_fraction_region' fraction of the allele fraction
        # points in the 'CN1_interval_candidate' falls within the range 'normal_range_af', then we call that the copy
        # number 2 segment. Otherwise, it will be 'CN2_interval_candidate'.
        normal_range_af = [self.__normal_minor_allele_fraction_threshold, 0.5]
        af__cn1_interval_candidate = []
        af__cn2_interval_candidate = []
        for i in range(len(self.__allele_fraction_median)):
            if cn2_interval_candidate[0] <= self.__copy_ratio_median[i] < cn2_interval_candidate[1]:
                af__cn2_interval_candidate.append([self.__allele_fraction_median[i], self.__weights[i][0]])
            elif cn1_interval_candidate[0] <= self.__copy_ratio_median[i] < cn1_interval_candidate[1]:
                af__cn1_interval_candidate.append([self.__allele_fraction_median[i], self.__weights[i][1]])
                []
        normal_range_af_cn1_interval_candidate = [af for af in af__cn1_interval_candidate
                                                  if (normal_range_af[0] <= af[0] <= normal_range_af[1])]

        if (np.sum(np.array(normal_range_af_cn1_interval_candidate)[:,1]) / np.sum(np.array(af__cn1_interval_candidate))
                > self.__min_fraction_of_points_in_normal_allele_fraction_region):
            cn2_interval = cn1_interval_candidate
            color_1 = 'k'
            color_2 = 'r'
            if self.__interactive and self.__log_filename != "":
                self.__logger.info("We chose the 1st interval to be of copy number 2.")
        else:
            cn2_interval = cn2_interval_candidate
            color_1 = 'r'
            color_2 = 'k'
            which_interval = 2
            if self.__interactive and self.__log_filename!="":
                self.__logger.info("We chose the 2nd interval to be of copy number 2.")

        # Plot histograms of allele fraction data falling in the intervals 'cn1_interval_candidate'
        # and 'cn2_interval_candidate'
        if self.__interactive:
            fig = plt.figure(2, dpi=400)
            plt.subplot(211)
            plt.hist(np.array(af__cn1_interval_candidate)[:,0],
                     weights=np.array(af__cn1_interval_candidate)[:,1],
                     bins=50, color=color_1)
            plt.xlim(0, 0.5)
            plt.xlabel("Allele fraction of the first interval")
            plt.subplot(212)
            plt.hist(np.array(af__cn2_interval_candidate)[:,0],
                     weights=np.array(af__cn2_interval_candidate)[:,1],
                     bins=50, color=color_2)
            plt.xlim(0, 0.5)
            plt.xlabel("Allele fraction of the second interval")
            plt.subplots_adjust(hspace=0.5)
            plt.savefig(self.__fig_allele_fraction_cn1_cn2_intervals)
            plt.close(fig)

        cn2_interval[0] = max(cn2_interval[0], normal_cr_segment_min)
        cn2_interval[1] = min(cn2_interval[1], normal_cr_segment_max)
        return cn2_interval

    def __indices_increasing_order(self, ls: List):
        """ Get the indices of the ordered version of ls.
        """
        ls_new = deepcopy(list(ls))
        ls_new.sort()
        ordering = [ls.index(i) for i in ls_new]
        return ordering

    def __fit_normal_1d(self, data: List, n_Gaussians: int, n_iterations: int):
        """ Fit Gaussians to the one dimensional data.
            n_Gaussians = the number of Gaussians used
            n_iterations: the number of scikit learn iterations used.
        """
        random.seed(RANDOM_SEED)
        gmix = mixture.BayesianGaussianMixture(n_components=n_Gaussians,
                                               weight_concentration_prior_type="dirichlet_distribution",
                                               weight_concentration_prior = 0.01 / n_Gaussians,
                                               covariance_prior = np.array([[0.01]]),
                                               max_iter = n_iterations)
        gmix.fit(np.asarray(data).reshape(-1,1))
        return [gmix.weights_, gmix.means_, gmix.covariances_]

    def __estimate_number_of_cr_clusters(self, data, max_n_Gaussians: int=10, alpha: float=1.,
                                         min_std_dev: float=0.03, max_std_dev: float=0.4):
        """ This function throws down max_n_Gaussians peaks on the one dimensional copy ratio data
            self.copy_ratio_sampled. Then, it merges those that have substantial overlap.
            The overlap is measured based on whether the peaks overlap up to alpha times their standard deviation.
            In order to avoid the need of dealing with too narrow or too wide peaks, we limit their minimum
            and maximum width using 'min_std_dev' and 'max_std_dev'.
        """

        def __overlapping_peaks(mu1: float, cov1: float, mu2: float, cov2: float, alpha0: float=alpha,
                                min_std_dev0: float=min_std_dev, max_std_dev0: float=max_std_dev):
            """ Function to determine whether we consider two Gaussian peaks to be overlapping.
                alpha is the scaling factor of their standard deviations.
            """
            std_dev1 = min(max_std_dev0, max(min_std_dev0, np.sqrt(cov1)))
            std_dev2 = min(max_std_dev0, max(min_std_dev0, np.sqrt(cov2)))
            if ((mu1 <= mu2 and (mu1 + alpha0 * std_dev1) >= (mu2 - alpha0 * std_dev2))
                or (mu2 <= mu1 and (mu2 + alpha0 * std_dev2) >= (mu1 - alpha0 * std_dev1))
                ):
                return True
            return False

        [w, mu, cov] = self.__fit_normal_1d(data=data,
                                            n_Gaussians=max_n_Gaussians,
                                            n_iterations=10000)
        weight_threshold = max(0.25 / max_n_Gaussians, 0.08)
        peak_indices = []
        for i in range(max_n_Gaussians):
            if w[i] > weight_threshold:
                peak_indices.append(i)
        peak_groups = [[peak_indices[-1]]]
        peak_indices.remove(peak_indices[-1])

        for i in peak_indices:
            # We first figure out which peak groups overlap with the peak i
            # The index of these groups in 'peak_groups' will be stored in
            # the variable 'groups_overlapping_with_i'
            groups_overlapping_with_i = []
            for k in range(len(peak_groups)):
                for j in peak_groups[k]:
                    if __overlapping_peaks(mu[i], cov[i][0][0], mu[j], cov[j][0][0]):
                        groups_overlapping_with_i.append(k)
                        break
            if len(groups_overlapping_with_i) == 0:
                peak_groups.append([i])
            elif len(groups_overlapping_with_i) == 1:
                peak_groups[groups_overlapping_with_i[0]].append(i)
            else:
                new_group_components = []
                new_group = []
                for g_ind in groups_overlapping_with_i:
                    new_group_components.append(peak_groups[g_ind])
                    new_group = new_group + peak_groups[g_ind]
                for g in new_group_components:
                    peak_groups.remove(g)
                peak_groups.append(new_group)

        n_peaks = len(peak_groups)
        mu_peaks = [0] * n_peaks
        sd_peaks = [0] * n_peaks
        w_peaks = [0] * n_peaks
        for i in range(n_peaks):
            for j in peak_groups[i]:
                w_peaks[i] = w_peaks[i] + w[j]
                mu_peaks[i] = mu_peaks[i] + w[j] * mu[j]
                sd_peaks[i] = sd_peaks[i] + w[j] * cov[j]
            mu_peaks[i] = mu_peaks[i] / w_peaks[i]
            sd_peaks[i] = np.sqrt(sd_peaks[i] / w_peaks[i])
        w_total = np.sum(w_peaks)
        w_peaks = [w_p / w_total for w_p in w_peaks]

        mu_peaks = np.array(mu_peaks).flatten()
        sd_peaks = np.array(sd_peaks).flatten()
        w_peaks  = np.array(w_peaks).flatten()

        return n_peaks, w_peaks, mu_peaks, sd_peaks

    def __plot_and_save_segments(self):
        """ This function plots the results of the calling scripts and also saves
            these results.
        """
        [avg_normal_cr, std_dev_normal_cr] = self.__average_and_std_dev_copy_ratio_normal_segments()

        x = []
        y_cr_median = []
        y_cr_10 = []
        y_cr_90 = []
        y_min_af_median = []
        y_min_af_10 = []
        y_min_af_90 = []
        y_maj_af_median = []
        y_maj_af_10 = []
        y_maj_af_90 = []
        y_responsibilities_normal = []
        y_responsibilities_normal_PHRED = []
        y_color = []

        site = 0
        for i in range(self.__n_segments):
            x.append([site, site + self.__segment_lengths[i]])
            site += self.__segment_lengths[i]
            y_cr_median.append([self.__copy_ratio_median[i]] * 2)
            y_cr_10.append([self.__copy_ratio_10th_perc[i]] * 2)
            y_cr_90.append([self.__copy_ratio_90th_perc[i]] * 2)
            y_min_af_median.append([self.__allele_fraction_median[i]] * 2)
            y_min_af_10.append([self.__allele_fraction_10th_perc[i]] * 2)
            y_min_af_90.append([self.__allele_fraction_90th_perc[i]] * 2)
            y_maj_af_median.append([1-self.__allele_fraction_median[i]] * 2)
            y_maj_af_10.append([1-self.__allele_fraction_10th_perc[i]] * 2)
            y_maj_af_90.append([1-self.__allele_fraction_90th_perc[i]] * 2)
            y_responsibilities_normal.append([self.__responsibilities_normal[i]] * 2)
            y_responsibilities_normal_PHRED.append([self.__get_phred_score(
                probability=self.__responsibilities_normal[i],
                max_phred_score=self.__max_PHRED_score
            )] * 2)

            n_d_a = self.__normal_del_ampl(self.__copy_ratio_median[i], avg_normal_cr,
                                           std_dev_normal_cr, self.__responsibilities_normal[i])
            if n_d_a == "+":
                y_color.append((1.0, 0.0, 1.0))
            elif n_d_a == "0":
                y_color.append((0.0, 1.0, 0.0))
            elif n_d_a == "CNLOH":
                y_color.append((0.0, 0.0, 1.0))
            else:
                y_color.append((0.4, 0.4, 0.4))

        contig_beginning_end = []
        current_contig = self.__contig[0]
        current_contig_beginning = 1
        site = 0
        for i in range(self.__n_segments):
            if self.__contig[i] != current_contig:
                contig_beginning_end.append([current_contig_beginning, site])
                current_contig_beginning = site + 1
                current_contig = self.__contig[i]
            site += self.__segment_lengths[i]
        contig_beginning_end.append([current_contig_beginning, site])

        fig1, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, dpi=400, figsize=(8, 8))
        self.__gray_background_contigs(contig_beginning_end, 0, 1.05 * self.__max_PHRED_score, ax1)

        lines_PHRED = []
        for i in range(len(x)):
            lines_PHRED.append([(x[i][0], y_responsibilities_normal_PHRED[i][0]),
                                (x[i][1], y_responsibilities_normal_PHRED[i][1])])
        lc_PHRED = LineCollection(lines_PHRED, color=(0.35, 0.35, 0.35))
        ax1.add_collection(lc_PHRED)
        ax1.set_xlim(contig_beginning_end[0][0], contig_beginning_end[-1][1])
        ax1.set_ylim([0, 1.05 * self.__max_PHRED_score])
        ax1.set_ylabel('PHRED(P(normal))')

        self.__gray_background_contigs(contig_beginning_end, -0.02, 3, ax2)
        lines_cr_10 = []
        lines_cr_median = []
        lines_cr_90 = []
        lines_cr_color = []
        for i in range(len(x)):
            lines_cr_10.append([(x[i][0], y_cr_10[i][0]), (x[i][1], y_cr_10[i][1])])
            lines_cr_median.append([(x[i][0], y_cr_median[i][0]), (x[i][1], y_cr_median[i][1])])
            lines_cr_90.append([(x[i][0], y_cr_90[i][0]), (x[i][1], y_cr_90[i][1])])
            lines_cr_color.append((1-self.__responsibilities_normal[i], 0, 0))
        lc_cr_10 = LineCollection(lines_cr_10, colors=lines_cr_color, linewidth=0.8)
        lc_cr_median = LineCollection(lines_cr_median, colors=lines_cr_color, linewidth=1.2)
        lc_cr_90 = LineCollection(lines_cr_90, colors=lines_cr_color, linewidth=0.8)
        ax2.add_collection(lc_cr_10)
        ax2.add_collection(lc_cr_median)
        ax2.add_collection(lc_cr_90)
        ax2.set_xlim(contig_beginning_end[0][0], contig_beginning_end[-1][1])
        ax2.set_ylim([-0.02, 3])
        ax2.set_ylabel('Copy ratio')
        ax2.plot([], color='r', label='Not normal')
        ax2.plot([], color='k', label='Normal')
        ax2.legend(loc="upper right", bbox_to_anchor=(1,1))

        self.__gray_background_contigs(contig_beginning_end, -0.02, 1.02, ax3)
        lines_min_af_10 = []
        lines_min_af_median = []
        lines_min_af_90 = []
        lines_maj_af_10 = []
        lines_maj_af_median = []
        lines_maj_af_90 = []
        lines_af_color = []
        for i in range(len(x)):
            lines_min_af_10.append([(x[i][0], y_min_af_10[i][0]), (x[i][1], y_min_af_10[i][1])])
            lines_min_af_median.append([(x[i][0], y_min_af_median[i][0]), (x[i][1], y_min_af_median[i][1])])
            lines_min_af_90.append([(x[i][0], y_min_af_90[i][0]), (x[i][1], y_min_af_90[i][1])])
            lines_af_color.append((1-self.__responsibilities_normal[i], 0, 0))
            lines_maj_af_10.append([(x[i][0], y_maj_af_10[i][0]), (x[i][1], y_maj_af_10[i][1])])
            lines_maj_af_median.append([(x[i][0], y_maj_af_median[i][0]), (x[i][1], y_maj_af_median[i][1])])
            lines_maj_af_90.append([(x[i][0], y_maj_af_90[i][0]), (x[i][1], y_maj_af_90[i][1])])
        lc_min_af_10 =LineCollection(lines_min_af_10, colors=lines_af_color, linewidth=0.8)
        lc_min_af_median = LineCollection(lines_min_af_median, colors=lines_af_color, linewidth=1.2)
        lc_min_af_90 = LineCollection(lines_min_af_90, colors=lines_af_color, linewidth=0.8)
        lc_maj_af_10 = LineCollection(lines_maj_af_10, colors=lines_af_color, linewidth=0.8)
        lc_maj_af_median = LineCollection(lines_maj_af_median, colors=lines_af_color, linewidth=1.2)
        lc_maj_af_90 = LineCollection(lines_maj_af_90, colors=lines_af_color, linewidth=0.8)
        ax3.add_collection(lc_min_af_10)
        ax3.add_collection(lc_min_af_median)
        ax3.add_collection(lc_min_af_90)
        ax3.add_collection(lc_maj_af_10)
        ax3.add_collection(lc_maj_af_median)
        ax3.add_collection(lc_maj_af_90)
        ax3.set_xlim(contig_beginning_end[0][0], contig_beginning_end[-1][1])
        ax3.set_ylim([-0.02, 1.02])
        ax3.set_ylabel('Allele fraction')
        ax3.plot([], color='r', label='Not normal')
        ax3.plot([], color='k', label='Normal')
        ax3.legend(loc="upper right", bbox_to_anchor=(1,1))

        plt.savefig(self.__fig_normal_segments_filename)
        plt.close(fig1)

        if self.__interactive:
            fig2, (ax1, ax2) = plt.subplots(2, sharex=True, dpi=400, figsize=(8, 6))
            self.__gray_background_contigs(contig_beginning_end, -0.02, 4, ax1)
            lines_prob = []
            for i in range(len(x)):
                lines_prob.append([(x[i][0], y_responsibilities_normal[i][0]),
                                   (x[i][1], y_responsibilities_normal[i][1])])
            lc_prob = LineCollection(lines_prob, colors=(1.0, 0.0, 0.0), linewidth=1.2)
            lc_cr_median = LineCollection(lines_cr_median, colors=y_color, linewidth=1.2)
            ax1.add_collection(lc_prob)
            ax1.add_collection(lc_cr_median)
            ax1.set_xlim(contig_beginning_end[0][0], contig_beginning_end[-1][1])
            ax1.set_ylim([-0.02, 3])
            ax1.set_ylabel('Copy ratio and Prob(normal)')
            ax1.plot([], color=(1.0, 0.0, 0.0), label='Prob(normal)')
            ax1.plot([], color=(0.0, 1.0, 0.0), label='Normal')
            ax1.plot([], color=(0.0, 0.0, 1.0), label='Deletion')
            ax1.plot([], color=(1.0, 0.0, 1.0), label='Amplification')
            ax1.plot([], color=(0.4, 0.4, 0.4), label='CNLOH')
            ax1.legend(loc="upper right", bbox_to_anchor=(1,1))

            self.__gray_background_contigs(contig_beginning_end, -0.02, 1.02, ax2)
            lc_prob = LineCollection(lines_prob, color=(1.0, 0.0, 0.0), linewidth=1.2)
            lc_min_af_median = LineCollection(lines_min_af_median, colors=y_color, linewidth=1.2)
            lc_maj_af_median = LineCollection(lines_maj_af_median, colors=y_color, linewidth=1.2)
            ax2.add_collection(lc_prob)
            ax2.add_collection(lc_min_af_median)
            ax2.add_collection(lc_maj_af_median)
            ax2.set_xlim(contig_beginning_end[0][0], contig_beginning_end[-1][1])
            ax2.set_ylim([-0.02, 1.02])
            ax2.set_ylabel('Allele fraction and Prob(normal)')
            ax2.set_xlabel('Genomic location (bp)')
            ax2.plot([], color=(1.0, 0.0, 0.0), label='Prob(normal)')
            ax2.plot([], color=(0.0, 1.0, 0.0), label='Normal')
            ax2.plot([], color=(0.0, 0.0, 1.0), label='Deletion')
            ax2.plot([], color=(1.0, 0.0, 1.0), label='Amplification')
            ax2.plot([], color=(0.4, 0.4, 0.4), label='CNLOH')
            ax2.legend(loc="upper right", bbox_to_anchor=(1,1))

            plt.savefig(self.__fig_del_ampl_filename)
            plt.close(fig2)

    def __plot_clustering(self):
        """ This function plots the details of how the fitting was done to the normal cluster
            in copy ratio and allele fraction space.
        """
        fig = plt.figure(1, figsize=(12, 6), dpi=400)
        plt.subplot(221)
        xedges = np.linspace(0, 5, 50)
        yedges = np.linspace(0, 0.5, 20)
        plt.hist2d(x=self.__copy_ratio_median, y=self.__allele_fraction_median,
                   weights=[2 / (1/self.__weights[i][0]**2 + 1/self.__weights[i][1]**2) for i in range(len(self.__weights))],
                   bins=(xedges, yedges), norm=LogNorm(), cmap="gray_r")
        plt.colorbar()
        plt.xlabel("Copy ratio")
        plt.ylabel("Allele fraction")
        plt.title("Data histogram")

        plt.subplot(222)
        self.__plot_classification()
        plt.xlabel("Copy ratio")
        plt.ylabel("Allele fraction")
        plt.title("Gaussian fit to the data")
        plt.xlim((0, 5))
        plt.ylim((0, 0.5))
        plt.xticks(np.arange(0, 6, 1))
        plt.yticks(np.arange(0, 0.6, 0.1))

        plt.subplot(223)
        self.__plot_Gaussian_mixture_fit()
        plt.xlabel("Copy ratio")
        plt.ylabel("Allele fraction")
        plt.title("Classification of segments")
        plt.xlim((0, 5))
        plt.ylim((0, 0.5))
        plt.xticks(np.arange(0, 6, 1))
        plt.yticks(np.arange(0, 0.6, 0.1))

        plt.subplot(224)
        plt.hist(self.__copy_ratio_sampled, bins = "auto", color = "k")
        plt.xlabel("Copy ratio")
        plt.xlim((0, 5))
        plt.xticks(np.arange(0, 6, 1))

        plt.subplots_adjust(hspace=0.5)
        plt.savefig(self.__fig_scatter_plot)
        plt.close(fig)

    def __plot_classification(self):
        """ Plot the medians of the distributions in copy ratio and allele fraction space
            and indicate the normal and not normal segments with different colors.
        """
        samples = np.asarray([self.__copy_ratio_median, self.__allele_fraction_median]).T
        colors = ["k" if i in self.__normal_segment_indices else "r" for i in range(len(samples))]
        plt.scatter(samples[:,0], samples[:,1], c=colors, alpha=0.8, s=10)

    def __plot_Gaussian_mixture_fit(self):
        """ Plot the fit of Gaussian clusters to the 2D copy ratio and allele fraction data.
        """
        samples = np.asarray([self.__copy_ratio_sampled, self.__allele_fraction_sampled]).T
        [pis, mu, cov] = self.__gaussian_mixture_fit
        n_Gaussians = len(pis)
        X, Y = np.mgrid[0:5:.05, 0:0.5:0.01]
        pos = np.empty(X.shape + (2,))
        pos[:, :, 0] = X
        pos[:, :, 1] = Y
        Z = 0 * X
        for i in range(n_Gaussians):
            rv = multivariate_normal(mu[i], cov[i])
            Z += pis[i] * rv.pdf(pos)
        Z = np.asarray([[(Z[i][j] + 0.0000000000001) for i in range(len(Z))] for j in range(len(Z[0]))]).T
        print("Z = ", Z)
        plt.contourf(X, Y, Z, cmap="gray_r", locator=ticker.LogLocator())
        plt.colorbar()

    def __gray_background_contigs(self, contig_beginning_end: List, ymin: float, ymax: float, ax):
        """Create a gray/white alternating background, in which the length of the individual stripes is
           proportional to the corresponding contig.
        """

        # The array contig_beginning_end should consist of pairs that contain the beginning and the end
        # of each chromosome
        n_contigs = len(contig_beginning_end)
        x = [x0 for contig in contig_beginning_end for x0 in contig]
        x = [x, x]
        y = [[ymin] * (2*n_contigs), [ymax] * (2*n_contigs)]
        z = []
        for i in range(n_contigs):
            z.append(i%2)
            z.append(i%2)
        z = [z]

        plt.xlim(contig_beginning_end[0][0], contig_beginning_end[-1][1])
        plt.ylim(ymin, ymax)
        ax.pcolor(x, y, z, cmap="gray", vmin=-20, vmax=1)

    def __get_phred_score(self, probability: float, max_phred_score: float=float("inf")):
        """ Get the phred score corresponding to the given probability.
        """
        if probability <= 0:
            return 0
        if probability >=1:
            return max_phred_score
        return min([-10*np.log10(1-probability), max_phred_score])

    def __save_calls_to_file(self):
        """ Save the calls and the corresponding PHRED scores into file.
        """
        input_filename = self.__cr_af_data.get_input_filename()
        [avg_normal_cr, std_dev_normal_cr] = self.__average_and_std_dev_copy_ratio_normal_segments()
        file = open(input_filename, "r")
        lines = file.readlines()

        file_header = ""
        file_data = ""
        i = 0
        for line in lines:
            values = line.strip().split()
            if not is_number(values[0]):
                file_header += line
            else:
                if self.__load_cr and self.__load_af:
                    if(len(values) >= 5
                       and not math.isnan(float(values[5]))
                       and not math.isnan(float(values[6]))
                       and not math.isnan(float(values[7]))
                       and not math.isnan(float(values[8]))
                       and not math.isnan(float(values[9]))
                       and not math.isnan(float(values[10]))
                       and 2**float(values[5]) <= self.__max_copy_ratio_possible
                       and 2**float(values[6]) <= self.__max_copy_ratio_possible
                       and 2**float(values[7]) <= self.__max_copy_ratio_possible
                       and 0. <= float(values[8]) <= 0.5
                       and 0. <= float(values[9]) <= 0.5
                       and 0. <= float(values[10]) <= 0.5
                       ):
                        file_data += line.strip() + "\t" + str(self.__normal_del_ampl(self.__copy_ratio_median[i],
                                                                                      avg_normal_cr, std_dev_normal_cr,
                                                                                      self.__responsibilities_normal[i]))
                        file_data += "\t"
                        file_data += str(self.__get_phred_score(probability=self.__responsibilities_normal[i],
                                                                max_phred_score=self.__max_PHRED_score))
                        file_data += "\n"
                        i += 1
                elif self.__load_cr:
                    if(len(values) >= 5
                       and not math.isnan(float(values[5]))
                       and not math.isnan(float(values[6]))
                       and not math.isnan(float(values[7]))
                       and 2**float(values[5]) <= self.__max_copy_ratio_possible
                       and 2**float(values[6]) <= self.__max_copy_ratio_possible
                       and 2**float(values[7]) <= self.__max_copy_ratio_possible
                       ):
                        file_data += line.strip() + "\t" + str(self.__normal_del_ampl(self.__copy_ratio_median[i],
                                                                                      avg_normal_cr, std_dev_normal_cr,
                                                                                      self.__responsibilities_normal[i]))
                        file_data += "\t"
                        file_data += str(self.__get_phred_score(probability=self.__responsibilities_normal[i],
                                                                max_phred_score=self.__max_PHRED_score))
                        file_data += "\n"
                        i += 1
                elif self.__load_af:
                    if(len(values) >= 5
                       and not math.isnan(float(values[8]))
                       and not math.isnan(float(values[9]))
                       and not math.isnan(float(values[10]))
                       and 0. <= float(values[8]) <= 0.5
                       and 0. <= float(values[9]) <= 0.5
                       and 0. <= float(values[10]) <= 0.5
                       ):
                        file_data += line.strip() + "\t" + str(self.__normal_del_ampl(self.__copy_ratio_median[i],
                                                                                      avg_normal_cr, std_dev_normal_cr,
                                                                                      self.__responsibilities_normal[i]))
                        file_data += "\t"
                        file_data += str(self.__get_phred_score(probability=self.__responsibilities_normal[i],
                                                                max_phred_score=self.__max_PHRED_score))
                        file_data += "\n"
                        i += 1
        file_header = file_header[0:-1] + "\tCALL\tPHRED\n"
        output_file = open(self.__output_calls_filename, "w")
        output_file.write(file_header + file_data)
        output_file.close()

    def __normal_del_ampl(self, cr: float, avg_normal_cr_: float, std_dev_normal_cr_: float,
                          responsibility_normal: float):
        """ Determines if the copy ratio value cr corresponds to a deletion (-),
            an amplification (+), a normal segment (0), a normal copy ratio
            with imbalanced allele fraction (CNLOH) or cannot be decided but not
            normal (/).
        """
        if responsibility_normal > self.__responsibility_threshold_normal:
            return "0"
        else:
            if self.__load_cr:
                if cr > avg_normal_cr_ + 2 * std_dev_normal_cr_:
                    return "+"
                elif cr < avg_normal_cr_ - 2 * std_dev_normal_cr_:
                    return "-"
                else:
                    return "CNLOH"
            else:
                return "/"

    def __average_and_std_dev_copy_ratio_normal_segments(self, responsibility_threshold: float=0.5):
        """ Auxiliary function to determine the averge and the standard deviation of the copy ratio
            value of all the normal segments.
        """
        total_length = 0
        total_cr = 0
        total_var_cr = 0
        for i in range(len(self.__copy_ratio_median)):
            if self.__responsibilities_normal[i] > responsibility_threshold:
                total_length += self.__n_points_cr[i]
                total_cr += self.__n_points_cr[i] * self.__copy_ratio_median[i]

        if total_length == 0:
            avg_normal_cr = 0
        else:
            avg_normal_cr = total_cr / total_length

        for i in range(len(self.__copy_ratio_median)):
            if self.__responsibilities_normal[i] > responsibility_threshold:
                total_var_cr += self.__n_points_cr[i] * (self.__copy_ratio_median[i] - avg_normal_cr)**2

        if total_length == 0:
            var_normal_cr = 0
        else:
            var_normal_cr = total_var_cr / total_length

        std_dev_normal_cr = var_normal_cr**0.5

        return [avg_normal_cr, std_dev_normal_cr]


