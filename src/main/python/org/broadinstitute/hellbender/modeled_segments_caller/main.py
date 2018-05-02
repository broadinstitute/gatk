import numpy as np
import scipy
from scipy.stats import multivariate_normal
from scipy.optimize import minimize
from sklearn import mixture, cluster
import math
import random
from matplotlib.collections import LineCollection
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import  ticker
from copy import deepcopy
from typing import List
import logging
import datetime

RANDOM_SEED = 123


def is_number(s):
    """ Decide if 's' is a number. """
    try:
        float(s)
        return True
    except ValueError:
        return False


class LoadAndSampleCrAndAf:
    """ Class to load the copy ratio and allele fraction data from .seg files.
        Both the copy ratio and the allele fraction posterior distribution is
        characterized by the distribution's median, 10th and 90th percentile.
        We read these data from the file on initialization and then sample
        the distribution using self.sample_data_points().
    """
    def __init__(self, filename: str, load_CR: bool = True, load_AF: bool=True, log_filename_prefix: str=""):
        """ Inputs:
            - filename: '.seg' file characterizing the posterior distribution of the
              copy ratio and the allele fraction
            - load_CR: whether to load copy ratio data
            - loasd_AF: whether to load the allele fraction data
        """

        # Start logging
        if not (log_filename_prefix==""):
            now = str(datetime.datetime.now())
            now2 = now.replace(" ", "_")
            now2 = now2.replace(":", "_")
            now2 = now2.replace(".", "_")
            self.log_filename = log_filename_prefix + now2 + ".log"
            logging.basicConfig(filename=self.log_filename,level=logging.DEBUG)
            logging.info("* ---- %s ---- *" % now)
            logging.info("Initializing class LoadAndSampleCrAndAf")
            logging.info("Input file: %s" % filename)
            logging.info("load_copy_ratio: %s" % load_CR)
            logging.info("load_allele_fraction: %s" % load_AF)
        else:
            self.log_filename=""

        # Initialize
        np.random.seed(RANDOM_SEED)
        self.__filename = filename
        self.__load_CR = load_CR
        self.__load_AF = load_AF

        # Load data from file
        if not (log_filename_prefix==""):
            logging.info("Loading data from file.")
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
         self.__n_points_CR,                 # number of copy ratio points in each segment
         self.__n_points_AF                  # number of allele fraction in each segment
         ] = self.__import_copy_ratio_and_allele_fraction()

        # Sample posterior distributions
        # Returns list of samples for each segment
        if not (log_filename_prefix==""):
            logging.info("Sampling data points.")
        [self.__copy_ratio_sampled_per_segment,       # sampled copy ratio values
         self.__allele_fraction_sampled_per_segment   # sampled allele fraction values
         ] = self.sample_data_points()

        # Flattened lists
        self.__copy_ratio_sampled = [i for sublist in self.__copy_ratio_sampled_per_segment for i in sublist]
        self.__allele_fraction_sampled = [i for sublist in self.__allele_fraction_sampled_per_segment for i in sublist]

        if not (log_filename_prefix==""):
            logging.info("Finished.\n\n\n")

    def get_data_from_file(self):
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
                self.__n_points_CR,
                self.__n_points_AF]

    def get_input_filename(self):
        """Return input filename."""
        return self.__filename

    def get_sampled_points(self):
        """Return sampled posterior distributions as flattened lists."""
        return [self.__copy_ratio_sampled, self.__allele_fraction_sampled]

    def get_sampled_points_per_segment(self):
        """Return samples from the posterior distributions of the copy ratio and the allele fraction."""
        return [self.__copy_ratio_sampled_per_segment, self.__allele_fraction_sampled_per_segment]

    def get_load_CR(self):
        """ Whether we loaded the copy ratio values."""
        return self.__load_CR

    def get_load_AF(self):
        """ Whether we loaded the allele fraction values."""
        return self.__load_AF

    def __import_copy_ratio_and_allele_fraction(self):
        """ Read data related to the copy ratio and allele fraction distributions from file.
        """
        copy_ratio_median = []         # median of the copy ratio distribution for the segment
        copy_ratio_10th_perc = []      # 10th percentile of the copy ratio distribution
        copy_ratio_90th_perc = []      # 90th percentile of the copy ratio distribution
        allele_fraction_median = []    # median of the allele fraction distribution
        allele_fraction_10th_perc = [] # 10th percentile of the allele fraction distribution
        allele_fraction_90th_perc = [] # 90th percentile of the allele fraction distribution
        contig = []                    # contig (which chromosome the segments are on)
        segment_start = []             # starting positions of the segments in base pair
        segment_end = []               # ending positions of the segments in base pair
        segment_lengths = []           # length of segments in base pair
        n_points_CR = []               # number of copy ratio data points taken from the segment
        n_points_AF = []               # number of allele fraction data points taken from the segment

        # Sample data points from each segment
        copy_ratio_sampled = []       # sampled copy ratio
        allele_fraction_sampled = []  # sampled allele fraction

        # Read data from file
        file = open(self.__filename, "r")
        lines = file.readlines()

        if self.__load_CR and self.__load_AF:
            for line in lines:
                values = line.strip().split()
                if is_number(values[0]):
                    if(len(values) >= 5):
                        if not float(values[4]) == 0 and not math.isnan(float(values[9])):
                            copy_ratio_median.append(2**float(values[6]))
                            copy_ratio_10th_perc.append(2**float(values[5]))
                            copy_ratio_90th_perc.append(2**float(values[7]))
                            allele_fraction_median.append(float(values[9]))
                            allele_fraction_10th_perc.append(float(values[8]))
                            allele_fraction_90th_perc.append(float(values[10].rstrip("\\")))
                            contig.append(int(values[0]))
                            segment_start.append(int(values[1]))
                            segment_end.append(int(values[2]))
                            segment_lengths.append(int(values[2]) - int(values[1]))
                            n_points_CR.append(int(values[3]))
                            n_points_AF.append(int(values[4]))
        elif self.__load_CR:
            for line in lines:
                values = line.strip().split()
                if is_number(values[0]):
                    if(len(values) >= 5):
                        copy_ratio_median.append(2**float(values[6]))
                        copy_ratio_10th_perc.append(2**float(values[5]))
                        copy_ratio_90th_perc.append(2**float(values[7]))
                        allele_fraction_median.append(float('nan'))
                        allele_fraction_10th_perc.append(float('nan'))
                        allele_fraction_90th_perc.append(float('nan'))
                        contig.append(int(values[0]))
                        segment_start.append(int(values[1]))
                        segment_end.append(int(values[2]))
                        segment_lengths.append(int(values[2]) - int(values[1]))
                        n_points_CR.append(int(values[3]))
                        n_points_AF.append(0)
        elif self.__load_AF:
            for line in lines:
                values = line.strip().split()
                if is_number(values[0]):
                    if(len(values) >= 5):
                        copy_ratio_median.append(float('nan'))
                        copy_ratio_10th_perc.append(float('nan'))
                        copy_ratio_90th_perc.append(float('nan'))
                        allele_fraction_median.append(float(values[9]))
                        allele_fraction_10th_perc.append(float(values[8]))
                        allele_fraction_90th_perc.append(float(values[10].rstrip("\\")))
                        contig.append(int(values[0]))
                        segment_start.append(int(values[1]))
                        segment_end.append(int(values[2]))
                        segment_lengths.append(int(values[2]) - int(values[1]))
                        n_points_CR.append(0)
                        n_points_AF.append(int(values[4]))

        return [copy_ratio_median, copy_ratio_10th_perc,
                copy_ratio_90th_perc, allele_fraction_median,
                allele_fraction_10th_perc, allele_fraction_90th_perc,
                contig, segment_start, segment_end, segment_lengths,
                n_points_CR, n_points_AF]

    def __fit_beta_distribution(self, af_median: float, af_10: float, af_90: float):
        """ Fit the beta distribution to the allele fraction distribution,
            when the 50th (median), 10th and 90th percentile are known.
        """

        def __find_beta_percentile(a: float, b: float, perc: float):
            # Get the value corresponding to percentile perc from
            # the beta distriubutionwith parameters a and b
            if perc < 0.:
                return 0.
            if perc > 1.:
                return 1.
            return scipy.special.betaincinv(a, b, perc)

        def __find_initial_ab(af_median: float, af_10: float, af_90: float):
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
            std_dev_estimate = 0.5 * (af_90 - af_10)
            var_estimate = std_dev_estimate ** 2
            a0 = (af_median**2 * (1 - af_median)) / var_estimate - af_median
            b0 = a0 * (1 - af_median) / af_median
            a0 = min(max(0.001, a0), 10000)
            b0 = min(max(0.001, b0), 10000)
            return [a0, b0]

        def __beta_err(ab):
            """ This function defines the fitting error. We require that the median
                and the width of the distribution (difference between the 10th and
                90th percentile) be as accurate as possible. """
            # Note, that we do not take the errors of the 10th and 90th percentile,
            # only the error of determining their difference.
            # The reason for this choice is that the errors in the 10th and 90th
            # percentile can be much larger than those of the median, as these are
            # at the tails of the distribution.

            a = ab[0]
            b = ab[1]
            af_median_prediction = __find_beta_percentile(a, b, 0.5)
            af_10_prediction = __find_beta_percentile(a, b, 0.1)
            af_90_prediction = __find_beta_percentile(a, b, 0.9)

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


    def __fit_gamma_distribution(self, cr_median: float, cr_10: float, cr_90: float):
        """ Fit the gamma distribution to the CR distribution data, given
            the median, the 10th and the 90th percentile of the distribution.
        """

        def __find_gamma_percentile(a, b, perc):
            # Get the copy ratio value corresponding to the percentile 'perc'
            if perc < 0:
                return 0
            if perc > 1:
                return 1
            return scipy.special.gammaincinv(a, perc) / b

        def __find_initial_ab(cr_median: float, cr_10: float, cr_90: float):
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
            std_dev_estimate = 0.5 * (cr_90 - cr_10)
            var_estimate = std_dev_estimate ** 2
            mean_estimate = 0.3 * cr_median + 0.7 * 0.5 * (cr_10 + cr_90)

            # The choice of mixing the known parameters is slightly arbitrary,
            # but it approximates the mean relatively well
            b0 = mean_estimate / var_estimate
            a0 = mean_estimate * b0
            return [a0, b0]

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
            a = ab[0]
            b = ab[1]
            cr_10_prediction = __find_gamma_percentile(a, b, 0.1)
            cr_90_prediction = __find_gamma_percentile(a, b, 0.9)
            cr_median_prediction = __find_gamma_percentile(a, b, 0.5)
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

    def sample_data_points(self):
        """ This function samples the posterior probability distribution of each segment,
            on average using 'avg_n_sampled_points_per_segment' points per segment
        """
        avg_n_sampled_points_per_segment = 12

        n_segments = len(self.__copy_ratio_median)
        total_length = np.sum(self.__segment_lengths)
        avg_segment_length = max(1, int(total_length / n_segments))

        copy_ratio_sampled_per_segment = []
        allele_fraction_sampled_per_segment = []

        for i in range(n_segments):
            n_sampled_points = max([1, int(avg_n_sampled_points_per_segment*self.__segment_lengths[i] / avg_segment_length)])
            [cr, af] = self.sample_data_points_from_one_segment(
                self.__copy_ratio_median[i],
                self.__copy_ratio_10th_perc[i],
                self.__copy_ratio_90th_perc[i],
                self.__allele_fraction_median[i],
                self.__allele_fraction_10th_perc[i],
                self.__allele_fraction_90th_perc[i],
                n_sampled_points
            )
            copy_ratio_sampled_per_segment.append(cr)
            allele_fraction_sampled_per_segment.append(af)

        return [copy_ratio_sampled_per_segment, allele_fraction_sampled_per_segment]

    def sample_data_points_from_one_segment(self, cr_median: float, cr_10: float, cr_90: float,
                                            af_median: float, af_10: float, af_90: float, n_points: int):
        """ Sample n_points data points from each segment, according to the gamma
            and beta distributions whose 10th, 50th and 90th percentiles are given.
        """

        # Notice that when we fit the beta distribution then the x axis needs
        # to be expanded two-fold, so that instead of the usual [0, 0.5] range for AF,
        # we get a [0,1] range of the beta distribution.
        if self.__load_CR:
            [cr_a, cr_b] = self.__fit_gamma_distribution(cr_median, cr_10, cr_90)
            cr = np.random.gamma(cr_a, 1 / cr_b, n_points)
        else:
            cr = [float('nan')] * n_points
        if self.__load_AF:
            [af_a, af_b] = self.__fit_beta_distribution(2 * af_median, 2 * af_10, 2 * af_90)
            af = 0.5 * np.random.beta(af_a, af_b, n_points)
        else:
            af = [float('nan')] * n_points
        return [cr, af]





class ModeledSegmentsCaller:
    """ This class takes an instance of the class LoadAndSampleCrAndAf as input
        and it identifies the normal segments. It outputs the plots of the segments
        (and optionally other plots about the fitting procedure) and saves the results
        into files.
    """
    def __init__(self, CR_AF_data: LoadAndSampleCrAndAf, interactive: bool=False,
                 output_image_dir: str= "",
                 output_calls_dir: str= "",
                 output_image_prefix:str= "",
                 output_calls_prefix: str= "",
                 output_image_suffix: str=".jpg",
                 output_calls_suffix: str=".called.seg",
                 interactive_output_del_ampl_image_suffix:str="_del_ampl.jpg",
                 interactive_output_scatter_plot_suffix:str="_scatter_plot.jpg",
                 interactive_output_allele_fraction_plot_suffix:str= "_allele_fraction_CN1_and_CN2_candidate_intervals.jpg",
                 interactive_output_copy_ratio_suffix:str="_copy_ratio_fit.jpg",
                 interactive_output_copy_ratio_clustering_suffix:str="_copy_ratio_clusters.jpg"):
        """ On initialization, the caller loads the copy ratio and allele fraction data from
            the LoadAndSampleCrAndAf class. It then identifies normal segments and saves the
            results, including the plots.
        """

        # Start logging
        self.log_filename=CR_AF_data.log_filename
        if not (self.log_filename==""):
            logging.info("* ---- %s ---- *" % str(datetime.datetime.now()))
            logging.info("Initializing class CNVCaller")

        # Copy input data
        self.__CR_AF_data = CR_AF_data
        self.__load_CR = CR_AF_data.get_load_CR()
        self.__load_AF = CR_AF_data.get_load_AF()
        self.__interactive = interactive
        self.__output_image_dir = output_image_dir
        self.__output_calls_dir = output_calls_dir
        self.__output_calls_prefix = output_calls_prefix
        self.__output_image_prefix = output_image_prefix
        self.__output_calls_suffix = output_calls_suffix
        self.__output_image_suffix = output_image_suffix
        self.__interactive_output_del_ampl_image_suffix = interactive_output_del_ampl_image_suffix
        self.__interactive_output_scatter_plot_suffix = interactive_output_scatter_plot_suffix
        self.__interactive_output_allele_fraction_plot_suffix = interactive_output_allele_fraction_plot_suffix
        self.__interactive_output_copy_ratio_suffix = interactive_output_copy_ratio_suffix
        self.__interactive_output_copy_ratio_clustering_suffix = interactive_output_copy_ratio_clustering_suffix

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
         self.__n_points_CR,
         self.__n_points_AF] = CR_AF_data.get_data_from_file()

        # Sample copy ratio and allele fraction points according to the distribution specified in
        # the input file. Number of points sampled is proportional to the segment lengths
        [self.__copy_ratio_sampled,
         self.__allele_fraction_sampled] = CR_AF_data.get_sampled_points()
        [self.__copy_ratio_sampled_per_segment,
         self.__allele_fraction_sampled_per_segment] = self.__CR_AF_data.get_sampled_points_per_segment()

        # Set the filenames for the output data
        if not self.log_filename=="":
            logging.info("Setting output filenames.")
        self.set_output_filenames()
        if not self.log_filename=="":
            logging.info("   Normal segments image file : %s" % self.fig_normal_segments_filename)
            logging.info("   Calls file : %s" % self.output_calls_filename)
            if self.__interactive:
                logging.info("   (interactive mode) Image of deletions and amplifications: %s"
                             % self.fig_del_ampl_filename)
                logging.info("   (interactive mode) Scatter plot of segments: %s" % self.fig_scatter_plot)
                logging.info("   (interactive mode) Plots of allele fraction data with copy ratios in the " +
                             "copy number 1 and 2 candidate intervals: %s"
                             % self.fig_allele_fraction_CN1_CN2_intervals)
                logging.info("   (interactive mode) Gaussian fit to the copy ratio data: %s" % self.fig_copy_ratio_fit)

        # Find normal segments
        if not self.log_filename=="":
            logging.info("Determining normal segments.")
        if self.__load_CR and self.__load_AF:
            [self.__responsibilities_normal,
             self.__normal_peak_indices,
             self.__gaussian_mixture_fit] = self.__choose_normal_segments__CR_AF_data()
        elif self.__load_CR:
            [self.__responsibilities_normal,
             self.__normal_segment_indices] = self.__choose_normal_segments__CR_data_only()

        # Save plots of the segments
        if not self.log_filename=="":
            logging.info("Plotting and saving segments.")
        self.__plot_and_save_segments()

        # Create auxiliary plots when needed
        if self.__load_CR and self.__load_AF and self.__interactive:
            if not self.log_filename=="":
                logging.info("Creating auxiliary plots in interactive mode.")
            self.__plot_clustering()

        # Save the results in a file
        if not self.log_filename=="":
            logging.info("Saving results.")
        self.__save_calls_to_file()
        if not self.log_filename=="":
            logging.info("Finished.\n\n\n")


    def set_output_filenames(self):
        """Set file names for all outputs."""
        self.fig_normal_segments_filename = self.__output_image_dir + self.__output_image_prefix + self.__output_image_suffix
        self.output_calls_filename = self.__output_calls_dir + self.__output_calls_prefix + self.__output_calls_suffix

        # Get input directory
        input_fname_ending = self.__CR_AF_data.get_input_filename().split("/")[-1]
        input_dir = self.__CR_AF_data.get_input_filename().split(input_fname_ending)[0]
        if input_dir[-1] != "/":
            input_dir = input_dir + "/"

        if self.__output_image_dir == "":
            self.__output_image_dir = input_dir

        if self.__output_calls_dir == "":
            self.__output_calls_dir = input_dir

        # Set file names for images produced in interactive mode
        if self.__output_image_prefix == "":
            image_fname_prefix = self.__CR_AF_data.get_input_filename()
            image_fname_extension = image_fname_prefix.split(".")[-1]
            image_fname_prefix = image_fname_prefix.split("." + image_fname_extension)[0]
        else:
            image_fname_prefix = self.__output_image_prefix

        # Set the names of the outputs used in interactive mode
        image_fname_prefix = self.__output_image_dir +  image_fname_prefix
        self.fig_del_ampl_filename = image_fname_prefix + self.__interactive_output_del_ampl_image_suffix
        self.fig_scatter_plot = image_fname_prefix + self.__interactive_output_scatter_plot_suffix
        self.fig_allele_fraction_CN1_CN2_intervals = image_fname_prefix + self.__interactive_output_allele_fraction_plot_suffix
        self.fig_copy_ratio_fit = image_fname_prefix + self.__interactive_output_copy_ratio_suffix
        self.fig_copy_ratio_clusters = image_fname_prefix + self.__interactive_output_copy_ratio_clustering_suffix

    def get_output_calls_filename(self):
        """Name of the output calls file."""
        return self.output_calls_filename

    def get_all_output_fig_filenames(self):
        """Name of all filenames output by the tool."""
        if not self.__interactive:
            return [self.fig_normal_segments_filename,
                    self.fig_del_ampl_filename,
                    self.fig_scatter_plot,
                    self.fig_allele_fraction_CN1_CN2_intervals,
                    self.fig_copy_ratio_fit]
        else:
            return [self.fig_normal_segments_filename]

    def get_output_fig_filename(self):
        """Name of the output fig file."""
        return self.fig_normal_segments_filename

    def get_interactive_figure_filenames(self):
        """Get the names of auxiliary files that we save in interactive mode."""
        if not self.__interactive:
            return None
        return [self.fig_del_ampl_filename,
                self.fig_scatter_plot,
                self.fig_allele_fraction_CN1_CN2_intervals,
                self.fig_copy_ratio_fit]

    def get_normal_segment_indices(self):
        """Output the indices of the normal segments."""
        return self.__normal_segment_indices

    def get_responsibilities_segments_are_normal(self):
        """Output the responsibility associated with each segment being normal."""
        return self.__responsibilities_normal

    def get_interactive(self):
        """Return whether the code was run in interactive mode."""
        return self.__interactive

    def __estimate_width_of_normal_cluster(self):
        """ This auxiliary function estimates the standard deviation of the tightest cluster
            in the region where the normal sample's peaks are expected to be.
            This result is used later in the caller.
        """
        n_Gaussian_components_to_look_for = 5
        weight_cutoff = 0.10   # Gaussians w/ weight below this will not be considered
        default_covariances = [0.01, 0.01] # default return (if no points in the region)

        # Heuristic parameters: choose points only in the region specified here:
        cr_region = [0.4, 2.5]
        af_region = [0.4, 0.5]
        cr_af_normal = [[self.__copy_ratio_sampled[i], self.__allele_fraction_sampled[i]]
                        for i in range(len(self.__copy_ratio_sampled))
                        if cr_region[0] <= self.__copy_ratio_sampled[i] <= cr_region[1]
                        and af_region[0] <= self.__allele_fraction_sampled[i] <= af_region[1]
                        ]
        if len(cr_af_normal) <= 1:
            return default_covariances
        [cr_normal, af_normal] = np.asarray(cr_af_normal).T
        samples = np.asarray([cr_normal, af_normal]).T
        gmix = mixture.GaussianMixture(n_components = min([n_Gaussian_components_to_look_for,
                                       len(cr_af_normal)]), covariance_type="diag")
        gmix.fit(samples)
        covariances = np.array([gmix.covariances_[i] for i in range(len(gmix.weights_))
                                if gmix.weights_[i] >= weight_cutoff])
        covariances = [min(covariances[:,0]), min(covariances[:,1])]
        return(covariances)

    def __choose_normal_segments__CR_AF_data(self, n_Gaussians: int=40):
        """ Choose those Gaussians from the fitted distribution that cover the
            normal segments. The Gaussians are fitted using Bayesian variational
            inference in the  two dimensional space of copy ratio and allele fraction axes.
        """

        np.random.seed(RANDOM_SEED)
        #[sigma_CR_normal, sigma_AF_normal] = self.__estimate_width_of_normal_cluster()
        sigma_CR_normal = 0.005
        sigma_AF_normal = 0.005

        # Fit mixtures using Gaussian variational inference
        samples = np.asarray([self.__copy_ratio_sampled, self.__allele_fraction_sampled]).T
        gmix = mixture.BayesianGaussianMixture(
            weight_concentration_prior_type="dirichlet_distribution",
            covariance_type="diag", n_components=n_Gaussians,
            init_params="random", tol=1e-4, n_init=2,
            max_iter=1500,
            covariance_prior=np.asarray([sigma_CR_normal, sigma_AF_normal]),
            weight_concentration_prior=1000.0/n_Gaussians
        )
        gmix.fit(samples)

        # We choose those peaks to be normal whose mean's copy ratio value is within the range specified
        # by '__choose_CN2_CR_cluster' and whose allele fraction value is within the range specified by 'normal_range_AF'
        normal_range_AF = [0.475, 0.500]
        normal_range_CR = self.__choose_CN2_CR_cluster()

        normal_peak_indices = []
        for i in range(len(gmix.weights_)):
            if ((normal_range_CR[0]-np.sqrt(gmix.covariances_[i, 0]) <= gmix.means_[i, 0]
                <= normal_range_CR[1]+np.sqrt(gmix.covariances_[i, 0]))
                and (normal_range_AF[0] <= gmix.means_[i, 1] + np.sqrt(gmix.covariances_[i, 1]))
                and (gmix.means_[i, 1] <= normal_range_AF[1])):
                normal_peak_indices.append(i)

        classification = gmix.predict(samples)
        responsibilities_normal = self.__responsibility_segment_is_normal(gmix.weights_,
                                                                          gmix.means_,
                                                                          gmix.covariances_,
                                                                          classification,
                                                                          normal_peak_indices)
        return [responsibilities_normal, normal_peak_indices, gmix]

    def __choose_normal_segments__CR_data_only(self):
        """ This function determines which peak in the 1D copy ratio data is the copy number 2 peak.
        """
        n_peaks_initial_fit, _, _, _ = self.__estimate_number_of_CR_clusters(max_n_Gaussians = 10, alpha = 0.1,
                                                                             min_std_dev = 0.05)
        n_peaks_initial_fit, _, _, _ = self.__estimate_number_of_CR_clusters(max_n_Gaussians = n_peaks_initial_fit + 2,
                                                                             alpha = 0.8, min_std_dev = 0.03)
        n_peaks_initial_fit, w_peaks_initial_fit, mu_peaks_initial_fit, _ = self.__estimate_number_of_CR_clusters(
            max_n_Gaussians = n_peaks_initial_fit + 2,
            alpha = 0.5, min_std_dev = 0.05
        )

        ordering = self.__indices_increasing_order(list(mu_peaks_initial_fit))
        w_peaks_initial_fit  = [w_peaks_initial_fit[i]  for i in ordering]
        mu_peaks_initial_fit = [mu_peaks_initial_fit[i] for i in ordering]

        if n_peaks_initial_fit <= 1:
            CN2_interval = [0, 100000]
        elif n_peaks_initial_fit == 2:
            if w_peaks_initial_fit[0] >= w_peaks_initial_fit[1]:
                CN2_interval = [0, 0.5 * (mu_peaks_initial_fit[0] + mu_peaks_initial_fit[1])]
            else:
                CN2_interval = [0.5 * (mu_peaks_initial_fit[0] + mu_peaks_initial_fit[1]), 100000]
        else:
            if w_peaks_initial_fit[0] >= w_peaks_initial_fit[1]:
                CN2_interval = [0, 0.5 * (mu_peaks_initial_fit[0] + mu_peaks_initial_fit[1])]
            else:
                CN2_interval = [0.5 * (mu_peaks_initial_fit[0] + mu_peaks_initial_fit[1]),
                                0.5 * (mu_peaks_initial_fit[1] + mu_peaks_initial_fit[2])]

        responsibilities_normal = [0] * len(self.__copy_ratio_median)
        normal_peak_indices = []
        for i in range(len(self.__copy_ratio_median)):
            if CN2_interval[0] <= self.__copy_ratio_median[i] <= CN2_interval[1]:
                responsibilities_normal[i] = 1
                normal_peak_indices.append(i)
            else:
                responsibilities_normal[i] = 0

        return [responsibilities_normal, normal_peak_indices]


    def __responsibility_segment_is_normal(self, gmix_weights, gmix_means,
                                           gmix_covariances, gmix_classification, normal_peak_indices):
        """ Based on a Gaussian fit to the clusters in copy ratio and allele fraction space,
            figure out the responsibilities of a segment being normal.
        """
        n_Gaussians = len(gmix_weights)
        n_segments = len(self.__copy_ratio_median)

        # The Gaussian mixture model determined the Gaussian mixture model that fits the data the best
        # Now, we assign the probability of each segment being Gaussian to each segment.
        # This probability is defined as the 'responsibility' assigned to the segment median
        responsibilities_normal = [0] * n_segments
        for i in range(n_segments):
            pts = np.array([self.__copy_ratio_sampled_per_segment[i],
                            self.__allele_fraction_sampled_per_segment[i]]).T
            normal_pdf = np.sum([gmix_weights[j] * scipy.stats.multivariate_normal.pdf(
                pts[k], gmix_means[j], gmix_covariances[j])
                                 for j in normal_peak_indices
                                 for k in range(len(pts))]) / len(pts)
            total_pdf = np.sum([gmix_weights[j] * scipy.stats.multivariate_normal.pdf(
                pts[k], gmix_means[j], gmix_covariances[j])
                                for j in range(n_Gaussians)
                                for k in range(len(pts))]) / len(pts)
            responsibilities_normal[i] = max(min(normal_pdf / total_pdf, 1), 0)
        return responsibilities_normal

    def __choose_CN2_CR_cluster(self):
        """ This function determines which peak in the 1D copy ratio data is the copy number 2 peak.
            It also makes use of the allele fraction data to make this decision.
        """
        n_peaks_initial_fit, _, _, _ = self.__estimate_number_of_CR_clusters(max_n_Gaussians = 10, alpha = 0.1,
                                                                             min_std_dev = 0.05)
        n_peaks_initial_fit, _, _, _ = self.__estimate_number_of_CR_clusters(max_n_Gaussians = n_peaks_initial_fit + 2,
                                                                             alpha = 0.8, min_std_dev = 0.03)
        n_peaks_initial_fit, _, mu_peaks_initial_fit, _ = self.__estimate_number_of_CR_clusters(
                                                                             max_n_Gaussians = n_peaks_initial_fit + 2,
                                                                             alpha = 0.5, min_std_dev = 0.05)

        ordering = self.__indices_increasing_order(list(mu_peaks_initial_fit))
        mu_peaks_initial_fit = [mu_peaks_initial_fit[i] for i in ordering]

        if n_peaks_initial_fit > 1:
            d_mean = mu_peaks_initial_fit[1] - mu_peaks_initial_fit[0]
            o_mean = mu_peaks_initial_fit[0]
        else:
            d_mean = 100000
            o_mean = 0
            logging.info("Only a single copy ratio peak found. %s was not created." % self.fig_allele_fraction_CN1_CN2_intervals)
            return [o_mean, d_mean]

        n_intervals_needed = n_peaks_initial_fit + 2

        def __CR_clusters(distance, offset, n_intervals_needed):
            n0 = int(np.ceil((offset - distance/2) / distance))
            n_intervals = n_intervals_needed + n0
            first_interval_beginning = offset - (n0 + 1/2) * distance
            return [first_interval_beginning, n_intervals]

        # Find those intervals which correspond to the CN = 1 and CN = 2 copy ratio values
        [first_interval_beginning, n_intervals] = __CR_clusters(d_mean, o_mean, n_intervals_needed)
        p_CR_segment = [0] * n_intervals
        max_CR = first_interval_beginning + n_intervals * d_mean
        for CR in self.__copy_ratio_sampled:
            if CR < max_CR:
                i = int(np.floor((CR - first_interval_beginning) / d_mean))
                p_CR_segment[i] += 1
        p_CR_segment = [pcr / len(self.__copy_ratio_sampled) for pcr in p_CR_segment]

        CN1_interval_candidate_index = -1
        CN2_interval_candidate_index = -1
        for i in range(len(p_CR_segment)):
            if p_CR_segment[i] > 0.03:
                if CN1_interval_candidate_index < 0:
                    CN1_interval_candidate_index = i
                elif CN2_interval_candidate_index < 0:
                    CN2_interval_candidate_index = i
                else:
                    break

        CN1_interval_candidate = [0, first_interval_beginning + (CN1_interval_candidate_index + 1) * d_mean]
        CN2_interval_candidate = [first_interval_beginning + CN2_interval_candidate_index * d_mean,
                                  first_interval_beginning + (CN2_interval_candidate_index + 1) * d_mean]
        if self.__interactive:
            fig = plt.figure(1, dpi=300)
            plt.hist(self.__copy_ratio_sampled, color = "b", bins = 100, density = True)
            plt.plot(CN1_interval_candidate[0], [0], "rv")
            plt.plot(CN1_interval_candidate[1], [0], "rv")
            plt.plot(CN2_interval_candidate[0], [0], "g*")
            plt.plot(CN2_interval_candidate[1], [0], "g*")
            plt.xlim(0,3)
            plt.savefig(self.fig_copy_ratio_clusters)
            plt.close(fig)

        # If at least 'normal_min_ratio' fraction of the AF points in the 'CN1_interval_candidate' falls within
        # the range 'normal_range_AF', then we call that the copy number 2 segment. Otherwise, it will be
        # 'CN2_interval_candidate'.
        normal_min_ratio = 0.15
        normal_range_AF = [0.475, 0.50]
        AF__CN1_interval_candidate = []
        AF__CN2_interval_candidate = []
        for i in range(len(self.__allele_fraction_sampled)):
            if CN2_interval_candidate[0] <= self.__copy_ratio_sampled[i] < CN2_interval_candidate[1]:
                AF__CN2_interval_candidate.append(self.__allele_fraction_sampled[i])
            elif CN1_interval_candidate[0] <= self.__copy_ratio_sampled[i] < CN1_interval_candidate[1]:
                AF__CN1_interval_candidate.append(self.__allele_fraction_sampled[i])
        n_points_in_normal_range__CN1_interval_candidate = len([af for af in AF__CN1_interval_candidate
                                                                if (af >= normal_range_AF[0] and af <= normal_range_AF[1])])
        if (n_points_in_normal_range__CN1_interval_candidate / len(AF__CN1_interval_candidate) > normal_min_ratio):
            CN2_interval = CN1_interval_candidate
            if self.__interactive and self.log_filename!="":
                logging.info("We chose the 1st interval to be of copy number 2.")
        else:
            CN2_interval = CN2_interval_candidate
            if self.__interactive and self.log_filename!="":
                logging.info("We chose the 2nd interval to be of copy number 2.")

        # Plot histograms of AF data falling in the intervals 'CN1_interval_candidate' and 'CN2_interval_candidate'
        if self.__interactive:
            fig = plt.figure(2, dpi=400)
            plt.subplot(211)
            plt.hist(AF__CN1_interval_candidate, bins=50)
            plt.xlim(0,0.5)
            plt.xlabel("Allele fraction of the first interval")
            plt.subplot(212)
            plt.hist(AF__CN2_interval_candidate, bins=50)
            plt.xlim(0,0.5)
            plt.xlabel("Allele fraction of the second interval")
            plt.savefig(self.fig_allele_fraction_CN1_CN2_intervals)
            plt.close(fig)

        return CN2_interval

    def __indices_increasing_order(self, ls: List):
        """ Get the indices of the ordered version of ls.
        """
        ls_new = deepcopy(list(ls))
        ls_new.sort()
        ordering = [ls.index(i) for i in ls_new]
        return ordering

    def __fit_normal_1D(self, data: List, n_Gaussians: int, n_iterations: int):
        """ Fit Gaussians to the one dimensional data.
            n_Gaussians = the number of Gaussians used
            n_iterations: the number of scikit learn iterations used.
        """
        random.seed(RANDOM_SEED)
        gmix = mixture.BayesianGaussianMixture(n_components=n_Gaussians,
                                               weight_concentration_prior_type="dirichlet_distribution",
                                               weight_concentration_prior = 0.001 / n_Gaussians,
                                               covariance_prior = np.array([[0.01]]),
                                               max_iter = n_iterations)
        gmix.fit(np.asarray(data).reshape(-1,1))
        return [gmix.weights_, gmix.means_, gmix.covariances_]

    def __estimate_number_of_CR_clusters(self, max_n_Gaussians: int=10, alpha: float=1.,
                                         min_std_dev: float=0.03, max_std_dev:float =0.4):
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

        [w, mu, cov] = self.__fit_normal_1D(self.__copy_ratio_sampled, max_n_Gaussians, 10000)
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

        if self.__interactive:
            fig = plt.figure(1, dpi=400)
            x = np.linspace(0, 5, 1000)
            cols = ["r", "m", "c", "k", "y"]
            plt.hist(self.__copy_ratio_sampled, density=1, bins=200)
            for i in range(len(peak_groups)):
                y = [0] * len(x)
                for j in peak_groups[i]:
                    rv = multivariate_normal(mu[j], cov[j])
                    y = y + w[j] * rv.pdf(x)
                plt.plot(x, y, cols[i % len(cols)])
            plt.xlim(0,2.5)
            plt.savefig(self.fig_copy_ratio_fit)
            plt.close(fig)

        return n_peaks, w_peaks, mu_peaks, sd_peaks

    def __plot_and_save_segments(self):
        """ This function plots the results of the calling scripts and also saves
            these results.
        """
        [avg_normal_CR, std_dev_normal_CR] = self.__average_and_std_dev_copy_ratio_normal_segments()
        n_segments = len(self.__copy_ratio_median)

        x = []
        y_CR_median = []
        y_CR_10 = []
        y_CR_90 = []
        y_min_AF_median = []
        y_min_AF_10 = []
        y_min_AF_90 = []
        y_maj_AF_median = []
        y_maj_AF_10 = []
        y_maj_AF_90 = []
        y_responsibilities_normal = []
        y_responsibilities_normal_PHRED = []
        y_color = []

        site = 0
        for i in range(n_segments):
            x.append([site, site + self.__segment_lengths[i]])
            site += self.__segment_lengths[i]
            y_CR_median.append([self.__copy_ratio_median[i]] * 2)
            y_CR_10.append([self.__copy_ratio_10th_perc[i]] * 2)
            y_CR_90.append([self.__copy_ratio_90th_perc[i]] * 2)
            y_min_AF_median.append([self.__allele_fraction_median[i]] * 2)
            y_min_AF_10.append([self.__allele_fraction_10th_perc[i]] * 2)
            y_min_AF_90.append([self.__allele_fraction_90th_perc[i]] * 2)
            y_maj_AF_median.append([1-self.__allele_fraction_median[i]] * 2)
            y_maj_AF_10.append([1-self.__allele_fraction_10th_perc[i]] * 2)
            y_maj_AF_90.append([1-self.__allele_fraction_90th_perc[i]] * 2)
            y_responsibilities_normal.append([self.__responsibilities_normal[i]] * 2)
            y_responsibilities_normal_PHRED.append([self.__get_phred_score(
                probability=self.__responsibilities_normal[i],
                max_phred_score=self.__max_PHRED_score
            )] * 2)

            n_d_a = self.__normal_del_ampl(self.__copy_ratio_median[i], avg_normal_CR,
                                           std_dev_normal_CR, self.__responsibilities_normal[i])
            if n_d_a == "+":
                y_color.append((1.0, 0.0, 1.0))
            elif n_d_a == "0":
                y_color.append((0.0, 1.0, 0.0))
            elif n_d_a == "-":
                y_color.append((0.0, 0.0, 1.0))
            else:
                y_color.append((0.4, 0.4, 0.4))

        contig_beginning_end = []
        current_contig = self.__contig[0]
        current_contig_beginning = 1
        site = 0
        for i in range(n_segments):
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
        lines_CR_10 = []
        lines_CR_median = []
        lines_CR_90 = []
        lines_CR_color = []
        for i in range(len(x)):
            lines_CR_10.append([(x[i][0], y_CR_10[i][0]), (x[i][1], y_CR_10[i][1])])
            lines_CR_median.append([(x[i][0], y_CR_median[i][0]), (x[i][1], y_CR_median[i][1])])
            lines_CR_90.append([(x[i][0], y_CR_90[i][0]), (x[i][1], y_CR_90[i][1])])
            lines_CR_color.append((1-self.__responsibilities_normal[i], 0, 0))
        lc_CR_10 = LineCollection(lines_CR_10, colors=lines_CR_color, linewidth=0.8)
        lc_CR_median = LineCollection(lines_CR_median, colors=lines_CR_color, linewidth=1.2)
        lc_CR_90 = LineCollection(lines_CR_90, colors=lines_CR_color, linewidth=0.8)
        ax2.add_collection(lc_CR_10)
        ax2.add_collection(lc_CR_median)
        ax2.add_collection(lc_CR_90)
        ax2.set_xlim(contig_beginning_end[0][0], contig_beginning_end[-1][1])
        ax2.set_ylim([-0.02, 3])
        ax2.set_ylabel('Copy ratio')
        ax2.plot([], color='r', label='Not normal')
        ax2.plot([], color='k', label='Normal')
        ax2.legend(loc="upper right", bbox_to_anchor=(1,1))

        self.__gray_background_contigs(contig_beginning_end, -0.02, 1.02, ax3)
        lines_min_AF_10 = []
        lines_min_AF_median = []
        lines_min_AF_90 = []
        lines_maj_AF_10 = []
        lines_maj_AF_median = []
        lines_maj_AF_90 = []
        lines_AF_color = []
        for i in range(len(x)):
            lines_min_AF_10.append([(x[i][0], y_min_AF_10[i][0]), (x[i][1], y_min_AF_10[i][1])])
            lines_min_AF_median.append([(x[i][0], y_min_AF_median[i][0]), (x[i][1], y_min_AF_median[i][1])])
            lines_min_AF_90.append([(x[i][0], y_min_AF_90[i][0]), (x[i][1], y_min_AF_90[i][1])])
            lines_AF_color.append((1-self.__responsibilities_normal[i], 0, 0))
            lines_maj_AF_10.append([(x[i][0], y_maj_AF_10[i][0]), (x[i][1], y_maj_AF_10[i][1])])
            lines_maj_AF_median.append([(x[i][0], y_maj_AF_median[i][0]), (x[i][1], y_maj_AF_median[i][1])])
            lines_maj_AF_90.append([(x[i][0], y_maj_AF_90[i][0]), (x[i][1], y_maj_AF_90[i][1])])
        lc_min_AF_10 =LineCollection(lines_min_AF_10, colors=lines_AF_color, linewidth=0.8)
        lc_min_AF_median = LineCollection(lines_min_AF_median, colors=lines_AF_color, linewidth=1.2)
        lc_min_AF_90 = LineCollection(lines_min_AF_90, colors=lines_AF_color, linewidth=0.8)
        lc_maj_AF_10 = LineCollection(lines_maj_AF_10, colors=lines_AF_color, linewidth=0.8)
        lc_maj_AF_median = LineCollection(lines_maj_AF_median, colors=lines_AF_color, linewidth=1.2)
        lc_maj_AF_90 = LineCollection(lines_maj_AF_90, colors=lines_AF_color, linewidth=0.8)
        ax3.add_collection(lc_min_AF_10)
        ax3.add_collection(lc_min_AF_median)
        ax3.add_collection(lc_min_AF_90)
        ax3.add_collection(lc_maj_AF_10)
        ax3.add_collection(lc_maj_AF_median)
        ax3.add_collection(lc_maj_AF_90)
        ax3.set_xlim(contig_beginning_end[0][0], contig_beginning_end[-1][1])
        ax3.set_ylim([-0.02, 1.02])
        ax3.set_ylabel('Allele fraction')
        ax3.plot([], color='r', label='Not normal')
        ax3.plot([], color='k', label='Normal')
        ax3.legend(loc="upper right", bbox_to_anchor=(1,1))

        plt.savefig(self.fig_normal_segments_filename)
        plt.close(fig1)

        if self.__interactive:
            fig2, (ax1, ax2) = plt.subplots(2, sharex=True, dpi=400, figsize=(8, 6))
            self.__gray_background_contigs(contig_beginning_end, -0.02, 4, ax1)
            lines_prob = []
            for i in range(len(x)):
                lines_prob.append([(x[i][0], y_responsibilities_normal[i][0]),
                                   (x[i][1], y_responsibilities_normal[i][1])])
            lc_prob = LineCollection(lines_prob, colors=(1.0, 0.0, 0.0), linewidth=1.2)
            lc_CR_median = LineCollection(lines_CR_median, colors=y_color, linewidth=1.2)
            ax1.add_collection(lc_prob)
            ax1.add_collection(lc_CR_median)
            ax1.set_xlim(contig_beginning_end[0][0], contig_beginning_end[-1][1])
            ax1.set_ylim([-0.02, 3])
            ax1.plot([], color=(1.0, 0.0, 0.0), label='Prob(normal)')
            ax1.plot([], color=(0.0, 1.0, 0.0), label='Normal')
            ax1.plot([], color=(0.0, 0.0, 1.0), label='Deletion')
            ax1.plot([], color=(1.0, 0.0, 1.0), label='Amplification')
            ax1.plot([], color=(0.4, 0.4, 0.4), label='CNLOH')
            ax1.legend(loc="upper right", bbox_to_anchor=(1,1))

            self.__gray_background_contigs(contig_beginning_end, -0.02, 1.02, ax2)
            lc_prob = LineCollection(lines_prob, color=(1.0, 0.0, 0.0), linewidth=1.2)
            lc_min_AF_median = LineCollection(lines_min_AF_median, colors=y_color, linewidth=1.2)
            lc_maj_AF_median = LineCollection(lines_maj_AF_median, colors=y_color, linewidth=1.2)
            ax2.add_collection(lc_prob)
            ax2.add_collection(lc_min_AF_median)
            ax2.add_collection(lc_maj_AF_median)
            ax2.set_xlim(contig_beginning_end[0][0], contig_beginning_end[-1][1])
            ax2.set_ylim([-0.02, 1.02])
            ax2.set_ylabel('Allele fraction')
            ax2.plot([], color=(1.0, 0.0, 0.0), label='Prob(normal)')
            ax2.plot([], color=(0.0, 1.0, 0.0), label='Normal')
            ax2.plot([], color=(0.0, 0.0, 1.0), label='Deletion')
            ax2.plot([], color=(1.0, 0.0, 1.0), label='Amplification')
            ax2.plot([], color=(0.4, 0.4, 0.4), label='CNLOH')
            ax2.legend(loc="upper right", bbox_to_anchor=(1,1))

            plt.savefig(self.fig_del_ampl_filename)
            plt.close(fig2)

    def __plot_clustering(self):
        """ This function plots the details of how the fitting was done to the normal cluster
            in copy ratio and allele fraction space.
        """
        fig = plt.figure(1, figsize=(12, 6), dpi=400)
        plt.subplot(221)
        xedges = np.linspace(0, 5, 50)
        yedges = np.linspace(0, 0.5, 20)
        plt.hist2d(self.__copy_ratio_sampled, self.__allele_fraction_sampled,
                   bins=(xedges, yedges), norm=LogNorm(), cmap="gray_r")
        plt.colorbar()
        plt.xlabel("Copy ratio samples")
        plt.ylabel("Allele fraction samples")

        plt.subplot(222)
        self.__plot_classification()
        plt.xlabel("Copy ratio samples")
        plt.ylabel("Allele fraction samples")
        plt.xlim((0, 5))
        plt.ylim((0, 0.5))
        plt.xticks(np.arange(0, 6, 1))
        plt.yticks(np.arange(0, 0.6, 0.1))

        plt.subplot(223)
        self.__plot_Gaussian_mixture_fit()
        plt.xlabel("Copy ratio samples")
        plt.ylabel("Allele fraction samples")
        plt.xlim((0, 5))
        plt.ylim((0, 0.5))
        plt.xticks(np.arange(0, 6, 1))
        plt.yticks(np.arange(0, 0.6, 0.1))

        plt.subplot(224)
        plt.hist(self.__copy_ratio_sampled, bins = "auto", color = "r")
        plt.xlabel("Copy ratio")
        plt.xlim((0, 5))
        plt.xticks(np.arange(0, 6, 1))

        plt.savefig(self.fig_scatter_plot)
        plt.close(fig)

    def __plot_classification(self):
        """ Plot the medians of the distributions in copy ratio and allele fraction space
            and indicate the normal and not normal segments with different colors.
        """
        samples = np.asarray([self.__copy_ratio_sampled, self.__allele_fraction_sampled]).T
        classification = self.__gaussian_mixture_fit.predict(samples)
        colors = ["k" if i in self.__normal_peak_indices else "r" for i in classification]
        plt.scatter(samples[:,0], samples[:,1], c=colors, alpha=0.8, s=10)

    def __plot_Gaussian_mixture_fit(self):
        """ Plot the fit of Gaussian clusters to the 2D copy ratio and allele fraction data.
        """
        samples = np.asarray([self.__copy_ratio_sampled, self.__allele_fraction_sampled]).T
        n_Gaussians = len(self.__gaussian_mixture_fit.weights_)
        X, Y = np.mgrid[0:5:.05, 0:0.5:0.01]
        pos = np.empty(X.shape + (2,))
        pos[:, :, 0] = X
        pos[:, :, 1] = Y
        Z = 0 * X
        for i in range(n_Gaussians):
            rv = multivariate_normal(self.__gaussian_mixture_fit.means_[i],
                                     self.__gaussian_mixture_fit.covariances_[i])
            Z += self.__gaussian_mixture_fit.weights_[i] * rv.pdf(pos)
        Z = np.asarray([[(Z[i][j] + 0.0000000000001) for i in range(len(Z))] for j in range(len(Z[0]))]).T
        plt.contourf(X, Y, Z, cmap="gray_r", locator=ticker.LogLocator())
        plt.colorbar()

    def __gray_background_contigs(self, contig_beginning_end: List, ymin: float, ymax: float, ax):
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
        """ Save the results.
        """
        input_filename = self.__CR_AF_data.get_input_filename()
        [avg_normal_CR, std_dev_normal_CR] = self.__average_and_std_dev_copy_ratio_normal_segments()
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
                if not math.isnan(float(values[9])):
                    file_data += line.strip() + "\t" + str(self.__normal_del_ampl(self.__copy_ratio_median[i],
                                                                                  avg_normal_CR, std_dev_normal_CR,
                                                                                  self.__responsibilities_normal[i]))
                    file_data += "\t"
                    file_data += str(self.__get_phred_score(probability=self.__responsibilities_normal[i],
                                                            max_phred_score=self.__max_PHRED_score))
                    file_data += "\n"
                    i += 1
        file_header = file_header[0:-1] + "\tCALL\tPHRED\n"

        output_file = open(self.output_calls_filename, "w")
        output_file.write(file_header + file_data)
        output_file.close()

    def __normal_del_ampl(self, CR: float, avg_normal_CR: float, std_dev_normal_CR: float,
                          responsibility_normal: float, responsibility_threshold: float=0.5):
        """ Determine if the copy ratio value CR corresponds to a deletion (-),
            an amplification (+), a normal segment (0) or a normal copy ratio
            with imbalanced allele fraction (/).
        """
        if responsibility_normal > responsibility_threshold:
            return "0"
        else:
            if CR > avg_normal_CR + 2 * std_dev_normal_CR:
                return "+"
            elif CR < avg_normal_CR - 2 * std_dev_normal_CR:
                return "-"
            else:
                return "/"

    def __average_and_std_dev_copy_ratio_normal_segments(self, responsibility_threshold: float=0.5):
        """ Auxiliary function to determine the averge and the standard deviation of the copy ratio
            value of all the normal segments.
        """
        total_length = 0
        total_CR = 0
        total_var_CR = 0
        for i in range(len(self.__copy_ratio_median)):
            if self.__responsibilities_normal[i] > responsibility_threshold:
                total_length += self.__n_points_CR[i]
                total_CR += self.__n_points_CR[i] * self.__copy_ratio_median[i]

        if total_length == 0:
            avg_normal_CR = 0
        else:
            avg_normal_CR = total_CR / total_length

        for i in range(len(self.__copy_ratio_median)):
            if self.__responsibilities_normal[i] > responsibility_threshold:
                total_var_CR += self.__n_points_CR[i] * (self.__copy_ratio_median[i] - avg_normal_CR)**2

        if total_length == 0:
            var_normal_CR = 0
        else:
            var_normal_CR = total_var_CR / total_length

        std_dev_normal_CR = var_normal_CR**0.5

        return [avg_normal_CR, std_dev_normal_CR]
