"""
Python script for meta-d' analysis based on SDT framework

Given data from an experiment where an observer discriminates between two stimulus alternatives
on every trial and provides confidence ratings, converts trial by trial experimental information for N trials
into response counts, and compute meta-d' using either MLE or SSE method.

see http://www.columbia.edu/~bsm2105/type2sdt/ for explanation.

== UPDATE HISTORY ==
2018-01-01 Created
"""
import numpy as np
import pandas as pd
from scipy.stats import norm


class MetaD:
    """The overall MetaD class instance on which the calculations are performed.
    """
    def __init__(self, stim_id, response, rating, num_ratings, pad_cells=0, pad_amount=None):
        """== To create an instance of MetaD, the following INPUTS must be supplied ==
        stim_id:     one-dimensional list or numpy array  stim_id[i - 1] = 0 --> stimulus on i'th trial was S1
                                                          stim_id[i - 1] = 1 --> stimulus on i'th trial was S2

        response:    one-dimensional list or numpy array  response[i - 1] = 0 --> response on i'th trial was 'S1'
                                                          response[i - 1] = 1 --> response on i'th trial was 'S2'

        rating:      one-dimensional list or numpy array  rating[i - 1] = X --> rating on i'th trial was X
                                                          X must be in the range 1 <= X <= num_ratings

        N.B. all trials where stim_id is not 0 or 1, response is not 0 or 1, or rating is
        not in the range [1, num_ratings] are omitted from the response count. Length of inputs must agree.

        num_ratings: total number of available subjective ratings for the subject.
                     e.g. if the subject can rate confidenceon a scale of 1 - 4, then num_ratings = 4

        == OPTIONAL INPUTS ==
        pad_cells:   if set to 1, each response count in the MetaD has the value of pad_amount added to it.
                     Padding cells is desirable if trial counts of 0 interfere with model fitting.
                     If set to 0, trial counts are not manipulated and 0s may be present in the data
                     contained in the MetaD instance.

        padAmount:   the value to add to each response count oif padCess is set to 1.
                     default value is 1/(2 * num_ratings)

        == CRITICAL ATTRIBUTES ==
        MetaD.nr_s1, MetaD.nr_s2
        These are numpy arrays containing the total number of responses in each response category,
        conditional on presentation of S1 and S2.

        e.g. if nr_s1 = [100, 50, 20, 10, 5, 1], then when stimulus S1 was
        presented, the subject had the following response counts:
        responded S1, rating=3 : 100 times
        responded S1, rating=2 : 50 times
        responded S1, rating=1 : 20 times
        responded S2, rating=1 : 10 times
        responded S2, rating=2 : 5 times
        responded S2, rating=3 : 1 time

        The ordering of response / rating counts for S2 should be the same as it
        is for S1. e.g. if nr_s2 = [3, 7, 8, 12, 27, 89], then when stimulus S2 was
        presented, the subject had the following response counts:
        responded S1, rating=3 : 3 times
        responded S1, rating=2 : 7 times
        responded S1, rating=1 : 8 times
        responded S2, rating=1 : 12 times
        responded S2, rating=2 : 27 times
        responded S2, rating=3 : 89 times

        Many other attributes are set to None when the instance is created,
        but will be filled up with numbers once certain estimation is performed.
        """

        # Set pad amount if padding cells needed
        if pad_cells == 1:
            if pad_amount is None:
                pad_amount = 1.0 / (2 * num_ratings)
        elif pad_cells == 0:
            pad_amount = 0

        # Check input type and if their dimensions agree
        if not isinstance(stim_id, np.ndarray):
            stim_id = np.array(stim_id, dtype=int)
        if not isinstance(response, np.ndarray):
            response = np.array(response, dtype=int)
        if not isinstance(rating, np.ndarray):
            rating = np.array(rating, dtype=int)

        if (len(stim_id) == len(response)) and (len(stim_id) == len(rating)):
            pass
        else:
            print "The dimension of your inputs do not agree!"
            quit()

        # Create a pandas DataFrame
        dic = {'stim_id': stim_id, 'response': response, 'rating': rating}
        d = pd.DataFrame(dic)

        # Filter bad trials
        d = d[(d['stim_id'] == 0) | (d['stim_id'] == 1)]
        d = d[(d['response'] == 0) | (d['response'] == 1)]
        d = d[(d['rating'] >= 1) & (d['rating'] <= num_ratings)]


        # Process data
        nr_s1 = np.zeros([num_ratings * 2, ], dtype=int)
        nr_s2 = np.zeros([num_ratings * 2, ], dtype=int)
        s1_only = d[d['response'] == 0]
        s2_only = d[d['response'] == 1]


        for j in range(num_ratings, 0, -1):
            nr_s1[num_ratings - j] = s1_only[(s1_only['stim_id'] == 0) & (s1_only['rating'] == j)].shape[0]
            nr_s2[num_ratings - j] = s1_only[(s1_only['stim_id'] == 1) & (s1_only['rating'] == j)].shape[0]


        for j in range(1, num_ratings + 1):
            nr_s1[num_ratings - 1 + j] = s2_only[(s2_only['stim_id'] == 0) & (s2_only['rating'] == j)].shape[0]
            nr_s2[num_ratings - 1 + j] = s2_only[(s2_only['stim_id'] == 1) & (s2_only['rating'] == j)].shape[0]


        if pad_cells:
            nr_s1 = nr_s1 + pad_amount
            nr_s2 = nr_s2 + pad_amount

        self.stim_id = stim_id
        self.response = response
        self.rating = rating
        self.data = d
        self.num_ratings = num_ratings
        self.pad_cells = pad_cells
        self.pad_amount = pad_amount
        self.nr_s1 = nr_s1
        self.nr_s2 = nr_s2
        self.type2_fit_sse = None
        self.d_a = None
        self.meta_d_a = None
        self.m_ratio = None
        self.m_diff = None
        self.s = None
        self.type2_fit_mle = None
        
    def type2_sdt_sse(self, d_min=-5, d_max=5, d_grain=0.01):
        """This method estimates meta-d' by minimizing the sum of squared errors (SSE)
        By calling it, the following attributes are added to the MetaD class instance
        or existing attributes are modified

        MetaD.d_a        : d_a for input data. If s = 1, d_a = d'
        MetaD.meta_d_a   : meta_d for input data
        MetaD.m_ratio    : meta_d_a / d_a; measure of metacognitive efficiency
        MetaD.m_diff     : meta_d_a - d_a; measure of metacognitive efficiency
        MetaD.s          : ratio of evidence distribution standard deviations assumed
                                 for the analysis
        MetaD.fit        : a dictionary that contains other information, see
                                 function fit_meta_d_sse for details.

        Given data from an experiment where an observer discriminates between two stimulus
        alternatives on every trial and provides confidence ratings, provides a type 2 SDT
        analysis of the data.

        The method does a standard type 1 SDT analysis on the raw behavioral data and then
        does a type 2 SDT analysis using the function fit_meta_d_sse with default values:
        d_min = -5
        d_grain = 0.01
        d_max = 5

        == UPDATE HISTORY ==
        2018-01-01 Created
        """

        # Check input validity, input arrays
        # (nr_s1 and nr_s2 must have an even number of elements)
        if (len(self.nr_s1) == len(self.nr_s2)) and len(self.nr_s1) % 2 == 0:
            pass
        else:
            print "nr_s1 and nr_s2 must be the same length and have an even number of elements.\n"
            print "Consider adding pad amounts."
            quit()

        num_ratings = self.num_ratings

        # Standard SDT analysis
        hr1 = sum(self.nr_s2[num_ratings:]) / float(sum(self.nr_s2))
        far1 = sum(self.nr_s1[num_ratings:]) / float(sum(self.nr_s1))


        rating_hrs = np.zeros([2 * num_ratings - 1, ], dtype=float)
        rating_fars = np.zeros([2 * num_ratings - 1, ], dtype=float)
        for i in range(2 * num_ratings - 2):
            rating_hrs[i] = sum(self.nr_s2[i + 1:]) / float(sum(self.nr_s2))
            rating_fars[i] = sum(self.nr_s1[i + 1:]) / float(sum(self.nr_s1))

        s = 1

        # d' and c in terms of S1 distribution standard deviation units

        d_1 = (1.0 / s) * norm.ppf(hr1) - norm.ppf(far1)
        c_1 = (-1.0 / (1 + s)) * (norm.ppf(hr1) + norm.ppf(far1))
        cprime = c_1 / d_1

        # Type 2 SDT analysis
        # Get type 2 HR and FAR for S1 responses
        hr2_rs1 = np.zeros([num_ratings - 1, ], dtype=float)
        far2_rs1 = np.zeros([num_ratings - 1, ], dtype=float)
        for i in range(num_ratings - 1):
            hr2_rs1[i] = sum(self.nr_s1[:i + 1]) / float(sum(self.nr_s1[:num_ratings]))
            far2_rs1[i] = sum(self.nr_s2[:i + 1]) / float(sum(self.nr_s2[:num_ratings]))

        # Get type 2 HR and FAR for S2 responses
        hr2_rs2 = np.zeros([num_ratings - 1, ], dtype=float)
        far2_rs2 = np.zeros([num_ratings - 1, ], dtype=float)
        for i in range(num_ratings + 1, 2 * num_ratings):
            hr2_rs2[i - (num_ratings + 1)] = sum(self.nr_s2[i:]) / float(sum(self.nr_s2[num_ratings:]))
            far2_rs2[i - (num_ratings + 1)] = sum(self.nr_s1[i:]) / float(sum(self.nr_s1[num_ratings:]))

        def fit_meta_d_sse(obs_hr2_rs1, obs_far2_rs1, obs_hr2_rs2, obs_far2_rs2, c_prime, s_ratio, dmin, dmax, dgrain):
            """Given response-conditional type 2 hit rates and type 2 false alarm rates, as well
            as the empirically estimated relative criterion cprime = c/d', use a signal detection
            theory model to estimate meta-d', the value of d' that would have been expected to
            generate the observed type 2 data.

            Estimation is done by testing different values for meta-d' to see what value gives the
            best fits to the observed data (best fit = minimize the sum squared error between observed
            and expected type 2 HR and type 2 FAR).

            == SOME INPUTS EXPLAINED ==
            s           : the ratio of the standard deviations of the evidence distribution for
                          stimulus classes S1 and S2. It can be estimated from rating data. If
                          unspecified, it will be set to 1.
            d_min       : the minimum value for meta-d' that will be tested, default to -5.
            d_max       : the maximum value for meta_d' that will be tested, default to 5.
            d_grain     : the step size used in testing values for meta-d', default to .01.

            == OUTPUT ==
            The output will be packaged in a dictionary, with key-value pairs as follows:
            ["meta_d"]       : meta-d' value that minimizes the SSE between observed and expected
                               type 2 data. If s is not equal to 1, meta-d' is specified in units
                               of the S1 distribution
            ["meta_c"]       : the value of type 1 criterion c used in conjunction with meta-d'.
                               meta_c / meta_d = cprime, the constant type 1 criterion specified
                               in the input. If s is not equal to 1, meta-c is specified in units
                               of the S1 distribution
            ["s_ratio"]      : the value of s used in the type 2 data fitting, where s = sd(S1) / sd(S2)
            ["t2c_rs1"]      : values for the type 2 criteria that, along with meta-d' and c', provide
                               the best fit for type 2 data for "S1" response
            ["t2c_rs2"]      : likewise, for "S2" responses
            ["sse"]          : sum of squared errors between observed and expected type 2 data
            ["est_hr2_rs1"]  : the type 2 hit rates for "S1" responses expected from meta_d, meta_c, s,
                               and t2c_rs1
            ["obs_hr2_rs1"]  : empirically observed type 2 hit rates for "S1" responses
            ["est_far2_rs1"]
            ["obs_far2_rs1"]
            ["est_hr2_rs2"]
            ...              : likewise as above
            """
            # Initialize analysis
            nratings = len(obs_hr2_rs1)
            ds = np.linspace(dmin, dmax, (dmax - dmin) / dgrain + 1)
            sse_min = float("inf")
            meta_d = 0
            meta_c = 0
            t2c_rs1 = np.array([], dtype=float).reshape([-1, ])
            t2c_rs2 = np.array([], dtype=float).reshape([-1, ])
            est_hr2_rs1 = np.array([], dtype=float).reshape([-1, ])
            est_far2_rs1 = np.array([], dtype=float).reshape([-1, ])
            est_hr2_rs2 = np.array([], dtype=float).reshape([-1, ])
            est_far2_rs2 = np.array([], dtype=float).reshape([-1, ])

            # Search for meta-d' that minimizes the type 2 SSE
            for k in range(len(ds)):

                # Initialize parameters for current level of meta-d'
                d = ds[k]
                c = c_prime * d
                s1mu = -d / 2
                s2mu = d / 2
                s1sd = 1
                s2sd = 1 / s_ratio

                x = np.linspace(s1mu - 5 * max([s1sd, s2sd]), s2mu + 5 * max([s1sd, s2sd]),
                                ((s1mu + 5 * max([s1sd, s2sd])) - (s2mu - 5 * max([s1sd, s2sd]))) / 0.001 + 1)
                c_ind = np.where(abs(x - c) == abs(x - c).min())[0][0]

                hrs = 1 - norm.cdf(x, loc=s2mu, scale=s2sd)
                fars = 1 - norm.cdf(x, loc=s1mu, scale=s1sd)

                # Fit type 2 data for S1 responses
                est_hr2s_rs1 = (1 - fars[:c_ind + 1]) / (1 - fars[c_ind])
                est_far2s_rs1 = (1 - hrs[:c_ind + 1]) / (1 - hrs[c_ind])

                sse_rs1 = np.array([], dtype=float).reshape([-1, ])
                rs1_ind = np.array([], dtype=int).reshape([-1, ])
                for n in range(nratings):
                    sse = (est_hr2s_rs1 - obs_hr2_rs1[n]) ** 2 + (est_far2s_rs1 - obs_far2_rs1[n]) ** 2
                    sse_rs1 = np.append(sse_rs1, sse.min())
                    rs1_ind = np.append(rs1_ind, np.where(sse == sse.min())[0][0])

                # Fit type 2 data for S2 responses
                est_hr2s_rs2 = hrs[c_ind:] / hrs[c_ind]
                est_far2s_rs2 = fars[c_ind:] / fars[c_ind]

                sse_rs2 = np.array([], dtype=float).reshape([-1, ])
                rs2_ind = np.array([], dtype=int).reshape([-1, ])
                for n in range(nratings):
                    sse = (est_hr2s_rs2 - obs_hr2_rs2[n]) ** 2 + (est_far2s_rs2 - obs_far2_rs2[n]) ** 2
                    sse_rs2 = np.append(sse_rs2, sse.min())
                    rs2_ind = np.append(rs2_ind, np.where(sse == sse.min())[0][0])

                # Update analysis
                sse_tot = sse_rs1.sum() + sse_rs2.sum()
                if sse_tot < sse_min:
                    sse_min = sse_tot
                    meta_d = d
                    meta_c = c
                    t2c_rs1 = x[rs1_ind]
                    t2c_rs2 = x[rs2_ind + c_ind - 1]
                    est_hr2_rs1 = est_hr2s_rs1[rs1_ind]
                    est_far2_rs1 = est_far2s_rs1[rs1_ind]
                    est_hr2_rs2 = est_hr2s_rs2[rs2_ind]
                    est_far2_rs2 = est_far2s_rs2[rs2_ind]

            # Package output
            fit = {
                "meta_d": meta_d,
                "meta_c": meta_c,
                "s_ratio": s_ratio,
                "t2c_rs1": t2c_rs1,
                "t2c_rs2": t2c_rs2,
                "sse": sse_min,
                "est_hr2_rs1": est_hr2_rs1,
                "obs_hr2_rs1": obs_hr2_rs1,
                "est_far2_rs1": est_far2_rs1,
                "obs_far2_rs1": obs_far2_rs1,
                "est_hr2_rs2": est_hr2_rs2,
                "obs_hr2_rs2": obs_hr2_rs2,
                "est_far2_rs2": est_far2_rs2,
                "obs_far2_rs2": obs_far2_rs2
            }

            return fit

        self.type2_fit_sse = fit_meta_d_sse(hr2_rs1, far2_rs1, hr2_rs2, far2_rs2, cprime, s, d_min, d_max, d_grain)
        self.d_a = d_1 * s * (2 / float(1 + s ** 2)) ** 0.5
        self.meta_d_a = self.type2_fit_sse['meta_d'] * s * (2 / float(1 + s ** 2)) ** 0.5
        self.m_ratio = self.meta_d_a / self.d_a
        self.m_diff = self.meta_d_a - self.d_a
        self.s = s

    # MLE method is under development
    def fit_meta_d_mle(self, s=1, fncdf=norm.cdf, fninv=norm.ppf):
        """Given data from an experiment where observer discriminates between two stimulus
        alternatives on every trial and provides confidence ratings,  this method
        provides a type 2 SDT analysis of the data.

        == INPUTS ==
        self.nr_s1 and self.nr_s2
        refer to the class MetaD

        *s
        The ratio of standard deviation for type 1 distributions, i.e.
        s = sd(S1) / sd(S2)
        if not specified, s is set to a default value of 1

        *fncdf
        A function handle for the CDF of the type 1 distribution.
        If not specified, default to scipy.stats.norm.cdf(loc=0, scale=1)

        *fninv
        A function handle for the inverse cdf of the type 1 distribution
        If not specified, default to scipy.stats.norm.ppf(loc=0, scale=1)

        == OUTPUTS ==
        The output is packaged in a dictionary called "type2_fit_mle", added as an
        attribute to MetaD instance, with the following key-value pairs:
        ["d_a"]          : mean(S2) - mean(S1), in room-mean-square (sd(S1), sd(S2)) units
        ["s"]            : sd(S1) / sd(S2)
        ["meta_d_a"]     : meta-d' in RMS units
        ["m_diff"]       : meta_d_a - d_a
        ["m_ratio"]      : meta_d_a / d_a
        ["meta_c_a"]     : type 1 criterion for meta-d' fit, RMS units
        ["t2ca_rs1"]     : type 2 criteria of "S1" responses for meta-d' fit, RMS units
        ["t2ca_rs2"]     : type 2 criteria of "S2" responses for meta-d' fit, RMS units
        ["s1units"]      : contains same parameters in sd(S1) units.
                           These may be of use since the data-fitting is conducted using
                           parameters specified in sd(S1) units.
        ["logl"]         : log likelihood of the data fit
        ["est_HR2_rS1"]  : estimated (from meta-d' fit) type 2 hit rates for S1 responses
        ["obs_HR2_rS1"]  : actual type 2 hit rates for S1 responses
        ["est_FAR2_rS1"] : estimated type 2 false alarm rates for S1 responses
        ["obs_FAR2_rS1"] : actual type 2 false alarm rates for S1 responses
        ["est_HR2_rS2"]  : estimated type 2 hit rates for S2 responses
        ["obs_HR2_rS2"]  : actual type 2 hit rates for S2 responses
        ["est_FAR2_rS2"] : estimated type 2 false alarm rates for S2 responses
        ["obs_FAR2_rS2"] : actual type 2 false alarm rates for S2 responses

        If there are N ratings, then there will be N - 1 type 2 hit rates and
        false alarm rates.
        """
