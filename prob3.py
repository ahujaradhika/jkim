import numpy as np

# trial parameters
size_of_distribution = 5
mu_target, sigma_target = 0.5, 1 # mean and standard deviation of target distribution
mu_non_target, sigma_non_target = 0.25, 1 # mean and standard deviation of non-target distribution
num_of_trials = 5


def run_trial(mu_target, sigma_target, mu_non_target, sigma_non_target, num_of_trials, size):

    # Create distributions
    target_dist = np.random.normal(mu_target, sigma_target, size_of_distribution)
    non_target_dist = np.random.normal(mu_non_target, sigma_non_target, size_of_distribution)

    # c is the criteria: midpoint of the two means
    c = (mu_target - mu_non_target)/2

    yes_count = 0
    no_count = 0

    for i in range(num_of_trials):
        chance = np.random.randint(0, 2)
        if chance == 0:
            val = np.random.choice(target_dist, 1, replace=False)
            if val > c:
                yes_count += 1
            else:
                no_count += 1
        else:
            val = np.random.choice(non_target_dist, 1, replace=False)
            if val > c:
                yes_count += 1
            else:
                no_count += 1
    
    print(yes_count)
    print(no_count)

    return (yes_count/no_count)

print(run_trial(mu_target, sigma_target, mu_non_target, sigma_non_target, num_of_trials, size_of_distribution))
