import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp


def two_compartment_model(t, k10, k12, k21, C1_0):
    """
    # Calculates concentration at a time point
    #   using the two compartment model equation.

    Parameters:
        time_points (array): Time points.
        k10: Elimination rate constant (plasma → elimination).
        k12: Rate constant for drug moving plasma → tissue.
        k21: Rate constant for drug moving tissue → plasma.
        C1_0 (float): Initial concentration in the compartment.

    Returns:
        optimized concentration
    """
    def odes(t, y):
        C1, C2 = y
        dC1_dt = -k10 * C1 - k12 * C1 + k21 * C2
        dC2_dt = k12 * C1 - k21 * C2
        return [dC1_dt, dC2_dt]

    #t = np.unique(t)  # Ensure time points are unique and sorted
    y0 = [C1_0, 0]  # Initial concentrations
    solution = solve_ivp(odes, [t[0], t[-1]], y0, t_eval=t)
    return solution.y[0]  # Return plasma concentration (C1)


def fit_replicate(time_points, concentrations, initial_guess, C1_0):
    """
    # Calculates goodness-of-fit metrics for the model.

    Parameters:
        time_points (array): Time points.
        cocentrations (array): Observed concentrations.
        initial_guess (list): Initial guess for the parameter [k].
        C1_0 (float): Initial concentration in the compartment.

    Returns:
        Optimized [k] parameters (list).
    """
    def fit_model(t, k10, k12, k21):
        return two_compartment_model(t, k10, k12, k21, C1_0)

    try:
        params_opt, _ = curve_fit(
            fit_model, time_points, concentrations,
            p0=initial_guess, bounds=([0, 0, 0], [np.inf, np.inf, np.inf])
        )
        return params_opt
    except Exception as e:
        print(f"Error fitting replicate: {e}")
        return [np.nan, np.nan, np.nan]


def calculate_goodness_of_fit(time_points, observed, predicted):
    """
    # Calculates goodness-of-fit metrics for the model.

    Parameters:
        time_points (array): Time points.
        observed (array): Observed concentrations.
        predicted (array): Model-predicted concentrations.

    Returns:
        tuple: Chi-squared statistic, RMSE, and K-S test statistic (ks_stat).
    """
    chi_squared = np.sum((observed - predicted) ** 2)  # Uniform variance assumed
    rmse = np.sqrt(np.mean((observed - predicted) ** 2))
    ks_stat, _ = ks_2samp(observed, predicted)
    return chi_squared, rmse, ks_stat


def model(file_path, data_raw, C1_0, initial_guess):
    """
    # Processes experimental data, Fits a two-compartment model to all replicates
    #   calculates goodness of fitfits the single-compartment model, and summarizes results.
    Parameters:
        file_path (str): Path to the input file.
        data_raw (DataFrame): Raw experimental data.
        C1_0 (float): Initial concentration in the compartment.
        initial_guess (list): Initial guess for the parameter [k].

    Returns:
        DataFrame: Summary of optimized parameters and goodness-of-fit metrics.
    """
    if 'Time' not in data_raw.columns:
        raise ValueError("The input file must contain a 'Time' column.")

    data = data_raw.dropna()
    if data.empty:
        raise ValueError("The input file contains no valid data after removing NaN values.")

    time_points = data['Time'].values #pulls data from dataframe
    replicates = data.drop(columns=['Time']).values

    optimized_parameters = []
    goodness_of_fit = []

    for i in range(replicates.shape[1]):
        replicate_data = replicates[:, i]
        params_opt = fit_replicate(time_points, replicate_data, initial_guess, C1_0)
        optimized_parameters.append(params_opt)

        if not np.isnan(params_opt).any():
            predicted_curve = two_compartment_model(time_points, *params_opt, C1_0)
            chi_squared, rmse, ks_stat = calculate_goodness_of_fit(time_points, replicate_data, predicted_curve)
            goodness_of_fit.append((chi_squared, rmse, ks_stat))
        else:
            goodness_of_fit.append((np.nan, np.nan, np.nan))

    param_df = pd.DataFrame(optimized_parameters, columns=['k10', 'k12', 'k21'])
    goodness_df = pd.DataFrame(goodness_of_fit, columns=['Chi-Squared', 'RMSE', "KS_STAT"])

    #Prints the results of the model
    results_df = pd.concat([param_df, goodness_df], axis=1)
    print(f"\nOptimized Parameters for '{file_path}' Across Replicates:")
    print(results_df)
    print("\nSummary Statistics:\n", results_df.describe())

    return results_df


def main():
    args = sys.argv[1:]
    for filename in args:
        try:
            data_raw = pd.read_excel(filename)
            print(f"Processing file: {filename}")
            model(filename, data_raw, C1_0=100, initial_guess=[0.1, 0.05, 0.05])
        except FileNotFoundError:
            print(f"Error: The file '{filename}' was not found.")
        except Exception as e:
            print(f"An error occurred while processing '{filename}': {e}")


if __name__ == '__main__':
    main()
