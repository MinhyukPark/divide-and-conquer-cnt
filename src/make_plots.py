import glob
import os

import click
import dendropy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

import CN3

@click.command()
@click.option("--error-cnp-dir", required=True, type=click.Path(exists=True), help="The input error cnp directory with the logs")
@click.option("--output-prefix", required=True, type=click.Path(), help="Output prefix")
def make_plots(error_cnp_dir, output_prefix):
    plt.rc("font", size=20)
    cumulative_rf_rate = 0.0
    k_arr = [4, 6, 8]
    n_arr = [15, 20, 30, 40]
    label_arr = []
    time_dict = {}
    rf_dict = {}
    for k in k_arr:
        for n in n_arr:
            label_arr.append("k" + str(k) + "_" + "n" + str(n))
    for label in label_arr:
        current_error_log = error_cnp_dir + "/" + label + "/main.err"
        current_count = 0
        with open(current_error_log, "r") as f:
            for line in f:
                if("Elapsed (wall clock)" in line):
                    time_dict[label] = (int(line.split()[7].split(":")[0]) * 60) / current_count
                if("rf rate so far" in line):
                    if(label not in rf_dict):
                        rf_dict[label] = []
                    rf_dict[label].append(float(line.split()[7]))
                    current_count += 1
    for label in rf_dict:
        print(label + ": " + str(rf_dict[label]) + " mean is " + str(np.mean(rf_dict[label])))
        print(label + ": " + str(time_dict[label]))
    print(time_dict)
    rf_arr = []
    time_arr = []
    for label in label_arr:
        rf_arr.append(np.array(rf_dict[label]))
        time_arr.append(time_dict[label])

    color_arr = ["#4169E1", "#7EC850", "#722F37", "#B19CD9"] * 3
    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
    parts = ax.violinplot(rf_arr, points=20, widths=0.3, showmeans=False, showextrema=False, showmedians=False)
    # ax.set_xticks(range(len(label_arr)))
    # ax.set_xticklabels(label_arr)
    ax.set_xticks([2.5, 6.5, 10.5])
    ax.set_xticklabels(["4", "6", "8"])
    for body_index,body in enumerate(parts["bodies"]):
        body.set_facecolor(color_arr[body_index])
        body.set_edgecolor("black")
        body.set_alpha(1)
    q1 = {}
    medians = {}
    q3 = {}
    for label in label_arr:
        current_q1,current_medians,current_q3 = np.percentile(rf_dict[label], [25, 50, 75])
        q1[label] = current_q1
        medians[label] = current_medians
        q3[label] = current_q3
    print(q1)
    print(medians)
    print(q3)
    whiskers = []
    for label in label_arr:
        whiskers.append(np.array(adjacent_values(rf_dict[label], q1[label], q3[label])))
    whiskers = np.stack(whiskers)
    print(whiskers)
    whiskers_min,whiskers_max = whiskers[:,0], whiskers[:,1]

    inds = np.arange(1, len(rf_arr) + 1)
    ax.scatter(inds, medians.values(), marker="o", color="white", s=30, zorder=3)
    ax.vlines(inds, q1.values(), q3.values(), color="k", linestyle="-", lw=5)
    ax.vlines(inds, whiskers_min, whiskers_max, color="k", linestyle="-", lw=1)
    ax.set_ylim((0.0, 1.0))
    ax.yaxis.grid("True", color="#ffffff")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_facecolor("#dddddd")
    ax.set_xlabel("k")
    ax.set_ylabel("RF rate")

    handles = []
    for color_index,color in enumerate(color_arr[:4]):
        handles.append(mpatches.Patch(color=color, label=str(n_arr[color_index])))
    ax.legend(handles=handles, prop={"size": 13})


    fig.savefig(output_prefix + "/fn_rate.png")
    plt.figure()


    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
    parts = ax.bar(label_arr, time_arr, color=color_arr)
    ax.set_xlabel("k")
    ax.set_ylabel("Time (s)")
    ax.set_xticks([1.5, 5.5, 9.5])
    ax.set_xticklabels(["4", "6", "8"])
    handles = []
    for color_index,color in enumerate(color_arr[:4]):
        handles.append(mpatches.Patch(color=color, label=str(n_arr[color_index])))
    ax.legend(handles=handles, prop={"size": 13})
    fig.savefig(output_prefix + "/time.png")



def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


if __name__ == "__main__":
    make_plots()

