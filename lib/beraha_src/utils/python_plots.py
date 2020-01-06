import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.lines import Line2D
from scipy.stats import norm, skewnorm, t

sns.set_context(rc={"lines.linewidth": 2.5})


def plot_yppd(y_ppd, y_obs=None, outfile=""):
    current_palette = sns.color_palette()

    numGroups = y_ppd.shape[1]
    fig, axes = plt.subplots(nrows=1, ncols=numGroups, figsize=(16, 6))
    for i in range(numGroups):
        sns.kdeplot(y_ppd[:, i], ax=axes[i])

    if y_obs is not None:
        for i in range(numGroups):
            sns.kdeplot(y_obs[i, :], ax=axes[i], linestyle="--")

    custom_lines = [
        Line2D([0], [0], color=current_palette[i], lw=4)
        for i in range(2)]
    fig.legend(
        custom_lines, ["Estimated", "Observed"], loc=8,
        ncol=2,  fontsize=14, frameon=False)

    for i in range(numGroups):
        axes[i].set_title("Group {0}".format(i + 1))

    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()


def plot_exp_1(y_ppd, y_obs=None, title="", outfile=""):
    # Group 1
    current_palette = sns.color_palette()

    x = np.linspace(-3, 3, num=100)
    fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(16, 6))
    if title:
        fig.suptitle(title, fontsize=18)
    sns.lineplot(x, norm.pdf(x), ax=axes[0])
    sns.kdeplot(y_ppd[:, 0], ax=axes[0])

    sns.lineplot(x, norm.pdf(x), ax=axes[1])
    sns.kdeplot(y_ppd[:, 1], ax=axes[1])

    sns.lineplot(x, norm.pdf(x), ax=axes[2])
    sns.kdeplot(y_ppd[:, 2], ax=axes[2])

    sns.lineplot(x, skewnorm.pdf(x, 1), ax=axes[3])
    sns.kdeplot(y_ppd[:, 3], ax=axes[3])

    for i in range(4):
        if y_obs is not None:
            sns.kdeplot(y_obs[i, :], ax=axes[i], linestyle="--")

    custom_lines = [
        Line2D([0], [0], color=current_palette[i], lw=4)
        for i in range(3)]
    fig.legend(
        custom_lines, ["True", "Estimated", "Observed"], loc=8,
        ncol=3,  fontsize=14, frameon=False)

    for i in range(4):
        axes[i].set_title("Group {0}".format(i + 1))
        axes[i].set_xlim((-5, 5))

    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()


def plot_exp_2(y_ppd, y_obs=None, title="", outfile=""):
    current_palette = sns.color_palette()

    x = np.linspace(-5, 5, num=100)
    fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(16, 6))
    if title:
        fig.suptitle(title, fontsize=18)
    sns.lineplot(x, norm.pdf(x), ax=axes[0])
    sns.kdeplot(y_ppd[:, 0], ax=axes[0])

    sns.lineplot(x, norm.pdf(x, scale=2.25), ax=axes[1])
    sns.kdeplot(y_ppd[:, 1], ax=axes[1])

    sns.lineplot(x, norm.pdf(x, scale=0.25), ax=axes[2])
    sns.kdeplot(y_ppd[:, 2], ax=axes[2])

    sns.lineplot(x, norm.pdf(x, scale=1), ax=axes[3])
    sns.kdeplot(y_ppd[:, 3], ax=axes[3])

    for i in range(4):
        if y_obs is not None:
            sns.kdeplot(y_obs[i, :], ax=axes[i], linestyle="--")

    custom_lines = [
        Line2D([0], [0], color=current_palette[i], lw=4)
        for i in range(3)]
    fig.legend(
        custom_lines, ["True", "Estimated", "Observed"], loc=8,
        ncol=3,  fontsize=14, frameon=False)

    for i in range(4):
        axes[i].set_title("Group {0}".format(i + 1))
        axes[i].set_xlim((-5, 5))

    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()


def plot_exp_3(y_ppd, y_obs=None, title="", outfile=""):
    current_palette = sns.color_palette()

    x = np.linspace(-5, 5, num=100)
    fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(16, 6))
    if title:
        fig.suptitle(title, fontsize=18)
    sns.lineplot(x, norm.pdf(x, scale=0.49), ax=axes[0])
    sns.kdeplot(y_ppd[:, 0], ax=axes[0])

    sns.lineplot(x, norm.pdf(x, loc=1, scale=1), ax=axes[1])
    sns.kdeplot(y_ppd[:, 1], ax=axes[1])

    sns.lineplot(
        x,
        0.5*norm.pdf(x, loc=-1.2, scale=0.5) + 0.5*norm.pdf(x, loc=1.2, scale=0.5),
        ax=axes[2])
    sns.kdeplot(y_ppd[:, 2], ax=axes[2])

    sns.lineplot(
        x,
        0.4*norm.pdf(x, loc=-1.2, scale=0.5) + 0.7*norm.pdf(x, loc=1.2, scale=0.5),
        ax=axes[3])
    sns.kdeplot(y_ppd[:, 3], ax=axes[3])

    for i in range(4):
        if y_obs is not None:
            sns.kdeplot(y_obs[i, :], ax=axes[i], linestyle="--")

    custom_lines = [
        Line2D([0], [0], color=current_palette[i], lw=4)
        for i in range(3)]
    fig.legend(
        custom_lines, ["True", "Estimated", "Observed"], loc=8,
        ncol=3,  fontsize=14, frameon=False)

    for i in range(4):
        axes[i].set_title("Group {0}".format(i + 1))
        axes[i].set_xlim((-5, 5))
        axes[i].set_ylim((0, 0.9))

    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()


# def plot_exp_4(y_ppd, title="", outfile=""):
#     current_palette = sns.color_palette()
#
#     x = np.linspace(-5, 5, num=100)
#     fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(16, 6))
#     if title:
#         fig.suptitle(title, fontsize=18)
#
#     for g in range(4):
#         sns.lineplot(
#             x,
#             0.5*norm.pdf(x, loc=-1.2, scale=0.5) + 0.5*norm.pdf(x, loc=1.2, scale=0.5),
#             ax=axes[g])
#         sns.kdeplot(y_ppd[:, g], ax=axes[g])
#
#     custom_lines = [
#         Line2D([0], [0], color=current_palette[i], lw=4)
#         for i in range(2)]
#     fig.legend(
#         custom_lines, ["True", "Estimated"], loc=8,
#         ncol=3,  fontsize=14, frameon=False)
#
#     for i in range(4):
#         axes[i].set_title("Group {0}".format(i + 1))
#         axes[i].set_xlim((-5, 5))
#
#     if outfile:
#         plt.savefig(outfile)
#     else:
#         plt.show()

def plot_exp_4(y_ppd, y_obs=None, title="", outfile=""):
    current_palette = sns.color_palette()

    x = np.linspace(-5, 5, num=100)
    fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(16, 6))
    if title:
        fig.suptitle(title, fontsize=18)

    for g in range(4):
        sns.lineplot(
            x,
            t.pdf(x, df=15),
            ax=axes[g])
        sns.kdeplot(y_ppd[:, g], ax=axes[g])

    for i in range(4):
        if y_obs is not None:
            sns.kdeplot(y_obs[i, :], ax=axes[i], linestyle="--")

    custom_lines = [
        Line2D([0], [0], color=current_palette[i], lw=4)
        for i in range(3)]
    fig.legend(
        custom_lines, ["True", "Estimated", "Observed"], loc=8,
        ncol=3,  fontsize=14, frameon=False)

    for i in range(4):
        axes[i].set_title("Group {0}".format(i + 1))
        axes[i].set_xlim((-5, 5))

    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()


def plot_exp_5(y_ppd, y_obs=None, title="", outfile=""):
    current_palette = sns.color_palette()

    x = np.linspace(-5, 5, num=100)
    fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(16, 6))
    if title:
        fig.suptitle(title, fontsize=18)
    sns.lineplot(x, norm.pdf(x, scale=1), ax=axes[0])
    sns.kdeplot(y_ppd[:, 0], ax=axes[0])

    sns.lineplot(x, norm.pdf(x, loc=0.2, scale=1), ax=axes[1])
    sns.kdeplot(y_ppd[:, 1], ax=axes[1])

    sns.lineplot(
        x,
        norm.pdf(x, loc=0, scale=1.2),
        ax=axes[2])
    sns.kdeplot(y_ppd[:, 2], ax=axes[2])

    sns.lineplot(
        x,
        norm.pdf(x, loc=0, scale=0.25),
        ax=axes[3])
    sns.kdeplot(y_ppd[:, 3], ax=axes[3])

    for i in range(4):
        if y_obs is not None:
            sns.kdeplot(y_obs[i, :], ax=axes[i], linestyle="--")

    custom_lines = [
        Line2D([0], [0], color=current_palette[i], lw=4)
        for i in range(3)]
    fig.legend(
        custom_lines, ["True", "Estimated", "Observed"], loc=8,
        ncol=3,  fontsize=14, frameon=False)

    for i in range(4):
        axes[i].set_title("Group {0}".format(i + 1))
        axes[i].set_xlim((-5, 5))

    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()
