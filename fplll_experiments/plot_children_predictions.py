import os 
import json
import profile
from random import gauss
from matplotlib import legend_handler
from sage.all import line, save, sqrt, log, prod, beta, list_plot, multi_graphics
from bkz_simulators.sim import Sim, LLLProfile

def gaussian_heuristic(v, n):
    """
        :param v:   lattice volume
        :param n:   lattice rank
    """
    from math import pi, sqrt, e
    return ((pi * n)**(1./(2.*n))) * sqrt(n/(2 * pi * e)) * (v ** (1./n))


def lwe_gh(lwe):
    n, q, m = lwe['n'], lwe['q'], lwe['m']
    vol = q**m
    dim = n+m+1
    return gaussian_heuristic(vol, dim)


def profile_gh(profile):
    # profile contains square norms
    dim = len(profile)
    vol = sqrt(prod(profile))
    return gaussian_heuristic(vol, dim)


def gen_plots(data, plots_dir, lwe, bkz, _tour=None, _block=None, dpi=150):
    n, q, m = lwe['n'], lwe['q'], lwe['m']
    sim = Sim(LLLProfile(n, q, m))

    if _tour:
        data_tuple = [data[int(_tour)]]
        tour = int(_tour)
        sim(bkz['beta'], tour)
    else:
        data_tuple = data
        tour = 0

    for tour_stats in data_tuple:
        if _block:
            tour_stats_tuple = [_block]
        else:
            tour_stats_tuple = tour_stats
        for index in tour_stats_tuple:
            print(f"tour: {tour}, index: {index}, beta: {bkz['beta']}")
            sim_block_prof = sim.profile[int(index.split('-')[0]):min(int(index.split('-')[0])+bkz['beta'], len(sim.profile))]
            sim_gh = profile_gh(sim_block_prof)
            sim_R = sqrt(bkz['ghfactor']) * sim_gh
            stats = tour_stats[index]

            stats_keys = set(stats.keys())
            gs2 = stats['gs'][int(index.split('-')[0]):min(int(index.split('-')[0])+bkz['beta'], len(sim.profile))]
            stats_keys.remove('pruning')
            stats_keys.remove('gs')
            stats_keys = list(map(str, sorted(map(int, list(stats_keys)))))

            measured_ln = [ (int(x)+1, stats[x]['mac']) for x in stats_keys ]
            theory_ln = [ (int(x)+1, -1 if stats[x]['tac'] <= 0 else stats[x]['tac']) for x in stats_keys ]
            # theory_ln_lower = [ (x, -1 if stats[x]['tac'] <= 0 else stats[x]['tac'] - sqrt(max(0, -1 if stats[x]['tvc'] <= 0  else stats[x]['tvc']))) for x in stats_keys ]
            # theory_ln_upper = [ (x, -1 if stats[x]['tac'] <= 0 else stats[x]['tac'] + sqrt(max(0, -1 if stats[x]['tvc'] <= 0  else stats[x]['tvc']))) for x in stats_keys ]
            theory_ln_upup = [ (int(x)+1, int(int(x) == 0) + (0.5 if int(x) == 0 else 1) * 2 * sim_R / sqrt(sim.profile[min(int(index.split('-')[0])+bkz['beta'], len(sim.profile))-int(x)-1]) ) for x in stats_keys ] # if k \approx 1, add a (1 if int(x) == 0 else 2) factor
            # theory_ln_upup_half = [ (x, (0.5 if int(x) == 0 else 1) * sim_R / sqrt(sim.profile[min(int(index.split('-')[0])+bkz['beta'], len(sim.profile))-int(x)-1]) ) for x in stats_keys ]
            # theory_ln_upup_pruned = [ (x, (1 if int(x) == 0 else 2)        * sim_R * sqrt(stats['pruning'][-int(x)-1]) / sqrt(sim.profile[min(int(index.split('-')[0])+bkz['beta'], len(sim.profile))-int(x)-1]) ) for x in stats_keys ]
            # theory_ln_upup_half_pruned = [ (x, (0.5 if int(x) == 0 else 1) * sim_R * sqrt(stats['pruning'][-int(x)-1]) / sqrt(sim.profile[min(int(index.split('-')[0])+bkz['beta'], len(sim.profile))-int(x)-1]) ) for x in stats_keys ]
            # no_prun_model = [(x, (0.5 if int(x) == 0 else 1) * 2 * sim_R * sqrt(2/(int(x)+2)) / sqrt(sim_block_prof[len(sim_block_prof) - int(x)-1])) for x in stats_keys]
            no_prun_model_2 = [(int(x)+1, (0.5 if int(x) == 0 else 1) * sim_R * beta(1/2, int(x)/2+1) / sqrt(sim_block_prof[len(sim_block_prof) - int(x)-1])) for x in stats_keys] # if k \approx 1, add a (.5 if int(x) == 0 else 1) factor
            no_prun_model_2_using_gs = [(int(x)+1, (0.5 if int(x) == 0 else 1) * sim_R * beta(1/2, int(x)/2+1) / sqrt(gs2[len(sim_block_prof) - int(x)-1])) for x in stats_keys] # if k \approx 1, add a (.5 if int(x) == 0 else 1) factor
            # no_prun_model_2 = [(x, (0.5 if int(x) == 0 else 1) * 2 * sim_R * beta(1/2, int(x)/2+1)/2 / sqrt(sim_block_prof[len(sim_block_prof) - int(x)-1])) for x in stats_keys]
            # no_prun_model = [(x, (0.5 if int(x) == 0 else 1) * 2 * sim_R * sqrt(2/(int(x)+1))) for x in stats_keys]
            g = line([])
            g += line(no_prun_model_2,          xmin=1, color='black', axes_labels=["$k$", ""], linestyle="--", legend_label="$\\mathbb{E}\\left[C\\left(g \in Z_k\\right)\\right]$ using [CN11]")
            g += line(no_prun_model_2_using_gs, xmin=1, color='red', axes_labels=["$k$", ""], linestyle="--", legend_label="$\\mathbb{E}\\left[C\\left(g \in Z_k\\right)\\right]$ using GS")
            g += list_plot(measured_ln,         xmin=1, color='green', marker="D", size=20, legend_label="Measured $C\\left(g \in Z_k\\right)$ average")
            g += line(theory_ln_upup,           xmin=1, color='blue', linestyle="--", legend_label="$C\\left(g \in Z_k\\right)$ upper bound")
                # + line(theory_ln, color='blue', legend_label="'theory' avg pruning") \
                # + line(theory_ln_lower, color='blue', linestyle='--') \
                # + line(theory_ln_upper, color='blue', linestyle='--') \
                # + line(theory_ln_upup_pruned, color='black', legend_label="upper bnd with pruning") \
                # + line(no_prun_model, color='green', thickness=2, legend_label="model avg, no pruning") \
                # + line(theory_ln_upup_half, color='red', linestyle='--') \
                # + line(theory_ln_upup_half_pruned, color='black', linestyle='--') \
            g.set_legend_options(loc='upper right')
            filename = f"{plots_dir}/{tour}-{index}"

            gs_line = line([(i, log(gs2[i],2)) for i in range(len(gs2))], xmin=1, transparent=False, frame=True, ticks=[[], []], title="basis profile:") # , legend_label="$\log_2(\|b^*_i\|^2$)"
            gs_line += line([(i, log(sim_block_prof[i],2)) for i in range(len(gs2))], xmin=1, linestyle="dashed", color="green", transparent=False, frame=True, ticks=[[], []]) # , legend_label="[CN11]"

            G = multi_graphics([g, (gs_line, (0.65, 0.4, 0.2, 0.2))])

            # save(G, filename+".png", dpi=dpi)
            save(G, filename+".pdf")


        sim(bkz['beta'], 1)
        tour += 1


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-tree-stats-filename', type=str, default=None)
    parser.add_argument('-plot-dir', type=str, default="./plots")
    parser.add_argument('-n', type=int, default=72, help="LWE secret dimension")
    parser.add_argument('-q', type=int, default=97, help="LWE modulo")
    parser.add_argument('-sd', type=float, default=1., help="LWE discrete gaussian standard deviation")
    parser.add_argument('-m', type=int, default=87, help="LWE samples")
    parser.add_argument('-beta', type=int, default=30, help="BKZ block size")
    parser.add_argument('-tours', type=int, default=20, help="Max BKZ tours")
    parser.add_argument('-ghfactor', type=float, default=1, help="BKZ GH factor")
    parser.add_argument('-tour', default=None)
    parser.add_argument('-block', default=None)
    args = parser.parse_args()

    data = json.load(open(args.tree_stats_filename))
    lwe = { 'n': args.n, 'q': args.q, 'sd': args.sd, 'm': args.m }
    bkz = { 'beta': args.beta, 'tours': args.tours, 'ghfactor': args.ghfactor }

    os.system('mkdir -p '+ args.plot_dir)
    gen_plots(data, args.plot_dir, lwe, bkz, _tour=args.tour, _block=args.block)

def sims(lwe, bkz):
    n, q, sd, m = lwe['n'], lwe['q'], lwe['sd'], lwe['m']
    beta, tours, ghfactor = bkz['beta'], bkz['tours'], bkz['ghfactor']
    print(n, q, sd, m, beta, tours, ghfactor)

    lll_prof = LLLProfile(n, q, m)
    sim = Sim(lll_prof)
    save(sim.plot(), "test_lll.png")
    sim(beta, tours)
    save(sim.plot(), "test_bkz.png")
