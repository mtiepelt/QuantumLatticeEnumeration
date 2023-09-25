from fpylll.fplll.bkz_param import BKZParam, Strategy, dump_strategies_json

b = 150

# params_fplll = BKZParam(
#     block_size=b,
#     strategies=None,
#     flags=0
#     | fpylll.BKZ.AUTO_ABORT
#     | fpylll.BKZ.MAX_LOOPS,
#     max_loops=20
# )
# strats = params_fplll.strategies


strats = list(map(Strategy, range(4)))
for beta in range(4, b):
    preproc = 45
    if beta <= 10:
        preproc = 3
    elif beta <= 30:
        preproc = 10
    elif beta <= 45:
        preproc = 30
    strats.append(Strategy(beta, [preproc]))

for strat in strats:
    print(strat)

dump_strategies_json("no-pruning.json", strats)

