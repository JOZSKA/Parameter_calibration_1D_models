from run import run
from run_seasons import run_seasons
from run_obs_reductions import run_obs_reductions
from run_obs_perturbations import run_obs_perturbations
from configurations import Boussole, BoussoleSatOnly, BATS, BATSSatOnly


def main():
    confs=[
        # Boussole(),
        BoussoleSatOnly(),
        # BATS(),
        # BATSSatOnly(),
        ]
    for conf in confs:
        # run(conf)
        # run_seasons(conf)
        for i in range(1,2):
            run_obs_reductions(conf, obs_ratio=0.5**i)
        # for i in range(1,2):
        run_obs_perturbations(conf)#, obs_ratio=0.5**i)

if __name__=="__main__":
    main()
     
