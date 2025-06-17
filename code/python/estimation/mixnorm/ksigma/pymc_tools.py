from pymc.backends.base import MultiTrace
from pymc.stats import rhat
from collections import namedtuple

ChainProgressTimers = namedtuple("ChainProgressTimers", ["tuning_elapsed", 
                                                         "tuning_draws", 
                                                         "elapsed", 
                                                         "draws",
                                                         "avg_acc_rate", 
                                                         "divergences"])

class ClusterFriendlyCallback:
    """print progress so far in NUTs in cluster environments
    Due to the interactive mechanics of pymc default sample behavior, in cluster environments
    it is not clear how many samples were generated and if the process crashed in interim the 
    progress made by the system is not visible. This method will print progress every x draw
    """
    def __init__(self, every=1000, max_rhat=1.05, chains=4):
        self.every = every
        self.max_rhat = max_rhat
        self.traces = {}
        self.timers = {}

    def __call__(self, trace, draw):
        chain_id = draw.chain
        if draw.tuning:
            perf_diffs = trace.get_sampler_stats('perf_counter_diff')
            total_time = sum(perf_diffs)
            
            if len(trace) % self.every == 0:
                # draw_dict = draw._asdict()
                # draw_idx = draw_dict.get("draw_idx", None)
                print(f"Chain {chain_id} tunning. {draw.draw_idx+1} samples so far in {total_time:.2f} seconds")
                
            return
        
        self.traces[draw.chain] = trace
        if len(trace) % self.every == 0:
            perf_diffs = trace.get_sampler_stats('perf_counter_diff')
            total_time = sum(perf_diffs)
            print(f"Chain {chain_id} sampling. {draw.draw_idx+1} samples so far in {total_time:.2f} seconds")
            # print(trace.report)
        
            multitrace = MultiTrace(list(self.traces.values()))
            if rhat(multitrace).to_array().max() < self.max_rhat:
                raise KeyboardInterrupt

