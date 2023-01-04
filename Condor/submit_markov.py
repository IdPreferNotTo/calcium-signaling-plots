import os

def write_submit():
    home = os.path.expanduser('~')
    script_name = home + "/SSH BCCN/ca_stochastic.submit"
    with open(script_name, "w") as file:
        file.write("universe = vanilla \n "
                   "notification = Error \n "
                   "Request_memory = 500 \n "
                   "initialdir=/home/lukasra/CLionProjects/calcium_spikes_markov/cmake-build-release \n "
                   "log = /home/lukasra/condor_out/condor.log "
                   "\n output = /home/lukasra/condor_out/condor.out \n error = /home/lukasra/condor_out/condor.err \n")
        N=0
        for i in range(41):
            for j in range(41):
                file.write("executable = /home/lukasra/CLionProjects/calcium_spikes_markov/cmake-build-release/calcium_spikes_markov \n arguments = {:d} \n queue \n # \n".format(N))
                N += 1
write_submit()