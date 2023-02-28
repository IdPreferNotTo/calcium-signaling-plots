import os

def write_submit():
    home = os.path.expanduser('~')
    script_name = home + "/SSH BCCN/ca_stochastic.submit"
    with open(script_name, "w") as file:
        file.write("universe = vanilla \n "
                   "notification = Error \n "
                   "Request_memory = 1000 \n "
                   "initialdir=/home/lukasra/CLionProjects/calcium_spikes_markov/cmake-build-release \n "
                   "log = /home/lukasra/condor/log/condor.log \n"
                   "output = /home/lukasra/condor/out/condor.out \n "
                   "error = /home/lukasra/condor/err/condor.err \n")
        N=0
        for i in range(200):
            for j in range(1):
                file.write("executable = /home/lukasra/CLionProjects/calcium_spikes_markov/cmake-build-release/calcium_spikes_markov \n "
                           "arguments = {:d} \n "
                           "requirements = TARGET.machine =!= \"lexa\" \n"
                           "queue \n # "
                           "\n".format(N))
                N += 1
write_submit()