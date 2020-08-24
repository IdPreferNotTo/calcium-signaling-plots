import os

def write_submit():
    home = os.path.expanduser('~')
    script_name = home + "/BCCN File system/lr_submit.sh"
    with open(script_name, "w") as file:
        file.write("universe = vanilla \n "
                   "notification = Error \n "
                   "Request_memory = 500 \n "
                   "initialdir=/home/lukasra/CLionProjects/li-rinzel-calcium-phd/cmake-build-release \n "
                   "log = /home/lukasra/condor_out/condor.log "
                   "\n output = /home/lukasra/condor_out/condor.out \n error = /home/lukasra/condor_out/condor.err \n")
        for i in range(10):
            for j in range(10):
                N = 2**(i+1)
                ip3 = 0.1*(j+1)
                file.write("executable = /home/lukasra/CLionProjects/li-rinzel-calcium-phd/cmake-build-release/li_rinzel_calcium_phd \n arguments = {:d} {:.2f} \n queue \n # \n".format(N, ip3))

write_submit()