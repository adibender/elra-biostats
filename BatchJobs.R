cluster.functions <- makeClusterFunctionsSSH(
  makeSSHWorker(nodename="bounty", max.jobs=20, max.load=63, ncpus=64),
  makeSSHWorker(nodename="snickers", max.jobs=20, max.load=63, ncpus=64)
  # makeSSHWorker(nodename="bolt", max.jobs=10, max.load=31, ncpus=32)
  #,makeSSHWorker(nodename="helios", max.jobs=10, max.load=47, ncpus=48)
  )