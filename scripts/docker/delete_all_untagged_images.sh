# USE AT YOUR OWN RISK
# from: https://forums.docker.com/t/command-to-remove-all-unused-images/20/2
docker rmi -f $(docker images | grep "<none>" | awk "{print \$3}")