# frenet path planning for f1tenth gym

## Install Gym
```bash
$ git clone https://gitlab.com/acrome-colab/riders-poc/f1tenth-riders-quickstart.git
$ cd f1tenth-riders-quickstart
$ pip install --user -e gym
$ cd pkg/src
$ python -m pkg.main
```
## Install Motion Planner
```bash
$ cd f1tenth-riders-quickstart/pkg/src/pkg
$ git clone https://github.com/Cho-Hyeongrae/f1tenth_overtake_FOT.git
```
## Run a Motion Planner
```bash
# main.py
# import your drivers here
from pkg.drivers import GapFollower
from pkg.drivers import DisparityExtender
from pkg.fot import FrenetPlaner
```

# Vedio
<img width="80%" src="https://github.com/Cho-Hyeongrae/Optimal_frenet_planning_for_f1tenth/blob/main/image/Screen%20Recording%202022-12-03%20at%209.07.49%20PM.gif"/>
