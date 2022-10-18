


def generate_target_course(x, y):
    csp = cubic_spline_planner.CubicSpline2D(x, y)
    s = np.arange(0, csp.s[-1], 0.1)

    rx, ry, ryaw, rk = [], [], [], []
    for i_s in s:
        ix, iy = csp.calc_position(i_s)
        rx.append(ix)
        ry.append(iy)
        ryaw.append(csp.calc_yaw(i_s))
        rk.append(csp.calc_curvature(i_s))

    return rx, ry, ryaw, rk, csp

class SimpleDriver:

    def __init__(self):
        # read waypoint

        # input waypoint -> 
        self.tx, self.ty, self.tyaw, self.tc, self.csp = generate_target_course(wx, wy)


        self.current_speed = 10.0 / 3.6 # 현재 속도 (c_speed)
        self.current_lateral_position = 2.0 # 현재 측면 위치 (c_d)
        self.current_lateral_speed = 0.0 # 현재 측면 속도 (c_d_d)
        self.current_lateral_acceleration = 0.0 # 현재 횡가속도 (c_d_dd)
        self.current_course_position = 0.0 # 현재 코스 위치 (s0)


    def frenet_optimal_planning(self, csp, s0, c_speed, c_d, c_d_d, c_d_dd, ob):
        fplist = calc_frenet_paths(c_speed, c_d, c_d_d, c_d_dd, s0)
        fplist = calc_global_paths(fplist, csp)
        fplist = check_paths(fplist, ob)

        # find minimum cost path
        min_cost = float("inf")
        best_path = None
        for fp in fplist:
            if min_cost >= fp.cf:
                min_cost = fp.cf
                best_path = fp

        return best_path
    def process_observation(self, ranges, ego_odom):

        # obstacle 위치 업데이트
        ## code
        
        # c_speed, c_d, c_d_d, c_d_dd, s0 업데이트
        path = self.frenet_optimal_planning(self.csp, 
        self.current_course_position, self.current_speed, 
        self.current_lateral_position, self.current_lateral_speed, 
        self.current_lateral_acceleration, ob)

        self.current_course_position = path.s[1] # traget wp x
        self.current_lateral_position = path.d[1] 
        self.current_lateral_speed = path.d_d[1]
        self.current_lateral_acceleration = path.d_dd[1]
        self.current_speed = path.s_d[1] # = speed
        # path.x[1] , path.y[1] : 현재 자기 위치
        # tx[-1], ty[-1] : 도착 지점 위치



        return speed, steering_angle